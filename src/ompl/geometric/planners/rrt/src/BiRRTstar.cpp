/*********************************************************************
* Software License Agreement (BSD License)
*
*  Copyright (c) 2011, Rice University
*  All rights reserved.
*
*  Redistribution and use in source and binary forms, with or without
*  modification, are permitted provided that the following conditions
*  are met:
*
*   * Redistributions of source code must retain the above copyright
*     notice, this list of conditions and the following disclaimer.
*   * Redistributions in binary form must reproduce the above
*     copyright notice, this list of conditions and the following
*     disclaimer in the documentation and/or other materials provided
*     with the distribution.
*   * Neither the name of the Rice University nor the names of its
*     contributors may be used to endorse or promote products derived
*     from this software without specific prior written permission.
*
*  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
*  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
*  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
*  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
*  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
*  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
*  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
*  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
*  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
*  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
*  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
*  POSSIBILITY OF SUCH DAMAGE.
*********************************************************************/

/* Authors: Shi Shenglei */

#include "ompl/geometric/planners/rrt/BiRRTstar.h"
#include "ompl/base/goals/GoalSampleableRegion.h"
#include "ompl/base/objectives/PathLengthOptimizationObjective.h"
#include "ompl/tools/config/SelfConfig.h"
#include <boost/math/constants/constants.hpp>

ompl::geometric::BiRRTstar::BiRRTstar(const base::SpaceInformationPtr &si) : base::Planner(si, "BiRRTstar")
{
    specs_.approximateSolutions = false;
    specs_.optimizingPaths = true;
    specs_.canReportIntermediateSolutions = true;

    Planner::declareParam<double>("range", this, &BiRRTstar::setRange, &BiRRTstar::getRange, "0.:1.:10000.");
    Planner::declareParam<double>("collision_range", this, &BiRRTstar::setMaxCollisionDistance, &BiRRTstar::getMaxCollisionDistance, "0.:1.:10000.");
    Planner::declareParam<double>("rewire_factor", this, &BiRRTstar::setRewireFactor, &BiRRTstar::getRewireFactor,
                                  "1.0:0.01:2.0");

    addPlannerProgressProperty("best cost REAL", [this] { return bestCostProperty(); });
}

ompl::geometric::BiRRTstar::~BiRRTstar()
{
    freeMemory();
}

void ompl::geometric::BiRRTstar::setup()
{
    Planner::setup();
    tools::SelfConfig sc(si_, getName());
    sc.configurePlannerRange(maxDistance_);
    sc.configurePlannerCollisionRange(maxCollisionDistance_);
    if (!si_->getStateSpace()->hasSymmetricDistance() || !si_->getStateSpace()->hasSymmetricInterpolate())
    {
        OMPL_WARN("%s requires a state space with symmetric distance and symmetric interpolation.", getName().c_str());
    }

    if (!tStart_)
        tStart_.reset(tools::SelfConfig::getDefaultNearestNeighbors<Motion *>(this));

    if (!tGoal_)
        tGoal_.reset(tools::SelfConfig::getDefaultNearestNeighbors<Motion *>(this));

    tStart_->setDistanceFunction([this](const Motion *a, const Motion *b)
                                 {
                                    return distanceFunction(a, b); 
                                 });
    tGoal_->setDistanceFunction([this](const Motion *a, const Motion *b)
                                {
                                    return distanceFunction(a, b); 
                                });

    // Setup optimization objective
    //
    // If no optimization objective was specified, then default to
    // optimizing path length as computed by the distance() function
    // in the state space.
    if (pdef_)
    {
        if (pdef_->hasOptimizationObjective())
            opt_ = pdef_->getOptimizationObjective();
        else
        {
            OMPL_INFORM("%s: No optimization objective specified. Defaulting to optimizing path length for the allowed "
                        "planning time.",
                        getName().c_str());
            opt_ = std::make_shared<base::PathLengthOptimizationObjective>(si_);

            // Store the new objective in the problem def'n
            pdef_->setOptimizationObjective(opt_);
        }

        // Set the bestCost_ and prunedCost_ as infinite
        bestCost_ = opt_->infiniteCost();
    }
    else
    {
        OMPL_INFORM("%s: problem definition is not set, deferring setup completion...", getName().c_str());
        setup_ = false;
    }

    // Calculate some constants:
    calculateRewiringLowerBounds();
}

void ompl::geometric::BiRRTstar::freeMemory()
{
    std::vector<Motion *> motions;
    if (tStart_)
    {
        tStart_->list(motions);
        for (auto &motion : motions)
        {
            if (motion->state != nullptr)
                si_->freeState(motion->state);
            delete motion;
        }
        motions.clear();
    }

    if (tGoal_)
    {
        tGoal_->list(motions);
        for (auto &motion : motions)
        {
            if (motion->state != nullptr)
                si_->freeState(motion->state);
            delete motion;
        }
        motions.clear();
    }
}

void ompl::geometric::BiRRTstar::clear()
{
    setup_ = false;
    Planner::clear();
    sampler_.reset();
    freeMemory();
    if (tStart_)
        tStart_->clear();
    if (tGoal_)
        tGoal_->clear();

    startMotions_.clear();
    goalMotions_.clear();

    bestCost_ = base::Cost(std::numeric_limits<double>::quiet_NaN());

    connectionPoint_.clear();
    bestStartMotion_ = nullptr;
    bestGoalMotion_ = nullptr;
}

ompl::base::PlannerStatus ompl::geometric::BiRRTstar::solve(const base::PlannerTerminationCondition &ptc)
{
    checkValidity();
    auto *goal = dynamic_cast<base::GoalSampleableRegion *>(pdef_->getGoal().get());

    if (goal == nullptr)
    {
        OMPL_ERROR("%s: Unknown type of goal", getName().c_str());
        return base::PlannerStatus::UNRECOGNIZED_GOAL_TYPE;
    }
    else if (!goal->couldSample())
    {
        OMPL_ERROR("%s: Insufficient states in sampleable goal region", getName().c_str());
        return base::PlannerStatus::INVALID_GOAL;
    }
    // Check if there are more starts
    if (pis_.haveMoreStartStates() == true)
    {
        // There are, add them
        while (const base::State *st = pis_.nextStart())
        {
            auto *motion = new Motion(si_);
            si_->copyState(motion->state, st);
            motion->root = motion->state;
            motion->cost = opt_->identityCost();
            tStart_->add(motion);
            startMotions_.push_back(motion);
        }
    }

    if (tStart_->size() == 0)
    {
        OMPL_ERROR("%s: There are no valid initial states!", getName().c_str());
        return base::PlannerStatus::INVALID_START;
    }

    sampler_ = si_->allocStateSampler();

    OMPL_INFORM("%s: Started planning with %u states. Seeking a solution better than %.5f.", getName().c_str(), tStart_->size(), opt_->getCostThreshold().value());

    const base::ReportIntermediateSolutionFn intermediateSolutionCallback = pdef_->getIntermediateSolutionCallback();

    TreeGrowingInfo tgi;
    tgi.xstate = si_->allocState();

    auto *rmotion = new Motion(si_);
    base::State *rstate = rmotion->state;

    bool optimal = false;
    bool startTree = true;

    OMPL_INFORM("%s: Initial k-nearest value of %u", getName().c_str(), (unsigned int)std::ceil(k_rrt_ * log((double)(tStart_->size() + 1u))));

    while (ptc == false)
    {
        if (tGoal_->size() == 0 || pis_.getSampledGoalsCount() < tGoal_->size() / 2)
        {
            const base::State *st = tGoal_->size() == 0 ? pis_.nextGoal(ptc) : pis_.nextGoal();
            if (st != nullptr)
            {
                auto *motion = new Motion(si_);
                si_->copyState(motion->state, st);
                motion->root = motion->state;
                motion->cost = opt_->identityCost();
                tGoal_->add(motion);
                goalMotions_.push_back(motion);
            }
            if (tGoal_->size() == 0)
            {
                OMPL_ERROR("%s: Unable to sample any valid states for goal tree", getName().c_str());
                break;
            }
        }

        TreeData &tree = startTree ? tStart_ : tGoal_;
        tgi.start = startTree;
        startTree = !startTree;
        TreeData &otherTree = startTree ? tStart_ : tGoal_;

        sampler_->sampleUniform(rstate);

        bool gs = growTree(tree, tgi, rmotion);
        /* remember which motion was just added */
        Motion *addedMotion = tgi.xmotion;
        if (!gs)
            si_->copyState(rstate, addedMotion->state);

        tgi.start = startTree;
        bool gsc = growTree(otherTree, tgi, rmotion);
        if (gsc)
        {
            Motion *startMotion = startTree ? tgi.xmotion : addedMotion;
            Motion *goalMotion = startTree ? addedMotion : tgi.xmotion;
            connectionPoint_.emplace_back(startMotion, goalMotion);
        }

        bool updatedSolution = false;
        for (auto & pair : connectionPoint_)
        {
            base::Cost temp = opt_->combineCosts(pair.first->cost, pair.second->cost);
            if (opt_->isCostBetterThan(temp, bestCost_))
            {
                updatedSolution = true;
                bestCost_ = temp;
                OMPL_INFORM("%s: Found a better solution with a cost of %.2f (%u vertices in the graph)",
                        getName().c_str(), bestCost_.value(), tStart_->size() + tGoal_->size());
                bestStartMotion_ = pair.first;
                bestGoalMotion_ = pair.second;
                if (opt_->isSatisfied(bestCost_))
                {
                    optimal = true;
                    break;
                }
            }
        }
        if (optimal)
            break;
        if (updatedSolution)
        {
            if (intermediateSolutionCallback)
            {
                std::vector<const base::State *> spath;

                bool del = false;
                const Motion *solution = bestStartMotion_;
                if (solution->parent != nullptr)
                {
                    del = true;
                    solution = solution->parent;
                }

                std::vector<const base::State *> mpath1;
                while (solution != nullptr)
                {
                    mpath1.push_back(solution->state);
                    solution = solution->parent;
                }

                for (std::size_t i = mpath1.size() - 1; i < mpath1.size(); --i)
                    spath.push_back(mpath1[i]);

                solution = bestGoalMotion_;
                if (!del)
                    solution = solution->parent;
                while (solution != nullptr)
                {
                    spath.push_back(solution->state);
                    solution = solution->parent;
                }

                intermediateSolutionCallback(this, spath, bestCost_);
            }
        }
    }

    if (bestStartMotion_)
    {
        std::vector<const base::State *> spath;

        bool del = false;
        const Motion *solution = bestStartMotion_;
        if (solution->parent != nullptr)
        {
            del = true;
            solution = solution->parent;
        }
        std::vector<const base::State *> mpath1;
        while (solution != nullptr)
        {
            mpath1.push_back(solution->state);
            solution = solution->parent;
        }

        for (std::size_t i = mpath1.size() - 1; i < mpath1.size(); i--)
            spath.push_back(mpath1[i]);

        solution = bestGoalMotion_;
        if (!del)
            solution = solution->parent;
        while (solution != nullptr)
        {
            spath.push_back(solution->state);
            solution = solution->parent;
        }

        auto path(std::make_shared<PathGeometric>(si_));
        path->getStates().reserve(spath.size());
        for (auto & i : spath)
            path->append(i);

        // Add the solution path.
        base::PlannerSolution psol(path);
        psol.setPlannerName(getName());

        psol.setOptimized(opt_, bestCost_, opt_->isSatisfied(bestCost_));
        pdef_->addSolutionPath(psol);
    }

    si_->freeState(tgi.xstate);
    si_->freeState(rstate);
    delete rmotion;

    OMPL_INFORM("%s: Created %u new states. Final solution cost %.3f", getName().c_str(), tStart_->size() + tGoal_->size(), bestCost_.value());
    return bestStartMotion_ ? base::PlannerStatus::EXACT_SOLUTION : base::PlannerStatus::TIMEOUT;
}

bool ompl::geometric::BiRRTstar::growTree(TreeData &tree, TreeGrowingInfo &tgi, Motion *rmotion)
{
    Motion *nmotion = tree->nearest(rmotion);
    tgi.xmotion = nmotion;
    if (si_->equalStates(nmotion->state, rmotion->state))
        return false;

    bool symCost = opt_->isSymmetric();

    std::vector<Motion *> nbh;
    std::vector<int> valid;
    std::vector<base::Cost> costs;
    std::vector<base::Cost> incCosts;
    std::vector<std::size_t> sortedCostIndices;
    CostIndexCompare compareFn(costs, *opt_);

    bool reach = false;
    while (!reach)
    {
        nmotion = tgi.xmotion;

        base::State *dstate = rmotion->state;
        double d = tgi.start ? si_->distance(nmotion->state, rmotion->state) : si_->distance(rmotion->state, nmotion->state);
        if (d > maxCollisionDistance_)
        {
            if (tgi.start)
                si_->getStateSpace()->interpolate(nmotion->state, rmotion->state, maxCollisionDistance_ / d, tgi.xstate);
            else 
                si_->getStateSpace()->interpolate(rmotion->state, nmotion->state, 1.0 - maxCollisionDistance_ / d, tgi.xstate);
            if (si_->equalStates(nmotion->state, tgi.xstate))
                break;
            dstate = tgi.xstate;
        }
        else 
            reach = true;

        if (!si_->checkMotion(nmotion->state, dstate))
        {
            reach = false;
            break;
        }

        Motion *motion = new Motion(si_);
        si_->copyState(motion->state, dstate);
        motion->parent = nmotion;
        motion->root = motion->parent->root;
        motion->incCost = opt_->motionCost(nmotion->state, motion->state);
        motion->cost = opt_->combineCosts(nmotion->cost, motion->incCost);
      
        getNeighbors(tree, motion, nbh);

        if (costs.size() < nbh.size())
        {
            costs.resize(nbh.size());
            incCosts.resize(nbh.size());
            sortedCostIndices.resize(nbh.size());
            valid.resize(nbh.size());
        }
        std::fill(valid.begin(), valid.begin() + nbh.size(), 0);

        for (std::size_t i = 0; i < nbh.size(); ++i)
        {
            incCosts[i] = opt_->motionCost(nbh[i]->state, motion->state);
            costs[i] = opt_->combineCosts(nbh[i]->cost, incCosts[i]);
        }

        for (std::size_t i = 0; i < nbh.size(); ++i)
            sortedCostIndices[i] = i;
        std::sort(sortedCostIndices.begin(), sortedCostIndices.begin() + nbh.size(), compareFn);

        for (std::vector<std::size_t>::const_iterator i = sortedCostIndices.begin(); i != sortedCostIndices.begin() + nbh.size(); ++i)
        {
            if (nbh[*i] == nmotion)
            {
                valid[*i] = 1;
                break;
            } 
            if (si_->distance(nbh[*i]->state, motion->state) < maxDistance_)
            {
                if (si_->checkMotion(nbh[*i]->state, motion->state, true))
                {
                    motion->incCost = incCosts[*i];
                    motion->cost = costs[*i];
                    motion->parent = nbh[*i];
                    motion->root = motion->parent->root;
                    break;
                }
                else
                    valid[*i] = -1;
            }
            else
                valid[*i] = -1;
        }

        motion->parent->children.push_back(motion);
        tree->add(motion);
        tgi.xmotion = motion;

        for (std::size_t i = 0; i < nbh.size(); ++i)
        {
            if (nbh[i] != motion->parent)
            {
                base::Cost nbhIncCost;
                if (symCost)
                    nbhIncCost = incCosts[i];
                else
                    nbhIncCost = opt_->motionCost(motion->state, nbh[i]->state);
                base::Cost nbhNewCost = opt_->combineCosts(motion->cost, nbhIncCost);
                if (opt_->isCostBetterThan(nbhNewCost, nbh[i]->cost))
                {
                    bool motionValid = false;
                    if (valid[i] == 0)
                    {
                        motionValid = si_->distance(motion->state, nbh[i]->state) < maxDistance_;
                        if (motionValid)
                            motionValid = si_->checkMotion(motion->state, nbh[i]->state, true);
                    }
                    else
                        motionValid = (valid[i] == 1);

                    if (motionValid)
                    {
                        removeFromParent(nbh[i]);

                        nbh[i]->parent = motion;
                        nbh[i]->root   = nbh[i]->parent->root;
                        nbh[i]->incCost = nbhIncCost;
                        nbh[i]->cost = nbhNewCost;
                        nbh[i]->parent->children.push_back(nbh[i]);

                        // Update the costs of the node's children
                        updateChildCosts(nbh[i]);
                    }
                }
            }
        }
    }

    return reach;
}

void ompl::geometric::BiRRTstar::getNeighbors(const TreeData tree, Motion *motion, std::vector<Motion *> &nbh) const
{
    auto cardDbl = static_cast<double>(tree->size() + 1u);
    unsigned int k = std::ceil(k_rrt_ * log(cardDbl));
    tree->nearestK(motion, k, nbh);
}

void ompl::geometric::BiRRTstar::removeFromParent(Motion *m)
{
    for (auto it = m->parent->children.begin(); it != m->parent->children.end(); ++it)
    {
        if (*it == m)
        {
            m->parent->children.erase(it);
            break;
        }
    }
}

void ompl::geometric::BiRRTstar::updateChildCosts(Motion *m)
{
    for (std::size_t i = 0; i < m->children.size(); ++i)
    {
        m->children[i]->cost = opt_->combineCosts(m->cost, m->children[i]->incCost);
        updateChildCosts(m->children[i]);
    }
}

void ompl::geometric::BiRRTstar::getPlannerData(base::PlannerData &data) const
{
    Planner::getPlannerData(data);
    std::vector<Motion *> motions;
    if (tStart_)
        tStart_->list(motions);
    for (auto &motion : motions)
    {
        if (!motion->parent)
            data.addStartVertex(base::PlannerDataVertex(motion->state, 1));
        else
        {
            data.addEdge(base::PlannerDataVertex(motion->parent->state, 1), base::PlannerDataVertex(motion->state, 1));
        }
    }

    motions.clear();
    if (tGoal_)
        tGoal_->list(motions);
    for (auto &motion : motions)
    {
        if (!motion->parent)
            data.addGoalVertex(base::PlannerDataVertex(motion->state, 2));
        else
        {
            // The edges in the goal tree are reversed to be consistent with start tree
            data.addEdge(base::PlannerDataVertex(motion->state, 2), base::PlannerDataVertex(motion->parent->state, 2));
        }
    }
}

void ompl::geometric::BiRRTstar::calculateRewiringLowerBounds()
{
    const auto dimDbl = static_cast<double>(si_->getStateDimension());
    // k_rrt > 2^(d + 1) * e * (1 + 1 / d).  K-nearest RRT*
    k_rrt_ = rewireFactor_ * (std::pow(2.0, dimDbl + 1.0) * boost::math::constants::e<double>() * (1.0 + 1.0 / dimDbl));
}
