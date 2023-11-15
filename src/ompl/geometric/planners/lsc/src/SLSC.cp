/*********************************************************************
 * Software License Agreement (BSD License)
 *
 *  Copyright (c) 2008, Willow Garage, Inc.
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
 *   * Neither the name of the Willow Garage nor the names of its
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

/* Author: Shi Shenglei */

#include "ompl/geometric/planners/lsc/SLSC.h"
#include "ompl/base/goals/GoalSampleableRegion.h"
#include "ompl/tools/config/SelfConfig.h"

#include <limits>
#include <cassert>
#include <vector>

ompl::geometric::SLSC::SLSC(const base::SpaceInformationPtr &si, const SafetyCertificateChecker &safetyCertificateChecker,
                            const CollisionCertificateChecker &collisionCertificateChecker, double confidence) :
    base::Planner(si, "SLSC"), geometric::CollisionChecking(si, safetyCertificateChecker, collisionCertificateChecker, confidence) 
{
    specs_.recognizedGoal = base::GOAL_SAMPLEABLE_REGION;

    Planner::declareParam<double>("range", this, &SLSC::setRange, &SLSC::getRange, "0.:1.:10000.");
	Planner::declareParam<unsigned int>("collision_confidence", this, &SLSC::setCollisionConfidence, &SLSC::getCollisionConfidence, "0:1:10");
	Planner::declareParam<unsigned int>("safety_confidence", this, &SLSC::setSafetyConfidence, &SLSC::getSafetyConfidence, "0:1:10");

    connectionPoint_ = std::make_pair<base::State *, base::State *>(nullptr, nullptr);
    distanceBetweenTrees_ = std::numeric_limits<double>::infinity();
}

ompl::geometric::SLSC::~SLSC()
{
    freeMemory();
}

void ompl::geometric::SLSC::setup()
{
    Planner::setup();
    tools::SelfConfig sc(si_, getName());

    symmetric_ = si_->getStateSpace()->hasSymmetricDistance() && si_->getStateSpace()->hasSymmetricInterpolate();

    sc.configureProjectionEvaluator(projectionEvaluator_);
    sc.configurePlannerRange(maxDistance_);

    tStart_.grid.setDimension(projectionEvaluator_->getDimension());
    tGoal_.grid.setDimension(projectionEvaluator_->getDimension());

    tStart_.grid.setDistanceFunction([this](const Motion *a, const Motion *b)
                                         {
                                            return distanceFunction(a, b); 
                                         });
    tGoal_.grid.setDistanceFunction([this](const Motion *a, const Motion *b)
                                        {
                                            return distanceFunction(a, b); 
                                        });

	if (!snn_)
        snn_.reset(tools::SelfConfig::getDefaultNearestNeighbors<base::SafetyCertificate *>(this));
    snn_->setDistanceFunction([this](const base::SafetyCertificate *a, const base::SafetyCertificate *b)
                             {
                                 return distanceFunction(a, b);
                             });
	
	if (!onn_)
        onn_.reset(tools::SelfConfig::getDefaultNearestNeighbors<base::SafetyCertificate *>(this));
    onn_->setDistanceFunction([this](const base::SafetyCertificate *a, const base::SafetyCertificate *b)
                             {
                                 return distanceFunction(a, b);
                             });
}

void ompl::geometric::SLSC::freeMemory()
{
    freeGridMotions(tStart_.grid);
    freeGridMotions(tGoal_.grid);
}

void ompl::geometric::SLSC::freeGridMotions(GridNearest<MotionInfo, Motion*> &grid)
{
    for (const auto &it : grid)
    {
        for (unsigned int i = 0; i < it.second->data.size(); ++i)
        {
            if (it.second->data[i]->state)
                si_->freeState(it.second->data[i]->state);
            delete it.second->data[i];
        }
    }
}

void ompl::geometric::SLSC::clear()
{
    Planner::clear();
    CollisionChecking::clear();

    sampler_.reset();

    freeMemory();

    tStart_.grid.clear();
    tStart_.size = 0;
    tStart_.pdf.clear();

    tGoal_.grid.clear();
    tGoal_.size = 0;
    tGoal_.pdf.clear();

    tStart_.minWeight = 1.0;
    tStart_.maxWeight = 0.0;

    tGoal_.minWeight = 1.0;
    tGoal_.maxWeight = 0.0;

    connectionPoint_ = std::make_pair<base::State *, base::State *>(nullptr, nullptr);
    distanceBetweenTrees_ = std::numeric_limits<double>::infinity();
}

ompl::base::PlannerStatus ompl::geometric::SLSC::solve(const base::PlannerTerminationCondition &ptc)
{
    checkValidity();

    auto *goal = dynamic_cast<base::GoalSampleableRegion *>(pdef_->getGoal().get());

    if (!goal)
    {
        OMPL_ERROR("%s: Unknown type of goal", getName().c_str());
        return base::PlannerStatus::UNRECOGNIZED_GOAL_TYPE;
    }
    else if (!goal->couldSample())
    {
        OMPL_ERROR("%s: Insufficient states in sampleable goal region", getName().c_str());
        return base::PlannerStatus::INVALID_GOAL;
    }

    if (pis_.haveMoreStartStates())
    {
        while (const base::State *st = pis_.nextStart())
        {
            auto *motion = new Motion(si_);
            si_->copyState(motion->state, st);
            motion->valid = true;
            motion->stateValid = true;
            motion->root = motion->state;
            addMotion(tStart_, motion);
        }
    }

    if (tStart_.size == 0)
    {
        OMPL_ERROR("%s: Motion planning start tree could not be initialized!", getName().c_str());
        return base::PlannerStatus::INVALID_START;
    }

    if (!sampler_)
        sampler_ = si_->allocStateSampler();

    OMPL_INFORM("%s: Starting planning with %d states already in datastructure", getName().c_str(),
                (int)(tStart_.size + tGoal_.size));

    TreeGrowingInfo tgi;
    tgi.xstate = si_->allocState();

    auto *rmotion = new Motion(si_);
    base::State *rstate = rmotion->state;

    bool startTree = true;
    bool solved = false;

    base::State *s1 = si_->allocState();
    base::State *s2 = si_->allocState() ; 

    Motion *startMotion = nullptr, *goalMotion = nullptr;

    while (!ptc)
    {
        TreeData &tree = startTree ? tStart_ : tGoal_;
        tgi.start = startTree;
        startTree = !startTree;
        TreeData &otherTree = startTree ? tStart_ : tGoal_;

        // if we have not sampled too many goals already
        if (tGoal_.size == 0 || pis_.getSampledGoalsCount() < tGoal_.size / 2)
        {
            const base::State *st = tGoal_.size == 0 ? pis_.nextGoal(ptc) : pis_.nextGoal();
            if (st != nullptr)
            {
                auto *motion = new Motion(si_);
                si_->copyState(motion->state, st);
                motion->valid = true;
                motion->stateValid = true;
                motion->root = motion->state;
                addMotion(tGoal_, motion);
            }

            if (tGoal_.size == 0)
            {
                OMPL_ERROR("%s: Unable to sample any valid states for goal tree", getName().c_str());
                break;
            }
        }

        /* sample a random state todo*/
        double r = rng_.uniform01();
        if (r > tree.maxWeight || r < tree.minWeight)
        {
            Motion *existing = selectMotion(tree); 
            assert(existing);
            sampler_->sampleUniformNear(rstate, existing->state, maxDistance_);
        }
        else if (!validStates_.empty() && rng_.uniform01() < 0.5)
        {
            si_->copyState(rstate, validStates_.back());
            si_->freeState(validStates_.back());
            validStates_.pop_back();
        }
        else 
            sampler_->sampleUniform(rstate);

        GrowState gs = growTree(tree, tgi, rmotion);

        if (gs != TRAPPED)
        {
            /* remember which motion was just added */
            Motion *addedMotion = tgi.xmotion;

            /* attempt to connect trees */

            /* if reached, it means we used rstate directly, no need to copy again */
            if (gs != REACHED)
                si_->copyState(rstate, addedMotion->state);

            GrowState gsc = ADVANCED;
            tgi.start = startTree;

            while (gsc == ADVANCED)
                gsc = growTree(otherTree, tgi, rmotion);

            /* update distance between trees */
            const double newDist = !startTree ? tree.grid.getDistanceFunction()(addedMotion, tgi.xmotion) :
                                                tree.grid.getDistanceFunction()(tgi.xmotion, addedMotion);
            if (newDist < distanceBetweenTrees_) 
                distanceBetweenTrees_ = newDist;

            /* if we connected the trees in a valid way (start and goal pair is valid)*/
            if (gsc == REACHED)
            {
                startMotion = startTree ? tgi.xmotion : addedMotion;
                goalMotion = startTree ? addedMotion : tgi.xmotion;

                if (startMotion->parent != nullptr)
                    startMotion = startMotion->parent;
                else
                    goalMotion = goalMotion->parent;

                if (goal->isStartGoalPairValid(startMotion->root, goalMotion->root))
                {
                    if (!startTree ? isPathValid(tree, otherTree, addedMotion, tgi.xmotion) : 
                                     isPathValid(otherTree, tree, tgi.xmotion, addedMotion))
                    {
                        connectionPoint_ = std::make_pair(startMotion->state, goalMotion->state);

                        /* construct the solution path */
                        Motion *solution = startMotion;
                        std::vector<Motion *> mpath1;
                        while (solution != nullptr)
                        {
                            mpath1.push_back(solution);
                            solution = solution->parent;
                        }

                        solution = goalMotion;
                        std::vector<Motion *> mpath2;
                        while (solution != nullptr)
                        {
                            mpath2.push_back(solution);
                            solution = solution->parent;
                        }

                        auto path(std::make_shared<PathGeometric>(si_));
                        path->getStates().reserve(mpath1.size() + mpath2.size());
                        for (int i = mpath1.size() - 1; i >= 0; --i)
                            path->append(mpath1[i]->state);
                        for (auto &i : mpath2)
                            path->append(i->state);

                        pdef_->addSolutionPath(path, false, 0.0, getName());
                        solved = true;
                        break;
                    }
                    else 
                        distanceBetweenTrees_ = std::numeric_limits<double>::infinity();
                }
            }
        }
 
        if (tGoal_.size == 0 || solved)
        {
            break;
        }
    }

    si_->freeState(tgi.xstate);
    si_->freeState(rstate);
    delete rmotion;

    si_->freeState(s1);
    si_->freeState(s2);


    OMPL_INFORM("%s: Created %u (%u start + %u goal) states in %u cells (%u start + %u goal)", getName().c_str(),
                tStart_.size + tGoal_.size, tStart_.size, tGoal_.size, tStart_.grid.size() + tGoal_.grid.size(),
                tStart_.grid.size(), tGoal_.grid.size());

    return solved ? base::PlannerStatus::EXACT_SOLUTION : base::PlannerStatus::TIMEOUT;
}

ompl::geometric::SLSC::GrowState ompl::geometric::SLSC::growTree(TreeData &tree, TreeGrowingInfo &tgi, Motion *rmotion)
{
    /* find closest state in the tree */
    Motion *nmotion = tree.grid.nearest(rmotion);

    /* assume we can reach the state we go towards */
    bool reach = true;

    /* find state to add */
    base::State *dstate = rmotion->state;
    double d = tgi.start ? si_->distance(nmotion->state, rmotion->state) : si_->distance(rmotion->state, nmotion->state);
    if (d > maxDistance_)
    {
        if (tgi.start)
            si_->getStateSpace()->interpolate(nmotion->state, rmotion->state, maxDistance_ / d, tgi.xstate);
        else 
            si_->getStateSpace()->interpolate(rmotion->state, nmotion->state, 1.0 - maxDistance_ / d, tgi.xstate);

        /* Check if we have moved at all. Due to some stranger state spaces (e.g., the constrained state spaces),
         * interpolate can fail and no progress is made. Without this check, the algorithm gets stuck in a loop as it
         * thinks it is making progress, when none is actually occurring. */
        if (si_->equalStates(nmotion->state, tgi.xstate))
        {
            tgi.xmotion = nmotion;
            return TRAPPED;
        }

        dstate = tgi.xstate;
        reach = false;
    }

    bool trueValid = false;
    bool trueValidState = false;

    if (isValid(dstate, trueValidState))
    {
        auto *motion = new Motion(si_);
        si_->copyState(motion->state, dstate);
        motion->parent = nmotion;
        motion->root = nmotion->root;
        motion->parent->children.push_back(motion);
        motion->stateValid = trueValidState;
        addMotion(tree, motion);

        tgi.xmotion = motion;
    }
    else 
    {
        tgi.xmotion = nmotion;
        return TRAPPED;
    }

    while(!reach)
    {
        nmotion = tgi.xmotion;

        dstate = rmotion->state;
        double d = tgi.start ? si_->distance(nmotion->state, rmotion->state) : si_->distance(rmotion->state, nmotion->state);
        if (d > maxDistance_)
        {
            if (tgi.start)
                si_->getStateSpace()->interpolate(nmotion->state, rmotion->state, maxDistance_ / d, tgi.xstate);
            else 
                si_->getStateSpace()->interpolate(rmotion->state, nmotion->state, 1.0 - maxDistance_ / d, tgi.xstate);

            /* Check if we have moved at all. Due to some stranger state spaces (e.g., the constrained state spaces),
             * interpolate can fail and no progress is made. Without this check, the algorithm gets stuck in a loop as it
             * thinks it is making progress, when none is actually occurring. */
            if (si_->equalStates(nmotion->state, tgi.xstate))
                break;

            dstate = tgi.xstate;
        }
        else 
            reach = true;

        trueValid = false;
        trueValidState = false;

        std::vector<int> stepValid;

        if (rng_.uniform01() < 0.5)
        { 
            if (tgi.start && checkStartMotion(nmotion->state, dstate, trueValid, trueValidState, stepValid)) 
            {
                auto *motion = new Motion(si_);
                si_->copyState(motion->state, dstate);
                motion->parent = nmotion;
                motion->root = nmotion->root;
                motion->valid= trueValid;
                motion->stateValid = trueValidState;
                if (!trueValid)
                    motion->stepValid.swap(stepValid);
                motion->parent->children.push_back(motion);
                addMotion(tree, motion);

                tgi.xmotion = motion;
            }
            else if (!tgi.start && checkGoalMotion(dstate, nmotion->state, trueValid, trueValidState, stepValid)) 
            {
                auto *motion = new Motion(si_);
                si_->copyState(motion->state, dstate);
                motion->parent = nmotion;
                motion->root = nmotion->root;
                motion->valid= trueValid;
                motion->stateValid = trueValidState;
                if (!trueValid)
                    motion->stepValid.swap(stepValid);
                motion->parent->children.push_back(motion);
                addMotion(tree, motion);

                tgi.xmotion = motion;
            }
            else
            {
                reach = false;
                break;
            }
        }
        else 
        {
            if (isValid(dstate, trueValidState))
            {
                auto *motion = new Motion(si_);
                si_->copyState(motion->state, dstate);
                motion->parent = nmotion;
                motion->root = nmotion->root;
                motion->parent->children.push_back(motion);
                motion->stateValid = trueValidState;
                addMotion(tree, motion);

                tgi.xmotion = motion;
            }
            else 
            {
                reach = false;
                break;
            }
        }
    }

    return reach ? REACHED : ADVANCED;
}

ompl::geometric::SLSC::Motion *ompl::geometric::SLSC::selectMotion(TreeData &tree)
{
    GridCell *cell = tree.pdf.sample(rng_.uniform01());
    return cell && !cell->data.empty() ? cell->data[rng_.uniformInt(0, cell->data.size() - 1)] : nullptr;
}

bool ompl::geometric::SLSC::isPathValid(TreeData &tree, TreeData &otherTree,
                                         Motion *motion, Motion *otherMotion)
{
    std::vector<Motion *> mpath1;
    std::vector<Motion *> mpath2;

    //const base::State *root = nullptr;

    bool start = true;

    bool stateValid = false;

    while (motion->parent != nullptr || otherMotion->parent != nullptr)
    {
        stateValid = false;
        if (start && motion->parent != nullptr)
        {
            if (!(motion->valid && motion->parent->stateValid))
            {
                if (checkStartMotion(motion->parent, motion, stateValid))
                {
                    mpath1.push_back(motion);

                    motion->valid = true;
                    motion->stateValid = true;
                    motion->parent->stateValid = true;

                    motion->stepValid.clear();
                }
                else
                {
                    if (stateValid)
                    {
                        motion->parent->stateValid = true;
                        removeMotion(tree, motion);
                    }
                    else 
                        removeMotion(tree, motion->parent);

                    mpath1.clear();
                    mpath2.clear();

                    return false;
                }
            }
            else 
                mpath1.push_back(motion);

            motion = motion->parent;
        }
        else if(otherMotion->parent != nullptr) 
        {
            if (!(otherMotion->valid && otherMotion->parent->stateValid))
            {
                if (checkGoalMotion(otherMotion, otherMotion->parent, stateValid))
                {
                    mpath2.push_back(otherMotion);

                    otherMotion->valid = true;
                    otherMotion->stateValid = true;
                    otherMotion->parent->stateValid = true;

                    otherMotion->stepValid.clear();
                }
                else
                {
                    if (stateValid)
                    {
                        otherMotion->parent->stateValid = true;
                        removeMotion(otherTree, otherMotion);
                    }
                    else
                        removeMotion(otherTree, otherMotion->parent);

                    mpath1.clear();
                    mpath2.clear();

                    return false;
                }
            }
            else 
                mpath2.push_back(otherMotion);

            otherMotion = otherMotion->parent;
        }

        start = !start;
    }

    mpath1.clear();
    mpath2.clear();

    return true;
}

void ompl::geometric::SLSC::removeMotion(TreeData &tree, Motion *motion)
{
    /* remove from grid */

    GridNearest<MotionInfo, Motion*>::Coord coord(projectionEvaluator_->getDimension());
    projectionEvaluator_->computeCoordinates(motion->state, coord);
    GridNearest<MotionInfo, Motion*>::Cell *cell = tree.grid.getCell(coord);
    if (cell)
    {
        for (unsigned int i = 0; i < cell->data.size(); ++i)
        {
            if (cell->data[i] == motion)
            {
                cell->data.erase(cell->data.begin() + i);
                tree.size--;
                break;
            }
        }

        double weight;

        if (cell->data.empty())
        {
            tree.pdf.remove(cell->data.elem_);
            tree.grid.remove(cell);
            tree.grid.destroyCell(cell);

            weight = 0.0;
            for (const auto &it : tree.grid)
            {
                if (weight < tree.pdf.getWeight(it.second->data.elem_))
                    weight = tree.pdf.getWeight(it.second->data.elem_);
            }
            tree.maxWeight = weight;
        }
        else
        {
            weight = 1.0 / cell->data.size();
            tree.pdf.update(cell->data.elem_, weight);
            if (tree.maxWeight < weight)
                tree.maxWeight = weight;
        }

        weight = 1.0;
        for (const auto &it : tree.grid)
        {
            if (weight > tree.pdf.getWeight(it.second->data.elem_))
                weight = tree.pdf.getWeight(it.second->data.elem_);
        }
        tree.minWeight = weight;
    }

    /* remove self from parent list */

    if (motion->parent)
    {
        for (unsigned int i = 0; i < motion->parent->children.size(); ++i)
        {
            if (motion->parent->children[i] == motion)
            {
                motion->parent->children.erase(motion->parent->children.begin() + i);
                break;
            }
        }
    }

    /* remove children */
    for (auto &i : motion->children)
    {
        i->parent = nullptr;
        removeMotion(tree, i);
    }

    if (motion->state)
        si_->freeState(motion->state);
    delete motion;
}

void ompl::geometric::SLSC::addMotion(TreeData &tree, Motion *motion)
{
    GridNearest<MotionInfo, Motion*>::Coord coord(projectionEvaluator_->getDimension());
//    Grid<MotionInfo>::Coord coord(projectionEvaluator_->getDimension());

    projectionEvaluator_->computeCoordinates(motion->state, coord);
    GridNearest<MotionInfo, Motion*>::Cell *cell = tree.grid.getCell(coord);

    if (cell)
    {
        cell->data.push_back(motion);
        double weight = 1.0 / cell->data.size(); 
        tree.pdf.update(cell->data.elem_, weight);

        if (tree.minWeight > weight)
            tree.minWeight = weight;

        weight = 0.0;
        for (const auto &it : tree.grid)
        {
            if (weight < tree.pdf.getWeight(it.second->data.elem_))
                weight = tree.pdf.getWeight(it.second->data.elem_);
        }
        tree.maxWeight = weight;
    }
    else
    {
        cell = tree.grid.createCell(coord);
        cell->data.push_back(motion);
        tree.grid.add(cell);
        cell->data.elem_ = tree.pdf.add(cell, 1.0);

        tree.maxWeight = 1.0;
    }
    tree.size++;
}

void ompl::geometric::SLSC::getPlannerData(base::PlannerData &data) const
{
    Planner::getPlannerData(data);

    std::vector<MotionInfo> motionInfo;
    tStart_.grid.getContent(motionInfo);

    for (auto &m : motionInfo)
        for (auto &motion : m.motions_)
            if (motion->parent == nullptr)
                data.addStartVertex(base::PlannerDataVertex(motion->state, 1));
            else
                data.addEdge(base::PlannerDataVertex(motion->parent->state, 1),
                             base::PlannerDataVertex(motion->state, 1));

    motionInfo.clear();
    tGoal_.grid.getContent(motionInfo);
    for (auto &m : motionInfo)
        for (auto &motion : m.motions_)
            if (motion->parent == nullptr)
                data.addGoalVertex(base::PlannerDataVertex(motion->state, 2));
            else
                // The edges in the goal tree are reversed so that they are in the same direction as start tree
                data.addEdge(base::PlannerDataVertex(motion->state, 2),
                             base::PlannerDataVertex(motion->parent->state, 2));

    data.addEdge(data.vertexIndex(connectionPoint_.first), data.vertexIndex(connectionPoint_.second));
}

