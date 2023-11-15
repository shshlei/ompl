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

/* Author: Shi Shenglei */

#include "ompl/geometric/planners/bispace/RRTBispace.h"
#include "ompl/base/goals/GoalSampleableRegion.h"
#include "ompl/base/objectives/PathLengthOptimizationObjective.h"
#include "ompl/tools/config/SelfConfig.h"

#include <boost/format.hpp> // to delete
#include <fstream>
#include "ompl/base/spaces/RealVectorStateSpace.h"

ompl::geometric::RRTBispace::RRTBispace(const base::SpaceInformationPtr &si) : base::Planner(si, "RRTBispace")
  , dStart_(0), dGoal_(0), mc_(opt_), bh_(mc_)
{
    specs_.recognizedGoal = base::GOAL_SAMPLEABLE_REGION;
    specs_.directed = true;
    specs_.approximateSolutions = false;

    Planner::declareParam<double>("range", this, &RRTBispace::setRange, &RRTBispace::getRange, "0.:0.1:10000.");
    Planner::declareParam<double>("pen_distance", this, &RRTBispace::setPenDistance, &RRTBispace::getPenDistance, "0.:0.1:10000.");
    Planner::declareParam<bool>("lazy_path", this, &RRTBispace::setLazyPath, &RRTBispace::getLazyPath, "0,1");
    Planner::declareParam<bool>("lazy_node", this, &RRTBispace::setLazyNode, &RRTBispace::getLazyNode, "0,1");
    Planner::declareParam<bool>("use_bispace", this, &RRTBispace::setUseBispace, &RRTBispace::getUseBispace, "0,1");
//    Planner::declareParam<bool>("use_biasgrow", this, &RRTBispace::setUseBiasGrow, &RRTBispace::getUseBiasGrow, "0,1");
//    Planner::declareParam<bool>("treated_as_multi_subapce", this, &RRTBispace::setTreatedAsMultiSubapce, &RRTBispace::getTreatedAsMultiSubapce, "0,1");
}

ompl::geometric::RRTBispace::~RRTBispace()
{
    freeMemory();
}

void ompl::geometric::RRTBispace::setup()
{
    Planner::setup();
    tools::SelfConfig sc(si_, getName());
    if (treatedAsMultiSubapce_)
    {
        maxDistance_ = penDistance_ = 0.0;
        sc.configurePlannerRange(maxDistance_);
        sc.configurePenetrationDistance(penDistance_);
        auto cs = si_->getStateSpace()->as<base::CompoundStateSpace>();
        maxDistanceV_.resize(cs->getSubspaceCount());
        penDistanceV_.resize(cs->getSubspaceCount());
        double ratio = maxDistance_ / cs->getMaximumExtent();
        for (std::size_t i = 0; i < maxDistanceV_.size(); i++)
            maxDistanceV_[i] = ratio * cs->getSubspace(i)->getMaximumExtent();
        ratio = penDistance_ / cs->getMaximumExtent();
        for (std::size_t i = 0; i < penDistanceV_.size(); i++)
            penDistanceV_[i] = ratio * cs->getSubspace(i)->getMaximumExtent();
    }
    else 
    {
        sc.configurePlannerRange(maxDistance_);
        sc.configurePenetrationDistance(penDistance_);
    }

    sc.configureProjectionEvaluator(projectionEvaluator_);
    if (rewire_)
    {
        dStart_.setDimension(projectionEvaluator_->getDimension());
        dGoal_.setDimension(projectionEvaluator_->getDimension());
        dStart_.setNeighborCell(2);
        dGoal_.setNeighborCell(2);
    }

    /*
    if (useBiasGrow_)
    {
        startBiasPdf_.grid.setDimension(projectionEvaluator_->getDimension());
        goalBiasPdf_.grid.setDimension(projectionEvaluator_->getDimension());
    }
    */

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
    if (pdef_)
    {
        if (pdef_->hasOptimizationObjective())
            opt_ = pdef_->getOptimizationObjective();
        else
        {
            OMPL_INFORM("%s: No optimization objective specified. Defaulting to optimizing path length for the allowed planning time.", getName().c_str());
            opt_ = std::make_shared<base::PathLengthOptimizationObjective>(si_);
            pdef_->setOptimizationObjective(opt_);
        }
        mc_ = MotionCompare(opt_);
        bh_ = BinaryHeap<Motion *, MotionCompare>(mc_);
    }
    else
    {
        OMPL_INFORM("%s: problem definition is not set, deferring setup completion...", getName().c_str());
        setup_ = false;
    }
    
    if (lazyNode_)
        lazyPath_ = true;
    if (!lazyPath_)
    {
        lazyNode_ = false;
        rewire_ = false;
    }
    if (!rewire_)
        rewireSort_ = false;
    symmetric_ = si_->getStateSpace()->hasSymmetricDistance() && si_->getStateSpace()->hasSymmetricInterpolate();
}

void ompl::geometric::RRTBispace::freeMemory()
{
    std::vector<Motion *> motions;
    if (tStart_)
    {
        tStart_->list(motions);
        for (auto & motion : motions)
            freeMotion(motion);
        motions.clear();
    }
    if (tGoal_)
    {
        tGoal_->list(motions);
        for (auto & motion : motions)
            freeMotion(motion);
        motions.clear();
    }

    if (rewire_)
    {
        std::vector<Cell *> cells;
        dStart_.getCells(cells);
        for (auto & cell : cells)
            delete cell->data;
        cells.clear();
        dGoal_.getCells(cells);
        for (auto & cell : cells)
            delete cell->data;
        dStart_.clear();
        dGoal_.clear();
    }
}

void ompl::geometric::RRTBispace::clear()
{
    setup_ = false;

    Planner::clear();
    sampler_.reset();

    freeMemory();

    if (tStart_)
    {
        tStart_->clear();
    }

    if (tGoal_)
    {
        tGoal_->clear();
    }

    startMotions_.clear();
    goalMotions_.clear();

    pnullStartMotions_.clear();
    pnullGoalMotions_.clear();

    invalidStartMotions_.clear();
    invalidGoalMotions_.clear();
    invalidStartNum_ = invalidGoalNum_ = 0;

    connectionPoint_.clear();

    /*
    if (useBiasGrow_)
    {
        startBiasPdf_.grid.clear();
        startBiasPdf_.size = 0;
        startBiasPdf_.pdf.clear();
        goalBiasPdf_.grid.clear();
        goalBiasPdf_.size = 0;
        goalBiasPdf_.pdf.clear();
        startBiasProb_ = goalBiasProb_ = 0.0;
    }
    */
}

ompl::base::PlannerStatus ompl::geometric::RRTBispace::solve(const base::PlannerTerminationCondition &ptc)
{
    base::PlannerStatus ps= prepareSolve(ptc);
    if (ps != base::PlannerStatus::PREPARE_SUCCESS)
        return ps;

    sampler_ = si_->allocStateSampler();

    OMPL_INFORM("%s: Starting planning with %u states already in datastructure.", getName().c_str(), (tStart_->size() + tGoal_->size()));

    bool startTree = true;
    bool solved = false;
    Motion *bestStartMotion = nullptr;
    Motion *bestGoalMotion = nullptr;
    while (!ptc)
    {
        if (pis_.getSampledGoalsCount() < tGoal_->size() / 2)
        {
            const base::State *st = pis_.nextGoal();
            if (st != nullptr)
            {
                auto *motion = new Motion(si_);
                si_->copyState(motion->state, st);
                motion->valid = true;
                motion->stateValid = Valid;
                motion->root = motion->state;
                motion->cost = opt_->identityCost();
                goalMotions_.push_back(motion);
                tGoal_->add(motion);
                if (rewire_)
                {
                    Coord xcoord(projectionEvaluator_->getDimension());
                    projectionEvaluator_->computeCoordinates(motion->state, xcoord);
                    addToDisc(dGoal_, motion, xcoord);
                    if (rewireSort_)
                    {
                        motion->cell->data->cmotion = motion;
                        motion->cell->data->mmotion = motion;
                    }
                }
            }
        }

        batchGrow(startTree);

        for (auto & pair : connectionPoint_)
        {
            if (opt_->isFinite(pair.first->cost) && opt_->isFinite(pair.second->cost))
            {
                bool valid = false;
                if (rewire_)
                    valid = isPathValid(pair.first, pair.second);
                else 
                    valid = isPathValid(pair.first, pair.second, !startTree);
                if (valid)
                {
                    solved = true;
                    bestStartMotion = pair.first;
                    bestGoalMotion = pair.second;
                    break;
                }
            }
        }
        if (solved)
            break;
        if (!rewire_)
            removeInvalidMotionsDirectly();
        else if (lazyNode_)
            removeInvalidMotions();
    }

    if (solved)
    {
        ptc.terminate();
        processSolution(bestStartMotion, bestGoalMotion);
    }

    OMPL_INFORM("%s: Created %u states (%u start + %u goal).", getName().c_str(), tStart_->size() + tGoal_->size(),
                tStart_->size(), tGoal_->size());
    return solved ? base::PlannerStatus::EXACT_SOLUTION : base::PlannerStatus::TIMEOUT;
}

ompl::base::PlannerStatus ompl::geometric::RRTBispace::prepareSolve(const base::PlannerTerminationCondition &ptc)
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

    if (pis_.haveMoreStartStates())
    {
        while (const base::State *st = pis_.nextStart())
        {
            auto *motion = new Motion(si_);
            si_->copyState(motion->state, st);
            motion->valid = true;
            motion->stateValid = Valid;
            motion->root = motion->state;
            motion->cost = opt_->identityCost();
            startMotions_.push_back(motion);
            tStart_->add(motion);
            if (rewire_)
            {
                Coord xcoord(projectionEvaluator_->getDimension());
                projectionEvaluator_->computeCoordinates(motion->state, xcoord);
                addToDisc(dStart_, motion, xcoord);
                if (rewireSort_)
                {
                    motion->cell->data->cmotion = motion;
                    motion->cell->data->mmotion = motion;
                }
            }
        }
    }

    if (tStart_->size() == 0)
    {
        OMPL_ERROR("%s: There are no valid initial states!", getName().c_str());
        return base::PlannerStatus::INVALID_START;
    }

    if (const base::State *st = pis_.nextGoal(ptc))
    {
        auto *motion = new Motion(si_);
        si_->copyState(motion->state, st);
        motion->valid = true;
        motion->stateValid = Valid;
        motion->root = motion->state;
        motion->cost = opt_->identityCost();
        goalMotions_.push_back(motion);
        tGoal_->add(motion);
        if (rewire_)
        {
            Coord xcoord(projectionEvaluator_->getDimension());
            projectionEvaluator_->computeCoordinates(motion->state, xcoord);
            addToDisc(dGoal_, motion, xcoord);
            if (rewireSort_)
            {
                motion->cell->data->cmotion = motion;
                motion->cell->data->mmotion = motion;
            }
        }
    }

    if (tGoal_->size() == 0)
    {
        OMPL_ERROR("%s: Unable to sample any valid states for goal tree", getName().c_str());
        return base::PlannerStatus::INVALID_GOAL;
    }

    if (useBispace_)
    {
        double dist = -1.0;
        for (auto & sm : startMotions_)
        {
            for (auto & gm : goalMotions_)
            {
                double d = si_->distance(sm->state, gm->state);
                if (dist < d)
                    dist = d;
            }
        }
        penDistance_ = std::min(penDistance_, 0.1*dist);

        if (treatedAsMultiSubapce_)
        {
            for (std::size_t i = 0; i < si_->getStateSpace()->getSubspaceCount(); i++)
            {
                double dist = -1.0;
                for (auto & sm : startMotions_)
                {
                    for (auto & gm : goalMotions_)
                    {
                        double d = si_->distance(sm->state, gm->state, i);
                        if (dist < d)
                            dist = d;
                    }
                }
                penDistanceV_[i] = std::min(penDistanceV_[i], 0.1*dist);
            }
        }
    }
    return base::PlannerStatus::PREPARE_SUCCESS;
}

void ompl::geometric::RRTBispace::processSolution(const Motion *bestStartMotion, const Motion *bestGoalMotion)
{
    std::vector<const base::State *> spath;

    bool del = false;
    const Motion *solution = bestStartMotion;
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

    solution = bestGoalMotion;
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

    pdef_->addSolutionPath(psol);
}

bool ompl::geometric::RRTBispace::batchGrow(bool &startTree)
{
    bool nconnect = false;
    TreeGrowingInfo tgi;
    tgi.xstate = si_->allocState();
    Motion *rmotion = new Motion(si_);
    base::State *rstate = rmotion->state;

    for (unsigned int i = 0; i < 1; i++)
    {
        tgi.start = startTree;
        startTree = !startTree;
        sampler_->sampleUniform(rstate);

        bool otherSide = false;
        bool change = false;
        bool gs = growTree(tgi, rmotion, otherSide, change);
        Motion *addedMotion = tgi.xmotion;
        if (!addedMotion)
            continue;
        if (!gs)
            si_->copyState(rstate, addedMotion->state);
        tgi.start = startTree;
        bool gsc = growTree(tgi, rmotion, otherSide, change);
        Motion *startMotion = nullptr, *goalMotion = nullptr;
        if (gsc && !change)
        {
            if (tgi.xmotion->stateValid == Valid || addedMotion->stateValid == Valid)
                tgi.xmotion->stateValid = addedMotion->stateValid = Valid;
            if (!addedMotion->inConnection)
            {
                startMotion = startTree ? tgi.xmotion : addedMotion;
                goalMotion = startTree ? addedMotion : tgi.xmotion;
            }
        }
        else if (otherSide) 
        {
            addedMotion = tgi.xmotion;
            if (!addedMotion)
                continue;
            si_->copyState(rstate, addedMotion->state);
            tgi.start = !startTree;
            gsc = growTree(tgi, rmotion, otherSide, change);
            if (gsc && !change)
            {
                if (tgi.xmotion->stateValid == Valid || addedMotion->stateValid == Valid)
                    tgi.xmotion->stateValid = addedMotion->stateValid = Valid;
                if (!addedMotion->inConnection)
                {
                    startMotion = startTree ? addedMotion : tgi.xmotion;
                    goalMotion = startTree ? tgi.xmotion : addedMotion;
                }
            }
        }

        /*
        if (useBiasGrow_ && !gsc)
        {
            tgi.start = !startTree;
            gsc = biasGrow(tree, tgi, addedMotion);
            if (gsc)
            {
                startMotion = startTree ? addedMotion : tgi.xmotion; 
                goalMotion = startTree ? tgi.xmotion : addedMotion;
            }
            else 
            {
                tgi.start = startTree;
                gsc = biasGrow(otherTree, tgi, addedMotion);
                if (gsc)
                {
                    startMotion = startTree ? tgi.xmotion : addedMotion;
                    goalMotion = startTree ? addedMotion : tgi.xmotion;
                }
            }
        }
        */

        if (startMotion && goalMotion && pdef_->getGoal()->isStartGoalPairValid(startMotion->root, goalMotion->root))
        {
            nconnect = true;
            startMotion->inConnection = true;
            goalMotion->inConnection = true;
            connectionPoint_.emplace_back(startMotion, goalMotion);
        }
    }

    si_->freeState(tgi.xstate);
    freeMotion(rmotion);
    return nconnect;
}

bool ompl::geometric::RRTBispace::growTree(TreeGrowingInfo &tgi, Motion *rmotion, bool &otherSide, bool &change)
{
    if (treatedAsMultiSubapce_)
        return growTreeMultiSpace(tgi, rmotion, otherSide, change);
    return growTreeSingleSpace(tgi, rmotion, otherSide, change);
}

bool ompl::geometric::RRTBispace::growTreeSingleSpace(TreeGrowingInfo &tgi, Motion *rmotion, bool &otherSide, bool &change)
{
    otherSide = false;
    change = false;
    TreeData &tree = tgi.start ? tStart_ : tGoal_;

    bool null = false;
    Motion *nmotion = selectNMotion(tree, rmotion, null);
    tgi.xmotion = nmotion;
    if (!nmotion)
        return false;
    if (si_->equalStates(nmotion->state, rmotion->state))
        return false;
    if (null)
    {
        change = true;
        sampler_->sampleUniformNear(rmotion->state, nmotion->state, maxDistance_);
    }
    CellDiscretizationData &disc = tgi.start ? dStart_ : dGoal_;
    bool currentTree = true;
    bool addpd = false;
    double pd = 0.0;
    if (useBispace_)
        currentTree = (tgi.start == growCurrentTree(rmotion->state));
    Motion *motion = nullptr;
    bool reach = false;
    while (!reach)
    {
        nmotion = tgi.xmotion;
        base::State *dstate = rmotion->state;
        double d = null ? 0.0 : (tgi.start ? si_->distance(nmotion->state, rmotion->state) : si_->distance(rmotion->state, nmotion->state));
        if (d > maxDistance_)
        {
            if (tgi.start)
                si_->getStateSpace()->interpolate(nmotion->state, rmotion->state, std::min(maxDistance_ / d, 0.5), tgi.xstate);
            else 
                si_->getStateSpace()->interpolate(rmotion->state, nmotion->state, std::max(1.0 - maxDistance_ / d, 0.5), tgi.xstate);
            if (si_->equalStates(nmotion->state, tgi.xstate) || si_->equalStates(rmotion->state, tgi.xstate))
                break;
            dstate = tgi.xstate;
        }
        else 
            reach = true;

        if (addpd || (!currentTree && tgi.start != growCurrentTree(dstate))) 
        {
            if (addpd)
                pd += tgi.start ? si_->distance(nmotion->state, dstate) : si_->distance(dstate, nmotion->state);
            else 
            {
                Motion *last = nmotion;
                while (last && tgi.start != growCurrentTree(last->state))
                    last = last->parent;
                if (last)
                    pd += penetrationDistance(last->state, dstate, tgi.start);
                else 
                {
                    otherSide = true;
                    reach = false;
                    break;
                }
                addpd = true;
            }
            if (pd > penDistance_)
            {
                reach = false;
                break;
            }
        }

        if (!lazyNode_ && !isValid(dstate))
        {
            reach = false;
            break;
        }

        motion = new Motion(si_);
        si_->copyState(motion->state, dstate);
        if (!lazyNode_)
            motion->stateValid = Valid;

        Motion *nb = nmotion;
        if (!lazyPath_)
        {
            if (!checkInterMotion(nb, motion, tgi.start))
            {
                removeFromInvalidNeighbor(motion);
                freeMotion(motion);
                reach = false;
                break;
            }
        }

        connectToPmotion(motion, nb, tgi.start);
        motion->parent->children.push_back(motion);
        tree->add(motion);
        tgi.xmotion = motion;

        if (addpd)
        {
            otherSide = true;
            motion->middle = true;
            /*
            if (useBiasGrow_ && motion->stateValid == Valid)
            {
                if (tgi.start)
                    addPdfMotion(startBiasPdf_, motion, true);
                else 
                    addPdfMotion(goalBiasPdf_, motion, false);
            }
            */
        }

        if (!lazyPath_)
            motion->valid = true;

        if (rewire_)
        {
            Coord xcoord(projectionEvaluator_->getDimension());
            projectionEvaluator_->computeCoordinates(motion->state, xcoord);
            addToDisc(disc, motion, xcoord);
            if (rewireSort_)
            {
                if (!motion->cell->data->cmotion || opt_->isCostBetterThan(motion->cost, motion->cell->data->cmotion->cost))
                    motion->cell->data->cmotion = motion;
                if (!motion->cell->data->mmotion || opt_->isCostBetterThan(motion->cell->data->mmotion->cost, motion->cost))
                    motion->cell->data->mmotion = motion;
            }
            if (tgi.start ? !pnullStartMotions_.empty() : !pnullGoalMotions_.empty())
                rewireTree(bh_, motion, tgi.start);
        }
    }
    return reach;
}

bool ompl::geometric::RRTBispace::growTreeMultiSpace(TreeGrowingInfo &tgi, Motion *rmotion, bool &otherSide, bool &change)
{
    otherSide = false;
    change = false;
    TreeData &tree = tgi.start ? tStart_ : tGoal_;

    bool null = false;
    Motion *nmotion = selectNMotion(tree, rmotion, null);
    tgi.xmotion = nmotion;
    if (!nmotion)
        return false;
    if (si_->equalStates(nmotion->state, rmotion->state))
        return false;
    if (null)
    {
        change = true;
        sampler_->sampleUniformNear(rmotion->state, nmotion->state, maxDistance_);
    }
    CellDiscretizationData &disc = tgi.start ? dStart_ : dGoal_;
    std::vector<bool> currentTreeV;
    std::vector<bool> addpdv;
    std::vector<double> pdv;
    std::vector<bool> pstop;
    bool currentTree = true;
    if (useBispace_)
    {
        currentTreeV.resize(si_->getStateSpace()->getSubspaceCount());
        for (std::size_t i = 0; i < currentTreeV.size(); i++)
            currentTreeV[i] = (tgi.start == growCurrentTree(rmotion->state, i));
        currentTree = std::accumulate(currentTreeV.begin(), currentTreeV.end(), true, [](bool a, bool b){return a&&b;});
        if (!currentTree)
        {
            addpdv.resize(currentTreeV.size());
            std::fill(addpdv.begin(), addpdv.end(), false);
            pdv.resize(currentTreeV.size());
            std::fill(pdv.begin(), pdv.end(), 0.0);
            pstop.resize(currentTreeV.size());
            std::fill(pstop.begin(), pstop.end(), false);
        }
    }

    Motion *motion = nullptr;
    bool reach = false;
    bool reachi = true;
    bool ppstop = false;
    while (!reach && !ppstop)
    {
        nmotion = tgi.xmotion;
        base::State *dstate = rmotion->state;
        if (!currentTree)
        {
            bool stop = false;
            bool advance = false;
            reach = true;
            ppstop= true;
            dstate = tgi.xstate;
            for (std::size_t i = 0; i < currentTreeV.size(); i++)
            {
                if (pstop[i])
                {
                    si_->copyState(tgi.xstate, nmotion->state, i);
                    continue;
                }
                double d = tgi.start ? si_->distance(nmotion->state, rmotion->state, i) : si_->distance(rmotion->state, nmotion->state, i);
                if (d > maxDistanceV_[i])
                {
                    if (tgi.start)
                        si_->getStateSpace()->interpolate(nmotion->state, rmotion->state, std::min(maxDistanceV_[i] / d, 0.5), tgi.xstate, i);
                    else 
                        si_->getStateSpace()->interpolate(rmotion->state, nmotion->state, std::max(1.0 - maxDistanceV_[i] / d, 0.5), tgi.xstate, i);
                    if (si_->equalStates(nmotion->state, tgi.xstate) || si_->equalStates(rmotion->state, tgi.xstate))
                    {
                        stop = true;
                        break;
                    }
                    reach = false;
                }
                else 
                {
                    si_->copyState(tgi.xstate, rmotion->state, i);
                    pstop[i] = true;
                }
                if (addpdv[i] || (!currentTreeV[i] && tgi.start != growCurrentTree(dstate, i)))
                {
                    if (addpdv[i])
                        pdv[i] += tgi.start ? si_->distance(nmotion->state, dstate, i) : si_->distance(dstate, nmotion->state, i);
                    else 
                    {
                        Motion *last = nmotion;
                        while (last != nullptr && tgi.start != growCurrentTree(last->state, i))
                            last = last->parent;
                        if (last != nullptr)
                            pdv[i] += penetrationDistance(last->state, dstate, tgi.start, i);
                        else 
                        {
                            otherSide = true;
                            stop = true;
                            break;
                        }
                    }
                    addpdv[i] = true;
                    otherSide = true;
                    if (pdv[i] > penDistanceV_[i])
                    {
                        reachi = false;
                        pstop[i] = true;
                        si_->copyState(tgi.xstate, nmotion->state, i);
                        continue;
                    }
                }
                advance = true;
                if (d > maxDistanceV_[i])
                    ppstop = false;
            }
            if (stop || !advance)
            {
                reach = false;
                break;
            }
        }
        else 
        {
            double d = null ? 0.0 : (tgi.start ? si_->distance(nmotion->state, rmotion->state) : si_->distance(rmotion->state, nmotion->state));
            if (d > maxDistance_)
            {
                if (tgi.start)
                    si_->getStateSpace()->interpolate(nmotion->state, rmotion->state, std::min(maxDistance_ / d, 0.5), tgi.xstate);
                else 
                    si_->getStateSpace()->interpolate(rmotion->state, nmotion->state, std::max(1.0 - maxDistance_ / d, 0.5), tgi.xstate);
                if (si_->equalStates(nmotion->state, tgi.xstate) || si_->equalStates(rmotion->state, tgi.xstate))
                    break;
                dstate = tgi.xstate;
            }
            else 
                reach = true;
        }

        if (!lazyNode_ && !isValid(dstate))
        {
            reach = false;
            break;
        }

        motion = new Motion(si_);
        si_->copyState(motion->state, dstate);
        if (!lazyNode_)
            motion->stateValid = Valid;

        Motion *nb = nmotion;
        if (!lazyPath_)
        {
            if (!checkInterMotion(nb, motion, tgi.start))
            {
                removeFromInvalidNeighbor(motion);
                freeMotion(motion);
                reach = false;
                break;
            }
        }

        connectToPmotion(motion, nb, tgi.start);
        motion->parent->children.push_back(motion);
        tree->add(motion);
        tgi.xmotion = motion;

        if (otherSide)
        {
            motion->middle = true;
            /*
            if (useBiasGrow_ && motion->stateValid == Valid)
            {
                if (tgi.start)
                    addPdfMotion(startBiasPdf_, motion, true);
                else 
                    addPdfMotion(goalBiasPdf_, motion, false);
            }
            */
        }

        if (!lazyPath_)
            motion->valid = true;

        if (rewire_)
        {
            Coord xcoord(projectionEvaluator_->getDimension());
            projectionEvaluator_->computeCoordinates(motion->state, xcoord);
            addToDisc(disc, motion, xcoord);
            if (rewireSort_)
            {
                if (!motion->cell->data->cmotion || opt_->isCostBetterThan(motion->cost, motion->cell->data->cmotion->cost))
                    motion->cell->data->cmotion = motion;
                if (!motion->cell->data->mmotion || opt_->isCostBetterThan(motion->cell->data->mmotion->cost, motion->cost))
                    motion->cell->data->mmotion = motion;
            }
            if (tgi.start ? !pnullStartMotions_.empty() : !pnullGoalMotions_.empty())
                rewireTree(bh_, motion, tgi.start);
        }
    }
    return reach && reachi;
}

ompl::geometric::RRTBispace::Motion *ompl::geometric::RRTBispace::selectNMotion(const TreeData &tree, Motion *rmotion, bool &null)
{
    null = false;
    Motion *nmotion = tree->nearest(rmotion);
    if (!opt_->isFinite(nmotion->cost))
    {
        null = true;
        Motion *last = nmotion;
        if (last->stateValid == InValid)
            return nullptr;
        bool found = false;
        while (!found)
        {
            for (auto & cv : last->cell->nbh)
            {
                if (!cv.empty())
                {
//                    Cell *c = cv.size() == 1 ? cv[0] : cv[rng_.uniformInt(0, cv.size()-1)]; // todo
                    Cell *c = cv[0];
                    nmotion = selectMotionInCell(c);
                    if (nmotion)
                    {
                        found = true;
                        break;
                    }
                }
            }
            if (!found)
            {
                if (last->parent)
                    last = last->parent;
                else
                {
                    Motion *pmotion = last->pmotion;
                    if (!opt_->isFinite(pmotion->cost))
                        last = pmotion;
                    else 
                    {
                        found = true;
                        nmotion = pmotion;
                    }
                }
            }
        }
    }
    return nmotion;
}

ompl::geometric::RRTBispace::Motion *ompl::geometric::RRTBispace::selectMotionInCell(Cell *cell)
{
    Motion *nmotion = nullptr;
    if (cell->data->motions.size() > cell->data->disabled)
    {
        for (std::size_t i = cell->data->motions.size() - 1; i < cell->data->motions.size(); i--)
        {
            Motion *motion = cell->data->motions[i];
            if (opt_->isFinite(motion->cost))
            {
                nmotion = motion;
                break;
            }
        }
    }
    return nmotion;
}

void ompl::geometric::RRTBispace::rewireTree(BinaryHeap<Motion *, MotionCompare> &bh, Motion *m, bool start)
{
    if (!rewire_)
        return;
    std::vector<Motion *> &pnullMotions = start ? pnullStartMotions_ : pnullGoalMotions_;
    OrderCellsByDisabled ocbd;
    updateQueue(bh, m);
    RewireSort rsort;
    while (!bh.empty())
    {
        Motion *motion = bh.top()->data;
        motion->handle = nullptr;
        bh.pop();
        if (!opt_->isFinite(motion->cost))
            break;
        for (auto & cv : motion->cell->nbh)
        {
            std::size_t numc = cv.size(), ic = numc - 1;
            while (ic < numc)
            {
                if (rewireSort_ && ic)
                    std::sort(cv.begin(), cv.begin() + ic + 1, ocbd);
                Cell *c = cv[ic];
                ic--;
                std::size_t pos = c->data->disabled;
                if (!pos)
                {
                    if (rewireSort_)
                        break;
                    else 
                        continue;
                }
                if (rewireSort_ && motion->cell != c && c->data->mmotion) // todo
                {
                    if (opt_->isCostBetterThan(c->data->mmotion->cost, motion->cell->data->cmotion->cost))
                        continue;
                    base::Cost cost1 = opt_->combineCosts(motion->cost, opt_->motionCost(motion->state, c->data->cmotion->state));
                    base::Cost cost2 = opt_->combineCosts(c->data->mmotion->cost, opt_->combineCosts(c->data->mmotion->cost, c->data->mmotion->cost));
                    if (opt_->isCostBetterThan(cost2, cost1))
                        continue;
                }
                std::vector<Motion *> nbh;
                nbh.reserve(pos);
                if (c->data->motions.size() == c->data->disabled)
                {
                    std::sort(c->data->motions.begin(), c->data->motions.end(), rsort);               
                    nbh.insert(nbh.end(), c->data->motions.begin(), c->data->motions.end());
                }
                else 
                {
                    for (std::size_t i = 0; i < c->data->motions.size() && nbh.size() < pos; i++)
                    {
                        if (!opt_->isFinite(c->data->motions[i]->cost))
                            nbh.push_back(c->data->motions[i]);
                    }
                    std::sort(nbh.begin(), nbh.end(), rsort);
                }
                for (auto & nb : nbh)
                {
                    if (opt_->isFinite(nb->cost))
                        continue;
                    if (isInvalidNeighbor(motion, nb))
                        continue;
                    if (!nb->parent)
                        removeFromPnull(pnullMotions, nb);
                    else
                        removeFromParent(nb);
                    nb->valid = isValidNeighbor(motion, nb);
                    connectToPmotion(nb, motion, start);
                    nb->parent->children.push_back(nb);
                    enableMotionInDisc(nb);
                    updateQueue(bh, nb);
                }
            }
        }
    }

    while (!bh.empty())
    {
        bh.top()->data->handle = nullptr;
        bh.pop();
    }
    bh.clear();
}

void ompl::geometric::RRTBispace::updateQueue(BinaryHeap<Motion *, MotionCompare> &bh, Motion *motion)
{
    if (motion->handle != nullptr)
        bh.update(motion->handle);
    else
        motion->handle = bh.insert(motion);
}

bool ompl::geometric::RRTBispace::growCurrentTree(const base::State *state) const
{
    return growStartTree(state);
}

bool ompl::geometric::RRTBispace::growStartTree(const base::State *state) const
{
    double dist1 = si_->distance(startMotions_[0]->state, state);
    double dist2 = si_->distance(state, goalMotions_[0]->state);

    for (std::size_t i = 1; i < startMotions_.size(); i++)
    {
        double d = si_->distance(startMotions_[i]->state, state);
        if (d < dist1)
            dist1 = d;
    }

    for (std::size_t i = 1; i < goalMotions_.size(); i++)
    {
        double d = si_->distance(state, goalMotions_[i]->state);
        if (d < dist2)
            dist2 = d;
    }

    return dist1 <= dist2;
}

double ompl::geometric::RRTBispace::penetrationDistance(const base::State *nstate, const base::State *state, bool start) const
{
    base::State *test = si_->allocState();
    base::State *state1 = si_->allocState();
    base::State *state2 = si_->allocState();

    si_->copyState(state1, nstate);
    si_->copyState(state2, state);

    if (start)
        si_->getStateSpace()->interpolate(state1, state2, 0.5, test);
    else 
        si_->getStateSpace()->interpolate(state2, state1, 0.5, test);

    double d = start ? si_->distance(state1, state2) : si_->distance(state2, state1);;

    while (d > 1.e-6)
    {
        if (start != growStartTree(test))
        {
            si_->copyState(state2, test);
        }
        else 
        {
            si_->copyState(state1, test);
        }

        if (start)
            si_->getStateSpace()->interpolate(state1, state2, 0.5, test);
        else 
            si_->getStateSpace()->interpolate(state2, state1, 0.5, test);

        d *= 0.5;
    }

    if (start)
        d = si_->distance(test, state);
    else 
        d = si_->distance(state, test);

    si_->freeState(test);
    si_->freeState(state1);
    si_->freeState(state2);

    return d;
}

bool ompl::geometric::RRTBispace::growCurrentTree(const base::State *state, std::size_t sub) const
{
    return growStartTree(state, sub);
}

bool ompl::geometric::RRTBispace::growStartTree(const base::State *state, std::size_t sub) const
{
    double dist1 = si_->distance(startMotions_[0]->state, state, sub);
    double dist2 = si_->distance(state, goalMotions_[0]->state, sub);

    for (std::size_t i = 1; i < startMotions_.size(); i++)
    {
        double d = si_->distance(startMotions_[i]->state, state, sub);
        if (d < dist1)
            dist1 = d;
    }

    for (std::size_t i = 1; i < goalMotions_.size(); i++)
    {
        double d = si_->distance(state, goalMotions_[i]->state, sub);
        if (d < dist2)
            dist2 = d;
    }

    return dist1 <= dist2;
}

double ompl::geometric::RRTBispace::penetrationDistance(const base::State *nstate, const base::State *state, bool start, std::size_t sub) const
{
    base::State *test = si_->allocState();
    base::State *state1 = si_->allocState();
    base::State *state2 = si_->allocState();

    si_->copyState(state1, nstate, sub);
    si_->copyState(state2, state, sub);

    if (start)
        si_->getStateSpace()->interpolate(state1, state2, 0.5, test, sub);
    else 
        si_->getStateSpace()->interpolate(state2, state1, 0.5, test, sub);

    double d = start ? si_->distance(state1, state2, sub) : si_->distance(state2, state1, sub);

    while (d > 1.e-6)
    {
        if (start != growStartTree(test, sub))
        {
            si_->copyState(state2, test, sub);
        }
        else 
        {
            si_->copyState(state1, test, sub);
        }

        if (start)
            si_->getStateSpace()->interpolate(state1, state2, 0.5, test, sub);
        else 
            si_->getStateSpace()->interpolate(state2, state1, 0.5, test, sub);

        d *= 0.5;
    }

    if (start)
        d = si_->distance(test, state, sub);
    else 
        d = si_->distance(state, test, sub);

    si_->freeState(test);
    si_->freeState(state1);
    si_->freeState(state2);

    return d;
}

/*
bool ompl::geometric::RRTBispace::biasGrow(TreeData &tree, TreeGrowingInfo &tgi, Motion *&rmotion)
{
    bool gs = false;

    if (tgi.start ? goalBiasPdf_.size > 0 && rng_.uniform01() < goalBiasProb_ : startBiasPdf_.size > 0 && rng_.uniform01() < startBiasProb_)
    {
        GridCell *cell = nullptr; 
        MotionPDF &motionPdf = tgi.start ? goalBiasPdf_ : startBiasPdf_;
        rmotion = selectPdfMotion(motionPdf, cell);

        if (rmotion != nullptr)
        {
            if (!rmotion->inConnection && opt_->isFinite(rmotion->cost))
            {
                bool null;
                Motion *nmotion = selectNMotion(tree, tgi.start, rmotion, null);
                if (opt_->isFinite(nmotion->cost) && isValid(nmotion, tgi.start))
                {
                    double maxDistance = maxDistance_;
                    if (treatedAsMultiSubapce_)
                        maxDistance *= si_->getStateSpace()->getSubspaceCount();
                    if (tgi.start ? si_->distance(nmotion->state, rmotion->state) < maxDistance : si_->distance(rmotion->state, nmotion->state) < maxDistance)
                    {
                        if (checkInterMotion(nmotion, rmotion, tgi.start))
                        {
                            Motion *motion = new Motion(si_);
                            si_->copyState(motion->state, rmotion->state);
                            connectToPmotion(motion, nmotion, tgi.start);
                            motion->stateValid = Valid;
                            motion->valid = true;
                            motion->parent->children.push_back(motion);

                            if (rewire_)
                            {
                                getNeighbors(motion, tgi.start);
                                unsigned int ind = 0;
                                if (!checkIfIn(motion->nbh, nmotion, ind))
                                    insertNeighbor(motion, nmotion);
                            }

                            for (auto it = motion->nbh.begin(); it != motion->nbh.end(); it++)
                                insertNeighbor(it->first, motion);
                            setMotionValid(nmotion, motion);

                            tree->add(motion);
                            tgi.xmotion = motion;

                            gs = REACHED;
                        }
                    }
                    else
                    {
                        bool checkConnection = false;
                        bool otherSide = false;
                        if (growTree(tree, tgi, rmotion, checkConnection, otherSide, optimal) == REACHED)
                        {
                            gs = REACHED;
                        }
                    } 
                }
            }

            if(gs == REACHED)
            {
                double w = motionPdf.pdf.getWeight(cell->data.elem_);
                if (w < 1.0)
                { 
                    w /= (1.0 - w);
                    if (treatedAsMultiSubapce_)
                        w *= static_cast<double>(si_->getStateSpace()->getSubspaceCount());
                    motionPdf.pdf.update(cell->data.elem_, w);
                }
            }
            else
            {
                double w = motionPdf.pdf.getWeight(cell->data.elem_);
                if (treatedAsMultiSubapce_)
                    w /= (static_cast<double>(si_->getStateSpace()->getSubspaceCount()) + w);
                else
                    w /= (1.0 + w);
                motionPdf.pdf.update(cell->data.elem_, w);

                if (treatedAsMultiSubapce_ ? w < 0.01 : w < 0.05)
                {
                    removePdfMotion(motionPdf, rmotion);
                    rmotion->middle = false;
                }
            }

            std::vector<double> weights;
            motionPdf.pdf.getWeights(weights);
            if (!weights.empty())
            {
                if (tgi.start)
                    goalBiasProb_ = *std::max_element(weights.begin(), weights.end());
                else 
                    startBiasProb_ = *std::max_element(weights.begin(), weights.end());
            }
        }
    }

    return gs;
}
*/

bool ompl::geometric::RRTBispace::isPathValid(Motion *motion, Motion *otherMotion, bool start)
{
    bool valid = true;
    if (start)
    {
        if (!isPathValid(motion, true))
            valid = false;
        else if (!isPathValid(otherMotion, false))
            valid = false;
    }
    else 
    {
        if (!isPathValid(otherMotion, false))
            valid = false;
        else if (!isPathValid(motion, true))
            valid = false;
    }
    return valid;
}

bool ompl::geometric::RRTBispace::isPathValid(Motion *motion, Motion *otherMotion)
{
    bool valid = true;
    if (!isPathValid(motion, true))
        valid = false;
    if (!isPathValid(otherMotion, false))
        valid = false;
    return valid;
}

bool ompl::geometric::RRTBispace::isPathValid(Motion *motion, bool start)
{
    if (!lazyPath_)
        return true;
    if (lazyNode_ && !isStateValid(motion, start))
        return false;
    bool tvalid = true;
    std::vector<Motion *> mpath;
    while (motion->parent)
    {
        mpath.push_back(motion);
        motion = motion->parent;
    }
    std::vector<Motion *> &pnullMotions= start ? pnullStartMotions_: pnullGoalMotions_;
    std::vector<Motion *> &invalidMotions = start ? invalidStartMotions_ : invalidGoalMotions_;
    for (std::size_t i = mpath.size() - 1; i < mpath.size(); i--)
    {
        motion = mpath[i];
        Motion *pmotion = motion->parent;
        if (!checkMotion(pmotion, motion, start))
        {
            removeFromParent(motion);
            motion->parent = nullptr;
            tvalid = false;
            if (rewire_)
            {
                Motion *ppmotion = nullptr;
                std::size_t epos = invalidMotions.size(), spos = epos;
                tvalid = backPathRewireMotion(motion, start, ppmotion);
                epos = invalidMotions.size();
                for (std::size_t i = spos; i < epos; i++)
                {
                    Motion *temp = invalidMotions[i];
                    for (auto & child : temp->children)
                        child->pmotion = pmotion;
                }
                if (tvalid)
                {
                    tvalid = isPathValidInter(ppmotion, start);
                    connectToPmotion(motion, ppmotion, start);
                    motion->parent->children.push_back(motion); 
                }
                else 
                {
                    pnullMotions.push_back(motion);
                    motion->valid = false;
                    motion->pmotion = pmotion;
                }
            }
            else 
                pnullMotions.push_back(motion);
            if (tvalid)
                enableMotionInDisc(motion);
            else 
                break;
        }
    }
    return tvalid;
}

bool ompl::geometric::RRTBispace::isStateValid(Motion *motion, bool start)  // todo
{
    if (!lazyNode_)
        return true;
    bool tvalid = true;
    std::vector<Motion *> mpath;
    while (motion->parent)
    {
        mpath.push_back(motion);
        motion = motion->parent;
    }
    std::vector<Motion *> &pnullMotions= start ? pnullStartMotions_: pnullGoalMotions_;
    std::vector<Motion *> &invalidMotions = start ? invalidStartMotions_ : invalidGoalMotions_;
    for (std::size_t i = mpath.size() - 1; i < mpath.size(); i--)
    {
        motion = mpath[i];
        Motion *pmotion = motion->parent;
        if (!isValid(motion, start))
        {
            tvalid = false;
            if (rewire_)
            {
                std::size_t epos = invalidMotions.size(), spos = epos - 1;
                Motion *last = nullptr;
                i--;
                while (i < mpath.size())
                {
                    if (isValid(mpath[i], start))
                    {
                        last = mpath[i];
                        break;
                    }
                    i--;
                }
                if (last != nullptr)
                {
                    removeFromParent(last);
                    last->parent = nullptr;
                    Motion *plast = nullptr;
                    tvalid = backPathRewireMotion(last, start, plast);
                    epos = invalidMotions.size();
                    if (tvalid)
                    {
                        tvalid = isPathValidInter(plast, start);
                        connectToPmotion(last, plast, start);
                        last->parent->children.push_back(last);
                        if (tvalid)
                            enableMotionInDisc(last);
                    }
                    else 
                    {
                        last->valid = false;
                        last->pmotion = pmotion;
                        pnullMotions.push_back(last);
                    }
                }
                else 
                    epos = invalidMotions.size();
                for (std::size_t i = spos; i < epos; i++)
                {
                    Motion *temp = invalidMotions[i];
                    for (auto & child : temp->children)
                        child->pmotion = pmotion;
                }
            }
            if (!tvalid)
                break;
        }
    }
    return tvalid;
}

bool ompl::geometric::RRTBispace::isPathValidInter(Motion *motion, bool start)
{
    if (lazyNode_ && !isStateValid(motion, start))
        return false;
    bool tvalid = true;
    std::vector<Motion *> mpath;
    while (motion->parent)
    {
        mpath.push_back(motion);
        motion = motion->parent;
    }
    std::vector<Motion *> &pnullMotions= start ? pnullStartMotions_: pnullGoalMotions_;
    std::vector<Motion *> &invalidMotions = start ? invalidStartMotions_ : invalidGoalMotions_;
    for (std::size_t i = mpath.size() - 1; i < mpath.size(); i--)
    {
        motion = mpath[i];
        Motion *pmotion = motion->parent;
        if (!checkMotion(pmotion, motion, start))
        {
            removeFromParent(motion);
            motion->parent = nullptr;
            tvalid = false;
            if (rewire_)
            {
                Motion *ppmotion = nullptr;
                std::size_t epos = invalidMotions.size(), spos = epos;
                tvalid = backPathRewireMotion(motion, start, ppmotion);
                epos = invalidMotions.size();
                for (std::size_t i = spos; i < epos; i++)
                {
                    Motion *temp = invalidMotions[i];
                    for (auto & child : temp->children)
                        child->pmotion = pmotion;
                }
                if (tvalid)
                {
                    tvalid = isPathValidInter(ppmotion, start);
                    connectToPmotion(motion, ppmotion, start);
                    motion->parent->children.push_back(motion); 
                }
                else 
                {
                    pnullMotions.push_back(motion);
                    motion->valid = false;
                    motion->pmotion = pmotion;
                }
            }
            else 
                pnullMotions.push_back(motion);
            if (tvalid)
                enableMotionInDisc(motion);
            else 
                break;
        }
    }
    return tvalid;
}

/*
bool ompl::geometric::RRTBispace::isPathValidInter(Motion *motion, bool start)
{
    if (lazyNode_ && !isStateValid(motion, start))
        return false;
    bool tvalid = true;
    std::vector<Motion *> &pnullMotions= start ? pnullStartMotions_: pnullGoalMotions_;
    std::vector<Motion *> &invalidMotions = start ? invalidStartMotions_ : invalidGoalMotions_;
    while (motion->parent)
    {
        Motion *pmotion = motion->parent;
        if (!checkMotion(pmotion, motion, start))
        {
            removeFromParent(motion);
            motion->parent = nullptr;
            tvalid = false;
            if (rewire_)
            {
                Motion *ppmotion = nullptr;
                std::size_t epos = invalidMotions.size(), spos = epos;
                tvalid = backPathRewireMotion(motion, start, ppmotion);
                epos = invalidMotions.size();
                for (std::size_t i = spos; i < epos; i++)
                {
                    Motion *temp = invalidMotions[i];
                    for (auto & child : temp->children)
                        child->pmotion = pmotion;
                }
                if (tvalid)
                {
                    tvalid = isPathValidInter(ppmotion, start);
                    connectToPmotion(motion, ppmotion, start);
                    motion->parent->children.push_back(motion);
                    if (tvalid)
                        enableMotionInDisc(motion);
                }
                else 
                {
                    pnullMotions.push_back(motion);
                    motion->valid = true;
                    motion->pmotion = pmotion;
                    motion->cell->data->root++;
                    for (auto & child : motion->children)
                    {
                        child->valid = false;
                        child->parent = nullptr;
                        child->pmotion = motion;
                        pnullMotions.push_back(child);
                        child->cell->data->root++;
                    }
                    motion->children.clear();
                }
            }
            else
                pnullMotions.push_back(motion);
            break;
        }
        motion = motion->parent;
    }
    return tvalid;
}
*/

bool ompl::geometric::RRTBispace::backPathRewireMotion(Motion *motion, bool start, Motion *&pmotion)
{
    bool valid = false;
    if (!rewire_)
        return valid;
    OrderCellsByCost ocbc(opt_);
    CellDiscretizationData &disc = start ? dStart_ : dGoal_;
    std::size_t cnbh = motion->cell->nbh.size();
    for (std::size_t cindex = 0; cindex < cnbh; cindex++)
    {
        std::size_t numc = motion->cell->nbh[cindex].size(), ic = numc - 1;
        while (ic < numc)
        {
            if (rewireSort_ && ic)
                std::sort(motion->cell->nbh[cindex].begin(), motion->cell->nbh[cindex].begin() + ic + 1, ocbc);
            Cell *c = motion->cell->nbh[cindex][ic];
            ic--;
            std::size_t num = c->data->motions.size() - c->data->disabled;
            if (!num)
            {
                if (rewireSort_)
                    break;
                else 
                    continue;
            }
            if (rewireSort_ && motion->cell != c && motion->cell->data->mmotion) // todo
            {
                if (opt_->isCostBetterThan(motion->cell->data->mmotion->cost, c->data->cmotion->cost))
                    continue;
                base::Cost cost1 = opt_->combineCosts(c->data->cmotion->cost, opt_->motionCost(c->data->cmotion->state, motion->state));
                base::Cost cost2 = opt_->combineCosts(motion->cell->data->mmotion->cost, opt_->combineCosts(motion->cell->data->mmotion->cost, motion->cell->data->mmotion->cost));
                if (opt_->isCostBetterThan(cost2, cost1))
                    continue;
            }
            std::vector<Motion *> nbh;
            nbh.reserve(num);
            if (!c->data->disabled)
                nbh.insert(nbh.end(), c->data->motions.begin(), c->data->motions.end());
            else 
            {
                for (std::size_t i = 0; i < c->data->motions.size() && nbh.size() < num; i++)
                {
                    if (opt_->isFinite(c->data->motions[i]->cost))
                        nbh.push_back(c->data->motions[i]);
                }
            }
            for (auto & nb : nbh)
            {
                if (!opt_->isFinite(nb->cost))
                    continue;
                if (!isValid(nb, start))
                    continue;
                if (isInvalidNeighbor(motion, nb))
                    continue;
                if (!isValidNeighbor(motion, nb))
                {
                    if (checkInterMotion(nb, motion, start))
                    {
                        insertNeighbor(nb, motion);
                        disc.updateNbh(nb->cell, motion->cell);
                    }
                    else
                    {
                        insertInvalidNeighbor(nb, motion);
                        continue;
                    }
                }
                pmotion = nb;
                motion->valid = true;
                valid = true;
                break;
            }
            if (valid)
                break;
        }
        if (valid)
            break;
    }

    return valid;
}

void ompl::geometric::RRTBispace::removeInvalidMotionsDirectly()
{
    for (auto & pair : connectionPoint_)
    {
        pair.first->inConnection = false;
        pair.second->inConnection = false;
    }
    connectionPoint_.clear();
    removeInvalidMotionsDirectlyTree();
}

void ompl::geometric::RRTBispace::removeInvalidMotionsDirectlyTree()
{
    if (!invalidStartMotions_.empty() || !pnullStartMotions_.empty())
    {
        tStart_->clear();
        for (auto & rootMotion : startMotions_)
            addToTree(tStart_, rootMotion);
        for (auto & pnull : invalidStartMotions_)
            removeFromTree(startBiasPdf_, pnull);
        for (auto & pnull : pnullStartMotions_)
            removeFromTree(startBiasPdf_, pnull);
        invalidStartMotions_.clear();
        pnullStartMotions_.clear();
    }

    if (!invalidGoalMotions_.empty() || !pnullGoalMotions_.empty())
    {
        tGoal_->clear();
        for (auto & rootMotion : goalMotions_)
            addToTree(tGoal_, rootMotion);
        for (auto & pnull : invalidGoalMotions_)
            removeFromTree(goalBiasPdf_, pnull);
        for (auto & pnull : pnullGoalMotions_)
            removeFromTree(goalBiasPdf_, pnull);
        invalidGoalMotions_.clear();
        pnullGoalMotions_.clear();
    }
}

void ompl::geometric::RRTBispace::addToTree(TreeData &tree, Motion *motion)
{
    tree->add(motion);
    for (auto & child : motion->children)
        addToTree(tree, child);
}

void ompl::geometric::RRTBispace::removeFromTree(MotionPDF &pdf, Motion *motion)
{
//    if (useBiasGrow_ && motion->middle && motion->stateValid == Valid)
//        removePdfMotion(pdf, motion);
    for (auto & child : motion->children)
    {
        child->parent = nullptr;
        removeFromTree(pdf, child);
    }
    freeMotion(motion);
}

void ompl::geometric::RRTBispace::removeInvalidMotions()
{
    std::size_t num = connectionPoint_.size(), i = num - 1;
    while (i < num)
    {
        auto pair = connectionPoint_[i];
        if (pair.first->stateValid == InValid || pair.second->stateValid == InValid)
        {
            pair.first->inConnection = false;
            pair.second->inConnection = false;
            std::iter_swap(connectionPoint_.begin() + i, connectionPoint_.end() - 1);
            connectionPoint_.pop_back();
        }
        i--;
    }
    removeInvalidMotionsTree(maxInvalidNodeRatio_);
}

void ompl::geometric::RRTBispace::removeInvalidMotionsTree(double ratio)
{
    for (std::size_t i = invalidStartNum_; i < invalidStartMotions_.size(); i++)
    {
        Motion *pnull = invalidStartMotions_[i];
        for (auto & child : pnull->children)
        {
            child->parent = nullptr;
            pnullStartMotions_.push_back(child);
        }
        pnull->children.clear();
    }
    invalidStartNum_ = invalidStartMotions_.size();
    if ((double)invalidStartMotions_.size() / (double) tStart_->size() > ratio)
    {
        tStart_->clear();
        for (auto & rootMotion : startMotions_)
            addToTree(tStart_, rootMotion);
        for (auto & pnull : invalidStartMotions_)
            freeMotion(pnull);
        invalidStartMotions_.clear();
        invalidStartNum_ = 0;
        for (auto & pnull : pnullStartMotions_)
            addToTree(tStart_, pnull);
    }

    for (std::size_t i = invalidGoalNum_; i < invalidGoalMotions_.size(); i++)
    {
        Motion *pnull = invalidGoalMotions_[i];
        for (auto & child : pnull->children)
        {
            child->parent = nullptr;
            pnullGoalMotions_.push_back(child);
        }
        pnull->children.clear();
    }
    invalidGoalNum_ = invalidGoalMotions_.size();
    if ((double)invalidGoalMotions_.size() / (double) tGoal_->size() > ratio)
    {
        tGoal_->clear();
        for (auto & rootMotion : goalMotions_)
            addToTree(tGoal_, rootMotion);
        for (auto & pnull : invalidGoalMotions_)
            freeMotion(pnull);
        invalidGoalMotions_.clear();
        invalidGoalNum_ = 0;
        for (auto & pnull : pnullGoalMotions_)
            addToTree(tGoal_, pnull);
    }
}

void ompl::geometric::RRTBispace::addPdfMotion(MotionPDF &pdf, Motion *motion, bool start)
{
    Grid<MotionInfo>::Coord coord(projectionEvaluator_->getDimension());
    projectionEvaluator_->computeCoordinates(motion->state, coord);

    Grid<MotionInfo>::Cell *cell = pdf.grid.getCell(coord);
    if (cell)
        cell->data.push_back(motion);
    else
    {
        cell = pdf.grid.createCell(coord);
        cell->data.push_back(motion);
        pdf.grid.add(cell);
        cell->data.elem_ = pdf.pdf.add(cell, 0.1);

        if (start)
            startBiasProb_ = 0.1;
        else 
            goalBiasProb_ = 0.1;
    }
    pdf.size++;

    /*
    if (savePd_)
    {
        base::PlannerData data(si_);
        Planner::getPlannerData(data);
        std::vector<Motion *> motions;
        tStart_->list(motions);
        for (auto & motion : motions)
        {
            if (motion->parent == nullptr)
                data.addStartVertex(base::PlannerDataVertex(motion->state, 1));
            else
            {
                data.addEdge(base::PlannerDataVertex(motion->parent->state, 1), base::PlannerDataVertex(motion->state, 1));
            }
        }
        pds_.store(data, boost::str(boost::format("pdstart_%i.txt") % pdi_).c_str());
    }

    if (savePd_)
    {
        base::PlannerData data(si_);
        Planner::getPlannerData(data);
        std::vector<Motion *> motions;
        tGoal_->list(motions);
        for (auto & motion : motions)
        {
            if (motion->parent == nullptr)
                data.addGoalVertex(base::PlannerDataVertex(motion->state, 2));
            else
            {
                // The edges in the goal tree are reversed to be consistent with start tree
                data.addEdge(base::PlannerDataVertex(motion->state, 2), base::PlannerDataVertex(motion->parent->state, 2));
            }
        }
        pds_.store(data, boost::str(boost::format("pdgoal_%i.txt") % pdi_).c_str());
    }

    if (savePd_)
    {
        MotionPDF &motionPdf = goalBiasPdf_;
        Grid<MotionInfo>::CellArray cells;
        motionPdf.grid.getCells(cells);
        std::ofstream ofs(boost::str(boost::format("pdstartbias_%i.txt") % pdi_).c_str());
        for (auto & cell : cells)
        {
            double w = motionPdf.pdf.getWeight(cell->data.elem_);
            for (auto & motion : cell->data)
            {
                auto state = motion->state->as<ompl::base::RealVectorStateSpace::StateType>();
                ofs << state->values[0] << " " << state->values[1] << " " << w << std::endl;
            }
        }
        ofs.close();
    }

    if (savePd_)
    {
        MotionPDF &motionPdf = startBiasPdf_;
        Grid<MotionInfo>::CellArray cells;
        motionPdf.grid.getCells(cells);
        std::ofstream ofs(boost::str(boost::format("pdgoalbias_%i.txt") % pdi_).c_str());
        for (auto & cell : cells)
        {
            double w = motionPdf.pdf.getWeight(cell->data.elem_);
            for (auto & motion : cell->data)
            {
                auto state = motion->state->as<ompl::base::RealVectorStateSpace::StateType>();
                ofs << state->values[0] << " " << state->values[1] << " " << w << std::endl;
            }
        }
        ofs.close();
    }

    pdi_++;
    */
}

ompl::geometric::RRTBispace::Motion *ompl::geometric::RRTBispace::selectPdfMotion(MotionPDF &pdf, GridCell *&cell)
{
    cell = pdf.pdf.sample(rng_.uniform01());
    if (cell && !cell->data.empty())
    {
        double w = pdf.pdf.getWeight(cell->data.elem_);
        if (treatedAsMultiSubapce_)
            w /= (static_cast<double>(si_->getStateSpace()->getSubspaceCount()) + w);
        else
            w /= (1.0 + w);
        pdf.pdf.update(cell->data.elem_, w);
        return cell->data[rng_.uniformInt(0, cell->data.size() - 1)];
    }
    else 
        return nullptr;
}

void ompl::geometric::RRTBispace::removePdfMotion(MotionPDF &pdf, Motion *motion)
{
    Grid<MotionInfo>::Coord coord(projectionEvaluator_->getDimension());
    projectionEvaluator_->computeCoordinates(motion->state, coord);

    Grid<MotionInfo>::Cell *cell = pdf.grid.getCell(coord);
    if (cell)
    {
        for (std::size_t i = 0; i < cell->data.size(); ++i)
        {
            if (cell->data[i] == motion)
            {
                std::iter_swap(cell->data.begin() + i, cell->data.end() - 1);
                cell->data.pop_back();
                pdf.size--;
                break;
            }
        }

        if (cell->data.empty())
        {
            pdf.pdf.remove(cell->data.elem_);
            pdf.grid.remove(cell);
            pdf.grid.destroyCell(cell);
        }
    }
}

void ompl::geometric::RRTBispace::enableMotionInDisc(Motion *motion)
{
    motion->cell->data->disabled--;
    if (rewireSort_)
    {
        if (!motion->cell->data->cmotion || opt_->isCostBetterThan(motion->cost, motion->cell->data->cmotion->cost))
            motion->cell->data->cmotion = motion;
        if (!motion->cell->data->mmotion || opt_->isCostBetterThan(motion->cell->data->mmotion->cost, motion->cost))
            motion->cell->data->mmotion = motion;
    }
    for (auto & child : motion->children)
        enableMotionInDisc(child);
}

void ompl::geometric::RRTBispace::connectToPmotion(Motion *motion, Motion *pmotion, bool start) const
{
    motion->parent = pmotion;
    motion->root = pmotion->root;
    motion->incCost = start ? opt_->motionCost(pmotion->state, motion->state) : opt_->motionCost(motion->state, pmotion->state);
    motion->cost = opt_->combineCosts(pmotion->cost, motion->incCost);
    if (opt_->isFinite(motion->cost))
        updateChildCosts(motion);
}

void ompl::geometric::RRTBispace::updateChildCosts(Motion *motion) const
{
    for (auto & child : motion->children)
    {
        child->cost = opt_->combineCosts(motion->cost, child->incCost);
        child->root = motion->root;
        updateChildCosts(child);
    }
}

void ompl::geometric::RRTBispace::insertNeighbor(Motion *pmotion, Motion *motion)
{
    pmotion->nbh.insert(motion);
    motion->nbh.insert(pmotion);
}

void ompl::geometric::RRTBispace::insertInvalidNeighbor(Motion *pmotion, Motion *motion)
{
    pmotion->invalidnbh.insert(motion);
    motion->invalidnbh.insert(pmotion);
}

void ompl::geometric::RRTBispace::removeFromNeighbor(Motion *motion)
{
    for (auto & pmotion : motion->nbh)
        pmotion->nbh.erase(motion);
}

void ompl::geometric::RRTBispace::removeFromInvalidNeighbor(Motion *motion)
{
    for (auto & pmotion : motion->invalidnbh)
        pmotion->invalidnbh.erase(motion);
}

bool ompl::geometric::RRTBispace::isValidNeighbor(Motion *motion, Motion *nb) const
{
    return motion->nbh.find(nb) != motion->nbh.end();
}

bool ompl::geometric::RRTBispace::isInvalidNeighbor(Motion *motion, Motion *nb) const
{
    return motion->invalidnbh.find(nb) != motion->invalidnbh.end();
}

void ompl::geometric::RRTBispace::addToDisc(CellDiscretizationData &disc, Motion *motion, const Coord& coord)
{
    Cell *cell = disc.getCell(coord);
    if (cell)
        cell->data->motions.push_back(motion);
    else
    {
        cell = static_cast<Cell *>(disc.createCell(coord));
        cell->data = new CellData();
        cell->data->motions.push_back(motion);
        disc.add(cell);
    }
    motion->cell = cell;
}

void ompl::geometric::RRTBispace::removeFromDisc(CellDiscretizationData &disc, Motion *motion)
{
    Cell *cell = motion->cell;
    if (!cell)
        return;
    if (removeFromVector(cell->data->motions, motion))
        motion->cell = nullptr;
    if (cell->data->motions.empty())
    {
        disc.remove(cell);
        delete cell->data;
        disc.destroyCell(cell);
    }
}

void ompl::geometric::RRTBispace::removeFromParent(Motion *motion)
{
    if (!motion->parent)
        return;
    removeFromVector(motion->parent->children, motion);
}

void ompl::geometric::RRTBispace::removeFromPnull(std::vector<Motion *> &pnullMotions, Motion *motion)
{
    if (removeFromVector(pnullMotions, motion))
        motion->pmotion = nullptr;
}

bool ompl::geometric::RRTBispace::removeFromVector(std::vector<Motion *> &motions, Motion *motion)
{
    bool found = false;
    for (std::size_t i = motions.size() - 1; i < motions.size(); i--)
    {
        if (motions[i] == motion)
        {
            found = true;
            std::iter_swap(motions.begin() + i, motions.end() - 1);
            motions.pop_back();
            break;
        }
    }
    return found;
}

void ompl::geometric::RRTBispace::setMotionInfinityCost(Motion *motion) const
{
    motion->cost = opt_->infiniteCost();
    for (auto & child : motion->children)
        setMotionInfinityCost(child);
}

void ompl::geometric::RRTBispace::setMotionInfinityCostWithDisable(Motion *motion, std::unordered_set<Cell *> &cells) const
{
    motion->cost = opt_->infiniteCost();
    motion->cell->data->disabled++;
    if (rewireSort_)
        cells.insert(motion->cell);
    for (auto & child : motion->children)
        setMotionInfinityCostWithDisable(child, cells);
}

// check 
bool ompl::geometric::RRTBispace::isValid(const base::State *state)
{
    bool valid = si_->isValid(state);
    return valid;
}

bool ompl::geometric::RRTBispace::isValid(Motion *motion, bool start)
{
    if (motion->stateValid == UnCkecked)
    {
        if (isValid(motion->state))
        {
            motion->stateValid = Valid;
            /*
            if (useBiasGrow_ && motion->middle)
            {
                if (start)
                    addPdfMotion(startBiasPdf_, motion, true);
                else 
                    addPdfMotion(goalBiasPdf_, motion, false);
            }
            */
            return true;
        }
        motion->stateValid = InValid;
        CellDiscretizationData &disc = start ? dStart_ : dGoal_;
        removeFromParent(motion);
        motion->parent = nullptr;
        if (start)
            invalidStartMotions_.push_back(motion);
        else 
            invalidGoalMotions_.push_back(motion);
        if (opt_->isFinite(motion->cost))
        {
            std::unordered_set<Cell *> cells;
            if (rewire_)
            {
                setMotionInfinityCostWithDisable(motion, cells);
                if (motion->cell->data->motions.size() == 1)
                    cells.erase(motion->cell);
                else 
                    motion->cell->data->disabled--;
                removeFromDisc(disc, motion);
            }
            else 
                setMotionInfinityCost(motion);
            if (rewireSort_)
            {
                for (auto & cell : cells)
                {
                    if (cell->data->motions.size() == cell->data->disabled)
                    {
                        cell->data->cmotion = nullptr;
                        cell->data->mmotion = nullptr;
                    }
                    else
                    {
                        bool e1 = !cell->data->cmotion || !opt_->isFinite(cell->data->cmotion->cost);
                        bool e2 = !cell->data->mmotion || !opt_->isFinite(cell->data->mmotion->cost);
                        if (e1 || e2)
                        {
                            EnableSort esort(opt_);
                            std::sort(cell->data->motions.begin(), cell->data->motions.end(), esort);
                            if (e1)
                                cell->data->cmotion = *std::min_element(cell->data->motions.begin(), cell->data->motions.end() - cell->data->disabled, mc_);
                            if (e2)
                                cell->data->mmotion = *std::max_element(cell->data->motions.begin(), cell->data->motions.end() - cell->data->disabled, mc_);
                        }
                    }
                }
            }
        }
        else if (rewire_)
        {
            motion->cell->data->disabled--;
            removeFromDisc(disc, motion);
        }
    }
    return motion->stateValid == Valid;
}

bool ompl::geometric::RRTBispace::checkMotion(Motion *pmotion, Motion *motion, bool start)
{
    if (!motion->valid)
    {
        if (!checkInterMotion(pmotion, motion, start))
        {
            if (opt_->isFinite(motion->cost))
            {
                if (rewire_)
                {
                    std::unordered_set<Cell *> cells;
                    setMotionInfinityCostWithDisable(motion, cells);
                    if (rewireSort_)
                    {
                        for (auto & cell : cells)
                        {
                            if (cell->data->motions.size() == cell->data->disabled)
                            {
                                cell->data->cmotion = nullptr;
                                cell->data->mmotion = nullptr;
                            }
                            else
                            {
                                bool e1 = !cell->data->cmotion || !opt_->isFinite(cell->data->cmotion->cost);
                                bool e2 = !cell->data->mmotion || !opt_->isFinite(cell->data->mmotion->cost);
                                if (e1 || e2)
                                {
                                    EnableSort esort(opt_);
                                    std::sort(cell->data->motions.begin(), cell->data->motions.end(), esort);
                                    if (e1)
                                        cell->data->cmotion = *std::min_element(cell->data->motions.begin(), cell->data->motions.end() - cell->data->disabled, mc_);
                                    if (e2)
                                        cell->data->mmotion = *std::max_element(cell->data->motions.begin(), cell->data->motions.end() - cell->data->disabled, mc_);
                                }
                            }
                        }
                    }
                }
                else 
                    setMotionInfinityCost(motion);
            }
            if (rewire_)
                insertInvalidNeighbor(pmotion, motion);
        }
        else 
        {
            motion->valid = true;
            if (rewire_)
            {
                insertNeighbor(pmotion, motion);
                CellDiscretizationData &disc = start ? dStart_ : dGoal_;
                disc.updateNbh(pmotion->cell, motion->cell);
            }
        }
    }
    return motion->valid;
}

bool ompl::geometric::RRTBispace::checkInterMotion(Motion *pmotion, Motion *motion, bool start)
{
    bool valid = false;
    if (start || symmetric_)
        valid = checkInterMotion1(pmotion, motion, start);
    else 
        valid = checkInterMotion2(motion, pmotion, start);
    return valid;
}

bool ompl::geometric::RRTBispace::checkInterMotion1(Motion *smotion, Motion *gmotion, bool start)
{
    bool valid = true;
    /*assume smotion, gmotion are valid*/
    const base::State *s1 = smotion->state, *s2 = gmotion->state;
    double dist = si_->distance(s1, s2);
    auto space = si_->getStateSpace();
    int nd = (int)ceil(dist / space->getLongestValidSegmentLength());
    if (nd >= 2)
    {
        base::State *state= si_->allocState();
        double ratio = 1.0 / (double)nd, delta = ratio;
        for (int i = 1; i < nd; i++)
        {
            si_->getStateSpace()->interpolate(s1, s2, ratio, state);
            if (!si_->isValid(state))
            {
                valid = false;
                break;
            }
            ratio += delta;
        }
        if (!valid && addIntermediateState_ && ratio > 0.5 && ratio * dist > 0.25 * maxDistance_)
        {
            ratio -= delta;
            Motion *last = new Motion(si_);
            si_->getStateSpace()->interpolate(s1, s2, ratio, last->state);
            addIntermediateMotion(smotion, gmotion, start, last);
        }
        si_->freeState(state);
    }
    return valid;
}

bool ompl::geometric::RRTBispace::checkInterMotion2(Motion *smotion, Motion *gmotion, bool start)
{
    bool valid = true;
    /*assume smotion, gmotion are valid*/
    const base::State *s1 = smotion->state, *s2 = gmotion->state;
    double dist = si_->distance(s1, s2);
    auto space = si_->getStateSpace();
    int nd = (int)ceil(dist / space->getLongestValidSegmentLength());
    if (nd >= 2)
    {
        base::State *state= si_->allocState();
        double delta = 1.0 / (double)nd, ratio = 1.0 - delta;
        for (int i = nd - 1; i > 0; i--)
        {
            si_->getStateSpace()->interpolate(s1, s2, ratio, state);
            if (!si_->isValid(state))
            {
                valid = false;
                break;
            }
            ratio -= delta;
        }
        if (!valid && addIntermediateState_ && ratio < 0.5 && (1.0 - ratio) * dist > 0.25 * maxDistance_)
        {
            ratio += delta;
            Motion *last = new Motion(si_);
            si_->getStateSpace()->interpolate(s1, s2, ratio, last->state);
            addIntermediateMotion(gmotion, smotion, start, last);
        }
        si_->freeState(state);
    }
    return valid;
}

void ompl::geometric::RRTBispace::addIntermediateMotion(Motion *pmotion, Motion *motion, bool start, Motion *last)
{
    TreeData &tree = start ? tStart_ : tGoal_;
    last->valid = true;
    last->stateValid = Valid;
    connectToPmotion(last, pmotion, start);
    last->parent->children.push_back(last);
    tree->add(last);
    if (rewire_)
    {
        CellDiscretizationData &disc = start ? dStart_ : dGoal_;
        Coord xcoord(projectionEvaluator_->getDimension());
        projectionEvaluator_->computeCoordinates(last->state, xcoord);
        addToDisc(disc, last, xcoord);
        if (!opt_->isFinite(last->cost))
            last->cell->data->disabled++;
        insertNeighbor(last, pmotion);
        disc.updateNbh(last->cell, pmotion->cell);
        if (motion->stateValid == Valid)
            insertInvalidNeighbor(last, motion);
        if (rewireSort_)
        {
            if (!last->cell->data->cmotion || opt_->isCostBetterThan(last->cost, last->cell->data->cmotion->cost))
                last->cell->data->cmotion = last;
            if (!last->cell->data->mmotion || opt_->isCostBetterThan(last->cell->data->mmotion->cost, last->cost))
                last->cell->data->mmotion = last;
        }
    }
}

void ompl::geometric::RRTBispace::getPlannerData(base::PlannerData &data) const
{
    Planner::getPlannerData(data);

    std::vector<Motion *> motions;
    if (tStart_)
        tStart_->list(motions);
    for (auto & motion : motions)
    {
        if (motion->parent == nullptr)
            data.addStartVertex(base::PlannerDataVertex(motion->state, 1));
        else
        {
            data.addEdge(base::PlannerDataVertex(motion->parent->state, 1), base::PlannerDataVertex(motion->state, 1));
        }
    }

    motions.clear();
    if (tGoal_)
        tGoal_->list(motions);
    for (auto & motion : motions)
    {
        if (motion->parent == nullptr)
            data.addGoalVertex(base::PlannerDataVertex(motion->state, 2));
        else
        {
            // The edges in the goal tree are reversed to be consistent with start tree
            data.addEdge(base::PlannerDataVertex(motion->state, 2), base::PlannerDataVertex(motion->parent->state, 2));
        }
    }
}

void ompl::geometric::RRTBispace::getPlannerData(base::PlannerData &data, int sub) const
{
    assert(sub >= 0);

    Planner::getPlannerData(data);

    std::vector<Motion *> motions;
    if (tStart_)
        tStart_->list(motions);

    for (auto & motion : motions)
    {
        if (motion->parent == nullptr)
            data.addStartVertex(base::PlannerDataVertex(motion->state->as<base::CompoundState>()->components[sub], 1));
        else
        {
            data.addEdge(base::PlannerDataVertex(motion->parent->state->as<base::CompoundState>()->components[sub], 1), base::PlannerDataVertex(motion->state->as<base::CompoundState>()->components[sub], 1));
        }
    }

    motions.clear();
    if (tGoal_)
        tGoal_->list(motions);

    for (auto & motion : motions)
    {
        if (motion->parent == nullptr)
            data.addGoalVertex(base::PlannerDataVertex(motion->state->as<base::CompoundState>()->components[sub], 2));
        else
        {
            // The edges in the goal tree are reversed to be consistent with start tree
            data.addEdge(base::PlannerDataVertex(motion->state->as<base::CompoundState>()->components[sub], 2), base::PlannerDataVertex(motion->parent->state->as<base::CompoundState>()->components[sub], 2));
        }
    }
}

void ompl::geometric::RRTBispace::getBiasData(base::PlannerData &data, bool start) const
{
    Planner::getPlannerData(data);

    if (start)
    {
        for (std::size_t i = 0; i < startBiasStates_.size(); i++)
        {
            if (i == 0)
                data.addStartVertex(base::PlannerDataVertex(startBiasStates_[i], 1));
            else 
                data.addEdge(base::PlannerDataVertex(startBiasStates_[i-1], 1), base::PlannerDataVertex(startBiasStates_[i], 1));
        }
    }
    else
    {
        for (std::size_t i = 0; i < goalBiasStates_.size(); i++)
        {
            if (i == 0)
                data.addStartVertex(base::PlannerDataVertex(goalBiasStates_[i], 2));
            else 
                data.addEdge(base::PlannerDataVertex(goalBiasStates_[i-1], 2), base::PlannerDataVertex(goalBiasStates_[i], 2));
        }
    }
}

/*
bool ompl::geometric::RRTBispace::selectCMotion(std::size_t &index, bool &reverse, double ratio1)
{
    index = 0;

    bool nconnect = false;

    std::size_t iter = 0;
    if (reverse)
        iter = connectionPoint_.size() - 1;

    while (iter < connectionPoint_.size())
    {
        auto pair = connectionPoint_[iter];

        base::Cost temp = opt_->combineCosts(pair.first->cost, pair.second->cost);

        if (opt_->isFinite(temp) && opt_->isCostBetterThan(temp, bestCost_))
        {
            double ratio = opt_->isFinite(bestCost_) ? std::abs((temp.value() - bestCost_.value()) / bestCost_.value()) : 1.0;

            if (ratio > ratio1 && pair.first != bestStartMotion_)
            {
                bool valid = true;
                if (pair.first->parent != nullptr)
                {
                    valid = isValid(pair.first->parent, true);
                }

                if (valid)
                {
                    valid = pair.first->valid;
                    if (!valid)
                    {
                        valid = checkMotion(pair.first->parent, pair.first, true);
                        if (!valid)
                        {
                            removeFromParent(pair.first);
                            pair.first->parent = nullptr;
                            valid = backPathRewireMotion(pair.first, true);
                            if (!valid)
                            {
                                pair.first->valid = true;
                                pnullStartMotions_.push_back(pair.first);
                            }
                        }
                    }
                }

                if (valid)
                {
                    if (pair.second->parent != nullptr)
                    {
                        valid = isValid(pair.second->parent, false);
                    }

                    if (valid)
                    {
                        valid = pair.second->valid;
                        if (!valid)
                        {
                            valid = checkMotion(pair.second->parent, pair.second, false);
                            if (!valid)
                            {
                                removeFromParent(pair.second);
                                pair.second->parent = nullptr;
                                valid = backPathRewireMotion(pair.second, false);
                                if (!valid)
                                {
                                    pair.second->valid = true;
                                    pnullGoalMotions_.push_back(pair.second);
                                }
                            }
                        }
                    }
                }

                if (valid)
                {
                    index = iter;
                    nconnect = true;
                    reverse = !reverse;
                    break;
                }
            }
        }

        if (reverse)
            iter--;
        else 
            iter++;
    }

    return nconnect;
}
*/
