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

#include "ompl/geometric/planners/bispace/RRTBispaceD.h"
#include "ompl/base/goals/GoalSampleableRegion.h"
#include "ompl/base/objectives/PathLengthOptimizationObjective.h"
#include "ompl/tools/config/SelfConfig.h"

ompl::geometric::RRTBispaceD::RRTBispaceD(const base::SpaceInformationPtr &si) : base::Planner(si, "RRTBispaceD")
  , dStart_(0), dGoal_(0), mc_(opt_), bh_(mc_)
{
    specs_.recognizedGoal = base::GOAL_SAMPLEABLE_REGION;
    specs_.directed = true;
    specs_.approximateSolutions = false;

    Planner::declareParam<double>("range", this, &RRTBispaceD::setRange, &RRTBispaceD::getRange, "0.:0.1:10000.");
    Planner::declareParam<double>("pen_distance", this, &RRTBispaceD::setPenDistance, &RRTBispaceD::getPenDistance, "0.:0.1:10000.");
    Planner::declareParam<bool>("lazy_path", this, &RRTBispaceD::setLazyPath, &RRTBispaceD::getLazyPath, "0,1");
    Planner::declareParam<bool>("lazy_node", this, &RRTBispaceD::setLazyNode, &RRTBispaceD::getLazyNode, "0,1");
    Planner::declareParam<bool>("add_intermediate_state", this, &RRTBispaceD::setAddIntermediateState, &RRTBispaceD::getAddIntermediateState, "0,1");
    Planner::declareParam<bool>("use_bispace", this, &RRTBispaceD::setUseBispace, &RRTBispaceD::getUseBispace, "0,1");
//    Planner::declareParam<bool>("use_biasgrow", this, &RRTBispaceD::setUseBiasGrow, &RRTBispaceD::getUseBiasGrow, "0,1");
//    Planner::declareParam<bool>("treated_as_multi_subapce", this, &RRTBispaceD::setTreatedAsMultiSubapce, &RRTBispaceD::getTreatedAsMultiSubapce, "0,1");
}

ompl::geometric::RRTBispaceD::~RRTBispaceD()
{
    freeMemory();
}

void ompl::geometric::RRTBispaceD::setup()
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

    if (minValidPathFraction_ < std::numeric_limits<double>::epsilon() || minValidPathFraction_ > 1.0)
        throw Exception("The minimum valid path fraction must be in the range (0,1]");

    sc.configureProjectionEvaluator(projectionEvaluator_);
    if (lazyPath_)
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
        lazyNode_ = false;
}

void ompl::geometric::RRTBispaceD::freeMemory()
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

    if (lazyPath_)
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

void ompl::geometric::RRTBispaceD::clear()
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

ompl::base::PlannerStatus ompl::geometric::RRTBispaceD::solve(const base::PlannerTerminationCondition &ptc)
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
                if (lazyPath_)
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
                if (isPathValid(pair.first, pair.second, !startTree))
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
        removeInvalidMotionsDirectly();
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

ompl::base::PlannerStatus ompl::geometric::RRTBispaceD::prepareSolve(const base::PlannerTerminationCondition &ptc)
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
            if (lazyPath_)
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
        if (lazyPath_)
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

void ompl::geometric::RRTBispaceD::processSolution(const Motion *bestStartMotion, const Motion *bestGoalMotion)
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

bool ompl::geometric::RRTBispaceD::batchGrow(bool &startTree)
{
    bool nconnect = false;
    TreeGrowingInfo tgi;
    tgi.xstate = si_->allocState();
    Motion *rmotion = new Motion(si_);
    base::State *rstate = rmotion->state;

    for (unsigned int i = 0; i < 10; i++)
    {
        tgi.start = startTree;
        startTree = !startTree;
        sampler_->sampleUniform(rstate);

        bool otherSide = false;
        bool change = false;
        bool gs = growTree(tgi, rmotion, otherSide, change);
        Motion *addedMotion = tgi.xmotion;
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

bool ompl::geometric::RRTBispaceD::growTree(TreeGrowingInfo &tgi, Motion *rmotion, bool &otherSide, bool &change)
{
    if (treatedAsMultiSubapce_)
        return growTreeMultiSpace(tgi, rmotion, otherSide, change);
    return growTreeSingleSpace(tgi, rmotion, otherSide, change);
}

bool ompl::geometric::RRTBispaceD::growTreeSingleSpace(TreeGrowingInfo &tgi, Motion *rmotion, bool &otherSide, bool &change)
{
    otherSide = false;
    change = false;
    TreeData &tree = tgi.start ? tStart_ : tGoal_;

    Motion *nmotion = tree->nearest(rmotion);
    tgi.xmotion = nmotion;
    if (si_->equalStates(nmotion->state, rmotion->state))
        return false;
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
        double d = tgi.start ? si_->distance(nmotion->state, rmotion->state) : si_->distance(rmotion->state, nmotion->state);
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
            if (tgi.start ? !checkInterMotion(nb, motion, true) : !checkInterMotion(motion, nb, false))
            {
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
        else
        {
            Coord xcoord(projectionEvaluator_->getDimension());
            projectionEvaluator_->computeCoordinates(motion->state, xcoord);
            addToDisc(disc, motion, xcoord);
            if (rewireSort_)
            {
                if (!motion->cell->data->cmotion || opt_->isCostBetterThan(motion->cost, motion->cell->data->cmotion->cost))
                    motion->cell->data->cmotion = motion;
                if (!motion->cell->data->mmotion || !opt_->isCostBetterThan(motion->cost, motion->cell->data->mmotion->cost))
                    motion->cell->data->mmotion = motion;
            }
        }
    }
    return reach;
}

bool ompl::geometric::RRTBispaceD::growTreeMultiSpace(TreeGrowingInfo &tgi, Motion *rmotion, bool &otherSide, bool &change)
{
    otherSide = false;
    change = false;
    TreeData &tree = tgi.start ? tStart_ : tGoal_;

    Motion *nmotion = tree->nearest(rmotion);
    tgi.xmotion = nmotion;
    if (si_->equalStates(nmotion->state, rmotion->state))
        return false;
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
            double d = tgi.start ? si_->distance(nmotion->state, rmotion->state) : si_->distance(rmotion->state, nmotion->state);
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
            if (tgi.start ? !checkInterMotion(nb, motion, true) : !checkInterMotion(motion, nb, false))
            {
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
        else
        {
            Coord xcoord(projectionEvaluator_->getDimension());
            projectionEvaluator_->computeCoordinates(motion->state, xcoord);
            addToDisc(disc, motion, xcoord);
            if (rewireSort_)
            {
                if (!motion->cell->data->cmotion || opt_->isCostBetterThan(motion->cost, motion->cell->data->cmotion->cost))
                    motion->cell->data->cmotion = motion;
                if (!motion->cell->data->mmotion || !opt_->isCostBetterThan(motion->cost, motion->cell->data->mmotion->cost))
                    motion->cell->data->mmotion = motion;
            }
        }
    }
    return reach && reachi;
}

bool ompl::geometric::RRTBispaceD::growCurrentTree(const base::State *state) const
{
    return growStartTree(state);
}

bool ompl::geometric::RRTBispaceD::growStartTree(const base::State *state) const
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

double ompl::geometric::RRTBispaceD::penetrationDistance(const base::State *nstate, const base::State *state, bool start) const
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

bool ompl::geometric::RRTBispaceD::growCurrentTree(const base::State *state, std::size_t sub) const
{
    return growStartTree(state, sub);
}

bool ompl::geometric::RRTBispaceD::growStartTree(const base::State *state, std::size_t sub) const
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

double ompl::geometric::RRTBispaceD::penetrationDistance(const base::State *nstate, const base::State *state, bool start, std::size_t sub) const
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

bool ompl::geometric::RRTBispaceD::isPathValid(Motion *motion, Motion *otherMotion, bool start)
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

bool ompl::geometric::RRTBispaceD::isPathValid(Motion *motion, Motion *otherMotion)
{
    bool valid = true;
    if (!isPathValid(motion, true))
        valid = false;
    if (!isPathValid(otherMotion, false))
        valid = false;
    return valid;
}

bool ompl::geometric::RRTBispaceD::isPathValid(Motion *motion, bool start)
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
    for (std::size_t i = mpath.size() - 1; i < mpath.size(); i--)
    {
        motion = mpath[i];
        Motion *pmotion = motion->parent;
        if (start ? !checkStartMotion(pmotion, motion) : !checkGoalMotion(motion, pmotion))
        {
            removeFromParent(motion);
            motion->parent = nullptr;
            tvalid = false;

            Motion *ppmotion = pmotion;
            tvalid = backPathRewireMotion(motion, start, ppmotion);
            if (tvalid)
            {
                tvalid = isPathValidInter(ppmotion, start);
                connectToPmotion(motion, ppmotion, start);
                motion->parent->children.push_back(motion); 
            }
            else 
                pnullMotions.push_back(motion);
            if (!tvalid)
                break;
        }
    }
    return tvalid;
}

bool ompl::geometric::RRTBispaceD::isStateValid(Motion *motion, bool start)  // todo
{
    bool tvalid = true;
    std::vector<Motion *> mpath;
    while (motion->parent)
    {
        mpath.push_back(motion);
        motion = motion->parent;
    }
    std::vector<Motion *> &pnullMotions= start ? pnullStartMotions_: pnullGoalMotions_;
    for (std::size_t i = mpath.size() - 1; i < mpath.size(); i--)
    {
        motion = mpath[i];
        if (!isValid(motion, start, true))
        {
            tvalid = false;
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
            if (last)
            {
                removeFromParent(last);
                last->parent = nullptr;
                Motion *plast = nullptr;
                tvalid = backPathRewireMotion(last, start, plast);
                if (tvalid)
                {
                    tvalid = isPathValidInter(plast, start);
                    connectToPmotion(last, plast, start);
                    last->parent->children.push_back(last);
                }
                else 
                    pnullMotions.push_back(last);
            }
            if (!tvalid)
                break;
        }
    }
    return tvalid;
}

bool ompl::geometric::RRTBispaceD::isPathValidInter(Motion *motion, bool start)
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
    for (std::size_t i = mpath.size() - 1; i < mpath.size(); i--)
    {
        motion = mpath[i];
        Motion *pmotion = motion->parent;
        if (start ? !checkStartMotion(pmotion, motion) : !checkGoalMotion(motion, pmotion))
        {
            removeFromParent(motion);
            motion->parent = nullptr;
            tvalid = false;

            Motion *ppmotion = pmotion;
            tvalid = backPathRewireMotion(motion, start, ppmotion);
            if (tvalid)
            {
                tvalid = isPathValidInter(ppmotion, start);
                connectToPmotion(motion, ppmotion, start);
                motion->parent->children.push_back(motion); 
            }
            else 
                pnullMotions.push_back(motion);
            if (!tvalid)
                break;
        }
    }
    return tvalid;
}

/*
bool ompl::geometric::RRTBispaceD::isPathValidInter(Motion *motion, bool start)
{
    if (lazyNode_ && !isStateValid(motion, start))
        return false;
    bool tvalid = true;
    std::vector<Motion *> &pnullMotions= start ? pnullStartMotions_: pnullGoalMotions_;
    std::vector<Motion *> &invalidMotions = start ? invalidStartMotions_ : invalidGoalMotions_;
    while (motion->parent)
    {
        Motion *pmotion = motion->parent;
        if (start ? !checkStartMotion(pmotion, motion) : !checkGoalMotion(motion, pmotion))
        {
            removeFromParent(motion);
            motion->parent = nullptr;
            tvalid = false;

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
            break;
        }
        motion = motion->parent;
    }
    return tvalid;
}
*/

bool ompl::geometric::RRTBispaceD::backPathRewireMotion(Motion *motion, bool start, Motion *&pmotion) // todo
{
    bool valid = false;
    CellDiscretizationData &disc = start ? dStart_ : dGoal_;
    for (auto & cv : motion->cell->nbh)
    {
        std::size_t numc = cv.size(), ic = numc - 1;
        while (ic < numc)
        {
            Cell *c = cv[ic];
            ic--;
            std::size_t num = c->data->motions.size();
            /*
            if (rewireSort_ && motion->cell != c && motion->cell->data->mmotion) // todo
            {
                if (opt_->isCostBetterThan(motion->cell->data->mmotion->cost, c->data->cmotion->cost))
                    continue;
                base::Cost cost1 = opt_->combineCosts(c->data->cmotion->cost, opt_->motionCost(c->data->cmotion->state, motion->state));
                base::Cost cost2 = opt_->combineCosts(motion->cell->data->mmotion->cost, opt_->combineCosts(motion->cell->data->mmotion->cost, motion->cell->data->mmotion->cost));
                if (opt_->isCostBetterThan(cost2, cost1))
                    continue;
            }
            */
            std::vector<Motion *> nbh;
            nbh.reserve(num);
            nbh.insert(nbh.end(), c->data->motions.begin(), c->data->motions.end());
            for (auto & nb : nbh)
            {
                if (!opt_->isFinite(nb->cost))
                    continue;
                if (!isValid(nb, start))
                    continue;
                if (nb == pmotion)
                    continue;
                if (start ? checkInterMotion(nb, motion, start) : checkInterMotion(motion, nb, start)) // todo
                    disc.updateNbh(nb->cell, motion->cell);
                else
                    continue;
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

void ompl::geometric::RRTBispaceD::removeInvalidMotionsDirectly()
{
    for (auto & pair : connectionPoint_)
    {
        pair.first->inConnection = false;
        pair.second->inConnection = false;
    }
    connectionPoint_.clear();
    removeInvalidMotionsDirectlyTree();
}

void ompl::geometric::RRTBispaceD::removeInvalidMotionsDirectlyTree()
{
    if (!invalidStartMotions_.empty() || !pnullStartMotions_.empty())
    {
        tStart_->clear();
        for (auto & rootMotion : startMotions_)
            addToTree(tStart_, rootMotion);
        for (auto & pnull : invalidStartMotions_)
            removeFromTree(dStart_, startBiasPdf_, pnull);
        for (auto & pnull : pnullStartMotions_)
            removeFromTree(dStart_, startBiasPdf_, pnull);
        invalidStartMotions_.clear();
        pnullStartMotions_.clear();
    }

    if (!invalidGoalMotions_.empty() || !pnullGoalMotions_.empty())
    {
        tGoal_->clear();
        for (auto & rootMotion : goalMotions_)
            addToTree(tGoal_, rootMotion);
        for (auto & pnull : invalidGoalMotions_)
            removeFromTree(dGoal_, goalBiasPdf_, pnull);
        for (auto & pnull : pnullGoalMotions_)
            removeFromTree(dGoal_, goalBiasPdf_, pnull);
        invalidGoalMotions_.clear();
        pnullGoalMotions_.clear();
    }
}

void ompl::geometric::RRTBispaceD::addToTree(TreeData &tree, Motion *motion)
{
    tree->add(motion);
    for (auto & child : motion->children)
        addToTree(tree, child);
}

void ompl::geometric::RRTBispaceD::removeFromTree(CellDiscretizationData &disc, MotionPDF &pdf, Motion *motion)
{
//    if (useBiasGrow_ && motion->middle && motion->stateValid == Valid)
//        removePdfMotion(pdf, motion);
    for (auto & child : motion->children)
    {
        child->parent = nullptr;
        removeFromTree(disc, pdf, child);
    }
    removeFromDisc(disc, motion);
    freeMotion(motion);
}

void ompl::geometric::RRTBispaceD::addPdfMotion(MotionPDF &pdf, Motion *motion, bool start)
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

ompl::geometric::RRTBispaceD::Motion *ompl::geometric::RRTBispaceD::selectPdfMotion(MotionPDF &pdf, GridCell *&cell)
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

void ompl::geometric::RRTBispaceD::removePdfMotion(MotionPDF &pdf, Motion *motion)
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

void ompl::geometric::RRTBispaceD::connectToPmotion(Motion *motion, Motion *pmotion, bool start) const
{
    motion->parent = pmotion;
    motion->root = pmotion->root;
    motion->incCost = start ? opt_->motionCost(pmotion->state, motion->state) : opt_->motionCost(motion->state, pmotion->state);
    motion->cost = opt_->combineCosts(pmotion->cost, motion->incCost);
    if (opt_->isFinite(motion->cost))
        updateChildCosts(motion);
}

void ompl::geometric::RRTBispaceD::updateChildCosts(Motion *motion) const
{
    for (auto & child : motion->children)
    {
        child->cost = opt_->combineCosts(motion->cost, child->incCost);
        child->root = motion->root;
        updateChildCosts(child);
    }
}

void ompl::geometric::RRTBispaceD::addToDisc(CellDiscretizationData &disc, Motion *motion, const Coord& coord)
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

void ompl::geometric::RRTBispaceD::removeFromDisc(CellDiscretizationData &disc, Motion *motion)
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

void ompl::geometric::RRTBispaceD::removeFromParent(Motion *motion)
{
    if (!motion->parent)
        return;
    removeFromVector(motion->parent->children, motion);
}

bool ompl::geometric::RRTBispaceD::removeFromVector(std::vector<Motion *> &motions, Motion *motion)
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

void ompl::geometric::RRTBispaceD::setMotionInfinityCost(Motion *motion) const
{
    motion->cost = opt_->infiniteCost();
    for (auto & child : motion->children)
        setMotionInfinityCost(child);
}

// check 
bool ompl::geometric::RRTBispaceD::isValid(const base::State *state)
{
    bool valid = si_->isValid(state);
    return valid;
}

bool ompl::geometric::RRTBispaceD::isValid(Motion *motion, bool start, bool add)
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

        if (add && addIntermediateState_ && motion->parent && motion->parent->stateValid == Valid)
        {
            bool intervalid = false;
            if (start)
                intervalid = checkInterMotion2(motion->parent, motion, start);
            else 
                intervalid = checkInterMotion2(motion, motion->parent, start);
            if (intervalid)
            {
                if (start)
                {
                    int nd = si_->getStateSpace()->validSegmentCount(motion->parent->state, motion->state);
                    if (nd >= 2)
                    {
                        Motion *last = new Motion(si_);
                        si_->getStateSpace()->interpolate(motion->parent->state, motion->state, (double)(nd-1) / (double)nd, last->state);
                        addIntermediateMotion(motion->parent, motion, start, last);
                    }
                }
                else 
                {
                    int nd = si_->getStateSpace()->validSegmentCount(motion->state, motion->parent->state);
                    if (nd >= 2)
                    {
                        Motion *last = new Motion(si_);
                        si_->getStateSpace()->interpolate(motion->state, motion->parent->state, 1.0 / (double)nd, last->state);
                        addIntermediateMotion(motion, motion->parent, start, last);
                    }
                }
            }
        }

        motion->stateValid = InValid;
        removeFromParent(motion);
        motion->parent = nullptr;
        if (start)
            invalidStartMotions_.push_back(motion);
        else 
            invalidGoalMotions_.push_back(motion);
        if (opt_->isFinite(motion->cost))
            setMotionInfinityCost(motion);
    }
    return motion->stateValid == Valid;
}

bool ompl::geometric::RRTBispaceD::checkStartMotion(Motion *smotion, Motion *gmotion)
{
    if (!gmotion->valid)
    {
        if (!checkInterMotion(smotion, gmotion, true))
        {
            if (opt_->isFinite(gmotion->cost))
                setMotionInfinityCost(gmotion);
        }
        else 
        {
            gmotion->valid = true;
            dStart_.updateNbh(gmotion->cell, smotion->cell);
        }
    }
    return gmotion->valid;
}

bool ompl::geometric::RRTBispaceD::checkGoalMotion(Motion *smotion, Motion *gmotion)
{
    if (!smotion->valid)
    {
        if (!checkInterMotion(smotion, gmotion, false))
        {
            if (opt_->isFinite(smotion->cost))
                setMotionInfinityCost(smotion);
        }
        else 
        {
            smotion->valid = true;
            dGoal_.updateNbh(gmotion->cell, smotion->cell);
        }
    }
    return smotion->valid;
}

bool ompl::geometric::RRTBispaceD::checkInterMotion(Motion *smotion, Motion *gmotion, bool start)
{
    bool valid = false;
    if (addIntermediateState_)
        valid = checkInterMotion2(smotion, gmotion, start);
    else
        valid = checkInterMotion1(smotion, gmotion, start);
    return valid;
}

bool ompl::geometric::RRTBispaceD::checkInterMotion1(Motion *smotion, Motion *gmotion, bool /*start*/)
{
    /*assume smotion, gmotion are valid*/
    bool valid = true;
    const base::State *s1 = smotion->state, *s2 = gmotion->state;
    int nd = si_->getStateSpace()->validSegmentCount(s1, s2);
    if (nd >= 2)
    {
        std::queue<std::pair<int, int>> pos;
        pos.emplace(1, nd - 1);

        base::State *state= si_->allocState();

        /* repeatedly subdivide the path segment in the middle (and check the middle) */
        while (!pos.empty())
        {
            std::pair<int, int> x = pos.front();
            int mid = (x.first + x.second) / 2;
            si_->getStateSpace()->interpolate(s1, s2, (double)mid / (double)nd, state);

            if (!si_->isValid(state))
            {
                valid = false;
                break;
            }

            pos.pop();

            if (x.first < mid)
                pos.emplace(x.first, mid - 1);
            if (x.second > mid)
                pos.emplace(mid + 1, x.second);
        }

        si_->freeState(state);
    }
    return valid;
}

bool ompl::geometric::RRTBispaceD::checkInterMotion2(Motion *smotion, Motion *gmotion, bool start)
{
    bool valid = true;
    /*assume smotion, gmotion are valid*/
    const base::State *s1 = smotion->state, *s2 = gmotion->state;
    int nd = si_->getStateSpace()->validSegmentCount(s1, s2);
    if (nd >= 2)
    {
        base::State *state= si_->allocState();

        int i = start ? 1 : nd - 1;
        while (start ? i < nd : i > 0)
        {
            si_->getStateSpace()->interpolate(s1, s2, (double)i / (double)nd, state);
            if (!si_->isValid(state))
            {
                valid = false;
                break;
            }
            i = start ? i + 1: i - 1;
        }

        if (!valid)
        {
            i = start ? i - 1: i + 1;
            double ratio = (double)i/(double)nd;
            if (start ? ratio > minValidPathFraction_ : ratio < 1.0 - minValidPathFraction_)
            {
                Motion *last = new Motion(si_);
                si_->getStateSpace()->interpolate(s1, s2, (double)i / (double)nd, last->state);
                addIntermediateMotion(smotion, gmotion, start, last);
            }
        }

        si_->freeState(state);
    }
    return valid;
}

void ompl::geometric::RRTBispaceD::addIntermediateMotion(Motion *smotion, Motion *gmotion, bool start, Motion *last)
{
    TreeData &tree = start ? tStart_ : tGoal_;
    last->valid = true;
    last->stateValid = Valid;
    Motion *v = start ? smotion : gmotion;
    connectToPmotion(last, v, start);
    last->parent->children.push_back(last);
    tree->add(last);

    if (lazyPath_)
    {
        CellDiscretizationData &disc = start ? dStart_ : dGoal_;
        Coord xcoord(projectionEvaluator_->getDimension());
        projectionEvaluator_->computeCoordinates(last->state, xcoord);
        addToDisc(disc, last, xcoord);
        disc.updateNbh(last->cell, v->cell);
        if (rewireSort_)
        {
            if (!last->cell->data->cmotion || opt_->isCostBetterThan(last->cost, last->cell->data->cmotion->cost))
                last->cell->data->cmotion = last;
            if (!last->cell->data->mmotion || !opt_->isCostBetterThan(last->cost, last->cell->data->mmotion->cost))
                last->cell->data->mmotion = last;
        }
    }
}

void ompl::geometric::RRTBispaceD::getPlannerData(base::PlannerData &data) const
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

void ompl::geometric::RRTBispaceD::getPlannerData(base::PlannerData &data, int sub) const
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

void ompl::geometric::RRTBispaceD::getBiasData(base::PlannerData &data, bool start) const
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
