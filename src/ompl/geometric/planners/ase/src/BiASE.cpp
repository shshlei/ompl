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

#include "ompl/geometric/planners/ase/BiASE.h"
#include "ompl/base/goals/GoalSampleableRegion.h"
#include "ompl/base/objectives/PathLengthOptimizationObjective.h"
#include "ompl/tools/config/SelfConfig.h"

#include "ompl/base/samplers/adinformed/RejectionAdInfSampler.h"
#include "ompl/base/samplers/adinformed/PathLengthDirectAdInfSampler.h"

ompl::geometric::BiASE::BiASE(const base::SpaceInformationPtr &si) : base::Planner(si, "BiASE")
  , dStart_(0), dGoal_(0), mc_(opt_), bh_(mc_)
{
    specs_.recognizedGoal = base::GOAL_SAMPLEABLE_REGION;
    specs_.directed = true;
    specs_.approximateSolutions = false;

    setInformedSampling(true);
    setSampleRejection(false);

    Planner::declareParam<double>("range", this, &BiASE::setRange, &BiASE::getRange, "0.:0.1:10000.");
    Planner::declareParam<double>("pen_distance", this, &BiASE::setPenDistance, &BiASE::getPenDistance, "0.:0.1:10000.");
    Planner::declareParam<bool>("use_bispace", this, &BiASE::setUseBispace, &BiASE::getUseBispace, "0,1");
//    Planner::declareParam<bool>("use_biasgrow", this, &BiASE::setUseBiasGrow, &BiASE::getUseBiasGrow, "0,1");
//    Planner::declareParam<bool>("treated_as_multi_subapce", this, &BiASE::setTreatedAsMultiSubapce, &BiASE::getTreatedAsMultiSubapce, "0,1");
    Planner::declareParam<bool>("update_nbhcell", this, &BiASE::setUpdateNbCell, &BiASE::getUpdateNbCell, "0,1");
    Planner::declareParam<bool>("add_inter_states", this, &BiASE::setAddIntermediateState, &BiASE::getAddIntermediateState, "0,1");

    Planner::declareParam<bool>("informed_sampling", this, &BiASE::setInformedSampling, &BiASE::getInformedSampling, "0,1");
    Planner::declareParam<bool>("sample_rejection", this, &BiASE::setSampleRejection, &BiASE::getSampleRejection, "0,1");
    Planner::declareParam<unsigned int>("number_sampling_attempts", this, &BiASE::setNumSamplingAttempts,
                                        &BiASE::getNumSamplingAttempts, "10:10:100000");
}

ompl::geometric::BiASE::~BiASE()
{
    freeMemory();
}

void ompl::geometric::BiASE::setup()
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
    dStart_.setDimension(projectionEvaluator_->getDimension());
    dGoal_.setDimension(projectionEvaluator_->getDimension());
    dStart_.setNeighborCell(2);
    dGoal_.setNeighborCell(2);

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
    symmetric_ = si_->getStateSpace()->hasSymmetricDistance() && si_->getStateSpace()->hasSymmetricInterpolate();
}

void ompl::geometric::BiASE::freeMemory()
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

void ompl::geometric::BiASE::clear()
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

    checkedStartPath_.clear();
    checkedGoalPath_.clear();

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

    clearStartAdInfSampler();
    clearGoalAdInfSampler();
    startAdInfProb_ = 0.0;
    tree_ = -1;
    localRatio_ = 0.75;

    ais_ = false;
}

ompl::base::PlannerStatus ompl::geometric::BiASE::solve(const base::PlannerTerminationCondition &ptc)
{
    base::PlannerStatus ps = prepareSolve(ptc);
    if (ps != base::PlannerStatus::PREPARE_SUCCESS)
        return ps;

    clearStartAdInfSampler();
    clearGoalAdInfSampler();
    sampler_ = si_->allocStateSampler();

    OMPL_INFORM("%s: Starting planning with %u states already in datastructure.", getName().c_str(), (tStart_->size() + tGoal_->size()));

    bool startTree = true;
    bool solved = false;
    Motion *bestStartMotion = nullptr;
    Motion *bestGoalMotion = nullptr;

    Motion *startAd = nullptr, *goalAd = nullptr;
    bool reverse = false;
    unsigned int adinfcount = 0;

    while (!ptc && !solved)
    {
        if (pis_.getSampledGoalsCount() < tGoal_->size() / 2)
        {
            const base::State *st = pis_.nextGoal();
            if (st != nullptr)
            {
                auto *motion = new Motion(si_);
                si_->copyState(motion->state, st);
                motion->valid = true;
                motion->root = motion->state;
                motion->cost = opt_->identityCost();
                goalMotions_.push_back(motion);
                tGoal_->add(motion);

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

        if (!batchGrow(startTree))
            continue;

        if (!startAd)
        {
            std::size_t index = 0;
            if (selectCMotion(index, reverse))
            {
                auto pair = connectionPoint_[index];
                startAd = pair.first;
                goalAd = pair.second;
            }
        }

        bool clearoradd = false;
        if (startAd && opt_->isFinite(startAd->cost) && opt_->isFinite(goalAd->cost))
        {
            if (isPathValid(startAd, goalAd))
            {
                bestStartMotion = startAd;
                bestGoalMotion = goalAd;
                solved = true;
            }
            else if (!ais_ && !checkedStartPath_.empty() && !checkedGoalPath_.empty())
            {
                bool locals = false, localg = false;
                localInfeasible(tree_, locals, localg);
                if (tree_ == 0)
                {
                    if (!goalAdInfSamplers_.empty())
                    {
                        clearoradd = true;
                        clearGoalAdInfSampler();
                    }
                    calculateInfSampler(locals, true, clearoradd);
                }
                else if (tree_ == 1)
                {
                    if (!startAdInfSamplers_.empty())
                    {
                        clearoradd = true;
                        clearStartAdInfSampler();
                    }
                    calculateInfSampler(localg, false, clearoradd);
                }
                else 
                {
                    calculateInfSampler(locals, true, clearoradd);
                    calculateInfSampler(localg, false, clearoradd);
                }
            }
        }

        if (!solved)
        {
            for (auto & pair : connectionPoint_)
            {
                if (opt_->isFinite(pair.first->cost) && opt_->isFinite(pair.second->cost))
                {
                    if (isPathValid(pair.first, pair.second))
                    {
                        bestStartMotion = pair.first;
                        bestGoalMotion = pair.second;
                        solved = true;
                        break;
                    }
                }
            }
        }

        if (solved)
            break;

        processAdEllipsoidRind(clearoradd, adinfcount);
        if (!ais_) 
            startAd = goalAd = nullptr;
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

ompl::base::PlannerStatus ompl::geometric::BiASE::prepareSolve(const base::PlannerTerminationCondition &ptc)
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
            motion->root = motion->state;
            motion->cost = opt_->identityCost();
            startMotions_.push_back(motion);
            tStart_->add(motion);

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
        motion->root = motion->state;
        motion->cost = opt_->identityCost();
        goalMotions_.push_back(motion);
        tGoal_->add(motion);

        Coord xcoord(projectionEvaluator_->getDimension());
        projectionEvaluator_->computeCoordinates(motion->state, xcoord);
        addToDisc(dGoal_, motion, xcoord);
        if (rewireSort_)
        {
            motion->cell->data->cmotion = motion;
            motion->cell->data->mmotion = motion;
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

void ompl::geometric::BiASE::processSolution(const Motion *bestStartMotion, const Motion *bestGoalMotion)
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

bool ompl::geometric::BiASE::growTree(TreeGrowingInfo &tgi, Motion *rmotion, bool &otherSide, bool &change, bool &add)
{
    if (treatedAsMultiSubapce_)
        return growTreeMultiSpace(tgi, rmotion, otherSide, change, add);
    return growTreeSingleSpace(tgi, rmotion, otherSide, change, add);
}

bool ompl::geometric::BiASE::growTreeSingleSpace(TreeGrowingInfo &tgi, Motion *rmotion, bool &otherSide, bool &change, bool &add)
{
    otherSide = false;
    change = false;
    TreeData &tree = tgi.start ? tStart_ : tGoal_;

    bool null = false;
    Motion *nmotion = selectNMotion(tree, rmotion, null);
    tgi.xmotion = nmotion;
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

        if (!isValid(dstate)) // todo
        {
            reach = false;
            break;
        }

        motion = new Motion(si_);
        si_->copyState(motion->state, dstate);

        Motion *nb = nmotion;
        connectToPmotion(motion, nb, tgi.start);
        motion->parent->children.push_back(motion);
        tree->add(motion);
        tgi.xmotion = motion;
        add = true;
        if (addpd)
        {
            otherSide = true;
            motion->middle = true;
            /*
            if (useBiasGrow_)
            {
                if (tgi.start)
                    addPdfMotion(startBiasPdf_, motion, true);
                else 
                    addPdfMotion(goalBiasPdf_, motion, false);
            }
            */
        }

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
    return reach;
}

bool ompl::geometric::BiASE::growTreeMultiSpace(TreeGrowingInfo &tgi, Motion *rmotion, bool &otherSide, bool &change, bool &add)
{
    otherSide = false;
    change = false;
    TreeData &tree = tgi.start ? tStart_ : tGoal_;

    bool null = false;
    Motion *nmotion = selectNMotion(tree, rmotion, null);
    tgi.xmotion = nmotion;
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

        if (!isValid(dstate))
        {
            reach = false;
            break;
        }

        motion = new Motion(si_);
        si_->copyState(motion->state, dstate);

        Motion *nb = nmotion;
        connectToPmotion(motion, nb, tgi.start);
        motion->parent->children.push_back(motion);
        tree->add(motion);
        tgi.xmotion = motion;
        add = true;
        if (otherSide)
        {
            motion->middle = true;
            /*
            if (useBiasGrow_)
            {
                if (tgi.start)
                    addPdfMotion(startBiasPdf_, motion, true);
                else 
                    addPdfMotion(goalBiasPdf_, motion, false);
            }
            */
        }

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
    return reach && reachi;
}

ompl::geometric::BiASE::Motion *ompl::geometric::BiASE::selectNMotion(const TreeData &tree, Motion *rmotion, bool &null)
{
    null = false;
    Motion *nmotion = tree->nearest(rmotion);
    if (!opt_->isFinite(nmotion->cost))
    {
        null = true;
        Motion *last = nmotion;
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

ompl::geometric::BiASE::Motion *ompl::geometric::BiASE::selectMotionInCell(Cell *cell)
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

void ompl::geometric::BiASE::rewireTree(BinaryHeap<Motion *, MotionCompare> &bh, Motion *m, bool start)
{
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
                if (ic)
                    std::sort(cv.begin(), cv.begin() + ic + 1, ocbd);
                Cell *c = cv[ic];
                ic--;
                std::size_t pos = c->data->disabled;
                if (!pos)
                    break;
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

void ompl::geometric::BiASE::updateQueue(BinaryHeap<Motion *, MotionCompare> &bh, Motion *motion)
{
    if (motion->handle != nullptr)
        bh.update(motion->handle);
    else
        motion->handle = bh.insert(motion);
}

bool ompl::geometric::BiASE::growCurrentTree(const base::State *state) const
{
    return growStartTree(state);
}

bool ompl::geometric::BiASE::growStartTree(const base::State *state) const
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

double ompl::geometric::BiASE::penetrationDistance(const base::State *nstate, const base::State *state, bool start) const
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

bool ompl::geometric::BiASE::growCurrentTree(const base::State *state, std::size_t sub) const
{
    return growStartTree(state, sub);
}

bool ompl::geometric::BiASE::growStartTree(const base::State *state, std::size_t sub) const
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

double ompl::geometric::BiASE::penetrationDistance(const base::State *nstate, const base::State *state, bool start, std::size_t sub) const
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

bool ompl::geometric::BiASE::backPathRewireMotion(Motion *motion, bool start, Motion *&pmotion)
{
    bool valid = false;
    OrderCellsByCost ocbc(opt_);
    CellDiscretizationData &disc = start ? dStart_ : dGoal_;
    std::size_t cnbh = motion->cell->nbh.size();
    for (std::size_t cindex = 0; cindex < cnbh; cindex++)
    {
        std::size_t numc = motion->cell->nbh[cindex].size(), ic = numc - 1;
        while (ic < numc)
        {
            if (ic)
                std::sort(motion->cell->nbh[cindex].begin(), motion->cell->nbh[cindex].begin() + ic + 1, ocbc);
            Cell *c = motion->cell->nbh[cindex][ic];
            ic--;
            std::size_t num = c->data->motions.size() - c->data->disabled;
            if (!num)
                break;
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
                if (isInvalidNeighbor(motion, nb))
                    continue;
                if (!isValidNeighbor(motion, nb))
                {
                    if (checkInterMotion(nb, motion, start))
                    {
                        insertNeighbor(nb, motion);
                        if (updateNbCell_)
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

void ompl::geometric::BiASE::addPdfMotion(MotionPDF &pdf, Motion *motion, bool start)
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
}

ompl::geometric::BiASE::Motion *ompl::geometric::BiASE::selectPdfMotion(MotionPDF &pdf, GridCell *&cell)
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

void ompl::geometric::BiASE::removePdfMotion(MotionPDF &pdf, Motion *motion)
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

void ompl::geometric::BiASE::enableMotionInDisc(Motion *motion)
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

void ompl::geometric::BiASE::connectToPmotion(Motion *motion, Motion *pmotion, bool start) const
{
    motion->parent = pmotion;
    motion->root = pmotion->root;
    motion->incCost = start ? opt_->motionCost(pmotion->state, motion->state) : opt_->motionCost(motion->state, pmotion->state);
    motion->cost = opt_->combineCosts(pmotion->cost, motion->incCost);
    if (opt_->isFinite(motion->cost))
        updateChildCosts(motion);
}

void ompl::geometric::BiASE::updateChildCosts(Motion *motion) const
{
    for (auto & child : motion->children)
    {
        child->cost = opt_->combineCosts(motion->cost, child->incCost);
        child->root = motion->root;
        updateChildCosts(child);
    }
}

void ompl::geometric::BiASE::insertNeighbor(Motion *pmotion, Motion *motion)
{
    pmotion->nbh.insert(motion);
    motion->nbh.insert(pmotion);
}

void ompl::geometric::BiASE::insertInvalidNeighbor(Motion *pmotion, Motion *motion)
{
    pmotion->invalidnbh.insert(motion);
    motion->invalidnbh.insert(pmotion);
}

void ompl::geometric::BiASE::removeFromNeighbor(Motion *motion)
{
    for (auto & pmotion : motion->nbh)
        pmotion->nbh.erase(motion);
}

void ompl::geometric::BiASE::removeFromInvalidNeighbor(Motion *motion)
{
    for (auto & pmotion : motion->invalidnbh)
        pmotion->invalidnbh.erase(motion);
}

bool ompl::geometric::BiASE::isValidNeighbor(Motion *motion, Motion *nb) const
{
    return motion->nbh.find(nb) != motion->nbh.end();
}

bool ompl::geometric::BiASE::isInvalidNeighbor(Motion *motion, Motion *nb) const
{
    return motion->invalidnbh.find(nb) != motion->invalidnbh.end();
}

void ompl::geometric::BiASE::addToDisc(CellDiscretizationData &disc, Motion *motion, const Coord& coord)
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

void ompl::geometric::BiASE::removeFromDisc(CellDiscretizationData &disc, Motion *motion)
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

void ompl::geometric::BiASE::removeFromParent(Motion *motion)
{
    if (!motion->parent)
        return;
    removeFromVector(motion->parent->children, motion);
}

void ompl::geometric::BiASE::removeFromPnull(std::vector<Motion *> &pnullMotions, Motion *motion)
{
    if (removeFromVector(pnullMotions, motion))
        motion->pmotion = nullptr;
}

bool ompl::geometric::BiASE::removeFromVector(std::vector<Motion *> &motions, Motion *motion)
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

void ompl::geometric::BiASE::setMotionInfinityCostWithDisable(Motion *motion, std::unordered_set<Cell *> &cells) const
{
    motion->cost = opt_->infiniteCost();
    motion->cell->data->disabled++;
    if (rewireSort_)
        cells.insert(motion->cell);
    for (auto & child : motion->children)
        setMotionInfinityCostWithDisable(child, cells);
}

// check 
bool ompl::geometric::BiASE::isValid(const base::State *state)
{
    bool valid = si_->isValid(state);
    return valid;
}

bool ompl::geometric::BiASE::checkMotion(Motion *pmotion, Motion *motion, bool start)
{
    if (!motion->valid)
    {
        if (!checkInterMotion(pmotion, motion, start))
        {
            if (opt_->isFinite(motion->cost))
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
            insertInvalidNeighbor(pmotion, motion);
        }
        else 
        {
            motion->valid = true;
            insertNeighbor(pmotion, motion);
            CellDiscretizationData &disc = start ? dStart_ : dGoal_;
            if (updateNbCell_)
                disc.updateNbh(pmotion->cell, motion->cell);
        }
    }
    return motion->valid;
}

bool ompl::geometric::BiASE::checkInterMotion(Motion *pmotion, Motion *motion, bool start)
{
    bool valid = false;
    if (start || symmetric_)
        valid = checkInterMotion1(pmotion, motion, start);
    else 
        valid = checkInterMotion2(motion, pmotion, start);
    return valid;
}

bool ompl::geometric::BiASE::checkInterMotion1(Motion *smotion, Motion *gmotion, bool start)
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

bool ompl::geometric::BiASE::checkInterMotion2(Motion *smotion, Motion *gmotion, bool start)
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

void ompl::geometric::BiASE::addIntermediateMotion(Motion *pmotion, Motion *motion, bool start, Motion *last)
{
    TreeData &tree = start ? tStart_ : tGoal_;
    last->valid = true;
    connectToPmotion(last, pmotion, start);
    last->parent->children.push_back(last);
    tree->add(last);

    CellDiscretizationData &disc = start ? dStart_ : dGoal_;
    Coord xcoord(projectionEvaluator_->getDimension());
    projectionEvaluator_->computeCoordinates(last->state, xcoord);
    addToDisc(disc, last, xcoord);
    if (!opt_->isFinite(last->cost))
        last->cell->data->disabled++;
    insertNeighbor(last, pmotion);
    if (updateNbCell_)
        disc.updateNbh(last->cell, pmotion->cell);
    insertInvalidNeighbor(last, motion);
    if (rewireSort_)
    {
        if (!last->cell->data->cmotion || opt_->isCostBetterThan(last->cost, last->cell->data->cmotion->cost))
            last->cell->data->cmotion = last;
        if (!last->cell->data->mmotion || opt_->isCostBetterThan(last->cell->data->mmotion->cost, last->cost))
            last->cell->data->mmotion = last;
    }
}

void ompl::geometric::BiASE::getPlannerData(base::PlannerData &data) const
{
    Planner::getPlannerData(data);

    std::vector<Motion *> motions;
    if (tStart_)
        tStart_->list(motions);
    for (auto & motion : motions)
    {
        if (!opt_->isFinite(motion->cost))
            continue;
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
        if (!opt_->isFinite(motion->cost))
            continue;
        if (motion->parent == nullptr)
            data.addGoalVertex(base::PlannerDataVertex(motion->state, 2));
        else
        {
            // The edges in the goal tree are reversed to be consistent with start tree
            data.addEdge(base::PlannerDataVertex(motion->state, 2), base::PlannerDataVertex(motion->parent->state, 2));
        }
    }
}

// adaptive informed sampling
bool ompl::geometric::BiASE::batchGrow(bool &startTree)
{
    bool add = false; // add a new random sample
    TreeGrowingInfo tgi;
    tgi.xstate = si_->allocState();
    Motion *rmotion = new Motion(si_);
    base::State *rstate = rmotion->state;

    for (unsigned int i = 0; i < 1; i++)
    {
        if (tree_ == 0)
            startTree = true;
        else if (tree_ == 1)
            startTree = false;
        else if (startAdInfProb_ > 0.0)
        {
            if (rng_.uniform01() < startAdInfProb_)
                startTree = true;
            else 
                startTree = false;
        }

        tgi.start = startTree;
        startTree = !startTree;
        if (!sampleUniformAd(rstate, tgi.start)) 
            continue;

        bool otherSide = false;
        bool change = false;
        bool gs = growTree(tgi, rmotion, otherSide, change, add);
        Motion *addedMotion = tgi.xmotion;
        if (!gs)
            si_->copyState(rstate, addedMotion->state);
        tgi.start = startTree;
        bool gsc = growTree(tgi, rmotion, otherSide, change, add);
        Motion *startMotion = nullptr, *goalMotion = nullptr;
        if (gsc && !change)
        {
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
            gsc = growTree(tgi, rmotion, otherSide, change, add);
            if (gsc && !change)
            {
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
            startMotion->inConnection = true;
            goalMotion->inConnection = true;
            connectionPoint_.emplace_back(startMotion, goalMotion);
        }
    }

    si_->freeState(tgi.xstate);
    freeMotion(rmotion);
    return add;
}

bool ompl::geometric::BiASE::selectCMotion(std::size_t &index, bool &reverse)
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
        if (opt_->isFinite(temp))
        {
            index = iter;
            nconnect = true;
            reverse = !reverse;
            break;
        }
        if (reverse)
            iter--;
        else 
            iter++;
    }
    return nconnect;
}

bool ompl::geometric::BiASE::isPathValid(Motion *motion, Motion *otherMotion)
{
    bool valid = true;
    checkedStartPath_.clear();
    checkedGoalPath_.clear();
    if (!isPathValid(motion, true))
        valid = false;
    if (!isPathValid(otherMotion, false))
        valid = false;
    return valid;
}

bool ompl::geometric::BiASE::isPathValid(Motion *motion, bool start)
{
    bool valid = true, tvalid = true;
    std::vector<Motion *> mpath;
    std::vector<Motion *> &checkedPath = start ? checkedStartPath_ : checkedGoalPath_;
    while (motion)
    {
        checkedPath.push_back(motion);
        mpath.push_back(motion);
        motion = motion->parent;
    }
    mpath.pop_back();
    std::vector<Motion *> &pnullMotions= start ? pnullStartMotions_: pnullGoalMotions_;
    for (std::size_t i = mpath.size() - 1; i < mpath.size(); i--)
    {
        motion = mpath[i];
        Motion *pmotion = motion->parent;
        if (!checkMotion(pmotion, motion, start))
        {
            removeFromParent(motion);
            motion->parent = nullptr;
            if (tvalid)
            {
                valid = false;
                if (backRewire_)
                {
                    Motion *ppmotion = nullptr;
                    valid = backPathRewireMotion(motion, start, ppmotion);
                    if (valid)
                    {
                        checkedPath.resize(i+1);
                        Motion *last = ppmotion;
                        while (last)
                        {
                            checkedPath.push_back(last);
                            last = last->parent;
                        }
                        valid = isPathValidInter(ppmotion, start);
                        connectToPmotion(motion, ppmotion, start);
                        motion->parent->children.push_back(motion); 
                    }
                }
                tvalid = valid;
            }
            if (tvalid)
                enableMotionInDisc(motion);
            else if (!motion->parent)
            {
                pnullMotions.push_back(motion);
                motion->valid = false;
                motion->pmotion = pmotion;
            }
        }
    }
    return tvalid;
}

bool ompl::geometric::BiASE::isPathValidInter(Motion *motion, bool start)
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
        Motion *pmotion = motion->parent;
        if (!checkMotion(pmotion, motion, start))
        {
            removeFromParent(motion);
            motion->parent = nullptr;
            tvalid = false;
            pnullMotions.push_back(motion);
            motion->valid = false;
            motion->pmotion = pmotion;
        }
    }
    return tvalid;
}

/*
bool ompl::geometric::BiASE::isPathValidInter(Motion *motion, bool start) // back rewire
{
    bool valid = true, tvalid = true;
    std::vector<Motion *> mpath;
    while (motion->parent)
    {
        mpath.push_back(motion);
        motion = motion->parent;
    }
    std::vector<Motion *> &pnullMotions= start ? pnullStartMotions_: pnullGoalMotions_;
    std::vector<Motion *> nullMotions;
    for (std::size_t i = mpath.size() - 1; i < mpath.size(); i--)
    {
        motion = mpath[i];
        Motion *pmotion = motion->parent;
        if (!checkMotion(pmotion, motion, start))
        {
            removeFromParent(motion);
            motion->parent = nullptr;
            if (tvalid)
            {
                valid = false;
                if (backRewire_)
                {
                    Motion *ppmotion = nullptr;
                    valid = backPathRewireMotion(motion, start, ppmotion);
                    if (valid)
                    {
                        valid = isPathValidInter(ppmotion, start);
                        connectToPmotion(motion, ppmotion, start);
                        motion->parent->children.push_back(motion); 
                    }
                }
                tvalid = valid;
            }
            if (tvalid)
                enableMotionInDisc(motion);
            else if (!motion->parent)
            {
                pnullMotions.push_back(motion);
                motion->valid = true;
                motion->pmotion = pmotion;
                nullMotions.push_back(motion);
                motion->cell->data->root++;
            }
        }
    }
    for (auto & nullm : nullMotions)
    {
        for (auto & child : nullm->children)
        {
            child->valid = false;
            child->parent = nullptr;
            child->pmotion = nullm;
            pnullMotions.push_back(child);
            child->cell->data->root++;
        }
        nullm->children.clear();
    }
    return tvalid;
}
*/

void ompl::geometric::BiASE::processAdEllipsoidRind(bool clearoradd, unsigned int &adinfcount)
{
    if (clearoradd)
    {
        startAdInfPdf_.clear();
        startAdElems_.clear();

        goalAdInfPdf_.clear();
        goalAdElems_.clear();

        calculateInfProb(false);
        ais_ = true;
        adinfcount = 0;
    }
    else if (ais_)
    {
        adinfcount++;
        if (adinfcount == 100) // todo
        {
            clearStartAdInfSampler();
            clearGoalAdInfSampler();
            ais_ = false;
            tree_ = -1;
            startAdInfProb_ = -1.0;
        }
        else if (adinfcount == 50) // todo
        {
            double measure = 0.0;
            for (auto & sampler : startAdInfSamplers_)
            {
                auto adsampler = sampler->as<base::PathLengthDirectAdInfSampler>();
                measure += adsampler->getDirectInformedMeasure() * adsampler->getDirectSamplingFraction();
            }
            for (auto & sampler : goalAdInfSamplers_)
            {
                auto adsampler = sampler->as<base::PathLengthDirectAdInfSampler>();
                measure += adsampler->getDirectInformedMeasure() * adsampler->getDirectSamplingFraction();
            }
            if (measure < 0.5 * si_->getSpaceMeasure())
            {
                if (useInformedSampling_)
                {
                    if (!startAdInfSamplers_.empty())
                    {
                        unsigned int dim = startAdInfSamplers_[0]->as<base::PathLengthDirectAdInfSampler>()->getInformedDimension();
                        factor_ = std::pow(2.0, 1.0 / static_cast<double>(dim));
                    }
                    else 
                    {
                        unsigned int dim = goalAdInfSamplers_[0]->as<base::PathLengthDirectAdInfSampler>()->getInformedDimension();
                        factor_ = std::pow(2.0, 1.0 / static_cast<double>(dim));
                    }
                }
                for (auto & sampler : startAdInfSamplers_)
                    sampler->update(factor_);
                for (auto & sampler : goalAdInfSamplers_)
                    sampler->update(factor_);
//                calculateInfProb(true);
            }
        }
    }
}

void ompl::geometric::BiASE::localInfeasible(int &tree, bool &locals, bool &localg)
{
    locals = localg = false;

    double valid = 0.0, num = 0.0;
    int invalid =  0;
    if (opt_->isFinite(checkedStartPath_.front()->cost))
    {
        tree = 1;
        for (std::size_t i = 0; i < checkedGoalPath_.size() - 1; i++)
        {
            Motion *motion = checkedGoalPath_[i];
            double d = si_->distance(motion->state, checkedGoalPath_[i+1]->state);
            if (motion->valid && motion->parent)
                valid += d;
            else
                invalid++;
            num += d;
        }
        if (invalid <= 2)
            localg = true;
        else
        {
            double ratio = valid/num;
            localg = ratio > localRatio_ ? true : false;
        }
    }
    else if (opt_->isFinite(checkedGoalPath_.front()->cost))
    {
        tree = 0;
        for (std::size_t i = 0; i < checkedStartPath_.size() - 1; i++)
        {
            Motion *motion = checkedStartPath_[i];
            double d = si_->distance(checkedStartPath_[i+1]->state, motion->state);
            if (motion->valid && motion->parent)
                valid += d;
            else 
                invalid++;
            num += d;
        }
        if (invalid <= 2)
            locals = true;
        else
        {
            double ratio = valid/num;
            locals = ratio > localRatio_ ? true : false;
        }
    }
    else 
    {
        tree = -1;
        for (std::size_t i = 0; i < checkedStartPath_.size() - 1; i++)
        {
            Motion *motion = checkedStartPath_[i];
            double d = si_->distance(checkedStartPath_[i+1]->state, motion->state);
            if (motion->valid && motion->parent)
                valid += d;
            else 
                invalid++;
            num += d;
        }
        if (invalid <= 2)
            locals = true;
        else 
        {
            double ratio = valid/num;
            locals = ratio > localRatio_ ? true : false;
        }

        valid = num = 0.0;
        invalid = 0;
        for (std::size_t i = 0; i < checkedGoalPath_.size() - 1; i++)
        {
            Motion *motion = checkedGoalPath_[i];
            double d = si_->distance(motion->state, checkedGoalPath_[i+1]->state);
            if (motion->valid && motion->parent)
                valid += d;
            else 
                invalid++;
            num += d;
        }
        if (invalid <= 2)
            localg = true;
        else 
        {
            double ratio = valid/num;
            localg = ratio > localRatio_ ? true : false;
        }
    }
}

void ompl::geometric::BiASE::calculateInfSampler(bool local, bool start, bool &clearoradd)
{
    if (start)
    {
        if (startAdInfSamplers_.empty() || local)
        {
            factor_ = 1.1;
            clearoradd = true;

            clearStartAdInfSampler();
            startInfSampler(local, startAdInfSamplers_);
        }
    }
    else 
    {
        if (goalAdInfSamplers_.empty() || local)
        {
            factor_ = 1.1;
            clearoradd = true;

            clearGoalAdInfSampler();
            goalInfSampler(local, goalAdInfSamplers_);
        }
    }
}

void ompl::geometric::BiASE::startInfSampler(bool local, std::vector<base::AdInformedSamplerPtr> &infSamplers)
{
    if (local)
        startLocalInfSampler(infSamplers);
    else 
    {
        Motion *motion1 = nullptr, *motion2 = nullptr;
        base::Cost cost = opt_->identityCost();
        base::Cost cost1= opt_->identityCost();
        motion1 = checkedStartPath_.back();
        if (false)//(checkedGoalPath_.size() > 8) // todo
        {
            motion2 = checkedGoalPath_.front();
            for (std::size_t i = 0; i < checkedGoalPath_.size() - 1 && i < 2; i++)
            {
                cost1 = opt_->combineCosts(cost1, opt_->motionCost(checkedGoalPath_[i]->state, checkedGoalPath_[i+1]->state));
                motion2 = checkedGoalPath_[i+1];
            }
        }
        else 
            motion2 = checkedStartPath_.front();

        cost = opt_->combineCosts(cost, cost1);
        for (std::size_t i = 0; i < checkedStartPath_.size() - 1; i++)
            cost = opt_->combineCosts(cost, opt_->motionCost(checkedStartPath_[i+1]->state, checkedStartPath_[i]->state));

        while (cost.value() <= si_->distance(motion1->state, motion2->state))
            cost = base::Cost(factor_ * cost.value());

        base::AdInformedSamplerPtr sampler = allocInfSampler(motion1->state, motion2->state, base::Cost(0.0), cost);
        double ratio = (double)tStart_->size() / (double)(tStart_->size() + tGoal_->size());
        if (useInformedSampling_ && sampler->getInformedMeasure() >= ratio * si_->getSpaceMeasure())
        {
            std::size_t localseg = std::ceil(0.2 * (double)(checkedStartPath_.size() - 1));
            localseg = localseg > 3 ? localseg : 3; // todo

            std::size_t i = 0;
            while (i < checkedStartPath_.size() - 1)
            {
                cost = opt_->identityCost();
                if (i == 0)
                    cost = opt_->combineCosts(cost, cost1);
                else 
                    motion2 = checkedStartPath_[i];

                if (i + localseg >= checkedStartPath_.size() - 3)
                    localseg = checkedStartPath_.size() - 1 - i;

                std::size_t j = i;
                while (j < i + localseg && j < checkedStartPath_.size() - 1)
                {
                    cost = opt_->combineCosts(cost, opt_->motionCost(checkedStartPath_[j+1]->state, checkedStartPath_[j]->state));
                    j++;
                }

                motion1 = checkedStartPath_[j];

                while (cost.value() <= si_->distance(motion1->state, motion2->state))
                    cost = base::Cost(factor_ * cost.value());

                sampler = allocInfSampler(motion1->state, motion2->state, base::Cost(0.0), cost);
                infSamplers.push_back(sampler);

                i = j;
            }
        }
        else 
            infSamplers.push_back(sampler);
    }
}

void ompl::geometric::BiASE::startLocalInfSampler(std::vector<base::AdInformedSamplerPtr> &infSamplers)
{
    Motion *motion1 = nullptr, *motion2 = nullptr;
    std::size_t localseg = std::ceil(0.1 * (double)(checkedStartPath_.size() - 1)); // todo
//    std::size_t localseg = 1;
    std::size_t j = 0;
    while (j < checkedStartPath_.size() - 1)
    {
        std::size_t i = j;
        bool feas = true;
        while (i < checkedStartPath_.size() - 1)
        {
            if (!checkedStartPath_[i]->valid || !checkedStartPath_[i]->parent)
            {
                feas = false;
                break;
            }
            i++;
        }
        if (feas)
            break;

        base::Cost cost = opt_->identityCost();

        std::size_t low = i > localseg ? (i - localseg) : 0;
        std::size_t high= std::min(i + localseg + 1, checkedStartPath_.size() - 1);

        if (false)//(j == 0 && low <= 2) // todo
        {
            if (false)//(checkedGoalPath_.size() > 8) // todo
            {
                motion2 = checkedGoalPath_.front();
                for (std::size_t i = 0; i < checkedGoalPath_.size() - 1 && i < 2; i++)
                {
                    cost = opt_->combineCosts(cost, opt_->motionCost(checkedGoalPath_[i]->state, checkedGoalPath_[i+1]->state));
                    motion2 = checkedGoalPath_[i+1];
                }
            }
            else 
                motion2 = checkedStartPath_.front();
        }
        else 
            motion2 = checkedStartPath_[low];

        while (high < checkedStartPath_.size() - 1)
        {
            if (checkedStartPath_[high]->valid && checkedStartPath_[high]->parent) 
                break;
            high++;
        }

        if (high == checkedStartPath_.size() - 2)
            high = checkedStartPath_.size() - 1;

        for (std::size_t i = low; i < high ; i++)
            cost = opt_->combineCosts(cost, opt_->motionCost(checkedStartPath_[i+1]->state, checkedStartPath_[i]->state));

        motion1 = checkedStartPath_[high];

        while (cost.value() <= si_->distance(motion1->state, motion2->state))
            cost = base::Cost(factor_ * cost.value());

        base::AdInformedSamplerPtr sampler = allocInfSampler(motion1->state, motion2->state, base::Cost(0.0), cost);
        infSamplers.push_back(sampler);
        j = high;
    }
}

void ompl::geometric::BiASE::goalInfSampler(bool local, std::vector<base::AdInformedSamplerPtr> &infSamplers)
{
    if (local)
        goalLocalInfSampler(infSamplers);
    else 
    {
        Motion *motion1 = nullptr, *motion2 = nullptr;
        base::Cost cost = opt_->identityCost();
        base::Cost cost1 = opt_->identityCost();

        if (false)//(checkedStartPath_.size() > 8) // todo
        {
            motion1 = checkedStartPath_.front();
            for (std::size_t i = 0; i < checkedStartPath_.size() - 1 && i < 2; i++)
            {
                cost1 = opt_->combineCosts(cost1, opt_->motionCost(checkedStartPath_[i+1]->state, checkedStartPath_[i]->state));
                motion1 = checkedStartPath_[i+1];
            }
        }
        else 
            motion1 = checkedGoalPath_.front();
        motion2 = checkedGoalPath_.back();

        cost = opt_->combineCosts(cost, cost1);
        for (std::size_t i = 0; i < checkedGoalPath_.size() - 1; i++)
            cost = opt_->combineCosts(cost, opt_->motionCost(checkedGoalPath_[i]->state, checkedGoalPath_[i+1]->state));

        while (cost.value() <= si_->distance(motion1->state, motion2->state))
            cost = base::Cost(factor_ * cost.value());

        base::AdInformedSamplerPtr sampler = allocInfSampler(motion1->state, motion2->state, base::Cost(0.0), cost);
        double ratio = (double)tGoal_->size() / (double)(tStart_->size() + tGoal_->size());
        if (useInformedSampling_ && sampler->getInformedMeasure() >= ratio * si_->getSpaceMeasure())
        {
            std::size_t localseg = std::ceil(0.2 * (double)(checkedGoalPath_.size() - 1));
            localseg = localseg > 3 ? localseg : 3; // todo
            
            std::size_t i = 0;
            while (i < checkedGoalPath_.size() - 1)
            {
                cost = opt_->identityCost();
                if (i == 0)
                    cost = opt_->combineCosts(cost, cost1);
                else 
                    motion1 = checkedGoalPath_[i];

                if (i + localseg >= checkedGoalPath_.size() - 3)
                    localseg = checkedGoalPath_.size() - 1 - i;

                std::size_t j = i;
                while (j < i + localseg && j < checkedGoalPath_.size() - 1)
                {
                    cost = opt_->combineCosts(cost, opt_->motionCost(checkedGoalPath_[j]->state, checkedGoalPath_[j+1]->state));
                    j++;
                }
                motion2 = checkedGoalPath_[j];

                while (cost.value() <= si_->distance(motion1->state, motion2->state))
                    cost = base::Cost(factor_ * cost.value());

                sampler = allocInfSampler(motion1->state, motion2->state, base::Cost(0.0), cost);
                infSamplers.push_back(sampler);
                i = j;
            }
        }
        else 
            infSamplers.push_back(sampler);
    }
}

void ompl::geometric::BiASE::goalLocalInfSampler(std::vector<base::AdInformedSamplerPtr> &infSamplers)
{
    Motion *motion1 = nullptr, *motion2 = nullptr;
    std::size_t localseg = std::ceil(0.1 * (double)(checkedGoalPath_.size() - 1)); // todo
//    std::size_t localseg = 1;
    std::size_t j = 0;
    while (j < checkedGoalPath_.size() - 1)
    {
        std::size_t i = j;
        bool feas = true;
        while (i < checkedGoalPath_.size() - 1)
        {
            if (!checkedGoalPath_[i]->valid || !checkedGoalPath_[i]->parent)
            {
                feas = false;
                break;
            }
            i++;
        }
        if (feas)
            break;

        std::size_t low = i > localseg ? (i - localseg) : 0;
        std::size_t high= std::min(i + localseg + 1, checkedGoalPath_.size() - 1);

        base::Cost cost = opt_->identityCost();
        if (false)//(j == 0 && low <= 2) // todo
        {
            if (false)//(checkedStartPath_.size() > 8) // todo
            {
                motion1 = checkedStartPath_.front();
                for (std::size_t i = 0; i < checkedStartPath_.size() - 1 && i < 2; i++)
                {
                    cost = opt_->combineCosts(cost, opt_->motionCost(checkedStartPath_[i+1]->state, checkedStartPath_[i]->state));
                    motion1 = checkedStartPath_[i+1];
                }
            }
            else 
                motion1 = checkedGoalPath_.front();
        }
        else 
            motion1 = checkedGoalPath_[low];

        while (high < checkedGoalPath_.size() - 1)
        {
            if (checkedGoalPath_[high]->valid && checkedGoalPath_[high]->parent)
                break;
            high++;
        }
        if (high == checkedGoalPath_.size() - 2)
            high = checkedGoalPath_.size() - 1;

        motion2 = checkedGoalPath_[high];
        for (std::size_t i = low; i < high ; i++)
            cost = opt_->combineCosts(cost, opt_->motionCost(checkedGoalPath_[i]->state, checkedGoalPath_[i+1]->state));
        while (cost.value() <= si_->distance(motion1->state, motion2->state))
            cost = base::Cost(factor_ * cost.value());

        base::AdInformedSamplerPtr sampler = allocInfSampler(motion1->state, motion2->state, base::Cost(0.0), cost);
        infSamplers.push_back(sampler);
        j = high;
    }
}

void ompl::geometric::BiASE::calculateInfProb(bool update)
{
    calculateInfProb(startAdInfSamplers_, goalAdInfSamplers_, startAdInfPdf_, goalAdInfPdf_, startAdElems_, goalAdElems_, update);
}

void ompl::geometric::BiASE::calculateInfProb(const std::vector<base::AdInformedSamplerPtr> &startInfSamplers, 
                                              const std::vector<base::AdInformedSamplerPtr> &goalInfSamplers,
                                              NumPdf &startInfPdf, NumPdf &goalInfPdf, 
                                              std::vector<NumElem *> &startelems, std::vector<NumElem *> &goalelems, bool update)
{
    double startm = 0.0, goalm = 0.0;
    if (!startInfSamplers.empty())
    {
        double measure = 0;
        for (auto & sampler : startInfSamplers)
            measure += sampler->getInformedMeasure();
        startm = measure;
        if (startInfSamplers.size() > 1)
        {
            if (!update)
            {
                for (std::size_t i = 0; i < startInfSamplers.size(); i++)
                    startelems.push_back(startInfPdf.add(i, startInfSamplers[i]->getInformedMeasure() / measure));
            }
            else 
            {
                for (std::size_t i = 0; i < startInfSamplers.size(); i++)
                    startInfPdf.update(startelems[i], startInfSamplers[i]->getInformedMeasure() / measure);
            }
        }
    }
    if (!goalInfSamplers.empty())
    {
        double measure = 0;
        for (auto & sampler : goalInfSamplers)
            measure += sampler->getInformedMeasure();
        goalm = measure;
        if (goalInfSamplers.size() > 1)
        {
            if (!update)
            {
                for (std::size_t i = 0; i < goalInfSamplers.size(); i++)
                    goalelems.push_back(goalInfPdf.add(i, goalInfSamplers[i]->getInformedMeasure() / measure));
            }
            else 
            {
                for (std::size_t i = 0; i < goalInfSamplers.size(); i++)
                    goalInfPdf.update(goalelems[i], goalInfSamplers[i]->getInformedMeasure() / measure);
            }
        }
    }

    if (startm != 0.0 && goalm != 0.0)
        startAdInfProb_ = startm / (startm + goalm);
}

ompl::base::AdInformedSamplerPtr ompl::geometric::BiASE::allocInfSampler(const base::State *s1, const base::State *s2,
                                                                       const base::Cost &minCost, const base::Cost &maxCost)
{
    if (useRejectionSampling_)
        return std::make_shared<base::RejectionAdInfSampler>(pdef_, s1, s2, minCost, maxCost, numSampleAttempts_);
    else 
        return std::make_shared<base::PathLengthDirectAdInfSampler>(pdef_, s1, s2, minCost, maxCost, numSampleAttempts_);
}

bool ompl::geometric::BiASE::sampleUniformAd(base::State *state, bool start)
{
    if (start)
        return sampleUniform(state, startAdInfSamplers_, startAdInfPdf_);
    return sampleUniform(state, goalAdInfSamplers_, goalAdInfPdf_);
}

bool ompl::geometric::BiASE::sampleUniform(base::State *state, const std::vector<base::AdInformedSamplerPtr> &infSamplers, const NumPdf &adInfPdf)
{
    if (!infSamplers.empty())
    {
        if (infSamplers.size() == 1)
            return infSamplers[0]->sampleUniform(state);
        else
        {
            std::size_t ii = adInfPdf.sample(rng_.uniform01());
            bool valid = infSamplers[ii]->sampleUniform(state);
            /*
            if (valid && useInformedSampling_) // todo
            {
                double reject = 0;
                for (std::size_t i = 0; i < infSamplers.size(); i++)
                {
                    if (i != ii && infSamplers[i]->as<base::PathLengthDirectAdInfSampler>()->isInPhs(state)) 
                        reject += 1.0;
                }
                if (reject != 0 && rng_.uniform01() < reject/(1.0 + reject))
                    return false;
            }
            */
            return valid;
        }
    }
    else
    {
        sampler_->sampleUniform(state);
        return true;
    }
}
