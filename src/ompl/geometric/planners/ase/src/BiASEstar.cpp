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

#include "ompl/geometric/planners/ase/BiASEstar.h"
#include "ompl/base/goals/GoalSampleableRegion.h"
#include "ompl/base/objectives/PathLengthOptimizationObjective.h"
#include "ompl/tools/config/SelfConfig.h"

#include "ompl/base/samplers/adinformed/RejectionAdInfSampler.h"
#include "ompl/base/samplers/adinformed/PathLengthDirectAdInfSampler.h"

ompl::geometric::BiASEstar::BiASEstar(const base::SpaceInformationPtr &si) : base::Planner(si, "BiASEstar")
  , dStart_(0), dGoal_(0), mc_(opt_), bh_(mc_)
{
    specs_.recognizedGoal = base::GOAL_SAMPLEABLE_REGION;
    specs_.directed = true;
    specs_.approximateSolutions = false;
    specs_.optimizingPaths = true;
    specs_.canReportIntermediateSolutions = true;

    setInformedSampling(true);
    setSampleRejection(false);

    Planner::declareParam<double>("range", this, &BiASEstar::setRange, &BiASEstar::getRange, "0.:0.1:10000.");
    Planner::declareParam<double>("pen_distance", this, &BiASEstar::setPenDistance, &BiASEstar::getPenDistance, "0.:0.1:10000.");
//    Planner::declareParam<bool>("use_bispace", this, &BiASEstar::setUseBispace, &BiASEstar::getUseBispace, "0,1");
    Planner::declareParam<double>("prune_threshold", this, &BiASEstar::setPruneThreshold, &BiASEstar::getPruneThreshold, "0.:.01:1.");
    Planner::declareParam<bool>("update_nbhcell", this, &BiASEstar::setUpdateNbCell, &BiASEstar::getUpdateNbCell, "0,1");

    addPlannerProgressProperty("iterations INTEGER", [this] { return numIterationsProperty(); });
    addPlannerProgressProperty("best cost REAL", [this] { return bestCostProperty(); });

    Planner::declareParam<bool>("informed_sampling", this, &BiASEstar::setInformedSampling, &BiASEstar::getInformedSampling, "0,1");
    Planner::declareParam<bool>("sample_rejection", this, &BiASEstar::setSampleRejection, &BiASEstar::getSampleRejection, "0,1");
    Planner::declareParam<unsigned int>("number_sampling_attempts", this, &BiASEstar::setNumSamplingAttempts,
                                        &BiASEstar::getNumSamplingAttempts, "10:10:100000");
}

ompl::geometric::BiASEstar::~BiASEstar()
{
    freeMemory();
}

void ompl::geometric::BiASEstar::setup()
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
    dStart_.setMaxNeighborCell(4);
    dGoal_.setMaxNeighborCell(4);

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
        currentStartCost_ = opt_->infiniteCost();
        currentGoalCost_ = opt_->infiniteCost();
        bestCost_ = opt_->infiniteCost();
        prunedCost_ = opt_->infiniteCost();
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

void ompl::geometric::BiASEstar::freeMemory()
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

void ompl::geometric::BiASEstar::clear()
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
    bestStartMotion_ = nullptr;
    bestGoalMotion_ = nullptr;

    bestCost_ = base::Cost(std::numeric_limits<double>::quiet_NaN());

    solved_ = false;

    iterations_ = 0;
    prunedCost_ = base::Cost(std::numeric_limits<double>::quiet_NaN());
    currentStartCost_ = base::Cost(std::numeric_limits<double>::quiet_NaN());
    currentGoalCost_ = base::Cost(std::numeric_limits<double>::quiet_NaN());

    clearStartAdInfSampler();
    clearGoalAdInfSampler();
    tree_ = -1;
    localRatio_ = 0.75;

    clearStartInfSampler();
    clearGoalInfSampler();
    guniform_ = true;

    adinfcount_ = infcount_ = 0;
}

ompl::base::PlannerStatus ompl::geometric::BiASEstar::solve(const base::PlannerTerminationCondition &ptc)
{
    base::PlannerStatus ps= prepareSolve(ptc);
    if (ps != base::PlannerStatus::PREPARE_SUCCESS)
        return ps;

    clearStartAdInfSampler();
    clearGoalAdInfSampler();
    clearStartInfSampler();
    clearGoalInfSampler();
    sampler_ = si_->allocStateSampler();

    OMPL_INFORM("%s: Starting planning with %u states already in datastructure. Seeking a solution better than %.5f.",
                getName().c_str(), (tStart_->size() + tGoal_->size()), opt_->getCostThreshold().value());

    const base::ReportIntermediateSolutionFn intermediateSolutionCallback = pdef_->getIntermediateSolutionCallback();

    bool startTree = true;
    bool optimal = false;

    Motion *startAd = nullptr, *goalAd = nullptr;
    bool ais = false;
    bool adinf = true;
    bool reverse = false;

    while (!ptc)
    {
        iterations_++;
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
        
        if (!batchGrow(startTree, ais, adinf))
            continue;

        if (!solved_ && !startAd)
        {
            std::size_t index = 0;
            if (selectCMotion(index, reverse))
            {
                auto pair = connectionPoint_[index];
                startAd = pair.first;
                goalAd = pair.second;
            }
        }

        bool ad = false, clearoradd = false;
        bool updatedSolution = findBetterSolution(startAd, goalAd, ad, ais, clearoradd, optimal);
        if (optimal)
            break;
        if (ad)
        {
            clearoradd = false;

            ais = false;
            startAd = nullptr;
            goalAd = nullptr;

            tree_ = -1;
            adinfcount_ = 0;
            clearStartAdInfSampler();
            clearGoalAdInfSampler();
        }
            
        if (updatedSolution)
        {
            if (startAd && !(keepCondition(startAd, bestCost_, true) && keepCondition(goalAd, bestCost_, false)))
            {
                clearoradd = false;

                ais = false;
                startAd = nullptr;
                goalAd = nullptr;

                tree_ = -1;
                adinfcount_ = 0;
                clearStartAdInfSampler();
                clearGoalAdInfSampler();
            }
            int numPruned = pruneTree(bestCost_);
            if (0)
                OMPL_INFORM("%s: %u states are pruned from the tree, %u states are left", getName().c_str(), numPruned, tStart_->size() + tGoal_->size());
            if (opt_->isFinite(bestStartMotion_->cost) && opt_->isFinite(bestGoalMotion_->cost))
            {
                reportBetterSolution(intermediateSolutionCallback);
            }
        }
        else if (solved_)
        {
            if (infcount_ >= 50u) // todo
            {
                infcount_ = 0;
                clearStartInfSampler();
                clearGoalInfSampler();
            }
            /*
            else if (infcount_ == 50u) // todo
            {
                double measure = 0.0;
                for (auto & sampler : startInfSamplers_)
                {
                    auto adsampler = sampler->as<base::PathLengthDirectAdInfSampler>();
                    measure += adsampler->getDirectInformedMeasure() * adsampler->getDirectSamplingFraction();
                }
                for (auto & sampler : goalInfSamplers_)
                {
                    auto adsampler = sampler->as<base::PathLengthDirectAdInfSampler>();
                    measure += adsampler->getDirectInformedMeasure() * adsampler->getDirectSamplingFraction();
                }
                if (measure < 0.5 * si_->getSpaceMeasure())
                {
                    if (useInformedSampling_)
                    {
                        if (!startInfSamplers_.empty())
                        {
                            unsigned int dim = startInfSamplers_[0]->as<base::PathLengthDirectAdInfSampler>()->getInformedDimension();
                            factor_ = std::pow(2.0, 1.0 / static_cast<double>(dim));
                        }
                        else 
                        {
                            unsigned int dim = goalInfSamplers_[0]->as<base::PathLengthDirectAdInfSampler>()->getInformedDimension();
                            factor_ = std::pow(2.0, 1.0 / static_cast<double>(dim));
                        }
                    }
                    for (auto & sampler : startInfSamplers_)
                        sampler->update(factor_);
                    for (auto & sampler : goalInfSamplers_)
                        sampler->update(factor_);
                }
            }
            */
        }

        processAdEllipsoidRind(clearoradd, ais);
        if (!ais)
            startAd = goalAd = nullptr;
    }

    if (solved_ || optimal)
    {
        ptc.terminate();
        if (!optimal)
            isPathValid(bestStartMotion_, bestGoalMotion_);
        processSolution(bestStartMotion_, bestGoalMotion_);
    }

    OMPL_INFORM("%s: Created %u states (%u start + %u goal). Final solution cost %.5f", getName().c_str(), tStart_->size() + tGoal_->size(),
                tStart_->size(), tGoal_->size(), bestCost_.value());

    return (solved_ || optimal) ? base::PlannerStatus::EXACT_SOLUTION : base::PlannerStatus::TIMEOUT;
}

// grow
ompl::base::PlannerStatus ompl::geometric::BiASEstar::prepareSolve(const base::PlannerTerminationCondition &ptc)
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

void ompl::geometric::BiASEstar::processSolution(const Motion *bestStartMotion, const Motion *bestGoalMotion)
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
    psol.setOptimized(opt_, bestCost_, opt_->isSatisfied(bestCost_));
    pdef_->addSolutionPath(psol);
}

bool ompl::geometric::BiASEstar::growTree(TreeGrowingInfo &tgi, Motion *rmotion, bool &otherSide, bool &change, bool &add)
{
    if (treatedAsMultiSubapce_)
        return growTreeMultiSpace(tgi, rmotion, otherSide, change, add);
    return growTreeSingleSpace(tgi, rmotion, otherSide, change, add);
}

bool ompl::geometric::BiASEstar::growTreeSingleSpace(TreeGrowingInfo &tgi, Motion *rmotion, bool &otherSide, bool &change, bool &add)
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

        if (!isValid(dstate))
        {
            reach = false;
            break;
        }

        Motion *motion = new Motion(si_);
        si_->copyState(motion->state, dstate);
        Coord xcoord(projectionEvaluator_->getDimension());
        projectionEvaluator_->computeCoordinates(motion->state, xcoord);
        addToDisc(disc, motion, xcoord);

        Motion *nb = nmotion;
        connectToPmotion(motion, nb, tgi.start);
        CostMotionCompare compareFn(motion, opt_, tgi.start);
        for (auto & cv : motion->cell->nbh)
        {
            for (auto & c : cv)
            {
                if (c->data->motions.size() > c->data->disabled)
                {
                    Motion *temp = *std::min_element(c->data->motions.begin(), c->data->motions.end(), compareFn);
                    if (temp != nb)
                    {
                        base::Cost incCost = tgi.start ? opt_->motionCost(temp->state, motion->state) : opt_->motionCost(motion->state, temp->state);
                        base::Cost cost = opt_->combineCosts(temp->cost, incCost);
                        if (opt_->isCostBetterThan(cost, motion->cost))
                        {
                            nb = temp;
                            motion->parent = temp;
                            motion->incCost = incCost;
                            motion->cost = cost;
                            motion->root = temp->root;
                        }
                    }
                }
            }
        }
        if (solved_ && guniform_ && !keepCondition2(motion, bestCost_, tgi.start))
        {
            removeFromDisc(disc, motion);
            freeMotion(motion);
            reach = false;
            break;
        }
        motion->parent->children.push_back(motion);
        tree->add(motion);
        tgi.xmotion = motion;
        add = true;
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
        if (rewireSort_)
        {
            if (!motion->cell->data->cmotion || opt_->isCostBetterThan(motion->cost, motion->cell->data->cmotion->cost))
                motion->cell->data->cmotion = motion;
            if (!motion->cell->data->mmotion || opt_->isCostBetterThan(motion->cell->data->mmotion->cost, motion->cost))
                motion->cell->data->mmotion = motion;
        }
        if (!solved_)
        {
            if (tgi.start ? !pnullStartMotions_.empty() : !pnullGoalMotions_.empty())
                rewireTree(bh_, motion, tgi.start);
        }
        else 
            optimalRewireTree(bh_, motion, tgi.start);
    }
    return reach;
}

bool ompl::geometric::BiASEstar::growTreeMultiSpace(TreeGrowingInfo &tgi, Motion *rmotion, bool &otherSide, bool &change, bool &add)
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
        Coord xcoord(projectionEvaluator_->getDimension());
        projectionEvaluator_->computeCoordinates(motion->state, xcoord);
        addToDisc(disc, motion, xcoord);

        Motion *nb = nmotion;
        connectToPmotion(motion, nb, tgi.start);
        CostMotionCompare compareFn(motion, opt_, tgi.start);
        for (auto & cv : motion->cell->nbh)
        {
            for (auto & c : cv)
            {
                if (c->data->motions.size() > c->data->disabled)
                {
                    Motion *temp = *std::min_element(c->data->motions.begin(), c->data->motions.end(), compareFn);
                    if (temp != nb)
                    {
                        base::Cost incCost = tgi.start ? opt_->motionCost(temp->state, motion->state) : opt_->motionCost(motion->state, temp->state);
                        base::Cost cost = opt_->combineCosts(temp->cost, incCost);
                        if (opt_->isCostBetterThan(cost, motion->cost))
                        {
                            nb = temp;
                            motion->parent = temp;
                            motion->incCost = incCost;
                            motion->cost = cost;
                            motion->root = temp->root;
                        }
                    }
                }
            }
        }

        if (solved_ && guniform_ && !keepCondition2(motion, bestCost_, tgi.start))
        {
            removeFromDisc(disc, motion);
            freeMotion(motion);
            reach = false;
            break;
        }

        motion->parent->children.push_back(motion);
        tree->add(motion);
        tgi.xmotion = motion;
        add = true;
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

        if (rewireSort_)
        {
            if (!motion->cell->data->cmotion || opt_->isCostBetterThan(motion->cost, motion->cell->data->cmotion->cost))
                motion->cell->data->cmotion = motion;
            if (!motion->cell->data->mmotion || opt_->isCostBetterThan(motion->cell->data->mmotion->cost, motion->cost))
                motion->cell->data->mmotion = motion;
        }
        if (!solved_)
        {
            if (tgi.start ? !pnullStartMotions_.empty() : !pnullGoalMotions_.empty())
                rewireTree(bh_, motion, tgi.start);
        }
        else 
            optimalRewireTree(bh_, motion, tgi.start);
    }
    return reach && reachi;
}

ompl::geometric::BiASEstar::Motion *ompl::geometric::BiASEstar::selectNMotion(const TreeData &tree, Motion *rmotion, bool &null)
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

ompl::geometric::BiASEstar::Motion *ompl::geometric::BiASEstar::selectMotionInCell(Cell *cell)
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

void ompl::geometric::BiASEstar::rewireTree(BinaryHeap<Motion *, MotionCompare> &bh, Motion *m, bool start)
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
                if (rewireSort_ && motion->cell != c && c->data->mmotion)
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

void ompl::geometric::BiASEstar::updateQueue(BinaryHeap<Motion *, MotionCompare> &bh, Motion *motion)
{
    if (motion->handle != nullptr)
        bh.update(motion->handle);
    else
        motion->handle = bh.insert(motion);
}

bool ompl::geometric::BiASEstar::growCurrentTree(const base::State *state) const
{
    return growStartTree(state);
}

bool ompl::geometric::BiASEstar::growStartTree(const base::State *state) const
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

double ompl::geometric::BiASEstar::penetrationDistance(const base::State *nstate, const base::State *state, bool start) const
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

bool ompl::geometric::BiASEstar::growCurrentTree(const base::State *state, std::size_t sub) const
{
    return growStartTree(state, sub);
}

bool ompl::geometric::BiASEstar::growStartTree(const base::State *state, std::size_t sub) const
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

double ompl::geometric::BiASEstar::penetrationDistance(const base::State *nstate, const base::State *state, bool start, std::size_t sub) const
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

bool ompl::geometric::BiASEstar::backPathRewireMotion(Motion *motion, bool start, Motion *&pmotion)
{
    bool valid = false;
    OrderCellsByCost ocbc(opt_);
    CostMotionCompare compareFn(motion, opt_, start);
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
            if (rewireSort_ && motion->cell != c && motion->cell->data->mmotion)
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
            {
                std::sort(c->data->motions.begin(), c->data->motions.end(), compareFn);
                nbh.insert(nbh.end(), c->data->motions.begin(), c->data->motions.end());
            }
            else 
            {
                for (std::size_t i = 0; i < c->data->motions.size() && nbh.size() < num; i++)
                {
                    if (opt_->isFinite(c->data->motions[i]->cost))
                        nbh.push_back(c->data->motions[i]);
                }
                std::sort(nbh.begin(), nbh.end(), compareFn);
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

void ompl::geometric::BiASEstar::addPdfMotion(MotionPDF &pdf, Motion *motion, bool start)
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

ompl::geometric::BiASEstar::Motion *ompl::geometric::BiASEstar::selectPdfMotion(MotionPDF &pdf, GridCell *&cell)
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

void ompl::geometric::BiASEstar::removePdfMotion(MotionPDF &pdf, Motion *motion)
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

void ompl::geometric::BiASEstar::enableMotionInDisc(Motion *motion)
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

void ompl::geometric::BiASEstar::connectToPmotion(Motion *motion, Motion *pmotion, bool start) const
{
    motion->parent = pmotion;
    motion->root = pmotion->root;
    motion->incCost = start ? opt_->motionCost(pmotion->state, motion->state) : opt_->motionCost(motion->state, pmotion->state);
    motion->cost = opt_->combineCosts(pmotion->cost, motion->incCost);
    if (opt_->isFinite(motion->cost))
        updateChildCosts(motion);
}

void ompl::geometric::BiASEstar::updateChildCosts(Motion *motion) const
{
    for (auto & child : motion->children)
    {
        child->cost = opt_->combineCosts(motion->cost, child->incCost);
        child->root = motion->root;
        updateChildCosts(child);
    }
}

void ompl::geometric::BiASEstar::insertNeighbor(Motion *pmotion, Motion *motion)
{
    pmotion->nbh.insert(motion);
    motion->nbh.insert(pmotion);
}

void ompl::geometric::BiASEstar::insertInvalidNeighbor(Motion *pmotion, Motion *motion)
{
    pmotion->invalidnbh.insert(motion);
    motion->invalidnbh.insert(pmotion);
}

void ompl::geometric::BiASEstar::removeFromNeighbor(Motion *motion)
{
    for (auto & pmotion : motion->nbh)
        pmotion->nbh.erase(motion);
}

void ompl::geometric::BiASEstar::removeFromInvalidNeighbor(Motion *motion)
{
    for (auto & pmotion : motion->invalidnbh)
        pmotion->invalidnbh.erase(motion);
}

bool ompl::geometric::BiASEstar::isValidNeighbor(Motion *motion, Motion *nb) const
{
    return motion->nbh.find(nb) != motion->nbh.end();
}

bool ompl::geometric::BiASEstar::isInvalidNeighbor(Motion *motion, Motion *nb) const
{
    return motion->invalidnbh.find(nb) != motion->invalidnbh.end();
}

void ompl::geometric::BiASEstar::addToDisc(CellDiscretizationData &disc, Motion *motion, const Coord& coord)
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

void ompl::geometric::BiASEstar::removeFromDisc(CellDiscretizationData &disc, Motion *motion)
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

void ompl::geometric::BiASEstar::removeFromParent(Motion *motion)
{
    if (!motion->parent)
        return;
    removeFromVector(motion->parent->children, motion);
}

void ompl::geometric::BiASEstar::removeFromPnull(std::vector<Motion *> &pnullMotions, Motion *motion)
{
    if (removeFromVector(pnullMotions, motion))
    {
        removeFromVector(motion->pmotion->pchildren, motion);
        motion->pmotion = nullptr;
    }
}

bool ompl::geometric::BiASEstar::removeFromVector(std::vector<Motion *> &motions, Motion *motion)
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

void ompl::geometric::BiASEstar::setMotionInfinityCostWithDisable(Motion *motion, std::unordered_set<Cell *> &cells) const
{
    motion->cost = opt_->infiniteCost();
    motion->cell->data->disabled++;
    if (rewireSort_)
        cells.insert(motion->cell);
    for (auto & child : motion->children)
        setMotionInfinityCostWithDisable(child, cells);
}

// check 
bool ompl::geometric::BiASEstar::isValid(const base::State *state)
{
    bool valid = si_->isValid(state);
    return valid;
}

bool ompl::geometric::BiASEstar::checkMotion(Motion *pmotion, Motion *motion, bool start)
{
    if (!motion->valid)
    {
        if (!checkInterMotion(pmotion, motion, start))
        {
            if (opt_->isFinite(motion->cost))
            {
                std::unordered_set<Cell *> cells;
                setMotionInfinityCostWithDisable(motion, cells);
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

bool ompl::geometric::BiASEstar::checkInterMotion(Motion *pmotion, Motion *motion, bool start)
{
    bool valid = false;
    if (start || symmetric_)
        valid = checkInterMotion1(pmotion, motion, start);
    else 
        valid = checkInterMotion2(motion, pmotion, start);
    return valid;
}

bool ompl::geometric::BiASEstar::checkInterMotion1(Motion *smotion, Motion *gmotion, bool start)
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

bool ompl::geometric::BiASEstar::checkInterMotion2(Motion *smotion, Motion *gmotion, bool start)
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

void ompl::geometric::BiASEstar::addIntermediateMotion(Motion *pmotion, Motion *motion, bool start, Motion *last)
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
    if (!last->cell->data->cmotion || opt_->isCostBetterThan(last->cost, last->cell->data->cmotion->cost))
        last->cell->data->cmotion = last;
    if (!last->cell->data->mmotion || opt_->isCostBetterThan(last->cell->data->mmotion->cost, last->cost))
        last->cell->data->mmotion = last;
}

void ompl::geometric::BiASEstar::getPlannerData(base::PlannerData &data) const
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

// optimal
bool ompl::geometric::BiASEstar::findBetterSolution(Motion *startAd, Motion *goalAd, bool &ad, bool ais, bool &clearoradd, bool &optimal)
{
    bool updatedSolution = false;
    if (startAd)
    {
        base::Cost temp = opt_->combineCosts(startAd->cost, goalAd->cost);
        if (opt_->isCostBetterThan(temp, bestCost_))
        {
            if (isPathValid(startAd, goalAd))
            {
                ad = true;
                temp = opt_->combineCosts(startAd->cost, goalAd->cost);
                if (opt_->isCostBetterThan(temp, bestCost_))
                {
                    bestCost_ = temp;

                    if (solved_)
                        OMPL_INFORM("%s: Found a better solution with a cost of %.2f in %u iterations (%u "
                                "vertices in the graph, %u start + %u goal)",
                                getName().c_str(), bestCost_.value(), iterations_, tStart_->size() + tGoal_->size(), tStart_->size(), tGoal_->size());
                    else 
                    {
                        OMPL_INFORM("%s: Found an initial solution with a cost of %.2f in %u iterations (%u "
                                "vertices in the graph)",
                                getName().c_str(), bestCost_.value(), iterations_, tStart_->size() + tGoal_->size());
                        for (auto & motion : startMotions_)
                            optimalRewireTree(bh_, motion, true);
                        for (auto & motion : goalMotions_)
                            optimalRewireTree(bh_, motion, false);
                        checkedStartPath_.clear();
                        checkedGoalPath_.clear();
                    }

                    solved_ = true;
                    updatedSolution = true;
                    bestStartMotion_ = startAd;
                    bestGoalMotion_ = goalAd;

                    if (opt_->isSatisfied(bestCost_))
                        optimal = true;
                    else 
                    {
                        infcount_ = 0;
                        clearStartInfSampler();
                        clearGoalInfSampler();
                        Motion *rtSt = nullptr, *rtG = nullptr;
                        for (auto & stMotion : startMotions_)
                        {
                            if (stMotion->root == bestStartMotion_->root)
                            {
                                rtSt = stMotion;
                                break;
                            }
                        }
                        for (auto & gMotion : goalMotions_)
                        {
                            if (gMotion->root == bestGoalMotion_->root)
                            {
                                rtG = gMotion;
                                break;
                            }
                        }
                        optimalInfSampler(rtSt, bestStartMotion_, true, startInfSamplers_);
                        optimalInfSampler(rtG,  bestGoalMotion_, false, goalInfSamplers_);
                        calculateOptimalInfProb();
                    }
                }
            }
            else if (!ais && !checkedStartPath_.empty() && !checkedGoalPath_.empty())
            {
                adinfcount_ = 0;
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
    }

    if (!optimal)
    {
        for (auto & pair : connectionPoint_)
        {
            base::Cost temp = opt_->combineCosts(pair.first->cost, pair.second->cost);
            if (opt_->isCostBetterThan(temp, bestCost_))
            {
                if (isPathValid(pair.first, pair.second))
                {
                    temp = opt_->combineCosts(pair.first->cost, pair.second->cost);
                    if (opt_->isCostBetterThan(temp, bestCost_))
                    {
                        bestCost_ = temp;
                        if (solved_)
                            OMPL_INFORM("%s: Found a better solution with a cost of %.2f in %u iterations (%u "
                                    "vertices in the graph, %u start + %u goal)",
                                    getName().c_str(), bestCost_.value(), iterations_, tStart_->size() + tGoal_->size(), tStart_->size(), tGoal_->size());
                        else 
                        {
                            OMPL_INFORM("%s: Found an initial solution with a cost of %.2f in %u iterations (%u "
                                    "vertices in the graph)",
                                    getName().c_str(), bestCost_.value(), iterations_, tStart_->size() + tGoal_->size());
                            for (auto & motion : startMotions_)
                                optimalRewireTree(bh_, motion, true);
                            for (auto & motion : goalMotions_)
                                optimalRewireTree(bh_, motion, false);
                            checkedStartPath_.clear();
                            checkedGoalPath_.clear();
                        }

                        solved_ = true;
                        updatedSolution = true;
                        bestStartMotion_ = pair.first;
                        bestGoalMotion_ = pair.second;

                        if (opt_->isSatisfied(bestCost_))
                        {
                            optimal = true;
                            break;
                        }
                        else 
                        {
                            infcount_ = 0;
                            clearStartInfSampler();
                            clearGoalInfSampler();
                            Motion *rtSt = nullptr, *rtG = nullptr;
                            for (auto & stMotion : startMotions_)
                            {
                                if (stMotion->root == bestStartMotion_->root)
                                {
                                    rtSt = stMotion;
                                    break;
                                }
                            }
                            for (auto & gMotion : goalMotions_)
                            {
                                if (gMotion->root == bestGoalMotion_->root)
                                {
                                    rtG = gMotion;
                                    break;
                                }
                            }
                            optimalInfSampler(rtSt, bestStartMotion_, true, startInfSamplers_);
                            optimalInfSampler(rtG,  bestGoalMotion_, false, goalInfSamplers_);
                            calculateOptimalInfProb();
                        }
                    }
                }
                //rewirePath(); // todo
            }
        }
    }

    return updatedSolution;
}

void ompl::geometric::BiASEstar::rewirePath() // todo update cmotion
{
    {
        bool updated = false;
        std::vector<Motion *> &checkedPath = checkedStartPath_;
        std::size_t i = checkedPath.size() - 1, j = i - 2, k = 0;
        if (j < checkedPath.size() && checkedPath[j]->valid && checkedPath[j + 1]->valid)
        {
            while (i < checkedPath.size())
            {
                j = i - 2;
                while (j < checkedPath.size() && checkedPath[j]->valid && opt_->isFinite(checkedPath[j]->cost))
                {
                    if (isInvalidNeighbor(checkedPath[i], checkedPath[j]))
                        break;
                    updated = true;
                    k = j;
                    j--;
                }
                if (updated)
                    break;
                i--;
            }
        }
        if (updated)
        {
            base::Cost incCost = opt_->motionCost(checkedPath[i]->state, checkedPath[k]->state);
            base::Cost cost = opt_->combineCosts(checkedPath[i]->cost, incCost);
            if (1.05 * cost.value() < checkedPath[k]->cost.value())
            {
                removeFromParent(checkedPath[k]);
                checkedPath[k]->parent = checkedPath[i];
                checkedPath[k]->incCost = incCost;
                checkedPath[k]->cost = cost;
                checkedPath[k]->root = checkedPath[i]->root;
                updateChildCosts(checkedPath[k]);
                checkedPath[k]->parent->children.push_back(checkedPath[k]);
                checkedPath[k]->valid = isValidNeighbor(checkedPath[k], checkedPath[i]);
            }
        }
    }
    {
        bool updated = false;
        std::vector<Motion *> &checkedPath = checkedGoalPath_;
        std::size_t i = checkedPath.size() - 1, j = i - 2, k = 0;
        if (j < checkedPath.size() && checkedPath[j]->valid && checkedPath[j + 1]->valid)
        {
            while (i < checkedPath.size())
            {
                j = i - 2;
                while (j < checkedPath.size() && checkedPath[j]->valid && opt_->isFinite(checkedPath[j]->cost))
                {
                    if (isInvalidNeighbor(checkedPath[i], checkedPath[j]))
                        break;
                    updated = true;
                    k = j;
                    j--;
                }
                if (updated)
                    break;
                i--;
            }
        }
        if (updated)
        {
            base::Cost incCost = opt_->motionCost(checkedPath[k]->state, checkedPath[i]->state);
            base::Cost cost = opt_->combineCosts(checkedPath[i]->cost, incCost);
            if (1.05 * cost.value() < checkedPath[k]->cost.value())
            {
                removeFromParent(checkedPath[k]);
                checkedPath[k]->parent = checkedPath[i];
                checkedPath[k]->incCost = incCost;
                checkedPath[k]->cost = cost;
                checkedPath[k]->root = checkedPath[i]->root;
                updateChildCosts(checkedPath[k]);
                checkedPath[k]->parent->children.push_back(checkedPath[k]);
                checkedPath[k]->valid = isValidNeighbor(checkedPath[k], checkedPath[i]);
            }
        }
    }
}

void ompl::geometric::BiASEstar::reportBetterSolution(const base::ReportIntermediateSolutionFn &intermediateSolutionCallback)
{
    if (!intermediateSolutionCallback)
        return;
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

double ompl::geometric::BiASEstar::improvementRatio(const base::Cost &temp, const base::State *sm, const base::State *gm) const
{
    base::Cost min = opt_->combineCosts(opt_->identityCost(), opt_->motionCost(sm, gm));
    double ratio = std::abs((temp.value() - min.value()) / temp.value());
    if (opt_->getCostThreshold().value() != 0.0)
        ratio = 0.7 * std::min(ratio, std::abs((temp.value() - opt_->getCostThreshold().value()) / temp.value()));
    else 
        ratio *= 0.5;
    return ratio;
}

void ompl::geometric::BiASEstar::optimalRewireTree(BinaryHeap<Motion *, MotionCompare> &bh, Motion *m, bool start)
{
    std::vector<Motion *> &pnullMotions = start ? pnullStartMotions_ : pnullGoalMotions_;
    updateQueue(bh, m);
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
                Cell *c = cv[ic];
                ic--;
                if (rewireSort_ && motion->cell != c && c->data->mmotion)
                {
                    if (opt_->isCostBetterThan(c->data->mmotion->cost, motion->cell->data->cmotion->cost))
                        continue;
                    base::Cost cost1 = opt_->combineCosts(motion->cost, opt_->motionCost(motion->state, c->data->cmotion->state));
                    base::Cost cost2 = opt_->combineCosts(c->data->mmotion->cost, opt_->combineCosts(c->data->mmotion->cost, c->data->mmotion->cost));
                    if (opt_->isCostBetterThan(cost2, cost1))
                        continue;
                }
                for (auto & nb : c->data->motions)
                {
                    if (motion == nb || isInvalidNeighbor(motion, nb))
                        continue;
                    Motion *pmotion = nullptr;
                    base::Cost nbhIncCost, nbhNewCost;
                    bool feas = false;
                    if (backRewire(motion, nb, start, pmotion, nbhIncCost, nbhNewCost, feas))
                    {
                        if (!opt_->isFinite(nb->cost))
                            enableMotionInDisc(nb);
                        if (!nb->parent)
                            removeFromPnull(pnullMotions, nb);
                        else
                            removeFromParent(nb);
                        nb->parent = pmotion;
                        nb->root = pmotion->root;
                        nb->incCost = nbhIncCost;
                        nb->cost = nbhNewCost;
                        updateChildCosts(nb);
                        nb->parent->children.push_back(nb);
                        nb->valid = feas;
                        updateQueue(bh, nb);
                        if (!nb->children.empty())
                            updateLeafQueue(bh, nb);
                    }
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

void ompl::geometric::BiASEstar::updateLeafQueue(BinaryHeap<Motion *, MotionCompare> &bh, Motion *motion)
{
    if (motion->children.empty())
        updateQueue(bh, motion);
    else 
    {
        for (auto & child : motion->children)
            updateLeafQueue(bh, child);
    }
}

bool ompl::geometric::BiASEstar::backRewire(Motion *motion, Motion *nb, bool start, Motion *&pmotion, base::Cost &nbhIncCost, base::Cost &nbhNewCost, bool &feas)
{
    bool valid = false;
    nbhIncCost = start ? opt_->motionCost(motion->state, nb->state) : opt_->motionCost(nb->state, motion->state);
    nbhNewCost = opt_->combineCosts(motion->cost, nbhIncCost);
    feas = isValidNeighbor(motion, nb);
    if (opt_->isCostBetterThan(nbhNewCost, nb->cost))
    {
        valid = true;
        pmotion = motion;
    }
    if (!feas)
        return valid;
    Motion *temp = nullptr;
    while (motion->valid && motion->parent)
    {
        if (isInvalidNeighbor(motion->parent, nb))
            break;
        base::Cost incCost = start ? opt_->motionCost(motion->parent->state, nb->state) : opt_->motionCost(nb->state, motion->parent->state);
        if (1.05 * incCost.value() < opt_->combineCosts(motion->incCost, nbhIncCost).value())
            break;
        base::Cost cost = opt_->combineCosts(motion->parent->cost, incCost);
        if (!opt_->isCostBetterThan(cost, nb->cost))
            break;
        temp = motion->parent;
        nbhIncCost = incCost;
        nbhNewCost = cost;
        motion = motion->parent;
    }
    if (temp)
    {
        valid = true;
        pmotion = temp;
        feas = isValidNeighbor(temp, nb);
    }
    return valid;
}

// optimal prune tree
std::size_t ompl::geometric::BiASEstar::pruneTree(const base::Cost &pruneTreeCost)
{
    double fracBetter;
    std::size_t numPruned = 0;
    if (opt_->isFinite(prunedCost_))
        fracBetter = std::abs((pruneTreeCost.value() - prunedCost_.value()) / prunedCost_.value());
    else
        fracBetter = 1.0;
    if (fracBetter > pruneThreshold_)
    {
        tStart_->clear();
        tGoal_->clear();

        bool invalid = false;
        std::queue<Motion *, std::deque<Motion *>> motionQueue;
        for (auto it = pnullStartMotions_.begin(); it != pnullStartMotions_.end();)
        {
            Motion *motion = *it;
            motionQueue.push(motion);
            numPruned += pruneTreeInternalDisabled(tStart_, dStart_, pruneTreeCost, true, motionQueue, invalid);
            if (invalid)
            {
                std::iter_swap(it, pnullStartMotions_.end() - 1);
                pnullStartMotions_.pop_back();
            }
            else 
                it++;
        }
        std::unordered_set<Cell *> cells;
        EnableSort esort(opt_);
        numPruned += pruneSingleTree(tStart_, dStart_, pruneTreeCost, true, startMotions_, cells);
        if (rewireSort_)
        {
            for (auto & cell : cells)
            {
                if (cell->data->disabled)
                    std::sort(cell->data->motions.begin(), cell->data->motions.end(), esort);
                if (!cell->data->cmotion)
                    cell->data->cmotion = *std::min_element(cell->data->motions.begin(), cell->data->motions.end() - cell->data->disabled, mc_);
                if (!cell->data->mmotion)
                    cell->data->mmotion = *std::max_element(cell->data->motions.begin(), cell->data->motions.end() - cell->data->disabled, mc_);
            }
        }
        
        for (auto it = pnullGoalMotions_.begin(); it != pnullGoalMotions_.end();)
        {
            Motion *motion = *it;
            motionQueue.push(motion);
            numPruned += pruneTreeInternalDisabled(tGoal_, dGoal_, pruneTreeCost, false, motionQueue, invalid);
            if (invalid)
            {
                std::iter_swap(it, pnullGoalMotions_.end() - 1);
                pnullGoalMotions_.pop_back();
            }
            else 
                it++;
        }
        cells.clear();
        numPruned += pruneSingleTree(tGoal_, dGoal_, pruneTreeCost, false, goalMotions_, cells);
        if (rewireSort_)
        {
            for (auto & cell : cells)
            {
                if (cell->data->disabled)
                    std::sort(cell->data->motions.begin(), cell->data->motions.end(), esort);
                if (!cell->data->cmotion)
                    cell->data->cmotion = *std::min_element(cell->data->motions.begin(), cell->data->motions.end() - cell->data->disabled, mc_);
                if (!cell->data->mmotion)
                    cell->data->mmotion = *std::max_element(cell->data->motions.begin(), cell->data->motions.end() - cell->data->disabled, mc_);
            }
        }
        prunedCost_ = pruneTreeCost;
    }
    return numPruned;
}

std::size_t ompl::geometric::BiASEstar::pruneSingleTree(TreeData &tree, CellDiscretizationData &disc, const base::Cost &pruneTreeCost, bool start,
                                                  const std::vector<Motion *> &rootMotions, std::unordered_set<Cell *> &cells)
{
    std::size_t numPruned = 0;
    std::queue<Motion *, std::deque<Motion *>> motionQueue;
    for (auto & rootMotion : rootMotions)
    {
        tree->add(rootMotion);
        addChildrenToList(&motionQueue, rootMotion);
    }
    numPruned += pruneTreeInternal(tree, disc, pruneTreeCost, start, motionQueue, cells);
    return numPruned;
}

std::size_t ompl::geometric::BiASEstar::pruneTreeInternal(TreeData &tree, CellDiscretizationData &disc, const base::Cost &pruneTreeCost, bool start,
        std::queue<Motion *, std::deque<Motion *>> &motionQueue, std::unordered_set<Cell *> &cells)
{
    std::size_t numPruned = 0;
    std::queue<Motion *, std::deque<Motion *>> leavesToPrune;
    std::list<Motion *> chainsToRecheck;
    toPrune(tree, motionQueue, leavesToPrune, chainsToRecheck, pruneTreeCost, start);
    while (!leavesToPrune.empty())
    {
        while (!leavesToPrune.empty())
        {
            pruneMotion(leavesToPrune.front(), disc, start, cells);
            leavesToPrune.pop();
            ++numPruned;
        }
        auto mIter = chainsToRecheck.begin();
        while (mIter != chainsToRecheck.end())
        {
            if ((*mIter)->children.empty())
            {
                leavesToPrune.push(*mIter);
                mIter = chainsToRecheck.erase(mIter);
            }
            else
                ++mIter;
        }
    }
    for (const auto & r : chainsToRecheck)
        tree->add(r);
    return numPruned;
}

std::size_t ompl::geometric::BiASEstar::pruneTreeInternalDisabled(TreeData &tree, CellDiscretizationData &disc, const base::Cost &pruneTreeCost, bool start,
        std::queue<Motion *, std::deque<Motion *>> &motionQueue, bool &invalid)
{
    invalid = false;
    std::size_t numPruned = 0;
    std::queue<Motion *, std::deque<Motion *>> leavesToPrune;
    std::list<Motion *> chainsToRecheck;
    toPrune(tree, motionQueue, leavesToPrune, chainsToRecheck, pruneTreeCost, start);
    while (!leavesToPrune.empty())
    {
        while (!leavesToPrune.empty())
        {
            if (!leavesToPrune.front()->parent)
            {
                invalid = true;
                removeFromVector(leavesToPrune.front()->pmotion->pchildren, leavesToPrune.front());
            }
            leavesToPrune.front()->cell->data->disabled--;
            pruneMotionDisabled(leavesToPrune.front(), disc, start);
            leavesToPrune.pop();
            ++numPruned;
        }
        auto mIter = chainsToRecheck.begin();
        while (mIter != chainsToRecheck.end())
        {
            if ((*mIter)->children.empty())
            {
                leavesToPrune.push(*mIter);
                mIter = chainsToRecheck.erase(mIter);
            }
            else
                ++mIter;
        }
    }
    for (const auto & r : chainsToRecheck)
        tree->add(r);
    return numPruned;
}

void ompl::geometric::BiASEstar::toPrune(TreeData &tree, std::queue<Motion *, std::deque<Motion *>> &motionQueue,
        std::queue<Motion *, std::deque<Motion *>> &leavesToPrune, std::list<Motion *> &chainsToRecheck, const base::Cost &pruneTreeCost, bool start)
{
    while (!motionQueue.empty())
    {
        if (keepCondition(motionQueue.front(), pruneTreeCost, start))
        {
            tree->add(motionQueue.front());
            addChildrenToList(&motionQueue, motionQueue.front());
        }
        else
        {
            if (!motionQueue.front()->children.empty())
            {
                bool keepAChild = false;
                for (unsigned int i = 0u; keepAChild == false && i < motionQueue.front()->children.size(); ++i)
                    keepAChild = keepCondition(motionQueue.front()->children.at(i), pruneTreeCost, start);
                if (keepAChild)
                    tree->add(motionQueue.front());
                else
                    chainsToRecheck.push_back(motionQueue.front());
                addChildrenToList(&motionQueue, motionQueue.front());
            }
            else
                leavesToPrune.push(motionQueue.front());
        }
        motionQueue.pop();
    }
}

void ompl::geometric::BiASEstar::pruneMotion(Motion *motion, CellDiscretizationData &disc, bool start, std::unordered_set<Cell *> &cells)
{
    if (rewireSort_)
    {
        if (motion->cell->data->cmotion == motion)
            motion->cell->data->cmotion = nullptr;
        if (motion->cell->data->mmotion == motion)
            motion->cell->data->mmotion = nullptr;
        if (motion->cell->data->motions.size() > 1)
        {
            if (!motion->cell->data->cmotion || !motion->cell->data->mmotion)
                cells.insert(motion->cell);
        }
        else 
            cells.erase(motion->cell);
    }
    pruneMotionDisabled(motion, disc, start);
}

void ompl::geometric::BiASEstar::pruneMotionDisabled(Motion *motion, CellDiscretizationData &disc, bool start)
{
    if (motion->inConnection)
    {
        for (auto it = connectionPoint_.begin(); it != connectionPoint_.end();)
        {
            if (start ? (*it).first == motion : (*it).second == motion)
            {
                (*it).first->inConnection = false;
                (*it).second->inConnection = false;
                std::iter_swap(it, connectionPoint_.end() - 1);
                connectionPoint_.pop_back();
            }
            else 
                ++it;
        }
    }
    removeFromNeighbor(motion);
    removeFromInvalidNeighbor(motion);
    removeFromParent(motion);
    removeFromDisc(disc, motion);
    freeMotion(motion);
}

void ompl::geometric::BiASEstar::addChildrenToList(std::queue<Motion *, std::deque<Motion *>> *motionList, Motion *motion)
{
    for (auto & child : motion->children)
        motionList->push(child);
}

bool ompl::geometric::BiASEstar::keepCondition(Motion *motion, const base::Cost &threshold, bool start) // todo
{
    if (start && motion == bestStartMotion_)
        return true;
    if (!start && motion == bestGoalMotion_)
        return true;
    if (!motion->pchildren.empty())
        return true;

    bool keep = !opt_->isCostBetterThan(threshold, solutionHeuristic2(motion, start));
    return keep;

    if (!useBispace_)
        return keep;
    if (keep)
    {
        if (start && opt_->isFinite(bestStartMotion_->cost))
        {
            if (opt_->isFinite(motion->cost))
                keep = !opt_->isCostBetterThan(bestStartMotion_->cost, base::Cost(0.5 * motion->cost.value()));
            if (keep)
                keep = !opt_->isCostBetterThan(bestStartMotion_->cost, base::Cost(0.75 * bordersolutionHeuristic(motion, start).value()));
        }
        else if (!start && opt_->isFinite(bestGoalMotion_->cost)) 
        {
            if (opt_->isFinite(motion->cost))
                keep = !opt_->isCostBetterThan(bestGoalMotion_->cost, base::Cost(0.5 * motion->cost.value()));
            if (keep)
                keep = !opt_->isCostBetterThan(bestGoalMotion_->cost, base::Cost(0.75 * bordersolutionHeuristic(motion, start).value()));
        }
    }

    return keep;
}

ompl::base::Cost ompl::geometric::BiASEstar::bordersolutionHeuristic(Motion *motion, bool start) const
{
    base::Cost costToCome = calculateCostToCome(motion, start);

    const base::State *rootSt = bestStartMotion_->root;
    const base::State *rootG = bestGoalMotion_->root;

    if (!treatedAsMultiSubapce_ || opt_->getDescription() != "Path Length")
    {
        double dist1 = si_->distance(rootSt, motion->state);
        double dist2 = si_->distance(motion->state, rootG);

        if (start != (dist1 <= dist2))
            return costToCome;

        base::State *test = si_->allocState();
        base::State *state1 = si_->allocState();
        base::State *state2 = si_->allocState();

        si_->copyState(state1, motion->state);

        double d = 0.;
        if (start)
        {
            si_->copyState(state2, rootG);
            si_->getStateSpace()->interpolate(state1, state2, 0.5, test);
            d = si_->distance(state1, state2);
        }
        else 
        {
            si_->copyState(state2, rootSt);
            si_->getStateSpace()->interpolate(state2, state1, 0.5, test);
            d = si_->distance(state2, state1);
        }

        while (d > 1.e-6)
        {
            dist1 = si_->distance(rootSt, test);
            dist2 = si_->distance(test, rootG);

            if (start != (dist1 <= dist2))
                si_->copyState(state2, test);
            else 
                si_->copyState(state1, test);

            if (start)
                si_->getStateSpace()->interpolate(state1, state2, 0.5, test);
            else 
                si_->getStateSpace()->interpolate(state2, state1, 0.5, test);

            d *= 0.5;
        }

        base::Cost costToGo = opt_->infiniteCost();
        if (start)
            costToGo = opt_->motionCost(motion->state, test);
        else 
            costToGo = opt_->motionCost(test, motion->state);

        si_->freeState(test);
        si_->freeState(state1);
        si_->freeState(state2);

        return opt_->combineCosts(costToCome, costToGo);
    }
    else 
    {
        for (std::size_t i = 0; i < si_->getStateSpace()->getSubspaceCount(); i++)
        {
            double dist1 = si_->distance(rootSt, motion->state, i);
            double dist2 = si_->distance(motion->state, rootG, i);

            if (start != (dist1 <= dist2))
                return costToCome;
        }

        base::State *test = si_->allocState();
        base::State *state1 = si_->allocState();
        base::State *state2 = si_->allocState();
        base::Cost costToGo = opt_->identityCost();

        for (std::size_t i = 0; i < si_->getStateSpace()->getSubspaceCount(); i++)
        {
            si_->copyState(state1, motion->state, i);

            double d = 0.;
            if (start)
            {
                si_->copyState(state2, rootG, i);
                si_->getStateSpace()->interpolate(state1, state2, 0.5, test, i);
                d = si_->distance(state1, state2, i);
            }
            else 
            {
                si_->copyState(state2, rootSt, i);
                si_->getStateSpace()->interpolate(state2, state1, 0.5, test, i);
                d = si_->distance(state2, state1, i);
            }

            while (d > 1.e-6)
            {
                double dist1 = si_->distance(rootSt, test, i);
                double dist2 = si_->distance(test, rootG, i);

                if (start != (dist1 <= dist2))
                    si_->copyState(state2, test, i);
                else 
                    si_->copyState(state1, test, i);

                if (start)
                    si_->getStateSpace()->interpolate(state1, state2, 0.5, test, i);
                else 
                    si_->getStateSpace()->interpolate(state2, state1, 0.5, test, i);

                d *= 0.5;
            }

            if (start)
                costToGo = opt_->combineCosts(costToGo, opt_->motionCost(motion->state, test, i));
            else 
                costToGo = opt_->combineCosts(costToGo, opt_->motionCost(test, motion->state, i));
        }

        si_->freeState(test);
        si_->freeState(state1);
        si_->freeState(state2);

        return opt_->combineCosts(costToCome, costToGo);
    }
}

bool ompl::geometric::BiASEstar::keepCondition2(Motion *motion, const base::Cost &threshold, bool start) const
{
    return !opt_->isCostBetterThan(threshold, solutionHeuristic2(motion, start));
}

ompl::base::Cost ompl::geometric::BiASEstar::solutionHeuristic2(Motion *motion, bool start) const
{
    base::Cost costToCome = calculateCostToCome(motion, start);
    base::Cost costToGo = calculateCostToGo(motion, start);
    return opt_->combineCosts(costToCome, costToGo);
}

ompl::base::Cost ompl::geometric::BiASEstar::calculateCostToCome(Motion *motion, bool start) const
{
    base::Cost costToCome = motion->cost;

    if (!opt_->isFinite(costToCome))
    {
        costToCome = opt_->infiniteCost();

        if (start)
        {
            for (auto & startMotion : startMotions_)
                costToCome = opt_->betterCost(costToCome, opt_->motionCost(startMotion->state, motion->state));
        }
        else 
        {
            for (auto & goalMotion : goalMotions_)
                costToCome = opt_->betterCost(costToCome, opt_->motionCost(motion->state, goalMotion->state));
        }
    }

    return costToCome;
}

ompl::base::Cost ompl::geometric::BiASEstar::calculateCostToGo(Motion *motion, bool start) const
{
    base::Cost costToGo = opt_->infiniteCost();

    if (start)
    {
        for (auto & goalMotion : goalMotions_)
            costToGo = opt_->betterCost(costToGo, opt_->motionCost(motion->state, goalMotion->state));
    }
    else 
    {
        for (auto & startMotion : startMotions_)
            costToGo = opt_->betterCost(costToGo, opt_->motionCost(startMotion->state, motion->state));
    }

    return costToGo;
}


// adaptive informed sampling
bool ompl::geometric::BiASEstar::batchGrow(bool &startTree, bool ais, bool &adinf)
{
    bool add = false;
    TreeGrowingInfo tgi;
    tgi.xstate = si_->allocState();
    Motion *rmotion = new Motion(si_);
    base::State *rstate = rmotion->state;
    for (unsigned int i = 0; i < 1; i++)
    {
        bool ad = false;
        if (!solved_ || (ais && adinf))
        {
            if (tree_ == 0)
                startTree = true;
            else if (tree_ == 1)
                startTree = false;
            adinf = !adinf;
            ad = true;
        }

        tgi.start = startTree;
        startTree = !startTree;
        guniform_ = true;
        if (ad)
        {
            if (!sampleUniformAd(rstate, tgi.start)) 
                continue;
        }
        else if (!sampleUniform(rstate, tgi.start))
            continue;
        if (!guniform_)
        {
            if (ad)
                adinfcount_++;
            else
                infcount_++;
        }

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

bool ompl::geometric::BiASEstar::isPathValid(Motion *motion, Motion *otherMotion)
{
    bool valid = true;
    currentStartCost_ = motion->cost;
    currentGoalCost_  = otherMotion->cost;
    checkedStartPath_.clear();
    checkedGoalPath_.clear();
    if (!isPathValid(motion, true))
        valid = false;
    if (!isPathValid(otherMotion, false))
        valid = false;
    return valid;
}

bool ompl::geometric::BiASEstar::isPathValid(Motion *motion, bool start)
{
    bool tvalid = true;
    std::vector<Motion *> mpath;
    std::vector<Motion *> &checkedPath = start ? checkedStartPath_ : checkedGoalPath_;
    while (motion)
    {
        checkedPath.push_back(motion);
        mpath.push_back(motion);
        motion = motion->parent;
    }
    mpath.pop_back();
    base::Cost &currentCost = start ? currentStartCost_ : currentGoalCost_;
    std::vector<Motion *> &pnullMotions= start ? pnullStartMotions_: pnullGoalMotions_;
    for (std::size_t i = mpath.size() - 1; i < mpath.size(); i--)
    {
        motion = mpath[i];
        Motion *pmotion = motion->parent;
        base::Cost cost = motion->cost;
        if (!checkMotion(pmotion, motion, start))
        {
            removeFromParent(motion);
            motion->parent = nullptr;
            bool stop = false;
            if (tvalid)
            {
                tvalid = false;
                if (backRewire_)
                {
                    Motion *ppmotion = nullptr;
                    tvalid = backPathRewireMotion(motion, start, ppmotion);
                    if (tvalid)
                    {
                        checkedPath.resize(i+1);
                        Motion *last = ppmotion;
                        while (last)
                        {
                            checkedPath.push_back(last);
                            last = last->parent;
                        }
                        base::Cost nbhIncCost = start ? opt_->motionCost(ppmotion->state, motion->state) : opt_->motionCost(motion->state, ppmotion->state);
                        base::Cost nbhNewCost = opt_->combineCosts(ppmotion->cost, nbhIncCost);
                        currentCost = opt_->combineCosts(currentCost, base::Cost(nbhNewCost.value() - cost.value()));
                        base::Cost temp = opt_->combineCosts(currentStartCost_, currentGoalCost_);
                        if (!opt_->isCostBetterThan(temp, bestCost_))
                            stop = true;
                        else 
                            tvalid = isPathValidInter(ppmotion, start);
                        connectToPmotion(motion, ppmotion, start);
                        motion->parent->children.push_back(motion); 
                    }
                }
            }
            if (tvalid)
                enableMotionInDisc(motion);
            else if (!motion->parent)
            {
                pnullMotions.push_back(motion);
                motion->valid = false;
                motion->pmotion = pmotion;
                motion->pmotion->pchildren.push_back(motion);
            }
            if (stop)
            {
                tvalid = false;
                checkedPath.clear();
                break;
            }
        }
    }
    return tvalid;
}

bool ompl::geometric::BiASEstar::isPathValidInter(Motion *motion, bool start)
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
            motion->pmotion->pchildren.push_back(motion);
        }
    }
    return tvalid;
}

bool ompl::geometric::BiASEstar::selectCMotion(std::size_t &index, bool &reverse)
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
            if (ratio > 0.2 && pair.first != bestStartMotion_) // todo
            {
                index = iter;
                nconnect = true;
                reverse = !reverse;
                break;
            }
        }
        if (reverse)
            iter--;
        else 
            iter++;
    }
    return nconnect;
}

void ompl::geometric::BiASEstar::processAdEllipsoidRind(bool clearoradd, bool &ais)
{
    if (clearoradd)
    {
        startAdInfPdf_.clear();
        startAdElems_.clear();

        goalAdInfPdf_.clear();
        goalAdElems_.clear();

        calculateInfProb(false);
        ais = true;
        adinfcount_ = 0;
    }
    else if (ais)
    {
        if (adinfcount_ >= 100u) // todo
        {
            clearStartAdInfSampler();
            clearGoalAdInfSampler();
            ais = false;
            tree_ = -1;
            adinfcount_ = 0;
        }
        else if (adinfcount_ == 50u) // todo
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

void ompl::geometric::BiASEstar::localInfeasible(int &tree, bool &locals, bool &localg)
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

void ompl::geometric::BiASEstar::calculateInfSampler(bool local, bool start, bool &clearoradd)
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

void ompl::geometric::BiASEstar::startInfSampler(bool local, std::vector<base::AdInformedSamplerPtr> &infSamplers)
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
//        double ratio = (double)tStart_->size() / (double)(tStart_->size() + tGoal_->size()); // todo
//        if (useInformedSampling_ && sampler->getInformedMeasure() >= ratio * si_->getSpaceMeasure())
        if (useInformedSampling_ && sampler->getInformedMeasure() >= si_->getSpaceMeasure())
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

void ompl::geometric::BiASEstar::startLocalInfSampler(std::vector<base::AdInformedSamplerPtr> &infSamplers)
{
    Motion *motion1 = nullptr, *motion2 = nullptr;
    std::size_t localseg = std::ceil(0.1 * (double)(checkedStartPath_.size() - 1)); // todo
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

void ompl::geometric::BiASEstar::goalInfSampler(bool local, std::vector<base::AdInformedSamplerPtr> &infSamplers)
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
//        double ratio = (double)tGoal_->size() / (double)(tStart_->size() + tGoal_->size()); // todo
//        if (useInformedSampling_ && sampler->getInformedMeasure() >= ratio * si_->getSpaceMeasure())
        if (useInformedSampling_ && sampler->getInformedMeasure() >= si_->getSpaceMeasure())
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

void ompl::geometric::BiASEstar::goalLocalInfSampler(std::vector<base::AdInformedSamplerPtr> &infSamplers)
{
    Motion *motion1 = nullptr, *motion2 = nullptr;
    std::size_t localseg = std::ceil(0.1 * (double)(checkedGoalPath_.size() - 1)); // todo
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

ompl::base::AdInformedSamplerPtr ompl::geometric::BiASEstar::allocInfSampler(const base::State *s1, const base::State *s2,
                                                                       const base::Cost &minCost, const base::Cost &maxCost)
{
    if (useRejectionSampling_)
        return std::make_shared<base::RejectionAdInfSampler>(pdef_, s1, s2, minCost, maxCost, numSampleAttempts_);
    else 
        return std::make_shared<base::PathLengthDirectAdInfSampler>(pdef_, s1, s2, minCost, maxCost, numSampleAttempts_);
}

void ompl::geometric::BiASEstar::calculateInfProb(bool update)
{
    calculateInfProb(startAdInfSamplers_, goalAdInfSamplers_, startAdInfPdf_, goalAdInfPdf_, startAdElems_, goalAdElems_, update);
}

void ompl::geometric::BiASEstar::calculateInfProb(const std::vector<base::AdInformedSamplerPtr> &startInfSamplers, 
                                              const std::vector<base::AdInformedSamplerPtr> &goalInfSamplers,
                                              NumPdf &startInfPdf, NumPdf &goalInfPdf, 
                                              std::vector<NumElem *> &startelems, std::vector<NumElem *> &goalelems, bool update)
{
    if (!startInfSamplers.empty())
    {
        double measure = 0;
        for (auto & sampler : startInfSamplers)
            measure += sampler->getInformedMeasure();
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
}

bool ompl::geometric::BiASEstar::sampleUniformAd(base::State *state, bool start)
{
    if (start)
        return sampleUniform(state, startAdInfSamplers_, startAdInfPdf_);
    return sampleUniform(state, goalAdInfSamplers_, goalAdInfPdf_);
}

bool ompl::geometric::BiASEstar::sampleUniform(base::State *state, const std::vector<base::AdInformedSamplerPtr> &infSamplers, const NumPdf &adInfPdf)
{
    if (!infSamplers.empty())
    {
        guniform_ = false;
        if (infSamplers.size() == 1)
            return infSamplers[0]->sampleUniform(state);
        else
        {
            std::size_t ii = adInfPdf.sample(rng_.uniform01());
            return infSamplers[ii]->sampleUniform(state);
        }
    }
    else
    {
        guniform_ = true;
        sampler_->sampleUniform(state);
        return true;
    }
}


// informed sampling
void ompl::geometric::BiASEstar::optimalInfSampler(const Motion *motion1, const Motion *motion2, bool start, std::vector<base::AdInformedSamplerPtr> &infSamplers)
{
    base::Cost cost = base::Cost(motion2->cost.value() - motion1->cost.value());
    while (start ? cost.value() <= si_->distance(motion1->state, motion2->state) : cost.value() <= si_->distance(motion2->state, motion1->state))
        cost = base::Cost(factor_ * cost.value());
    base::AdInformedSamplerPtr sampler = allocInfSampler(motion1->state, motion2->state, base::Cost(0.0), cost);
    if (useInformedSampling_ && sampler->getInformedMeasure() >= si_->getSpaceMeasure())
    {
        unsigned int segment = 0;
        const Motion *temp = motion2;
        while (temp->parent)
        {
            temp = temp->parent;
            segment++;
            if (temp == motion1)
                break;
        }
        if (segment <= 1)
            return;
        unsigned int n = 0;
        temp = motion2;
        while (temp->parent && n < segment/2)
        {
            temp = temp->parent;
            n++;
        }
        optimalInfSampler(motion1, temp, start, infSamplers);
        optimalInfSampler(temp, motion2, start, infSamplers);
    }
    else
        infSamplers.push_back(sampler);
}

void ompl::geometric::BiASEstar::calculateOptimalInfProb()
{
    calculateInfProb(startInfSamplers_, goalInfSamplers_, startInfPdf_, goalInfPdf_, startElems_, goalElems_, false);
}

bool ompl::geometric::BiASEstar::sampleUniform(base::State *state, bool start)
{
    if (start)
        return sampleUniform(state, startInfSamplers_, startInfPdf_);
    return sampleUniform(state, goalInfSamplers_, goalInfPdf_);
}
