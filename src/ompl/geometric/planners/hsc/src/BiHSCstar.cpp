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

#include "ompl/geometric/planners/hsc/BiHSCstar.h"
#include "ompl/base/goals/GoalSampleableRegion.h"
#include "ompl/base/objectives/PathLengthOptimizationObjective.h"
#include "ompl/tools/config/SelfConfig.h"

ompl::geometric::BiHSCstar::BiHSCstar(const base::SpaceInformationPtr &si) : base::Planner(si, "BiHSCstar")
  , dStart_(0), dGoal_(0), mc_(opt_), bh_(mc_)
{
    specs_.recognizedGoal = base::GOAL_SAMPLEABLE_REGION;
    specs_.directed = true;
    specs_.approximateSolutions = false;
    specs_.optimizingPaths = true;
    specs_.canReportIntermediateSolutions = true;

    Planner::declareParam<double>("range", this, &BiHSCstar::setRange, &BiHSCstar::getRange, "0.:0.1:10000.");
    Planner::declareParam<double>("pen_distance", this, &BiHSCstar::setPenDistance, &BiHSCstar::getPenDistance, "0.:0.1:10000.");
    Planner::declareParam<bool>("lazy_path", this, &BiHSCstar::setLazyPath, &BiHSCstar::getLazyPath, "0,1");
    Planner::declareParam<bool>("lazy_node", this, &BiHSCstar::setLazyNode, &BiHSCstar::getLazyNode, "0,1");
    Planner::declareParam<bool>("add_intermediate_state", this, &BiHSCstar::setAddIntermediateState, &BiHSCstar::getAddIntermediateState, "0,1");
    Planner::declareParam<bool>("use_bispace", this, &BiHSCstar::setUseBispace, &BiHSCstar::getUseBispace, "0,1");
//    Planner::declareParam<bool>("use_biasgrow", this, &BiHSCstar::setUseBiasGrow, &BiHSCstar::getUseBiasGrow, "0,1");
//    Planner::declareParam<bool>("treated_as_multi_subapce", this, &BiHSCstar::setTreatedAsMultiSubapce, &BiHSCstar::getTreatedAsMultiSubapce, "0,1");
    Planner::declareParam<double>("prune_threshold", this, &BiHSCstar::setPruneThreshold, &BiHSCstar::getPruneThreshold, "0.:.01:1.");

    Planner::declareParam<bool>("use_collision_sc", this, &BiHSCstar::setUseCollisionCertificateChecker, &BiHSCstar::getUseCollisionCertificateChecker, "0,1");

    addPlannerProgressProperty("iterations INTEGER", [this] { return numIterationsProperty(); });
    addPlannerProgressProperty("best cost REAL", [this] { return bestCostProperty(); });

    addPlannerProgressProperty("collision check time REAL", [this] { return collisionCheckTimeProperty(); });
}

ompl::geometric::BiHSCstar::~BiHSCstar()
{
    freeMemory();
}

void ompl::geometric::BiHSCstar::setup()
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
    dStart_.setDimension(projectionEvaluator_->getDimension());
    dGoal_.setDimension(projectionEvaluator_->getDimension());
    dStart_.useNeighbor(false);
    dGoal_.useNeighbor(false);
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
    
    if (lazyNode_)
        lazyPath_ = true;
    if (!lazyPath_)
      lazyNode_ = false;
    if (useCollisionCertificateChecker_) 
    {
        onn_.reset(new SimpleGridSC());
        onn_->setDimension(projectionEvaluator_->getDimension());
    }

    symmetric_ = si_->getStateSpace()->hasSymmetricDistance() && si_->getStateSpace()->hasSymmetricInterpolate();
}

void ompl::geometric::BiHSCstar::freeMemory()
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

    if (lazyNode_)
    {
        for (auto & sc : ssnne_)
        {
            freeCertificate(sc->sc);
            delete sc;
        }
        ssnne_.clear();
        for (auto & sc : gsnne_)
        {
            freeCertificate(sc->sc);
            delete sc;
        }
        gsnne_.clear();
    }

    if (onn_)
    {
        std::vector<CellSC *> cells;
        onn_->getGrid().getCells(cells);
        for (auto & cell : cells)
        {
            for (auto &sc : cell->data->datas)
                freeCertificate(sc);
        }
        onn_->freeMemory();
    }
}

void ompl::geometric::BiHSCstar::clear()
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

    invalidStartMotions_.clear();
    invalidGoalMotions_.clear();
    invalidStartNum_ = invalidGoalNum_ = 0;

    connectionPoint_.clear();
    currentStartCost_ = base::Cost(std::numeric_limits<double>::quiet_NaN());
    currentGoalCost_ = base::Cost(std::numeric_limits<double>::quiet_NaN());

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

    oTime_ = 0.0;
    if (onn_)
    {
        onn_->clear();
    }

    ssnne_.clear();
    gsnne_.clear();
}

ompl::base::PlannerStatus ompl::geometric::BiHSCstar::solve(const base::PlannerTerminationCondition &ptc)
{
    base::PlannerStatus ps= prepareSolve(ptc);
    if (ps != base::PlannerStatus::PREPARE_SUCCESS)
        return ps;   

    sampler_ = si_->allocStateSampler();

    OMPL_INFORM("%s: Starting planning with %u states already in datastructure. Seeking a solution better than %.5f.", getName().c_str(), (tStart_->size() + tGoal_->size()),
            opt_->getCostThreshold().value());

    const base::ReportIntermediateSolutionFn intermediateSolutionCallback = pdef_->getIntermediateSolutionCallback();

    bool startTree = true;
    bool optimal = false;

    unsigned int connect1 = 0;
    double ratio1 = 0.5, maxratio1 = 0.0, connectTresh1 = 10.0;

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
                motion->valid = ValidP;
                motion->stateValid = Valid;
                motion->root = motion->state;
                motion->cost = opt_->identityCost();
                goalMotions_.push_back(motion);
                tGoal_->add(motion);
                Coord xcoord(projectionEvaluator_->getDimension());
                projectionEvaluator_->computeCoordinates(motion->state, xcoord);
                addToDisc(dGoal_, motion, xcoord);
                motion->cell->data->cmotion = motion;
                motion->cell->data->mmotion = motion;
                if (lazyNode_)
                {
                    base::SafetyCertificate *sc = new base::SafetyCertificate(si_, certificateDim_);
                    si_->copyState(sc->state, st);
                    delete sc->contact;
                    sc->contact = nullptr;

                    SafetyCertificateWithElems *sce = new SafetyCertificateWithElems();
                    sce->sc = sc;
                    sce->objects.push_back(motion);
                    gsnne_.push_back(sce);

                    motion->sce = sce;
                    motion->scd.resize(certificateDim_);
                    std::fill(motion->scd.begin(), motion->scd.end(), 0.0);
                }
            }
        }

        batchGrow(startTree);
        bool updatedSolution = findBetterSolution(optimal, ratio1, maxratio1, connect1, connectTresh1);
        if (optimal)
            break;
        if (lazyNode_)
            removeInvalidMotions();
        if (updatedSolution)
        {
            reportBetterSolution(intermediateSolutionCallback);
            int numPruned = pruneTree(bestCost_);
            if (false)
                OMPL_INFORM("%s: %u states are pruned from the tree, %u states are left", getName().c_str(), numPruned, tStart_->size() + tGoal_->size());
        }
        if (lazyNode_)
            removeInvalidMotionsTree(maxInvalidNodeRatio_);
    }

    if (solved_ || optimal)
    {
        ptc.terminate();
        if (!optimal)
            isPathValid(bestStartMotion_, bestGoalMotion_);
        processSolution(bestStartMotion_, bestGoalMotion_);
    }

    OMPL_INFORM("%s: Created %u states (%u start + %u goal). Final solution cost %.5f.", getName().c_str(), tStart_->size() + tGoal_->size(),
                tStart_->size(), tGoal_->size(), bestCost_.value());
    return (solved_ || optimal) ? base::PlannerStatus::EXACT_SOLUTION : base::PlannerStatus::TIMEOUT;
}

ompl::base::PlannerStatus ompl::geometric::BiHSCstar::prepareSolve(const base::PlannerTerminationCondition &ptc)
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
            motion->valid = ValidP;
            motion->stateValid = Valid;
            motion->root = motion->state;
            motion->cost = opt_->identityCost();
            startMotions_.push_back(motion);
            tStart_->add(motion);
            Coord xcoord(projectionEvaluator_->getDimension());
            projectionEvaluator_->computeCoordinates(motion->state, xcoord);
            addToDisc(dStart_, motion, xcoord);
            motion->cell->data->cmotion = motion;
            motion->cell->data->mmotion = motion;
            if (lazyNode_)
            {
                base::SafetyCertificate *sc = new base::SafetyCertificate(si_, certificateDim_);
                si_->copyState(sc->state, st);
                delete sc->contact;
                sc->contact = nullptr;

                SafetyCertificateWithElems *sce = new SafetyCertificateWithElems();
                sce->sc = sc;
                sce->objects.push_back(motion);
                ssnne_.push_back(sce);

                motion->sce = sce;
                motion->scd.resize(certificateDim_);
                std::fill(motion->scd.begin(), motion->scd.end(), 0.0);
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
        motion->valid = ValidP;
        motion->stateValid = Valid;
        motion->root = motion->state;
        motion->cost = opt_->identityCost();
        goalMotions_.push_back(motion);
        tGoal_->add(motion);
        Coord xcoord(projectionEvaluator_->getDimension());
        projectionEvaluator_->computeCoordinates(motion->state, xcoord);
        addToDisc(dGoal_, motion, xcoord);
        motion->cell->data->cmotion = motion;
        motion->cell->data->mmotion = motion;
        if (lazyNode_)
        {
            base::SafetyCertificate *sc = new base::SafetyCertificate(si_, certificateDim_);
            si_->copyState(sc->state, st);
            delete sc->contact;
            sc->contact = nullptr;

            SafetyCertificateWithElems *sce = new SafetyCertificateWithElems();
            sce->sc = sc;
            sce->objects.push_back(motion);
            gsnne_.push_back(sce);

            motion->sce = sce;
            motion->scd.resize(certificateDim_);
            std::fill(motion->scd.begin(), motion->scd.end(), 0.0);
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

void ompl::geometric::BiHSCstar::processSolution(const Motion *bestStartMotion, const Motion *bestGoalMotion)
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

bool ompl::geometric::BiHSCstar::batchGrow(bool &startTree)
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

bool ompl::geometric::BiHSCstar::growTree(TreeGrowingInfo &tgi, Motion *rmotion, bool &otherSide, bool &change)
{
    if (treatedAsMultiSubapce_)
        return growTreeMultiSpace(tgi, rmotion, otherSide, change);
    else 
        return growTreeSingleSpace(tgi, rmotion, otherSide, change);
}

bool ompl::geometric::BiHSCstar::growTreeSingleSpace(TreeGrowingInfo &tgi, Motion *rmotion, bool &otherSide, bool &change)
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

    std::vector<SafetyCertificateWithElems *> &snne = tgi.start ? ssnne_ : gsnne_;
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

        base::SafetyCertificate *sc = new base::SafetyCertificate(si_, certificateDim_);
        si_->copyState(sc->state, dstate);
        bool lazy = false;
        bool cvalid = true;
        std::vector<double> cdist;
        time::point starto = time::now();
        if (useCollisionCertificateChecker_)
        {
            CoordSC xcoord(projectionEvaluator_->getDimension());
            projectionEvaluator_->computeCoordinates(sc->state, xcoord);
            CellSC *cell = onn_->getCell(xcoord);
            if (cell && collisionCertificateChecker_(sc->state, cell->data->datas))
            {
                cvalid = false;
                freeCertificate(sc);
            }
        }
        if (cvalid)
        {
            if (!lazyNode_)
            {
                if (!isValid(sc))
                    cvalid = false;
            }
            // the invalid motion does not have a safety certificate node
            else if (nmotion->sce != nullptr && safetyCertificateChecker_(sc->state, nmotion->sce->sc, cdist))
                lazy = true;
            else if (!isValid(sc))
                cvalid = false;
        }
        oTime_ += time::seconds(time::now() - starto);
        if (!cvalid)
        {
            reach = false;
            break;
        }
        delete sc->contact;
        sc->contact = nullptr;

        motion = new Motion(si_);
        si_->copyState(motion->state, dstate);
        if (!lazy)
            motion->stateValid = Valid;
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
        if (!lazyPath_ && !checkInterMotion(nb, motion, tgi.start))
        {
            removeFromInvalidNeighbor(motion);
            removeFromDisc(disc, motion);
            freeMotion(motion);
            freeCertificate(sc);
            reach = false;
            break;
        }

        if (solved_ && !keepCondition2(motion, bestCost_, tgi.start))
        {
            removeFromDisc(disc, motion);
            freeMotion(motion);
            freeCertificate(sc);
            reach = false;
            break;
        }

        if (lazyNode_)
        {
            if (!lazy)
            {
                SafetyCertificateWithElems *sce = new SafetyCertificateWithElems();
                sce->sc = sc;
                snne.push_back(sce);

                motion->sce = sce;
                motion->sce->sc->confidence_ = distanceCertificate_(nb->state, motion->state);
                motion->sce->objects.push_back(motion);
                motion->scd.resize(certificateDim_);
                std::fill(motion->scd.begin(), motion->scd.end(), 0.0);
            }
            else 
            {
                freeCertificate(sc);
                motion->sce = nmotion->sce;
                motion->scd = cdist;
                motion->sce->objects.push_back(motion);
            }
        }
        else 
            freeCertificate(sc);

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
        {
            motion->valid = ValidP;
            insertNeighbor(nb, motion);
            disc.updateNbh(nb->cell, motion->cell);
        }

        if (!motion->cell->data->cmotion || opt_->isCostBetterThan(motion->cost, motion->cell->data->cmotion->cost))
            motion->cell->data->cmotion = motion;
        if (!motion->cell->data->mmotion || !opt_->isCostBetterThan(motion->cost, motion->cell->data->mmotion->cost))
            motion->cell->data->mmotion = motion;
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

bool ompl::geometric::BiHSCstar::growTreeMultiSpace(TreeGrowingInfo &tgi, Motion *rmotion, bool &otherSide, bool &change)
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

    std::vector<SafetyCertificateWithElems *> &snne = tgi.start ? ssnne_ : gsnne_;
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
                
        base::SafetyCertificate *sc = new base::SafetyCertificate(si_, certificateDim_);
        si_->copyState(sc->state, dstate);
        bool lazy = false;
        bool cvalid = true;
        std::vector<double> cdist;
        time::point starto = time::now();
        if (useCollisionCertificateChecker_)
        {
            CoordSC xcoord(projectionEvaluator_->getDimension());
            projectionEvaluator_->computeCoordinates(sc->state, xcoord);
            CellSC *cell = onn_->getCell(xcoord);
            if (cell && collisionCertificateChecker_(sc->state, cell->data->datas))
            {
                cvalid = false;
                freeCertificate(sc);
            }
        }
        if (cvalid)
        {
            if (!lazyNode_)
            {
                if (!isValid(sc))
                    cvalid = false;
            }
            // the invalid motion does not have a safety certificate node // the invalid motion does not have a safety certificate node
            else if (nmotion->sce != nullptr && safetyCertificateChecker_(sc->state, nmotion->sce->sc, cdist))
                lazy = true;
            else if (!isValid(sc))
                cvalid = false;
        }
        oTime_ += time::seconds(time::now() - starto);
        if (!cvalid)
        {
            reach = false;
            break;
        }
        delete sc->contact;
        sc->contact = nullptr;

        motion = new Motion(si_);
        si_->copyState(motion->state, dstate);
        if (!lazy)
            motion->stateValid = Valid;
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
        if (!lazyPath_ && !checkInterMotion(nb, motion, tgi.start))
        {
            removeFromInvalidNeighbor(motion);
            removeFromDisc(disc, motion);
            freeMotion(motion);
            freeCertificate(sc);
            reach = false;
            break;
        }

        if (solved_ && !keepCondition2(motion, bestCost_, tgi.start))
        {
            removeFromDisc(disc, motion);
            freeMotion(motion);
            freeCertificate(sc);
            reach = false;
            break;
        }

        if (lazyNode_)
        {
            if (!lazy)
            {
                SafetyCertificateWithElems *sce = new SafetyCertificateWithElems();
                sce->sc = sc;
                snne.push_back(sce);

                motion->sce = sce;
                motion->sce->sc->confidence_ = distanceCertificate_(nb->state, motion->state);
                motion->sce->objects.push_back(motion);
                motion->scd.resize(certificateDim_);
                std::fill(motion->scd.begin(), motion->scd.end(), 0.0);
            }
            else 
            {
                freeCertificate(sc);
                motion->sce = nmotion->sce;
                motion->scd = cdist;
                motion->sce->objects.push_back(motion);
            }
        }
        else 
            freeCertificate(sc);

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
        {
            motion->valid = ValidP;
            insertNeighbor(nb, motion);
            disc.updateNbh(nb->cell, motion->cell);
        }

        if (!motion->cell->data->cmotion || opt_->isCostBetterThan(motion->cost, motion->cell->data->cmotion->cost))
            motion->cell->data->cmotion = motion;
        if (!motion->cell->data->mmotion || !opt_->isCostBetterThan(motion->cost, motion->cell->data->mmotion->cost))
            motion->cell->data->mmotion = motion;
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

ompl::geometric::BiHSCstar::Motion *ompl::geometric::BiHSCstar::selectNMotion(const TreeData &tree, Motion *rmotion, bool &null)
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

ompl::geometric::BiHSCstar::Motion *ompl::geometric::BiHSCstar::selectMotionInCell(Cell *cell)
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

void ompl::geometric::BiHSCstar::rewireTree(BinaryHeap<Motion *, MotionCompare> &bh, Motion *m, bool start)
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
                if (motion->cell != c && c->data->mmotion)
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
                    if (isValidNeighbor(motion, nb))
                        nb->valid = ValidP;
                    else 
                        nb->valid = UnCkeckedP;
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

void ompl::geometric::BiHSCstar::updateQueue(BinaryHeap<Motion *, MotionCompare> &bh, Motion *motion)
{
    if (motion->handle != nullptr)
        bh.update(motion->handle);
    else
        motion->handle = bh.insert(motion);
}

bool ompl::geometric::BiHSCstar::growCurrentTree(const base::State *state) const
{
    return growStartTree(state);
}

bool ompl::geometric::BiHSCstar::growStartTree(const base::State *state) const
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

double ompl::geometric::BiHSCstar::penetrationDistance(const base::State *nstate, const base::State *state, bool start) const
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

bool ompl::geometric::BiHSCstar::growCurrentTree(const base::State *state, std::size_t sub) const
{
    return growStartTree(state, sub);
}

bool ompl::geometric::BiHSCstar::growStartTree(const base::State *state, std::size_t sub) const
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

double ompl::geometric::BiHSCstar::penetrationDistance(const base::State *nstate, const base::State *state, bool start, std::size_t sub) const
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

bool ompl::geometric::BiHSCstar::isPathValid(Motion *motion, Motion *otherMotion)
{
    bool valid = true;
    currentStartCost_ = motion->cost;
    currentGoalCost_  = otherMotion->cost;
    checkedStartPath_.clear();
    checkedGoalPath_.clear();
    if (!isPathValid(motion, true))
        valid = false;
//    currentStartCost_ = motion->cost;
    if (!isPathValid(otherMotion, false))
        valid = false;
    return valid;
}

bool ompl::geometric::BiHSCstar::isPathValid(Motion *motion, bool start)
{
    if (!lazyPath_)
        return true;
    bool stop = false;
    if (!isPathValidLazy(motion, start, stop) || stop)
        return false;
    bool tvalid = true;
    std::vector<Motion *> mpath;
    std::vector<Motion *> &checkedPath = start ? checkedStartPath_ : checkedGoalPath_;
    while (motion != nullptr)
    {
        checkedPath.push_back(motion);
        mpath.push_back(motion);
        motion = motion->parent;
    }
    mpath.pop_back();
    base::Cost &currentCost = start ? currentStartCost_ : currentGoalCost_;
    std::vector<Motion *> &pnullMotions= start ? pnullStartMotions_: pnullGoalMotions_;
    std::vector<Motion *> &invalidMotions = start ? invalidStartMotions_ : invalidGoalMotions_;
    std::vector<Motion *> nullMotions;
    for (std::size_t i = mpath.size() - 1; i < mpath.size(); i--)
    {
        motion = mpath[i];
        Motion *pmotion = motion->parent;
        base::Cost cost = motion->cost;
        std::size_t spos = invalidMotions.size();
        if (!checkMotion(pmotion, motion, start))
        {
            removeFromParent(motion);
            motion->parent = nullptr;
            tvalid = false;
            Motion *ppmotion = nullptr;
            tvalid = backPathRewireMotion(motion, start, ppmotion);
            std::size_t epos = invalidMotions.size();
            for (std::size_t i = spos; i < epos; i++)
            {
                Motion *temp = invalidMotions[i];
                for (auto & child : temp->children)
                {
                    child->pmotion = pmotion;
                    child->pmotion->pchildren.push_back(child);
                }
            }
            bool stop = false;
            if (tvalid)
            {
                base::Cost nbhIncCost = start ? opt_->motionCost(ppmotion->state, motion->state) : opt_->motionCost(motion->state, ppmotion->state);
                base::Cost nbhNewCost = opt_->combineCosts(ppmotion->cost, nbhIncCost);
                currentCost = opt_->combineCosts(currentCost, base::Cost(nbhNewCost.value() - cost.value()));
                base::Cost temp = opt_->combineCosts(currentStartCost_, currentGoalCost_);
                if (!opt_->isCostBetterThan(temp, bestCost_))
                    stop = true;
                else 
                    tvalid = isPathValidInter(ppmotion, start, stop);
                connectToPmotion(motion, ppmotion, start);
                motion->parent->children.push_back(motion); 
            }
            else 
            {
                pnullMotions.push_back(motion);
                motion->valid = ValidP;
                motion->pmotion = pmotion;
                motion->pmotion->pchildren.push_back(motion);
                nullMotions.push_back(motion);
                motion->cell->data->root++;
            }
            if (tvalid)
            {
                enableMotionInDisc(motion);
                checkedPath.resize(i+1);
                Motion *last = motion->parent;
                while (last != nullptr)
                {
                    checkedPath.push_back(last);
                    last = last->parent;
                }
            }
            else 
                break;
            if (stop)
            {
                tvalid = false;
                break;
            }
        }
    }
    for (auto & nullm : nullMotions)
    {
        for (auto & child : nullm->children)
        {
            child->valid = UnCkeckedP;
            child->parent = nullptr;
            child->pmotion = nullm;
            child->pmotion->pchildren.push_back(child);
            pnullMotions.push_back(child);
            child->cell->data->root++;
        }
        nullm->children.clear();
    }
    return tvalid;
}

bool ompl::geometric::BiHSCstar::isPathValidInter(Motion *motion, bool start, bool &stop)
{
    if (!isPathValidLazy(motion, start, stop))
        return false;
    if (stop)
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
    std::vector<Motion *> nullMotions;
    for (std::size_t i = mpath.size() - 1; i < mpath.size(); i--)
    {
        motion = mpath[i];
        Motion *pmotion = motion->parent;
        std::size_t spos = invalidMotions.size();
        if (!checkMotion(pmotion, motion, start))
        {
            removeFromParent(motion);
            motion->parent = nullptr;
            tvalid = false;
            std::size_t epos = invalidMotions.size();
            for (std::size_t i = spos; i < epos; i++)
            {
                Motion *temp = invalidMotions[i];
                for (auto & child : temp->children)
                {
                    child->pmotion = pmotion;
                    child->pmotion->pchildren.push_back(child);
                }
            }

            pnullMotions.push_back(motion);
            motion->valid = ValidP;
            motion->pmotion = pmotion;
            motion->pmotion->pchildren.push_back(motion);
            nullMotions.push_back(motion);
            motion->cell->data->root++;
            break;
        }
    }
    for (auto & nullm : nullMotions)
    {
        for (auto & child : nullm->children)
        {
            child->valid = UnCkeckedP;
            child->parent = nullptr;
            child->pmotion = nullm;
            child->pmotion->pchildren.push_back(child);
            pnullMotions.push_back(child);
            child->cell->data->root++;
        }
        nullm->children.clear();
    }
    return tvalid;
}

/*
bool ompl::geometric::BiHSCstar::isPathValidInter(Motion *motion, bool start, bool &stop) // back rewire
{
    if (!isPathValidLazy(motion, start, stop))
        return false;
    if (stop)
        return true;
    bool tvalid = true;
    std::vector<Motion *> mpath;
    while (motion->parent)
    {
        mpath.push_back(motion);
        motion = motion->parent;
    }
    base::Cost &currentCost = start ? currentStartCost_ : currentGoalCost_;
    std::vector<Motion *> &pnullMotions= start ? pnullStartMotions_: pnullGoalMotions_;
    std::vector<Motion *> &invalidMotions = start ? invalidStartMotions_ : invalidGoalMotions_;
    std::vector<Motion *> nullMotions;
    for (std::size_t i = mpath.size() - 1; i < mpath.size(); i--)
    {
        motion = mpath[i];
        Motion *pmotion = motion->parent;
        base::Cost cost = motion->cost;
        std::size_t spos = invalidMotions.size();
        if (!checkMotion(pmotion, motion, start))
        {
            removeFromParent(motion);
            motion->parent = nullptr;
            tvalid = false;
            Motion *ppmotion = nullptr;
            tvalid = backPathRewireMotion(motion, start, ppmotion);
            std::size_t epos = invalidMotions.size();
            for (std::size_t i = spos; i < epos; i++)
            {
                Motion *temp = invalidMotions[i];
                for (auto & child : temp->children)
                {
                    child->pmotion = pmotion;
                    child->pmotion->pchildren.push_back(child);
                }
            }
            if (tvalid)
            {
                base::Cost nbhIncCost = start ? opt_->motionCost(ppmotion->state, motion->state) : opt_->motionCost(motion->state, ppmotion->state);
                base::Cost nbhNewCost = opt_->combineCosts(ppmotion->cost, nbhIncCost);
                currentCost = opt_->combineCosts(currentCost, base::Cost(nbhNewCost.value() - cost.value()));
                base::Cost temp = opt_->combineCosts(currentStartCost_, currentGoalCost_);
                if (!opt_->isCostBetterThan(temp, bestCost_))
                    stop = true;
                else 
                    tvalid = isPathValidInter(ppmotion, start, stop);
                connectToPmotion(motion, ppmotion, start);
                motion->parent->children.push_back(motion); 
            }
            else 
            {
                pnullMotions.push_back(motion);
                motion->valid = ValidP;
                motion->pmotion = pmotion;
                motion->pmotion->pchildren.push_back(motion);
                nullMotions.push_back(motion);
                motion->cell->data->root++;
            }
            if (tvalid)
                enableMotionInDisc(motion);
            else 
                break;
            if (stop)
                break;
        }
    }
    for (auto & nullm : nullMotions)
    {
        for (auto & child : nullm->children)
        {
            child->valid = UnCkeckedP;
            child->parent = nullptr;
            child->pmotion = nullm;
            child->pmotion->pchildren.push_back(child);
            pnullMotions.push_back(child);
            child->cell->data->root++;
        }
        nullm->children.clear();
    }
    return tvalid;
}
*/

/*
bool ompl::geometric::BiHSCstar::isPathValidLazy(Motion *motion, bool start, bool &stop)
{
    if (lazyNode_ && !isStateValid(motion, start, stop))
        return false;
    if (stop)
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
    std::vector<Motion *> nullMotions;
    for (std::size_t i = mpath.size() - 1; i < mpath.size(); i--)
    {
        motion = mpath[i];
        Motion *pmotion = motion->parent;
        std::size_t spos = invalidMotions.size();
        if (!checkMotionLazy(pmotion, motion, start))
        {
            removeFromParent(motion);
            motion->parent = nullptr;
            tvalid = false;

            std::size_t epos = invalidMotions.size();
            for (std::size_t i = spos; i < epos; i++)
            {
                Motion *temp = invalidMotions[i];
                for (auto & child : temp->children)
                {
                    child->pmotion = pmotion;
                    child->pmotion->pchildren.push_back(child);
                }
            }

            pnullMotions.push_back(motion);
            motion->valid = ValidP;
            motion->pmotion = pmotion;
            motion->pmotion->pchildren.push_back(motion);
            nullMotions.push_back(motion);
            motion->cell->data->root++;
            break;
        }
    }
    for (auto & nullm : nullMotions)
    {
        for (auto & child : nullm->children)
        {
            child->valid = UnCkeckedP;
            child->parent = nullptr;
            child->pmotion = nullm;
            child->pmotion->pchildren.push_back(child);
            pnullMotions.push_back(child);
            child->cell->data->root++;
        }
        nullm->children.clear();
    }
    return tvalid;
}
*/

bool ompl::geometric::BiHSCstar::isPathValidLazy(Motion *motion, bool start, bool &stop) // back rewire
{
    if (lazyNode_ && !isStateValid(motion, start, stop))
        return false;
    if (stop)
        return true;
    /*
    bool tvalid = true;
    std::vector<Motion *> mpath;
    while (motion->parent)
    {
        mpath.push_back(motion);
        motion = motion->parent;
    }
    base::Cost &currentCost = start ? currentStartCost_ : currentGoalCost_;
    std::vector<Motion *> &pnullMotions= start ? pnullStartMotions_: pnullGoalMotions_;
    std::vector<Motion *> &invalidMotions = start ? invalidStartMotions_ : invalidGoalMotions_;
    std::vector<Motion *> nullMotions;
    for (std::size_t i = mpath.size() - 1; i < mpath.size(); i--)
    {
        motion = mpath[i];
        Motion *pmotion = motion->parent;
        base::Cost cost = motion->cost;
        std::size_t spos = invalidMotions.size();
        if (!checkMotionLazy(pmotion, motion, start))
        {
            removeFromParent(motion);
            motion->parent = nullptr;
            tvalid = false;
            Motion *ppmotion = nullptr;
            tvalid = backPathRewireMotionLazy(motion, start, ppmotion);
            if (motion->parent)
            {
                ppmotion = motion->parent;
                removeFromParent(motion);
                motion->parent = nullptr;
            }
            std::size_t epos = invalidMotions.size();
            for (std::size_t i = spos; i < epos; i++)
            {
                Motion *temp = invalidMotions[i];
                for (auto & child : temp->children)
                {
                    child->pmotion = pmotion;
                    child->pmotion->pchildren.push_back(child);
                }
            }
            if (tvalid)
            {
                base::Cost nbhIncCost = start ? opt_->motionCost(ppmotion->state, motion->state) : opt_->motionCost(motion->state, ppmotion->state);
                base::Cost nbhNewCost = opt_->combineCosts(ppmotion->cost, nbhIncCost);
                currentCost = opt_->combineCosts(currentCost, base::Cost(nbhNewCost.value() - cost.value()));
                base::Cost temp = opt_->combineCosts(currentStartCost_, currentGoalCost_);
                if (!opt_->isCostBetterThan(temp, bestCost_))
                    stop = true;
                else 
                    tvalid = isPathValidLazy(ppmotion, start, stop);
                connectToPmotion(motion, ppmotion, start);
                motion->parent->children.push_back(motion); 
            }
            else 
            {
                pnullMotions.push_back(motion);
                motion->valid = ValidP;
                motion->pmotion = pmotion;
                motion->pmotion->pchildren.push_back(motion);
                nullMotions.push_back(motion);
                motion->cell->data->root++;
            }
            if (tvalid)
                enableMotionInDisc(motion);
            else 
                break;
            if (stop)
                break;
        }
    }
    for (auto & nullm : nullMotions)
    {
        for (auto & child : nullm->children)
        {
            child->valid = UnCkeckedP;
            child->parent = nullptr;
            child->pmotion = nullm;
            child->pmotion->pchildren.push_back(child);
            pnullMotions.push_back(child);
            child->cell->data->root++;
        }
        nullm->children.clear();
    }
    return tvalid;
    */
    return true;
}

void ompl::geometric::BiHSCstar::removeInvalidMotions()
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

    for (std::size_t i = invalidStartNum_; i < invalidStartMotions_.size(); i++)
    {
        Motion *pnull = invalidStartMotions_[i];
        for (auto & child : pnull->children)
        {
            child->parent = nullptr;
            pnullStartMotions_.push_back(child);
            child->cell->data->root++;
        }
        pnull->children.clear();
    }
    invalidStartNum_ = invalidStartMotions_.size();

    for (std::size_t i = invalidGoalNum_; i < invalidGoalMotions_.size(); i++)
    {
        Motion *pnull = invalidGoalMotions_[i];
        for (auto & child : pnull->children)
        {
            child->parent = nullptr;
            pnullGoalMotions_.push_back(child);
            child->cell->data->root++;
        }
        pnull->children.clear();
    }
    invalidGoalNum_ = invalidGoalMotions_.size();
}

void ompl::geometric::BiHSCstar::removeInvalidMotionsTree(double ratio)
{
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

void ompl::geometric::BiHSCstar::addToTree(TreeData &tree, Motion *motion)
{
    tree->add(motion);
    for (auto & child : motion->children)
        addToTree(tree, child);
}

void ompl::geometric::BiHSCstar::addPdfMotion(MotionPDF &pdf, Motion *motion, bool start)
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

ompl::geometric::BiHSCstar::Motion *ompl::geometric::BiHSCstar::selectPdfMotion(MotionPDF &pdf, GridCell *&cell)
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

void ompl::geometric::BiHSCstar::removePdfMotion(MotionPDF &pdf, Motion *motion)
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

void ompl::geometric::BiHSCstar::enableMotionInDisc(Motion *motion)
{
    motion->cell->data->disabled--;
    if (!motion->cell->data->cmotion || opt_->isCostBetterThan(motion->cost, motion->cell->data->cmotion->cost))
        motion->cell->data->cmotion = motion;
    if (!motion->cell->data->mmotion || !opt_->isCostBetterThan(motion->cost, motion->cell->data->mmotion->cost))
        motion->cell->data->mmotion = motion;
    for (auto & child : motion->children)
        enableMotionInDisc(child);
}

void ompl::geometric::BiHSCstar::connectToPmotion(Motion *motion, Motion *pmotion, bool start) const
{
    motion->parent = pmotion;
    motion->root = pmotion->root;
    motion->incCost = start ? opt_->motionCost(pmotion->state, motion->state) : opt_->motionCost(motion->state, pmotion->state);
    motion->cost = opt_->combineCosts(pmotion->cost, motion->incCost);
    if (opt_->isFinite(motion->cost))
        updateChildCosts(motion);
}

void ompl::geometric::BiHSCstar::updateChildCosts(Motion *motion) const
{
    for (auto & child : motion->children)
    {
        child->cost = opt_->combineCosts(motion->cost, child->incCost);
        child->root = motion->root;
        updateChildCosts(child);
    }
}

void ompl::geometric::BiHSCstar::insertNeighbor(Motion *pmotion, Motion *motion)
{
    pmotion->nbh.insert(motion);
    motion->nbh.insert(pmotion);
}

void ompl::geometric::BiHSCstar::insertInvalidNeighbor(Motion *pmotion, Motion *motion)
{
    pmotion->invalidnbh.insert(motion);
    motion->invalidnbh.insert(pmotion);
}

void ompl::geometric::BiHSCstar::removeFromNeighbor(Motion *motion)
{
    for (auto & pmotion : motion->nbh)
        pmotion->nbh.erase(motion);
}

void ompl::geometric::BiHSCstar::removeFromInvalidNeighbor(Motion *motion)
{
    for (auto & pmotion : motion->invalidnbh)
        pmotion->invalidnbh.erase(motion);
}

bool ompl::geometric::BiHSCstar::isValidNeighbor(Motion *motion, Motion *nb) const
{
    return motion->nbh.find(nb) != motion->nbh.end();
}

bool ompl::geometric::BiHSCstar::isInvalidNeighbor(Motion *motion, Motion *nb) const
{
    return motion->invalidnbh.find(nb) != motion->invalidnbh.end();
}

void ompl::geometric::BiHSCstar::addToDisc(CellDiscretizationData &disc, Motion *motion, const Coord& coord)
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

void ompl::geometric::BiHSCstar::removeFromDisc(CellDiscretizationData &disc, Motion *motion)
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

void ompl::geometric::BiHSCstar::removeFromParent(Motion *motion)
{
    if (!motion->parent)
        return;
    removeFromVector(motion->parent->children, motion);
}

void ompl::geometric::BiHSCstar::removeFromPnull(std::vector<Motion *> &pnullMotions, Motion *motion)
{
    if (removeFromVector(pnullMotions, motion))
    {
        removeFromVector(motion->pmotion->pchildren, motion);
        motion->pmotion = nullptr;
        motion->cell->data->root--;
    }
}

bool ompl::geometric::BiHSCstar::removeFromVector(std::vector<Motion *> &motions, Motion *motion)
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

void ompl::geometric::BiHSCstar::setMotionInfinityCostWithDisable(Motion *motion, std::unordered_set<Cell *> &cells) const
{
    motion->cost = opt_->infiniteCost();
    motion->cell->data->disabled++;
    cells.insert(motion->cell);
    for (auto & child : motion->children)
        setMotionInfinityCostWithDisable(child, cells);
}

// check 
bool ompl::geometric::BiHSCstar::checkMotion(Motion *pmotion, Motion *motion, bool start)
{
    if (motion->valid != ValidP)
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
            motion->valid = ValidP;
            insertNeighbor(motion, pmotion);
            CellDiscretizationData &disc = start ? dStart_ : dGoal_;
            disc.updateNbh(motion->cell, pmotion->cell);
        }
    }
    return motion->valid == ValidP;
}

bool ompl::geometric::BiHSCstar::checkMotionLazy(Motion *pmotion, Motion *motion, bool start)
{
    if (motion->valid == UnCkeckedP && !checkInterMotionLazy(pmotion, motion, start))
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
    return motion->valid >= LazyValid;
}

bool ompl::geometric::BiHSCstar::checkInterMotion(Motion *pmotion, Motion *motion, bool start)
{
    bool valid = false;
    if (start || symmetric_)
        valid = checkInterMotion1(pmotion, motion, start);
    else 
        valid = checkInterMotion2(motion, pmotion, start);
    return valid;
}

bool ompl::geometric::BiHSCstar::checkInterMotionLazy(Motion *pmotion, Motion *motion, bool start)
{
    bool valid = true;
    if (pmotion->sce == motion->sce)
    {
        motion->valid = LazyValid;
        return valid; 
    }
    std::vector<double> scd1 = distanceCertificate_(pmotion->sce->sc->state, motion->state);
    if (!certificateOutside(scd1, pmotion->sce->sc->confidence_))
    {
        motion->valid = LazyValid;
        return valid; 
    }
    std::vector<double> scd2 = distanceCertificate_(motion->sce->sc->state, pmotion->state);
    if (!certificateOutside(scd2, motion->sce->sc->confidence_))
    {
        motion->valid = LazyValid;
        return valid; 
    }
    double low = 0.0, high = 1.0;
    base::State *state = si_->allocState();
    while (true)
    {
        double middle = 0.5 * (low + high);
        si_->getStateSpace()->interpolate(pmotion->state, motion->state, middle, state);
        std::vector<double> scd1 = distanceCertificate_(pmotion->sce->sc->state, state);
        std::vector<double> scd2 = distanceCertificate_(motion->sce->sc->state, state);
        if (!certificateOutside(scd1, pmotion->sce->sc->confidence_))
        {
            if (!certificateOutside(scd2, motion->sce->sc->confidence_))
            {
                if (isValid(state))
                {
                    Motion *last = new Motion(si_);
                    si_->copyState(last->state, state);
                    addIntermediateMotionLazy(pmotion, start, last);
                    removeFromParent(motion);
                    motion->parent = last;
                    motion->parent->children.push_back(motion);
                    motion->incCost = start ? opt_->motionCost(last->state, motion->state) : opt_->motionCost(motion->state, last->state);
                    motion->valid = LazyValid;
                }
                else 
                {
                    valid = false;
                    removeInvalidCertificateInter(pmotion->sce, scd1, start);
                    removeInvalidCertificateInter(motion->sce, scd2, start);
                }
                break;
            }
            else 
                low = middle;
        }
        else if (!certificateOutside(scd2, motion->sce->sc->confidence_))
            high = middle;
        else 
        {
            Motion *motion1 = nullptr, *motion2 = nullptr;
            if (low == 0.0)
                motion1 = pmotion;
            else 
            {
                motion1 = new Motion(si_);
                si_->getStateSpace()->interpolate(pmotion->state, motion->state, low, motion1->state);
                if (!isValid(motion1->state))
                {
                    valid = false;
                    std::vector<double> scd1 = distanceCertificate_(pmotion->sce->sc->state, motion1->state);
                    removeInvalidCertificateInter(pmotion->sce, scd1, start);
                    freeMotion(motion1);
                    break;
                }
                else
                    addIntermediateMotionLazy(pmotion, start, motion1);
            }
            if (high == 1.0)
                motion2 = motion;
            else 
            {
                motion2 = new Motion(si_);
                si_->getStateSpace()->interpolate(pmotion->state, motion->state, high, motion2->state);
                if (!isValid(motion2->state))
                {
                    valid = false;
                    std::vector<double> scd2 = distanceCertificate_(motion->sce->sc->state, motion2->state);
                    removeInvalidCertificateInter(motion->sce, scd2, start);
                    freeMotion(motion2);
                    break;
                }
            }
            if (checkInterMotion(motion1, motion2, start))
            {
                CellDiscretizationData &disc = start ? dStart_ : dGoal_;
                if (low == 0.0 && high == 1.0)
                {
                    motion->valid = ValidP;
                    insertNeighbor(motion, pmotion);
                    disc.updateNbh(motion->cell, pmotion->cell);
                }
                else if (low == 0.0)
                {
                    addIntermediateMotionLazy(pmotion, start, motion2);
                    motion2->valid = ValidP;
                    insertNeighbor(motion2, pmotion);
                    disc.updateNbh(motion2->cell, pmotion->cell);
                    removeFromParent(motion);
                    motion->parent = motion2;
                    motion->parent->children.push_back(motion);
                    motion->incCost = start ? opt_->motionCost(motion2->state, motion->state) : opt_->motionCost(motion->state, motion2->state);
                    motion->valid = LazyValid;
                }
                else if (high == 1.0)
                {
                    removeFromParent(motion);
                    motion->parent = motion1;
                    motion->parent->children.push_back(motion);
                    motion->incCost = start ? opt_->motionCost(motion1->state, motion->state) : opt_->motionCost(motion->state, motion1->state);
                    motion->valid = ValidP;
                    insertNeighbor(motion1, motion);
                    disc.updateNbh(motion1->cell, motion->cell);
                }
                else 
                {
                    addIntermediateMotionLazy(motion1, start, motion2);
                    removeFromParent(motion);
                    motion->parent = motion2;
                    motion->parent->children.push_back(motion);
                    motion->incCost = start ? opt_->motionCost(motion2->state, motion->state) : opt_->motionCost(motion->state, motion2->state);
                    motion->valid = LazyValid;
                    motion2->valid = ValidP;
                    insertNeighbor(motion1, motion2);
                    disc.updateNbh(motion1->cell, motion2->cell);
                }
            }
            else 
            {
                valid = false;
                if (high != 1.0)
                    freeMotion(motion2);
            }
            break;
        }
        if (high - low < 0.03)
        {
            motion->valid = LazyValid;
            break;
        }
    }
    si_->freeState(state);
    return valid;
}

void ompl::geometric::BiHSCstar::getPlannerData(base::PlannerData &data) const
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

// optimal
bool ompl::geometric::BiHSCstar::findBetterSolution(bool &optimal, double &ratio1, double &maxratio1, unsigned int &connect1, double &connectTresh1) // todo
{
    bool updatedSolution = false;
    for (auto & pair : connectionPoint_)
    {
        base::Cost temp = opt_->combineCosts(pair.first->cost, pair.second->cost);
        if (opt_->isFinite(temp) && opt_->isCostBetterThan(temp, bestCost_))
        {
            if (checkPath(temp, bestCost_, ratio1, maxratio1, connect1, connectTresh1))
            {
                if (isPathValid(pair.first, pair.second))
                {
                    temp = opt_->combineCosts(pair.first->cost, pair.second->cost);
                    if (opt_->isCostBetterThan(temp, bestCost_))
                    {
                        bestCost_ = temp;

                        double ratio = improvementRatio(temp, pair.first->root, pair.second->root);
                        ratio1 = std::min(ratio, ratio1);

                        if (solved_)
                            OMPL_INFORM("%s: Found a better solution with a cost of %.2f in %u iterations (%u "
                                    "vertices in the graph)",
                                    getName().c_str(), bestCost_.value(), iterations_, tStart_->size() + tGoal_->size());
                        else 
                            OMPL_INFORM("%s: Found an initial solution with a cost of %.2f in %u iterations (%u "
                                    "vertices in the graph)",
                                    getName().c_str(), bestCost_.value(), iterations_, tStart_->size() + tGoal_->size());

                        solved_ = true;
                        updatedSolution = true;
                        bestStartMotion_ = pair.first;
                        bestGoalMotion_ = pair.second;

                        if (opt_->isSatisfied(bestCost_))
                        {
                            optimal = true;
                            break;
                        }
                    }
                }
                rewirePath();
            }
        }
    }

    return updatedSolution;
}

void ompl::geometric::BiHSCstar::rewirePath()
{
    {
        bool updated = false;
        std::vector<Motion *> &checkedPath = checkedStartPath_;
        std::size_t i = checkedPath.size() - 1, j = i - 2, k = 0;
        while (i < checkedPath.size() && checkedPath[i]->valid && opt_->isFinite(checkedPath[i]->cost))
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
        if (updated)
        {
            removeFromParent(checkedPath[k]);
            connectToPmotion(checkedPath[k], checkedPath[i], true);
            checkedPath[k]->parent->children.push_back(checkedPath[k]);
            if (isValidNeighbor(checkedPath[k], checkedPath[i]))
                checkedPath[k]->valid = ValidP;
            else
                checkedPath[k]->valid = UnCkeckedP;
        }
    }
    {
        bool updated = false;
        std::vector<Motion *> &checkedPath = checkedGoalPath_;
        std::size_t i = checkedPath.size() - 1, j = i - 2, k = 0;
        while (i < checkedPath.size() && checkedPath[i]->valid && opt_->isFinite(checkedPath[i]->cost))
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
        if (updated)
        {
            removeFromParent(checkedPath[k]);
            connectToPmotion(checkedPath[k], checkedPath[i], false);
            checkedPath[k]->parent->children.push_back(checkedPath[k]);
            if (isValidNeighbor(checkedPath[k], checkedPath[i]))
                checkedPath[k]->valid = ValidP;
            else 
                checkedPath[k]->valid = UnCkeckedP;
        }
    }
}

void ompl::geometric::BiHSCstar::reportBetterSolution(const base::ReportIntermediateSolutionFn &intermediateSolutionCallback)
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

double ompl::geometric::BiHSCstar::improvementRatio(const base::Cost &temp, const base::State *sm, const base::State *gm) const
{
    base::Cost min = opt_->combineCosts(opt_->identityCost(), opt_->motionCost(sm, gm));
    double ratio = std::abs((temp.value() - min.value()) / temp.value());
    if (opt_->getCostThreshold().value() != 0.0)
        ratio = 0.7 * std::min(ratio, std::abs((temp.value() - opt_->getCostThreshold().value()) / temp.value()));
    else 
        ratio *= 0.5;
    return ratio;
}

bool ompl::geometric::BiHSCstar::checkPath(const base::Cost &temp, const base::Cost &best,
                                           double &ratio1, double &maxratio1, unsigned int &connect1, double &connectTresh1) const
{
    if (!solved_)
        return true;
    if (!lazyPath_)
        return true;
    bool check = false;
    double ratio = opt_->isFinite(best) ? std::abs((temp.value() - best.value()) / best.value()) : 2.0 * ratio1;
    if (ratio >= ratio1) 
    {
        connect1 = 0;
        maxratio1 = 0.0;
        ratio1 = 0.35 * (ratio1 + ratio);
        check = true;
    }
    else if (ratio > maxratio1)
    {
        maxratio1 = ratio;
        connect1++;
        if ((double)connect1 >= connectTresh1)
        {
            connect1 = 0;
            connectTresh1 *= 0.9;
            if (maxratio1 > 0.5 * ratio1)
            {
                ratio1 = 0.35 * (maxratio1 + ratio1);
                maxratio1 = 0.0;
                check = true;
            }
            else 
                ratio1 *= 0.9;
        }
    }
    return check;
}

void ompl::geometric::BiHSCstar::optimalRewireTree(BinaryHeap<Motion *, MotionCompare> &bh, Motion *m, bool start)
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
                if (motion->cell != c && c->data->mmotion)
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
                    if (isInvalidNeighbor(motion, nb))
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
                        if (feas)
                            nb->valid = ValidP;
                        else 
                            nb->valid = UnCkeckedP;
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

void ompl::geometric::BiHSCstar::updateLeafQueue(BinaryHeap<Motion *, MotionCompare> &bh, Motion *motion)
{
    if (motion->children.empty())
        updateQueue(bh, motion);
    else 
    {
        for (auto & child : motion->children)
            updateLeafQueue(bh, child);
    }
}

bool ompl::geometric::BiHSCstar::backRewire(Motion *motion, Motion *nb, bool start, Motion *&pmotion, base::Cost &nbhIncCost, base::Cost &nbhNewCost, bool &feas)
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
        nbhIncCost = start ? opt_->motionCost(motion->parent->state, nb->state) : opt_->motionCost(nb->state, motion->parent->state);
        nbhNewCost = opt_->combineCosts(motion->parent->cost, nbhIncCost);
        if (opt_->isCostBetterThan(nbhNewCost, nb->cost))
            temp = motion->parent;
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

std::size_t ompl::geometric::BiHSCstar::pruneTree(const base::Cost &pruneTreeCost)
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
        for (auto & pnull : invalidStartMotions_)
            freeMotion(pnull);
        invalidStartMotions_.clear();
        invalidStartNum_ = 0;
        for (auto & pnull : invalidGoalMotions_)
            freeMotion(pnull);
        invalidGoalMotions_.clear();
        invalidGoalNum_ = 0;

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
        for (auto & cell : cells)
        {
            if (cell->data->disabled)
                std::sort(cell->data->motions.begin(), cell->data->motions.end(), esort);
            if (!cell->data->cmotion)
                cell->data->cmotion = *std::min_element(cell->data->motions.begin(), cell->data->motions.end() - cell->data->disabled, mc_);
            if (!cell->data->mmotion)
                cell->data->mmotion = *std::max_element(cell->data->motions.begin(), cell->data->motions.end() - cell->data->disabled, mc_);
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
        for (auto & cell : cells)
        {
            if (cell->data->disabled)
                std::sort(cell->data->motions.begin(), cell->data->motions.end(), esort);
            if (!cell->data->cmotion)
                cell->data->cmotion = *std::min_element(cell->data->motions.begin(), cell->data->motions.end() - cell->data->disabled, mc_);
            if (!cell->data->mmotion)
                cell->data->mmotion = *std::max_element(cell->data->motions.begin(), cell->data->motions.end() - cell->data->disabled, mc_);
        }
        prunedCost_ = pruneTreeCost;
    }
    return numPruned;
}

std::size_t ompl::geometric::BiHSCstar::pruneSingleTree(TreeData &tree, CellDiscretizationData &disc, const base::Cost &pruneTreeCost, bool start,
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

std::size_t ompl::geometric::BiHSCstar::pruneTreeInternal(TreeData &tree, CellDiscretizationData &disc, const base::Cost &pruneTreeCost, bool start,
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

std::size_t ompl::geometric::BiHSCstar::pruneTreeInternalDisabled(TreeData &tree, CellDiscretizationData &disc, const base::Cost &pruneTreeCost, bool start,
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
                leavesToPrune.front()->cell->data->root--;
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

void ompl::geometric::BiHSCstar::toPrune(TreeData &tree, std::queue<Motion *, std::deque<Motion *>> &motionQueue,
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

void ompl::geometric::BiHSCstar::pruneMotion(Motion *motion, CellDiscretizationData &disc, bool start, std::unordered_set<Cell *> &cells)
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
    pruneMotionDisabled(motion, disc, start);
}

void ompl::geometric::BiHSCstar::pruneMotionDisabled(Motion *motion, CellDiscretizationData &disc, bool start)
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
    removeFromSafetyCerficate(motion);
    freeMotion(motion);
}

void ompl::geometric::BiHSCstar::addChildrenToList(std::queue<Motion *, std::deque<Motion *>> *motionList, Motion *motion)
{
    for (auto & child : motion->children)
        motionList->push(child);
}

bool ompl::geometric::BiHSCstar::keepCondition(Motion *motion, const base::Cost &threshold, bool start) // todo
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

ompl::base::Cost ompl::geometric::BiHSCstar::bordersolutionHeuristic(Motion *motion, bool start) const
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

bool ompl::geometric::BiHSCstar::keepCondition2(Motion *motion, const base::Cost &threshold, bool start) const
{
    return !opt_->isCostBetterThan(threshold, solutionHeuristic2(motion, start));
}

ompl::base::Cost ompl::geometric::BiHSCstar::solutionHeuristic2(Motion *motion, bool start) const
{
    base::Cost costToCome = calculateCostToCome(motion, start);
    base::Cost costToGo = calculateCostToGo(motion, start);
    return opt_->combineCosts(costToCome, costToGo);
}

ompl::base::Cost ompl::geometric::BiHSCstar::calculateCostToCome(Motion *motion, bool start) const
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

ompl::base::Cost ompl::geometric::BiHSCstar::calculateCostToGo(Motion *motion, bool start) const
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

/// safety certificate
bool ompl::geometric::BiHSCstar::isValid(base::SafetyCertificate *sc)
{
    if (useCollisionCertificateChecker_)
    {
        double dist = 0.0;
        if (!si_->isValid(sc->state, *sc->contact, dist))
        {
            CoordSC xcoord(projectionEvaluator_->getDimension());
            projectionEvaluator_->computeCoordinates(sc->state, xcoord);
            onn_->add(sc, xcoord);
            return false;
        }
    }
    else if (!si_->isValid(sc->state)) 
    {
        freeCertificate(sc);
        return false;
    }
    return true;
}

bool ompl::geometric::BiHSCstar::isValid(const base::State *state)
{
    time::point starto = time::now();
    base::SafetyCertificate *sc = new base::SafetyCertificate(si_, certificateDim_);
    si_->copyState(sc->state, state);
    bool valid = true;
    if (!isValid(sc))
        valid = false;
    else 
        freeCertificate(sc);
    oTime_ += time::seconds(time::now() - starto);
    return valid;
}

bool ompl::geometric::BiHSCstar::isValid(Motion *motion, bool start, bool add)
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
        time::point starto = time::now();
        motion->stateValid = InValid;
        if (add && addIntermediateState_ && motion->parent && motion->parent->stateValid == Valid)
        {
            if (checkInterMotion(motion->parent, motion, start))
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
                        addIntermediateMotion(motion->parent, motion, start, last);
                    }
                }
            }
        }
        
        if (!motion->parent)
        {
            if (start)
                removeFromPnull(pnullStartMotions_, motion);
            else 
                removeFromPnull(pnullGoalMotions_, motion);
        }
        else 
        {
            removeFromVector(motion->parent->children, motion);
            motion->parent = nullptr;
        }
        if (start)
            invalidStartMotions_.push_back(motion);
        else 
            invalidGoalMotions_.push_back(motion);
        CellDiscretizationData &disc = start ? dStart_ : dGoal_;
        if (opt_->isFinite(motion->cost))
        {
            std::unordered_set<Cell *> cells;
            setMotionInfinityCostWithDisable(motion, cells);
            if (motion->cell->data->motions.size() == 1)
                cells.erase(motion->cell);
            else 
                motion->cell->data->disabled--;
            removeFromDisc(disc, motion);
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
        else
        {
            motion->cell->data->disabled--;
            removeFromDisc(disc, motion);
        }
        oTime_ += time::seconds(time::now() - starto);
//        std::queue<SafetyCertificatePairs, std::deque<SafetyCertificatePairs>> scQueue;
        removeInvalidCertificate(motion->sce, motion->scd, start);
//        removeInvalidCertificate(scQueue, start);
        return false;
    }
    return motion->stateValid == Valid;
}

std::vector<ompl::geometric::BiHSCstar::Motion *> ompl::geometric::BiHSCstar::removeInvalidCertificate(SafetyCertificateWithElems *msce, const std::vector<double> &scd, bool start)
{
    std::vector<Motion *> invalid;
    if (!msce || certificateOutside(scd, msce->sc->confidence_))
        return invalid;
    return removeInvalidCertificateInter(msce, scd, start);
}

std::vector<ompl::geometric::BiHSCstar::Motion *> ompl::geometric::BiHSCstar::removeInvalidCertificateInter(SafetyCertificateWithElems *msce, const std::vector<double> &scd, bool start)
{
    std::transform(scd.begin(), scd.end(), msce->sc->confidence_.begin(), [](double a){return 0.9*a;});
    std::vector<Motion *> valid, invalid;
    for (auto it = msce->objects.begin(); it != msce->objects.end();)
    {
        Motion *scm = *it;
        if (certificateOutside(scm->scd, msce->sc->confidence_))
        {
            scm->sce = nullptr;
            if (isValid(scm, start))
                valid.push_back(scm);
            else 
                invalid.push_back(scm);
            std::iter_swap(it, msce->objects.end() - 1);
            msce->objects.pop_back();
            continue;
        }
        it++;
    }
    std::vector<SafetyCertificateWithElems *> &snne = start ? ssnne_ : gsnne_;
    for (auto &vmotion : valid)
    {
        base::SafetyCertificate *sc = new base::SafetyCertificate(si_, certificateDim_);
        si_->copyState(sc->state, vmotion->state);
        delete sc->contact;
        sc->contact = nullptr;

        SafetyCertificateWithElems *sce = new SafetyCertificateWithElems();
        sce->sc = sc;
        sce->objects.push_back(vmotion);
        snne.push_back(sce);

        vmotion->sce = sce;
        std::fill(vmotion->scd.begin(), vmotion->scd.end(), 0.0);
        vmotion->sce->sc->confidence_ = distanceCertificate_(msce->sc->state, vmotion->state);
        for (auto &invalidm : invalid)
        {
            std::vector<double> temp = distanceCertificate_(invalidm->state, vmotion->state);
            std::transform(temp.begin(), temp.end(), vmotion->sce->sc->confidence_.begin(), vmotion->sce->sc->confidence_.begin(), [](double a, double b){return std::min(0.9*a, b);} );
        }
    }
    return invalid;
}

void ompl::geometric::BiHSCstar::removeInvalidCertificate(SafetyCertificateWithElems *msce, const std::vector<double> &scd, bool start,
        std::queue<SafetyCertificatePairs, std::deque<SafetyCertificatePairs>> &scQueue)
{
    std::vector<Motion *> invalid = removeInvalidCertificate(msce, scd, start);
//    std::unordered_set<base::SafetyCertificate *> temscs;
    for (auto &invmotion : invalid) // todo
    {
        for (auto &child : invmotion->children)
        {
            if (child->sce && child->sce != msce)// && temscs.find(child->sce) == temscs.end())
            {
                std::vector<double> scd = distanceCertificate_(child->sce->sc->state, invmotion->state);
                if (!certificateOutside(scd, child->sce->sc->confidence_))
                {
                    SafetyCertificatePairs scp;
                    scp.sce = child->sce;
                    scp.scd = scd;
                    scQueue.push(scp);
                }
            }
        }
    }
}

void ompl::geometric::BiHSCstar::removeInvalidCertificate(std::queue<SafetyCertificatePairs, std::deque<SafetyCertificatePairs>> &scQueue, bool start)
{
    while (!scQueue.empty())
    {
        SafetyCertificatePairs scp = scQueue.front();
        scQueue.pop();
        removeInvalidCertificate(scp.sce, scp.scd, start);
    }
}

bool ompl::geometric::BiHSCstar::certificateOutside(const std::vector<double> &scd, const std::vector<double> &confidence) const
{
    bool outside = false;
    for (std::size_t i = 0; i < scd.size(); i++)
    {
        if (scd[i] > confidence[i])
        {
            outside = true;
            break;
        }
    }
    return outside;
}

bool ompl::geometric::BiHSCstar::checkInterMotion1(Motion *smotion, Motion *gmotion, bool start)
{
    /*assume smotion, gmotion are valid*/
    bool valid = true;
    const base::State *s1 = smotion->state, *s2 = gmotion->state;
    int nd = si_->getStateSpace()->validSegmentCount(s1, s2);
    if (nd >= 2)
    {
        base::SafetyCertificate *sc = new base::SafetyCertificate(si_, certificateDim_);
        base::State *st = nullptr;
        int i = 1;
        while (i < nd)
        {
            si_->getStateSpace()->interpolate(s1, s2, (double)i / (double)nd, sc->state);
            if (useCollisionCertificateChecker_)
            {
                CoordSC xcoord(projectionEvaluator_->getDimension());
                projectionEvaluator_->computeCoordinates(sc->state, xcoord);
                CellSC *cell = onn_->getCell(xcoord);
                if (cell && collisionCertificateChecker_(sc->state, cell->data->datas))
                {	
                    st = si_->allocState();
                    si_->copyState(st, sc->state);
                    valid = false;
                    freeCertificate(sc);
                    break;
                }
            }
            if (!isValid(sc))
            {
                st = si_->allocState();
                if (useCollisionCertificateChecker_)
                    si_->copyState(st, sc->state);
                else 
                    si_->getStateSpace()->interpolate(s1, s2, (double)i / (double)nd, st);
                valid = false;
                break;
            }
            i++;
        }
        if (!valid)
        {
            i--;
            double ratio = (double)i/(double)nd;
            if (ratio > minValidPathFraction_)
            {
                Motion *last = new Motion(si_);
                si_->getStateSpace()->interpolate(s1, s2, ratio, last->state);
                addIntermediateMotion(smotion, gmotion, start, last);
            }
            if (smotion->sce) // todo 
            {
                std::vector<double> scd = distanceCertificate_(smotion->sce->sc->state, st);
//                std::queue<SafetyCertificatePairs, std::deque<SafetyCertificatePairs>> scQueue;
                removeInvalidCertificate(smotion->sce, scd, start);
//                removeInvalidCertificate(scQueue, start);
            }
            if (gmotion->sce && gmotion->sce != smotion->sce)
            {
                std::vector<double> scd = distanceCertificate_(gmotion->sce->sc->state, st);
//                std::queue<SafetyCertificatePairs, std::deque<SafetyCertificatePairs>> scQueue;
                removeInvalidCertificate(gmotion->sce, scd, start);
//                removeInvalidCertificate(scQueue, start);
            }
            si_->freeState(st);
        }
        else
            freeCertificate(sc);
    }
    return valid;
}

bool ompl::geometric::BiHSCstar::checkInterMotion2(Motion *smotion, Motion *gmotion, bool start)
{
    /*assume smotion, gmotion are valid*/
    bool valid = true;
    const base::State *s1 = smotion->state, *s2 = gmotion->state;
    int nd = si_->getStateSpace()->validSegmentCount(s1, s2);
    if (nd >= 2)
    {
        base::SafetyCertificate *sc = new base::SafetyCertificate(si_, certificateDim_);
        base::State *st = nullptr;
        int i = nd - 1;
        while (i > 0)
        {
            si_->getStateSpace()->interpolate(s1, s2, (double)i / (double)nd, sc->state);
            if (useCollisionCertificateChecker_)
            {
                CoordSC xcoord(projectionEvaluator_->getDimension());
                projectionEvaluator_->computeCoordinates(sc->state, xcoord);
                CellSC *cell = onn_->getCell(xcoord);
                if (cell && collisionCertificateChecker_(sc->state, cell->data->datas))
                {	
                    st = si_->allocState();
                    si_->copyState(st, sc->state);
                    valid = false;
                    freeCertificate(sc);
                    break;
                }
            }
            if (!isValid(sc))
            {
                st = si_->allocState();
                if (useCollisionCertificateChecker_)
                    si_->copyState(st, sc->state);
                else 
                    si_->getStateSpace()->interpolate(s1, s2, (double)i / (double)nd, st);
                valid = false;
                break;
            }
            i--;
        }
        if (!valid)
        {
            i++;
            double ratio = (double)i/(double)nd;
            if (ratio < 1.0 - minValidPathFraction_)
            {
                Motion *last = new Motion(si_);
                si_->getStateSpace()->interpolate(s1, s2, ratio, last->state);
                addIntermediateMotion(gmotion, smotion, start, last);
            }
            if (smotion->sce) // todo 
            {
                std::vector<double> scd = distanceCertificate_(smotion->sce->sc->state, st);
//                std::queue<SafetyCertificatePairs, std::deque<SafetyCertificatePairs>> scQueue;
                removeInvalidCertificate(smotion->sce, scd, start);
//                removeInvalidCertificate(scQueue, start);
            }
            if (gmotion->sce && gmotion->sce != smotion->sce)
            {
                std::vector<double> scd = distanceCertificate_(gmotion->sce->sc->state, st);
//                std::queue<SafetyCertificatePairs, std::deque<SafetyCertificatePairs>> scQueue;
                removeInvalidCertificate(gmotion->sce, scd, start);
//                removeInvalidCertificate(scQueue, start);
            }
            si_->freeState(st);
        }
        else
            freeCertificate(sc);
    }
    return valid;
}

void ompl::geometric::BiHSCstar::removeFromSafetyCerficate(Motion *motion)
{
    if (motion->sce)
        removeFromVector(motion->sce->objects, motion);
}

void ompl::geometric::BiHSCstar::addIntermediateMotion(Motion *pmotion, Motion *motion, bool start, Motion *last)
{
    TreeData &tree = start ? tStart_ : tGoal_;
    last->valid = ValidP;
    last->stateValid = Valid;
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
    disc.updateNbh(last->cell, pmotion->cell);
    if (motion->stateValid == Valid)
        insertInvalidNeighbor(last, motion);
    if (!last->cell->data->cmotion || opt_->isCostBetterThan(last->cost, last->cell->data->cmotion->cost))
        last->cell->data->cmotion = last;
    if (!last->cell->data->mmotion || !opt_->isCostBetterThan(last->cost, last->cell->data->mmotion->cost))
        last->cell->data->mmotion = last;
    if (lazyNode_)
    {
        if (motion->sce)
            last->sce = motion->sce;
        else if (pmotion->sce)
            last->sce = pmotion->sce;
        last->scd = distanceCertificate_(last->sce->sc->state, last->state);
        last->sce->objects.push_back(last);
    }
}

void ompl::geometric::BiHSCstar::addIntermediateMotionLazy(Motion *pmotion, bool start, Motion *last)
{
    TreeData &tree = start ? tStart_ : tGoal_;
    last->valid = LazyValid;
    last->stateValid = Valid;
    connectToPmotion(last, pmotion, start);
    last->parent->children.push_back(last);
    tree->add(last);
    CellDiscretizationData &disc = start ? dStart_ : dGoal_;
    Coord xcoord(projectionEvaluator_->getDimension());
    projectionEvaluator_->computeCoordinates(last->state, xcoord);
    addToDisc(disc, last, xcoord);
    if (!opt_->isFinite(last->cost))
        last->cell->data->disabled++;
    if (!last->cell->data->cmotion || opt_->isCostBetterThan(last->cost, last->cell->data->cmotion->cost))
        last->cell->data->cmotion = last;
    if (!last->cell->data->mmotion || !opt_->isCostBetterThan(last->cost, last->cell->data->mmotion->cost))
        last->cell->data->mmotion = last;
    if (lazyNode_ && pmotion->sce)
    {
        last->sce = pmotion->sce;
        last->scd = distanceCertificate_(last->sce->sc->state, last->state);
        last->sce->objects.push_back(last);
    }
}

bool ompl::geometric::BiHSCstar::isStateValid(Motion *motion, bool start, bool &stop) // todo inter
{
    stop = false;
    bool tvalid = true;
    std::vector<Motion *> mpath;
    while (motion)
    {
        mpath.push_back(motion);
        motion = motion->parent;
    }
    base::Cost &currentCost = start ? currentStartCost_ : currentGoalCost_;
    std::vector<Motion *> &pnullMotions= start ? pnullStartMotions_: pnullGoalMotions_;
    std::vector<Motion *> &invalidMotions = start ? invalidStartMotions_ : invalidGoalMotions_;
    std::vector<Motion *> nullMotions;
    for (std::size_t i = mpath.size() - 2; i < mpath.size(); i--)
    {
        motion = mpath[i];
        Motion *pmotion = mpath[i+1];
        base::Cost cost = motion->cost;
        std::size_t epos = invalidMotions.size(), spos = epos;
        if (!isValid(motion, start, true))
        {
            tvalid = false;
            Motion *last = nullptr;
            i--;
            while (i < mpath.size())
            {
                cost = opt_->combineCosts(cost, mpath[i]->incCost);
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
                    base::Cost nbhIncCost = start ? opt_->motionCost(plast->state, last->state) : opt_->motionCost(last->state, plast->state);
                    base::Cost nbhNewCost = opt_->combineCosts(plast->cost, nbhIncCost);
                    currentCost = opt_->combineCosts(currentCost, base::Cost(nbhNewCost.value() - cost.value()));
                    base::Cost temp = opt_->combineCosts(currentStartCost_, currentGoalCost_);
                    if (!opt_->isCostBetterThan(temp, bestCost_))
                        stop = true;
                    else
                        tvalid = isPathValidInter(plast, start, stop); // todo
                    connectToPmotion(last, plast, start);
                    last->parent->children.push_back(last);
                    if (tvalid)
                        enableMotionInDisc(last);
                }
                else 
                {
                    last->valid = ValidP;
                    last->pmotion = pmotion;
                    last->pmotion->pchildren.push_back(last);
                    pnullMotions.push_back(last);
                    nullMotions.push_back(last);
                    last->cell->data->root++;
                }
            }
            else 
                epos = invalidMotions.size();
            for (std::size_t i = spos; i < epos; i++)
            {
                Motion *temp = invalidMotions[i];
                for (auto & child : temp->children)
                {
                    child->pmotion = pmotion;
                    child->pmotion->pchildren.push_back(child);
                }
            }
            if (!tvalid)
                break;
            if (stop)
                break;
        }
    }
    for (auto & nullm : nullMotions)
    {
        for (auto & child : nullm->children)
        {
            child->valid = UnCkeckedP;
            child->parent = nullptr;
            child->pmotion = nullm;
            child->pmotion->pchildren.push_back(child);
            pnullMotions.push_back(child);
            child->cell->data->root++;
        }
        nullm->children.clear();
    }
    return tvalid;
}

bool ompl::geometric::BiHSCstar::backPathRewireMotion(Motion *motion, bool start, Motion *&pmotion)
{
    bool valid = false;
    OrderCellsByCost ocbc(opt_);
    CostMotionCompare compareFn(motion, opt_, start);
    CellDiscretizationData &disc = start ? dStart_ : dGoal_;
    for (auto & cv : motion->cell->nbh)
    {
        std::size_t numc = cv.size(), ic = numc - 1;
        while (ic < numc)
        {
            if (ic)
                std::sort(cv.begin(), cv.begin() + ic + 1, ocbc);
            Cell *c = cv[ic];
            ic--;
            std::size_t num = c->data->motions.size() - c->data->disabled;
            if (!num)
                break;
            if (motion->cell != c && motion->cell->data->mmotion)
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
            std::unordered_set<Cell *> ctemp;
            std::size_t ictemp = ic;
            while (ictemp < numc)
            {
                ctemp.insert(cv[ictemp]);
                ictemp--;
            }
            bool invalid = false;
            for (auto & nb : nbh)
            {
                if (!opt_->isFinite(nb->cost))
                    continue;
                if (!isValid(nb, start))
                {
                    invalid = true;
                    continue;
                }
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
                motion->valid = ValidP;
                valid = true;
                break;
            }
            if (valid)
                break;
            if (invalid && !ctemp.empty())
            {
                numc = cv.size(), ic = numc - 1;
                while (ic < numc)
                {
                    if (ctemp.find(cv[ic]) != ctemp.end())
                        break;
                    ic--;
                }
            }
        }
        if (valid)
            break;
    }

    return valid;
}

bool ompl::geometric::BiHSCstar::backPathRewireMotionLazy(Motion *motion, bool start, Motion *&pmotion)
{
    bool valid = false;
    OrderCellsByCost ocbc(opt_);
    CostMotionCompare compareFn(motion, opt_, start);
    for (auto & cv : motion->cell->nbh)
    {
        std::size_t numc = cv.size(), ic = numc - 1;
        while (ic < numc)
        {
            if (ic)
                std::sort(cv.begin(), cv.begin() + ic + 1, ocbc);
            Cell *c = cv[ic];
            ic--;
            std::size_t num = c->data->motions.size() - c->data->disabled;
            if (!num)
                break;
            if (motion->cell != c && motion->cell->data->mmotion)
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
            std::unordered_set<Cell *> ctemp;
            std::size_t ictemp = ic;
            while (ictemp < numc)
            {
                ctemp.insert(cv[ictemp]);
                ictemp--;
            }
            bool invalid = false;
            for (auto & nb : nbh)
            {
                if (!opt_->isFinite(nb->cost))
                    continue;
                if (!isValid(nb, start))
                {
                    invalid = true;
                    continue;
                }
                if (isInvalidNeighbor(motion, nb))
                    continue;
                if (!isValidNeighbor(motion, nb))
                {
                    if (!checkInterMotionLazy(nb, motion, start))
                    {
                        insertInvalidNeighbor(nb, motion);
                        continue;
                    }
                }
                else 
                    motion->valid = ValidP;
                pmotion = nb;
                valid = true;
                break;
            }
            if (valid)
                break;
            if (invalid && !ctemp.empty())
            {
                numc = cv.size(), ic = numc - 1;
                while (ic < numc)
                {
                    if (ctemp.find(cv[ic]) != ctemp.end())
                        break;
                    ic--;
                }
            }
        }
        if (valid)
            break;
    }

    return valid;
}
