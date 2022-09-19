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

#include "ompl/geometric/planners/hsc/BiHSCCellstar.h"
#include "ompl/base/goals/GoalSampleableRegion.h"
#include "ompl/base/objectives/PathLengthOptimizationObjective.h"
#include "ompl/tools/config/SelfConfig.h"

ompl::geometric::BiHSCCellstar::BiHSCCellstar(const base::SpaceInformationPtr &si) : base::Planner(si, "BiHSCCellstar")
  ,mc_(opt_), bh_(mc_)
{
    specs_.recognizedGoal = base::GOAL_SAMPLEABLE_REGION;
    specs_.directed = true;
    specs_.approximateSolutions = false;
    specs_.optimizingPaths = true;
    specs_.canReportIntermediateSolutions = true;

    Planner::declareParam<double>("range", this, &BiHSCCellstar::setRange, &BiHSCCellstar::getRange, "0.:0.1:10000.");
    Planner::declareParam<double>("pen_distance", this, &BiHSCCellstar::setPenDistance, &BiHSCCellstar::getPenDistance, "0.:0.1:10000.");
    Planner::declareParam<double>("min_valid_path_fraction", this, &BiHSCCellstar::setMinValidPathFraction, &BiHSCCellstar::getMinValidPathFraction, "0.:0.05:1.");
    Planner::declareParam<bool>("lazy_path", this, &BiHSCCellstar::setLazyPath, &BiHSCCellstar::getLazyPath, "0,1");
    Planner::declareParam<bool>("lazy_node", this, &BiHSCCellstar::setLazyNode, &BiHSCCellstar::getLazyNode, "0,1");
    Planner::declareParam<bool>("add_intermediate_state", this, &BiHSCCellstar::setAddIntermediateState, &BiHSCCellstar::getAddIntermediateState, "0,1");
    Planner::declareParam<bool>("use_bispace", this, &BiHSCCellstar::setUseBispace, &BiHSCCellstar::getUseBispace, "0,1");
//    Planner::declareParam<bool>("treated_as_multi_subapce", this, &BiHSCCellstar::setTreatedAsMultiSubapce, &BiHSCCellstar::getTreatedAsMultiSubapce, "0,1");
    Planner::declareParam<double>("prune_threshold", this, &BiHSCCellstar::setPruneThreshold, &BiHSCCellstar::getPruneThreshold, "0.:.01:1.");

    addPlannerProgressProperty("iterations INTEGER", [this] { return numIterationsProperty(); });
    addPlannerProgressProperty("best cost REAL", [this] { return bestCostProperty(); });

    Planner::declareParam<bool>("strice_certificate", this, &BiHSCCellstar::setStrictCertificate, &BiHSCCellstar::getStrictCertificate, "0,1");
}

ompl::geometric::BiHSCCellstar::~BiHSCCellstar()
{
    freeMemory();
}

void ompl::geometric::BiHSCCellstar::setup()
{
    Planner::setup();
    tools::SelfConfig sc(si_, getName());
    sc.configurePlannerRange(maxDistance_);
    sc.configurePenetrationDistance(penDistance_);
    if (treatedAsMultiSubapce_)
    {
        maxDistance_ /= si_->getStateSpace()->getSubspaceCount();
        penDistance_ /= si_->getStateSpace()->getSubspaceCount();
    }

    if (minValidPathFraction_ < std::numeric_limits<double>::epsilon() || minValidPathFraction_ > 1.0)
        throw Exception("The minimum valid path fraction must be in the range (0,1]");

    sc.configureProjectionEvaluator(projectionEvaluator_);

    const std::vector<double> &csizes = projectionEvaluator_->getCellSizes();
    double maxSize = std::accumulate(csizes.begin(), csizes.end(), 0.0, [](double a, double b){return a + b*b;});
    maxSize = std::sqrt(maxSize);
    if (maxDistance_ < maxSize)
    {
        OMPL_WARN("%s: The range is smaller than the cell extent (%.4f VS %.4f)", getName().c_str(), maxDistance_, maxSize);
        maxDistance_ = 2.0 * maxSize;
    }

    dStart_.setDimension(projectionEvaluator_->getDimension());
    dGoal_.setDimension(projectionEvaluator_->getDimension());
    dStart_.setFreeMotionFn([this](Motion *m) { freeMotion(m); });
    dGoal_.setFreeMotionFn([this](Motion *m) { freeMotion(m); });
    dStart_.setNeighborCell(2);
    dGoal_.setNeighborCell(2);

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
        bestCost_ = opt_->infiniteCost();
        prunedCost_ = opt_->infiniteCost();
        mc_ = MotionCompare(opt_);
        bh_ = BinaryHeap<Motion *, MotionCompare>(mc_);
        dStart_.setIsEnabledMotionFn([this](Motion *m) { return opt_->isFinite(m->cost); });
        dGoal_.setIsEnabledMotionFn([this](Motion *m) { return opt_->isFinite(m->cost); });
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
        if (!onn_)
            onn_.reset(tools::SelfConfig::getDefaultNearestNeighbors<base::SafetyCertificate *>(this));
        onn_->setDistanceFunction([this](const base::SafetyCertificate *a, const base::SafetyCertificate *b)
                {
                return distanceFunction(a, b);
                });
    }
    const base::StateSpacePtr &space = si_->getStateSpace();
    const base::PlannerSpecs &specs = this->getSpecs();
    dStart_.getGrid().useThread(specs.multithreaded);
    dStart_.getGrid().metricSpace(space->isMetricSpace());
    dStart_.getGrid().setDistanceFunction([this](const Motion *a, const Motion *b){ return distanceFunction(a, b); });
    dGoal_.getGrid().useThread(specs.multithreaded);
    dGoal_.getGrid().metricSpace(space->isMetricSpace());
    dGoal_.getGrid().setDistanceFunction([this](const Motion *a, const Motion *b){ return distanceFunction(a, b); });
}

void ompl::geometric::BiHSCCellstar::freeMemory()
{
    if (lazyNode_)
    {
        for (auto & sc : snne_)
        {
            freeCertificate(sc->sc);
            delete sc;
        }
        snne_.clear();
    }

    if (onn_)
    {
        std::vector<base::SafetyCertificate *> safetycertificates;
        onn_->list(safetycertificates);
        for (auto & sc : safetycertificates)
            freeCertificate(sc);
        safetycertificates.clear();
    }
}

void ompl::geometric::BiHSCCellstar::clear()
{
    setup_ = false;

    Planner::clear();
    sampler_.reset();

    freeMemory();

    dStart_.clear();
    dGoal_.clear();

    startMotions_.clear();
    goalMotions_.clear();

    pnullStartMotions_.clear();
    pnullGoalMotions_.clear();

    checkedStartPath_.clear();
    checkedGoalPath_.clear();

    invalidStartMotions_.clear();
    invalidGoalMotions_.clear();

    connectionPoint_.clear();

    bestStartMotion_ = nullptr;
    bestGoalMotion_ = nullptr;

    solved_ = false;
    iterations_ = 0;
    bestCost_ = base::Cost(std::numeric_limits<double>::quiet_NaN());
    prunedCost_ = base::Cost(std::numeric_limits<double>::quiet_NaN());

    if (onn_)
    {
        onn_->clear();
    }
    snne_.clear();
}

ompl::base::PlannerStatus ompl::geometric::BiHSCCellstar::solve(const base::PlannerTerminationCondition &ptc)
{
    base::PlannerStatus ps= prepareSolve(ptc);
    if (ps != base::PlannerStatus::PREPARE_SUCCESS)
        return ps;   

    sampler_ = si_->allocStateSampler();

    OMPL_INFORM("%s: Starting planning with %u states already in datastructure. Seeking a solution better than %.5f.",
            getName().c_str(), (dStart_.getMotionCount() + dGoal_.getMotionCount()), opt_->getCostThreshold().value());

    if (!si_->getStateSpace()->isMetricSpace())
    {
        OMPL_WARN("%s: The state space (%s) is not metric and as a result the optimization objective may not satisfy "
                  "the triangle inequality. "
                  "You may need to disable pruning or rejection.",
                  getName().c_str(), si_->getStateSpace()->getName().c_str());
    }

    const base::ReportIntermediateSolutionFn intermediateSolutionCallback = pdef_->getIntermediateSolutionCallback();

    bool startTree = true;
    bool optimal = false;

    unsigned int connect1 = 0;
    double ratio1 = 0.5, maxratio1 = 0.0, connectTresh1 = 10.0;

    while (!ptc)
    {
        iterations_++;

        if (pis_.getSampledGoalsCount() < dGoal_.getMotionCount() / 2)
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
                Coord xcoord(projectionEvaluator_->getDimension());
                projectionEvaluator_->computeCoordinates(motion->state, xcoord);
                dGoal_.addMotion(motion, xcoord);
                dGoal_.addRootCell(motion->cell);
                motion->cell->auxData->cmotion = motion;
                motion->cell->auxData->mmotion = motion;
                if (motion->cell->auxData->disabled) 
                    dGoal_.updateCell(motion->cell);
                if (lazyNode_)
                {
                    base::SafetyCertificate *sc = new base::SafetyCertificate(si_, certificateDim_);
                    si_->copyState(sc->state, st);
                    delete sc->contact;
                    sc->contact = nullptr;

                    SafetyCertificateWithElems *sce = new SafetyCertificateWithElems();
                    sce->sc = sc;
                    sce->objects.push_back(motion);
                    snne_.push_back(sce);

                    motion->sce = sce;
                    motion->scd.resize(certificateDim_);
                    std::fill(motion->scd.begin(), motion->scd.end(), 0.0);
                }
            }
        }

        batchGrow(startTree);

        bool updatedSolution = findBetterSolution(optimal, ratio1, maxratio1, connect1, connectTresh1);
        if (lazyNode_)
            removeInvalidMotions();
        if (optimal)
            break;
        if (updatedSolution)
        {
            reportBetterSolution(intermediateSolutionCallback);
            int numPruned = pruneTree(bestCost_);
            if (false)
                OMPL_INFORM("%s: %u states are pruned from the tree, %u states are left", getName().c_str(), numPruned, dStart_.getMotionCount() + dGoal_.getMotionCount());
        }
    }

    if (solved_ || optimal)
    {
        ptc.terminate();
        if (!optimal)
        {
            isPathValid(bestStartMotion_, bestGoalMotion_);
            if (lazyNode_)
                removeInvalidMotions();
        }
        processSolution(bestStartMotion_, bestGoalMotion_);
    }

    OMPL_INFORM("%s: Created %u (%u start + %u goal) states in %u cells (%u start + %u goal). Final solution cost %.5f.",
                getName().c_str(), dStart_.getMotionCount() + dGoal_.getMotionCount(), dStart_.getMotionCount(),
                dGoal_.getMotionCount(), dStart_.getCellCount() + dGoal_.getCellCount(), dStart_.getCellCount(),
                dGoal_.getCellCount(), bestCost_.value());
    return (solved_ || optimal) ? base::PlannerStatus::EXACT_SOLUTION : base::PlannerStatus::TIMEOUT;
}

ompl::base::PlannerStatus ompl::geometric::BiHSCCellstar::prepareSolve(const base::PlannerTerminationCondition &ptc)
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
            Coord xcoord(projectionEvaluator_->getDimension());
            projectionEvaluator_->computeCoordinates(motion->state, xcoord);
            dStart_.addMotion(motion, xcoord);
            dStart_.addRootCell(motion->cell);
            motion->cell->auxData->cmotion = motion;
            motion->cell->auxData->mmotion = motion;
            if (lazyNode_)
            {
                base::SafetyCertificate *sc = new base::SafetyCertificate(si_, certificateDim_);
                si_->copyState(sc->state, st);
                delete sc->contact;
                sc->contact = nullptr;

                SafetyCertificateWithElems *sce = new SafetyCertificateWithElems();
                sce->sc = sc;
                sce->objects.push_back(motion);
                snne_.push_back(sce);

                motion->sce = sce;
                motion->scd.resize(certificateDim_);
                std::fill(motion->scd.begin(), motion->scd.end(), 0.0);
            }
        }
    }

    if (dStart_.getMotionCount() == 0)
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
        Coord xcoord(projectionEvaluator_->getDimension());
        projectionEvaluator_->computeCoordinates(motion->state, xcoord);
        dGoal_.addMotion(motion, xcoord);
        dGoal_.addRootCell(motion->cell);
        motion->cell->auxData->cmotion = motion;
        motion->cell->auxData->mmotion = motion;
        if (lazyNode_)
        {
            base::SafetyCertificate *sc = new base::SafetyCertificate(si_, certificateDim_);
            si_->copyState(sc->state, st);
            delete sc->contact;
            sc->contact = nullptr;

            SafetyCertificateWithElems *sce = new SafetyCertificateWithElems();
            sce->sc = sc;
            sce->objects.push_back(motion);
            snne_.push_back(sce);

            motion->sce = sce;
            motion->scd.resize(certificateDim_);
            std::fill(motion->scd.begin(), motion->scd.end(), 0.0);
        }
    }

    if (dGoal_.getMotionCount() == 0)
    {
        OMPL_ERROR("%s: Unable to sample any valid states for goal tree", getName().c_str());
        return base::PlannerStatus::INVALID_GOAL;
    }

    double dist = -1.0;
    for (auto & sm : startMotions_)
    {
        for (auto & gm : goalMotions_)
        {
            double d = si_->distance(sm->state, gm->state);
            if (dist < d)
            {
                dist = d;
            }
        }
    }

    penDistance_ = std::min(penDistance_, 0.1*dist);

    return base::PlannerStatus::PREPARE_SUCCESS;
}

void ompl::geometric::BiHSCCellstar::processSolution(const Motion *bestStartMotion, const Motion *bestGoalMotion)
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

    // Does the solution satisfy the optimization objective?
    psol.setOptimized(opt_, bestCost_, opt_->isSatisfied(bestCost_));
    pdef_->addSolutionPath(psol);
}

bool ompl::geometric::BiHSCCellstar::batchGrow(bool &startTree)
{
    bool nconnect = false;
    TreeGrowingInfo tgi;
    tgi.xstate = si_->allocState();
    for (unsigned int i = 0; i < 1; i++)
    {
        tgi.start = startTree;
        startTree = !startTree;
        Motion *startMotion = nullptr, *goalMotion = nullptr;
        growTree(tgi, startMotion, goalMotion);
        if (startMotion && goalMotion && pdef_->getGoal()->isStartGoalPairValid(startMotion->root, goalMotion->root))
        {
            nconnect = true;
            startMotion->inConnection = true;
            goalMotion->inConnection = true;
            connectionPoint_.emplace_back(startMotion, goalMotion);
        }
    }
    si_->freeState(tgi.xstate);
    return nconnect;
}

void ompl::geometric::BiHSCCellstar::growTree(TreeGrowingInfo &tgi, Motion *&startMotion, Motion *&goalMotion)
{
    CellDiscretizationData &disc = tgi.start ? dStart_ : dGoal_;
    Cell *ncell = nullptr;
    Motion *nmotion = nullptr;
    disc.selectMotion(nmotion, ncell);
    disc.updateCell(ncell);
    if (!opt_->isFinite(nmotion->cost))
        return;

    disc.countIteration();

    double maxDistance = maxDistance_;
    sampler_->sampleUniformNear(tgi.xstate, nmotion->state, maxDistance);
    if (useBispace_)
    {
        bool currentTree = (tgi.start == growCurrentTree(tgi.xstate));
        if (!currentTree) 
        {
            Motion *last = nmotion;
            while (last && tgi.start != growCurrentTree(last->state))
                last = last->parent;
            if (!last || penetrationDistance(last->state, tgi.xstate, tgi.start) > penDistance_)
                return;
        }
    }

    Motion *motion = new Motion(si_);
    si_->copyState(motion->state, tgi.xstate);
    Coord xcoord(projectionEvaluator_->getDimension());
    projectionEvaluator_->computeCoordinates(motion->state, xcoord);
    Motion *scm = lazyNode_ ? disc.getGrid().nearest(motion, xcoord) : nullptr;

    base::SafetyCertificate *sc = new base::SafetyCertificate(si_, certificateDim_);
    si_->copyState(sc->state, tgi.xstate);

    bool lazy = false;
    bool cvalid = true;
    std::vector<double> cdist;
    if (useCollisionCertificateChecker_ && onn_->size() > oscNum_)
    {
        std::vector<base::SafetyCertificate *> nsc;
        onn_->nearestK(sc, 5, nsc);
        if (collisionCertificateChecker_(sc->state, nsc))
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
        else if (scm->sce != nullptr && safetyCertificateChecker_(sc->state, scm->sce->sc, cdist))
            lazy = true;
        else if (!isValid(sc))
            cvalid = false;
    }
    if (!cvalid)
    {
        freeMotion(motion);
        return;
    }
    delete sc->contact;
    sc->contact = nullptr;

    if (!lazy)
        motion->stateValid = Valid;
    disc.addMotion(motion, xcoord);

    Motion *nb = nmotion;
    connectToPmotion(motion, nb, tgi.start);

    CostMotionCompare compareFn(motion, opt_, tgi.start);
    for (auto & cv : motion->cell->nbh)
    {
        for (auto & c : cv)
        {
            if (c->data.size() > c->auxData->disabled)
            {
                Motion *temp = *std::min_element(c->data.begin(), c->data.end(), compareFn);
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

    if (!lazyPath_)
    {
        if (tgi.start ? !checkInterMotion(nb, motion, true) : !checkInterMotion(motion, nb, false))
        {
            disc.removeMotion(motion);
            removeFromInvalidNeighbor(motion);
            freeMotion(motion);
            freeCertificate(sc);
            return;
        }
    }

    if (solved_ && !keepCondition2(motion, bestCost_, tgi.start))
    {
        disc.removeMotion(motion);
        freeMotion(motion);
        freeCertificate(sc);
        return;
    }
    
    if (lazyNode_)
    {
        if (!lazy)
        {
            SafetyCertificateWithElems *sce = new SafetyCertificateWithElems();
            sce->sc = sc;
            snne_.push_back(sce);

            motion->sce = sce;
            motion->sce->sc->confidence_ = distanceCertificate_(scm->state, motion->state);
            motion->sce->objects.push_back(motion);
            motion->scd.resize(certificateDim_);
            std::fill(motion->scd.begin(), motion->scd.end(), 0.0);

            if (strictCertificate_ && useCollisionCertificateChecker_ && onn_->size() > oscNum_)
            {
                base::SafetyCertificate *nsc = onn_->nearest(sc);
                std::vector<double> temp = distanceCertificate_(nsc->state, motion->state);
                std::transform(temp.begin(), temp.end(), motion->sce->sc->confidence_.begin(), motion->sce->sc->confidence_.begin(), [](double a, double b){return std::min(0.9*a, b);} );
            }
        }
        else 
        {
            freeCertificate(sc);
            motion->sce = scm->sce;
            motion->scd = cdist;
            motion->sce->objects.push_back(motion);
        }
    }
    else 
        freeCertificate(sc);

    if (motion->cell->auxData->disabled)
        disc.updateCell(motion->cell);
    motion->parent->children.push_back(motion);

    if (!lazyPath_)
    {
        motion->valid = true;
        insertNeighbor(nb, motion);
        disc.getGrid().updateNbh(nb->cell, motion->cell);
    }

    if (!motion->cell->auxData->cmotion || opt_->isCostBetterThan(motion->cost, motion->cell->auxData->cmotion->cost))
        motion->cell->auxData->cmotion = motion;
    if (!motion->cell->auxData->mmotion || !opt_->isCostBetterThan(motion->cost, motion->cell->auxData->mmotion->cost))
        motion->cell->auxData->mmotion = motion;
    if (!solved_)
    {
        if (tgi.start ? !pnullStartMotions_.empty() : !pnullGoalMotions_.empty())
            rewireTree(bh_, motion, tgi.start);
    }
    else 
        optimalRewireTree(bh_, motion, tgi.start);

    CellDiscretizationData &otherDisc = !tgi.start ? dStart_ : dGoal_;
    Cell *ocell = otherDisc.getGrid().getCell(xcoord);
    bool nconnect = false;
    Motion *connectOther = nullptr;
    std::size_t inConnection = 0;
    while (!nconnect && ocell && !ocell->data.empty())
    {
        std::size_t osize = ocell->data.size() - inConnection;
        std::size_t index = rng_.uniformInt(0, osize - 1);
        connectOther = ocell->data[index];
        if (!connectOther->inConnection)
        {
            if (true)
                nconnect = true;
            else if (osize == 1)
                break;
        }
        else if (osize == 1) 
            break;
        else 
        {
            inConnection++;
            std::iter_swap(ocell->data.begin() + index, ocell->data.end() - inConnection);
        }
    }

    if (nconnect)
    {
        Motion *connect = new Motion(si_);
        si_->copyState(connect->state, connectOther->state);
        connectToPmotion(connect, motion, tgi.start);
        connect->parent->children.push_back(connect);
        disc.addMotion(connect, xcoord);
        connect->stateValid = connectOther->stateValid;
        if (connect->cell->auxData->disabled)
            disc.updateCell(connect->cell);
        if (lazyNode_)
        {
            connect->sce = motion->sce;
            connect->scd = distanceCertificate_(connect->sce->sc->state, connect->state);
            connect->sce->objects.push_back(connect);
        }
        if (!connect->cell->auxData->cmotion || opt_->isCostBetterThan(connect->cost, connect->cell->auxData->cmotion->cost))
            connect->cell->auxData->cmotion = connect;
        if (!connect->cell->auxData->mmotion || !opt_->isCostBetterThan(connect->cost, connect->cell->auxData->mmotion->cost))
            connect->cell->auxData->mmotion = connect;
        if (!solved_)
        {
            if (tgi.start ? !pnullStartMotions_.empty() : !pnullGoalMotions_.empty())
                rewireTree(bh_, connect, tgi.start);
        }
        else
            optimalRewireTree(bh_, connect, tgi.start);

        if (tgi.start)
        {
            startMotion = connect;
            goalMotion = connectOther;
        }
        else 
        {
            startMotion = connectOther;
            goalMotion = connect;
        }
    }
}

void ompl::geometric::BiHSCCellstar::rewireTree(BinaryHeap<Motion *, MotionCompare> &bh, Motion *m, bool start)
{
    CellDiscretizationData &disc = start ? dStart_ : dGoal_;
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
                std::size_t pos = c->auxData->disabled;
                if (!pos)
                    break;
                if (motion->cell != c && c->auxData->mmotion)
                {
                    if (opt_->isCostBetterThan(c->auxData->mmotion->cost, motion->cell->auxData->cmotion->cost))
                        continue;
                    base::Cost cost1 = opt_->combineCosts(motion->cost, opt_->motionCost(motion->state, c->auxData->cmotion->state));
                    base::Cost cost2 = opt_->combineCosts(c->auxData->mmotion->cost, opt_->combineCosts(c->auxData->mmotion->cost, c->auxData->mmotion->cost));
                    if (opt_->isCostBetterThan(cost2, cost1))
                        continue;
                }
                std::vector<Motion *> nbh;
                nbh.reserve(pos);
                if (c->data.size() == c->auxData->disabled)
                {
                    std::sort(c->data.begin(), c->data.end(), rsort);               
                    nbh.insert(nbh.end(), c->data.begin(), c->data.end());
                }
                else 
                {
                    for (std::size_t i = 0; i < c->data.size() && nbh.size() < pos; i++)
                    {
                        if (!opt_->isFinite(c->data[i]->cost))
                            nbh.push_back(c->data[i]);
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
                    enableMotionInDisc(disc, nb);
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

void ompl::geometric::BiHSCCellstar::updateQueue(BinaryHeap<Motion *, MotionCompare> &bh, Motion *motion)
{
    if (motion->handle != nullptr)
        bh.update(motion->handle);
    else
        motion->handle = bh.insert(motion);
}

void ompl::geometric::BiHSCCellstar::updateLeafQueue(BinaryHeap<Motion *, MotionCompare> &bh, Motion *motion)
{
    if (motion->children.empty())
        updateQueue(bh, motion);
    else 
    {
        for (auto & child : motion->children)
            updateLeafQueue(bh, child);
    }
}

bool ompl::geometric::BiHSCCellstar::growCurrentTree(const base::State *state) const
{
    return growStartTree(state);
}

bool ompl::geometric::BiHSCCellstar::growStartTree(const base::State *state) const
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

double ompl::geometric::BiHSCCellstar::penetrationDistance(const base::State *nstate, const base::State *state, bool start) const
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

bool ompl::geometric::BiHSCCellstar::growCurrentTree(const base::State *state, std::size_t sub) const
{
    return growStartTree(state, sub);
}

bool ompl::geometric::BiHSCCellstar::growStartTree(const base::State *state, std::size_t sub) const
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

double ompl::geometric::BiHSCCellstar::penetrationDistance(const base::State *nstate, const base::State *state, bool start, std::size_t sub) const
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

bool ompl::geometric::BiHSCCellstar::isPathValid(Motion *motion, Motion *otherMotion)
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

bool ompl::geometric::BiHSCCellstar::isPathValid(Motion *motion, bool start) // todo
{
    if (!lazyPath_)
        return true;
    if (lazyNode_ && !isStateValid(motion, start))
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
    std::vector<Motion *> &pnullMotions= start ? pnullStartMotions_: pnullGoalMotions_;
    CellDiscretizationData &disc = start ? dStart_ : dGoal_;
    std::vector<Motion *> nullMotions;
    for (std::size_t i = mpath.size() - 1; i < mpath.size(); i--)
    {
        motion = mpath[i];
        if (start ? !checkStartMotion(motion->parent, motion) : !checkGoalMotion(motion, motion->parent))
        {
            removeFromParent(motion);
            motion->parent = nullptr;
            tvalid = false;
            Motion *pmotion = nullptr;
            if (backPathRewireMotion(motion, start, pmotion))
            {
                checkedPath.resize(i+1);
                Motion *last = pmotion;
                while (last)
                {
                    checkedPath.push_back(last);
                    last = last->parent;
                }
                tvalid = isPathValidInter(pmotion, start);
                connectToPmotion(motion, pmotion, start);
                motion->parent->children.push_back(motion); 
            }
            else 
            {
                motion->valid = true;
                pnullMotions.push_back(motion);
                nullMotions.push_back(motion);
                motion->cell->auxData->root++;
            }
            if (tvalid)
                enableMotionInDisc(disc, motion);
            else
                break;
        }
    }
    for (auto & nullm : nullMotions)
    {
        for (auto & child : nullm->children)
        {
            child->valid = false;
            child->parent = nullptr;
            pnullMotions.push_back(child);
            child->cell->auxData->root++;
        }
        nullm->children.clear();
    }
    return tvalid;
}

bool ompl::geometric::BiHSCCellstar::isStateValid(Motion *motion, bool start) // todo
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
    CellDiscretizationData &disc = start ? dStart_ : dGoal_;
    std::vector<Motion *> nullMotions;
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
            if (last != nullptr)
            {
                removeFromParent(last);
                last->parent = nullptr;
                Motion *plast = nullptr;
                if (backPathRewireMotion(last, start, plast))
                {
                    tvalid = isPathValidInter(plast, start);
                    connectToPmotion(last, plast, start);
                    last->parent->children.push_back(last);
                }
                else 
                {
                    last->valid = true;
                    pnullMotions.push_back(last);
                    nullMotions.push_back(last);
                    last->cell->auxData->root++;
                }
            }
            if (tvalid)
                enableMotionInDisc(disc, last);
            else
                break;
        }
    }
    for (auto & nullm : nullMotions)
    {
        for (auto & child : nullm->children)
        {
            child->valid = false;
            child->parent = nullptr;
            pnullMotions.push_back(child);
            child->cell->auxData->root++;
        }
        nullm->children.clear();
    }
    return tvalid;
}

bool ompl::geometric::BiHSCCellstar::isPathValidInter(Motion *motion, bool start) // todo
{
    if (lazyNode_ && !isStateValid(motion, start))
        return false;
    bool tvalid = true;
    std::vector<Motion *> &pnullMotions= start ? pnullStartMotions_: pnullGoalMotions_;
    while (motion->parent)
    {
        if (start ? !checkStartMotion(motion->parent, motion) : !checkGoalMotion(motion, motion->parent))
        {
            removeFromParent(motion);
            motion->parent = nullptr;
            tvalid = false;

            pnullMotions.push_back(motion);
            motion->valid = true;
            motion->cell->auxData->root++;
            for (auto & child : motion->children)
            {
                child->valid = false;
                child->parent = nullptr;
                pnullMotions.push_back(child);
                child->cell->auxData->root++;
            }
            motion->children.clear();
            break;
        }
        motion = motion->parent;
    }
    return tvalid;
}

/*
bool ompl::geometric::BiHSCCellstar::isPathValidInter(Motion *motion, bool start) // back rewire
{
    if (lazyNode_ && !isStateValid(motion, start))
        return false;
    bool tvalid = true;
    std::vector<Motion *> &pnullMotions= start ? pnullStartMotions_: pnullGoalMotions_;
    CellDiscretizationData &disc = start ? dStart_ : dGoal_;
    while (motion->parent)
    {
        if (start ? !checkStartMotion(motion->parent, motion) : !checkGoalMotion(motion, motion->parent))
        {
            removeFromParent(motion);
            motion->parent = nullptr;
            tvalid = false;
            Motion *pmotion = nullptr;
            if (backPathRewireMotion(motion, start, pmotion))
            {
                tvalid = isPathValidInter(pmotion, start);
                connectToPmotion(motion, pmotion, start);
                motion->parent->children.push_back(motion);
                if (tvalid)
                    enableMotionInDisc(disc, motion);
            }
            else 
            {
                pnullMotions.push_back(motion);
                motion->valid = true;
                motion->cell->auxData->root++;
                for (auto & child : motion->children)
                {
                    child->valid = false;
                    child->parent = nullptr;
                    pnullMotions.push_back(child);
                    child->cell->auxData->root++;
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

void ompl::geometric::BiHSCCellstar::removeInvalidMotions()
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
    removeInvalidMotionsDisc();
}

void ompl::geometric::BiHSCCellstar::removeInvalidMotionsDisc()
{
    for (auto & pnull : invalidStartMotions_)
    {
        for (auto & child : pnull->children)
        {
            child->parent = nullptr;
            pnullStartMotions_.push_back(child);
            child->cell->auxData->root++;
        }
        freeMotion(pnull);
    }
    invalidStartMotions_.clear();

    for (auto & pnull : invalidGoalMotions_)
    {
        for (auto & child : pnull->children)
        {
            child->parent = nullptr;
            pnullGoalMotions_.push_back(child);
            child->cell->auxData->root++;
        }
        freeMotion(pnull);
    }
    invalidGoalMotions_.clear();
}

void ompl::geometric::BiHSCCellstar::enableMotionInDisc(CellDiscretizationData &disc, Motion *motion)
{
    std::unordered_set<Cell *> cells;
    enableMotionInDisc2(motion, cells);
    for (auto & cell : cells)
    {
        if (cell->auxData->score <= std::numeric_limits<double>::epsilon())
        {
            cell->auxData->score = 1.0 + log((double)(cell->auxData->iteration));
            disc.updateCell(cell);
        }
    }
}

void ompl::geometric::BiHSCCellstar::enableMotionInDisc2(Motion *motion, std::unordered_set<Cell *> &cells)
{
    motion->cell->auxData->disabled--;
    cells.insert(motion->cell);
    if (!motion->cell->auxData->cmotion || opt_->isCostBetterThan(motion->cost, motion->cell->auxData->cmotion->cost))
        motion->cell->auxData->cmotion = motion;
    if (!motion->cell->auxData->mmotion || !opt_->isCostBetterThan(motion->cost, motion->cell->auxData->mmotion->cost))
        motion->cell->auxData->mmotion = motion;
    for (auto & child : motion->children)
        enableMotionInDisc2(child, cells);
}

void ompl::geometric::BiHSCCellstar::connectToPmotion(Motion *motion, Motion *pmotion, bool start) const
{
    motion->parent = pmotion;
    motion->root = pmotion->root;
    motion->incCost = start ? opt_->motionCost(pmotion->state, motion->state) : opt_->motionCost(motion->state, pmotion->state);
    motion->cost = opt_->combineCosts(pmotion->cost, motion->incCost);
    if (opt_->isFinite(motion->cost))
        updateChildCosts(motion);
}

void ompl::geometric::BiHSCCellstar::updateChildCosts(Motion *motion) const
{
    for (auto & child : motion->children)
    {
        child->cost = opt_->combineCosts(motion->cost, child->incCost);
        child->root = motion->root;
        updateChildCosts(child);
    }
}

void ompl::geometric::BiHSCCellstar::insertNeighbor(Motion *pmotion, Motion *motion)
{
    pmotion->nbh.insert(motion);
    motion->nbh.insert(pmotion);
}

void ompl::geometric::BiHSCCellstar::insertInvalidNeighbor(Motion *pmotion, Motion *motion)
{
    pmotion->invalidnbh.insert(motion);
    motion->invalidnbh.insert(pmotion);
}

void ompl::geometric::BiHSCCellstar::removeFromNeighbor(Motion *motion)
{
    for (auto & pmotion : motion->nbh)
        pmotion->nbh.erase(motion);
}

void ompl::geometric::BiHSCCellstar::removeFromInvalidNeighbor(Motion *motion)
{
    for (auto & pmotion : motion->invalidnbh)
        pmotion->invalidnbh.erase(motion);
}

bool ompl::geometric::BiHSCCellstar::isValidNeighbor(Motion *motion, Motion *nb) const
{
    return motion->nbh.find(nb) != motion->nbh.end();
}

bool ompl::geometric::BiHSCCellstar::isInvalidNeighbor(Motion *motion, Motion *nb) const
{
    return motion->invalidnbh.find(nb) != motion->invalidnbh.end();
}

void ompl::geometric::BiHSCCellstar::removeFromParent(Motion *motion)
{
    if (!motion->parent)
        return;
    removeFromVector(motion->parent->children, motion);
}

void ompl::geometric::BiHSCCellstar::removeFromPnull(std::vector<Motion *> &pnullMotions, Motion *motion)
{
    if (removeFromVector(pnullMotions, motion))
        motion->cell->auxData->root--;
}

bool ompl::geometric::BiHSCCellstar::removeFromVector(std::vector<Motion *> &motions, Motion *motion)
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

void ompl::geometric::BiHSCCellstar::setMotionInfinityCostWithDisable(Motion *motion, std::unordered_set<Cell *> &cells) const
{
    motion->cost = opt_->infiniteCost();
    motion->cell->auxData->disabled++;
    cells.insert(motion->cell);
    for (auto & child : motion->children)
        setMotionInfinityCostWithDisable(child, cells);
}

// check 
bool ompl::geometric::BiHSCCellstar::checkStartMotion(Motion *smotion, Motion *gmotion)
{
    if (!gmotion->valid)
    {
        if (!checkInterMotion(smotion, gmotion, true))
        {
            if (opt_->isFinite(gmotion->cost))
            {
                std::unordered_set<Cell *> cells;
                setMotionInfinityCostWithDisable(gmotion, cells);
                for (auto & cell : cells)
                {
                    dStart_.updateCell(cell);
                    if (cell->data.size() == cell->auxData->disabled)
                    {
                        cell->auxData->cmotion = nullptr;
                        cell->auxData->mmotion = nullptr;
                    }
                    else
                    {
                        bool e1 = !cell->auxData->cmotion || !opt_->isFinite(cell->auxData->cmotion->cost);
                        bool e2 = !cell->auxData->mmotion || !opt_->isFinite(cell->auxData->mmotion->cost);
                        if (e1 || e2)
                        {
                            dStart_.enableSort(cell);
                            if (e1)
                                cell->auxData->cmotion = *std::min_element(cell->data.begin(), cell->data.end() - cell->auxData->disabled, mc_);
                            if (e2)
                                cell->auxData->mmotion = *std::max_element(cell->data.begin(), cell->data.end() - cell->auxData->disabled, mc_);
                        }
                    }
                }
            }
            insertInvalidNeighbor(smotion, gmotion);
        }
        else 
        {
            gmotion->valid = true;
            insertNeighbor(gmotion, smotion);
            dStart_.getGrid().updateNbh(gmotion->cell, smotion->cell);
        }
    }
    return gmotion->valid;
}

bool ompl::geometric::BiHSCCellstar::checkGoalMotion(Motion *smotion, Motion *gmotion)
{
    if (!smotion->valid)
    {
        if (!checkInterMotion(smotion, gmotion, false))
        {
            if (opt_->isFinite(smotion->cost))
            {
                std::unordered_set<Cell *> cells;
                setMotionInfinityCostWithDisable(smotion, cells);
                for (auto & cell : cells)
                {
                    dGoal_.updateCell(cell);
                    if (cell->data.size() == cell->auxData->disabled)
                    {
                        cell->auxData->cmotion = nullptr;
                        cell->auxData->mmotion = nullptr;
                    }
                    else 
                    {
                        bool e1 = !cell->auxData->cmotion || !opt_->isFinite(cell->auxData->cmotion->cost);
                        bool e2 = !cell->auxData->mmotion || !opt_->isFinite(cell->auxData->mmotion->cost);
                        if (e1 || e2)
                        {
                            dGoal_.enableSort(cell);
                            if (e1)
                                cell->auxData->cmotion = *std::min_element(cell->data.begin(), cell->data.end() - cell->auxData->disabled, mc_);
                            if (e2)
                                cell->auxData->mmotion = *std::max_element(cell->data.begin(), cell->data.end() - cell->auxData->disabled, mc_);
                        }
                    }
                }
            }
            insertInvalidNeighbor(smotion, gmotion);
        }
        else 
        {
            smotion->valid = true;
            insertNeighbor(gmotion, smotion);
            dGoal_.getGrid().updateNbh(gmotion->cell, smotion->cell);
        }
    }
    return smotion->valid;
}

bool ompl::geometric::BiHSCCellstar::checkInterMotion(Motion *smotion, Motion *gmotion, bool start)
{
    bool valid = false;
    if (addIntermediateState_)
        valid = checkInterMotion2(smotion, gmotion, start);
    else
        valid = checkInterMotion1(smotion, gmotion, start);
    return valid;
}

void ompl::geometric::BiHSCCellstar::getPlannerData(base::PlannerData &data) const
{
    Planner::getPlannerData(data);
    for (auto & motion : startMotions_)
    {
        data.addStartVertex(base::PlannerDataVertex(motion->state, 1));
        for (auto & child : motion->children)
            getPlannerData(data, child, 1);
    }

    for (auto & motion : goalMotions_)
    {
        data.addGoalVertex(base::PlannerDataVertex(motion->state, 2));
        for (auto & child : motion->children)
            getPlannerData(data, child, 2);
    }
}

void ompl::geometric::BiHSCCellstar::getPlannerData(base::PlannerData &data, Motion *motion, int tag) const
{
    data.addEdge(base::PlannerDataVertex(motion->parent->state, tag), base::PlannerDataVertex(motion->state, tag));
    for (auto & child : motion->children)
    {
        getPlannerData(data, child, tag);
    }
}

// optimal
bool ompl::geometric::BiHSCCellstar::findBetterSolution(bool &optimal, double &ratio1, double &maxratio1, unsigned int &connect1, double &connectTresh1)
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
                                    getName().c_str(), bestCost_.value(), iterations_, dStart_.getMotionCount() + dGoal_.getMotionCount());
                        else 
                            OMPL_INFORM("%s: Found an initial solution with a cost of %.2f in %u iterations (%u "
                                    "vertices in the graph)",
                                    getName().c_str(), bestCost_.value(), iterations_, dStart_.getMotionCount() + dGoal_.getMotionCount());

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

void ompl::geometric::BiHSCCellstar::rewirePath()
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
            checkedPath[k]->valid = isValidNeighbor(checkedPath[k], checkedPath[i]);
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
            checkedPath[k]->valid = isValidNeighbor(checkedPath[k], checkedPath[i]);
        }
    }
}

void ompl::geometric::BiHSCCellstar::reportBetterSolution(const base::ReportIntermediateSolutionFn &intermediateSolutionCallback)
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

double ompl::geometric::BiHSCCellstar::improvementRatio(const base::Cost &temp, const base::State *sm, const base::State *gm) const
{
    base::Cost min = opt_->combineCosts(opt_->identityCost(), opt_->motionCost(sm, gm));
    double ratio = std::abs((temp.value() - min.value()) / temp.value());
    if (opt_->getCostThreshold().value() != 0.0)
        ratio = 0.7 * std::min(ratio, std::abs((temp.value() - opt_->getCostThreshold().value()) / temp.value()));
    else 
        ratio *= 0.5;
    return ratio;
}

bool ompl::geometric::BiHSCCellstar::checkPath(const base::Cost &temp, const base::Cost &best,
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

void ompl::geometric::BiHSCCellstar::optimalRewireTree(BinaryHeap<Motion *, MotionCompare> &bh, Motion *m, bool start)
{
    CellDiscretizationData &disc = start ? dStart_ : dGoal_;
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
                if (motion->cell != c && c->auxData->mmotion)
                {
                    if (opt_->isCostBetterThan(c->auxData->mmotion->cost, motion->cell->auxData->cmotion->cost))
                        continue;
                    base::Cost cost1 = opt_->combineCosts(motion->cost, opt_->motionCost(motion->state, c->auxData->cmotion->state));
                    base::Cost cost2 = opt_->combineCosts(c->auxData->mmotion->cost, opt_->combineCosts(c->auxData->mmotion->cost, c->auxData->mmotion->cost));
                    if (opt_->isCostBetterThan(cost2, cost1))
                        continue;
                }
                for (auto & nb : c->data)
                {
                    if (isInvalidNeighbor(motion, nb))
                        continue;
                    Motion *pmotion = nullptr;
                    base::Cost nbhIncCost, nbhNewCost;
                    bool feas = false;
                    if (backRewire(motion, nb, start, pmotion, nbhIncCost, nbhNewCost, feas))
                    {
                        if (!opt_->isFinite(nb->cost))
                            enableMotionInDisc(disc, nb);
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

bool ompl::geometric::BiHSCCellstar::backRewire(Motion *motion, Motion *nb, bool start, Motion *&pmotion, base::Cost &nbhIncCost, base::Cost &nbhNewCost, bool &feas)
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

std::size_t ompl::geometric::BiHSCCellstar::pruneTree(const base::Cost &pruneTreeCost)
{
    double fracBetter;
    std::size_t numPruned = 0;
    if (opt_->isFinite(prunedCost_))
        fracBetter = std::abs((pruneTreeCost.value() - prunedCost_.value()) / prunedCost_.value());
    else
        fracBetter = 1.0;
    if (fracBetter > pruneThreshold_)
    {
        std::unordered_set<Cell *> cells;
        numPruned += pruneSingleTree(dStart_, pruneTreeCost, true, startMotions_, cells);
        bool invalid = false;
        std::queue<Motion *, std::deque<Motion *>> motionQueue;
        for (auto it = pnullStartMotions_.begin(); it != pnullStartMotions_.end();)
        {
            Motion *motion = *it;
            motionQueue.push(motion);
            numPruned += pruneTreeInternalDisabled(dStart_, pruneTreeCost, true, motionQueue, invalid, cells);
            if (invalid)
            {
                std::iter_swap(it, pnullStartMotions_.end() - 1);
                pnullStartMotions_.pop_back();
            }
            else 
                it++;
        }
        for (auto & cell : cells)
        {
            if (!cell->removed)
            {
                dStart_.updateCell(cell);
                if ((!cell->auxData->cmotion || !cell->auxData->mmotion) && cell->data.size() > cell->auxData->disabled)
                {
                    if (cell->auxData->disabled)
                        dStart_.enableSort(cell);
                    if (!cell->auxData->cmotion)
                        cell->auxData->cmotion = *std::min_element(cell->data.begin(), cell->data.end() - cell->auxData->disabled, mc_);
                    if (!cell->auxData->mmotion)
                        cell->auxData->mmotion = *std::max_element(cell->data.begin(), cell->data.end() - cell->auxData->disabled, mc_);
                }
            }
        }

        cells.clear();
        numPruned += pruneSingleTree(dGoal_, pruneTreeCost, false, goalMotions_, cells);
        for (auto it = pnullGoalMotions_.begin(); it != pnullGoalMotions_.end();)
        {
            Motion *motion = *it;
            motionQueue.push(motion);
            numPruned += pruneTreeInternalDisabled(dGoal_, pruneTreeCost, false, motionQueue, invalid, cells);
            if (invalid)
            {
                std::iter_swap(it, pnullGoalMotions_.end() - 1);
                pnullGoalMotions_.pop_back();
            }
            else 
                it++;
        }
        for (auto & cell : cells)
        {
            if (!cell->removed)
            {
                dGoal_.updateCell(cell);
                if ((!cell->auxData->cmotion || !cell->auxData->mmotion) && cell->data.size() > cell->auxData->disabled)
                {
                    if (cell->auxData->disabled)
                        dGoal_.enableSort(cell);
                    if (!cell->auxData->cmotion)
                        cell->auxData->cmotion = *std::min_element(cell->data.begin(), cell->data.end() - cell->auxData->disabled, mc_);
                    if (!cell->auxData->mmotion)
                        cell->auxData->mmotion = *std::max_element(cell->data.begin(), cell->data.end() - cell->auxData->disabled, mc_);
                }
            }
        }

        prunedCost_ = pruneTreeCost;
    }
    return numPruned;
}

std::size_t ompl::geometric::BiHSCCellstar::pruneSingleTree(CellDiscretizationData &disc, const base::Cost &pruneTreeCost, bool start,
                                                  const std::vector<Motion *> &rootMotions, std::unordered_set<Cell *> &cells)
{
    std::size_t numPruned = 0;
    std::queue<Motion *, std::deque<Motion *>> motionQueue;
    for (auto & rootMotion : rootMotions)
        addChildrenToList(&motionQueue, rootMotion);
    numPruned += pruneTreeInternal(disc, pruneTreeCost, start, motionQueue, cells);
    return numPruned;
}

std::size_t ompl::geometric::BiHSCCellstar::pruneTreeInternal(CellDiscretizationData &disc, const base::Cost &pruneTreeCost, bool start,
        std::queue<Motion *, std::deque<Motion *>> &motionQueue, std::unordered_set<Cell *> &cells)
{
    std::size_t numPruned = 0;
    std::queue<Motion *, std::deque<Motion *>> leavesToPrune;
    std::list<Motion *> chainsToRecheck;
    toPrune(motionQueue, leavesToPrune, chainsToRecheck, pruneTreeCost, start);
    while (!leavesToPrune.empty())
    {
        while (!leavesToPrune.empty())
        {
            if (leavesToPrune.front()->cell->auxData->cmotion == leavesToPrune.front())
                leavesToPrune.front()->cell->auxData->cmotion = nullptr;
            if (leavesToPrune.front()->cell->auxData->mmotion == leavesToPrune.front())
                leavesToPrune.front()->cell->auxData->mmotion = nullptr;
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

    return numPruned;
}

std::size_t ompl::geometric::BiHSCCellstar::pruneTreeInternalDisabled(CellDiscretizationData &disc, const base::Cost &pruneTreeCost, bool start,
        std::queue<Motion *, std::deque<Motion *>> &motionQueue, bool &invalid, std::unordered_set<Cell *> &cells)
{
    invalid = false;

    std::size_t numPruned = 0;
    std::queue<Motion *, std::deque<Motion *>> leavesToPrune;
    std::list<Motion *> chainsToRecheck;
    toPrune(motionQueue, leavesToPrune, chainsToRecheck, pruneTreeCost, start);
    while (!leavesToPrune.empty())
    {
        while (!leavesToPrune.empty())
        {
            if (leavesToPrune.front()->parent == nullptr)
            {
                invalid = true;
                leavesToPrune.front()->cell->auxData->root--;
            }
            leavesToPrune.front()->cell->auxData->disabled--;
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

    return numPruned;
}

void ompl::geometric::BiHSCCellstar::toPrune(std::queue<Motion *, std::deque<Motion *>> &motionQueue,
        std::queue<Motion *, std::deque<Motion *>> &leavesToPrune, std::list<Motion *> &chainsToRecheck, const base::Cost &pruneTreeCost, bool start)
{
    while (!motionQueue.empty())
    {
        if (keepCondition(motionQueue.front(), pruneTreeCost, start))
            addChildrenToList(&motionQueue, motionQueue.front());
        else
        {
            if (!motionQueue.front()->children.empty())
            {
                bool keepAChild = false;
                for (unsigned int i = 0u; keepAChild == false && i < motionQueue.front()->children.size(); ++i)
                    keepAChild = keepCondition(motionQueue.front()->children.at(i), pruneTreeCost, start);
                if (!keepAChild)
                    chainsToRecheck.push_back(motionQueue.front());
                addChildrenToList(&motionQueue, motionQueue.front());
            }
            else
                leavesToPrune.push(motionQueue.front());
        }
        motionQueue.pop();
    }
}

void ompl::geometric::BiHSCCellstar::pruneMotion(Motion *motion, CellDiscretizationData &disc, bool start, std::unordered_set<Cell *> &cells)
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
    cells.insert(motion->cell);
    disc.removeMotion(motion);
    removeFromSafetyCerficate(motion);
    freeMotion(motion);
}

void ompl::geometric::BiHSCCellstar::addChildrenToList(std::queue<Motion *, std::deque<Motion *>> *motionList, Motion *motion)
{
    for (auto & child : motion->children)
        motionList->push(child);
}

bool ompl::geometric::BiHSCCellstar::keepCondition(Motion *motion, const base::Cost &threshold, bool start) // todo
{
    if (start && motion == bestStartMotion_)
        return true;
    if (!start && motion == bestGoalMotion_)
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

ompl::base::Cost ompl::geometric::BiHSCCellstar::bordersolutionHeuristic(Motion *motion, bool start) const
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

bool ompl::geometric::BiHSCCellstar::keepCondition2(Motion *motion, const base::Cost &threshold, bool start) const
{
    return !opt_->isCostBetterThan(threshold, solutionHeuristic2(motion, start));
}

ompl::base::Cost ompl::geometric::BiHSCCellstar::solutionHeuristic2(Motion *motion, bool start) const
{
    base::Cost costToCome = calculateCostToCome(motion, start);
    base::Cost costToGo = calculateCostToGo(motion, start);
    return opt_->combineCosts(costToCome, costToGo);
}

ompl::base::Cost ompl::geometric::BiHSCCellstar::calculateCostToCome(Motion *motion, bool start) const
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

ompl::base::Cost ompl::geometric::BiHSCCellstar::calculateCostToGo(Motion *motion, bool start) const
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

// safety certificate 
bool ompl::geometric::BiHSCCellstar::isValid(base::SafetyCertificate *sc)
{
    if (useCollisionCertificateChecker_)
    {
        double dist = 0.0;
        if (!si_->isValid(sc->state, *sc->contact, dist))
        {
            onn_->add(sc);
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

bool ompl::geometric::BiHSCCellstar::isValid(const base::State *state)
{
    base::SafetyCertificate *sc = new base::SafetyCertificate(si_, certificateDim_);
    si_->copyState(sc->state, state);
    bool valid = true;
    if (!isValid(sc))
        valid = false;
    else 
        freeCertificate(sc);
    return valid;
}

bool ompl::geometric::BiHSCCellstar::isValid(Motion *motion, bool start, bool add)
{
    if (!lazyNode_)
        return true;
    if (motion->stateValid == UnCkecked)
    {
        bool valid = true;
        bool intervalid = false;
        valid = isValid(motion->state);
        if (!valid && add && addIntermediateState_ && motion->parent && motion->parent->stateValid == Valid)
        {
            if (start)
                intervalid = checkInterMotion2(motion->parent, motion, start);
            else 
                intervalid = checkInterMotion2(motion, motion->parent, start);
        }

        CellDiscretizationData &disc = start ? dStart_ : dGoal_;
        if (valid)
            motion->stateValid = Valid;
        else 
        {
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
            motion->stateValid = InValid;
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
            if (opt_->isFinite(motion->cost))
            {
                std::unordered_set<Cell *> cells;
                setMotionInfinityCostWithDisable(motion, cells);
                if (motion->cell->data.size() == 1)
                    cells.erase(motion->cell);
                else 
                    motion->cell->auxData->disabled--;
                disc.removeMotion(motion);
                for (auto & cell : cells)
                {
                    disc.updateCell(cell);
                    if (cell->data.size() == cell->auxData->disabled)
                    {
                        cell->auxData->cmotion = nullptr;
                        cell->auxData->mmotion = nullptr;
                    }
                    else
                    {
                        bool e1 = !cell->auxData->cmotion || !opt_->isFinite(cell->auxData->cmotion->cost);
                        bool e2 = !cell->auxData->mmotion || !opt_->isFinite(cell->auxData->mmotion->cost);
                        if (e1 || e2)
                        {
                            disc.enableSort(cell);
                            if (e1)
                                cell->auxData->cmotion = *std::min_element(cell->data.begin(), cell->data.end() - cell->auxData->disabled, mc_);
                            if (e2)
                                cell->auxData->mmotion = *std::max_element(cell->data.begin(), cell->data.end() - cell->auxData->disabled, mc_);
                        }
                    }
                }
            }
            else
            {
                Cell *cell = motion->cell;
                motion->cell->auxData->disabled--;
                disc.removeMotion(motion);
                if (!cell->removed)
                    disc.updateCell(cell);
            }
            removeInvalidCertificate(motion, start);
        }
    }
    return motion->stateValid == Valid;
}

void ompl::geometric::BiHSCCellstar::removeInvalidCertificate(Motion *motion, bool start)
{
    if (motion->sce != nullptr)
    {
        SafetyCertificateWithElems *msce = motion->sce;
        std::transform(motion->scd.begin(), motion->scd.end(), msce->sc->confidence_.begin(), [](double a){return 0.9*a;});
        std::vector<Motion *> valid, invalid;
        for (auto it = msce->objects.begin(); it != msce->objects.end();)
        {
            Motion *scm = *it;
            bool outside = false;
            for (std::size_t i = 0; i < scm->scd.size(); i++)
            {
                if (scm->scd[i] > msce->sc->confidence_[i])
                {
                    outside = true;
                    break;
                }
            }
            if (outside)
            {
                scm->sce = nullptr;
                std::fill(scm->scd.begin(), scm->scd.end(), 0.0);
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
        for (auto &vmotion : valid)
        {
            base::SafetyCertificate *sc = new base::SafetyCertificate(si_, certificateDim_);
            si_->copyState(sc->state, vmotion->state);
            delete sc->contact;
            sc->contact = nullptr;

            SafetyCertificateWithElems *sce = new SafetyCertificateWithElems();
            sce->sc = sc;
            sce->objects.push_back(vmotion);
            snne_.push_back(sce);

            vmotion->sce = sce;
            vmotion->sce->sc->confidence_ = distanceCertificate_(msce->sc->state, vmotion->state);
            if (strictCertificate_ && useCollisionCertificateChecker_)
            {
                base::SafetyCertificate *nsc = onn_->nearest(sc);
                std::vector<double> temp = distanceCertificate_(nsc->state, vmotion->state);
                std::transform(temp.begin(), temp.end(), vmotion->sce->sc->confidence_.begin(), vmotion->sce->sc->confidence_.begin(), [](double a, double b){return std::min(0.9*a, b);} );
            }
            else 
            {
                for (auto &invalidm : invalid)
                {
                    std::vector<double> temp = distanceCertificate_(invalidm->state, vmotion->state);
                    std::transform(temp.begin(), temp.end(), vmotion->sce->sc->confidence_.begin(), vmotion->sce->sc->confidence_.begin(), [](double a, double b){return std::min(0.9*a, b);} );
                }
            }
        }
    }
}

bool ompl::geometric::BiHSCCellstar::checkInterMotion1(Motion *smotion, Motion *gmotion, bool /*start*/)
{
    /*assume smotion, gmotion are valid*/
    bool result = true;
    const base::State *s1 = smotion->state, *s2 = gmotion->state;
    int nd = si_->getStateSpace()->validSegmentCount(s1, s2);
    if (nd >= 2)
    {
        std::queue<std::pair<int, int>> pos;
        pos.emplace(1, nd - 1);

        base::SafetyCertificate *sc = new base::SafetyCertificate(si_, certificateDim_);
        /* repeatedly subdivide the path segment in the middle (and check the middle) */
        while (!pos.empty())
        {
            std::pair<int, int> x = pos.front();

            int mid = (x.first + x.second) / 2;
            si_->getStateSpace()->interpolate(s1, s2, (double)mid / (double)nd, sc->state);

            if (useCollisionCertificateChecker_ && onn_->size() > oscNum_)
            {
                std::vector<base::SafetyCertificate *> nsc;
                onn_->nearestK(sc, 5, nsc);
                if (collisionCertificateChecker_(sc->state, nsc))
                {	
                    result = false;
                    freeCertificate(sc);
                    break;
                }
            }

            if (!isValid(sc))
            {
                result = false;
                break;
            }

            pos.pop();
            if (x.first < mid)
                pos.emplace(x.first, mid - 1);
            if (x.second > mid)
                pos.emplace(mid + 1, x.second);
        }

        if (result)
            freeCertificate(sc);
    }
    return result;
}

bool ompl::geometric::BiHSCCellstar::checkInterMotion2(Motion *smotion, Motion *gmotion, bool start)
{
    /*assume smotion, gmotion are valid*/
    bool valid = true;
    const base::State *s1 = smotion->state, *s2 = gmotion->state;
    int nd = si_->getStateSpace()->validSegmentCount(s1, s2);
    if (nd >= 2)
    {
        base::SafetyCertificate *sc = new base::SafetyCertificate(si_, certificateDim_);
        int i = start ? 1 : nd - 1;
        while (start ? i < nd : i > 0)
        {
            si_->getStateSpace()->interpolate(s1, s2, (double)i / (double)nd, sc->state);
            if (useCollisionCertificateChecker_ && onn_->size() > oscNum_)
            {
                std::vector<base::SafetyCertificate *> nsc;
                onn_->nearestK(sc, 5, nsc);
                if (collisionCertificateChecker_(sc->state, nsc))
                {	
                    valid = false;
                    freeCertificate(sc);
                    break;
                }
            }

            if (!isValid(sc))
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
        else
            freeCertificate(sc);
    }
    return valid;
}

void ompl::geometric::BiHSCCellstar::addIntermediateMotion(Motion *smotion, Motion *gmotion, bool start, Motion *last)
{
    CellDiscretizationData &disc = start ? dStart_ : dGoal_;
    Coord xcoord(projectionEvaluator_->getDimension());
    projectionEvaluator_->computeCoordinates(last->state, xcoord);
    disc.addMotion(last, xcoord);
    if (!opt_->isFinite(last->cost))
        last->cell->auxData->disabled++;
    if (last->cell->auxData->disabled)
        disc.updateCell(last->cell);
    last->valid = true;
    last->stateValid = Valid;
    Motion *v = start ? smotion : gmotion;
    Motion *inv = start ? gmotion : smotion;
    connectToPmotion(last, v, start);
    last->parent->children.push_back(last);
    insertNeighbor(last, v);
    disc.getGrid().updateNbh(last->cell, v->cell);
    if (inv->stateValid == Valid)
        insertInvalidNeighbor(last, inv);

    if (!last->cell->auxData->cmotion || opt_->isCostBetterThan(last->cost, last->cell->auxData->cmotion->cost))
        last->cell->auxData->cmotion = last;
    if (!last->cell->auxData->mmotion || !opt_->isCostBetterThan(last->cost, last->cell->auxData->mmotion->cost))
        last->cell->auxData->mmotion = last;
    if (lazyNode_)
    {
        last->sce = inv->sce;
        last->scd = distanceCertificate_(last->sce->sc->state, last->state);
        last->sce->objects.push_back(last);
    }
}

bool ompl::geometric::BiHSCCellstar::backPathRewireMotion(Motion *motion, bool start, Motion *&pmotion)
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
            std::size_t num = c->data.size() - c->auxData->disabled;
            if (!num)
                break;
            if (motion->cell != c && motion->cell->auxData->mmotion)
            {
                if (opt_->isCostBetterThan(motion->cell->auxData->mmotion->cost, c->auxData->cmotion->cost))
                    continue;
                base::Cost cost1 = opt_->combineCosts(c->auxData->cmotion->cost, opt_->motionCost(c->auxData->cmotion->state, motion->state));
                base::Cost cost2 = opt_->combineCosts(motion->cell->auxData->mmotion->cost, opt_->combineCosts(motion->cell->auxData->mmotion->cost, motion->cell->auxData->mmotion->cost));
                if (opt_->isCostBetterThan(cost2, cost1))
                    continue;
            }
            std::vector<Motion *> nbh;
            nbh.reserve(num);
            if (!c->auxData->disabled)
            {
                std::sort(c->data.begin(), c->data.end(), compareFn);
                nbh.insert(nbh.end(), c->data.begin(), c->data.end());
            }
            else 
            {
                for (std::size_t i = 0; i < c->data.size() && nbh.size() < num; i++)
                {
                    if (opt_->isFinite(c->data[i]->cost))
                        nbh.push_back(c->data[i]);
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
                    if (start ? checkInterMotion(nb, motion, start) : checkInterMotion(motion, nb, start))
                    {
                        insertNeighbor(nb, motion);
                        disc.getGrid().updateNbh(nb->cell, motion->cell);
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

void ompl::geometric::BiHSCCellstar::removeFromSafetyCerficate(Motion *motion)
{
    if (motion->sce)
        removeFromVector(motion->sce->objects, motion);
}
