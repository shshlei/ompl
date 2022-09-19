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

#include "ompl/geometric/planners/hsc/BiHSCCell.h"
#include "ompl/base/goals/GoalSampleableRegion.h"
#include "ompl/base/objectives/PathLengthOptimizationObjective.h"
#include "ompl/tools/config/SelfConfig.h"

ompl::geometric::BiHSCCell::BiHSCCell(const base::SpaceInformationPtr &si) : base::Planner(si, "BiHSCCell")
  ,mc_(opt_), bh_(mc_)
{
    specs_.recognizedGoal = base::GOAL_SAMPLEABLE_REGION;
    specs_.directed = true;
    specs_.approximateSolutions = false;

    Planner::declareParam<double>("range", this, &BiHSCCell::setRange, &BiHSCCell::getRange, "0.:0.1:10000.");
    Planner::declareParam<double>("pen_distance", this, &BiHSCCell::setPenDistance, &BiHSCCell::getPenDistance, "0.:0.1:10000.");
    Planner::declareParam<double>("min_valid_path_fraction", this, &BiHSCCell::setMinValidPathFraction, &BiHSCCell::getMinValidPathFraction, "0.:0.05:1.");
    Planner::declareParam<bool>("lazy_path", this, &BiHSCCell::setLazyPath, &BiHSCCell::getLazyPath, "0,1");
    Planner::declareParam<bool>("lazy_node", this, &BiHSCCell::setLazyNode, &BiHSCCell::getLazyNode, "0,1");
    Planner::declareParam<bool>("add_intermediate_state", this, &BiHSCCell::setAddIntermediateState, &BiHSCCell::getAddIntermediateState, "0,1");
    Planner::declareParam<bool>("use_bispace", this, &BiHSCCell::setUseBispace, &BiHSCCell::getUseBispace, "0,1");
//    Planner::declareParam<bool>("treated_as_multi_subapce", this, &BiHSCCell::setTreatedAsMultiSubapce, &BiHSCCell::getTreatedAsMultiSubapce, "0,1");
    Planner::declareParam<bool>("strice_certificate", this, &BiHSCCell::setStrictCertificate, &BiHSCCell::getStrictCertificate, "0,1");
}

ompl::geometric::BiHSCCell::~BiHSCCell()
{
    freeMemory();
}

void ompl::geometric::BiHSCCell::setup()
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
    if (rewire_)
    {
        dStart_.setNeighborCell(2);
        dGoal_.setNeighborCell(2);
    }
    else 
        batch_ = 1;

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
    {
        lazyNode_ = false;
        rewire_ = false;
    }
    if (!rewire_)
        rewireSort_ = false;
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

void ompl::geometric::BiHSCCell::freeMemory()
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

void ompl::geometric::BiHSCCell::clear()
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

    invalidStartMotions_.clear();
    invalidGoalMotions_.clear();

    connectionPoint_.clear();

    if (onn_)
    {
        onn_->clear();
    }

    snne_.clear();
}

ompl::base::PlannerStatus ompl::geometric::BiHSCCell::solve(const base::PlannerTerminationCondition &ptc)
{
    base::PlannerStatus ps= prepareSolve(ptc);
    if (ps != base::PlannerStatus::PREPARE_SUCCESS)
        return ps;

    sampler_ = si_->allocStateSampler();

    OMPL_INFORM("%s: Starting planning with %u states already in datastructure.", getName().c_str(), (dStart_.getMotionCount() + dGoal_.getMotionCount()));

    bool startTree = true;
    bool solved = false;
    Motion *bestStartMotion = nullptr;
    Motion *bestGoalMotion = nullptr;
    while (!ptc && !solved)
    {
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
                if (rewire_ && rewireSort_)
                {
                    motion->cell->auxData->cmotion = motion;
                    motion->cell->auxData->mmotion = motion;
                }
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

    OMPL_INFORM("%s: Created %u (%u start + %u goal) states in %u cells (%u start (%u on boundary) + %u goal (%u on "
                "boundary))",
                getName().c_str(), dStart_.getMotionCount() + dGoal_.getMotionCount(), dStart_.getMotionCount(),
                dGoal_.getMotionCount(), dStart_.getCellCount() + dGoal_.getCellCount(), dStart_.getCellCount(),
                dStart_.getGrid().countExternal(), dGoal_.getCellCount(), dGoal_.getGrid().countExternal());
    return solved ? base::PlannerStatus::EXACT_SOLUTION : base::PlannerStatus::TIMEOUT;
}

ompl::base::PlannerStatus ompl::geometric::BiHSCCell::prepareSolve(const base::PlannerTerminationCondition &ptc)
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
            if (rewire_ && rewireSort_)
            {
                motion->cell->auxData->cmotion = motion;
                motion->cell->auxData->mmotion = motion;
            }
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
        if (rewire_ && rewireSort_)
        {
            motion->cell->auxData->cmotion = motion;
            motion->cell->auxData->mmotion = motion;
        }
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

void ompl::geometric::BiHSCCell::processSolution(const Motion *bestStartMotion, const Motion *bestGoalMotion)
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

bool ompl::geometric::BiHSCCell::batchGrow(bool &startTree)
{
    bool nconnect = false;
    TreeGrowingInfo tgi;
    tgi.xstate = si_->allocState();
    for (unsigned int i = 0; i < batch_; i++)
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

void ompl::geometric::BiHSCCell::growTree(TreeGrowingInfo &tgi, Motion *&startMotion, Motion *&goalMotion)
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
    connectToPmotion(motion, nb, tgi.start);
    motion->parent->children.push_back(motion);

    if (!lazyPath_)
        motion->valid = true;

    if (rewire_)
    {
        if (rewireSort_)
        {
            if (!motion->cell->auxData->cmotion || opt_->isCostBetterThan(motion->cost, motion->cell->auxData->cmotion->cost))
                motion->cell->auxData->cmotion = motion;
            if (!motion->cell->auxData->mmotion || !opt_->isCostBetterThan(motion->cost, motion->cell->auxData->mmotion->cost))
                motion->cell->auxData->mmotion = motion;
        }
        if (tgi.start ? !pnullStartMotions_.empty() : !pnullGoalMotions_.empty())
            rewireTree(bh_, motion, tgi.start);
    }

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

        if (rewire_)
        {
            if (rewireSort_)
            {
                if (!connect->cell->auxData->cmotion || opt_->isCostBetterThan(connect->cost, connect->cell->auxData->cmotion->cost))
                    connect->cell->auxData->cmotion = connect;
                if (!connect->cell->auxData->mmotion || !opt_->isCostBetterThan(connect->cost, connect->cell->auxData->mmotion->cost))
                    connect->cell->auxData->mmotion = connect;
            }
            if (tgi.start ? !pnullStartMotions_.empty() : !pnullGoalMotions_.empty())
                rewireTree(bh_, connect, tgi.start);
        }

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

void ompl::geometric::BiHSCCell::rewireTree(BinaryHeap<Motion *, MotionCompare> &bh, Motion *m, bool start)
{
    if (!rewire_)
        return;
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
                if (rewireSort_ && ic > 0)
                    std::sort(cv.begin(), cv.begin() + ic + 1, ocbd);
                Cell *c = cv[ic];
                ic--;
                std::size_t pos = c->auxData->disabled;
                if (!pos)
                {
                    if (rewireSort_)
                        break;
                    else 
                        continue;
                }
                if (rewireSort_ && motion->cell != c && c->auxData->mmotion)
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

void ompl::geometric::BiHSCCell::updateQueue(BinaryHeap<Motion *, MotionCompare> &bh, Motion *motion)
{
    if (motion->handle != nullptr)
        bh.update(motion->handle);
    else
        motion->handle = bh.insert(motion);
}

bool ompl::geometric::BiHSCCell::growCurrentTree(const base::State *state) const
{
    return growStartTree(state);
}

bool ompl::geometric::BiHSCCell::growStartTree(const base::State *state) const
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

double ompl::geometric::BiHSCCell::penetrationDistance(const base::State *nstate, const base::State *state, bool start) const
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

bool ompl::geometric::BiHSCCell::growCurrentTree(const base::State *state, std::size_t sub) const
{
    return growStartTree(state, sub);
}

bool ompl::geometric::BiHSCCell::growStartTree(const base::State *state, std::size_t sub) const
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

double ompl::geometric::BiHSCCell::penetrationDistance(const base::State *nstate, const base::State *state, bool start, std::size_t sub) const
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

bool ompl::geometric::BiHSCCell::isPathValid(Motion *motion, Motion *otherMotion, bool start)
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

bool ompl::geometric::BiHSCCell::isPathValid(Motion *motion, Motion *otherMotion)
{
    bool valid = true;
    if (!isPathValid(motion, true))
        valid = false;
    if (!isPathValid(otherMotion, false))
        valid = false;
    return valid;
}

bool ompl::geometric::BiHSCCell::isPathValid(Motion *motion, bool start)
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
            if (rewire_)
            {
                Motion *pmotion = nullptr;
                if (backPathRewireMotion(motion, start, pmotion))
                {
                    tvalid = isPathValidInter(pmotion, start);
                    connectToPmotion(motion, pmotion, start);
                    motion->parent->children.push_back(motion); 
                }
                else 
                {
                    pnullMotions.push_back(motion);
                    motion->valid = true;
                    nullMotions.push_back(motion);
                    motion->cell->auxData->root++;
                }
            }
            else 
                pnullMotions.push_back(motion);
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

bool ompl::geometric::BiHSCCell::isStateValid(Motion *motion, bool start)
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
            if (rewire_)
            {
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
                        if (tvalid)
                            enableMotionInDisc(disc, last);
                    }
                    else 
                    {
                        last->valid = true;
                        pnullMotions.push_back(last);
                        nullMotions.push_back(last);
                        last->cell->auxData->root++;
                    }
                }
            }
            if (!tvalid)
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

bool ompl::geometric::BiHSCCell::isPathValidInter(Motion *motion, bool start)
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
            if (rewire_)
            {
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
            }
            else 
                pnullMotions.push_back(motion);
            break;
        }
        motion = motion->parent;
    }
    return tvalid;
}

void ompl::geometric::BiHSCCell::removeInvalidMotionsDirectly()
{
    for (auto & pair : connectionPoint_)
    {
        pair.first->inConnection = false;
        pair.second->inConnection = false;
    }
    connectionPoint_.clear();
    removeInvalidMotionsDirectlyDisc();
}

void ompl::geometric::BiHSCCell::removeInvalidMotionsDirectlyDisc()
{
    if (!invalidStartMotions_.empty() || !pnullStartMotions_.empty())
    {
        for (auto & pnull : invalidStartMotions_)
            removeFromDisc(dStart_, pnull);
        for (auto & pnull : pnullStartMotions_)
            removeFromDisc(dStart_, pnull);
        invalidStartMotions_.clear();
        pnullStartMotions_.clear();
    }

    if (!invalidGoalMotions_.empty() || !pnullGoalMotions_.empty())
    {
        for (auto & pnull : invalidGoalMotions_)
            removeFromDisc(dGoal_, pnull);
        for (auto & pnull : pnullGoalMotions_)
            removeFromDisc(dGoal_, pnull);
        invalidGoalMotions_.clear();
        pnullGoalMotions_.clear();
    }
}

void ompl::geometric::BiHSCCell::removeFromDisc(CellDiscretizationData &disc, Motion *motion)
{
    disc.removeMotion(motion);
    for (auto & child : motion->children)
    {
        child->parent = nullptr;
        removeFromDisc(disc, child);
    }
    freeMotion(motion);
}

void ompl::geometric::BiHSCCell::removeInvalidMotions()
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

void ompl::geometric::BiHSCCell::removeInvalidMotionsDisc()
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

void ompl::geometric::BiHSCCell::enableMotionInDisc(CellDiscretizationData &disc, Motion *motion)
{
    std::unordered_set<Cell *> cells;
    enableMotionInDisc2(motion, cells);
    for (auto & cell : cells)
        disc.updateCell(cell);
}

void ompl::geometric::BiHSCCell::enableMotionInDisc2(Motion *motion, std::unordered_set<Cell *> &cells)
{
    motion->cell->auxData->disabled--;
    cells.insert(motion->cell);
    if (rewireSort_)
    {
        if (!motion->cell->auxData->cmotion || opt_->isCostBetterThan(motion->cost, motion->cell->auxData->cmotion->cost))
            motion->cell->auxData->cmotion = motion;
        if (!motion->cell->auxData->mmotion || !opt_->isCostBetterThan(motion->cost, motion->cell->auxData->mmotion->cost))
            motion->cell->auxData->mmotion = motion;
    }
    for (auto & child : motion->children)
        enableMotionInDisc2(child, cells);
}

void ompl::geometric::BiHSCCell::disableCheck(Motion *motion, const std::string &context) const
{
    std::size_t pos = 0;
    for (auto & m : motion->cell->data)
    {
        if (!opt_->isFinite(m->cost))
            pos++;
    }
    if (pos != motion->cell->auxData->disabled)
        OMPL_ERROR("%s: Error in %s!", getName().c_str(), context.c_str());
    for (auto & child : motion->children)
        disableCheck(child, context);
}

void ompl::geometric::BiHSCCell::connectToPmotion(Motion *motion, Motion *pmotion, bool start) const
{
    motion->parent = pmotion;
    motion->root = pmotion->root;
    motion->incCost = start ? opt_->motionCost(pmotion->state, motion->state) : opt_->motionCost(motion->state, pmotion->state);
    motion->cost = opt_->combineCosts(pmotion->cost, motion->incCost);
    if (opt_->isFinite(motion->cost))
        updateChildCosts(motion);
}

void ompl::geometric::BiHSCCell::updateChildCosts(Motion *motion) const
{
    for (auto & child : motion->children)
    {
        child->cost = opt_->combineCosts(motion->cost, child->incCost);
        child->root = motion->root;
        updateChildCosts(child);
    }
}

void ompl::geometric::BiHSCCell::insertNeighbor(Motion *pmotion, Motion *motion)
{
    pmotion->nbh.insert(motion);
    motion->nbh.insert(pmotion);
}

void ompl::geometric::BiHSCCell::insertInvalidNeighbor(Motion *pmotion, Motion *motion)
{
    pmotion->invalidnbh.insert(motion);
    motion->invalidnbh.insert(pmotion);
}

void ompl::geometric::BiHSCCell::removeFromNeighbor(Motion *motion)
{
    for (auto & pmotion : motion->nbh)
        pmotion->nbh.erase(motion);
}

void ompl::geometric::BiHSCCell::removeFromInvalidNeighbor(Motion *motion)
{
    for (auto & pmotion : motion->invalidnbh)
        pmotion->invalidnbh.erase(motion);
}

bool ompl::geometric::BiHSCCell::isValidNeighbor(Motion *motion, Motion *nb) const
{
    return motion->nbh.find(nb) != motion->nbh.end();
}

bool ompl::geometric::BiHSCCell::isInvalidNeighbor(Motion *motion, Motion *nb) const
{
    return motion->invalidnbh.find(nb) != motion->invalidnbh.end();
}

void ompl::geometric::BiHSCCell::removeFromParent(Motion *motion)
{
    if (!motion->parent)
        return;
    removeFromVector(motion->parent->children, motion);
}

void ompl::geometric::BiHSCCell::removeFromPnull(std::vector<Motion *> &pnullMotions, Motion *motion)
{
    if (removeFromVector(pnullMotions, motion))
        motion->cell->auxData->root--;
}

bool ompl::geometric::BiHSCCell::removeFromVector(std::vector<Motion *> &motions, Motion *motion)
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

void ompl::geometric::BiHSCCell::setMotionInfinityCost(Motion *motion) const
{
    motion->cost = opt_->infiniteCost();
    for (auto & child : motion->children)
        setMotionInfinityCost(child);
}

void ompl::geometric::BiHSCCell::setMotionInfinityCostWithDisable(Motion *motion, std::unordered_set<Cell *> &cells) const
{
    motion->cost = opt_->infiniteCost();
    motion->cell->auxData->disabled++;
    cells.insert(motion->cell);
    for (auto & child : motion->children)
        setMotionInfinityCostWithDisable(child, cells);
}

// check 
bool ompl::geometric::BiHSCCell::checkStartMotion(Motion *smotion, Motion *gmotion)
{
    if (!gmotion->valid)
    {
        if (!checkInterMotion(smotion, gmotion, true))
        {
            if (opt_->isFinite(gmotion->cost))
            {
                if (rewire_)
                {
                    std::unordered_set<Cell *> cells;
                    setMotionInfinityCostWithDisable(gmotion, cells);
                    for (auto & cell : cells)
                    {
                        dStart_.updateCell(cell);
                        if (rewireSort_)
                        {
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
                }
                else 
                    setMotionInfinityCost(gmotion);
            }
            if (rewire_)
                insertInvalidNeighbor(smotion, gmotion);
        }
        else 
        {
            gmotion->valid = true;
            if (rewire_)
            {
                insertNeighbor(gmotion, smotion);
                dStart_.getGrid().updateNbh(gmotion->cell, smotion->cell);
            }
        }
    }
    return gmotion->valid;
}

bool ompl::geometric::BiHSCCell::checkGoalMotion(Motion *smotion, Motion *gmotion)
{
    if (!smotion->valid)
    {
        if (!checkInterMotion(smotion, gmotion, false))
        {
            if (opt_->isFinite(smotion->cost))
            {
                if (rewire_)
                {
                    std::unordered_set<Cell *> cells;
                    setMotionInfinityCostWithDisable(smotion, cells);
                    for (auto & cell : cells)
                    {
                        dGoal_.updateCell(cell);
                        if (rewireSort_)
                        {
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
                }
                else 
                    setMotionInfinityCost(smotion);
            }
            if (rewire_)
                insertInvalidNeighbor(smotion, gmotion);
        }
        else 
        {
            smotion->valid = true;
            if (rewire_)
            {
                insertNeighbor(gmotion, smotion);
                dGoal_.getGrid().updateNbh(gmotion->cell, smotion->cell);
            }
        }
    }
    return smotion->valid;
}

bool ompl::geometric::BiHSCCell::checkInterMotion(Motion *smotion, Motion *gmotion, bool start)
{
    bool valid = false;
    if (addIntermediateState_)
        valid = checkInterMotion2(smotion, gmotion, start);
    else
        valid = checkInterMotion1(smotion, gmotion, start);
    return valid;
}

void ompl::geometric::BiHSCCell::getPlannerData(base::PlannerData &data) const
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

void ompl::geometric::BiHSCCell::getPlannerData(base::PlannerData &data, Motion *motion, int tag) const
{
    data.addEdge(base::PlannerDataVertex(motion->parent->state, tag), base::PlannerDataVertex(motion->state, tag));
    for (auto & child : motion->children)
    {
        getPlannerData(data, child, tag);
    }
}

// safety certificate 
bool ompl::geometric::BiHSCCell::isValid(base::SafetyCertificate *sc)
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

bool ompl::geometric::BiHSCCell::isValid(const base::State *state)
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

bool ompl::geometric::BiHSCCell::isValid(Motion *motion, bool start, bool add)
{
    if (!lazyNode_)
        return true;
    if (motion->stateValid == UnCkecked)
    {
        bool valid = true;
        bool intervalid = false;
        if (rewire_)
        {
            valid = isValid(motion->state);
            if (!valid && add && addIntermediateState_ && motion->parent && motion->parent->stateValid == Valid)
            {
                if (start)
                    intervalid = checkInterMotion2(motion->parent, motion, start);
                else 
                    intervalid = checkInterMotion2(motion, motion->parent, start);
            }
        }
        else
        {
            if (add && addIntermediateState_ && motion->parent && motion->parent->stateValid == Valid)
            {
                if (start)
                    intervalid = checkInterMotion2(motion->parent, motion, start);
                else 
                    intervalid = checkInterMotion2(motion, motion->parent, start);
                valid = intervalid;
            }
            if (valid)
                valid = isValid(motion->state);
        }

        CellDiscretizationData &disc = start ? dStart_ : dGoal_;
        if (valid)
        {
            motion->stateValid = Valid;
            if (intervalid)
            {
                motion->valid = true;
                if (rewire_)
                {
                    insertNeighbor(motion, motion->parent);
                    disc.getGrid().updateNbh(motion->cell, motion->parent->cell);
                }
            }
        }
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
                if (rewire_)
                {
                    setMotionInfinityCostWithDisable(motion, cells);
                    if (motion->cell->data.size() == 1)
                        cells.erase(motion->cell);
                    else 
                        motion->cell->auxData->disabled--;
                }
                else 
                    setMotionInfinityCost(motion);
                disc.removeMotion(motion);
                for (auto & cell : cells)
                {
                    disc.updateCell(cell);
                    if (rewireSort_)
                    {
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
            }
            else if (rewire_)
            {
                Cell *cell = motion->cell;
                cell->auxData->disabled--;
                disc.removeMotion(motion);
                if (!cell->removed)
                    disc.updateCell(cell);
            }
            removeInvalidCertificate(motion, start);
        }
    }
    return motion->stateValid == Valid;
}

void ompl::geometric::BiHSCCell::removeInvalidCertificate(Motion *motion, bool start)
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

bool ompl::geometric::BiHSCCell::checkInterMotion1(Motion *smotion, Motion *gmotion, bool /*start*/)
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

bool ompl::geometric::BiHSCCell::checkInterMotion2(Motion *smotion, Motion *gmotion, bool start)
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

void ompl::geometric::BiHSCCell::addIntermediateMotion(Motion *smotion, Motion *gmotion, bool start, Motion *last)
{
    CellDiscretizationData &disc = start ? dStart_ : dGoal_;
    Coord xcoord(projectionEvaluator_->getDimension());
    projectionEvaluator_->computeCoordinates(last->state, xcoord);
    disc.addMotion(last, xcoord);
    if (rewire_ && !opt_->isFinite(last->cost))
        last->cell->auxData->disabled++;
    if (last->cell->auxData->disabled) 
        disc.updateCell(last->cell);

    last->valid = true;
    last->stateValid = Valid;
    Motion *v = start ? smotion : gmotion;
    Motion *inv = start ? gmotion : smotion;
    connectToPmotion(last, v, start);
    last->parent->children.push_back(last);
    if (rewire_)
    {
        insertNeighbor(last, v);
        disc.getGrid().updateNbh(last->cell, v->cell);
        if (inv->stateValid == Valid)
            insertInvalidNeighbor(last, inv);
        if (rewireSort_)
        {
            if (!last->cell->auxData->cmotion || opt_->isCostBetterThan(last->cost, last->cell->auxData->cmotion->cost))
                last->cell->auxData->cmotion = last;
            if (!last->cell->auxData->mmotion || !opt_->isCostBetterThan(last->cost, last->cell->auxData->mmotion->cost))
                last->cell->auxData->mmotion = last;
        }
    }
    if (lazyNode_)
    {
        last->sce = inv->sce;
        last->scd = distanceCertificate_(last->sce->sc->state, last->state);
        last->sce->objects.push_back(last);
    }
}

bool ompl::geometric::BiHSCCell::backPathRewireMotion(Motion *motion, bool start, Motion *&pmotion)
{
    bool valid = false;
    if (!rewire_)
        return valid;
    OrderCellsByCost ocbc(opt_);
    CostMotionCompare compareFn(motion, opt_, start);
    CellDiscretizationData &disc = start ? dStart_ : dGoal_;
    for (auto & cv : motion->cell->nbh)
    {
        std::size_t numc = cv.size(), ic = numc - 1;
        while (ic < numc)
        {
            if (rewireSort_ && ic)
                std::sort(cv.begin(), cv.begin() + ic + 1, ocbc);
            Cell *c = cv[ic];
            ic--;
            std::size_t num = c->data.size() - c->auxData->disabled;
            if (!num)
            {
                if (rewireSort_)
                    break;
                else 
                    continue;
            }
            if (rewireSort_ && motion->cell != c && motion->cell->auxData->mmotion)
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
