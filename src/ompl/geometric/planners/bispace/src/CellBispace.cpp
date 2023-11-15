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

#include "ompl/geometric/planners/bispace/CellBispace.h"
#include "ompl/base/goals/GoalSampleableRegion.h"
#include "ompl/base/objectives/PathLengthOptimizationObjective.h"
#include "ompl/tools/config/SelfConfig.h"

ompl::geometric::CellBispace::CellBispace(const base::SpaceInformationPtr &si) : base::Planner(si, "CellBispace")
  ,mc_(opt_), bh_(mc_)
{
    specs_.recognizedGoal = base::GOAL_SAMPLEABLE_REGION;
    specs_.directed = true;
    specs_.approximateSolutions = false;

    Planner::declareParam<double>("range", this, &CellBispace::setRange, &CellBispace::getRange, "0.:0.1:10000.");
    Planner::declareParam<bool>("lazy_node", this, &CellBispace::setLazyNode, &CellBispace::getLazyNode, "0,1");
    Planner::declareParam<bool>("use_bispace", this, &CellBispace::setUseBispace, &CellBispace::getUseBispace, "0,1");
}

void ompl::geometric::CellBispace::setup()
{
    Planner::setup();
    tools::SelfConfig sc(si_, getName());
    sc.configurePlannerRange(maxDistance_);
    sc.configureProjectionEvaluator(projectionEvaluator_);

    projectionEvaluator_->mulCellSizes(2.0); // todo delete

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
    
    symmetric_ = si_->getStateSpace()->hasSymmetricDistance() && si_->getStateSpace()->hasSymmetricInterpolate();
}

void ompl::geometric::CellBispace::clear()
{
    setup_ = false;

    Planner::clear();
    sampler_.reset();

    dStart_.clear();
    dGoal_.clear();

    startMotions_.clear();
    goalMotions_.clear();

    pnullStartMotions_.clear();
    pnullGoalMotions_.clear();

    invalidStartMotions_.clear();
    invalidGoalMotions_.clear();

    connectionPoint_.clear();
}

ompl::base::PlannerStatus ompl::geometric::CellBispace::solve(const base::PlannerTerminationCondition &ptc)
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
    while (!ptc)
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
                if (rewireSort_)
                {
                    motion->cell->data->cmotion = motion;
                    motion->cell->data->mmotion = motion;
                }
            }
        }

        batchGrow(startTree);

        for (auto & pair : connectionPoint_)
        {
            if (opt_->isFinite(pair.first->cost) && opt_->isFinite(pair.second->cost))
            {
                if (isPathValid(pair.first, pair.second))
                {
                    solved = true;
                    bestStartMotion = pair.first;
                    bestGoalMotion = pair.second;
                    break;
                }
            }
        }
        //if (solved) // todo delete
        //    break;
        if (lazyNode_)
            removeInvalidMotions();
    }

    if (solved)
    {
        ptc.terminate();
        processSolution(bestStartMotion, bestGoalMotion);
    }

    OMPL_INFORM("%s: Created %u (%u start + %u goal) states in %u cells (%u start + %u goal)",
                getName().c_str(), dStart_.getMotionCount() + dGoal_.getMotionCount(), dStart_.getMotionCount(),
                dGoal_.getMotionCount(), dStart_.getCellCount() + dGoal_.getCellCount(), dStart_.getCellCount(),
                dGoal_.getCellCount());
    return solved ? base::PlannerStatus::EXACT_SOLUTION : base::PlannerStatus::TIMEOUT;
}

ompl::base::PlannerStatus ompl::geometric::CellBispace::prepareSolve(const base::PlannerTerminationCondition &ptc)
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
            if (rewireSort_)
            {
                motion->cell->data->cmotion = motion;
                motion->cell->data->mmotion = motion;
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
        if (rewireSort_)
        {
            motion->cell->data->cmotion = motion;
            motion->cell->data->mmotion = motion;
        }
    }

    if (dGoal_.getMotionCount() == 0)
    {
        OMPL_ERROR("%s: Unable to sample any valid states for goal tree", getName().c_str());
        return base::PlannerStatus::INVALID_GOAL;
    }

    if (useBispace_)
    {
        scoord_ = Coord(projectionEvaluator_->getDimension());
        projectionEvaluator_->computeCoordinates(startMotions_[0]->state, scoord_);
        Coord gcoord(projectionEvaluator_->getDimension());
        projectionEvaluator_->computeCoordinates(goalMotions_[0]->state, gcoord);
        dcoord_ = (gcoord - scoord_).cast<double>();
        dc_ = dcoord_.norm();
        dcoord_ /= dc_;
        dc_ *= 0.5;
    }

    return base::PlannerStatus::PREPARE_SUCCESS;
}

void ompl::geometric::CellBispace::processSolution(const Motion *bestStartMotion, const Motion *bestGoalMotion)
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

bool ompl::geometric::CellBispace::batchGrow(bool &startTree)
{
    bool nconnect = false;
    TreeGrowingInfo tgi;
    tgi.xstate = si_->allocState();
    for (unsigned int i = 0; i < 1; i++)
    {
        tgi.start = startTree;
        startTree = !startTree;
        nconnect = nconnect || growTree(tgi);
    }
    si_->freeState(tgi.xstate);
    return nconnect;
}

bool ompl::geometric::CellBispace::growTree(TreeGrowingInfo &tgi)
{
    bool nconnect = false;
    sampler_->sampleUniform(tgi.xstate);
    if (useBispace_ && !growCurrentTree(tgi.xstate, tgi.start))
        return nconnect;
    Motion *root = tgi.start ? startMotions_[0] : goalMotions_[0];
    Motion *motion = backGrowTree(root, tgi.xstate, tgi.start);
    if (!motion)
        return nconnect;
    CellDiscretizationData &disc = tgi.start ? dStart_ : dGoal_;
    CellDiscretizationData &otherDisc = !tgi.start ? dStart_ : dGoal_;
    std::vector<Motion *> &pnullMotions = tgi.start ? pnullStartMotions_ : pnullGoalMotions_;
    Cell *ocell = otherDisc.getGrid().getCell(motion->cell->coord);
    if (!ocell)
        return nconnect;
    Motion *connectOther = nullptr;
    std::size_t inConnection = 0, osize = ocell->data->motions.size();
    while (osize)
    {
        std::size_t index = rng_.uniformInt(0, osize - 1);
        connectOther = ocell->data->motions[index];
        if (!connectOther->inConnection)
        {
            nconnect = true;
            break;
        }
        else 
        {
            inConnection++;
            std::iter_swap(ocell->data->motions.begin() + index, ocell->data->motions.end() - inConnection);
        }
        osize--;
    }
    if (nconnect)
    {
        Motion *connect = new Motion(si_);
        si_->copyState(connect->state, connectOther->state);
        connectToPmotion(connect, motion, tgi.start);
        connect->parent->children.push_back(connect);
        disc.addMotion(connect, motion->cell->coord);
        connect->stateValid = connectOther->stateValid;
        if (rewireSort_)
        {
            if (!connect->cell->data->cmotion || opt_->isCostBetterThan(connect->cost, connect->cell->data->cmotion->cost))
                connect->cell->data->cmotion = connect;
            if (!connect->cell->data->mmotion || opt_->isCostBetterThan(connect->cell->data->mmotion->cost, connect->cost))
                connect->cell->data->mmotion = connect;
        }
        if (!pnullMotions.empty())
            rewireTree(bh_, connect, tgi.start);
        Motion *startMotion = connect, *goalMotion = connectOther;
        if (!tgi.start)
        {
            startMotion = connectOther;
            goalMotion = connect;
        }
        startMotion->inConnection = true;
        goalMotion->inConnection = true;
        connectionPoint_.emplace_back(startMotion, goalMotion);
    }
    return nconnect;
}

ompl::geometric::CellBispace::Motion* ompl::geometric::CellBispace::backGrowTree(Motion *root, base::State *dstate, bool start)
{
    Motion *child = nullptr, *leaf = nullptr;
    CellDiscretizationData &disc = start ? dStart_ : dGoal_;
    std::vector<Motion *> &pnullMotions = start ? pnullStartMotions_ : pnullGoalMotions_;
    while (true)
    {
        if (!lazyNode_ && !isValid(dstate))
        {
            if (child)
                pnullMotions.push_back(child);
            break;
        }
        Motion *motion = new Motion(si_);
        si_->copyState(motion->state, dstate);
        if (!lazyNode_)
            motion->stateValid = Valid;
        Coord xcoord(projectionEvaluator_->getDimension());
        projectionEvaluator_->computeCoordinates(motion->state, xcoord);
        disc.addMotion(motion, xcoord);
        motion->cost = opt_->infiniteCost();
        motion->cell->data->disabled++;
        if (child)
        {
            connectToPmotion(child, motion, start);
            child->parent->children.push_back(child);
        }
        if (!leaf)
            leaf = motion;
        Motion *pmotion = nullptr;
        if (!backPathRewireMotionSimple(motion, start, pmotion))
        {
            double d = start ? si_->distance(root->state, motion->state) : si_->distance(motion->state, root->state);
            if (d > maxDistance_)
            {
                if (start)
                    si_->getStateSpace()->interpolate(root->state, motion->state, std::max(1.0 - maxDistance_ / d, 0.5), dstate);
                else 
                    si_->getStateSpace()->interpolate(motion->state, root->state, std::min(maxDistance_ / d, 0.5), dstate);
                if (si_->equalStates(root->state, dstate) || si_->equalStates(motion->state, dstate))
                {
                    pnullMotions.push_back(motion);
                    break;
                }
            }
            else 
                pmotion = root;
        }
        if (pmotion)
        {
            connectToPmotion(motion, pmotion, start);
            motion->parent->children.push_back(motion);
            if (rewireSort_)
                enableMotionInDiscWithRewireSort(motion);
            else
                enableMotionInDisc(motion);
            if (!pnullMotions.empty())
                rewireTree(bh_, motion, start);
            if (!pnullMotions.empty() && leaf && leaf != motion)
                rewireTree(bh_, leaf, start);
            break;
        }
        child = motion;
    }
    return leaf;
}

void ompl::geometric::CellBispace::rewireTree(BinaryHeap<Motion *, MotionCompare> &bh, Motion *m, bool start)
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
                    if (rewireSort_)
                        enableMotionInDiscWithRewireSort(nb);
                    else
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

void ompl::geometric::CellBispace::updateQueue(BinaryHeap<Motion *, MotionCompare> &bh, Motion *motion)
{
    if (motion->handle != nullptr)
        bh.update(motion->handle);
    else
        motion->handle = bh.insert(motion);
}

bool ompl::geometric::CellBispace::growCurrentTree(const base::State *state, bool start) const
{
    double dist1 = si_->distance(startMotions_[0]->state, state);
    double dist2 = si_->distance(state, goalMotions_[0]->state);
    if (start == (dist1 <= dist2))
        return true;
    return false; // todo delete
    Coord coord(projectionEvaluator_->getDimension());
    projectionEvaluator_->computeCoordinates(state, coord);
    Coord dist = (coord - scoord_).array() + 1;;
    if (start)
        dist = dist.array() - 2;
    double d = dist.cast<double>().dot(dcoord_);
    return start == (d <= dc_);
}

bool ompl::geometric::CellBispace::backPathRewireMotionSimple(Motion *motion, bool start, Motion *&pmotion)
{
    bool valid = false;
    OrderCellsByCost ocbc(opt_);
    CostMotionCompare compareFn(motion, opt_, start);
    for (auto & cv : motion->cell->nbh)
    {
        std::size_t numc = cv.size(), ic = numc - 1;
        while (ic < numc)
        {
            if (rewireSort_ && ic)
                std::sort(cv.begin(), cv.begin() + ic + 1, ocbc);
            Cell *c = cv[ic];
            ic--;
            std::size_t num = c->data->motions.size() - c->data->disabled;
            if (!num)
            {
                if (rewireSort_)
                    break;
                else 
                    continue;
            }
            pmotion = *std::min_element(c->data->motions.begin(), c->data->motions.end(), compareFn);
            valid = true;
            break;
        }
        if (valid)
            break;
    }
    return valid;
}

bool ompl::geometric::CellBispace::isPathValid(Motion *motion, Motion *otherMotion)
{
    bool valid = true;
    if (!isPathValid(motion, true))
        valid = false;
    if (!isPathValid(otherMotion, false))
        valid = false;
    return valid;
}

bool ompl::geometric::CellBispace::isPathValid(Motion *motion, bool start)
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
        if (!checkMotion(motion->parent, motion, start))
        {
            removeFromParent(motion);
            motion->parent = nullptr;
            tvalid = false;
            if (backRewire_)
            {
                Motion *ppmotion = nullptr;
                if (backPathRewireMotion(motion, start, ppmotion))
                {
                    tvalid = isPathValidInter(ppmotion, start);
                    if (!motion->parent)
                    {
                        connectToPmotion(motion, ppmotion, start);
                        motion->parent->children.push_back(motion); 
                        if (opt_->isFinite(motion->cost))
                        {
                            if (rewireSort_)
                                enableMotionInDiscWithRewireSort(motion);
                            else
                                enableMotionInDisc(motion);
                        }
                    }
                }
                else 
                    pnullMotions.push_back(motion);
            }
            else 
                pnullMotions.push_back(motion);
            if (!tvalid)
            {
                rewireTree(bh_, pmotion, start);
                break;
            }
        }
    }
    return tvalid;
}

bool ompl::geometric::CellBispace::isPathValidInter(Motion *motion, bool start)
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
        if (!checkMotion(motion->parent, motion, start))
        {
            removeFromParent(motion);
            motion->parent = nullptr;
            tvalid = false;
            pnullMotions.push_back(motion);
            rewireTree(bh_, pmotion, start); // todo
            break;
        }
    }
    return tvalid;
}

bool ompl::geometric::CellBispace::isStateValid(Motion *motion, bool start)
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
        if (!isValid(motion, start))
        {
            tvalid = false;
            if (backRewire_)
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
                        if (!last->parent)
                        {
                            connectToPmotion(last, plast, start);
                            last->parent->children.push_back(last);
                            if (opt_->isFinite(last->cost))
                            {
                                if (rewireSort_)
                                    enableMotionInDiscWithRewireSort(last);
                                else
                                    enableMotionInDisc(last);
                            }
                        }
                    }
                    else 
                        pnullMotions.push_back(last);
                }
            }
            if (!tvalid)
            {
                rewireTree(bh_, pmotion, start); // todo
                break;
            }
        }
    }
    return tvalid;
}

bool ompl::geometric::CellBispace::backPathRewireMotion(Motion *motion, bool start, Motion *&pmotion)
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
        }
        if (valid)
            break;
    }
    return valid;
}

void ompl::geometric::CellBispace::removeInvalidMotions()
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

void ompl::geometric::CellBispace::removeInvalidMotionsDisc()
{
    for (auto & pnull : invalidStartMotions_)
    {
        for (auto & child : pnull->children)
        {
            child->parent = nullptr;
            pnullStartMotions_.push_back(child);
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
        }
        freeMotion(pnull);
    }
    invalidGoalMotions_.clear();
}

void ompl::geometric::CellBispace::enableMotionInDisc(Motion *motion)
{
    motion->cell->data->disabled--;
    for (auto & child : motion->children)
        enableMotionInDisc(child);
}

void ompl::geometric::CellBispace::enableMotionInDiscWithRewireSort(Motion *motion)
{
    motion->cell->data->disabled--;
    if (!motion->cell->data->cmotion || opt_->isCostBetterThan(motion->cost, motion->cell->data->cmotion->cost))
        motion->cell->data->cmotion = motion;
    if (!motion->cell->data->mmotion || opt_->isCostBetterThan(motion->cell->data->mmotion->cost, motion->cost))
        motion->cell->data->mmotion = motion;
    for (auto & child : motion->children)
        enableMotionInDiscWithRewireSort(child);
}

void ompl::geometric::CellBispace::disableCheck(Motion *motion, const std::string &context) const
{
    std::size_t pos = 0;
    for (auto & m : motion->cell->data->motions)
    {
        if (!opt_->isFinite(m->cost))
            pos++;
    }
    if (pos != motion->cell->data->disabled)
        OMPL_ERROR("%s: Error in %s!", getName().c_str(), context.c_str());
    for (auto & child : motion->children)
        disableCheck(child, context);
}

void ompl::geometric::CellBispace::connectToPmotion(Motion *motion, Motion *pmotion, bool start) const
{
    motion->parent = pmotion;
    motion->root = pmotion->root;
    motion->incCost = start ? opt_->motionCost(pmotion->state, motion->state) : opt_->motionCost(motion->state, pmotion->state);
    motion->cost = opt_->combineCosts(pmotion->cost, motion->incCost);
    if (opt_->isFinite(motion->cost))
        updateChildCosts(motion);
}

void ompl::geometric::CellBispace::updateChildCosts(Motion *motion) const
{
    for (auto & child : motion->children)
    {
        child->cost = opt_->combineCosts(motion->cost, child->incCost);
        child->root = motion->root;
        updateChildCosts(child);
    }
}

void ompl::geometric::CellBispace::insertNeighbor(Motion *pmotion, Motion *motion)
{
    pmotion->nbh.insert(motion);
    motion->nbh.insert(pmotion);
}

void ompl::geometric::CellBispace::insertInvalidNeighbor(Motion *pmotion, Motion *motion)
{
    pmotion->invalidnbh.insert(motion);
    motion->invalidnbh.insert(pmotion);
}

void ompl::geometric::CellBispace::removeFromNeighbor(Motion *motion)
{
    for (auto & pmotion : motion->nbh)
        pmotion->nbh.erase(motion);
}

void ompl::geometric::CellBispace::removeFromInvalidNeighbor(Motion *motion)
{
    for (auto & pmotion : motion->invalidnbh)
        pmotion->invalidnbh.erase(motion);
}

bool ompl::geometric::CellBispace::isValidNeighbor(Motion *motion, Motion *nb) const
{
    return motion->nbh.find(nb) != motion->nbh.end();
}

bool ompl::geometric::CellBispace::isInvalidNeighbor(Motion *motion, Motion *nb) const
{
    return motion->invalidnbh.find(nb) != motion->invalidnbh.end();
}

void ompl::geometric::CellBispace::removeFromParent(Motion *motion)
{
    if (!motion->parent)
        return;
    removeFromVector(motion->parent->children, motion);
}

void ompl::geometric::CellBispace::removeFromPnull(std::vector<Motion *> &pnullMotions, Motion *motion)
{
    removeFromVector(pnullMotions, motion);
}

bool ompl::geometric::CellBispace::removeFromVector(std::vector<Motion *> &motions, Motion *motion)
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

void ompl::geometric::CellBispace::setMotionInfinityCost(Motion *motion) const
{
    motion->cost = opt_->infiniteCost();
    motion->cell->data->disabled++;
    for (auto & child : motion->children)
        setMotionInfinityCost(child);
}

void ompl::geometric::CellBispace::setMotionInfinityCostWithRewireSort(Motion *motion, std::unordered_set<Cell *> &cells) const
{
    motion->cost = opt_->infiniteCost();
    motion->cell->data->disabled++;
    cells.insert(motion->cell);
    for (auto & child : motion->children)
        setMotionInfinityCostWithRewireSort(child, cells);
}

// check 
bool ompl::geometric::CellBispace::isValid(const base::State *state)
{
    return si_->isValid(state);
}

bool ompl::geometric::CellBispace::isValid(Motion *motion, bool start)
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
        removeFromParent(motion);
        motion->parent = nullptr;
        if (start)
            invalidStartMotions_.push_back(motion);
        else 
            invalidGoalMotions_.push_back(motion);
        std::unordered_set<Cell *> cells;
        if (opt_->isFinite(motion->cost))
        {
            if (rewireSort_)
                setMotionInfinityCostWithRewireSort(motion, cells);
            else
                setMotionInfinityCost(motion);
        }
        CellDiscretizationData &disc = start ? dStart_ : dGoal_;
        motion->cell->data->disabled--;
        disc.removeMotion(motion);
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
                    disc.enableSort(cell);
                    if (e1)
                        cell->data->cmotion = *std::min_element(cell->data->motions.begin(), cell->data->motions.end() - cell->data->disabled, mc_);
                    if (e2)
                        cell->data->mmotion = *std::max_element(cell->data->motions.begin(), cell->data->motions.end() - cell->data->disabled, mc_);
                }
            }
        }
    }
    return motion->stateValid == Valid;
}

bool ompl::geometric::CellBispace::checkMotion(Motion *pmotion, Motion *motion, bool start)
{
    if (!motion->valid)
    {
        CellDiscretizationData &disc = start ? dStart_ : dGoal_;
        if (!checkInterMotion(pmotion, motion, start))
        {
            if (opt_->isFinite(motion->cost))
            {
                if (rewireSort_)
                {
                    std::unordered_set<Cell *> cells;
                    setMotionInfinityCostWithRewireSort(motion, cells);
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
                                disc.enableSort(cell);
                                if (e1)
                                    cell->data->cmotion = *std::min_element(cell->data->motions.begin(), cell->data->motions.end() - cell->data->disabled, mc_);
                                if (e2)
                                    cell->data->mmotion = *std::max_element(cell->data->motions.begin(), cell->data->motions.end() - cell->data->disabled, mc_);
                            }
                        }
                    }
                }
                else
                    setMotionInfinityCost(motion);
            }
            insertInvalidNeighbor(pmotion, motion);
        }
        else 
        {
            motion->valid = true;
            insertNeighbor(pmotion, motion);
            disc.getGrid().updateNbh(pmotion->cell, motion->cell);
        }
    }
    return motion->valid;
}

bool ompl::geometric::CellBispace::checkInterMotion(Motion *pmotion, Motion *motion, bool start)
{
    bool valid = false;
    if (start || symmetric_)
        valid = checkInterMotion1(pmotion, motion, start);
    else 
        valid = checkInterMotion2(motion, pmotion, start);
    return valid;
}

bool ompl::geometric::CellBispace::checkInterMotion1(Motion *smotion, Motion *gmotion, bool start)
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

bool ompl::geometric::CellBispace::checkInterMotion2(Motion *smotion, Motion *gmotion, bool start)
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

void ompl::geometric::CellBispace::addIntermediateMotion(Motion *pmotion, Motion *motion, bool start, Motion *last)
{
    CellDiscretizationData &disc = start ? dStart_ : dGoal_;
    Coord xcoord(projectionEvaluator_->getDimension());
    projectionEvaluator_->computeCoordinates(last->state, xcoord);
    disc.addMotion(last, xcoord);
    if (!opt_->isFinite(last->cost))
        last->cell->data->disabled++;
    last->valid = true;
    last->stateValid = Valid;
    connectToPmotion(last, pmotion, start);
    last->parent->children.push_back(last);
    insertNeighbor(last, pmotion);
    disc.getGrid().updateNbh(last->cell, pmotion->cell);
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

void ompl::geometric::CellBispace::getPlannerData(base::PlannerData &data) const
{
    Planner::getPlannerData(data);
    dStart_.getPlannerData(data, 1, true);
    dGoal_.getPlannerData(data, 2, false);
}
