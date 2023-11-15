/*********************************************************************
* Software License Agreement (BSD License)
*
*  Copyright (c) 2014, University of Toronto
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
*   * Neither the name of the University of Toronto nor the names of its
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

#include "ompl/base/samplers/adinformed/PathLengthDirectAdInfSampler.h"
#include "ompl/util/Exception.h"
#include "ompl/base/goals/GoalSampleableRegion.h"
#include "ompl/base/StateSpace.h"
#include "ompl/base/spaces/RealVectorStateSpace.h"

#include <memory>
#include <vector>

ompl::base::PathLengthDirectAdInfSampler::PathLengthDirectAdInfSampler(const ProblemDefinitionPtr &probDefn,
                                                                       const State *s1, const State *s2,
                                                                       const Cost &minCost, const Cost &maxCost,
                                                                       unsigned int maxNumberCalls)
  : AdInformedSampler(probDefn, s1, s2, minCost, maxCost, maxNumberCalls),
    informedIdx_(0u), uninformedIdx_(0u)
{
    baseSampler_ = AdInformedSampler::space_->allocStateSampler();

    // Check that the provided statespace is compatible and extract the necessary indices.
    // The statespace must either be R^n or SE(2) or SE(3).
    // If it is UNKNOWN, warn and treat it as R^n
    if (!AdInformedSampler::space_->isCompound())
    {
        if (AdInformedSampler::space_->getType() == STATE_SPACE_REAL_VECTOR)
        {
            informedIdx_ = 0u;
            uninformedIdx_ = 0u;
        }
        else if (AdInformedSampler::space_->getType() == STATE_SPACE_UNKNOWN)
        {
            OMPL_WARN("PathLengthDirectAdInfSampler: Treating the StateSpace of type \"STATE_SPACE_UNKNOWN\" as type \"STATE_SPACE_REAL_VECTOR\".");
            informedIdx_ = 0u;
            uninformedIdx_ = 0u;
        }
        else
        {
            throw Exception("PathLengthDirectAdInfSampler only supports Unknown, RealVector, SE2, and SE3 StateSpaces.");
        }
        informedSubSpace_ = AdInformedSampler::space_;
        uninformedSubSpace_ = StateSpacePtr();
        uninformedSubSampler_ = StateSamplerPtr();
    }
    else if (AdInformedSampler::space_->getType() == STATE_SPACE_SE2 || AdInformedSampler::space_->getType() == STATE_SPACE_SE3)
    {
        const CompoundStateSpace *compoundSpace = AdInformedSampler::space_->as<CompoundStateSpace>();
        if (compoundSpace->getSubspaceCount() != 2u)
        {
            throw Exception("The provided compound StateSpace is SE(2) or SE(3) but does not have exactly "
                            "2 subspaces.");
        }
        for (unsigned int idx = 0u; idx < compoundSpace->getSubspaceCount(); ++idx)
        {
            if (compoundSpace->getSubspace(idx)->getType() == STATE_SPACE_REAL_VECTOR)
            {
                informedIdx_ = idx;
            }
            else if (compoundSpace->getSubspace(idx)->getType() == STATE_SPACE_SO2)
            {
                uninformedIdx_ = idx;
            }
            else if (compoundSpace->getSubspace(idx)->getType() == STATE_SPACE_SO3)
            {
                uninformedIdx_ = idx;
            }
            else
            {
                throw Exception("The provided compound StateSpace is SE(2) or SE(3) but contains a "
                                "subspace that is not R^2, R^3, SO(2), or SO(3).");
            }
        }
        informedSubSpace_ = compoundSpace->getSubspace(informedIdx_);
        uninformedSubSpace_ = compoundSpace->getSubspace(uninformedIdx_);
        uninformedSubSampler_ = uninformedSubSpace_->allocStateSampler();
    }
    else
    {
        const CompoundStateSpace *compoundSpace = AdInformedSampler::space_->as<CompoundStateSpace>();
        if (compoundSpace->getSubspaceCount() != 1u)
            throw Exception("PathLengthDirectAdInfSampler only supports RealVector, SE2 and SE3 statespaces.");
        compound_ = true;
        if (compoundSpace->getSubspace(0)->getType() == STATE_SPACE_REAL_VECTOR)       
        {
            informedIdx_ = 0u;
            uninformedIdx_ = 0u;
            informedSubSpace_ = compoundSpace->getSubspace(0);
            uninformedSubSpace_ = StateSpacePtr();
            uninformedSubSampler_ = StateSamplerPtr();
        }
        else if (compoundSpace->getSubspace(0)->getType() == STATE_SPACE_SE2 || compoundSpace->getSubspace(0)->getType() == STATE_SPACE_SE3)
        {
            compoundCompound_ = true;
            const CompoundStateSpace *compoundCompoundSpace = compoundSpace->getSubspace(0)->as<CompoundStateSpace>();
            if (compoundCompoundSpace->getSubspaceCount() != 2u)
                throw Exception("The provided compound compound StateSpace is SE(2) or SE(3) but does not have exactly 2 subspaces.");
            for (unsigned int idx = 0u; idx < compoundCompoundSpace->getSubspaceCount(); ++idx)
            {
                if (compoundCompoundSpace->getSubspace(idx)->getType() == STATE_SPACE_REAL_VECTOR)
                {
                    informedIdx_ = idx;
                }
                else if (compoundCompoundSpace->getSubspace(idx)->getType() == STATE_SPACE_SO2)
                {
                    uninformedIdx_ = idx;
                }
                else if (compoundCompoundSpace->getSubspace(idx)->getType() == STATE_SPACE_SO3)
                {
                    uninformedIdx_ = idx;
                }
                else
                    throw Exception("The provided compound compound StateSpace is SE(2) or SE(3) but contains a subspace that is not R^2, R^3, SO(2), or SO(3).");
            }
            informedSubSpace_ = compoundCompoundSpace->getSubspace(informedIdx_);
            uninformedSubSpace_ = compoundCompoundSpace->getSubspace(uninformedIdx_);
            uninformedSubSampler_ = uninformedSubSpace_->allocStateSampler();           
        }
        else 
            throw Exception("PathLengthDirectAdInfSampler only supports RealVector, SE2 and SE3 statespaces.");
    }

    std::vector<double> startFocusVector = getInformedSubstate(this->s1_);
    std::vector<double> goalFocusVector  = getInformedSubstate(this->s2_);
    phsPtr_ = std::make_shared<ProlateHyperspheroidRing>(informedSubSpace_->getDimension(), 
                                                         &startFocusVector[0], &goalFocusVector[0]);

    if (minCost.value() > 0)
    {
        phsPtr_->setTransverseDiameter(minCost.value());
        phsPtr_->setExpansionFactor(maxCost.value()/minCost.value());
    }
    else 
        phsPtr_->setTransverseDiameter(maxCost.value());

    resetDirectCounter();
}

void ompl::base::PathLengthDirectAdInfSampler::update(const Cost &minCost, const Cost &maxCost)
{
    AdInformedSampler::update(minCost, maxCost);

    if (minCost.value() > 0)
    {
        phsPtr_->setTransverseDiameter(minCost.value());
        phsPtr_->setExpansionFactor(maxCost.value()/minCost.value());
    }
    else 
        phsPtr_->setTransverseDiameter(maxCost.value());

    resetDirectCounter();
}

void ompl::base::PathLengthDirectAdInfSampler::update(const Cost &maxCost)
{
    AdInformedSampler::update(this->minCost_, maxCost);

    if (phsPtr_->isProlateHyperspheroid())
        phsPtr_->setTransverseDiameter(maxCost.value());
    else 
        update(this->minCost_, maxCost);

    resetDirectCounter();
}

void ompl::base::PathLengthDirectAdInfSampler::update(const double factor)
{
    AdInformedSampler::update(this->maxCost_, Cost(factor * this->maxCost_.value()));
    phsPtr_->setExpansionFactor(factor * phsPtr_->getExternalExpansionFactor());

    resetDirectCounter();
}

bool ompl::base::PathLengthDirectAdInfSampler::sampleUniform(State *state)
{
    unsigned int iter = 0u;

    return sampleUniform(state, &iter);
}

bool ompl::base::PathLengthDirectAdInfSampler::hasInformedMeasure() const
{
    return true;
}

double ompl::base::PathLengthDirectAdInfSampler::getDirectInformedMeasure() const
{
    double informedMeasure = phsPtr_->getPhsrMeasure();
    if ((AdInformedSampler::space_->isCompound() && !compound_) || compoundCompound_)
        informedMeasure *= uninformedSubSpace_->getMeasure();
    return direct_ == 0 ? std::min(AdInformedSampler::space_->getMeasure(), informedMeasure) : informedMeasure;
}

double ompl::base::PathLengthDirectAdInfSampler::getInformedMeasure(double cost) const
{
    double informedMeasure = phsPtr_->getPhsMeasure(cost);
    if ((AdInformedSampler::space_->isCompound() && !compound_) || compoundCompound_)
        informedMeasure *= uninformedSubSpace_->getMeasure();
    return std::min(AdInformedSampler::space_->getMeasure(), informedMeasure);
}

double ompl::base::PathLengthDirectAdInfSampler::getInformedMeasure() const
{
    double informedMeasure = getDirectInformedMeasure();
    return std::min(AdInformedSampler::space_->getMeasure(), informedMeasure);
}

ompl::base::Cost ompl::base::PathLengthDirectAdInfSampler::heuristicSolnCost(const State *state) const
{
    std::vector<double> rawData = getInformedSubstate(state);

    return Cost(phsPtr_->getPathLength(&rawData[0]));
}

bool ompl::base::PathLengthDirectAdInfSampler::isInPhs(const State *state) const
{
    std::vector<double> data = getInformedSubstate(state);

    return phsPtr_->isInPhs(&data[0]);
}

bool ompl::base::PathLengthDirectAdInfSampler::isOutsideMiddlePhs(const State *state) const
{
    std::vector<double> data = getInformedSubstate(state);

    return phsPtr_->isOutsideMiddlePhs(&data[0]);
}

double ompl::base::PathLengthDirectAdInfSampler::getMinTransverseDiameter() const
{
    return phsPtr_->getMinTransverseDiameter();
}

bool ompl::base::PathLengthDirectAdInfSampler::sampleUniform(State *state, unsigned int *iters)
{
    bool foundSample = false;
    std::vector<double> informedVector(informedSubSpace_->getDimension());
    while (!foundSample && *iters < AdInformedSampler::numIters_)
    {
        rng_.uniformProlateHyperspheroidRing(phsPtr_, &informedVector[0]);
        if (!AdInformedSampler::space_->isCompound())
        {
            informedSubSpace_->copyFromReals(state, informedVector);
            if (informedSubSpace_->satisfiesBounds(state))
            {
                direct_++;
                foundSample = true;
            }
            else 
            {
                indirect_++;
                informedSubSpace_->enforceBoundsRandom(state);
                informedSubSpace_->copyToReals(informedVector, state);
                foundSample = isInPhs(informedVector);
            }
        }
        else if (!compoundCompound_)
        {
            informedSubSpace_->copyFromReals(state->as<CompoundState>()->components[informedIdx_], informedVector);
            if (informedSubSpace_->satisfiesBounds(state->as<CompoundState>()->components[informedIdx_]))
            {
                direct_++;
                foundSample = true;
            }
            else 
            {
                indirect_++;
                informedSubSpace_->enforceBoundsRandom(state->as<CompoundState>()->components[informedIdx_]);
                informedSubSpace_->copyToReals(informedVector, state->as<CompoundState>()->components[informedIdx_]);
                foundSample = isInPhs(informedVector);
            }
        }
        else
        {
            informedSubSpace_->copyFromReals(state->as<CompoundState>()->components[0]->as<CompoundState>()->components[informedIdx_], informedVector);
            if (informedSubSpace_->satisfiesBounds(state->as<CompoundState>()->components[0]->as<CompoundState>()->components[informedIdx_]))
            {
                direct_++;
                foundSample = true;
            }
            else 
            {
                indirect_++;
                informedSubSpace_->enforceBoundsRandom(state->as<CompoundState>()->components[0]->as<CompoundState>()->components[informedIdx_]);
                informedSubSpace_->copyToReals(informedVector, state->as<CompoundState>()->components[0]->as<CompoundState>()->components[informedIdx_]);
                foundSample = isInPhs(informedVector);
            }
        }

        if (foundSample)
        {
            if (AdInformedSampler::space_->isCompound() && !compound_)
            {
                State *uninformedState = uninformedSubSpace_->allocState();
                uninformedSubSampler_->sampleUniform(uninformedState);
                uninformedSubSpace_->copyState(state->as<CompoundState>()->components[uninformedIdx_], uninformedState);
                uninformedSubSpace_->freeState(uninformedState);
            }
            else if (compoundCompound_)
            {
                State *uninformedState = uninformedSubSpace_->allocState();
                uninformedSubSampler_->sampleUniform(uninformedState);
                uninformedSubSpace_->copyState(state->as<CompoundState>()->components[0]->as<CompoundState>()->components[uninformedIdx_], uninformedState);
                uninformedSubSpace_->freeState(uninformedState);
            }
        }
        ++(*iters);
    }
    return foundSample;
}

std::vector<double> ompl::base::PathLengthDirectAdInfSampler::getInformedSubstate(const State *state) const
{
    std::vector<double> rawData(informedSubSpace_->getDimension());
    if (!AdInformedSampler::space_->isCompound())
        informedSubSpace_->copyToReals(rawData, state);
    else if (!compoundCompound_)
        informedSubSpace_->copyToReals(rawData, state->as<CompoundState>()->components[informedIdx_]);
    else
        informedSubSpace_->copyToReals(rawData, state->as<CompoundState>()->components[0]->as<CompoundState>()->components[informedIdx_]);
    return rawData;
}

bool ompl::base::PathLengthDirectAdInfSampler::isInPhs(const std::vector<double> &informedVector) const
{
    return phsPtr_->isInPhs(&informedVector[0]);
}
