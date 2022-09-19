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

#include "ompl/base/samplers/AdInformedStateSampler.h"
#include "ompl/base/samplers/adinformed/RejectionAdInfSampler.h"
#include "ompl/util/Exception.h"
#include "ompl/base/objectives/PathLengthOptimizationObjective.h"
#include "ompl/base/Goal.h"

ompl::base::AdInformedSampler::AdInformedSampler(const ProblemDefinitionPtr &probDefn, const State *s1, const State *s2,
                                                 const Cost &minCost, const Cost &maxCost, unsigned int maxNumberCalls)
  : probDefn_(probDefn), space_(probDefn->getSpaceInformation()->getStateSpace()), 
    numIters_(maxNumberCalls), minCost_(minCost), maxCost_(maxCost)
{
    if (!probDefn_->hasOptimizationObjective())
        opt_ = std::make_shared<PathLengthOptimizationObjective>(probDefn->getSpaceInformation());
    else 
        opt_ = probDefn_->getOptimizationObjective();

    s1_ = space_->cloneState(s1);
    s2_ = space_->cloneState(s2);
}

void ompl::base::AdInformedSampler::update(const Cost &minCost, const Cost &maxCost)
{
    minCost_ = minCost;
    maxCost_ = maxCost;
}

void ompl::base::AdInformedSampler::update(const Cost &maxCost)
{
    minCost_ = maxCost_;
    maxCost_ = maxCost; 
}

void ompl::base::AdInformedSampler::update(double factor)
{
    minCost_ = maxCost_;
    maxCost_ = base::Cost(factor * maxCost_.value());
}

void ompl::base::AdInformedSampler::clear()
{
    space_->freeState(s1_);
    space_->freeState(s2_);
}

ompl::base::Cost ompl::base::AdInformedSampler::heuristicSolnCost(const State *state) const
{
    return opt_->combineCosts(opt_->motionCostHeuristic(s1_, state), opt_->motionCostHeuristic(state, s2_));
}

ompl::base::ProblemDefinitionPtr ompl::base::AdInformedSampler::getProblemDefn() const
{
    return probDefn_;
}

unsigned int ompl::base::AdInformedSampler::getMaxNumberOfIters() const
{
    return numIters_;
}

ompl::base::AdInformedStateSampler::AdInformedStateSampler(const ProblemDefinitionPtr &probDefn, 
                                                           const State *s1, const State *s2,
                                                           const Cost &minCost, const Cost &maxCost,
                                                           unsigned int maxNumberCalls)
  : StateSampler(probDefn->getSpaceInformation()->getStateSpace().get())
{
    infSampler_ = std::make_shared<RejectionAdInfSampler>(probDefn, s1, s2, minCost, maxCost, maxNumberCalls);

    baseSampler_ = StateSampler::space_->allocStateSampler();
}

void ompl::base::AdInformedStateSampler::update(const Cost &minCost, const Cost &maxCost)
{
    infSampler_->update(minCost, maxCost);
}

void ompl::base::AdInformedStateSampler::update(const Cost &maxCost)
{
    infSampler_->update(maxCost);
}

void ompl::base::AdInformedStateSampler::update(double factor)
{
    infSampler_->update(factor);
}

void ompl::base::AdInformedStateSampler::sampleUniform(State *state)
{
    bool informedSuccess;

    informedSuccess = infSampler_->sampleUniform(state);

    if (!informedSuccess)
    {
        baseSampler_->sampleUniform(state);
    }
}

void ompl::base::AdInformedStateSampler::sampleUniformNear(State *state, const State *near, const double distance)
{
    OMPL_WARN("sampleUniformNear is not informed.");
    return baseSampler_->sampleUniformNear(state, near, distance);
}

void ompl::base::AdInformedStateSampler::sampleGaussian(State *state, const State *mean, const double stdDev)
{
    OMPL_WARN("sampleGaussian is not informed.");
    return baseSampler_->sampleGaussian(state, mean, stdDev);
}
