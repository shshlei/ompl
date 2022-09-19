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

#include "ompl/base/samplers/adinformed/RejectionAdInfSampler.h"
#include "ompl/base/OptimizationObjective.h"

ompl::base::RejectionAdInfSampler::RejectionAdInfSampler(const ProblemDefinitionPtr &probDefn,
                                                         const State *s1, const State *s2,
                                                         const Cost &minCost, const Cost &maxCost,
                                                         unsigned int maxNumberCalls)
  : AdInformedSampler(probDefn, s1, s2, minCost, maxCost, maxNumberCalls)
{
    baseSampler_ = AdInformedSampler::space_->allocStateSampler();
}

bool ompl::base::RejectionAdInfSampler::sampleUniform(State *state)
{
    return sampleUniform(state, this->minCost_, this->maxCost_);
}

bool ompl::base::RejectionAdInfSampler::sampleUniform(State *state, const Cost &minCost, const Cost &maxCost)
{
    bool foundSample = false;

    for (unsigned int i = 0u; i < AdInformedSampler::numIters_ && !foundSample; ++i)
    {
        foundSample = sampleUniform(state, maxCost, &i);

        if (foundSample)
        {
            Cost sampledCost = AdInformedSampler::heuristicSolnCost(state);

            foundSample = this->opt_->isCostBetterThan(minCost, sampledCost);
        }
    }

    return foundSample;
}

bool ompl::base::RejectionAdInfSampler::hasInformedMeasure() const
{
    return false;
}

double ompl::base::RejectionAdInfSampler::getInformedMeasure() const
{
    return AdInformedSampler::space_->getMeasure();
}

bool ompl::base::RejectionAdInfSampler::sampleUniform(State *state, const Cost &maxCost, unsigned int *iterPtr)
{
    bool foundSample = false;

    for (; *iterPtr < AdInformedSampler::numIters_ && !foundSample; ++(*iterPtr))
    {
        baseSampler_->sampleUniform(state);

        foundSample = this->opt_->isCostBetterThan(AdInformedSampler::heuristicSolnCost(state), maxCost);
    }

    return foundSample;
}
