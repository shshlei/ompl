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

#ifndef OMPL_BASE_SAMPLERS_ADINFORMED_GENERAL_REJECTION_INFORMED_SAMPLER_
#define OMPL_BASE_SAMPLERS_ADINFORMED_GENERAL_REJECTION_INFORMED_SAMPLER_

#include "ompl/base/samplers/AdInformedStateSampler.h"

namespace ompl
{
    namespace base
    {
        class RejectionAdInfSampler : public AdInformedSampler
        {
        public:
            RejectionAdInfSampler(const ProblemDefinitionPtr &probDefn, const State *s1, const State *s2,
                                  const Cost &minCost, const Cost &maxCost, unsigned int maxNumberCalls);
            ~RejectionAdInfSampler() = default;

            bool sampleUniform(State *state) override;

            bool sampleUniform(State *state, const Cost &minCost, const Cost &maxCost);

            bool hasInformedMeasure() const override;

            double getInformedMeasure() const override;

        private:
            StateSamplerPtr baseSampler_;

            bool sampleUniform(State *state, const Cost &maxCost, unsigned int *iterPtr);
        };
    }
}

#endif  // OMPL_BASE_SAMPLERS_ADINFORMED_REJECTION_INFORMED_SAMPLER_
