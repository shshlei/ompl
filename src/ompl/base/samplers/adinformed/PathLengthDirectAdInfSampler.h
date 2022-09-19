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

#ifndef OMPL_BASE_SAMPLERS_ADINFORMED_PATH_LENGTH_DIRECT_INFORMED_SAMPLER_
#define OMPL_BASE_SAMPLERS_ADINFORMED_PATH_LENGTH_DIRECT_INFORMED_SAMPLER_

#include "ompl/base/samplers/AdInformedStateSampler.h"

namespace ompl
{
    namespace base
    {
        class PathLengthDirectAdInfSampler : public AdInformedSampler
        {
        public:
            PathLengthDirectAdInfSampler(const ProblemDefinitionPtr &probDefn, const State *s1, const State *s2,
                                         const Cost &minCost, const Cost &maxCost, unsigned int maxNumberCalls);

            ~PathLengthDirectAdInfSampler() = default;

            void update(const Cost &minCost, const Cost &maxCost) override;

            void update(const Cost &maxCost) override;

            void update(double factor) override;

            bool sampleUniform(State *state) override;

            bool hasInformedMeasure() const override;

            double getDirectInformedMeasure() const;

            double getInformedMeasure(double cost) const;

            double getInformedMeasure() const override;

            Cost heuristicSolnCost(const State *state) const override;

            unsigned int getDirectSampingCount() const 
            {
                return direct_;
            }

            unsigned int getInDirectSamplingCount() const 
            {
                return indirect_;
            }

            double getDirectSamplingFraction() const
            {
                return direct_ == 0 ? 1.0 : (double)direct_ / (double)(direct_ + indirect_);
            }

            void resetDirectCounter() 
            {
                direct_ = indirect_ = 0;
            }

            unsigned int getInformedDimension() const
            {
                return informedSubSpace_->getDimension();
            } 

            bool isInPhs(const State *state) const;

            bool isOutsideMiddlePhs(const State *state) const;

            double getMinTransverseDiameter() const;

        private:

            using ProlateHyperspheroidRingCPtr = std::shared_ptr<const ompl::ProlateHyperspheroidRing>;

            bool sampleUniform(State *state, unsigned int *iters);

            std::vector<double> getInformedSubstate(const State *state) const;

            void createFullState(State *state, const std::vector<double> &informedVector);

            bool isInPhs(const std::vector<double> &informedVector) const;

            void enforceBoundsRandom(std::vector<double> &informedVector);

            ompl::ProlateHyperspheroidRingPtr phsPtr_;

            unsigned int informedIdx_;

            StateSpacePtr informedSubSpace_;

            unsigned int uninformedIdx_;

            StateSpacePtr uninformedSubSpace_;

            StateSamplerPtr baseSampler_;

            StateSamplerPtr uninformedSubSampler_;

            RNG rng_;

            mutable unsigned int direct_;

            mutable unsigned int indirect_;
        };
    }
}

#endif  // OMPL_BASE_SAMPLERS_ADINFORMED_DIRECT_PATH_LENGTH_INFORMED_SAMPLER_
