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

#ifndef OMPL_BASE_SAMPLERS_ADINFORMED_SAMPLER_
#define OMPL_BASE_SAMPLERS_ADINFORMED_SAMPLER_

#include "ompl/base/StateSampler.h"
#include "ompl/base/Cost.h"
#include "ompl/base/ProblemDefinition.h"

namespace ompl
{
    namespace base
    {
        OMPL_CLASS_FORWARD(AdInformedSampler);
        OMPL_CLASS_FORWARD(AdInformedStateSampler);

        /** \brief An abstract class for the concept of using information about the state space
        and the current infeasible solution cost to limit future search to a planning
        subproblem that contains all possibly better solutions. */
        class AdInformedSampler
        {
        public:
            AdInformedSampler(const ProblemDefinitionPtr &probDefn, const State *s1, const State *s2,
                              const Cost &minCost, const Cost &maxCost, unsigned int maxNumberCalls);

            ~AdInformedSampler()
            {
                clear();
            }

            // non-copyable
            AdInformedSampler(const AdInformedSampler &) = delete;
            AdInformedSampler &operator=(const AdInformedSampler &) = delete;

            virtual void update(const Cost &minCost, const Cost &maxCost);

            virtual void update(const Cost &maxCost);

            virtual void update(double factor);

            void clear();

            virtual bool sampleUniform(State *state) = 0;

            /** \brief Whether the sampler can provide a measure of the informed subset */
            virtual bool hasInformedMeasure() const = 0;

            /** \brief The measure of the subset of the state space defined by the current solution cost that is being
             * searched. Does not consider problem boundaries but returns the measure of the entire space if no solution
             * has been found or if a closed form expression for the measure does not exist. */
            virtual double getInformedMeasure() const = 0;

            /** \brief A helper function to calculate the heuristic estimate of the solution cost for a given state
             * using the optimization objective stored in the problem definition. */
            /** \todo With the future invention of a heuristic class, this should move.  */
            virtual Cost heuristicSolnCost(const State *state) const;

            /** Helper for the OrderedInfSampler wrapper */
            ProblemDefinitionPtr getProblemDefn() const;

            /** Helper for the OrderedInfSampler wrapper */
            unsigned int getMaxNumberOfIters() const;

            Cost getMinCost() const 
            {
                return minCost_;
            }

            Cost getMaxCost() const 
            {
                return maxCost_;
            }

            /** \brief Cast this instance to a desired type. */
            template <class T>
            const T *as() const
            {
                BOOST_CONCEPT_ASSERT((boost::Convertible<T *, AdInformedSampler *>));

                return static_cast<const T *>(this);
            }

            /** \brief Cast this instance to a desired type. */
            template <class T>
            T *as()
            {
                BOOST_CONCEPT_ASSERT((boost::Convertible<T *, AdInformedSampler *>));

                return static_cast<T *>(this);
            }

        protected:
            /** \brief A copy of the problem definition */
            ProblemDefinitionPtr probDefn_;
            /** \brief A copy of the optimization objective */
            OptimizationObjectivePtr opt_;
            /** \brief A copy of the state space*/
            StateSpacePtr space_;
            /** \brief The number of iterations I'm allowed to attempt */
            unsigned int numIters_;
            /** \brief The start and end states of the planning path*/
            State *s1_, *s2_;
            Cost minCost_, maxCost_;
        };

        /** \brief A wrapper class that allows an AdInformedSampler to be used as a StateSampler. */
        class AdInformedStateSampler : public StateSampler
        {
        public:
            AdInformedStateSampler(const ProblemDefinitionPtr &probDefn, const State *s1, const State *s2,
                                   const Cost &minCost, const Cost &maxCost, unsigned int maxNumberCalls);

            ~AdInformedStateSampler() override = default;

            void update(const Cost &minCost, const Cost &maxCost);

            void update(const Cost &maxCost);

            void update(double factor);

            void sampleUniform(State *state) override;

            void sampleUniformNear(State *state, const State *near, double distance) override;

            void sampleGaussian(State *state, const State *mean, double stdDev) override;

        private:
            /** \brief A basic sampler */
            StateSamplerPtr baseSampler_;
            /** \brief The wrapped informed sampler */
            AdInformedSamplerPtr infSampler_;
        };
    }
}

#endif  // OMPL_BASE_SAMPLERS_INFORMED_SAMPLER_
