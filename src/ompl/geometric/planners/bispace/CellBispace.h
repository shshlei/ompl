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

/* Author: Shi Shenglei*/

#ifndef OMPL_GEOMETRIC_PLANNERS_BISPACE_CELLBISPACE_
#define OMPL_GEOMETRIC_PLANNERS_BISPACE_CELLBISPACE_

#include "ompl/geometric/planners/PlannerIncludes.h"
#include "ompl/geometric/planners/bispace/CellDiscretization.h"

#include "ompl/base/State.h"
#include "ompl/base/OptimizationObjective.h"
#include "ompl/base/ProjectionEvaluator.h"
#include "ompl/datastructures/BinaryHeap.h"

#include <unordered_set>

namespace ompl
{
    namespace geometric
    {
        /** \brief Cell Bispace */
        class CellBispace : public base::Planner
        {
        public:

            CellBispace(const base::SpaceInformationPtr &si);

            ~CellBispace() override = default;

            void setup() override;

            base::PlannerStatus solve(const base::PlannerTerminationCondition &ptc) override;

            void clear() override;

            void getPlannerData(base::PlannerData &data) const override;

            void setRange(double distance)
            {
                maxDistance_ = distance;
            }

            double getRange() const
            {
                return maxDistance_;
            }
            
            void setLazyNode(bool lazy)
            {
                lazyNode_ = lazy;
            }

            bool getLazyNode() const
            {
                return lazyNode_;
            }

            void setAddIntermediateState(bool add)
            {
                addIntermediateState_ = add;
            }

            bool getAddIntermediateState() const
            {
                return addIntermediateState_;
            }

            void setUseBispace(bool bispace)
            {
                useBispace_ = bispace;
            }

            bool getUseBispace() const
            {
                return useBispace_;
            }

            void setBackRewire(bool rewire)
            {
                backRewire_ = rewire;
            }

            bool getBackRewire() const
            {
                return backRewire_;
            }

            void setRewireSort(bool rewire)
            {
                rewireSort_ = rewire;
            }

            bool getRewireSort() const
            {
                return rewireSort_;
            }

        protected:

            ////////////////////////////////////////////////////////////////////////////////////////////////////
            // general functions
            enum StateValid
            {
                UnCkecked,
                InValid,
                Valid 
            };

            class Motion;
            using CellDiscretizationData = CellDiscretization<Motion>; 
            using Cell = CellDiscretizationData::Cell;
            using Coord = CellDiscretizationData::Coord;

            struct MotionCompare
            {
                MotionCompare(base::OptimizationObjectivePtr opt) : opt_(std::move(opt))
                {
                }

                /** \brief Ordering of motions */
                inline bool operator()(const Motion *motion1, const Motion *motion2) const
                {
                    return opt_->isCostBetterThan(motion1->cost, motion2->cost);
                }

                /** \brief Pointer to the Optimization Objective */
                base::OptimizationObjectivePtr opt_;
            };

            struct CostMotionCompare
            {
                CostMotionCompare(Motion *motion, base::OptimizationObjectivePtr opt, bool start)
                  : motion_(motion), opt_(std::move(opt)), start_(start)
                {
                }
                inline bool operator()(Motion * motion1, Motion * motion2)
                {
                    base::Cost cost1 = motion1->cost;
                    base::Cost cost2 = motion2->cost;

                    if (start_)
                    {
                        cost1 = opt_->combineCosts(cost1, opt_->motionCost(motion1->state, motion_->state));
                        cost2 = opt_->combineCosts(cost2, opt_->motionCost(motion2->state, motion_->state));
                    }
                    else 
                    {
                        cost1 = opt_->combineCosts(cost1, opt_->motionCost(motion_->state, motion1->state));
                        cost2 = opt_->combineCosts(cost2, opt_->motionCost(motion_->state, motion2->state));
                    }

                    return opt_->isCostBetterThan(cost1, cost2);
                }
                Motion *motion_;
                base::OptimizationObjectivePtr opt_;
                bool start_;
            };

            struct RewireSort 
            {
                RewireSort()
                {
                }
                inline bool operator()(Motion * motion1, Motion * motion2)
                {
                    if (motion1->valid && motion2->valid)
                    {
                        std::size_t len1 = 0, len2 = 0;
                        Motion *motion = motion1;
                        while (motion)
                        {
                            len1++;
                            motion = motion->parent;
                        }
                        motion = motion2;
                        while (motion)
                        {
                            len2++;
                            motion = motion->parent;
                        }
                        return len1 < len2;
                    }
                    else
                        return motion2->valid;
                }
            };

            struct OrderCellsByDisabled
            {
                /** \brief Order function */
                bool operator()(const Cell *const a, const Cell *const b) const
                {
                    return a->data->disabled < b->data->disabled;
                }
            };

            struct OrderCellsByCost
            {
                OrderCellsByCost(base::OptimizationObjectivePtr opt) : opt_(std::move(opt))
                {
                }

                /** \brief Order function */
                bool operator()(const Cell *const a, const Cell *const b) const
                {
                    if (a->data->cmotion && b->data->cmotion)
                        return opt_->isCostBetterThan(b->data->cmotion->cost, a->data->cmotion->cost);
                    else
                        return b->data->cmotion;
                }

                /** \brief Pointer to the Optimization Objective */
                base::OptimizationObjectivePtr opt_;
            };

            /** \brief Representation of a motion */
            class Motion
            {
            public:
                Motion(const base::SpaceInformationPtr &si) : state(si->allocState())
                {
                }

                ~Motion() = default;

                /** \brief The state contained by the motion */
                base::State *state{nullptr};

                const base::State *root{nullptr};

                /** \brief The parent motion in the exploration tree */
                Motion *parent{nullptr};

                Cell *cell{nullptr};

                StateValid stateValid{UnCkecked};

                bool valid{false};

                bool inConnection{false};

                /** \brief The cost up to this motion */
                base::Cost cost;

                /** \brief The incremental cost of this motion's parent to this motion (this is stored to save distance
                 * computations in the updateChildCosts() method) */
                base::Cost incCost;

                /** \brief The set of motions descending from the current motion */
                std::vector<Motion *> children;

                /** \brief The valid neighborhood */
                std::unordered_set<Motion *> nbh;

                /** \brief The invalid neighborhood */
                std::unordered_set<Motion *> invalidnbh;

                /** \brief Handle to identify the motion in the queue */
                BinaryHeap<Motion *, MotionCompare>::Element *handle{nullptr};
            };

            /** \brief Information attached to growing a tree of motions (used internally) */
            struct TreeGrowingInfo
            {
                base::State *xstate;
                Motion *xmotion;
                bool start;
            };

            base::PlannerStatus prepareSolve(const base::PlannerTerminationCondition &ptc);

            void processSolution(const Motion *bestStartMotion, const Motion *bestGoalMotion);

            bool batchGrow(bool &startTree);

            /** feasible */
            bool growTree(TreeGrowingInfo &tgi);

            Motion* backGrowTree(Motion *root, base::State *dstate, bool start);

            void rewireTree(BinaryHeap<Motion *, MotionCompare> &bh, Motion *m, bool start);

            void updateQueue(BinaryHeap<Motion *, MotionCompare> &bh, Motion *m);

            bool growCurrentTree(const base::State *state, bool start) const;

            bool backPathRewireMotionSimple(Motion *motion, bool start, Motion *&pmotion);

            /** \brief Check if the connected path is valid */
            bool isPathValid(Motion *motion, Motion *otherMotion);

            /** \brief Check from the root to the connect point, stop immediately if an invalid path is found  */
            bool isPathValid(Motion *motion, bool start);

            /** \brief Check from the connect point to the root */
            bool isPathValidInter(Motion *motion, bool start);

            /** \brief Check from the connect point to the root */
            bool isStateValid(Motion *motion, bool start);

            bool backPathRewireMotion(Motion *motion, bool start, Motion *&pmotion);

            void removeInvalidMotions();

            void removeInvalidMotionsDisc();

            void enableMotionInDisc(Motion *motion);

            void enableMotionInDiscWithRewireSort(Motion *motion);

            void disableCheck(Motion *motion, const std::string &context) const;

            void connectToPmotion(Motion *motion, Motion *pmotion, bool start) const;

            /** \brief Updates the cost of the children of this node if the cost up to this node has changed */
            void updateChildCosts(Motion *motion) const;

            void insertNeighbor(Motion *pmotion, Motion *motion);

            void insertInvalidNeighbor(Motion *pmotion, Motion *motion);

            void removeFromNeighbor(Motion *motion);

            void removeFromInvalidNeighbor(Motion *motion);

            bool isValidNeighbor(Motion *motion, Motion *nb) const;

            bool isInvalidNeighbor(Motion *motion, Motion *nb) const;

            /** \brief Removes the given motion from the parent's child list */
            void removeFromParent(Motion *motion);

            void removeFromPnull(std::vector<Motion *> &pnullMotions, Motion *motion);

            bool removeFromVector(std::vector<Motion *> &motions, Motion *motion);

            void setMotionInfinityCost(Motion *motion) const;

            void setMotionInfinityCostWithRewireSort(Motion *motion, std::unordered_set<Cell *> &cells) const;

            // check 
            bool isValid(const base::State *state);

            // check
            bool isValid(Motion *motion, bool start);

            bool checkMotion(Motion *pmotion, Motion *motion, bool start);

            bool checkInterMotion(Motion *pmotion, Motion *motion, bool start);

            bool checkInterMotion1(Motion *smotion, Motion *gmotion, bool start);

            bool checkInterMotion2(Motion *smotion, Motion *gmotion, bool start);

            void addIntermediateMotion(Motion *pmotion, Motion *motion, bool start, Motion *last);

            base::ProjectionEvaluatorPtr projectionEvaluator_;

            /** \brief The start tree */
            CellDiscretizationData dStart_;

            /** \brief The goal tree */
            CellDiscretizationData dGoal_;

            std::vector<Motion *> pnullStartMotions_;

            std::vector<Motion *> pnullGoalMotions_;

            /** \brief Stores the start states as Motions. */
            std::vector<Motion *> startMotions_;

            /** \brief A list of states in the tree that satisfy the goal condition */
            std::vector<Motion *> goalMotions_;

            std::vector<Motion *> invalidStartMotions_;

            std::vector<Motion *> invalidGoalMotions_;

            /** \brief The pair of motions in each tree connected during planning. */
            std::vector<std::pair<Motion *, Motion *>> connectionPoint_;

            Coord scoord_;
            Eigen::VectorXd dcoord_;
            double dc_;

            double maxDistance_{0.0};

            bool lazyNode_{true};

            bool addIntermediateState_{true};

            bool backRewire_{true};

            bool rewireSort_{true};

            bool symmetric_{true};

            /** \brief Objective we're optimizing */
            base::OptimizationObjectivePtr opt_;

            MotionCompare mc_;

            BinaryHeap<Motion *, MotionCompare> bh_;

            void freeMotion(Motion *motion)
            {
                if (motion->state != nullptr)
                    si_->freeState(motion->state);
                delete motion;
            }

            base::StateSamplerPtr sampler_;

            RNG rng_;

            ////////////////////////////////////////////////////////////////////////////////////////////////////
            // bispace 
            bool useBispace_{true};
        };
    }
}

#endif
