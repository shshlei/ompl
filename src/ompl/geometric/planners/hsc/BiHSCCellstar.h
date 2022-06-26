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

#ifndef OMPL_GEOMETRIC_PLANNERS_HSC_BIHSCCELLSTAR_
#define OMPL_GEOMETRIC_PLANNERS_HSC_BIHSCCELLSTAR_

#include "ompl/geometric/planners/PlannerIncludes.h"
#include "ompl/geometric/planners/hsc/CellDiscretizationN.h"

#include "ompl/datastructures/NearestNeighbors.h"

#include "ompl/base/State.h"
#include "ompl/base/SafetyCertificate.h"
#include "ompl/base/OptimizationObjective.h"
#include "ompl/base/ProjectionEvaluator.h"

#include <unordered_set>
#include <queue>
#include <deque>
#include <utility>
#include <list>

namespace ompl
{
    namespace geometric
    {
        /** \brief Optimal Bidirectional Hybrid Safety Certificate Using Cell Exploration */
        class BiHSCCellstar : public base::Planner
        {
        public:
            /** \brief Check if the state is certificated as collision-free by sc, and return the certificate distances */
            using SafetyCertificateChecker = std::function<bool(const base::State *state, const base::SafetyCertificate *sc, std::vector<double> &dist)>;

            using CollisionCertificateChecker = std::function<bool(const base::State *state, const std::vector<base::SafetyCertificate *> &ocv)>;

            /** \brief Calculate the certificate distances between two states */
            using DistanceCertificate = std::function<std::vector<double>(const base::State *a, const base::State *b)>;

            BiHSCCellstar(const base::SpaceInformationPtr &si);

            ~BiHSCCellstar() override;

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
            
            void setPenDistance(double distance)
            {
                penDistance_ = distance;
            }

            double getPenDistance() const
            {
                return penDistance_;
            }

            void setLazyPath(bool lazy)
            {
                lazyPath_ = lazy;
                if (!lazyPath_)
                    lazyNode_ = false;
            }

            bool getLazyPath() const
            {
                return lazyPath_;
            }

            void setLazyNode(bool lazy)
            {
                if (lazy)
                    lazyPath_ = true;
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

            /** \brief When extending a motion, the planner can decide
                to keep the first valid part of it, even if invalid
                states are found, as long as the valid part represents
                a sufficiently large fraction from the original
                motion. This function sets the minimum acceptable
                fraction. */
            void setMinValidPathFraction(double fraction)
            {
                minValidPathFraction_ = fraction;
            }

            /** \brief Get the value of the fraction set by setMinValidPathFraction() */
            double getMinValidPathFraction() const
            {
                return minValidPathFraction_;
            }

            void setUseBispace(bool bispace)
            {
                if (!bispace)
                {
                    treatedAsMultiSubapce_ = false;
                }
                useBispace_ = bispace;
            }

            bool getUseBispace() const
            {
                return useBispace_;
            }

            void setTreatedAsMultiSubapce(bool sub)
            {
                treatedAsMultiSubapce_ = sub;
                if (!useBispace_ && sub)
                {
                    OMPL_WARN("%s: The treate as multisubspace feature is not usable since the bispace feature is not set.", getName().c_str());
                    treatedAsMultiSubapce_ = false;
                }
                if (treatedAsMultiSubapce_)
                {
                    OMPL_WARN("%s: The treate as multisubspace feature is implemented in combination with the feature cell CellDiscretization.", getName().c_str());
                    treatedAsMultiSubapce_ = false;
                }
            }

            bool getTreatedAsMultiSubapce() const
            {
                return treatedAsMultiSubapce_;
            }

            /** \brief Set the projection evaluator. This class is
                able to compute the projection of a given state. */
            void setProjectionEvaluator(const base::ProjectionEvaluatorPtr &projectionEvaluator)
            {
                projectionEvaluator_ = projectionEvaluator;
            }

            /** \brief Set the projection evaluator (select one from
                the ones registered with the state space). */
            void setProjectionEvaluator(const std::string &name)
            {
                projectionEvaluator_ = si_->getStateSpace()->getProjection(name);
            }

            /** \brief Get the projection evaluator. */
            const base::ProjectionEvaluatorPtr &getProjectionEvaluator() const
            {
                return projectionEvaluator_;
            }

            /** \brief Set the fraction of time for focusing on the
                border (between 0 and 1). This is the minimum fraction
                used to select cells that are exterior (minimum
                because if 95% of cells are on the border, they will
                be selected with 95% chance, even if this fraction is
                set to 90%)*/
            void setBorderFraction(double bp)
            {
                dStart_.setBorderFraction(bp);
                dGoal_.setBorderFraction(bp);
            }

            /** \brief Get the fraction of time to focus exploration
                on boundary */
            double getBorderFraction() const
            {
                return dStart_.getBorderFraction();
            }

            void setPruneThreshold(const double pp)
            {
                pruneThreshold_ = pp;
            }

            /** \brief Get the current prune states percentage threshold parameter. */
            double getPruneThreshold() const
            {
                return pruneThreshold_;
            }

            unsigned int numIterations() const
            {
                return iterations_;
            }

            ompl::base::Cost bestCost() const
            {
                return bestCost_;
            }

            void setStrictCertificate(bool strictCertificate)
            {
                strictCertificate_ = strictCertificate;
            }

            bool getStrictCertificate() const
            {
                return strictCertificate_;
            }

            void setUseCollisionCertificateChecker(bool use)
            {
                useCollisionCertificateChecker_ = use;
            }

            bool getUseCollisionCertificateChecker() const
            {
                return useCollisionCertificateChecker_;
            }

            void setSafetyCertificateChecker(const SafetyCertificateChecker &safetyCertificateChecker)
            {
                safetyCertificateChecker_ = safetyCertificateChecker;
            }

            void setCollisionCertificateChecker(const CollisionCertificateChecker &collisionCertificateChecker)
            {
                collisionCertificateChecker_ = collisionCertificateChecker;
            }

            void setDistanceCertificate(const DistanceCertificate &distanceCertificate)
            {
                distanceCertificate_ = distanceCertificate;
            }

            void setCertificateDim(unsigned int certificateDim)
            {
                certificateDim_ = certificateDim;
            }

            unsigned int getCertificateDim() const 
            {
                return certificateDim_;
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
            using CellDiscretizationData = CellDiscretizationN<Motion>; 
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
                    return 100 * a->auxData->root + a->auxData->disabled < 100 * b->auxData->root + b->auxData->disabled;
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
                    if (a->auxData->cmotion && b->auxData->cmotion)
                        return opt_->isCostBetterThan(b->auxData->cmotion->cost, a->auxData->cmotion->cost);
                    else
                        return b->auxData->cmotion;
                }

                /** \brief Pointer to the Optimization Objective */
                base::OptimizationObjectivePtr opt_;
            };

            struct SafetyCertificateWithElems 
            {
                base::SafetyCertificate *sc;
                std::vector<Motion *> objects;
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

                SafetyCertificateWithElems *sce{nullptr};

                /** \brief The certificate distance between the motion and its certificate node */
                std::vector<double> scd;
            };

            /** \brief Information attached to growing a tree of motions (used internally) */
            struct TreeGrowingInfo
            {
                base::State *xstate;
                Motion *xmotion;
                bool start;
            };

            void processSolution(const Motion *bestStartMotion, const Motion *bestGoalMotion);

            bool batchGrow(bool &startTree);

            bool growCurrentTree(const base::State *state) const;

            bool growStartTree(const base::State *state) const;

            double penetrationDistance(const base::State *nstate, const base::State *state, bool start) const;

            bool growCurrentTree(const base::State *state, std::size_t sub) const;

            bool growStartTree(const base::State *state, std::size_t sub) const;

            double penetrationDistance(const base::State *nstate, const base::State *state, bool start, std::size_t sub) const;

            void rewireTree(BinaryHeap<Motion *, MotionCompare> &bh, Motion *m, bool start);

            void optimalRewireTree(BinaryHeap<Motion *, MotionCompare> &bh, Motion *m, bool start);

            /** \brief If successively connected local paths are all valid, regard the straight line connecting 
                the endpoints is valid*/
            bool backRewire(Motion *motion, Motion *nb, bool start, Motion *&pmotion, base::Cost &nbhIncCost, base::Cost &nbhNewCost, bool &feas);

            void updateQueue(BinaryHeap<Motion *, MotionCompare> &bh, Motion *m);

            void updateLeafQueue(BinaryHeap<Motion *, MotionCompare> &bh, Motion *motion);

            /** \brief Check if the connected path is valid */
            bool isPathValid(Motion *motion, Motion *otherMotion);

            /** \brief Check from the root to the connect point, stop immediately if an invalid path is found  */
            bool isPathValid(Motion *motion, bool start);

            /** \brief Check from the connect point to the root */
            bool isStateValid(Motion *motion, bool start);

            /** \brief Check from the connect point to the root */
            bool isPathValidInter(Motion *motion, bool start);

            void removeInvalidMotions();

            void removeInvalidMotionsDisc();

            void enableMotionInDisc(CellDiscretizationData &disc, Motion *motion);

            void enableMotionInDisc2(Motion *motion, std::unordered_set<Cell *> &cells);

            void connectToPmotion(Motion *motion, Motion *pmotion, bool start) const;

            /** \brief Updates the cost of the children of this node if the cost up to this node has changed */
            void updateChildCosts(Motion *motion) const;

            bool checkStartMotion(Motion *smotion, Motion *gmotion);

            bool checkGoalMotion(Motion *smotion, Motion *gmotion);

            bool checkInterMotion(Motion *smotion, Motion *gmotion, bool start);

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

            void setMotionInfinityCostWithDisable(Motion *motion, std::unordered_set<Cell *> &cells) const;

            void getPlannerData(base::PlannerData &data, Motion *motion, int tag) const;

            base::ProjectionEvaluatorPtr projectionEvaluator_;

            /** \brief The start tree */
            CellDiscretizationData dStart_;

            /** \brief The goal tree */
            CellDiscretizationData dGoal_;

            std::vector<Motion *> pnullStartMotions_;

            std::vector<Motion *> pnullGoalMotions_;

            std::vector<Motion *> checkedStartPath_;

            std::vector<Motion *> checkedGoalPath_;

            /** \brief Stores the start states as Motions. */
            std::vector<Motion *> startMotions_;

            /** \brief A list of states in the tree that satisfy the goal condition */
            std::vector<Motion *> goalMotions_;

            std::vector<Motion *> invalidStartMotions_;

            std::vector<Motion *> invalidGoalMotions_;

            /** \brief The pair of motions in each tree connected during planning. */
            std::vector<std::pair<Motion *, Motion *>> connectionPoint_;

            double maxDistance_{0.0};

            double penDistance_{0.0};

            /** \brief When extending a motion, the planner can decide
                to keep the first valid part of it, even if invalid
                states are found, as long as the valid part represents
                a sufficiently large fraction from the original
                motion */
            double minValidPathFraction_{0.5};

            bool lazyPath_{true};

            bool lazyNode_{true};

            bool addIntermediateState_{true};

            bool solved_{false};

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
            // optimal 
            /** \brief Check if a better path is found */
            bool findBetterSolution(bool &optimal, double &ratio1, double &maxratio1, unsigned int &connect1, double &connectTresh1);

            void rewirePath();

            void reportBetterSolution(const base::ReportIntermediateSolutionFn &intermediateSolutionCallback);

            /** \brief Determine an improvement ratio that a found path's feasibility is rechecked */
            double improvementRatio(const base::Cost &temp, const base::State *sm, const base::State *gm) const;

            /** \brief If lazy strategy is used, the feasibility of a new found path is only checked
                if it is ratio1 better than the stored best path, or after connectTresh1 connection times */
            bool checkPath(const base::Cost &temp, const base::Cost &best,
                           double &ratio1, double &maxratio1, unsigned int &connect1, double &connectTresh1) const;

            /** \brief Prunes all those states which estimated total cost is higher than pruneTreeCost.
                Returns the number of motions pruned. Depends on the parameter set by
               setPruneStatesImprovementThreshold() */
            std::size_t pruneTree(const base::Cost &pruneTreeCost);

            std::size_t pruneSingleTree(CellDiscretizationData &disc, const base::Cost &pruneTreeCost, bool start, const std::vector<Motion *> &rootMotions, std::unordered_set<Cell *> &cells);

            std::size_t pruneTreeInternal(CellDiscretizationData &disc, const base::Cost &pruneTreeCost, bool start,
                                std::queue<Motion *, std::deque<Motion *>> &motionQueue, std::unordered_set<Cell *> &cells);

            std::size_t pruneTreeInternalDisabled(CellDiscretizationData &disc, const base::Cost &pruneTreeCost, bool start,
                                std::queue<Motion *, std::deque<Motion *>> &motionQueue, bool &invalid, std::unordered_set<Cell *> &cells);

            void toPrune(std::queue<Motion *, std::deque<Motion *>> &motionQueue,
                    std::queue<Motion *, std::deque<Motion *>> &leavesToPrune, std::list<Motion *> &chainsToRecheck, const base::Cost &pruneTreeCost, bool start);

            void pruneMotion(Motion *motion, CellDiscretizationData &disc, bool start, std::unordered_set<Cell *> &cells);

            /** \brief Add the children of a vertex to the given list. */
            void addChildrenToList(std::queue<Motion *, std::deque<Motion *>> *motionList, Motion *motion);

            /** \brief Check whether the given motion passes the specified cost threshold, meaning it will be \e kept
             * during pruning */
            bool keepCondition(Motion *motion, const base::Cost &threshold, bool start);

            /** \brief Calculate the path cost to reach the bispace border */
            base::Cost bordersolutionHeuristic(Motion *motion, bool start) const;

            bool keepCondition2(Motion *motion, const base::Cost &threshold, bool start) const;

            base::Cost solutionHeuristic2(Motion *motion, bool start) const;

            base::Cost calculateCostToCome(Motion *motion, bool start) const;

            base::Cost calculateCostToGo(Motion *motion, bool start) const;

            Motion *bestStartMotion_{nullptr};

            Motion *bestGoalMotion_{nullptr};

            /** \brief Best cost found so far by algorithm */
            base::Cost bestCost_{std::numeric_limits<double>::quiet_NaN()};

            /** \brief The cost at which the graph was last pruned */
            base::Cost prunedCost_{std::numeric_limits<double>::quiet_NaN()};

            /** \brief The tree is pruned when the change in solution cost is greater than this fraction. */
            double pruneThreshold_{0.05};

            unsigned int iterations_{0u};

            std::string numIterationsProperty() const
            {
                return std::to_string(numIterations());
            }
            std::string bestCostProperty() const
            {
                return std::to_string(bestCost().value());
            }


            ////////////////////////////////////////////////////////////////////////////////////////////////////
            // bispace 
            bool useBispace_{true};

            bool treatedAsMultiSubapce_{false};


            ////////////////////////////////////////////////////////////////////////////////////////////////////
            // safety certificate 
            /** \brief Compute distance between motions (actually distance between contained states) */
            double distanceFunction(const Motion *a, const Motion *b) const
            {
                return si_->distance(a->state, b->state);
            }

            double distanceFunction(const base::SafetyCertificate *a, const base::SafetyCertificate *b) const
            {
                std::vector<double> cdist = distanceCertificate_(a->state, b->state);
                double dist = std::accumulate(cdist.begin(), cdist.end(), 0.0);
                return dist;
            }

            base::PlannerStatus prepareSolve(const base::PlannerTerminationCondition &ptc);

            /** feasible */
            void growTree(TreeGrowingInfo &tgi, Motion *& startMotion, Motion *& goalMotion);

//            GrowState biasGrow(TreeData &tree, TreeGrowingInfo &tgi, Motion *&rmotion);

            // check 
            bool isValid(base::SafetyCertificate *sc);

            bool isValid(const base::State *state);

            // check
            bool isValid(Motion *motion, bool start, bool add = false);

            void removeInvalidCertificate(Motion *motion, bool start);

            bool checkInterMotion1(Motion *smotion, Motion *gmotion, bool start);

            bool checkInterMotion2(Motion *smotion, Motion *gmotion, bool start);

            void addIntermediateMotion(Motion *smotion, Motion *gmotion, bool start, Motion *last);

            bool backPathRewireMotion(Motion *motion, bool start, Motion *&pmotion);

            /** \brief Free the memory allocated by this planner */
            void freeMemory();

            void removeFromSafetyCerficate(Motion *motion);

            void freeCertificate(base::SafetyCertificate *sc)
            {
                if (sc->state)
                    si_->freeState(sc->state);
                if (sc->contact)
                    delete sc->contact;
                delete sc;
            }

            std::vector<SafetyCertificateWithElems *> snne_;

            std::shared_ptr<NearestNeighbors<base::SafetyCertificate *>> onn_;

            SafetyCertificateChecker safetyCertificateChecker_;

            CollisionCertificateChecker collisionCertificateChecker_;

            DistanceCertificate distanceCertificate_;

            unsigned int oscNum_{20};

            unsigned int certificateDim_{1u};

            /** \brief Whether or not calculate the exact certificate distance according the graph data structure. */
            bool strictCertificate_{false};

            bool useCollisionCertificateChecker_{false};
        };
    }
}

#endif
