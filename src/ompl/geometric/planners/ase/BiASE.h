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

#ifndef OMPL_GEOMETRIC_PLANNERS_ASE_BIASE_
#define OMPL_GEOMETRIC_PLANNERS_ASE_BIASE_

#include "ompl/geometric/planners/PlannerIncludes.h"

#include "ompl/datastructures/NearestNeighbors.h"
#include "ompl/datastructures/BinaryHeap.h"
#include "ompl/datastructures/GridNR.h"
#include "ompl/datastructures/Grid.h"
#include "ompl/datastructures/PDF.h"

#include "ompl/base/State.h"
#include "ompl/base/OptimizationObjective.h"
#include "ompl/base/ProjectionEvaluator.h"

#include "ompl/base/samplers/AdInformedStateSampler.h"

#include <unordered_set>

namespace ompl
{
    namespace geometric
    {
        /** \brief Bidirectional Adaptive Space Expansion */
        class BiASE : public base::Planner
        {
        public:

            BiASE(const base::SpaceInformationPtr &si);

            ~BiASE() override;

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
                if (!bispace)
                    treatedAsMultiSubapce_ = false;
                useBispace_ = bispace;
            }

            bool getUseBispace() const
            {
                return useBispace_;
            }

            void setUseBiasGrow(bool use)
            {
                useBiasGrow_ = use;
                if (!useBispace_ && use)
                {
                    OMPL_WARN("%s: The bias grow feature is not usable since the bispace feature is not set.", getName().c_str());
                    useBiasGrow_ = false;
                }
            }

            bool getUseBiasGrow() const
            {
                return useBiasGrow_;
            }

            void setTreatedAsMultiSubapce(bool sub)
            {
                treatedAsMultiSubapce_ = sub;
                if (!useBispace_ && sub)
                {
                    OMPL_WARN("%s: The treate as multisubspace feature is not usable since the bispace feature is not set.", getName().c_str());
                    treatedAsMultiSubapce_ = false;
                }
            }

            bool getTreatedAsMultiSubapce() const
            {
                return treatedAsMultiSubapce_;
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

            void setUpdateNbCell(bool update)
            {
                updateNbCell_ = update;
            }

            bool getUpdateNbCell() const
            {
                return updateNbCell_;
            }

            void setInformedSampling(bool informedSampling)
            {
                useInformedSampling_ = informedSampling;
            }

            bool getInformedSampling() const
            {
                return useInformedSampling_;
            }

            void setSampleRejection(bool reject)
            {
                useRejectionSampling_ = reject;
            }

            bool getSampleRejection() const
            {
                return useRejectionSampling_;
            }

            void setNumSamplingAttempts(unsigned int numAttempts)
            {
                numSampleAttempts_ = numAttempts;
            }

            unsigned int getNumSamplingAttempts() const
            {
                return numSampleAttempts_;
            }

        protected:

            ////////////////////////////////////////////////////////////////////////////////////////////////////
            // general functions
            class Motion;

            struct CellData
            {
                CellData() = default;

                ~CellData() = default;

                /** \brief The set of motions contained in this grid cell */
                std::vector<Motion *> motions;

                Motion *cmotion{nullptr};

                Motion *mmotion{nullptr};

                /** \brief The disabled motion number in this cell */
                std::size_t disabled{0};
            };

            using CellDiscretizationData = GridNR<CellData *>;
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

            struct EnableSort
            {
                EnableSort(base::OptimizationObjectivePtr opt) : opt_(std::move(opt))
                {
                }
                inline bool operator()(Motion *motion1, Motion *motion2)
                {
                    if (opt_->isFinite(motion1->cost) && opt_->isFinite(motion2->cost))
                        return false;
                    else
                        return opt_->isFinite(motion1->cost);
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

                /** \brief The parent motion in the exploration tree */
                Motion *pmotion{nullptr};

                Cell *cell{nullptr};

                bool valid{false};

                bool inConnection{false};

                bool middle{false};

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

            /** \brief Compute distance between motions (actually distance between contained states) */
            double distanceFunction(const Motion *a, const Motion *b) const
            {
                return si_->distance(a->state, b->state);
            }

            /** \brief Information attached to growing a tree of motions (used internally) */
            struct TreeGrowingInfo
            {
                base::State *xstate;
                Motion *xmotion;
                bool start;
            };

            /** \brief A nearest-neighbor datastructure representing a tree of motions */
            using TreeData = std::shared_ptr<NearestNeighbors<Motion *>>;

            base::PlannerStatus prepareSolve(const base::PlannerTerminationCondition &ptc);

            void processSolution(const Motion *bestStartMotion, const Motion *bestGoalMotion);

            /** feasible */
            bool growTree(TreeGrowingInfo &tgi, Motion *rmotion, bool &otherSide, bool &change, bool &add);

            bool growTreeSingleSpace(TreeGrowingInfo &tgi, Motion *rmotion, bool &otherSide, bool &change, bool &add);

            bool growTreeMultiSpace(TreeGrowingInfo &tgi, Motion *rmotion, bool &otherSide, bool &change, bool &add);

            Motion *selectNMotion(const TreeData &tree, Motion *rmotion, bool &null);

            Motion *selectMotionInCell(Cell *cell);

            bool growCurrentTree(const base::State *state) const;

            bool growStartTree(const base::State *state) const;

            double penetrationDistance(const base::State *nstate, const base::State *state, bool start) const;

            bool growCurrentTree(const base::State *state, std::size_t sub) const;

            bool growStartTree(const base::State *state, std::size_t sub) const;

            double penetrationDistance(const base::State *nstate, const base::State *state, bool start, std::size_t sub) const;

            void rewireTree(BinaryHeap<Motion *, MotionCompare> &bh, Motion *m, bool start);

            void updateQueue(BinaryHeap<Motion *, MotionCompare> &bh, Motion *m);

//            bool biasGrow(TreeData &tree, TreeGrowingInfo &tgi, Motion *&rmotion, bool optimal = false);

            bool backPathRewireMotion(Motion *motion, bool start, Motion *&pmotion);

            void addToDisc(CellDiscretizationData &disc, Motion *motion, const Coord& coord);

            void removeFromDisc(CellDiscretizationData &disc, Motion *motion);

            struct MotionPDF;

            void enableMotionInDisc(Motion *motion);

            void connectToPmotion(Motion *motion, Motion *pmotion, bool start) const;

            /** \brief Updates the cost of the children of this node if the cost up to this node has changed */
            void updateChildCosts(Motion *motion) const;

            // check 
            bool isValid(const base::State *state);

            // check
            bool checkMotion(Motion *pmotion, Motion *motion, bool start);

            bool checkInterMotion(Motion *pmotion, Motion *motion, bool start);

            bool checkInterMotion1(Motion *smotion, Motion *gmotion, bool start);

            bool checkInterMotion2(Motion *smotion, Motion *gmotion, bool start);

            void addIntermediateMotion(Motion *pmotion, Motion *motion, bool start, Motion *last);

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

            /** \brief Free the memory allocated by this planner */
            void freeMemory();

            void freeMotion(Motion *motion)
            {
                if (motion->state != nullptr)
                    si_->freeState(motion->state);
                delete motion;
            }

            base::ProjectionEvaluatorPtr projectionEvaluator_;

            /** \brief The start tree */
            TreeData tStart_;

            /** \brief The goal tree */
            TreeData tGoal_;

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

            /** \brief The pair of motions in each tree connected during planning. */
            std::vector<std::pair<Motion *, Motion *>> connectionPoint_;

            double maxDistance_{0.0};

            double penDistance_{0.0};

            std::vector<double> maxDistanceV_, penDistanceV_;

            bool addIntermediateState_{true};

            bool backRewire_{true};

            bool rewireSort_{true};

            bool symmetric_{true};

            bool updateNbCell_{false};

            /** \brief Objective we're optimizing */
            base::OptimizationObjectivePtr opt_;

            MotionCompare mc_;

            BinaryHeap<Motion *, MotionCompare> bh_;

            base::StateSamplerPtr sampler_;

            RNG rng_;


            ////////////////////////////////////////////////////////////////////////////////////////////////////
            // bispace 
            struct MotionInfo;

            using GridCell = Grid<MotionInfo>::Cell;

            using CellPDF = PDF<GridCell *>;

            struct MotionInfo
            {
                Motion *operator[](unsigned int i)
                {
                    return motions_[i];
                }
                std::vector<Motion *>::iterator begin()
                {
                    return motions_.begin();
                }
                std::vector<Motion *>::iterator end()
                {
                    return motions_.end();
                }
                std::vector<Motion *>::iterator erase(std::vector<Motion *>::iterator iter)
                {
                    std::vector<Motion *>::iterator it = motions_.erase(iter);
                    return it;
                }
                void push_back(Motion *m)
                {
                    motions_.push_back(m);
                }
                void pop_back()
                {
                    motions_.pop_back();
                }
                std::size_t size() const
                {
                    return motions_.size();
                }
                bool empty() const
                {
                    return motions_.empty();
                }

                std::vector<Motion *> motions_;
                CellPDF::Element *elem_;
            };

            struct MotionPDF
            {
                MotionPDF() = default;

                Grid<MotionInfo> grid{0};

                std::size_t size{0};

                CellPDF pdf;
            };

            MotionPDF startBiasPdf_, goalBiasPdf_;

            void addPdfMotion(MotionPDF &pdf, Motion *motion, bool start);

            Motion *selectPdfMotion(MotionPDF &pdf, GridCell *&cell);

            void removePdfMotion(MotionPDF &pdf, Motion *motion);

            double startBiasProb_{0.0}, goalBiasProb_{0.0};

            bool useBispace_{true};

            bool useBiasGrow_{false};

            bool treatedAsMultiSubapce_{false};

            std::vector<ompl::base::State *> startBiasStates_, goalBiasStates_;


            ////////////////////////////////////////////////////////////////////////////////////////////////////
            // adaptive informed sampling 
            bool batchGrow(bool &startTree);

            bool selectCMotion(std::size_t &index, bool &reverse);

            /** \brief Check if the connected path is valid */
            bool isPathValid(Motion *motion, Motion *otherMotion);

            /** \brief Check from the root to the connect point, stop immediately if an invalid path is found  */
            bool isPathValid(Motion *motion, bool start);

            /** \brief Check from the connect point to the root */
            bool isPathValidInter(Motion *motion, bool start);


            ////////////////////////////////////////////////////////////////////////////////////////////////////
            // adaptive informed
            using NumPdf = PDF<std::size_t>;
            using NumElem = NumPdf::Element;
            
            void processAdEllipsoidRind(bool clearoradd, unsigned int &adinfcount);

            /** \brief Create the samplers */
            base::AdInformedSamplerPtr allocInfSampler(const base::State *s1, const base::State *s2,
                                                       const base::Cost &minCost, const base::Cost &maxCost);

            void clearStartAdInfSampler() 
            {
                if (!startAdInfSamplers_.empty())
                {
                    for (auto & sampler : startAdInfSamplers_)
                        sampler.reset();
                    startAdInfSamplers_.clear();
                    startAdInfPdf_.clear();
                    startAdElems_.clear();
                }
            }

            void clearGoalAdInfSampler() 
            {
                if (!goalAdInfSamplers_.empty())
                {
                    for (auto & sampler : goalAdInfSamplers_)
                        sampler.reset();
                    goalAdInfSamplers_.clear();
                    goalAdInfPdf_.clear();
                    goalAdElems_.clear();
                }
            }

            void localInfeasible(int &tree, bool &locals, bool &localg);

            void calculateInfSampler(bool local, bool start, bool &updated);

            void startInfSampler(bool local, std::vector<base::AdInformedSamplerPtr> &infSamplers);

            void startLocalInfSampler(std::vector<base::AdInformedSamplerPtr> &infSamplers);

            void goalInfSampler(bool local, std::vector<base::AdInformedSamplerPtr> &infSamplers);

            void goalLocalInfSampler(std::vector<base::AdInformedSamplerPtr> &infSamplers);

            void calculateInfProb(bool update);

            void calculateInfProb(const std::vector<base::AdInformedSamplerPtr> &startInfSamplers, 
                                  const std::vector<base::AdInformedSamplerPtr> &goalInfSamplers,
                                  NumPdf &startInfPdf, NumPdf &goalInfPdf, 
                                  std::vector<NumElem *> &startelems, std::vector<NumElem *> &goalelems, bool update);

            bool sampleUniformAd(base::State *state, bool start);

            virtual bool sampleUniform(base::State *state, const std::vector<base::AdInformedSamplerPtr> &infSamplers, const NumPdf &adInfPdf);

            std::vector<base::AdInformedSamplerPtr> startAdInfSamplers_;

            std::vector<base::AdInformedSamplerPtr> goalAdInfSamplers_;

            NumPdf startAdInfPdf_;
            std::vector<NumElem *> startAdElems_;

            NumPdf goalAdInfPdf_;
            std::vector<NumElem *> goalAdElems_;

            double startAdInfProb_{0.0};

            double localRatio_{0.75};

            double factor_{1.1};

            int tree_{-1};

            /** \brief Option to use informed sampling */
            bool useInformedSampling_{false};

            /** \brief The status of the sample rejection parameter. */
            bool useRejectionSampling_{false};

            /** \brief The number of attempts to make at informed sampling */
            unsigned int numSampleAttempts_{100u};

            bool ais_{false};
        };
    }
}

#endif
