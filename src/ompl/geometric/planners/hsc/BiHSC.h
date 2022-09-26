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

#ifndef OMPL_GEOMETRIC_PLANNERS_HSC_BIHSC_
#define OMPL_GEOMETRIC_PLANNERS_HSC_BIHSC_

#include "ompl/geometric/planners/PlannerIncludes.h"
#include "ompl/geometric/planners/hsc/SimpleGrid.h"

#include "ompl/datastructures/NearestNeighbors.h"
#include "ompl/datastructures/BinaryHeap.h"
#include "ompl/datastructures/GridN.h"
#include "ompl/datastructures/Grid.h"
#include "ompl/datastructures/PDF.h"

#include "ompl/base/State.h"
#include "ompl/base/SafetyCertificate.h"
#include "ompl/base/OptimizationObjective.h"
#include "ompl/base/ProjectionEvaluator.h"

#include <unordered_set>
#include <queue>
#include <deque>
#include <utility>

namespace ompl
{
    namespace geometric
    {
        /** \brief Bidirectional Hybrid Safety Certificate */
        class BiHSC : public base::Planner
        {
        public:

            /** \brief Check if the state is certificated as collision-free by sc, and return the certificate distances */
            using SafetyCertificateChecker = std::function<bool(const base::State *state, const base::SafetyCertificate *sc, std::vector<double> &dist)>;

            using CollisionCertificateChecker = std::function<bool(const base::State *state, const std::vector<base::SafetyCertificate *> &ocv)>;

            /** \brief Calculate the certificate distances between two states */
            using DistanceCertificate = std::function<std::vector<double>(const base::State *a, const base::State *b)>;

            BiHSC(const base::SpaceInformationPtr &si);

            ~BiHSC() override;

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

            void setMaxInvalidNodeRatio(double ratio)
            {
                maxInvalidNodeRatio_ = ratio;
            }

            double getMaxInvalidNodeRatio() const
            {
                return maxInvalidNodeRatio_;
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

            void setRewire(bool rewire)
            {
                rewire_ = rewire;
            }

            bool getRewire() const
            {
                return rewire_;
            }

            void setRewireSort(bool rewire)
            {
                rewireSort_ = rewire;
            }

            bool getRewireSort() const
            {
                return rewireSort_;
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

            enum PathValid
            {
                UnCkeckedP,
                LazyValid,
                ValidP 
            };

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

            using CellDiscretizationData = GridN<CellData *>;
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
                    if (motion1->valid == ValidP && motion2->valid == ValidP)
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
                        return motion2->valid == ValidP;
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

            struct SafetyCertificateWithElems 
            {
                base::SafetyCertificate *sc;
                std::vector<Motion *> objects;
            };

            struct SafetyCertificatePairs
            {
                SafetyCertificateWithElems *sce;
                std::vector<double> scd;
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

                StateValid stateValid{UnCkecked};

                PathValid valid{UnCkeckedP};

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

                SafetyCertificateWithElems *sce{nullptr};

                /** \brief The certificate distance between the motion and its certificate node */
                std::vector<double> scd;
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

            void processSolution(const Motion *bestStartMotion, const Motion *bestGoalMotion);

            bool batchGrow(bool &startTree);

            /** feasible */
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

            bool isPathValid(Motion *motion, Motion *otherMotion, bool start);

            /** \brief Check if the connected path is valid */
            bool isPathValid(Motion *motion, Motion *otherMotion);

            /** \brief Check from the root to the connect point, stop immediately if an invalid path is found  */
            bool isPathValid(Motion *motion, bool start);

            /** \brief Check from the connect point to the root */
            bool isPathValidInter(Motion *motion, bool start);

            bool isPathValidLazy(Motion *motion, bool start);

            void addToDisc(CellDiscretizationData &disc, Motion *motion, const Coord& coord);

            void removeFromDisc(CellDiscretizationData &disc, Motion *motion);

            void removeInvalidMotionsDirectly();

            void removeInvalidMotionsDirectlyTree();

            void addToTree(TreeData &tree, Motion *motion);

            struct MotionPDF;
            void removeFromTree(MotionPDF &pdf, Motion *motion);

            void removeInvalidMotions();

            void removeInvalidMotionsTree(double ratio);

            void enableMotionInDisc(Motion *motion);

            void connectToPmotion(Motion *motion, Motion *pmotion, bool start) const;

            /** \brief Updates the cost of the children of this node if the cost up to this node has changed */
            void updateChildCosts(Motion *motion) const;

            bool checkMotion(Motion *pmotion, Motion *motion, bool start);

            bool checkMotionLazy(Motion *pmotion, Motion *motion, bool start);

            bool checkInterMotion(Motion *pmotion, Motion *motion, bool start);

            bool checkInterMotionLazy(Motion *pmotion, Motion *motion, bool start);

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

            void setMotionInfinityCostWithDisable(Motion *motion, std::unordered_set<Cell *> &cells) const;

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

            /** \brief Stores the start states as Motions. */
            std::vector<Motion *> startMotions_;

            /** \brief A list of states in the tree that satisfy the goal condition */
            std::vector<Motion *> goalMotions_;

            std::vector<Motion *> invalidStartMotions_;

            std::vector<Motion *> invalidGoalMotions_;

            std::size_t invalidStartNum_{0}, invalidGoalNum_{0};

            /** \brief The pair of motions in each tree connected during planning. */
            std::vector<std::pair<Motion *, Motion *>> connectionPoint_;

            double maxDistance_{0.0};

            double penDistance_{0.0};

            std::vector<double> maxDistanceV_, penDistanceV_;

            double maxInvalidNodeRatio_{0.1};

            bool lazyPath_{true};

            bool lazyNode_{true};

            bool addIntermediateState_{true};

            bool rewire_{true};

            bool rewireSort_{true};

            bool symmetric_{true};

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
            // safety certificate 
            double distanceFunction(const base::SafetyCertificate *a, const base::SafetyCertificate *b) const
            {
                std::vector<double> cdist = distanceCertificate_(a->state, b->state);
                double dist = std::accumulate(cdist.begin(), cdist.end(), 0.0);
                return dist;
            }

            base::PlannerStatus prepareSolve(const base::PlannerTerminationCondition &ptc);

            /** feasible */
            bool growTree(TreeGrowingInfo &tgi, Motion *rmotion, bool &otherSide, bool &change);

            bool growTreeSingleSpace(TreeGrowingInfo &tgi, Motion *rmotion, bool &otherSide, bool &change);

            bool growTreeMultiSpace(TreeGrowingInfo &tgi, Motion *rmotion, bool &otherSide, bool &change);

//            GrowState biasGrow(TreeData &tree, TreeGrowingInfo &tgi, Motion *&rmotion);

            // check 
            bool isValid(base::SafetyCertificate *sc);

            bool isValid(const base::State *state);

            // check
            bool isValid(Motion *motion, bool start);

            std::vector<Motion *> removeInvalidCertificate(SafetyCertificateWithElems *msce, const std::vector<double> &scd, bool start);

            std::vector<Motion *> removeInvalidCertificateInter(SafetyCertificateWithElems *msce, const std::vector<double> &scd, bool start);

            void removeInvalidCertificate(SafetyCertificateWithElems *msce, const std::vector<double> &scd, bool start,
                    std::queue<SafetyCertificatePairs, std::deque<SafetyCertificatePairs>> &scQueue);

            void removeInvalidCertificate(std::queue<SafetyCertificatePairs, std::deque<SafetyCertificatePairs>> &scQueue, bool start);

            bool certificateOutside(const std::vector<double> &scd, const std::vector<double> &confidence) const;

            /** \brief The metric space is symmetric */
            bool checkInterMotion1(Motion *smotion, Motion *gmotion, bool start);

            /** \brief The metric space is not symmetric */
            bool checkInterMotion2(Motion *smotion, Motion *gmotion, bool start);

            void addIntermediateMotion(Motion *pmotion, Motion *motion, bool start, Motion *last);

            void addIntermediateMotionLazy(Motion *pmotion, bool start, Motion *last);

            /** \brief Check from the connect point to the root */
            bool isStateValid(Motion *motion, bool start);

            bool backPathRewireMotion(Motion *motion, bool start, Motion *&pmotion);

            bool backPathRewireMotionLazy(Motion *motion, bool start, Motion *&pmotion);

            void removeFromSafetyCerficate(Motion *motion);

            /** \brief Free the memory allocated by this planner */
            void freeMemory();

            void freeCertificate(base::SafetyCertificate *sc)
            {
                if (sc->state)
                    si_->freeState(sc->state);
                if (sc->contact)
                    delete sc->contact;
                delete sc;
            }

            std::vector<SafetyCertificateWithElems *> ssnne_, gsnne_;

            using SimpleGridSC = SimpleGrid<ompl::base::SafetyCertificate>;
            using CellSC = SimpleGridSC::Cell;
            using CoordSC = SimpleGridSC::Coord;

            std::shared_ptr<SimpleGridSC> onn_;

            SafetyCertificateChecker safetyCertificateChecker_;

            CollisionCertificateChecker collisionCertificateChecker_;

            DistanceCertificate distanceCertificate_;

            unsigned int certificateDim_{1u};

            bool useCollisionCertificateChecker_{false};
        };
    }
}

#endif
