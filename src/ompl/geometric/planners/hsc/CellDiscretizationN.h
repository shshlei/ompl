/*********************************************************************
* Software License Agreement (BSD License)
*
*  Copyright (c) 2008, Rice University
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

/* Author: Shi Shenglei */

#ifndef OMPL_GEOMETRIC_PLANNERS_HSC_CELLDISCRETIZATIONN_
#define OMPL_GEOMETRIC_PLANNERS_HSC_CELLDISCRETIZATIONN_

#include "ompl/base/Planner.h"
#include "ompl/datastructures/GridNearestB.h"
#include "ompl/util/Exception.h"
#include <functional>
#include <utility>
#include <vector>
#include <limits>
#include <cassert>
#include <cstdlib>
#include <cmath>

namespace ompl
{
    namespace geometric
    {
        /** \brief One-level CellDiscretization */
        template <typename Motion>
        class CellDiscretizationN
        {
        public:
            /** \brief The data held by a cell in the grid of motions */
            struct CellData
            {
                CellData() = default;

                ~CellData() = default;

                Motion *cmotion{nullptr};

                Motion *mmotion{nullptr};

                /** \brief The number of times this cell has been
                    selected for expansion */
                std::size_t selections{1};

                /** \brief The iteration at which this cell was created */
                std::size_t iteration{0};

                /** \brief The disabled motion number in this cell */
                std::size_t disabled{0};

                /** \brief The root motion number in this cell */
                std::size_t root{0};

                std::size_t lpd{1};

                /** \brief A measure of coverage for this cell. For
                    this implementation, this is the sum of motion
                    lengths */
                double coverage{0.0};

                /** \brief A heuristic score computed based on
                    distance to goal (if available), successes and
                    failures at expanding from this cell. */
                double score{1.0};

                /** \brief The computed importance (based on other class members) */
                double importance{0.0};
            };

            /** \brief Definintion of an operator passed to the Grid
                structure, to order cells by importance */
            struct OrderCellsByImportance
            {
                /** \brief Order function */
                bool operator()(const CellData *const a, const CellData *const b) const
                {
                    return a->importance > b->importance;
                }
            };

            /** \brief The datatype for the maintained grid datastructure */
            using Grid = GridNearestB<Motion *, CellData *, OrderCellsByImportance>;

            /** \brief The datatype for the maintained grid cells */
            using Cell = typename Grid::Cell;

            /** \brief The datatype for the maintained grid coordinates */
            using Coord = typename Grid::Coord;

            /** \brief The signature of a function that frees the memory for a motion */
            using FreeMotionFn = typename std::function<void(Motion *)>;

            using IsEnabledMotionFn = std::function<bool(Motion *)>;

            CellDiscretizationN()
              : grid_(0), size_(0), iteration_(1), recentCell_(nullptr)
            {
                grid_.onCellUpdate(computeImportance, nullptr);
                selectBorderFraction_ = 0.9;
            }

            ~CellDiscretizationN()
            {
                freeMemory();
            }

            void setFreeMotionFn(FreeMotionFn freeMotion)
            {
                freeMotion_ = std::move(freeMotion);
            }

            void setIsEnabledMotionFn(IsEnabledMotionFn func)
            {
                isEnabledMotion_ = std::move(func);
            }

            /** \brief Set the fraction of time for focusing on the
                border (between 0 and 1). This is the minimum fraction
                used to select cells that are exterior (minimum
                because if 95% of cells are on the border, they will
                be selected with 95% chance, even if this fraction is
                set to 90%)*/
            void setBorderFraction(double bp)
            {
                if (bp < std::numeric_limits<double>::epsilon() || bp > 1.0)
                    throw Exception("The fraction of time spent selecting border cells must be in the range (0,1]");
                selectBorderFraction_ = bp;
            }

            /** \brief Set the fraction of time for focusing on the
                border (between 0 and 1). */
            double getBorderFraction() const
            {
                return selectBorderFraction_;
            }

            /** \brief Set the dimension of the grid to be maintained */
            void setDimension(unsigned int dim)
            {
                grid_.setDimension(dim);
            }

            void setNeighborCell(int cell)
            {
                grid_.setNeighborCell(cell);
            }

            /** \brief Restore the CellDiscretizationN to its original form */
            void clear()
            {
                freeMemory();
                size_ = 0;
                iteration_ = 1;
                recentCell_ = nullptr;
                rootCells_.clear();
            }

            void countIteration()
            {
                ++iteration_;
            }

            std::size_t getMotionCount() const
            {
                return size_;
            }

            std::size_t getCellCount() const
            {
                return grid_.size();
            }

            /** \brief Free the memory for the motions contained in a grid */
            void freeMemory()
            {
                for (auto it = grid_.begin(); it != grid_.end(); ++it)
                    freeCellData(static_cast<Cell *>(it->second));
                grid_.clear();
            }

            void setRootCell(Cell *cell)
            {
                rootCells_.clear();
                rootCells_.push_back(cell);
            }

            void addRootCell(Cell *cell)
            {
                rootCells_.push_back(cell);
            }

            /** \brief Add a motion to the grid containing motions. As
                a hint, \e dist specifies the distance to the goal
                from the state of the motion being added. The function
                returns the number of cells created to accommodate the
                new motion (0 or 1). The CellDiscretizationN takes
                ownership of the motion passed as argument, and the
                memory for the motion is freed by calling the function
                passed to the constructor. */
            unsigned int addMotion(Motion *motion, const Coord &coord, double dist = 0.0)
            {
                Cell *cell = grid_.getCell(coord);

                unsigned int created = 0;
                if (cell)
                {
                    grid_.add(motion, cell);
                    cell->auxData->coverage += 1.0;
                    if (cell->removed)
                    {
                        cell->auxData->iteration = iteration_;
                        cell->auxData->score = (1.0 + log((double)(iteration_))) / (1.0 + dist);
                        recentCell_ = cell;
                        grid_.lazyAdd(cell);
                    }
                }
                else
                {
                    cell = grid_.createCell(coord);
                    grid_.add(motion, cell);
                    cell->auxData = new CellData();
                    cell->auxData->coverage = 1.0;
                    cell->auxData->iteration = iteration_;
                    cell->auxData->selections = 1;
                    cell->auxData->score = (1.0 + log((double)(iteration_))) / (1.0 + dist);
                    grid_.add(cell);
                    recentCell_ = cell;
                    created = 1;
                    for (auto & c : rootCells_)
                    {
                        std::size_t d = (std::size_t)((cell->coord - c->coord).template lpNorm<Eigen::Infinity>() + 1);
                        if (cell->auxData->lpd < d)
                            cell->auxData->lpd = d;
                    }
                }
                ++size_;
                motion->cell = cell;
                return created;
            }

            /** \brief Select a motion and the cell it is part of from
                the grid of motions. This is where preference is given
                to cells on the boundary of the grid.*/
            void selectMotion(Motion *&smotion, Cell *&scell)
            {
                scell = rng_.uniform01() < std::max(selectBorderFraction_, grid_.fracExternal()) ? grid_.topExternal() :
                                                                                                   grid_.topInternal();

                ++scell->auxData->selections;
                smotion = scell->data[rng_.halfNormalInt(0, scell->data.size() - 1)];
                if (!isEnabledMotion_(smotion) && scell->data.size() > scell->auxData->disabled)
                {
                    std::size_t index = rng_.halfNormalInt(0, scell->data.size() - 1 - scell->auxData->disabled);
                    smotion = scell->data[index];
                    if (!isEnabledMotion_(smotion))
                    {
                        enableSort(scell);
                        smotion = scell->data[index];
                    }
                }
            }

            /** \brief Remove a motion from the cell.
             * The user code must update the disabled motion number, and the cell's score
             */
            bool removeMotion(Motion *motion)
            {
                Cell *cell = motion->cell;
                if (!cell)
                    return false;
                bool found = grid_.remove(motion, cell);
                if (found)
                {
                    --size_;
                    motion->cell = nullptr;
                }
                if (cell->data.empty())
                {
                    cell->auxData->disabled = 0;
                    cell->auxData->root = 0;
                    cell->auxData->cmotion = nullptr;
                    cell->auxData->mmotion = nullptr;
                    grid_.lazyRemove(cell);
                }
                return found;
            }

            void enableSort(Cell *cell)
            {
                EnableSort esort(isEnabledMotion_);
                std::sort(cell->data.begin(), cell->data.end(), esort);
            }

            void updateCell(Cell *cell)
            {
                grid_.update(cell);
            }

            void updateAll()
            {
                grid_.updateAll();
            }

            Grid &getGrid()
            {
                return grid_;
            }

            const Grid &getGrid() const
            {
                return grid_;
            }

            void getPlannerData(base::PlannerData &data, int tag, bool start, const Motion *lastGoalMotion) const
            {
                std::vector<Cell *> cells;
                grid_.getCells(cells);

                if (lastGoalMotion)
                    data.addGoalVertex(base::PlannerDataVertex(lastGoalMotion->state, tag));

                for (auto & cell : cells)
                    for (auto & motion : cell->data)
                    {
                        if (motion->parent == nullptr)
                        {
                            if (start)
                                data.addStartVertex(base::PlannerDataVertex(motion->state, tag));
                            else
                                data.addGoalVertex(base::PlannerDataVertex(motion->state, tag));
                        }
                        else
                        {
                            if (start)
                                data.addEdge(base::PlannerDataVertex(motion->parent->state, tag),
                                             base::PlannerDataVertex(motion->state, tag));
                            else
                                data.addEdge(base::PlannerDataVertex(motion->state, tag),
                                             base::PlannerDataVertex(motion->parent->state, tag));
                        }
                    }
            }

        private:
            /** \brief Free the memory for the data contained in a grid cell */
            void freeCellData(Cell *cell)
            {
                for (std::size_t i = 0; i < cell->data.size(); ++i)
                    freeMotion_(cell->data[i]);
                delete cell->auxData;
            }

            /** \brief This function is provided as a callback to the
                grid datastructure to update the importance of a
                cell */
            static void computeImportance(Cell *cell, void * /*unused*/)
            {
                CellData &cd = *(cell->auxData);
                cd.importance = cd.score * cd.lpd / ((cell->neighbors + 1) * cd.coverage * cd.selections);
                if (cd.disabled)
                    cd.importance *= (double)(cell->data.size() - cd.disabled)/(double)cell->data.size();
            }

            struct EnableSort
            {
                EnableSort(IsEnabledMotionFn isEnabledMotion) : isEnabledMotion_(isEnabledMotion)
                {
                }
                inline bool operator()(Motion *motion1, Motion *motion2)
                {
                    if (isEnabledMotion_(motion1) && isEnabledMotion_(motion2))
                        return false;
                    else
                        return isEnabledMotion_(motion1);
                }
                IsEnabledMotionFn isEnabledMotion_;
            };

            /** \brief A grid containing motions, imposed on a
                projection of the state space */
            Grid grid_;

            /** \brief The total number of motions (there can be
                multiple per cell) in the grid */
            std::size_t size_;

            /** \brief The number of iterations performed on this tree */
            unsigned int iteration_;

            /** \brief The most recently created cell */
            Cell *recentCell_;

            std::vector<Cell *> rootCells_;

            /** \brief Method that can free the memory for a stored motion */
            FreeMotionFn freeMotion_;

            IsEnabledMotionFn isEnabledMotion_;

            /** \brief The fraction of time to focus exploration on
                the border of the grid. */
            double selectBorderFraction_;

            /** \brief The random number generator */
            RNG rng_;
        };
    }
}

#endif
