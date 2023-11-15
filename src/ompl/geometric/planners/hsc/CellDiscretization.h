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

#ifndef OMPL_GEOMETRIC_PLANNERS_HSC_CELLDISCRETIZATION_
#define OMPL_GEOMETRIC_PLANNERS_HSC_CELLDISCRETIZATION_

#include "ompl/base/Planner.h"
#include "ompl/datastructures/GridNR.h"
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
        class CellDiscretization
        {
        public:
            /** \brief The data held by a cell in the grid of motions */
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

            /** \brief The datatype for the maintained grid datastructure */
            using Grid = GridNR<CellData *>;

            /** \brief The datatype for the maintained grid cells */
            using Cell = typename Grid::Cell;

            /** \brief The datatype for the maintained grid coordinates */
            using Coord = typename Grid::Coord;

            /** \brief The signature of a function that frees the memory for a motion */
            using FreeMotionFn = typename std::function<void(Motion *)>;

            using IsEnabledMotionFn = std::function<bool(Motion *)>;

            CellDiscretization()
              : grid_(0), size_(0)
            {
            }

            ~CellDiscretization()
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

            /** \brief Set the dimension of the grid to be maintained */
            void setDimension(unsigned int dim)
            {
                grid_.setDimension(dim);
            }

            void setNeighborCell(int cell)
            {
                grid_.setNeighborCell(cell);
            }

            void setMaxNeighborCell(int cell)
            {
                grid_.setMaxNeighborCell(cell);
            }

            /** \brief Restore the CellDiscretization to its original form */
            void clear()
            {
                freeMemory();
                size_ = 0;
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
                    freeCellData(it->second->data);
                grid_.clear();
            }

            /** \brief Add a motion to the grid containing motions. As
                a hint, \e dist specifies the distance to the goal
                from the state of the motion being added. The function
                returns the number of cells created to accommodate the
                new motion (0 or 1). The CellDiscretization takes
                ownership of the motion passed as argument, and the
                memory for the motion is freed by calling the function
                passed to the constructor. */
            void addMotion(Motion *motion, const Coord &coord)
            {
                Cell *cell = grid_.getCell(coord);
                if (!cell)
                {
                    cell = grid_.createCell(coord);
                    cell->data = new CellData();
                    grid_.add(cell);
                }
                cell->data->motions.push_back(motion);
                ++size_;
                motion->cell = cell;
            }

            /** \brief Remove a motion from the cell.
             * The user code must update the disabled motion number, and the cell's score
             */
            bool removeMotion(Motion *motion)
            {
                Cell *cell = motion->cell;
                if (!cell)
                    return false;
                bool found = false;
                for (std::size_t i = cell->data->motions.size() - 1; i < cell->data->motions.size(); i--)
                    if (cell->data->motions[i] == motion)
                    {
                        std::iter_swap(cell->data->motions.begin() + i, cell->data->motions.end() - 1);
                        cell->data->motions.pop_back();
                        found = true;
                        --size_;
                        motion->cell = nullptr;
                        break;
                    }
                return found;
            }

            void enableSort(Cell *cell)
            {
                EnableSort esort(isEnabledMotion_);
                std::sort(cell->data->motions.begin(), cell->data->motions.end(), esort);
            }

            const Grid &getGrid() const
            {
                return grid_;
            }

            Grid &getGrid()
            {
                return grid_;
            }

            void getPlannerData(base::PlannerData &data, int tag, bool start) const
            {
                std::vector<CellData *> cdata;
                grid_.getContent(cdata);
                for (std::size_t i = 0; i < cdata.size(); ++i)
                    for (std::size_t j = 0; j < cdata[i]->motions.size(); ++j)
                    {
                        if (cdata[i]->motions[j]->parent == nullptr)
                        {
                            if (start)
                                data.addStartVertex(base::PlannerDataVertex(cdata[i]->motions[j]->state, tag));
                            else
                                data.addGoalVertex(base::PlannerDataVertex(cdata[i]->motions[j]->state, tag));
                        }
                        else
                        {
                            if (start)
                                data.addEdge(base::PlannerDataVertex(cdata[i]->motions[j]->parent->state, tag),
                                             base::PlannerDataVertex(cdata[i]->motions[j]->state, tag));
                            else
                                data.addEdge(base::PlannerDataVertex(cdata[i]->motions[j]->state, tag),
                                             base::PlannerDataVertex(cdata[i]->motions[j]->parent->state, tag));
                        }
                    }
            }

        private:
            /** \brief Free the memory for the data contained in a grid cell */
            void freeCellData(CellData *cdata)
            {
                for (std::size_t i = 0; i < cdata->motions.size(); ++i)
                    freeMotion_(cdata->motions[i]);
                delete cdata;
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

            /** \brief Method that can free the memory for a stored motion */
            FreeMotionFn freeMotion_;

            IsEnabledMotionFn isEnabledMotion_;
        };
    }
}

#endif
