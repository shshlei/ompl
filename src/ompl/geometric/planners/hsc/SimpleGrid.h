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

#ifndef OMPL_GEOMETRIC_PLANNERS_HSC_SIMPLEGRID_
#define OMPL_GEOMETRIC_PLANNERS_HSC_SIMPLEGRID_

#include "ompl/datastructures/Grid.h"

namespace ompl
{
    namespace geometric
    {
        template <typename T>
        class SimpleGrid
        {
        public:
            struct CellData
            {
                CellData() = default;
                ~CellData() = default;
                std::vector<T *> datas;
            };

            using SGrid = Grid<CellData *>;
            using Cell = typename SGrid::Cell;
            using Coord = typename SGrid::Coord;

            SimpleGrid() : grid_(0), size_(0)
            {
            }

            ~SimpleGrid()
            {
                freeMemory();
            }

            void setDimension(unsigned int dim)
            {
                grid_.setDimension(dim);
            }

            unsigned int add(T *data, const Coord &coord)
            {
                Cell *cell = grid_.getCell(coord);
                unsigned int created = 0;
                if (cell)
                    cell->data->datas.push_back(data);
                else
                {
                    cell = grid_.createCell(coord);
                    cell->data = new CellData();
                    cell->data->datas.push_back(data);
                    grid_.add(cell);
                }
                size_++;
                return created;
            }

            Cell *getCell(const Coord &coord)
            {
                return grid_.getCell(coord);
            }

            const Cell *getCell(const Coord &coord) const
            {
                return grid_.getCell(coord);
            }

            const SGrid &getGrid() const
            {
                return grid_;
            }

            SGrid &getGrid()
            {
                return grid_;
            }

            void freeMemory()
            {
                for (auto it = grid_.begin(); it != grid_.end(); ++it)
                    freeCellData(it->second->data);
                grid_.clear();
            }

            std::size_t size() const
            {
                return size_;
            }

            void clear()
            {
                grid_.clear();
            }

        private:
            void freeCellData(CellData *cdata)
            {
                delete cdata;
            }

            SGrid grid_;
            std::size_t size_{0};
        };
    }
}
#endif
