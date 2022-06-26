/*********************************************************************
* Software License Agreement (BSD License)
*
*  Copyright (c) 2010, Rice University
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

#ifndef OMPL_DATASTRUCTURES_GRIDNEAREST_N_
#define OMPL_DATASTRUCTURES_GRIDNEAREST_N_

#include "ompl/datastructures/GridNearest.h"

namespace ompl
{
    /** \brief Representation of a grid where cells keep track of how many neighbors they have */
    template <typename _T, typename _TT = void>
    class GridNearestN : public GridNearest<_T, _TT>
    {
    public:
        /// Datatype for cell in base class
        using BaseCell = typename GridNearest<_T, _TT>::Cell;

        /// Datatype for array of cells in base class
        using BaseCellArray = typename GridNearest<_T, _TT>::CellArray;

        /// Datatype for cell coordinates
        using Coord = typename GridNearest<_T, _TT>::Coord;

        /// Definition of a cell in this grid
        struct Cell : public BaseCell
        {
            /// The number of neighbors
            unsigned int neighbors{0};

            /// A flag indicating whether this cell is on the border or not
            bool border{true};

            /// A flag indicating whether this cell is lazily removed 
            bool removed{false};

            std::vector<Cell *> nnbh;

            std::vector<std::vector<Cell *>> nbh;

            Cell() = default;

            ~Cell() override = default;

            EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        };

        /// The datatype for arrays of cells
        using CellArray = std::vector<Cell *>;

        /// The constructor takes the dimension of the grid as argument
        explicit GridNearestN(unsigned int dimension) : GridNearest<_T, _TT>(dimension)
        {
            hasBounds_ = false;
            overrideCellNeighborsLimit_ = false;
            setDimension(dimension);
        }

        ~GridNearestN() override = default;

        /// Update the dimension of the grid; this should not be done
        /// unless the grid is empty
        void setDimension(unsigned int dimension)
        {
            assert((GridNearest<_T, _TT>::empty()));
            GridNearest<_T, _TT>::dimension_ = dimension;
            GridNearest<_T, _TT>::maxNeighbors_ = 2 * dimension;
            if (!overrideCellNeighborsLimit_)
                interiorCellNeighborsLimit_ = GridNearest<_T, _TT>::maxNeighbors_;
        }

        void useNeighbor(bool use)
        {
            useNeighbor_ = use;
        }

        /// Store the neighbors of a cell in (2*n - 1) grid 
        void setNeighborCell(int n)
        {
            neighborCell_ = n;
        }

        /// If bounds for the grid need to be considered, we can set them here.
        /// When the number of neighbors are counted, whether the
        /// Space is bounded matters, in the sense that if a cell is on
        /// the boundary, we know some of its neighbors cannot exist.
        /// In order to allow such a cell to reflect the fact it has
        /// Achieved its maximal number of neighbors, the boundary is
        /// counted as the number of neighbors it prevents from
        /// existing.
        void setBounds(const Coord &low, const Coord &up)
        {
            lowBound_ = low;
            upBound_ = up;
            hasBounds_ = true;
        }

        /// Set the limit of neighboring cells to determine when a cell becomes interior
        /// by default, this is 2 * dimension of grid
        void setInteriorCellNeighborLimit(unsigned int count)
        {
            interiorCellNeighborsLimit_ = count;
            assert(interiorCellNeighborsLimit_ > 0);
            overrideCellNeighborsLimit_ = true;
        }

        /// Get the cell at a specified coordinate
        Cell *getCell(const Coord &coord) const
        {
            return static_cast<Cell *>(GridNearest<_T, _TT>::getCell(coord));
        }

        /** \brief Get the nearest neighbor of a point */
        _T nearest(const _T &data, const Coord &coord) const override
        {
            Cell *cell = getCell(coord);
            if (cell && !cell->removed)
                return cell->nn->nearest(data);
            else if (GridNearest<_T, _TT>::cnn_->size() > 0) 
            {
                cell = getCell(GridNearest<_T, _TT>::cnn_->nearest(coord));
                return cell->nn->nearest(data);
            }
            else 
                throw Exception("No elements found in gridn nearest neighbors data structure");
        }

        /// Get the list of neighbors for a given cell
        void neighbors(const Cell *cell, CellArray &list) const
        {
            Coord test = cell->coord;
            neighbors(test, list);
        }

        /// Get the list of neighbors for a given coordinate
        void neighbors(const Coord &coord, CellArray &list) const
        {
            Coord test = coord;
            neighbors(test, list);
        }

        /// Get the list of neighbors for a given coordinate
        void neighbors(Coord &coord, CellArray &list) const
        {
            BaseCellArray baselist;
            GridNearest<_T, _TT>::neighbors(coord, baselist);
            list.reserve(list.size() + baselist.size());
            for (const auto &c : baselist)
                list.push_back(static_cast<Cell *>(c));
        }

        /// Get the list of neighbors for a given cell
        void neighbors(const Cell *cell, int n, std::vector<CellArray> &list, bool ring = false) const
        {
            if (n == 1)
            {
                list.resize(n);
                list[0].push_back(const_cast<Cell *>(cell));
                return;
            }
            Coord test = cell->coord;
            neighbors(test, n, list, ring);
        }

        /// Get the list of neighbors for a given coordinate
        void neighbors(const Coord &coord, int n, std::vector<CellArray> &list, bool ring = false) const
        {
            Coord test = coord;
            neighbors(test, n, list, ring);
        }

        /// Get the list of neighbors for a given coordinate
        void neighbors(Coord &coord, int n, std::vector<CellArray> &list, bool ring = false) const
        {
            std::vector<BaseCellArray> baselist;
            GridNearest<_T, _TT>::neighbors(coord, n, baselist, ring);
            list.resize(baselist.size());
            if (!ring)
            {
                int index =0;
                for (const auto &cv : baselist)
                {
                    list[index].reserve(cv.size());
                    for (const auto &c : cv)
                        list[index].push_back(static_cast<Cell *>(c));
                    index++;
                }
            }
            else 
            {
                int index = n-1;
                list[index].reserve(baselist[index].size());
                for (const auto &c : baselist[index])
                    list[index].push_back(static_cast<Cell *>(c));
            }
        }

        /// Instantiate a new cell at given coordinates;
        /// Optionally return the list of future neighbors.  Note:
        /// this call only creates the cell, but does not add it to
        /// the grid.  It however updates the neighbor count for
        /// neighboring cells
        BaseCell *createCell(const Coord &coord, BaseCellArray *nbh = nullptr) override
        {
            Cell *cell = new Cell();
            cell->coord = coord;
            if (GridNearest<_T, _TT>::metric_)
            {
                if (GridNearest<_T, _TT>::useThread_)
                    cell->nn.reset(new NearestNeighborsGNAT<_T>());
                else
                    cell->nn.reset(new NearestNeighborsGNATNoThreadSafety<_T>());
            }
            else 
                cell->nn.reset(new NearestNeighborsSqrtApprox<_T>());
            cell->nn->setDistanceFunction(GridNearest<_T, _TT>::distFun_);

            if (useNeighbor_)
            {
                BaseCellArray *list = nbh ? nbh : new BaseCellArray();
                GridNearest<_T, _TT>::neighbors(cell->coord, *list);

                cell->nnbh.reserve(list->size());
                for (auto cl = list->begin(); cl != list->end(); ++cl)
                {
                    auto *c = static_cast<Cell *>(*cl);
                    c->neighbors++;
                    if (c->border && c->neighbors >= interiorCellNeighborsLimit_)
                        c->border = false;
                    cell->nnbh.push_back(c);
                    c->nnbh.push_back(cell);
                }

                cell->neighbors = numberOfBoundaryDimensions(cell->coord);
                if (removed_)
                    cell->neighbors += std::accumulate(list->begin(), list->end(), 0, [](int a, BaseCell *b){ int bb = static_cast<Cell *>(b)->removed ? 0 : 1; return a + bb;});
                else
                    cell->neighbors += list->size();
                if (cell->border && cell->neighbors >= interiorCellNeighborsLimit_)
                    cell->border = false;

                if (!nbh)
                    delete list;
            }
            else if (nbh)
                GridNearest<_T, _TT>::neighbors(cell->coord, *nbh);
            
            return cell;
        }

        /// Remove a cell from the grid. If the cell has not been
        /// Added to the grid, only update the neighbor list
        bool remove(BaseCell *cell) override
        {
            if (!cell)
                return false;
            Cell *c = static_cast<Cell *>(cell);
            if (useNeighbor_)
            {
                for (auto & cl : c->nnbh)
                {
                    cl->neighbors--;
                    if (!cl->border && cl->neighbors < interiorCellNeighborsLimit_)
                        cl->border = true;
                    for (std::size_t i = cl->nnbh.size() - 1; i < cl->nnbh.size(); i--)
                    {
                        if (cl->nnbh[i] == c)
                        {
                            std::iter_swap(cl->nnbh.begin() + i, cl->nnbh.end() - 1);
                            cl->nnbh.pop_back();
                            break;
                        }
                    }
                }
            }

            for (std::size_t i = 1; i < c->nbh.size(); i++)
            {
                for (std::size_t j = 0; j < c->nbh[i].size(); j++)
                {
                    Cell* cl = c->nbh[i][j];
                    for (std::size_t k = 0; k < cl->nbh[i].size(); k++)
                    {
                        if (cl->nbh[i][k] == c)
                        {
                            cl->nbh[i].erase(cl->nbh[i].begin() + k);
                            break;
                        }
                    }
                }
            }

            auto pos = GridNearest<_T, _TT>::hash_.find(&cell->coord);
            if (pos != GridNearest<_T, _TT>::hash_.end())
            {
                GridNearest<_T, _TT>::hash_.erase(pos);
                GridNearest<_T, _TT>::cnn_->remove(cell->coord);
                return true;
            }
            return false;
        }

        void add(BaseCell *cell) override
        {
            GridNearest<_T, _TT>::add(cell); 
            if (neighborCell_ > 0)
            {
                std::vector<CellArray> list;
                if (neighborCell_ == 1)
                {
                    list.resize(1);
                    list[0].push_back(static_cast<Cell *>(cell));
                }
                else 
                {
                    Coord test = cell->coord;
                    neighbors(test, neighborCell_, list);
                }

                Cell *c = static_cast<Cell *>(cell);
                for (std::size_t i = 1; i < list.size(); i++)
                {
                    for (std::size_t j = 0; j < list[i].size(); j++)
                    {
                        Cell *cl = list[i][j];
                        cl->nbh[i].push_back(c);
                    }
                }
                c->nbh.swap(list);
            }
        }

        void add(const _T &data, const Coord &coord) override
        {
            BaseCell *cell = GridNearest<_T, _TT>::getCell(coord);
            if (!cell)
            {
                cell = createCell(coord);
                add(cell);
            }
            GridNearest<_T, _TT>::add(data, cell);
        }

        /// Lazily remove a cell from the grid. That is, only update the neighbor list
        virtual bool lazyRemove(BaseCell *cell)
        {
            if (!cell)
                return false;
            Cell *c = static_cast<Cell *>(cell);
            if (c->removed)
                return false;
            for (auto & cl : c->nnbh)
            {
                cl->neighbors--;
                if (!cl->border && cl->neighbors < interiorCellNeighborsLimit_)
                    cl->border = true;
            }

            auto pos = GridNearest<_T, _TT>::hash_.find(&cell->coord);
            if (pos != GridNearest<_T, _TT>::hash_.end())
            {
                c->removed = true;
                removed_++;
                GridNearest<_T, _TT>::cnn_->remove(cell->coord);
                return true;
            }
            return false;
        }

        /// Lazily add a cell to the grid. That is, only update the neighbor list
        virtual bool lazyAdd(BaseCell *cell)
        {
            if (!cell)
                return false;
            Cell *c = static_cast<Cell *>(cell);
            if (!c->removed)
                return false;
            for (auto & cl : c->nnbh)
            {
                cl->neighbors++;
                if (cl->border && cl->neighbors >= interiorCellNeighborsLimit_)
                    cl->border = false;
            }

            auto pos = GridNearest<_T, _TT>::hash_.find(&cell->coord);
            if (pos != GridNearest<_T, _TT>::hash_.end())
            {
                c->removed = false;
                removed_--;
                GridNearest<_T, _TT>::cnn_->add(cell->coord);
                return true;
            }
            return false;
        }

        virtual void updateNbh(Cell *cell1, Cell *cell2)
        {
            if (cell1 == cell2)
                return;
            std::size_t n = (std::size_t)((cell1->coord - cell2->coord).template lpNorm<Eigen::Infinity>() + 1);
            if (n > cell1->nbh.size())
                cell1->nbh.resize(n);
            if (n > cell2->nbh.size())
                cell2->nbh.resize(n);
            bool found = false;
            if (!cell1->nbh[n-1].empty() && !cell2->nbh[n-1].empty())
            {
                for (auto & c : cell1->nbh[n-1])
                {
                    if (c == cell2)
                    {
                        found = true;
                        break;
                    }
                }
            }
            if (!found)
            {
                cell1->nbh[n-1].push_back(cell2);
                cell2->nbh[n-1].push_back(cell1);
            }
        }

        /// Get the set of instantiated cells in the grid
        void getCells(CellArray &cells) const
        {
            for (const auto &h : GridNearest<_T, _TT>::hash_)
                cells.push_back(static_cast<Cell *>(h.second));
        }

    protected:
        /// Compute how many sides of a coordinate touch the boundaries of the grid
        unsigned int numberOfBoundaryDimensions(const Coord &coord) const
        {
            unsigned int result = 0;
            if (hasBounds_)
            {
                for (unsigned int i = 0; i < GridNearest<_T, _TT>::dimension_; ++i)
                    if (coord[i] <= lowBound_[i] || coord[i] >= upBound_[i])
                        result++;
            }
            return result;
        }

        /// Flag indicating whether storing the neighboring cells
        bool useNeighbor_{true};

        /// Flag indicating whether bounds are in effect for this grid
        bool hasBounds_;

        /// If bounds are set, this defines the lower corner cell
        Coord lowBound_;

        /// If bounds are set, this defines the upper corner cell
        Coord upBound_;

        /// By default, cells are considered on the border if 2n
        /// neighbors are created, for a space of dimension n.
        /// this value is overridden and set in this member variable
        unsigned int interiorCellNeighborsLimit_;

        /// The removed cell number in this grid
        unsigned int removed_{0};

        int neighborCell_{-1};

        /// Flag indicating whether the neighbor count used to determine whether
        /// a cell is on the border or not
        bool overrideCellNeighborsLimit_;

    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    };
}

#endif
