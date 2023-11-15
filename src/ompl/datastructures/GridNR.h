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

/* Authors: Shi Shenglei */

#ifndef OMPL_DATASTRUCTURES_GRID_NR_
#define OMPL_DATASTRUCTURES_GRID_NR_

#include "ompl/datastructures/Grid.h"
#include <numeric>

namespace ompl
{
    /** \brief Representation of a grid where cells keep track of how many neighbors they have in a ring */
    template <typename _T>
    class GridNR : public Grid<_T>
    {
    public:
        /// Datatype for cell in base class
        using BaseCell = typename Grid<_T>::Cell;

        /// Datatype for array of cells in base class
        using BaseCellArray = typename Grid<_T>::CellArray;

        /// Datatype for cell coordinates
        using Coord = typename Grid<_T>::Coord;

        /// Definition of a cell in this grid
        struct Cell : public BaseCell
        {
            std::vector<std::vector<Cell *>> nbh;

            Cell() = default;

            ~Cell() override = default;

            EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        };

        /// The datatype for arrays of cells
        using CellArray = std::vector<Cell *>;

        /// The constructor takes the dimension of the grid as argument
        explicit GridNR(unsigned int dimension) : Grid<_T>(dimension)
        {
            Grid<_T>::setDimension(dimension);
        }

        ~GridNR() override = default;

        /// Store the neighbors of a cell in (2*n - 1) grid 
        void setNeighborCell(int n)
        {
            neighborCell_ = n;
            maxneighborCell_ = 2 * n;
        }

        /// Store the neighbors of a cell in (2*n - 1) grid 
        void setMaxNeighborCell(int n)
        {
            maxneighborCell_ = n;
        }

        /// Get the cell at a specified coordinate
        Cell *getCell(const Coord &coord) const
        {
            return static_cast<Cell *>(Grid<_T>::getCell(coord));
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
            Grid<_T>::neighbors(coord, baselist);
            list.reserve(list.size() + baselist.size());
            for (const auto &c : baselist)
                list.push_back(static_cast<Cell *>(c));
        }

        /// Get the list of neighbors for a given cell in (2*n - 1) cell neighborhood 
        void neighbors(const Cell *cell, int n, std::vector<CellArray> &list, bool ring = false) const
        {
            if (n == 1)
            {
                list.resize(n);
                list[0].push_back(cell);
                return;
            }
            Coord test = cell->coord;
            neighbors(test, n, list, ring);
        }

        /// Get the list of neighbors for a given cell in (2*n - 1) cell neighborhood
        void neighbors(const Coord &coord, int n, std::vector<CellArray> &list, bool ring = false) const
        {
            Coord test = coord;
            neighbors(test, n, list, ring);
        }

        /// Get the list of neighbors for a given cell in (2*n - 1) cell neighborhood
        void neighbors(Coord &coord, int n, std::vector<CellArray> &list, bool ring = false) const
        {
            list.resize(n);
            if (n == 1)
            {
                Cell *cell = getCell(coord);
                if (cell)
                    list[0].push_back(cell);
                return;
            }
            assert(Grid<_T>::dimension_ <= 9);
            if (!ring)
            {
                int s = 0;
                for (int i = 0; i < n; i++)
                {
                    int e = std::pow(2*i + 1, Grid<_T>::dimension_);
                    list[i].reserve(e - s);
                    s = e;
                }
            }
            else
            {
                int i = n - 2;
                int s = std::pow(2*i + 1, Grid<_T>::dimension_);
                i = n - 1;
                int e = std::pow(2*i + 1, Grid<_T>::dimension_);
                list[i].reserve(e - s);
            }

            if (Grid<_T>::dimension_ == 1)
                neighbors1(coord, n, list, ring);
            else if (Grid<_T>::dimension_ == 2)
                neighbors2(coord, n, list, ring);
            else if (Grid<_T>::dimension_ == 3)
                neighbors3(coord, n, list, ring);
            else if (Grid<_T>::dimension_ == 4)
                neighbors4(coord, n, list, ring);
            else if (Grid<_T>::dimension_ == 5)
                neighbors5(coord, n, list, ring);
            else if (Grid<_T>::dimension_ == 6)
                neighbors6(coord, n, list, ring);
            else if (Grid<_T>::dimension_ == 7)
                neighbors7(coord, n, list, ring);
            else if (Grid<_T>::dimension_ == 8)
                neighbors8(coord, n, list, ring);
            else if (Grid<_T>::dimension_ == 9)
                neighbors9(coord, n, list, ring);
            else 
                throw Exception("Grid: Neighbor function is not implemented yet for that high dimension!");
        }

        virtual Cell *createCell(const Coord &coord)
        {
            Cell *cell = new Cell();
            cell->coord = coord;
            return cell;
        }

        /// Remove a cell from the grid. If the cell has not been
        /// Added to the grid, only update the neighbor list
        virtual bool remove(Cell *cell)
        {
            if (!cell)
                return false;
            for (std::size_t i = 1; i < cell->nbh.size(); i++)
            {
                for (std::size_t j = 0; j < cell->nbh[i].size(); j++)
                {
                    Cell* cl = cell->nbh[i][j];
                    for (std::size_t k = 0; k < cl->nbh[i].size(); k++)
                    {
                        if (cl->nbh[i][k] == cell)
                        {
                            std::iter_swap(cl->nbh[i].begin() + k, cl->nbh[i].end() - 1);
                            cl->nbh[i].pop_back();
                            break;
                        }
                    }
                }
            }
            auto pos = Grid<_T>::hash_.find(&cell->coord);
            if (pos != Grid<_T>::hash_.end())
            {
                Grid<_T>::hash_.erase(pos);
                return true;
            }
            return false;
        }

        virtual void add(Cell *cell)
        {
            Grid<_T>::add(cell); 
            std::vector<CellArray> list;
            if (neighborCell_ == 1)
            {
                list.resize(1);
                list[0].push_back(cell);
            }
            else 
            {
                Coord test = cell->coord;
                neighbors(test, neighborCell_, list);
            }
            for (std::size_t i = 1; i < list.size(); i++)
            {
                for (std::size_t j = 0; j < list[i].size(); j++)
                {
                    Cell *cl = list[i][j];
                    cl->nbh[i].push_back(cell);
                }
            }
            cell->nbh.swap(list);
        }

        virtual void updateNbh(Cell *cell1, Cell *cell2)
        {
            if (cell1 == cell2)
                return;
            std::size_t n = (std::size_t)((cell1->coord - cell2->coord).template lpNorm<Eigen::Infinity>() + 1);
            if (n > (std::size_t)maxneighborCell_)
                return;
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
            for (const auto &h : Grid<_T>::hash_)
                cells.push_back(static_cast<Cell *>(h.second));
        }

    protected:
        void neighbors1(Coord &coord, int n, std::vector<CellArray> &list, bool ring = false, int ii = 0, int pabs = 0) const
        {
            if (!ring || pabs == n-1)
            {
                Cell *cell = getCell(coord);
                if (cell)
                    list[pabs].push_back(cell);

                for (int i = 1; i < n; i++)
                {
                    int ppabs = std::max(i, pabs);

                    coord[ii] += i;
                    Cell *cell = getCell(coord);
                    if (cell)
                        list[ppabs].push_back(cell);

                    coord[ii] -= 2*i;
                    cell = getCell(coord);
                    if (cell)
                        list[ppabs].push_back(cell);

                    coord[ii] += i;
                }
            }
            else
            {
                int i = n - 1;
                int ppabs = i;

                coord[ii] += i;
                Cell *cell = getCell(coord);
                if (cell)
                    list[ppabs].push_back(cell);

                coord[ii] -= 2*i;
                cell = getCell(coord);
                if (cell)
                    list[ppabs].push_back(cell);

                coord[ii] += i;
            }
        }

        void neighbors2(Coord &coord, int n, std::vector<CellArray> &list, bool ring = false, int ii = 0, int jj = 1, int pabs = 0) const
        {
            neighbors1(coord, n, list, ring, jj, pabs);

            for (int i = 1; i < n; i++)
            {
                int ppabs = std::max(i, pabs);

                coord[ii] += i;
                neighbors1(coord, n, list, ring, jj, ppabs);

                coord[ii] -= 2*i;
                neighbors1(coord, n, list, ring, jj, ppabs);

                coord[ii] += i;
            }
        }

        void neighbors3(Coord &coord, int n, std::vector<CellArray> &list, bool ring = false, int ii = 0, int jj = 1, int kk = 2, int pabs = 0) const
        {
            neighbors2(coord, n, list, ring, jj, kk, pabs);

            for (int i = 1; i < n; i++)
            {
                int ppabs = std::max(i, pabs);

                coord[ii] += i;
                neighbors2(coord, n, list, ring, jj, kk, ppabs);

                coord[ii] -= 2*i;
                neighbors2(coord, n, list, ring, jj, kk, ppabs);

                coord[ii] += i;
            }
        }

        void neighbors4(Coord &coord, int n, std::vector<CellArray> &list, bool ring = false, int ii = 0, int jj = 1, int kk = 2, int ll = 3, int pabs = 0) const
        {
            neighbors3(coord, n, list, ring, jj, kk, ll, pabs);

            for (int i = 1; i < n; i++)
            {
                int ppabs = std::max(i, pabs);

                coord[ii] += i;
                neighbors3(coord, n, list, ring, jj, kk, ll, ppabs);

                coord[ii] -= 2*i;
                neighbors3(coord, n, list, ring, jj, kk, ll, ppabs);

                coord[ii] += i;
            }
        }

        void neighbors5(Coord &coord, int n, std::vector<CellArray> &list, bool ring = false, int ii = 0, int jj = 1, int kk = 2, int ll = 3, int mm = 4, int pabs = 0) const
        {
            neighbors4(coord, n, list, ring, jj, kk, ll, mm, pabs);

            for (int i = 1; i < n; i++)
            {
                int ppabs = std::max(i, pabs);

                coord[ii] += i;
                neighbors4(coord, n, list, ring, jj, kk, ll, mm, ppabs);

                coord[ii] -= 2*i;
                neighbors4(coord, n, list, ring, jj, kk, ll, mm, ppabs);

                coord[ii] += i;
            }
        }

        void neighbors6(Coord &coord, int n, std::vector<CellArray> &list, bool ring = false, int ii = 0, int jj = 1, int kk = 2, int ll = 3, int mm = 4, int nn = 5, int pabs = 0) const
        {
            neighbors5(coord, n, list, ring, jj, kk, ll, mm, nn, pabs);

            for (int i = 1; i < n; i++)
            {
                int ppabs = std::max(i, pabs);

                coord[ii] += i;
                neighbors5(coord, n, list, ring, jj, kk, ll, mm, nn, ppabs);

                coord[ii] -= 2*i;
                neighbors5(coord, n, list, ring, jj, kk, ll, mm, nn, ppabs);

                coord[ii] += i;
            }
        }

        void neighbors7(Coord &coord, int n, std::vector<CellArray> &list, bool ring = false, int ii = 0, int jj = 1, int kk = 2, int ll = 3, int mm = 4, int nn = 5, int oo = 6, int pabs = 0) const
        {
            neighbors6(coord, n, list, ring, jj, kk, ll, mm, nn, oo, pabs);

            for (int i = 1; i < n; i++)
            {
                int ppabs = std::max(i, pabs);

                coord[ii] += i;
                neighbors6(coord, n, list, ring, jj, kk, ll, mm, nn, oo, ppabs);

                coord[ii] -= 2*i;
                neighbors6(coord, n, list, ring, jj, kk, ll, mm, nn, oo, ppabs);

                coord[ii] += i;
            }
        }

        void neighbors8(Coord &coord, int n, std::vector<CellArray> &list, bool ring = false, int ii = 0, int jj = 1, int kk = 2, int ll = 3, int mm = 4, int nn = 5, int oo = 6, int pp = 7, int pabs = 0) const
        {
            neighbors7(coord, n, list, ring, jj, kk, ll, mm, nn, oo, pp, pabs);

            for (int i = 1; i < n; i++)
            {
                int ppabs = std::max(i, pabs);

                coord[ii] += i;
                neighbors7(coord, n, list, ring, jj, kk, ll, mm, nn, oo, pp, ppabs);

                coord[ii] -= 2*i;
                neighbors7(coord, n, list, ring, jj, kk, ll, mm, nn, oo, pp, ppabs);

                coord[ii] += i;
            }
        }

        void neighbors9(Coord &coord, int n, std::vector<CellArray> &list, bool ring = false, int ii = 0, int jj = 1, int kk = 2, int ll = 3, int mm = 4, int nn = 5, int oo = 6, int pp = 7, int qq = 8, int pabs = 0) const
        {
            neighbors8(coord, n, list, ring, jj, kk, ll, mm, nn, oo, pp, qq, pabs);

            for (int i = 1; i < n; i++)
            {
                int ppabs = std::max(i, pabs);

                coord[ii] += i;
                neighbors8(coord, n, list, ring, jj, kk, ll, mm, nn, oo, pp, qq, ppabs);

                coord[ii] -= 2*i;
                neighbors8(coord, n, list, ring, jj, kk, ll, mm, nn, oo, pp, qq, ppabs);

                coord[ii] += i;
            }
        }

    protected:

        int neighborCell_{1};

        int maxneighborCell_{2};

    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    };
}

#endif
