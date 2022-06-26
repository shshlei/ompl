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

/* Authors: Ioan Sucan, Shi Shenglei */

#ifndef OMPL_DATASTRUCTURES_GRID_
#define OMPL_DATASTRUCTURES_GRID_

#include <Eigen/Core>
#include <vector>
#include <iostream>
#include <cstdlib>
#include <unordered_map>
#include <algorithm>
#include "ompl/util/Exception.h"

namespace ompl
{
    /** \brief Representation of a simple grid */
    template <typename _T>
    class Grid
    {
    public:
        /// Definition of a coordinate within this grid
        using Coord = Eigen::VectorXi;

        /// Definition of a cell in this grid
        struct Cell
        {
            /// The data we store in the cell
            _T data;

            /// The coordinate of the cell
            Coord coord;

            Cell() = default;

            virtual ~Cell() = default;

            EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        };

        /// The datatype for arrays of cells
        using CellArray = std::vector<Cell *>;

        /// The constructor takes the dimension of the grid as argument
        explicit Grid(unsigned int dimension)
        {
            setDimension(dimension);
        }

        /// Destructor
        virtual ~Grid()
        {
            freeMemory();
        }

        /// Clear all cells in the grid
        virtual void clear()
        {
            freeMemory();
        }

        /// Return the dimension of the grid
        unsigned int getDimension() const
        {
            return dimension_;
        }

        /// Update the dimension of the grid; this should not be done
        /// unless the grid is empty
        void setDimension(unsigned int dimension)
        {
            if (!empty())
                throw;
            dimension_ = dimension;
            maxNeighbors_ = 2 * dimension_;
        }

        /// Check if a cell exists at the specified coordinate
        bool has(const Coord &coord) const
        {
            return getCell(coord) != nullptr;
        }

        /// Get the cell at a specified coordinate
        Cell *getCell(const Coord &coord) const
        {
            auto pos = hash_.find(const_cast<Coord *>(&coord));
            Cell *c = (pos != hash_.end()) ? pos->second : nullptr;
            return c;
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
            list.reserve(list.size() + maxNeighbors_);

            for (int i = dimension_ - 1; i >= 0; --i)
            {
                coord[i]--;

                auto pos = hash_.find(&coord);
                Cell *cell = (pos != hash_.end()) ? pos->second : nullptr;

                if (cell)
                    list.push_back(cell);
                coord[i] += 2;

                pos = hash_.find(&coord);
                cell = (pos != hash_.end()) ? pos->second : nullptr;

                if (cell)
                    list.push_back(cell);
                coord[i]--;
            }
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
            assert(dimension_ <= 9);
            if (!ring)
            {
                int s = 0;
                for (int i = 0; i < n; i++)
                {
                    int e = std::pow(2*i + 1, dimension_);
                    list[i].reserve(e - s);
                    s = e;
                }
            }
            else
            {
                int i = n - 2;
                int s = std::pow(2*i + 1, dimension_);
                i = n - 1;
                int e = std::pow(2*i + 1, dimension_);
                list[i].reserve(e - s);
            }

            if (dimension_ == 1)
                neighbors1(coord, n, list, ring);
            else if (dimension_ == 2)
                neighbors2(coord, n, list, ring);
            else if (dimension_ == 3)
                neighbors3(coord, n, list, ring);
            else if (dimension_ == 4)
                neighbors4(coord, n, list, ring);
            else if (dimension_ == 5)
                neighbors5(coord, n, list, ring);
            else if (dimension_ == 6)
                neighbors6(coord, n, list, ring);
            else if (dimension_ == 7)
                neighbors7(coord, n, list, ring);
            else if (dimension_ == 8)
                neighbors8(coord, n, list, ring);
            else if (dimension_ == 9)
                neighbors9(coord, n, list, ring);
            else 
                throw Exception("Grid: Neighbor function is not implemented yet for that high dimension!");
        }

        /// Get the connected components formed by the cells in this grid (based on neighboring relation)
        std::vector<std::vector<Cell *>> components() const
        {
            using ComponentHash = std::unordered_map<Coord *, int, HashFunCoordPtr, EqualCoordPtr>;

            int components = 0;
            ComponentHash ch;
            std::vector<std::vector<Cell *>> res;

            for (auto & i: hash_)
            {
                Cell *c0 = i.second;
                auto pos = ch.find(&c0->coord);
                int comp = (pos != ch.end()) ? pos->second : -1;

                if (comp < 0)
                {
                    res.resize(res.size() + 1);
                    std::vector<Cell *> &q = res.back();
                    q.push_back(c0);
                    std::size_t index = 0;
                    while (index < q.size())
                    {
                        Cell *c = q[index++];
                        pos = ch.find(&c->coord);
                        comp = (pos != ch.end()) ? pos->second : -1;

                        if (comp < 0)
                        {
                            ch.insert(std::make_pair(&c->coord, components));
                            std::vector<Cell *> nbh;
                            neighbors(c, nbh);
                            for (const auto &n : nbh)
                            {
                                pos = ch.find(&n->coord);
                                comp = (pos != ch.end()) ? pos->second : -1;
                                if (comp < 0)
                                    q.push_back(n);
                            }
                        }
                        else
                        {
                            --index;
                            q.erase(q.begin() + index);
                        }
                    }
                    ++components;
                }
            }
            std::sort(res.begin(), res.end(), SortComponents());
            return res;
        }

        /// Instantiate a new cell at given coordinates; optionally
        /// Return the list of future neighbors.
        /// Note: this call only creates the cell, but does not add it to the grid.
        /// It however updates the neighbor count for neighboring cells
        virtual Cell *createCell(const Coord &coord, CellArray *nbh = nullptr)
        {
            auto *cell = new Cell();
            cell->coord = coord;
            if (nbh)
                neighbors(cell->coord, *nbh);
            return cell;
        }

        /// Remove a cell from the grid. If the cell has not been
        /// Added to the grid, only update the neighbor list
        virtual bool remove(Cell *cell)
        {
            if (cell)
            {
                auto pos = hash_.find(&cell->coord);
                if (pos != hash_.end())
                {
                    hash_.erase(pos);
                    return true;
                }
            }
            return false;
        }

        /// Add an instantiated cell to the grid
        virtual void add(Cell *cell)
        {
            hash_.insert(std::make_pair(&cell->coord, cell));
        }

        /// Clear the memory occupied by a cell; do not call this function unless remove() was called first
        virtual void destroyCell(Cell *cell) const
        {
            delete cell;
        }

        /// Get the data stored in the cells we are aware of
        void getContent(std::vector<_T> &content) const
        {
            for (const auto &h : hash_)
                content.push_back(h.second->data);
        }

        /// Get the set of coordinates where there are cells
        void getCoordinates(std::vector<Coord *> &coords) const
        {
            for (const auto &h : hash_)
                coords.push_back(h.first);
        }

        /// Get the set of instantiated cells in the grid
        void getCells(CellArray &cells) const
        {
            for (const auto &h : hash_)
                cells.push_back(h.second);
        }

        /// Print the value of a coordinate to a stream
        void printCoord(Coord &coord, std::ostream &out = std::cout) const
        {
            out << "[ ";
            for (unsigned int i = 0; i < dimension_; ++i)
                out << coord[i] << " ";
            out << "]" << std::endl;
        }

        /// Check if the grid is empty
        bool empty() const
        {
            return hash_.empty();
        }

        /// Check the size of the grid
        unsigned int size() const
        {
            return hash_.size();
        }

        /// Print information about the data in this grid structure
        virtual void status(std::ostream &out = std::cout) const
        {
            out << size() << " total cells " << std::endl;
            const std::vector<std::vector<Cell *>> &comp = components();
            out << comp.size() << " connected components: ";
            for (const auto &c : comp)
                out << c.size() << " ";
            out << std::endl;
        }

    protected:
        /// Free the allocated memory
        void freeMemory()
        {
            CellArray content;
            getCells(content);
            hash_.clear();

            for (auto &c : content)
                delete c;
        }

        /// Hash function for coordinates; see
        /// http://www.cs.hmc.edu/~geoff/classes/hmc.cs070.200101/homework10/hashfuncs.html
        struct HashFunCoordPtr
        {
            /// Hash function for coordinates
            std::size_t operator()(const Coord *const s) const
            {
                unsigned long h = 0;
                for (int i = s->size() - 1; i >= 0; --i)
                {
                    int high = h & 0xf8000000;
                    h = h << 5;
                    h = h ^ (high >> 27);
                    h = h ^ std::abs((*s)[i]);
                }
                return (std::size_t)h;
            }
        };

        /// Equality operator for coordinate pointers
        struct EqualCoordPtr
        {
            /// Equality operator for coordinate pointers
            bool operator()(const Coord *const c1, const Coord *const c2) const
            {
                return *c1 == *c2;
            }
        };

        /// Define the datatype for the used hash structure
        using CoordHash = std::unordered_map<Coord *, Cell *, HashFunCoordPtr, EqualCoordPtr>;

        /// Helper to sort components by size
        struct SortComponents
        {
            /// Helper to sort components by size
            bool operator()(const std::vector<Cell *> &a, const std::vector<Cell *> &b) const
            {
                return a.size() > b.size();
            }
        };

    public:
        /// We only allow const iterators
        using iterator = typename CoordHash::const_iterator;

        /// Return the begin() iterator for the grid
        iterator begin() const
        {
            return hash_.begin();
        }

        /// Return the end() iterator for the grid
        iterator end() const
        {
            return hash_.end();
        }

        iterator erase(iterator it)
        {
            it = hash_.erase(it);
            return it;
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
        /// The dimension of the grid
        unsigned int dimension_;

        /// The maximum number of neighbors a cell can have (2 * dimension)
        unsigned int maxNeighbors_;

        /// The hash holding the cells
        CoordHash hash_;
    };
}  // namespace ompl

#endif
