/*********************************************************************
* Software License Agreement (BSD License)
*
*  Copyright (c) 2008, Willow Garage, Inc.
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
*   * Neither the name of the Willow Garage nor the names of its
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

#ifndef OMPL_BASE_CONTACT_
#define OMPL_BASE_CONTACT_

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <vector>
#include <memory>
#include <map>
#include <unordered_map>

namespace ompl
{
    namespace base
    {
        struct ContactResult
        {
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            std::string link_names[2];

            Eigen::Vector3d nearest_points[2];

            Eigen::Vector3d nearest_points_local[2];

            Eigen::Vector3d nearest_points_local2[2];

            Eigen::Isometry3d transform[2];

            Eigen::Vector3d normal;

            Eigen::Vector3d center[2];

            Eigen::Vector3d local_center[2];

            double radius[2];

            double distance;

            int type_id[2];

            int shape_id[2];

            int subshape_id[2];

            bool has_sphere[2];

            ContactResult() { clear(); }

            void clear()
            {
                distance = std::numeric_limits<double>::infinity();
                radius[0] = radius[1] = -1.0;
                center[0].setZero();
                center[1].setZero();
                local_center[0].setZero();
                local_center[1].setZero();
                has_sphere[0] = has_sphere[1] = false;
                nearest_points[0].setZero();
                nearest_points[1].setZero();
                nearest_points_local[0].setZero();
                nearest_points_local[1].setZero();
                nearest_points_local2[0].setZero();
                nearest_points_local2[1].setZero();
                transform[0] = Eigen::Isometry3d::Identity();
                transform[1] = Eigen::Isometry3d::Identity();
                link_names[0] = "";
                link_names[1] = "";
                shape_id[0] = -1;
                shape_id[1] = -1;
                subshape_id[0] = -1;
                subshape_id[1] = -1;
                type_id[0] = 0;
                type_id[1] = 0;
                normal.setZero();
            }
        };

        template <typename T>
        using AlignedVector = std::vector<T, Eigen::aligned_allocator<T>>;

        template <typename Key, typename Value>
        using AlignedMap = std::map<Key, Value, std::less<Key>, Eigen::aligned_allocator<std::pair<const Key, Value>>>;

        template <typename Key, typename Value>
        using AlignedUnorderedMap = std::unordered_map<Key, Value, std::hash<Key>, std::equal_to<Key>, Eigen::aligned_allocator<std::pair<const Key, Value>>>;

        using ContactResultVector = AlignedVector<ContactResult>;
        using ContactResultMap = AlignedMap<std::pair<std::string, std::string>, ContactResultVector>;

        inline std::size_t flattenMoveResults(ContactResultMap&& m, ContactResultVector& v)
        {
            v.clear();
            v.reserve(m.size());
            for (const auto& mv : m)
                std::move(mv.second.begin(), mv.second.end(), std::back_inserter(v));

            return v.size();
        }

        inline std::size_t flattenCopyResults(const ContactResultMap& m, ContactResultVector& v)
        {
            v.clear();
            v.reserve(m.size());
            for (const auto& mv : m)
                std::copy(mv.second.begin(), mv.second.end(), std::back_inserter(v));

            return v.size();
        }
    }
}
#endif
