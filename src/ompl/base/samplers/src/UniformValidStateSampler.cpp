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

/* Author: Ioan Sucan */

#include "ompl/base/samplers/UniformValidStateSampler.h"
#include "ompl/base/SpaceInformation.h"

ompl::base::UniformValidStateSampler::UniformValidStateSampler(const SpaceInformation *si)
  : ValidStateSampler(si), sampler_(si->allocStateSampler())
{
    name_ = "uniform";
}

bool ompl::base::UniformValidStateSampler::sample(State *state)
{
    oTime_ = 0;
    unsigned int attempts = 0;
    bool valid = false;
    do
    {
        sampler_->sampleUniform(state);
        time::point starto = time::now();
        valid = si_->isValid(state);
        oTime_ += time::seconds(time::now() - starto);
        ++attempts;
    } while (!valid && attempts < attempts_);
    return valid;
}

bool ompl::base::UniformValidStateSampler::sampleNear(State *state, const State *near, const double distance)
{
    oTime_ = 0;
    unsigned int attempts = 0;
    bool valid = false;
    do
    {
        sampler_->sampleUniformNear(state, near, distance);
        time::point starto = time::now();
        valid = si_->isValid(state);
        oTime_ += time::seconds(time::now() - starto);
        ++attempts;
    } while (!valid && attempts < attempts_);
    return valid;
}
