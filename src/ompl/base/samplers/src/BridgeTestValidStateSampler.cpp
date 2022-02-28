#include "ompl/base/samplers/BridgeTestValidStateSampler.h"
#include "ompl/base/SpaceInformation.h"
#include "ompl/tools/config/MagicConstants.h"

ompl::base::BridgeTestValidStateSampler::BridgeTestValidStateSampler(const SpaceInformation *si)
  : ValidStateSampler(si)
  , sampler_(si->allocStateSampler())
  , stddev_(si->getMaximumExtent() * magic::STD_DEV_AS_SPACE_EXTENT_FRACTION)
{
    name_ = "bridge_test";
    params_.declareParam<double>("standard_deviation", [this](double stddev) { setStdDev(stddev); },
                                 [this] { return getStdDev(); });
}

bool ompl::base::BridgeTestValidStateSampler::sample(State *state)
{
    oTime_ = 0;
    unsigned int attempts = 0;
    bool valid = false;
    State *endpoint = si_->allocState();
    do
    {
        sampler_->sampleUniform(state);
        time::point starto = time::now();
        bool v1 = si_->isValid(state);
        oTime_ += time::seconds(time::now() - starto);
        if (!v1)
        {
            sampler_->sampleGaussian(endpoint, state, stddev_);
            starto = time::now();
            bool v2 = si_->isValid(endpoint);
            oTime_ += time::seconds(time::now() - starto);
            if (!v2)
            {
                si_->getStateSpace()->interpolate(endpoint, state, 0.5, state);
                starto = time::now();
                valid = si_->isValid(state);
                oTime_ += time::seconds(time::now() - starto);
            }
        }
        ++attempts;
    } while (!valid && attempts < attempts_);

    si_->freeState(endpoint);
    return valid;
}

bool ompl::base::BridgeTestValidStateSampler::sampleNear(State *state, const State *near, const double distance)
{
    oTime_ = 0;
    unsigned int attempts = 0;
    bool valid = false;
    State *endpoint = si_->allocState();
    do
    {
        sampler_->sampleUniformNear(state, near, distance);
        time::point starto = time::now();
        bool v1 = si_->isValid(state);
        oTime_ += time::seconds(time::now() - starto);
        if (!v1)
        {
            sampler_->sampleGaussian(endpoint, state, distance);
            starto = time::now();
            bool v2 = si_->isValid(endpoint);
            oTime_ += time::seconds(time::now() - starto);
            if (!v2)
            {
                si_->getStateSpace()->interpolate(endpoint, state, 0.5, state);
                starto = time::now();
                valid = si_->isValid(state);
                oTime_ += time::seconds(time::now() - starto);
            }
        }
        ++attempts;
    } while (!valid && attempts < attempts_);

    si_->freeState(endpoint);
    return valid;
}
