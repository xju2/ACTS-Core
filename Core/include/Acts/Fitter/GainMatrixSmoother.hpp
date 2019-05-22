// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <boost/range/adaptors.hpp>
#include <memory>
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/TrackParameters.hpp"

namespace Acts {

/// @brief Kalman smoother implementation based on Gain matrix formalism
///
/// @tparam parameters_t Type of the track parameters
/// @tparam jacobian_t Type of the Jacobian
template <typename parameters_t>
class GainMatrixSmoother {
  using jacobian_t = typename parameters_t::CovMatrix_t;

 public:
  /// @brief Gain Matrix smoother implementation
  ///

  template <typename source_link_t>
  parameters_t operator()(const GeometryContext& gctx,
                          MultiTrajectory<source_link_t>& trajectory,
                          size_t entryIndex) const {
    using namespace boost::adaptors;

    // using ParVector_t = typename parameters_t::ParVector_t;
    using CovMatrix_t = typename parameters_t::CovMatrix_t;
    using gain_matrix_t = CovMatrix_t;

    // For the last state: smoothed is filtered - also: switch to next
    auto prev_ts = trajectory.getTrackState(entryIndex);

    prev_ts.smoothed() = prev_ts.filtered();
    prev_ts.smoothedCovariance() = prev_ts.filteredCovariance();

    // Smoothing gain matrix
    gain_matrix_t G;

    trajectory.applyBackwards(prev_ts.previous(), [&prev_ts, &G](auto ts) {
      // should have filtered and predicted, this should also include the
      // covariances.
      assert(ts.hasFiltered());
      assert(ts.hasPredicted());
      assert(ts.hasJacobian());

      // previous trackstate should have smoothed and predicted
      assert(prev_ts.hasSmoothed());
      assert(prev_ts.hasPredicted());

      // Gain smoothing matrix
      G = ts.filteredCovariance() * ts.jacobian().transpose() *
          prev_ts.predictedCovariance().inverse();

      // Calculate the smoothed parameters
      ts.smoothed() =
          ts.filtered() + G * (prev_ts.smoothed() - prev_ts.predicted());

      // And the smoothed covariance
      ts.smoothedCovariance() =
          ts.filteredCovariance() -
          G * (prev_ts.predictedCovariance() - prev_ts.smoothedCovariance()) *
              G.transpose();

      prev_ts = ts;
    });

    // construct parameters from last track state
    parameters_t lastSmoothed(
        gctx, std::make_unique<CovMatrix_t>(prev_ts.smoothedCovariance()),
        prev_ts.smoothed(), prev_ts.referenceSurface().getSharedPtr());

    return lastSmoothed;
  }
};
}  // namespace Acts
