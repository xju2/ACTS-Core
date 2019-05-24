// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Definitions.hpp"

namespace Acts {

/// @enum TargetStatus
///
/// These are used to identify the current
/// status of the navigation status
enum TargetStatus : int {
  undefined = -1,
  unreachable = 0,
  missed = 0,
  approaching = 1,
  overstepped = 2,
  hit = 3
};

/// Tests if the state reached a surface for single component
///
/// @tparam state_t is the stepper state type
///
/// @param [in,out] state State of the Stepper, may be changed
/// @param [in] surface Surface that is tested
/// @param [in] surface position that is used for testing
/// @param [in] surface direction that is used for testing
/// @param [in] bcheck is the BoundaryCheck for this directive
///
/// @return TargetStatus enum
template <typename state_t>
TargetStatus surfaceSatus(state_t& state, const Surface& surface,
                          const Vector3D& position, const Vector3D& direction,
                          const BoundaryCheck& bcheck,
                          bool updateStepSize = true) const {
  /// intersect, and allow a on surface tolerance
  auto intersect = surface.intersectionEstimate(state.geoContext, pposition,
                                                direction, navDir, bcheck);
  // The surface is reached within tolerance
  if (intersect and intersect.pathLength * intersect.pathLength <
                        s_onSurfaceTolerance * s_onSurfaceTolerance) {
    return SurfaceTarget::onSurface;
  } else if (intersect) {
    // The surface intersection has one solution
    state.stepSize.update(intersect.pathLength, cstep::actor, true);
    return intersect.pathLength > 0 ? SurfaceTarget::onApproach
                                    : SurfaceTarget::overstepped;
  }
  return SurfaceTarget::missed;
}

}  // namespace Acts