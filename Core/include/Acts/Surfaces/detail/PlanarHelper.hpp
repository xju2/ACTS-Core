// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/////////////////////////////////////////////////////////////////
// PlanarHelper.ipp, Acts project
///////////////////////////////////////////////////////////////////
#pragma once

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/BoundaryCheck.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Intersection.hpp"

namespace Acts {

namespace detail {

/// Straight line intersection from position and momentum to planar surfaces
///
/// @param surface The surface in question to be intersected
/// @param gctx The current geometry context object, e.g. alignment
/// @param gpos global 3D position - considered to be on surface but not
///        inside bounds (check is done)
/// @param gdir 3D direction representation - expected to be normalized
///        (no check done)
/// @param bwdTolerance a tolerance for which an intersection is accepted
///        in opposite direction
/// @param bcheck boundary check directive for this operation
/// @param corr is a correction function on position and momentum to do
///        a more appropriate intersection
///
/// @return Intersection object
static Intersection planarIntersectionEstimate(
    const Surface& surface, const GeometryContext& gctx, const Vector3D& gpos,
    const Vector3D& gdir, const BoundaryCheck& bcheck, double bwdTolerance,
    CorrFnc correct = nullptr) {
  // minimize the call to transform()
  const auto& tMatrix = surface.transform(gctx).matrix();
  const Vector3D pnormal = tMatrix.block<3, 1>(0, 2).transpose();
  const Vector3D pcenter = tMatrix.block<3, 1>(0, 3).transpose();
  // return solution and path
  Vector3D solution(0., 0., 0.);
  double path = std::numeric_limits<double>::infinity();
  // lemma : the solver -> should catch current values
  auto solve = [&solution, &path, &pnormal, &pcenter, &bwdTolerance](
                   const Vector3D& lpos,
                   const Vector3D& ldir) -> IntersectionStatus {
    double denom = ldir.dot(pnormal);
    if (denom != 0.0) {
      path = (pnormal.dot((pcenter - lpos))) / (denom);
      solution = (lpos + path * ldir);
    }
    if (path > 0.) {
      return IntersectionStatus::reachable;
    }
    if (path * path < s_onSurfaceTolerance * s_onSurfaceTolerance) {
      return IntersectionStatus::onSurface;
    }
    // is valid if it goes into the right direction, or is within
    return (path * path < bwdTolerance * bwdTolerance
                ? IntersectionStatus::overstepped
                : IntersectionStatus::unreachable);
  };
  // solve first without corrector
  auto istatus = solve(gpos, gdir);
  // if configured to correct, do it and solve again
  if (correct) {
    // copy as the corrector may change them
    Vector3D lposc = gpos;
    Vector3D ldirc = gdir;
    if (correct(lposc, ldirc, path)) {
      istatus = solve(lposc, ldirc);
    }
  }
  // Evaluate boundaries if necessary, since the solution was done in global
  // frame a global to local is necessary for the boundary check
  if (bcheck && istatus != IntersectionStatus::unreachable) {
    istatus = surface.isOnSurface(gctx, solution, gdir, bcheck)
                  ? istatus
                  : IntersectionStatus::missed;
  }
  // return the result
  return Intersection(solution, path, istatus);
}

}  // namespace detail
}  // namespace Acts