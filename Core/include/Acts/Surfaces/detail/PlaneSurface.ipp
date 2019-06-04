// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// PlaneSurface.ipp, Acts project
///////////////////////////////////////////////////////////////////

inline const Vector3D PlaneSurface::normal(const GeometryContext& gctx,
                                           const Vector2D& /*lpos*/) const {
  // fast access via tranform matrix (and not rotation())
  const auto& tMatrix = transform(gctx).matrix();
  return Vector3D(tMatrix(0, 2), tMatrix(1, 2), tMatrix(2, 2));
}

inline const Vector3D PlaneSurface::binningPosition(
    const GeometryContext& gctx, BinningValue /*bValue*/) const {
  return center(gctx);
}

inline double PlaneSurface::pathCorrection(const GeometryContext& gctx,
                                           const Vector3D& gpos,
                                           const Vector3D& gmom) const {
  /// we can ignore the global position here
  return 1. / std::abs(Surface::normal(gctx, gpos).dot(gmom.normalized()));
}

inline Intersection PlaneSurface::intersectionEstimate(
    const GeometryContext& gctx, const Vector3D& gpos, const Vector3D& gdir,
    const BoundaryCheck& bcheck, double bwdTolerance, CorrFnc correct) const {
  // Call the common method for planar surfaces
  return detail::planarIntersectionEstimate(*this, gctx, gpos, gdir, bcheck,
                                            bwdTolerance, correct);
}
