// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/GenericApproachDescriptor.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Intersection.hpp"

void Acts::GenericApproachDescriptor::registerLayer(const Layer& lay) {
  // go through the surfaces
  for (auto& sf : m_surfaceCache) {
    auto mutableSf = const_cast<Surface*>(sf);
    mutableSf->associateLayer(lay);
  }
}

Acts::ObjectIntersection<Acts::Surface>
Acts::GenericApproachDescriptor::approachSurface(
    const GeometryContext& gctx, const Vector3D& gpos, const Vector3D& gdir,
    const BoundaryCheck& bcheck, double bwdTolerance, CorrFnc corrfnc) const {
  // the intersection estimates
  std::vector<ObjectIntersection<Surface>> sIntersections;
  sIntersections.reserve(m_surfaceCache.size());
  for (auto& sf : m_surfaceCache) {
    // intersect
    auto intersection = sf->intersectionEstimate(gctx, gpos, gdir, bcheck,
                                                 bwdTolerance, corrfnc);
    sIntersections.push_back(ObjectIntersection<Surface>(intersection, sf));
  }
  // sort the solution , bwdTolerance should be respected
  std::sort(sIntersections.begin(), sIntersections.end());
  // return what you have
  return (*sIntersections.begin());
}

const std::vector<const Acts::Surface*>&
Acts::GenericApproachDescriptor::containedSurfaces() const {
  return m_surfaceCache;
}

std::vector<const Acts::Surface*>&
Acts::GenericApproachDescriptor::containedSurfaces() {
  return m_surfaceCache;
}
