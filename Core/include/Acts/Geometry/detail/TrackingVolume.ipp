// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// TrackingVolume.ipp, Acts project
///////////////////////////////////////////////////////////////////

inline const Acts::Layer* TrackingVolume::associatedLayer(
    const GeometryContext& /*gctx*/, const Vector3D& gp) const {
  // confined static layers - highest hierarchy
  if (m_confinedLayers) {
    return (m_confinedLayers->object(gp).get());
  }
  // return the null pointer
  return nullptr;
}

template <typename parameters_t, typename options_t, typename corrector_t>
std::vector<LayerIntersection> TrackingVolume::compatibleLayers(
    const GeometryContext& gctx, const parameters_t& parameters,
    const options_t& options, const corrector_t& corrfnc) const {
  // the layer intersections which are valid
  std::vector<LayerIntersection> lIntersections;

  // the confinedLayers
  if (m_confinedLayers) {
    // start layer given or not - test layer
    const Layer* tLayer = options.startObject
                              ? options.startObject
                              : associatedLayer(gctx, parameters.position());
    while (tLayer != nullptr) {
      // check if the layer needs resolving
      // - resolveSensitive -> always take layer if it has a surface array
      // - resolveMaterial -> always take layer if it has material
      // - resolvePassive -> always take, unless it's a navigation layer
      // skip the start object
      if (tLayer != options.startObject && tLayer->resolve(options)) {
        // if it's a resolveable start layer, you are by definition on it
        // layer on approach intersection
        auto atIntersection =
            tLayer->surfaceOnApproach(gctx, parameters, options, corrfnc);
        // Intersection is ok - take it (move to surface on appraoch)
        if (atIntersection &&
            (atIntersection.object != options.targetSurface)) {
          // create a layer intersection
          lIntersections.push_back(LayerIntersection(
              atIntersection.intersection, tLayer, atIntersection.object));
        }
      }
      // move to next one or break because you reached the end layer
      tLayer = (tLayer == options.endObject)
                   ? nullptr
                   : tLayer->nextLayer(gctx, parameters, options);
    }
    // sort them accordingly to the navigation direction
    if (options.navDir == forward) {
      std::sort(lIntersections.begin(), lIntersections.end());
    } else {
      std::sort(lIntersections.begin(), lIntersections.end(), std::greater<>());
    }
  }
  // and return
  return lIntersections;
}

// Returns the boundary surfaces ordered in probability to hit them based on
template <typename parameters_t, typename options_t, typename corrector_t>
std::vector<BoundaryIntersection> TrackingVolume::compatibleBoundaries(
    const GeometryContext& gctx, const parameters_t& parameters,
    const options_t& options, const corrector_t& corrfnc) const {
  // Loop over boundarySurfaces and calculate the intersection
  auto& bSurfaces = boundarySurfaces();
  std::vector<BoundaryIntersection> cIntersections;
  cIntersections.reserve(bSurfaces.size());

  // Collect the boundary surfaces both directions
  for (const auto& boundary : bSurfaces) {
    // get the surface representation
    const Surface& surface = boundary->surfaceRepresentation();
    // don't intersect the exlude surface
    if (&surface != options.startObject) {
      // intersect the surface
      SurfaceIntersection sIntersection = surface.surfaceIntersectionEstimate(
          gctx, parameters, options, corrfnc);
      // surface intersection can be: reachable || on surface || overstepped
      cIntersections.push_back(BoundaryIntersection(sIntersection.intersection,
                                                    boundary.get(), &surface));
    }
  }
  // sort them accordingly to the navigation direction
  if (options.navDir == forward) {
    std::sort(cIntersections.begin(), cIntersections.end());
  } else {
    std::sort(cIntersections.begin(), cIntersections.end(), std::greater<>());
  }
  // return the compatible boundary intersections
  return cIntersections;
}

template <typename parameters_t, typename options_t, typename corrector_t>
std::vector<SurfaceIntersection>
TrackingVolume::compatibleSurfacesFromHierarchy(
    const GeometryContext& gctx, const parameters_t& parameters, double angle,
    const options_t& options, const corrector_t& corrfnc) const {
  std::vector<SurfaceIntersection> sIntersections;
  sIntersections.reserve(20);  // arbitrary

  if (m_bvhTop == nullptr || !options.navDir) {
    return sIntersections;
  }

  Vector3D dir = parameters.direction();
  if (options.navDir == backward) {
    dir *= -1;
  }

  std::vector<const Volume*> hits;
  if (angle == 0) {
    // use ray
    Ray3D obj(parameters.position(), dir);
    hits = intersectSearchHierarchy(std::move(obj), m_bvhTop);
  } else {
    Acts::Frustum<double, 3, 4> obj(parameters.position(), dir, angle);
    hits = intersectSearchHierarchy(std::move(obj), m_bvhTop);
  }

  // have cells, decompose to surfaces
  for (const Volume* vol : hits) {
    const AbstractVolume* avol = dynamic_cast<const AbstractVolume*>(vol);
    const std::vector<std::shared_ptr<const BoundarySurfaceT<AbstractVolume>>>&
        boundarySurfaces = avol->boundarySurfaces();
    for (const auto& bs : boundarySurfaces) {
      const Surface& srf = bs->surfaceRepresentation();
      SurfaceIntersection sfi =
          srf.surfaceIntersectionEstimate(gctx, parameters, options, corrfnc);
      if (sfi) {
        sIntersections.push_back(std::move(sfi));
      }
    }
  }

  // sort according to the path length
  if (options.navDir == forward) {
    std::sort(sIntersections.begin(), sIntersections.end());
  } else {
    std::sort(sIntersections.begin(), sIntersections.end(), std::greater<>());
  }

  return sIntersections;
}

template <typename T>
std::vector<const Volume*> TrackingVolume::intersectSearchHierarchy(
    const T obj, const Volume::BoundingBox* lnode) {
  std::vector<const Volume*> hits;
  hits.reserve(20);  // arbitrary
  do {
    if (lnode->intersect(obj)) {
      if (lnode->hasEntity()) {
        // found primitive
        // check obb to limit false positivies
        const Volume* vol = lnode->entity();
        const auto& obb = vol->orientedBoundingBox();
        if (obb.intersect(obj.transformed(vol->itransform()))) {
          hits.push_back(vol);
        }
        // we skip in any case, whether we actually hit the OBB or not
        lnode = lnode->getSkip();
      } else {
        // go over children
        lnode = lnode->getLeftChild();
      }
    } else {
      lnode = lnode->getSkip();
    }
  } while (lnode != nullptr);

  return hits;
}
