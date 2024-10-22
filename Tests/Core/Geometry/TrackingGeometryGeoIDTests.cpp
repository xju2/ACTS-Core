// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///  Boost include(s)
#define BOOST_TEST_MODULE GeometryID Tests
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Utilities/Units.hpp"
#include "TrackingVolumeCreation.hpp"

using namespace Acts::UnitLiterals;

namespace Acts {
namespace Test {

// Create a test context
GeometryContext tgContext = GeometryContext();

///  create three cylinder surfaces
///  the surface radius (will also be the layer radius)
double iVsurfaceHalfLengthZ = 50_mm;
double iVsurfaceRadius = 25_mm;
double iVsurfaceRstagger = 5_mm;
double iVsurfaceZoverlap = 10_mm;
double iVlayerEnvelope = 0.5_mm;
double iVvolumeEnvelope = 10_mm;
double iVvolumeRadius = iVsurfaceRadius + 0.5 * iVsurfaceRstagger +
                        iVlayerEnvelope + iVvolumeEnvelope;

///  the surface radius (will also be the layer radius)
double oVsurfaceHalfLengthZ = 50_mm;
double oVsurfaceRadius = 100_mm;
double oVsurfaceRstagger = 5_mm;
double oVsurfaceZoverlap = 10_mm;
double oVlayerEnvelope = 0.5_mm;
double oVvolumeEnvelope = 10_mm;
double oVvolumeRadius = oVsurfaceRadius + 0.5 * oVsurfaceRstagger +
                        oVlayerEnvelope + oVvolumeEnvelope;

///  inner volume
auto iVolume = constructCylinderVolume(
    tgContext, iVsurfaceHalfLengthZ, iVsurfaceRadius, iVsurfaceRstagger,
    iVsurfaceZoverlap, iVlayerEnvelope, iVvolumeEnvelope, 0., iVvolumeRadius,
    "InnerVolume");

BOOST_AUTO_TEST_CASE(GeometryID_innervolume_test) {
  BOOST_CHECK_EQUAL(0ul, iVolume->geoID().value());
  // check the boundary surfaces
  for (auto bSf : iVolume->boundarySurfaces()) {
    BOOST_CHECK_EQUAL(0ul, bSf->surfaceRepresentation().geoID().value());
    for (auto lay : iVolume->confinedLayers()->arrayObjects()) {
      BOOST_CHECK_EQUAL(0ul, lay->geoID().value());
      // check the approach surfaces
      for (auto asf : lay->approachDescriptor()->containedSurfaces()) {
        BOOST_CHECK_EQUAL(0ul, asf->geoID().value());
      }
      // check the layer surface array
      for (auto ssf : lay->surfaceArray()->surfaces()) {
        BOOST_CHECK_EQUAL(0ul, ssf->geoID().value());
      }
    }
  }
}

///  outer volume
auto oVolume = constructCylinderVolume(
    tgContext, oVsurfaceHalfLengthZ, oVsurfaceRadius, oVsurfaceRstagger,
    oVsurfaceZoverlap, oVlayerEnvelope, oVvolumeEnvelope, iVvolumeRadius,
    oVvolumeRadius, "OuterVolume");

BOOST_AUTO_TEST_CASE(GeometryID_outervolume_test) {
  BOOST_CHECK_EQUAL(0ul, oVolume->geoID().value());
  // check the boundary surfaces
  for (auto bSf : iVolume->boundarySurfaces()) {
    BOOST_CHECK_EQUAL(0ul, bSf->surfaceRepresentation().geoID().value());
    for (auto lay : oVolume->confinedLayers()->arrayObjects()) {
      BOOST_CHECK_EQUAL(0ul, lay->geoID().value());
      // check the approach surfaces
      for (auto asf : lay->approachDescriptor()->containedSurfaces()) {
        BOOST_CHECK_EQUAL(0ul, asf->geoID().value());
      }
      // check the layer surface array
      for (auto ssf : lay->surfaceArray()->surfaces()) {
        BOOST_CHECK_EQUAL(0ul, ssf->geoID().value());
      }
    }
  }
}
//
double oVvolumeHalfZ =
    (4 * oVsurfaceHalfLengthZ - oVsurfaceZoverlap) + oVvolumeEnvelope;
// now create the container volume
auto hVolume = constructContainerVolume(
    tgContext, iVolume, oVolume, oVvolumeRadius, oVvolumeHalfZ, "Container");

///  pre-check on GeometryID
BOOST_AUTO_TEST_CASE(GeometryID_containervolume_test) {
  ///  let's check that the geometry ID values are all 0
  BOOST_CHECK_EQUAL(0ul, hVolume->geoID().value());
  /// check the boundaries of the hVolume, should also be 0
  for (auto hbsf : hVolume->boundarySurfaces()) {
    BOOST_CHECK_EQUAL(0ul, hbsf->surfaceRepresentation().geoID().value());
  }
  for (auto cVol : hVolume->confinedVolumes()->arrayObjects()) {
    /// let's check everything is set to 0
    BOOST_CHECK_EQUAL(0ul, cVol->geoID().value());
    // check the boundary surfaces
    for (auto bSf : cVol->boundarySurfaces()) {
      BOOST_CHECK_EQUAL(0ul, bSf->surfaceRepresentation().geoID().value());
    }
    for (auto lay : cVol->confinedLayers()->arrayObjects()) {
      BOOST_CHECK_EQUAL(0ul, lay->geoID().value());
      // check the approach surfaces
      for (auto asf : lay->approachDescriptor()->containedSurfaces()) {
        BOOST_CHECK_EQUAL(0ul, asf->geoID().value());
      }
      // check the layer surface array
      for (auto ssf : lay->surfaceArray()->surfaces()) {
        BOOST_CHECK_EQUAL(0ul, ssf->geoID().value());
      }
    }
  }
}

}  //  end of namespace Test
}  //  end of namespace Acts
