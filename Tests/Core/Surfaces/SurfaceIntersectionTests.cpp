// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#define BOOST_TEST_MODULE Surface Bounds Tests

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
// leave blank line

#include <boost/test/data/test_case.hpp>
// leave blank line

#include <boost/test/output_test_stream.hpp>
// leave blank line

#include "Acts/Surfaces/ConeSurface.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/StrawSurface.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Intersection.hpp"

namespace Acts {
namespace Test {

// Create a test context
GeometryContext tgContext = GeometryContext();

// Test for planar surfaces
template <typename planar_surface_t>
void testPlanarSurface(const planar_surface_t& bound,
                       const planar_surface_t& unbound, double zDist,
                       double xDelta) {
  // Simple tests
  Vector3D start{0., 0., 0.};
  Vector3D onSurface{0., 0., zDist};
  Vector3D onSurfaceDelta{0., 0., zDist + 0.5 * s_onSurfaceTolerance};
  Vector3D oversteppedWithin{0., 0., zDist + 1. * units::_mm};
  Vector3D oversteppedToomuch{0., 0., zDist + 10. * units::_mm};
  double overstepTolerance = 3. * units::_mm;

  Vector3D alongZ{0., 0., 1.};
  Vector3D oppositeZ{0., 0., -1.};
  Vector3D missing = Vector3D(xDelta, 0., zDist).normalized();

  // check reachable mode
  auto intersection =
      bound.intersectionEstimate(tgContext, start, alongZ, true);

  bool isreachable = (intersection.status == IntersectionStatus::reachable);
  BOOST_TEST(isreachable);
  CHECK_CLOSE_ABS(intersection.pathLength, zDist, s_onSurfaceTolerance);

  // check unreachable mode
  intersection = bound.intersectionEstimate(tgContext, start, oppositeZ, true);
  bool isnotreachable =
      (intersection.status == IntersectionStatus::unreachable);
  BOOST_TEST(isnotreachable);
  CHECK_CLOSE_ABS(intersection.pathLength, -zDist, s_onSurfaceTolerance);

  // check onSurface mode
  intersection = bound.intersectionEstimate(tgContext, onSurface, alongZ, true);
  bool isonsurface = (intersection.status == IntersectionStatus::onSurface);
  BOOST_TEST(isonsurface);
  CHECK_CLOSE_ABS(intersection.pathLength, 0., s_onSurfaceTolerance);

  // check onSurface mode + delta
  intersection =
      bound.intersectionEstimate(tgContext, onSurfaceDelta, alongZ, true);
  isonsurface = (intersection.status == IntersectionStatus::onSurface);
  BOOST_TEST(isonsurface);
  CHECK_CLOSE_ABS(intersection.pathLength, 0., s_onSurfaceTolerance);

  // check for overstepped mode, within tolerqance
  intersection = bound.intersectionEstimate(tgContext, oversteppedWithin,
                                            alongZ, true, overstepTolerance);
  bool isoverstepped = (intersection.status == IntersectionStatus::overstepped);
  BOOST_TEST(isoverstepped);
  CHECK_CLOSE_ABS(intersection.pathLength, 0., overstepTolerance);

  // check for overstepped mode, too much
  intersection = bound.intersectionEstimate(tgContext, oversteppedToomuch,
                                            alongZ, true, overstepTolerance);
  isnotreachable = (intersection.status == IntersectionStatus::unreachable);
  BOOST_TEST(isnotreachable);

  // check boundary testing: unreachable
  double inclinedPath = std::sqrt(xDelta * xDelta + zDist * zDist);

  intersection = bound.intersectionEstimate(tgContext, start, missing, true);

  isreachable = (intersection.status == IntersectionStatus::reachable);
  BOOST_TEST(!isreachable);

  CHECK_CLOSE_ABS(intersection.pathLength, inclinedPath, s_onSurfaceTolerance);

  // check boundary testing: reachable again
  intersection = bound.intersectionEstimate(tgContext, start, missing, false);

  isreachable = (intersection.status == IntersectionStatus::reachable);
  BOOST_TEST(isreachable);
  CHECK_CLOSE_ABS(intersection.pathLength, inclinedPath, s_onSurfaceTolerance);

  // check boundary testing: unbound surface is always reachable - independent
  // of boundary check
  intersection = unbound.intersectionEstimate(tgContext, start, missing, true);

  isreachable = (intersection.status == IntersectionStatus::reachable);
  BOOST_TEST(isreachable);
  CHECK_CLOSE_ABS(intersection.pathLength, inclinedPath, s_onSurfaceTolerance);

  // second check
  intersection = unbound.intersectionEstimate(tgContext, start, missing, false);

  isreachable = (intersection.status == IntersectionStatus::reachable);
  BOOST_TEST(isreachable);
  CHECK_CLOSE_ABS(intersection.pathLength, inclinedPath, s_onSurfaceTolerance);
}

/// Open the test suite
///
/// @brief these tests check the validity of the Intersection and
/// and IntersectionStatus handling
BOOST_AUTO_TEST_SUITE(Surfaces)

/// Unit test for intersecting plane surfaces
BOOST_AUTO_TEST_CASE(PlaneSurfaceIntersection) {
  // place a plane along z
  double zPos = 10. * units::_cm;
  Translation3D translation{0., 0., zPos};
  auto pTransform = std::make_shared<const Transform3D>(translation);
  auto rBounds = std::make_shared<const RectangleBounds>(10. * units::_cm,
                                                         10. * units::_cm);
  auto boundPlane = Surface::makeShared<PlaneSurface>(pTransform, rBounds);
  auto unboundPlane = Surface::makeShared<PlaneSurface>(pTransform, nullptr);

  // Test the plane surface
  testPlanarSurface<PlaneSurface>(*boundPlane, *unboundPlane, zPos,
                                  zPos + 1. * units::_cm);
}

/// Unit test for intersecting plane surfaces
BOOST_AUTO_TEST_CASE(DiscSurfaceIntersection) {
  // place a plane along z
  double zPos = 10. * units::_cm;
  Translation3D translation{0., 0., zPos};
  auto pTransform = std::make_shared<const Transform3D>(translation);
  auto rBounds = std::make_shared<const RadialBounds>(0., 10. * units::_cm);
  auto boundDisc = Surface::makeShared<DiscSurface>(pTransform, rBounds);
  auto unboundDisc = Surface::makeShared<DiscSurface>(pTransform, nullptr);

  // Test the plane surface
  testPlanarSurface<DiscSurface>(*boundDisc, *unboundDisc, zPos,
                                 zPos + 1. * units::_cm);
}

/// Unit test for intersecting perigee surfaces
BOOST_AUTO_TEST_CASE(PerigeeSurfaceTest) {
  // place a plane along z
  Vector3D ip{0., 0., 0.};
  auto perigee = Surface::makeShared<PerigeeSurface>(ip);

  Vector3D start{20. * units::_mm, 1. * units::_mm, 0.};
  Vector3D reach{-1., 0., 0.};
  Vector3D unreach{1., 0., 0.};

  // check reachable mode
  auto intersection =
      perigee->intersectionEstimate(tgContext, start, reach, true);

  bool isreachable = (intersection.status == IntersectionStatus::reachable);
  BOOST_TEST(isreachable);
  CHECK_CLOSE_ABS(intersection.pathLength, 20. * units::_mm,
                  s_onSurfaceTolerance);

  // check unreachable mode
  intersection = perigee->intersectionEstimate(tgContext, start, unreach, true);

  bool isnotreachable =
      (intersection.status == IntersectionStatus::unreachable);
  BOOST_TEST(isnotreachable);
  CHECK_CLOSE_ABS(intersection.pathLength, -20. * units::_mm,
                  s_onSurfaceTolerance);
}

/// Unit test for intersecting cylinder surfaces
BOOST_AUTO_TEST_CASE(CylinderSurfaceTest) {
  // place a plane along z
  double R = 10. * units::_cm;
  double hZ = 100. * units::_cm;
  double alpha = 0.1;
  Translation3D translation{0., 0., 0.};
  auto cTransform = std::make_shared<const Transform3D>(translation);
  auto cBounds = std::make_shared<const CylinderBounds>(R, hZ);
  auto cylinder = Surface::makeShared<CylinderSurface>(cTransform, cBounds);

  /// the start position inside
  Vector3D inside(R * cos(alpha), 0., 0.);
  Vector3D outside(cos(alpha), R, 0.);
  Vector3D onSurfaceAlong(R * cos(alpha), R * sin(alpha), 0.);
  Vector3D onSurfaceOpposite(R * cos(alpha), -R * sin(alpha), 0.);
  Vector3D oversteppedWithin(R * cos(alpha), R * sin(alpha) + 3 * units::_mm,
                             0.);
  Vector3D oversteppedToomuch(R * cos(alpha), R * sin(alpha) + 7 * units::_mm,
                              0.);
  double overstepTolerance = 5. * units::_mm;

  Vector3D direction(0., 1., 0.);
  Vector3D outsideDirection = Vector3D(0., 2., 1000.).normalized();

  double deltaY = R * sin(alpha);

  // check reachable mode
  auto intersection =
      cylinder->intersectionEstimate(tgContext, inside, direction, true);
  bool isreachable = (intersection.status == IntersectionStatus::reachable);
  BOOST_TEST(isreachable);
  bool isalong = onSurfaceAlong.isApprox(intersection.position);
  BOOST_TEST(isalong);
  CHECK_CLOSE_ABS(intersection.pathLength, deltaY, s_onSurfaceTolerance);

  // check reachable mode backwards, will have two solutions so this is fine as
  // well
  intersection =
      cylinder->intersectionEstimate(tgContext, inside, -direction, true);
  isreachable = (intersection.status == IntersectionStatus::reachable);
  BOOST_TEST(isreachable);

  // check overstepping mode - inside tolerance
  intersection = cylinder->intersectionEstimate(
      tgContext, oversteppedWithin, direction, true, overstepTolerance);

  bool isoverstepped = (intersection.status == IntersectionStatus::overstepped);
  BOOST_TEST(isoverstepped);

  CHECK_CLOSE_ABS(intersection.pathLength, 0., overstepTolerance);

  // check overstepping mode - outside tolerance
  intersection = cylinder->intersectionEstimate(
      tgContext, oversteppedToomuch, direction, true, overstepTolerance);

  // not overstepped and not reachable
  isoverstepped = (intersection.status == IntersectionStatus::overstepped);
  BOOST_TEST(!isoverstepped);
  isreachable = (intersection.status == IntersectionStatus::reachable);
  BOOST_TEST(!isreachable);

  // check reachable mode - without bounds
  intersection = cylinder->intersectionEstimate(tgContext, inside,
                                                outsideDirection, false);
  isreachable = (intersection.status == IntersectionStatus::reachable);
  BOOST_TEST(isreachable);

  // check reachable mode - with bounds
  intersection =
      cylinder->intersectionEstimate(tgContext, inside, outsideDirection, true);
  isreachable = (intersection.status == IntersectionStatus::reachable);
  BOOST_TEST(!isreachable);

  // check on surface - on surface / opposite
  intersection = cylinder->intersectionEstimate(tgContext, onSurfaceOpposite,
                                                direction, true);
  bool isonsurface = (intersection.status == IntersectionStatus::onSurface);
  BOOST_TEST(isonsurface);
  bool isnototherside = !onSurfaceAlong.isApprox(intersection.position);
  BOOST_TEST(isnototherside);

  // check on surface - on surface / opposite
  intersection = cylinder->intersectionEstimate(tgContext, onSurfaceAlong,
                                                direction, true);
  isonsurface = (intersection.status == IntersectionStatus::onSurface);
  BOOST_TEST(isonsurface);

  // now intersect the other side
  intersection = cylinder->intersectionEstimate(tgContext, onSurfaceOpposite,
                                                direction, true, -1.);
  isreachable = (intersection.status == IntersectionStatus::reachable);
  BOOST_TEST(isreachable);
  bool isotherside = onSurfaceAlong.isApprox(intersection.position);
  BOOST_TEST(isotherside);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Test

}  // namespace Acts