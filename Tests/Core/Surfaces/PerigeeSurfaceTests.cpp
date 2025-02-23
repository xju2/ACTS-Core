// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// clang-format off
#define BOOST_TEST_MODULE PerigeeSurface Tests
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/output_test_stream.hpp>
// clang-format on

#include <limits>

#include "Acts/Geometry/PlaneLayer.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Surfaces/InfiniteBounds.hpp"  //to get s_noBounds
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"  //to get s_noBounds
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Definitions.hpp"

using boost::test_tools::output_test_stream;
namespace utf = boost::unit_test;

namespace Acts {

namespace Test {

// Create a test context
GeometryContext tgContext = GeometryContext();

BOOST_AUTO_TEST_SUITE(PerigeeSurfaces)
/// Unit test for creating compliant/non-compliant PerigeeSurface object
BOOST_AUTO_TEST_CASE(PerigeeSurfaceConstruction) {
  // PerigeeSurface default constructor is deleted
  //
  /// Constructor with Vector3D
  Vector3D unitXYZ{1., 1., 1.};
  auto perigeeSurfaceObject = Surface::makeShared<PerigeeSurface>(unitXYZ);
  BOOST_CHECK_EQUAL(Surface::makeShared<PerigeeSurface>(unitXYZ)->type(),
                    Surface::Perigee);
  //
  /// Constructor with transform pointer, null or valid
  Translation3D translation{0., 1., 2.};
  auto pTransform = std::make_shared<const Transform3D>(translation);
  auto pNullTransform = std::make_shared<const Transform3D>();
  BOOST_CHECK_EQUAL(Surface::makeShared<PerigeeSurface>(pNullTransform)->type(),
                    Surface::Perigee);
  BOOST_CHECK_EQUAL(Surface::makeShared<PerigeeSurface>(pTransform)->type(),
                    Surface::Perigee);
  //
  /// Copy constructor
  auto copiedPerigeeSurface =
      Surface::makeShared<PerigeeSurface>(*perigeeSurfaceObject);
  BOOST_CHECK_EQUAL(copiedPerigeeSurface->type(), Surface::Perigee);
  BOOST_CHECK(*copiedPerigeeSurface == *perigeeSurfaceObject);
  //
  /// Copied and transformed
  auto copiedTransformedPerigeeSurface = Surface::makeShared<PerigeeSurface>(
      tgContext, *perigeeSurfaceObject, *pTransform);
  BOOST_CHECK_EQUAL(copiedTransformedPerigeeSurface->type(), Surface::Perigee);
}
//
/// Unit test for testing PerigeeSurface properties
BOOST_AUTO_TEST_CASE(PerigeeSurfaceProperties) {
  /// Test clone method
  Vector3D unitXYZ{1., 1., 1.};
  auto perigeeSurfaceObject = Surface::makeShared<PerigeeSurface>(unitXYZ);
  auto pClonedPerigeeSurface =
      perigeeSurfaceObject->clone(tgContext, Transform3D::Identity());
  BOOST_CHECK_EQUAL(pClonedPerigeeSurface->type(), Surface::Perigee);
  //
  /// Test type (redundant)
  BOOST_CHECK_EQUAL(perigeeSurfaceObject->type(), Surface::Perigee);
  //
  /// Test name
  BOOST_CHECK_EQUAL(perigeeSurfaceObject->name(),
                    std::string("Acts::PerigeeSurface"));
  //
  /// Test dump
  boost::test_tools::output_test_stream dumpOuput;
  perigeeSurfaceObject->toStream(tgContext, dumpOuput);
  BOOST_CHECK(
      dumpOuput.is_equal("Acts::PerigeeSurface:\n\
     Center position  (x, y, z) = (1.0000000, 1.0000000, 1.0000000)"));
}

BOOST_AUTO_TEST_CASE(EqualityOperators) {
  Vector3D unitXYZ{1., 1., 1.};
  Vector3D invalidPosition{0.0, 0.0, 0.0};
  auto perigeeSurfaceObject = Surface::makeShared<PerigeeSurface>(unitXYZ);
  auto perigeeSurfaceObject2 = Surface::makeShared<PerigeeSurface>(unitXYZ);
  auto assignedPerigeeSurface =
      Surface::makeShared<PerigeeSurface>(invalidPosition);
  /// Test equality operator
  BOOST_CHECK(*perigeeSurfaceObject == *perigeeSurfaceObject2);
  /// Test assignment
  *assignedPerigeeSurface = *perigeeSurfaceObject;
  /// Test equality of assigned to original
  BOOST_CHECK(*assignedPerigeeSurface == *perigeeSurfaceObject);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Test

}  // namespace Acts
