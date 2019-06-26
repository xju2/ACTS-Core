// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// clang-format off
#define BOOST_TEST_MODULE BoundParameters Tests
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
// clang-format on

#include "Acts/EventData/NeutralParameters.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/StrawSurface.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Units.hpp"

#include "ParametersTestHelper.hpp"

namespace bdata = boost::unit_test::data;
namespace tt    = boost::test_tools;

namespace Acts {

namespace Test {

  // Create a test context
  GeometryContext tgContext = GeometryContext();

  /// @brief Unit test for parameters at a plane
  ///
  BOOST_DATA_TEST_CASE(
      bound_to_plane_test,
      bdata::random((bdata::seed = 1240,
                     bdata::distribution
                     = std::uniform_real_distribution<>(-1000., 1000.)))
          ^ bdata::random((bdata::seed = 2351,
                           bdata::distribution
                           = std::uniform_real_distribution<>(-1000., 1000.)))
          ^ bdata::random((bdata::seed = 3412,
                           bdata::distribution
                           = std::uniform_real_distribution<>(-1000., 1000.)))
          ^ bdata::random((bdata::seed = 5732,
                           bdata::distribution
                           = std::uniform_real_distribution<>(0., M_PI)))
          ^ bdata::random((bdata::seed = 8941,
                           bdata::distribution
                           = std::uniform_real_distribution<>(0., M_PI)))
          ^ bdata::random((bdata::seed = 1295,
                           bdata::distribution
                           = std::uniform_real_distribution<>(0., M_PI)))
          ^ bdata::xrange(1),
      x,
      y,
      z,
      a,
      b,
      c,
      index)
  {
    (void)index;
    Vector3D center{x, y, z};
    auto     transform = std::make_shared<Transform3D>();
    transform->setIdentity();
    RotationMatrix3D rot;
    rot = AngleAxis3D(a, Vector3D::UnitX()) * AngleAxis3D(b, Vector3D::UnitY())
        * AngleAxis3D(c, Vector3D::UnitZ());
    transform->prerotate(rot);
    transform->pretranslate(center);
    // create the surfacex
    auto bounds   = std::make_shared<RectangleBounds>(100., 100.);
    auto pSurface = Surface::makeShared<PlaneSurface>(
        transform, bounds);  // surface use +1

    // now create parameters on this surface
    // l_x, l_y, phi, theta, q/p (1/p)
    std::array<double, 5> pars_array = {{-0.1234, 9.8765, 0.45, 0.888, 0.001}};
    TrackParametersBase::ParVector_t pars;
    pars << pars_array[0], pars_array[1], pars_array[2], pars_array[3],
        pars_array[4];

    BoundParameters* ataPlane_from_pars_0
        = new BoundParameters(tgContext, nullptr, pars, pSurface);  //+2
    BoundParameters* ataPlane_from_pars_1
        = new BoundParameters(tgContext, nullptr, pars, pSurface);  //+3

    // make multi bound par
    // MultipleBoundParameters multi_ataPlane_from_pars( 0.6,
    // ataPlane_from_pars_0, pSurface);  //+4
    MultipleBoundParameters multi_ataPlane_from_pars(pSurface);  //+4
    multi_ataPlane_from_pars.append(0.6, ataPlane_from_pars_0);
    multi_ataPlane_from_pars.append(0.7, ataPlane_from_pars_1);

    // test the append method
    MultipleBoundParameters::TrackParMapConstIter it
        = multi_ataPlane_from_pars.getTrackList().begin();
    BOOST_CHECK_EQUAL((*it).first, 0.7);
    ++it;
    BOOST_CHECK_EQUAL((*it).first, 0.6);

    // check shared ownership of same surface
    BOOST_CHECK_EQUAL(&multi_ataPlane_from_pars.referenceSurface(),
                      pSurface.get());
    BOOST_CHECK_EQUAL(pSurface.use_count(), 4);

    // check that the reference frame is the rotation matrix
    CHECK_CLOSE_REL(
        multi_ataPlane_from_pars.referenceFrame(tgContext), rot, 1e-6);
  }
}
}
