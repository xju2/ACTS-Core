// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// clang-format off
#define BOOST_TEST_MODULE MultiCurvilinearParameters Tests
#define BOOST_TEST_DYN_LINK
#include <boost/test/included/unit_test.hpp>
// clang-format on

#include "Acts/EventData/NeutralParameters.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Definitions.hpp"

#include "ParametersTestHelper.hpp"

namespace Acts {
using VectorHelpers::phi;
using VectorHelpers::theta;
namespace Test {

  /// @brief Unit test for Curvilinear parameters
  ///
  BOOST_AUTO_TEST_CASE(multi_curvilinear_initialization)
  {
    // Create a test context
    GeometryContext tgContext = GeometryContext();

    // some position and momentum
    Vector3D pos0(1., 2., 3.);
    Vector3D pos1(2.01, 2.01, 3.01);
    Vector3D pos2(3.02, 2.02, 3.02);
    Vector3D mom0(1000., 1000., -0.100);
    Vector3D mom1(1000.01, 1000., -0.100);
    Vector3D mom2(1000.02, 1000., -0.100);
    Vector3D dir0(mom0.normalized());
    Vector3D dir1(mom1.normalized());
    Vector3D dir2(mom2.normalized());
    Vector3D z_axis_global(0., 0., 1.);
    /// create curvilinear parameters without covariance +1/-1 charge
    CurvilinearParameters* curvilinear_pos_0
        = new CurvilinearParameters(nullptr, pos0, mom0, 1.);
    CurvilinearParameters* curvilinear_pos_1
        = new CurvilinearParameters(nullptr, pos1, mom1, 1.);
    CurvilinearParameters* curvilinear_pos_2
        = new CurvilinearParameters(nullptr, pos2, mom2, 1.);

    // MultipleCurvilinearParameters multi_curvilinear_pos(0.4,
    // curvilinear_pos_0);
    MultipleCurvilinearParameters multi_curvilinear_pos;
    multi_curvilinear_pos.append(0.4, curvilinear_pos_0);
    multi_curvilinear_pos.append(0.7, curvilinear_pos_1);
    multi_curvilinear_pos.append(0.5, curvilinear_pos_2);
    multi_curvilinear_pos.makeReferenceSurface();
    BOOST_CHECK_EQUAL(multi_curvilinear_pos.size(), 3);

    // test sort method
    MultipleCurvilinearParameters::TrackParMapConstIter it
        = multi_curvilinear_pos.getTrackList().begin();
    BOOST_CHECK_EQUAL((*it).first, 0.7);
    ++it;
    BOOST_CHECK_EQUAL((*it).first, 0.5);
    ++it;
    BOOST_CHECK_EQUAL((*it).first, 0.4);
  }

}  // namespace Test
}  // namespace Acts
