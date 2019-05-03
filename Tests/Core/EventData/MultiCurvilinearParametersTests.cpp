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
    Vector3D pos_combine = 0.4 * pos0 + 0.7 * pos1 + 0.5 * pos2;
    Vector3D dir_combine = (0.4 * mom0 + 0.7 * mom1 + 0.5 * mom2).normalized();
    Vector3D z_axis_global(0., 0., 1.);
    /// create curvilinear parameters without covariance +1/-1 charge
    CurvilinearParameters* curvilinear_pos_0
        = new CurvilinearParameters(nullptr, pos0, mom0, 1.);
    CurvilinearParameters* curvilinear_pos_1
        = new CurvilinearParameters(nullptr, pos1, mom1, 1.);
    CurvilinearParameters* curvilinear_pos_2
        = new CurvilinearParameters(nullptr, pos2, mom2, 1.);
    CurvilinearParameters* curvilinear_neg
        = new CurvilinearParameters(nullptr, pos0, mom0, -1.);

    MultipleCurvilinearParameters multi_curvilinear_pos(0.4, curvilinear_pos_0);
    // MultipleCurvilinearParameters multi_curvilinear_pos;
    // multi_curvilinear_pos.append(0.4,curvilinear_pos_0);
    multi_curvilinear_pos.append(0.7, curvilinear_pos_1);
    multi_curvilinear_pos.append(0.5, curvilinear_pos_2);
    multi_curvilinear_pos.makeReferenceSurface();

    // test sort method
    MultipleCurvilinearParameters::TrackParMapConstIter it
        = multi_curvilinear_pos.getTrackList().begin();
    BOOST_CHECK_EQUAL((*it).first.first, 0.7);
    ++it;
    BOOST_CHECK_EQUAL((*it).first.first, 0.5);
    ++it;
    BOOST_CHECK_EQUAL((*it).first.first, 0.4);

    // test position/momentum(id) method
    BOOST_CHECK_EQUAL(multi_curvilinear_pos.position(0), pos0);
    BOOST_CHECK_EQUAL(multi_curvilinear_pos.position(1), pos1);
    BOOST_CHECK_EQUAL(multi_curvilinear_pos.position(2), pos2);
    BOOST_CHECK_EQUAL(multi_curvilinear_pos.momentum(0), mom0);
    BOOST_CHECK_EQUAL(multi_curvilinear_pos.momentum(1), mom1);
    BOOST_CHECK_EQUAL(multi_curvilinear_pos.momentum(2), mom2);

    /// check local coordinates
    const auto   fphi   = phi(mom0);
    const auto   ftheta = theta(mom0);
    const double oOp    = 1. / mom0.norm();
    MultiConsistencyCheck(multi_curvilinear_pos,
                          0,
                          pos0,
                          mom0,
                          +1.,
                          {{0., 0., fphi, ftheta, oOp}});

    // check the created reference surface is the combination
    CHECK_CLOSE_REL(multi_curvilinear_pos.referenceSurface().center(tgContext),
                    pos_combine,
                    1e-6);
    CHECK_CLOSE_REL(multi_curvilinear_pos.referenceSurface().normal(tgContext),
                    dir_combine,
                    1e-6);

    // check the reference frame of curvilinear parameters
    // it is the x-y frame of the created surface
    RotationMatrix3D mFrame = RotationMatrix3D::Zero();
    Vector3D         tAxis  = dir_combine;
    Vector3D         uAxis  = (z_axis_global.cross(tAxis)).normalized();
    Vector3D         vAxis  = tAxis.cross(uAxis);
    mFrame.col(0)           = uAxis;
    mFrame.col(1)           = vAxis;
    mFrame.col(2)           = tAxis;
    CHECK_CLOSE_OR_SMALL(
        mFrame, multi_curvilinear_pos.referenceFrame(tgContext), 1e-6, 1e-9);

    /// modification test with set methods
    double ux = 0.1;
    double uy = 0.5;
    multi_curvilinear_pos.set<eLOC_0>(tgContext, ux, 0);
    multi_curvilinear_pos.set<eLOC_1>(tgContext, uy, 0);
    // the local parameter should still be (0,0) for Curvilinear
    BOOST_CHECK_EQUAL(multi_curvilinear_pos.parameters(0)[eLOC_0], 0);
    BOOST_CHECK_EQUAL(multi_curvilinear_pos.parameters(0)[eLOC_1], 0);
    // the position should be updated though
    Vector3D uposition
        = curvilinear_neg->referenceSurface().transform(tgContext)
        * Vector3D(ux, uy, 0.);  // the position should be updated
    // the first (0st) curvilinear should be updated
    CHECK_CLOSE_REL(multi_curvilinear_pos.position(0), uposition, 1e-6);
    // it should be the position of the surface
    CHECK_CLOSE_REL(multi_curvilinear_pos.referenceSurface(0).center(tgContext),
                    uposition,
                    1e-6);

    /// modification test with set methods
    double uphi   = 1.2;
    double utheta = 0.2;
    double uqop   = 0.025;
    multi_curvilinear_pos.set<ePHI>(tgContext, uphi, 0);
    multi_curvilinear_pos.set<eTHETA>(tgContext, utheta, 0);
    multi_curvilinear_pos.set<eQOP>(tgContext, uqop, 0);
    // we should have a new updated momentum
    Vector3D umomentum = 40. * Vector3D(cos(uphi) * sin(utheta),
                                        sin(uphi) * sin(utheta),
                                        cos(utheta));
    CHECK_CLOSE_REL(umomentum, multi_curvilinear_pos.momentum(0), 1e-6);
    // the updated momentum should be the col(2) of the transform
    CHECK_CLOSE_REL(umomentum.normalized(),
                    multi_curvilinear_pos.referenceSurface(0)
                        .transform(tgContext)
                        .rotation()
                        .col(2),
                    1e-6);

    // set<> method should also update the combined reference surface
    // this is to test if set<> correctly update the combined reference surface
    CurvilinearParameters* curvilinear_pos_update
        = new CurvilinearParameters(nullptr, pos0, mom0, 1.);
    curvilinear_pos_update->set<eLOC_0>(tgContext, ux);
    curvilinear_pos_update->set<eLOC_1>(tgContext, uy);
    curvilinear_pos_update->set<ePHI>(tgContext, uphi);
    curvilinear_pos_update->set<eTHETA>(tgContext, utheta);
    curvilinear_pos_update->set<eQOP>(tgContext, uqop);

    MultipleCurvilinearParameters multi_curvilinear_pos_update(
        0.4, curvilinear_pos_update);
    CurvilinearParameters* curvilinear_pos_update_1
        = new CurvilinearParameters(nullptr, pos1, mom1, 1.);
    CurvilinearParameters* curvilinear_pos_update_2
        = new CurvilinearParameters(nullptr, pos2, mom2, 1.);
    multi_curvilinear_pos_update.append(0.7, curvilinear_pos_update_1);
    multi_curvilinear_pos_update.append(0.5, curvilinear_pos_update_2);
    multi_curvilinear_pos_update.makeReferenceSurface();

    CHECK_CLOSE_REL(
        multi_curvilinear_pos_update.referenceSurface().center(tgContext),
        multi_curvilinear_pos.referenceSurface().center(tgContext),
        1e-6);
  }

}  // namespace Test
}  // namespace Acts
