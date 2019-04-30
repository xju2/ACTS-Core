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
    Vector3D pos(1.34, 2.34, 3.45);
    Vector3D mom(1000., 1000., -0.100);
    Vector3D dir(mom.normalized());
    Vector3D z_axis_global(0., 0., 1.);
    /// create curvilinear parameters without covariance +1/-1 charge
    CurvilinearParameters* curvilinear_pos         = new CurvilinearParameters(nullptr, pos, mom, 1.);
    CurvilinearParameters* curvilinear_neg         = new CurvilinearParameters(nullptr, pos, mom, -1);
    NeutralCurvilinearParameters* curvilinear_neut = new NeutralCurvilinearParameters(nullptr, pos, mom);
	MultipleCurvilinearParameters multi_curvilinear_pos(0.7,curvilinear_pos);
	MultipleCurvilinearParameters multi_curvilinear_neg(0.1,curvilinear_neg);
	MultipleNeutralCurvilinearParameters multi_curvilinear_neut(0.1,curvilinear_neut);


    /// check local coordinates
    const auto   fphi   = phi(mom);
    const auto   ftheta = theta(mom);
    const double oOp    = 1. / mom.norm();

    // check parameters
    MultiConsistencyCheck(
        multi_curvilinear_pos,0, pos, mom, +1., {{0., 0., fphi, ftheta, oOp}});
    MultiConsistencyCheck(
        multi_curvilinear_neg,0, pos, mom, -1., {{0., 0., fphi, ftheta, -oOp}});
    MultiConsistencyCheck(
        multi_curvilinear_neut,0, pos, mom, 0., {{0., 0., fphi, ftheta, oOp}});

    // check that the created surface is at the position
    CHECK_CLOSE_REL(
        multi_curvilinear_pos.referenceSurface().center(tgContext), pos, 1e-6);
    CHECK_CLOSE_REL(
        multi_curvilinear_neg.referenceSurface().center(tgContext), pos, 1e-6);
    CHECK_CLOSE_REL(
        multi_curvilinear_neut.referenceSurface().center(tgContext), pos, 1e-6);
    // check that the z-axis of the created surface is along momentum direction
    CHECK_CLOSE_REL(
        multi_curvilinear_pos.referenceSurface().normal(tgContext, pos), dir, 1e-6);
    CHECK_CLOSE_REL(
        multi_curvilinear_neg.referenceSurface().normal(tgContext, pos), dir, 1e-6);
    CHECK_CLOSE_REL(
        multi_curvilinear_neut.referenceSurface().normal(tgContext, pos), dir, 1e-6);

    // check the reference frame of curvilinear parameters
    // it is the x-y frame of the created surface
    RotationMatrix3D mFrame = RotationMatrix3D::Zero();
    Vector3D         tAxis  = dir;
    Vector3D         uAxis  = (z_axis_global.cross(tAxis)).normalized();
    Vector3D         vAxis  = tAxis.cross(uAxis);
    mFrame.col(0)           = uAxis;
    mFrame.col(1)           = vAxis;
    mFrame.col(2)           = tAxis;
    CHECK_CLOSE_OR_SMALL(mFrame, multi_curvilinear_pos.referenceFrame(tgContext), 1e-6, 1e-9);
    CHECK_CLOSE_OR_SMALL(mFrame, multi_curvilinear_neg.referenceFrame(tgContext), 1e-6, 1e-9);
    CHECK_CLOSE_OR_SMALL(mFrame, multi_curvilinear_neut.referenceFrame(tgContext), 1e-6, 1e-9);

    /// modification test with set methods
    double ux = 0.1;
    double uy = 0.5;
    multi_curvilinear_pos.set<eLOC_0>(tgContext,ux,0);
    multi_curvilinear_pos.set<eLOC_1>(tgContext,uy,0);
    // the local parameter should still be (0,0) for Curvilinear
    BOOST_CHECK_EQUAL(multi_curvilinear_pos.parameters(0)[eLOC_0], 0);
    BOOST_CHECK_EQUAL(multi_curvilinear_pos.parameters(0)[eLOC_1], 0);
    // the position should be updated though
    Vector3D uposition = multi_curvilinear_neg.referenceSurface().transform(tgContext)
        * Vector3D(ux, uy, 0.);
    // the position should be updated
    CHECK_CLOSE_REL(multi_curvilinear_pos.position(0), uposition, 1e-6);
    // it should be the position of the surface
    CHECK_CLOSE_REL(
        multi_curvilinear_pos.referenceSurface().center(tgContext), uposition, 1e-6);

	
    double uphi   = 1.2;
    double utheta = 0.2;
    double uqop   = 0.025;

    multi_curvilinear_pos.set<ePHI>(tgContext,uphi,0);
    multi_curvilinear_pos.set<eTHETA>(tgContext,utheta,0);
    multi_curvilinear_pos.set<eQOP>(tgContext,uqop,0);
    // we should have a new updated momentum
    Vector3D umomentum = 40. * Vector3D(cos(uphi) * sin(utheta),
                                        sin(uphi) * sin(utheta),
                                        cos(utheta));
    CHECK_CLOSE_REL(umomentum, multi_curvilinear_pos.momentum(0), 1e-6);
    // the updated momentum should be the col(2) of the transform
    CHECK_CLOSE_REL(
        umomentum.normalized(),
        multi_curvilinear_pos.referenceSurface().transform(tgContext).rotation().col(2),
        1e-6);

	/// append a component 
    Vector3D pos_append_1(1.34, 2.34, 3.45);
    Vector3D mom_append_1(1000., 1000., -0.100);
    Vector3D dir_append_1(mom.normalized());
    Vector3D pos_append_2(1.34, 2.34, 3.45);
    Vector3D mom_append_2(1000., 1000., -0.100);
    Vector3D dir_append_2(mom.normalized());
    /// create curvilinear parameters without covariance +1/-1 charge
    CurvilinearParameters* curvilinear_pos_append_1  = new CurvilinearParameters(nullptr, pos_append_1, mom_append_1, 1.);
    CurvilinearParameters* curvilinear_pos_append_2  = new CurvilinearParameters(nullptr, pos_append_2, mom_append_2, 1.);
	multi_curvilinear_pos.append(0.1,curvilinear_pos_append_1);
	multi_curvilinear_pos.append(0.2,curvilinear_pos_append_2);
  }


}  // namespace Test
}  // namespace Acts
