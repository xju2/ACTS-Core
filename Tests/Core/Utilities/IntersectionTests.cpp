// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// Boost include(s)
#define BOOST_TEST_MODULE GeometryID Tests
#include <algorithm>
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Utilities/Intersection.hpp"

namespace Acts {
namespace Test {

/// test of the intersection class
BOOST_AUTO_TEST_CASE(IntersectionTest) {
  // a few valid intersections
  // all positively sortered
  Intersection fIp(Vector3D(0., 1., 0.), 1., IntersectionStatus::reachable);
  Intersection sIp(Vector3D(0., 2., 0.), 2., IntersectionStatus::reachable);
  Intersection tIp(Vector3D(0., 3., 0.), 3., IntersectionStatus::reachable);

  // a non-valid intersection
  Intersection nIp(Vector3D(3., 3., 0.), 3., IntersectionStatus::unreachable);
  BOOST_CHECK(!nIp);

  std::vector<Intersection> fstpIntersections = {fIp, sIp, tIp};
  std::vector<Intersection> tsfpIntersections = {tIp, sIp, fIp};

  // let's sort the tsf intersection, it should give fst
  std::sort(tsfpIntersections.begin(), tsfpIntersections.end());
  BOOST_CHECK_EQUAL(fstpIntersections[0].pathLength,
                    tsfpIntersections[0].pathLength);
  BOOST_CHECK_EQUAL(fstpIntersections[1].pathLength,
                    tsfpIntersections[1].pathLength);
  BOOST_CHECK_EQUAL(fstpIntersections[2].pathLength,
                    tsfpIntersections[2].pathLength);

  // let's sort them with greater
  std::sort(tsfpIntersections.begin(), tsfpIntersections.end(),
            std::greater<>());
  BOOST_CHECK_EQUAL(fstpIntersections[0].pathLength,
                    tsfpIntersections[2].pathLength);
  BOOST_CHECK_EQUAL(fstpIntersections[1].pathLength,
                    tsfpIntersections[1].pathLength);
  BOOST_CHECK_EQUAL(fstpIntersections[2].pathLength,
                    tsfpIntersections[0].pathLength);

  // now let's create one with a non-valid intersection, it should be shuffled
  // last
  std::vector<Intersection> ntfspIntersections = {nIp, tIp, fIp, sIp};
  std::vector<Intersection> tfnsnpIntersections = {tIp, fIp, nIp, sIp, nIp};

  // shuffle the intersections
  std::sort(ntfspIntersections.begin(), ntfspIntersections.end());
  BOOST_CHECK_EQUAL(fstpIntersections[0].pathLength,
                    ntfspIntersections[0].pathLength);
  BOOST_CHECK_EQUAL(fstpIntersections[1].pathLength,
                    ntfspIntersections[1].pathLength);
  BOOST_CHECK_EQUAL(fstpIntersections[2].pathLength,
                    ntfspIntersections[2].pathLength);
  BOOST_CHECK_EQUAL(bool(ntfspIntersections[3]), false);

  std::sort(tfnsnpIntersections.begin(), tfnsnpIntersections.end());
  BOOST_CHECK_EQUAL(fstpIntersections[0].pathLength,
                    tfnsnpIntersections[0].pathLength);
  BOOST_CHECK_EQUAL(fstpIntersections[1].pathLength,
                    tfnsnpIntersections[1].pathLength);
  BOOST_CHECK_EQUAL(fstpIntersections[2].pathLength,
                    tfnsnpIntersections[2].pathLength);
  BOOST_CHECK_EQUAL(bool(tfnsnpIntersections[3]), false);
  BOOST_CHECK_EQUAL(bool(tfnsnpIntersections[4]), false);

  /// let's make a bunch of negative solution
  Intersection fIn(Vector3D(0., -1., 0.), -1., IntersectionStatus::reachable);
  Intersection sIn(Vector3D(0., -2., 0.), -2., IntersectionStatus::reachable);
  Intersection tIn(Vector3D(0., -3., 0.), -3., IntersectionStatus::reachable);

  std::vector<Intersection> tsfnIntersections = {tIn, sIn, fIn};
  std::vector<Intersection> fstnIntersections = {fIn, sIn, tIn};

  // this time around, sort the f-s-t-n to match the t-s-f-n
  std::sort(fstnIntersections.begin(), fstnIntersections.end());
  BOOST_CHECK_EQUAL(fstnIntersections[0].pathLength,
                    tsfnIntersections[0].pathLength);
  BOOST_CHECK_EQUAL(fstnIntersections[1].pathLength,
                    tsfnIntersections[1].pathLength);
  BOOST_CHECK_EQUAL(fstnIntersections[2].pathLength,
                    tsfnIntersections[2].pathLength);

  // let's sort them with greater
  std::sort(fstnIntersections.begin(), fstnIntersections.end(),
            std::greater<>());
  BOOST_CHECK_EQUAL(fstnIntersections[0].pathLength,
                    tsfnIntersections[2].pathLength);
  BOOST_CHECK_EQUAL(fstnIntersections[1].pathLength,
                    tsfnIntersections[1].pathLength);
  BOOST_CHECK_EQUAL(fstnIntersections[2].pathLength,
                    tsfnIntersections[0].pathLength);

  // shuffle negative and positive solutions
  std::vector<Intersection> pnsolutions = {tIp, sIn, sIp, fIn, tIn, fIp};
  std::sort(pnsolutions.begin(), pnsolutions.end());

  BOOST_CHECK_EQUAL(pnsolutions[0].pathLength, -3.);
  BOOST_CHECK_EQUAL(pnsolutions[1].pathLength, -2.);
  BOOST_CHECK_EQUAL(pnsolutions[2].pathLength, -1.);
  BOOST_CHECK_EQUAL(pnsolutions[3].pathLength, 1.);
  BOOST_CHECK_EQUAL(pnsolutions[4].pathLength, 2.);
  BOOST_CHECK_EQUAL(pnsolutions[5].pathLength, 3.);

  // sort intersections with zero path length
  Intersection zI(Vector3D(0., 0., 0.), 0., IntersectionStatus::reachable);
  std::vector<Intersection> tszfpIntersections = {tIp, sIp, zI, fIp};

  std::sort(tszfpIntersections.begin(), tszfpIntersections.end());
  BOOST_CHECK_EQUAL(tszfpIntersections[0].pathLength, 0.);
  BOOST_CHECK_EQUAL(tszfpIntersections[1].pathLength, 1.);
  BOOST_CHECK_EQUAL(tszfpIntersections[2].pathLength, 2.);
  BOOST_CHECK_EQUAL(tszfpIntersections[3].pathLength, 3.);

  std::vector<Intersection> tfsznIntersections = {tIn, fIn, sIn, zI};
  std::vector<Intersection> ztfsnIntersections = {zI, tIn, fIn, sIn};

  std::sort(tfsznIntersections.begin(), tfsznIntersections.end(),
            std::greater<>());
  BOOST_CHECK_EQUAL(tfsznIntersections[0].pathLength, 0.);
  BOOST_CHECK_EQUAL(tfsznIntersections[1].pathLength, -1.);
  BOOST_CHECK_EQUAL(tfsznIntersections[2].pathLength, -2.);
  BOOST_CHECK_EQUAL(tfsznIntersections[3].pathLength, -3.);

  std::sort(ztfsnIntersections.begin(), ztfsnIntersections.end(),
            std::greater<>());
  BOOST_CHECK_EQUAL(ztfsnIntersections[0].pathLength, 0.);
  BOOST_CHECK_EQUAL(ztfsnIntersections[1].pathLength, -1.);
  BOOST_CHECK_EQUAL(ztfsnIntersections[2].pathLength, -2.);
  BOOST_CHECK_EQUAL(ztfsnIntersections[3].pathLength, -3.);
}

/// test of the object intersection class
BOOST_AUTO_TEST_CASE(ObjectIntersectionTest) {
  auto psf6 = Surface::makeShared<PlaneSurface>(Vector3D(6., 0., 0.),
                                                Vector3D(1., 0., 0.));
  auto psf7 = Surface::makeShared<PlaneSurface>(Vector3D(7., 0., 0.),
                                                Vector3D(1., 0., 0.));
  auto psf8 = Surface::makeShared<PlaneSurface>(Vector3D(8., 0., 0.),
                                                Vector3D(1., 0., 0.));
  auto psf9 = Surface::makeShared<PlaneSurface>(Vector3D(9., 0., 0.),
                                                Vector3D(1., 0., 0.));
  auto psf10 = Surface::makeShared<PlaneSurface>(Vector3D(10., 0., 0.),
                                                 Vector3D(1., 0., 0.));

  using PlaneIntersection = ObjectIntersection<PlaneSurface>;

  PlaneIntersection int6(
      Intersection(Vector3D(6., 0., 0.), 6., IntersectionStatus::reachable),
      psf6.get());
  PlaneIntersection int7(
      Intersection(Vector3D(7., 0., 0.), 7., IntersectionStatus::reachable),
      psf7.get());
  PlaneIntersection int8(
      Intersection(Vector3D(8., 0., 0.), 8., IntersectionStatus::reachable),
      psf8.get());
  PlaneIntersection int9a(
      Intersection(Vector3D(9., 0., 0.), 9., IntersectionStatus::reachable),
      psf9.get());
  PlaneIntersection int9b(
      Intersection(Vector3D(9., 1., 0.), std::sqrt(9. * 9. + 1.),
                   IntersectionStatus::reachable),
      psf9.get());
  PlaneIntersection int10(
      Intersection(Vector3D(10., 0., 0.), 10., IntersectionStatus::reachable),
      psf10.get());
}

}  // namespace Test
}  // namespace Acts
