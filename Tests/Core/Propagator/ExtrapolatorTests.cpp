// This file is part of the Acts project.
//
// Copyright (C) 2018-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// clang-format off
#define BOOST_TEST_MODULE Extrapolator Tests
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/output_test_stream.hpp>
// clang-format on

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Propagator/MaterialInteractor.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/SurfaceCollector.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Propagator/ActionList.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/detail/DebugOutputActor.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Tests/CommonHelpers/CylindricalTrackingGeometry.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Units.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"

namespace bdata = boost::unit_test::data;
namespace tt = boost::test_tools;
using namespace Acts::UnitLiterals;

namespace Acts {
namespace Test {

// Create a test context
GeometryContext tgContext = GeometryContext();
MagneticFieldContext mfContext = MagneticFieldContext();

// Global definitions
// The path limit abort
using path_limit = detail::PathLimitReached;

std::vector<std::unique_ptr<const Surface>> stepState;

CylindricalTrackingGeometry cGeometry(tgContext);
auto tGeometry = cGeometry();

// Get the navigator and provide the TrackingGeometry
Navigator navigator(tGeometry);

using BFieldType = ConstantBField;
using EigenStepperType = EigenStepper<BFieldType>;
using EigenPropagatorType = Propagator<EigenStepperType, Navigator>;
using Covariance = BoundSymMatrix;

BFieldType bField(0, 0, 2_T);
EigenStepperType estepper(bField);
EigenPropagatorType epropagator(std::move(estepper), std::move(navigator));

const int ntests = 100;
bool debugMode = false;

// A plane selector for the SurfaceCollector
struct PlaneSelector {
  /// Call operator
  /// @param sf The input surface to be checked
  bool operator()(const Surface& sf) const {
    return (sf.type() == Surface::Plane);
  }
};

// This test case checks that no segmentation fault appears
// - simple extrapolation test
BOOST_DATA_TEST_CASE(
    test_extrapolation_,
    bdata::random((bdata::seed = 0,
                   bdata::distribution =
                       std::uniform_real_distribution<>(0.4_GeV, 10_GeV))) ^
        bdata::random((bdata::seed = 1,
                       bdata::distribution =
                           std::uniform_real_distribution<>(-M_PI, M_PI))) ^
        bdata::random((bdata::seed = 2,
                       bdata::distribution =
                           std::uniform_real_distribution<>(1.0, M_PI - 1.0))) ^
        bdata::random(
            (bdata::seed = 3,
             bdata::distribution = std::uniform_int_distribution<>(0, 1))) ^
        bdata::random(
            (bdata::seed = 4,
             bdata::distribution = std::uniform_int_distribution<>(0, 100))) ^
        bdata::xrange(ntests),
    pT, phi, theta, charge, time, index) {
  double dcharge = -1 + 2 * charge;
  (void)index;

  // define start parameters
  double x = 0;
  double y = 0;
  double z = 0;
  double px = pT * cos(phi);
  double py = pT * sin(phi);
  double pz = pT / tan(theta);
  double q = dcharge;
  Vector3D pos(x, y, z);
  Vector3D mom(px, py, pz);
  /// a covariance matrix to transport
  Covariance cov;
  // take some major correlations (off-diagonals)
  cov << 10_mm, 0, 0.123, 0, 0.5, 0, 0, 10_mm, 0, 0.162, 0, 0, 0.123, 0, 0.1, 0,
      0, 0, 0, 0.162, 0, 0.1, 0, 0, 0.5, 0, 0, 0, 1. / (10_GeV), 0, 0, 0, 0, 0,
      0, 0;
  CurvilinearParameters start(cov, pos, mom, q, time);

  PropagatorOptions<> options(tgContext, mfContext);
  options.maxStepSize = 10_cm;
  options.pathLimit = 25_cm;

  BOOST_CHECK(epropagator.propagate(start, options).value().endParameters !=
              nullptr);
}

// This test case checks that no segmentation fault appears
// - this tests the collection of surfaces
BOOST_DATA_TEST_CASE(
    test_surface_collection_,
    bdata::random((bdata::seed = 10,
                   bdata::distribution =
                       std::uniform_real_distribution<>(0.4_GeV, 10_GeV))) ^
        bdata::random((bdata::seed = 11,
                       bdata::distribution =
                           std::uniform_real_distribution<>(-M_PI, M_PI))) ^
        bdata::random((bdata::seed = 12,
                       bdata::distribution =
                           std::uniform_real_distribution<>(1.0, M_PI - 1.0))) ^
        bdata::random(
            (bdata::seed = 13,
             bdata::distribution = std::uniform_int_distribution<>(0, 1))) ^
        bdata::random(
            (bdata::seed = 14,
             bdata::distribution = std::uniform_int_distribution<>(0, 100))) ^
        bdata::xrange(ntests),
    pT, phi, theta, charge, time, index) {
  double dcharge = -1 + 2 * charge;
  (void)index;

  // define start parameters
  double x = 0;
  double y = 0;
  double z = 0;
  double px = pT * cos(phi);
  double py = pT * sin(phi);
  double pz = pT / tan(theta);
  double q = dcharge;
  Vector3D pos(x, y, z);
  Vector3D mom(px, py, pz);
  /// a covariance matrix to transport
  Covariance cov;
  // take some major correlations (off-diagonals)
  cov << 10_mm, 0, 0.123, 0, 0.5, 0, 0, 10_mm, 0, 0.162, 0, 0, 0.123, 0, 0.1, 0,
      0, 0, 0, 0.162, 0, 0.1, 0, 0, 0.5, 0, 0, 0, 1. / (10_GeV), 0, 0, 0, 0, 0,
      0, 0;
  CurvilinearParameters start(cov, pos, mom, q, time);

  // A PlaneSelector for the SurfaceCollector
  using PlaneCollector = SurfaceCollector<PlaneSelector>;

  PropagatorOptions<ActionList<PlaneCollector>> options(tgContext, mfContext);

  options.maxStepSize = 10_cm;
  options.pathLimit = 25_cm;
  options.debug = debugMode;

  const auto& result = epropagator.propagate(start, options).value();
  auto collector_result = result.get<PlaneCollector::result_type>();

  // step through the surfaces and go step by step
  PropagatorOptions<> optionsEmpty(tgContext, mfContext);

  optionsEmpty.maxStepSize = 25_cm;
  optionsEmpty.debug = true;
  // Try propagation from start to each surface
  for (const auto& colsf : collector_result.collected) {
    const auto& csurface = colsf.surface;
    // Avoid going to the same surface
    // @todo: decide on strategy and write unit test for this
    if (csurface == &start.referenceSurface()) {
      continue;
    }
    // Extrapolate & check
    const auto& cresult = epropagator.propagate(start, *csurface, optionsEmpty)
                              .value()
                              .endParameters;
    BOOST_CHECK(cresult != nullptr);
  }
}

// This test case checks that no segmentation fault appears
// - this tests the collection of surfaces
BOOST_DATA_TEST_CASE(
    test_material_interactor_,
    bdata::random((bdata::seed = 20,
                   bdata::distribution =
                       std::uniform_real_distribution<>(0.4_GeV, 10_GeV))) ^
        bdata::random((bdata::seed = 21,
                       bdata::distribution =
                           std::uniform_real_distribution<>(-M_PI, M_PI))) ^
        bdata::random((bdata::seed = 22,
                       bdata::distribution =
                           std::uniform_real_distribution<>(1.0, M_PI - 1.0))) ^
        bdata::random(
            (bdata::seed = 23,
             bdata::distribution = std::uniform_int_distribution<>(0, 1))) ^
        bdata::random(
            (bdata::seed = 24,
             bdata::distribution = std::uniform_int_distribution<>(0, 100))) ^
        bdata::xrange(ntests),
    pT, phi, theta, charge, time, index) {
  double dcharge = -1 + 2 * charge;
  (void)index;

  // define start parameters
  double x = 0;
  double y = 0;
  double z = 0;
  double px = pT * cos(phi);
  double py = pT * sin(phi);
  double pz = pT / tan(theta);
  double q = dcharge;
  Vector3D pos(x, y, z);
  Vector3D mom(px, py, pz);
  /// a covariance matrix to transport
  Covariance cov;
  // take some major correlations (off-diagonals)
  cov << 10_mm, 0, 0.123, 0, 0.5, 0, 0, 10_mm, 0, 0.162, 0, 0, 0.123, 0, 0.1, 0,
      0, 0, 0, 0.162, 0, 0.1, 0, 0, 0.5, 0, 0, 0, 1. / (10_GeV), 0, 0, 0, 0, 0,
      0, 0;
  CurvilinearParameters start(cov, pos, mom, q, time);

  using DebugOutput = detail::DebugOutputActor;

  PropagatorOptions<ActionList<MaterialInteractor, DebugOutput>> options(
      tgContext, mfContext);
  options.debug = debugMode;
  options.maxStepSize = 25_cm;
  options.pathLimit = 25_cm;

  const auto& result = epropagator.propagate(start, options).value();
  if (result.endParameters) {
    // test that you actually lost some energy
    BOOST_CHECK_LT(result.endParameters->momentum().norm(),
                   start.momentum().norm());
  }

  if (debugMode) {
    const auto& output = result.get<DebugOutput::result_type>();
    std::cout << ">>> Extrapolation output " << std::endl;
    std::cout << output.debugString << std::endl;
  }
}

// This test case checks that no segmentation fault appears
// - this tests the loop protection
BOOST_DATA_TEST_CASE(
    loop_protection_test,
    bdata::random((bdata::seed = 20,
                   bdata::distribution =
                       std::uniform_real_distribution<>(0.1_GeV, 0.5_GeV))) ^
        bdata::random((bdata::seed = 21,
                       bdata::distribution =
                           std::uniform_real_distribution<>(-M_PI, M_PI))) ^
        bdata::random((bdata::seed = 22,
                       bdata::distribution =
                           std::uniform_real_distribution<>(1.0, M_PI - 1.0))) ^
        bdata::random(
            (bdata::seed = 23,
             bdata::distribution = std::uniform_int_distribution<>(0, 1))) ^
        bdata::random(
            (bdata::seed = 24,
             bdata::distribution = std::uniform_int_distribution<>(0, 100))) ^
        bdata::xrange(ntests),
    pT, phi, theta, charge, time, index) {
  double dcharge = -1 + 2 * charge;
  (void)index;

  // define start parameters
  double x = 0;
  double y = 0;
  double z = 0;
  double px = pT * cos(phi);
  double py = pT * sin(phi);
  double pz = pT / tan(theta);
  double q = dcharge;
  Vector3D pos(x, y, z);
  Vector3D mom(px, py, pz);
  /// a covariance matrix to transport
  Covariance cov;
  // take some major correlations (off-diagonals)
  cov << 10_mm, 0, 0.123, 0, 0.5, 0, 0, 10_mm, 0, 0.162, 0, 0, 0.123, 0, 0.1, 0,
      0, 0, 0, 0.162, 0, 0.1, 0, 0, 0.5, 0, 0, 0, 1. / (10_GeV), 0, 0, 0, 0, 0,
      0, 0;
  CurvilinearParameters start(cov, pos, mom, q, time);

  // Action list and abort list
  using DebugOutput = detail::DebugOutputActor;

  PropagatorOptions<ActionList<MaterialInteractor, DebugOutput>> options(
      tgContext, mfContext);
  options.maxStepSize = 25_cm;
  options.pathLimit = 1500_mm;

  const auto& status = epropagator.propagate(start, options).value();
  // this test assumes state.options.loopFraction = 0.5
  // maximum momentum allowed
  double pmax = options.pathLimit * bField.getField(pos).norm() / M_PI;
  if (mom.norm() < pmax) {
    BOOST_CHECK_LT(status.pathLength, options.pathLimit);
  } else {
    CHECK_CLOSE_REL(status.pathLength, options.pathLimit, 1e-3);
  }
}

}  // namespace Test
}  // namespace Acts
