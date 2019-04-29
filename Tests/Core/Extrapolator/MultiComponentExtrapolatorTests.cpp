// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
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
#include "Acts/Extrapolator/MaterialInteractor.hpp"
#include "Acts/Extrapolator/MultiMaterialInteractor.hpp"
#include "Acts/Extrapolator/Navigator.hpp"
#include "Acts/Extrapolator/SurfaceCollector.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Propagator/ActionList.hpp"
#include "Acts/Propagator/MultiEigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/detail/DebugOutputActor.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Tests/CommonHelpers/CylindricalTrackingGeometry.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Units.hpp"
#include "Acts/Utilities/GeometryContext.hpp"
#include "Acts/Utilities/MagneticFieldContext.hpp"

namespace bdata = boost::unit_test::data;
namespace tt    = boost::test_tools;

namespace Acts {

namespace Test {

  // Create a test context
  GeometryContext      tgContext = GeometryContext();
  MagneticFieldContext mfContext = MagneticFieldContext();

  // Global definitions
  // The path limit abort
  using path_limit = detail::PathLimitReached;

  std::vector<std::unique_ptr<const Surface>> stepState;

  CylindricalTrackingGeometry cGeometry(tgContext);
  auto                        tGeometry = cGeometry();

  // Get the navigator and provide the TrackingGeometry
  Navigator navigator(tGeometry);
  Navigator multi_navigator(tGeometry);

  using BFieldType          = ConstantBField;
  using EigenStepperType    = MultiEigenStepper<BFieldType>;
  using EigenPropagatorType = Propagator<EigenStepperType, Navigator>;
  using MultiEigenStepperType    = MultiEigenStepper<BFieldType>;
  using MultiEigenPropagatorType = Propagator<MultiEigenStepperType, Navigator>;

  const double        Bz = 2. * units::_T;
  BFieldType          bField(0, 0, Bz);

  // define the scs
  EigenStepperType    estepper(bField);
  EigenPropagatorType epropagator(std::move(estepper), std::move(navigator));

  // define the mcs
  MultiEigenStepperType    multi_estepper(bField);
  MultiEigenPropagatorType multi_epropagator(std::move(multi_estepper),
                                             std::move(multi_navigator));

  const int ntests    = 10;
  bool      debugMode = true;

  // A plane selector for the SurfaceCollector
  struct PlaneSelector
  {
    /// Call operator
    /// @param sf The input surface to be checked
    bool
    operator()(const Surface& sf) const
    {
      return (sf.type() == Surface::Plane);
    }
  };

  // This test case checks that no segmentation fault appears in the mcs
  // the basic multi stepper propagate
  BOOST_DATA_TEST_CASE(
      test_mcs_extrapolation_,
      bdata::random((bdata::seed = 0,
                     bdata::distribution
                     = std::uniform_real_distribution<>(0.4 * units::_GeV,
                                                        10. * units::_GeV)))
          ^ bdata::random((bdata::seed = 1,
                           bdata::distribution
                           = std::uniform_real_distribution<>(-M_PI, M_PI)))
          ^ bdata::random((bdata::seed = 2,
                           bdata::distribution
                           = std::uniform_real_distribution<>(1.0, M_PI - 1.0)))
          ^ bdata::random((bdata::seed = 3,
                           bdata::distribution
                           = std::uniform_int_distribution<>(0, 1)))
          ^ bdata::xrange(ntests),
      pT,
      phi,
      theta,
      charge,
      index)
  {

    double dcharge = -1 + 2 * charge;
    (void)index;

    // define start parameters
    double   x  = 0;
    double   y  = 0;
    double   z  = 0;
    double   px = pT * cos(phi);
    double   py = pT * sin(phi);
    double   pz = pT / tan(theta);
    double   q  = dcharge;
    Vector3D pos(x, y, z);
    Vector3D mom(px, py, pz);
    /// a covariance matrix to transport
    ActsSymMatrixD<5> cov;
    // take some major correlations (off-diagonals)
    cov << 10 * units::_mm, 0, 0.123, 0, 0.5, 0, 10 * units::_mm, 0, 0.162, 0,
        0.123, 0, 0.1, 0, 0, 0, 0.162, 0, 0.1, 0, 0.5, 0, 0, 0,
        1. / (10 * units::_GeV);
    auto covPtr = std::make_unique<const ActsSymMatrixD<5>>(cov);
    CurvilinearParameters start(std::move(covPtr), pos, mom, q);

    using DebugOutput = detail::DebugOutputActor;

    PropagatorOptions<ActionList<DebugOutput>> options(
        tgContext, mfContext);
    options.debug       = debugMode;
    options.maxStepSize = 10. * units::_cm;
    options.pathLimit   = 25 * units::_cm;

    const auto& result = multi_epropagator.propagate(start, options).value();
    if (debugMode) {
      const auto& output = result.get<DebugOutput::result_type>();
      std::cout << ">>> Extrapolation output " << std::endl;
      std::cout << output.debugString << std::endl;
    }
  }
  /*
  // This test case checks that no segmentation fault appears
  // - this tests the same surfaceHit of different stepper
  BOOST_DATA_TEST_CASE(
      test_equal_scs_mcs_collection_,
      bdata::random((bdata::seed = 0,
                     bdata::distribution
                     = std::uniform_real_distribution<>(0.4 * units::_GeV,
                                                        10. * units::_GeV)))
          ^ bdata::random((bdata::seed = 1,
                           bdata::distribution
                           = std::uniform_real_distribution<>(-M_PI, M_PI)))
          ^ bdata::random((bdata::seed = 2,
                           bdata::distribution
                           = std::uniform_real_distribution<>(1.0, M_PI - 1.0)))
          ^ bdata::random((bdata::seed = 3,
                           bdata::distribution
                           = std::uniform_int_distribution<>(0, 1)))
          ^ bdata::xrange(ntests),
      pT,
      phi,
      theta,
      charge,
      index)
  {

    double dcharge = -1 + 2 * charge;
    (void)index;

    // define start parameters
    double   x  = 0;
    double   y  = 0;
    double   z  = 0;
    double   px = pT * cos(phi);
    double   py = pT * sin(phi);
    double   pz = pT / tan(theta);
    double   q  = dcharge;
    Vector3D pos(x, y, z);
    Vector3D mom(px, py, pz);
    /// a covariance matrix to transport
    ActsSymMatrixD<5> cov;
    // take some major correlations (off-diagonals)
    cov << 10 * units::_mm, 0, 0.123, 0, 0.5, 0, 10 * units::_mm, 0, 0.162, 0,
        0.123, 0, 0.1, 0, 0, 0, 0.162, 0, 0.1, 0, 0.5, 0, 0, 0,
        1. / (10 * units::_GeV);
    auto covPtr = std::make_unique<const ActsSymMatrixD<5>>(cov);
    CurvilinearParameters start(std::move(covPtr), pos, mom, q);

    // A PlaneSelector for the SurfaceCollector
    using PlaneCollector = SurfaceCollector<PlaneSelector>;

    PropagatorOptions<ActionList<PlaneCollector>> options(
        tgContext, mfContext);
    options.maxStepSize       = 10. * units::_cm;
    options.pathLimit         = 25 * units::_cm;
    options.debug             = debugMode;

    PropagatorOptions<ActionList<PlaneCollector>> multi_options(
        tgContext, mfContext);
    multi_options.maxStepSize = 10. * units::_cm;
    multi_options.pathLimit   = 25 * units::_cm;
    multi_options.debug       = debugMode;

    using DebugOutput = detail::DebugOutputActor;
    PropagatorOptions<ActionList<DebugOutput,PlaneCollector,MultiMaterialInteractor>> multi_material_options(
        tgContext, mfContext);
    multi_material_options.maxStepSize = 10. * units::_cm;
    multi_material_options.pathLimit   = 25 * units::_cm;
    multi_material_options.debug       = debugMode;

    // sigle component
    const auto& result           = epropagator.propagate(start, options).value();
    auto        collector_result = result.get<PlaneCollector::result_type>();

    // multi component
    const auto& multi_result 	 = multi_epropagator.propagate(start, multi_options).value();
    auto multi_collector_result  = multi_result.get<PlaneCollector::result_type>();

	// multi_material component
    const auto& multi_material_result 	 = multi_epropagator.propagate(start, multi_material_options).value();
    auto multi_material_collector_result  = multi_material_result.get<PlaneCollector::result_type>();

    BOOST_CHECK_EQUAL(collector_result.collected.size(),
                      multi_collector_result.collected.size());
    BOOST_CHECK(collector_result.collected == multi_collector_result.collected);

    if (debugMode) {
      const auto& output = multi_material_result.get<DebugOutput::result_type>();
      std::cout << ">>> Extrapolation output " << std::endl;
      std::cout << output.debugString << std::endl;
    }
  }
  */

}  // namespace Test
}  // namespace Acts
