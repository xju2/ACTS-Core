// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// clang-format off
#define BOOST_TEST_MODULE Overstepping Tests
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/output_test_stream.hpp>
// clang-format on

#include <memory>

#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Tests/CommonHelpers/CylindricalTrackingGeometry.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Units.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "Acts/Propagator/detail/DebugOutputActor.hpp"

namespace bdata = boost::unit_test::data;
namespace tt = boost::test_tools;

/// This test emulates over-stepping to a surface, under-stepping is not
/// a problem at all, as an additional step will just be inserted to compensate
/// overstepping, however, could result in missing the surface and within
/// tolerance should be reverted.
namespace Acts {
namespace Test {

// Create a test context
GeometryContext tgContext = GeometryContext();
MagneticFieldContext mfContext = MagneticFieldContext();

// Global definitions
// The path limit abort
using path_limit = detail::PathLimitReached;

CylindricalTrackingGeometry cGeometry(tgContext);
auto tGeometry = cGeometry();

const int ntests = 1;
const int skip = 0;

bool debugMode = true;

/// A Navigator that oversteps when targetting for a delta
class OversteppingNavigator : public Navigator {
 public:
  /// Constructor with shared tracking geometry
  ///
  /// @param tGeometry The tracking geometry for the navigator
  OversteppingNavigator(std::shared_ptr<const TrackingGeometry> trkGeo,
                        double deltaStep = 0. * units::_um)
      : Navigator(std::move(trkGeo)), m_deltaStep(deltaStep) {}

  /// @brief Navigator target call
  ///
  /// Overwrites teh Navigator::target call and adds the delta
  ///
  /// @tparam propagator_state_t is the type of Propagatgor state
  /// @tparam stepper_t is the used type of the Stepper by the Propagator
  ///
  /// @param [in,out] state is the mutable propagator state object
  /// @param [in] stepper Stepper in use
  template <typename propagator_state_t, typename stepper_t>
  void target(propagator_state_t& state, const stepper_t& stepper) const {
    // Call the BaseClass::target
    Navigator::target(state, stepper);
    // Modify the stepSize target
    double actorStep = state.stepping.stepSize.value(Cstep::actor);
    state.stepping.stepSize.update(actorStep + m_deltaStep, Cstep::actor, true);

    debugLog(state, [&] {
      std::stringstream dstream;
      dstream << "Navigation forced to overstep by " << m_deltaStep << " to ";
      dstream << state.stepping.stepSize.toString();
      return dstream.str();
    });

    // And now return
    return;
  }

 private:
  double m_deltaStep{0.};
};

// create a navigator for this tracking geometry
Navigator referenceNavigator(tGeometry);
OversteppingNavigator oversteppingNavigator(tGeometry, 50. * units::_um);

using BField = ConstantBField;
const double Bz = 2. * units::_T;
BField bField(0, 0, Bz);

// the actual test nethod that runs the test
/// can be used with several propagator types
/// @tparam stepper_t is the stepper type for the test
///
/// @param prop is the propagator instance
/// @param pT the transverse momentum
/// @param phi the azimuthal angle of the track at creation
/// @param theta the polar angle of the track at creation
/// @param charge is the charge of the particle
/// @param index is the run index from the test
template <typename stepper_t>
void runTest(const stepper_t& stepper, double pT, double phi, double theta,
             int charge, int index) {
  double dcharge = -1 + 2 * charge;

  if (index < skip) {
    return;
  }

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
  CurvilinearParameters start(nullptr, pos, mom, q);

  using DebugOutput = detail::DebugOutputActor;

  // Action list and abort list
  using ActionListType = ActionList<DebugOutput>;
  using AbortListType = AbortList<>;

  // declare the overstepper and standard propagator
  using Referencer = Propagator<stepper_t, Navigator>;
  using Overstepper = Propagator<stepper_t, OversteppingNavigator>;
  using Options = PropagatorOptions<ActionListType, AbortListType>;

  Referencer referencePropagator(stepper, referenceNavigator);
  Overstepper overstepPropagator(stepper, oversteppingNavigator);

  Options options(tgContext, mfContext);
  options.debug = debugMode;

  // propagation with the standard propagator
  if (debugMode) {
    std::cout << ">>> Propagation with standard Navigator: start." << std::endl;
  }

  auto referenceResult = referencePropagator.propagate(start, options);
  auto overstepResult = overstepPropagator.propagate(start, options);

  // propagation with the standard propagator
  if (debugMode && referenceResult.ok()) {
    std::cout << ">>> Result achieved with " << (*referenceResult).steps
              << " steps." << std::endl;
    const auto& referenceOutput =
        (*referenceResult).template get<DebugOutput::result_type>();
    std::cout << referenceOutput.debugString << std::endl;
  }

  // propagation with the overstepping propagator
  if (debugMode) {
    std::cout << ">>> Propagation with overstepping Navigator: start."
              << std::endl;
  }

  // propagation with the standard propagator
  if (debugMode && overstepResult.ok()) {
    std::cout << ">>> Result achieved with " << (*overstepResult).steps
              << " steps." << std::endl;
    const auto& oversteppingOutput =
        (*overstepResult).template get<DebugOutput::result_type>();
    std::cout << oversteppingOutput.debugString << std::endl;
  }
}

// This test case checks that no segmentation fault appears
// - this tests the collection of surfaces
BOOST_DATA_TEST_CASE(
    overstepping_test_straightline,
    bdata::random((bdata::seed = 20,
                   bdata::distribution = std::uniform_real_distribution<>(
                       0.5 * units::_GeV, 10. * units::_GeV))) ^
        bdata::random((bdata::seed = 21,
                       bdata::distribution =
                           std::uniform_real_distribution<>(-M_PI, M_PI))) ^
        bdata::random((bdata::seed = 22,
                       bdata::distribution =
                           std::uniform_real_distribution<>(1.0, M_PI - 1.0))) ^
        bdata::random(
            (bdata::seed = 23,
             bdata::distribution = std::uniform_int_distribution<>(0, 1))) ^
        bdata::xrange(ntests),
    pT, phi, theta, charge, index) {
  // The straight line stepper
  StraightLineStepper slstepper;
  runTest(slstepper, pT, phi, theta, charge, index);
}

}  // namespace Test
}  // namespace Acts
