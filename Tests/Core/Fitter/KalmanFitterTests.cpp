// This file is part of the Acts project.
//
// Copyright (C) 2016-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// clang-format off
#define BOOST_TEST_MODULE KalmanFitter Tests
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
// clang-format on

#include <algorithm>
#include <cmath>
#include <random>
#include <vector>

#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/EventData/Measurement.hpp"
#include "Acts/EventData/MeasurementHelpers.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/TrackState.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/SurfaceCollector.hpp"
#include "Acts/Fitter/GainMatrixSmoother.hpp"
#include "Acts/Fitter/GainMatrixUpdater.hpp"
#include "Acts/Fitter/KalmanFitter.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "Acts/Propagator/detail/DebugOutputActor.hpp"
#include "Acts/Propagator/detail/StandardAborters.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Tests/CommonHelpers/CubicTrackingGeometry.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Geometry/GeometryID.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Utilities/CalibrationContext.hpp"

using namespace Acts::UnitLiterals;

namespace Acts {
namespace Test {

// A few initialisations and definitionas
using SourceLink = MinimalSourceLink;
using Jacobian = BoundParameters::CovMatrix_t;
using Covariance = BoundSymMatrix;

using TrackState = TrackState<SourceLink, BoundParameters>;
using Resolution = std::pair<ParID_t, double>;
using ElementResolution = std::vector<Resolution>;
using VolumeResolution = std::map<geo_id_value, ElementResolution>;
using DetectorResolution = std::map<geo_id_value, VolumeResolution>;

using DebugOutput = detail::DebugOutputActor;

std::normal_distribution<double> gauss(0., 1.);
std::default_random_engine generator(42);

ActsSymMatrixD<1> cov1D;
ActsSymMatrixD<2> cov2D;

bool debugMode = false;

// Create a test context
GeometryContext tgContext = GeometryContext();
MagneticFieldContext mfContext = MagneticFieldContext();
CalibrationContext calContext = CalibrationContext();

template <ParID_t... params>
using MeasurementType = Measurement<SourceLink, params...>;

/// @brief This struct creates FittableMeasurements on the
/// detector surfaces, according to the given smearing xxparameters
///
struct MeasurementCreator {
  /// @brief Constructor
  MeasurementCreator() = default;

  /// The detector resolution
  DetectorResolution detectorResolution;

  using result_type = std::vector<FittableMeasurement<SourceLink>>;

  /// @brief Operater that is callable by an ActionList. The function collects
  /// the surfaces
  ///
  /// @tparam propagator_state_t Type of the propagator state
  /// @tparam stepper_t Type of the stepper
  /// @param [in] state State of the propagator
  /// @param [out] result Vector of matching surfaces
  template <typename propagator_state_t, typename stepper_t>
  void operator()(propagator_state_t& state, const stepper_t& stepper,
                  result_type& result) const {
    // monitor the current surface
    auto surface = state.navigation.currentSurface;
    if (surface and surface->associatedDetectorElement()) {
      auto geoID = surface->geoID();
      geo_id_value volumeID = geoID.value(GeometryID::volume_mask);
      geo_id_value layerID = geoID.value(GeometryID::layer_mask);
      // find volume and layer information for this
      auto vResolution = detectorResolution.find(volumeID);
      if (vResolution != detectorResolution.end()) {
        // find layer resolutions
        auto lResolution = vResolution->second.find(layerID);
        if (lResolution != vResolution->second.end()) {
          // Apply global to local
          Acts::Vector2D lPos;
          surface->globalToLocal(state.geoContext,
                                 stepper.position(state.stepping),
                                 stepper.direction(state.stepping), lPos);
          if (lResolution->second.size() == 1) {
            double sp = lResolution->second[0].second;
            cov1D << sp * sp;
            double dp = sp * gauss(generator);
            if (lResolution->second[0].first == eLOC_0) {
              // push back & move a LOC_0 measurement
              MeasurementType<eLOC_0> m0(surface->getSharedPtr(), {}, cov1D,
                                         lPos[eLOC_0] + dp);
              result.push_back(std::move(m0));
            } else {
              // push back & move a LOC_1 measurement
              MeasurementType<eLOC_1> m1(surface->getSharedPtr(), {}, cov1D,
                                         lPos[eLOC_1] + dp);
              result.push_back(std::move(m1));
            }
          } else if (lResolution->second.size() == 2) {
            // Create the measurment and move it
            double sx = lResolution->second[eLOC_0].second;
            double sy = lResolution->second[eLOC_1].second;
            cov2D << sx * sx, 0., 0., sy * sy;
            double dx = sx * gauss(generator);
            double dy = sy * gauss(generator);
            // push back & move a LOC_0, LOC_1 measurement
            MeasurementType<eLOC_0, eLOC_1> m01(surface->getSharedPtr(), {},
                                                cov2D, lPos[eLOC_0] + dx,
                                                lPos[eLOC_1] + dy);
            result.push_back(std::move(m01));
          }
        }
      }
    }
  }
};

double dX, dY;
Vector3D pos;
const Surface* sur;

///
/// @brief Simplified material interaction effect by pure gaussian
/// deflection
///
struct MaterialScattering {
  /// @brief Constructor
  MaterialScattering() = default;

  /// @brief Main action list call operator for the scattering on material
  ///
  /// @todo deal momentum in a gaussian way properly
  ///
  /// @tparam propagator_state_t State of the propagator
  /// @param stepper_t Type of the stepper
  /// @param [in] state State of the propagation
  /// @param [in] stepper Stepper of the propagation
  template <typename propagator_state_t, typename stepper_t>
  void operator()(propagator_state_t& state, const stepper_t& stepper) const {
    // Check if there is a surface with material and a covariance is existing
    if (state.navigation.currentSurface &&
        state.navigation.currentSurface->surfaceMaterial() &&
        state.stepping.cov != Covariance::Zero()) {
      // Sample angles
      std::normal_distribution<double> scatterAngle(
          0., 0.017);  //< \approx 1 degree
      double dPhi = scatterAngle(generator), dTheta = scatterAngle(generator);

      // Update the covariance
      state.stepping.cov(ePHI, ePHI) += dPhi * dPhi;
      state.stepping.cov(eTHETA, eTHETA) += dTheta * dTheta;

      // Update the angles
      auto direction = stepper.direction(state.stepping);
      double theta = std::acos(direction.z());
      double phi = std::atan2(direction.y(), direction.x());

      state.stepping.update(
          stepper.position(state.stepping),
          {std::sin(theta + dTheta) * std::cos(phi + dPhi),
           std::sin(theta + dTheta) * std::sin(phi + dPhi),
           std::cos(theta + dTheta)},
          std::max(stepper.momentum(state.stepping) -
                       std::abs(gauss(generator)) * UnitConstants::MeV,
                   0.));
    }
  }
};

///
/// @brief Unit test for Kalman fitter with measurements along the x-axis
///
BOOST_AUTO_TEST_CASE(kalman_fitter_zero_field) {
  // Build detector
  CubicTrackingGeometry cGeometry(tgContext);
  auto detector = cGeometry();

  // Build navigator for the measurement creatoin
  Navigator mNavigator(detector);
  mNavigator.resolvePassive = false;
  mNavigator.resolveMaterial = true;
  mNavigator.resolveSensitive = true;

  // Use straingt line stepper to create the measurements
  StraightLineStepper mStepper;

  // Define the measurement propagator
  using MeasurementPropagator = Propagator<StraightLineStepper, Navigator>;

  // Build propagator for the measurement creation
  MeasurementPropagator mPropagator(mStepper, mNavigator);
  Vector3D mPos(-3_m, 0., 0.), mMom(1_GeV, 0., 0);
  SingleCurvilinearTrackParameters<NeutralPolicy> mStart(std::nullopt, mPos,
                                                         mMom, 42_ns);

  // Create action list for the measurement creation
  using MeasurementActions = ActionList<MeasurementCreator, DebugOutput>;
  using MeasurementAborters = AbortList<detail::EndOfWorldReached>;

  auto pixelResX = Resolution(eLOC_0, 25_um);
  auto pixelResY = Resolution(eLOC_1, 50_um);
  auto stripResX = Resolution(eLOC_0, 100_um);
  auto stripResY = Resolution(eLOC_1, 150_um);

  ElementResolution pixelElementRes = {pixelResX, pixelResY};
  ElementResolution stripElementResI = {stripResX};
  ElementResolution stripElementResO = {stripResY};

  VolumeResolution pixelVolumeRes;
  pixelVolumeRes[2] = pixelElementRes;
  pixelVolumeRes[4] = pixelElementRes;

  VolumeResolution stripVolumeRes;
  stripVolumeRes[2] = stripElementResI;
  stripVolumeRes[4] = stripElementResO;
  stripVolumeRes[6] = stripElementResI;
  stripVolumeRes[8] = stripElementResO;

  DetectorResolution detRes;
  detRes[2] = pixelVolumeRes;
  detRes[3] = stripVolumeRes;

  // Set options for propagator
  PropagatorOptions<MeasurementActions, MeasurementAborters> mOptions(
      tgContext, mfContext);
  mOptions.debug = debugMode;
  auto& mCreator = mOptions.actionList.get<MeasurementCreator>();
  mCreator.detectorResolution = detRes;

  // Launch and collect - the measurements
  auto mResult = mPropagator.propagate(mStart, mOptions).value();
  if (debugMode) {
    const auto debugString =
        mResult.template get<DebugOutput::result_type>().debugString;
    std::cout << ">>>> Measurement creation: " << std::endl;
    std::cout << debugString;
  }

  // Extract measurements from result of propagation.
  // This vector owns the measurements
  std::vector<FittableMeasurement<SourceLink>> measurements =
      std::move(mResult.template get<MeasurementCreator::result_type>());
  BOOST_CHECK_EQUAL(measurements.size(), 6);

  // Make a vector of source links as input to the KF
  std::vector<SourceLink> sourcelinks;
  std::transform(measurements.begin(), measurements.end(),
                 std::back_inserter(sourcelinks),
                 [](const auto& m) { return SourceLink{&m}; });

  // The KalmanFitter - we use the eigen stepper for covariance transport
  // Build navigator for the measurement creatoin
  Navigator rNavigator(detector);
  rNavigator.resolvePassive = false;
  rNavigator.resolveMaterial = true;
  rNavigator.resolveSensitive = true;

  // Configure propagation with deactivated B-field
  ConstantBField bField(Vector3D(0., 0., 0.));
  using RecoStepper = EigenStepper<ConstantBField>;
  RecoStepper rStepper(bField);
  using RecoPropagator = Propagator<RecoStepper, Navigator>;
  RecoPropagator rPropagator(rStepper, rNavigator);

  // Set initial parameters for the particle track
  Covariance cov;
  cov << 1000_um, 0., 0., 0., 0., 0., 0., 1000_um, 0., 0., 0., 0., 0., 0., 0.05,
      0., 0., 0., 0., 0., 0., 0.05, 0., 0., 0., 0., 0., 0., 0.01, 0., 0., 0.,
      0., 0., 0., 1.;

  Vector3D rPos(-3_m, 10_um * gauss(generator), 100_um * gauss(generator));
  Vector3D rMom(1_GeV, 0.025_GeV * gauss(generator),
                0.025_GeV * gauss(generator));

  SingleCurvilinearTrackParameters<ChargedPolicy> rStart(cov, rPos, rMom, 1.,
                                                         42.);

  const Surface* rSurface = &rStart.referenceSurface();

  using Updater = GainMatrixUpdater<BoundParameters>;
  using Smoother = GainMatrixSmoother<BoundParameters>;
  using KalmanFitter = KalmanFitter<RecoPropagator, Updater, Smoother>;

  KalmanFitter kFitter(rPropagator,
                       getDefaultLogger("KalmanFilter", Logging::VERBOSE));

  KalmanFitterOptions kfOptions(tgContext, mfContext, calContext, rSurface);

  // Fit the track
  auto fittedTrack = kFitter.fit(sourcelinks, rStart, kfOptions);
  auto fittedParameters = fittedTrack.fittedParameters.get();

  // Make sure it is deterministic
  auto fittedAgainTrack = kFitter.fit(sourcelinks, rStart, kfOptions);
  auto fittedAgainParameters = fittedAgainTrack.fittedParameters.get();

  CHECK_CLOSE_REL(fittedParameters.parameters().template head<5>(),
                  fittedAgainParameters.parameters().template head<5>(), 1e-5);
  CHECK_CLOSE_ABS(fittedParameters.parameters().template tail<1>(),
                  fittedAgainParameters.parameters().template tail<1>(), 1e-5);

  // Change the order of the sourcelinks
  std::vector<SourceLink> shuffledMeasurements = {
      sourcelinks[3], sourcelinks[2], sourcelinks[1],
      sourcelinks[4], sourcelinks[5], sourcelinks[0]};

  // Make sure it works for shuffled measurements as well
  auto fittedShuffledTrack =
      kFitter.fit(shuffledMeasurements, rStart, kfOptions);
  auto fittedShuffledParameters = fittedShuffledTrack.fittedParameters.get();

  CHECK_CLOSE_REL(fittedParameters.parameters().template head<5>(),
                  fittedShuffledParameters.parameters().template head<5>(),
                  1e-5);
  CHECK_CLOSE_ABS(fittedParameters.parameters().template tail<1>(),
                  fittedShuffledParameters.parameters().template tail<1>(),
                  1e-5);

  // Remove one measurement and find a hole
  std::vector<SourceLink> measurementsWithHole = {
      sourcelinks[0], sourcelinks[1], sourcelinks[2], sourcelinks[4],
      sourcelinks[5]};

  // Make sure it works for shuffled measurements as well
  auto fittedWithHoleTrack =
      kFitter.fit(measurementsWithHole, rStart, kfOptions);
  auto fittedWithHoleParameters = fittedWithHoleTrack.fittedParameters.get();

  // Count one hole
  BOOST_CHECK_EQUAL(fittedWithHoleTrack.missedActiveSurfaces.size(), 1);
  // And the parameters should be different
  //~
  // BOOST_CHECK(!Acts::Test::checkCloseRel(fittedParameters.parameters().template
  // head<5>(), ~ fittedWithHoleParameters.parameters().template head<5>(), ~
  // 1e-6));
  BOOST_CHECK(!Acts::Test::checkCloseRel(fittedParameters.parameters(),
                                         fittedWithHoleParameters.parameters(),
                                         1e-6));
}

}  // namespace Test
}  // namespace Acts
