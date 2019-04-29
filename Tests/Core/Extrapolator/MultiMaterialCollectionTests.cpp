// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
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

#include "Acts/Detector/TrackingGeometry.hpp"
#include "Acts/EventData/Measurement.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/TrackState.hpp"
#include "Acts/Extrapolator/Navigator.hpp"
#include "Acts/Extrapolator/SurfaceCollector.hpp"
#include "Acts/Fitter/GainMatrixSmoother.hpp"
#include "Acts/Fitter/GainMatrixUpdator.hpp"
#include "Acts/Fitter/KalmanFitter.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/MultiEigenStepper.hpp"
#include "Acts/Extrapolator/MultiMaterialInteractor.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "Acts/Propagator/detail/DebugOutputActor.hpp"
#include "Acts/Propagator/detail/StandardAborters.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Tests/CommonHelpers/CubicTrackingGeometry.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/GeometryID.hpp"
#include "Acts/Utilities/GeometryContext.hpp"
#include "Acts/Utilities/MagneticFieldContext.hpp"
#include "Acts/Utilities/CalibrationContext.hpp"

namespace Acts {
namespace Test {

  using DebugOutput = detail::DebugOutputActor;

  std::normal_distribution<double> gauss(0., 1.);
  std::default_random_engine       generator(42);

  bool debugMode = true;

  // Create a test context
  GeometryContext      tgContext  = GeometryContext();
  MagneticFieldContext mfContext  = MagneticFieldContext();
  CalibrationContext   calContext = CalibrationContext();

  ///
  /// @brief Unit test for Kalman fitter with measurements along the x-axis
  ///
  BOOST_AUTO_TEST_CASE(kalman_fitter_zero_field)
  {
    // Build detector
    CubicTrackingGeometry cGeometry(tgContext);
    auto                  detector = cGeometry();

    // Build navigator for the measurement creatoin
    Navigator mNavigator(detector);
    mNavigator.resolvePassive   = false;
    mNavigator.resolveMaterial  = true;
    mNavigator.resolveSensitive = true;

	/*
    // Use straingt line stepper to create the measurements
    StraightLineStepper mStepper;

    // Define the measurement propagator
    using MeasurementPropagator = Propagator<StraightLineStepper, Navigator>;

    // Build propagator for the measurement creation
    MeasurementPropagator mPropagator(mStepper, mNavigator);
    Vector3D mPos(-3. * units::_m, 0., 0.), mMom(1. * units::_GeV, 0., 0);
    SingleCurvilinearTrackParameters<NeutralPolicy> mStart(nullptr, mPos, mMom);

    // Set options for propagator
    PropagatorOptions<ActionList<DebugOutput> > mOptions(
        tgContext, mfContext);
    mOptions.debug              = debugMode;
    //mOptions.maxStepSize = 10. * units::_cm;
    //mOptions.pathLimit   = 25 * units::_cm;

    // Launch and collect - the measurements
    auto mResult = mPropagator.propagate(mStart, mOptions).value();
    if (debugMode) {
      const auto debugString
          = mResult.template get<DebugOutput::result_type>().debugString;
      std::cout << ">>>> Measurement creation: " << std::endl;
      std::cout << debugString;
    } 
	*/

    ConstantBField bField(Vector3D(0., 0., 0.));
    //using RecoStepper = EigenStepper<ConstantBField>;
    //RecoStepper rStepper(bField);
    //using RecoPropagator = Propagator<RecoStepper, Navigator>;
    //RecoPropagator rPropagator(rStepper, mNavigator);

    using RecoStepper = MultiEigenStepper<ConstantBField>;
    RecoStepper rStepper(bField);
    using RecoPropagator = Propagator<RecoStepper, Navigator>;
    RecoPropagator rPropagator(rStepper, mNavigator);
	//
    // Set initial parameters for the particle track
    ActsSymMatrixD<5> cov;
    cov << 1000. * units::_um, 0., 0., 0., 0., 0., 1000. * units::_um, 0., 0.,
        0., 0., 0., 0.05, 0., 0., 0., 0., 0., 0.05, 0., 0., 0., 0., 0., 0.01;

    auto covPtr = std::make_unique<const ActsSymMatrixD<5>>(cov);

    Vector3D rPos(-3. * units::_m,
                  10. * units::_um * gauss(generator),
                  100. * units::_um * gauss(generator));
    Vector3D rMom(1. * units::_GeV,
                  0.025 * units::_GeV * gauss(generator),
                  0.025 * units::_GeV * gauss(generator));

    SingleCurvilinearTrackParameters<ChargedPolicy> rStart(
        std::move(covPtr), rPos, rMom, 1.);

    PropagatorOptions<ActionList<DebugOutput,MultiMaterialInteractor>,AbortList<detail::EndOfWorldReached> > rOptions(
        tgContext, mfContext);
    rOptions.debug              = debugMode;

    // Launch and collect - the measurements
    auto mResult = rPropagator.propagate(rStart, rOptions).value();
	auto numOfComponents = mResult.template get<MultiMaterialInteractor::result_type>().numComponents;
    if (debugMode) {
      const auto debugString
          = mResult.template get<DebugOutput::result_type>().debugString;
      std::cout << ">>>> Measurement creation: " << std::endl;
      std::cout << debugString;
      std::cout << " There collects "<<numOfComponents<<" components.";
    } 

	// Test if the number of components split into 128 in the interactions of 6 surfaces
	BOOST_CHECK( numOfComponents == 128 );

	// Test for each surface all material interaction recorded are the same 
	// because the component split don't change anything
	const auto& material_interactions_result = mResult.template get<MultiMaterialInteractor::result_type>().multi_materialInteractions;
	for(const auto& materialInteractionPair: material_interactions_result)
	{
	  const auto& materialInteractionVec = materialInteractionPair.second;
	 BOOST_CHECK( std::all_of(materialInteractionVec.begin()+1,materialInteractionVec.end(),
		            [&](const Acts::MaterialInteractionVec::value_type& r) {return r == materialInteractionVec.front();}) );

	}

  }
}
}

