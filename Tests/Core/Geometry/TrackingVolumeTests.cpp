// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// clang-format off
#define BOOST_TEST_MODULE TrackingVolume Tests
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
// clang-format on

#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Utilities/BoundingBox.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/Ray.hpp"
#include "Acts/Geometry/AbstractVolume.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/Volume.hpp"

#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Propagator/ActionList.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/detail/ConstrainedStep.hpp"
#include "Acts/Propagator/detail/DebugOutputActor.hpp"
#include "Acts/Propagator/detail/StandardAborters.hpp"
#include "Acts/Propagator/detail/SteppingLogger.hpp"

#include "Acts/Geometry/GeometryContext.hpp"

namespace Acts {
namespace Test {

BOOST_AUTO_TEST_SUITE(Geometry)

#define NBOXES 10
#define NTESTS 20
#include "BVHDataTestCase.hpp"

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Test
}  // namespace Acts
