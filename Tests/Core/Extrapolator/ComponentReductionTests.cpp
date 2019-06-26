// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// clang-format off
#define BOOST_TEST_MODULE BoundParameters Tests
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
// clang-format on

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/StrawSurface.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Units.hpp"

#include "Acts/Extrapolator/detail/ComponentDistance.hpp"
#include "Acts/Extrapolator/detail/ComponentCombiner.hpp"
#include "Acts/Extrapolator/detail/component_reduction_impl.hpp"

namespace bdata = boost::unit_test::data;
namespace tt    = boost::test_tools;

namespace Acts {

namespace Test {

  // Create a test context
  GeometryContext tgContext = GeometryContext();

  /// @brief Unit test for parameters at a plane
  ///
  BOOST_DATA_TEST_CASE(
      component_reduction_test,
      bdata::random((bdata::seed = 1240,
                     bdata::distribution
                     = std::uniform_real_distribution<>(-1000., 1000.)))
          ^ bdata::random((bdata::seed = 2351,
                           bdata::distribution
                           = std::uniform_real_distribution<>(-1000., 1000.)))
          ^ bdata::random((bdata::seed = 3412,
                           bdata::distribution
                           = std::uniform_real_distribution<>(-1000., 1000.)))
          ^ bdata::random((bdata::seed = 5732,
                           bdata::distribution
                           = std::uniform_real_distribution<>(0., M_PI)))
          ^ bdata::random((bdata::seed = 8941,
                           bdata::distribution
                           = std::uniform_real_distribution<>(0., M_PI)))
          ^ bdata::random((bdata::seed = 1295,
                           bdata::distribution
                           = std::uniform_real_distribution<>(0., M_PI)))
          ^ bdata::xrange(1),
      x,
      y,
      z,
      a,
      b,
      c,
      index)
  {
    (void)index;
    Vector3D center{x, y, z};
    auto     transform = std::make_shared<Transform3D>();
    transform->setIdentity();
    RotationMatrix3D rot;
    rot = AngleAxis3D(a, Vector3D::UnitX()) * AngleAxis3D(b, Vector3D::UnitY())
        * AngleAxis3D(c, Vector3D::UnitZ());
    transform->prerotate(rot);
    transform->pretranslate(center);
    // create the surface
    auto bounds   = std::make_shared<RectangleBounds>(100., 100.);
    auto pSurface = Surface::makeShared<PlaneSurface>(
        transform, bounds);  // surface use +1

    // now create parameters on this surface
    // l_x, l_y, phi, theta, q/p (1/p)
    std::array<double, 5> pars_array_1
        = {{-0.1234, 9.8765, 0.45, 0.888, 0.001}};
    std::array<double, 5> pars_array_2
        = {{-0.1234, 9.8765, 0.45, 0.888, 0.002}};
    std::array<double, 5> pars_array_3
        = {{-0.1234, 9.8765, 0.45, 0.888, 0.004}};
    TrackParametersBase::ParVector_t pars1, pars2, pars3, pars4;
    pars1 << pars_array_1[0], pars_array_1[1], pars_array_1[2], pars_array_1[3],
        pars_array_1[4];
    pars2 << pars_array_2[0], pars_array_2[1], pars_array_2[2], pars_array_2[3],
        pars_array_2[4];
    pars3 << pars_array_3[0], pars_array_3[1], pars_array_3[2], pars_array_3[3],
        pars_array_3[4];

    ActsSymMatrixD<5> cov_1;
    ActsSymMatrixD<5> cov_2;
    ActsSymMatrixD<5> cov_3;
    cov_1 << 10 * units::_mm, 0, 0.123, 0, 0.5, 0, 10 * units::_mm, 0, 0.162, 0,
        0.123, 0, 0.1, 0, 0, 0, 0.162, 0, 0.1, 0, 0.5, 0, 0, 0,
        1. / (10 * units::_GeV);
    cov_2 = cov_1;
    cov_3 = cov_1;

    auto covPtr_1 = std::make_unique<const ActsSymMatrixD<5>>(cov_1);
    auto covPtr_2 = std::make_unique<const ActsSymMatrixD<5>>(cov_2);
    auto covPtr_3 = std::make_unique<const ActsSymMatrixD<5>>(cov_3);

    BoundParameters* ataPlane_from_pars_0
        = new BoundParameters(tgContext, std::move(covPtr_1), pars1, pSurface);
    BoundParameters* ataPlane_from_pars_1
        = new BoundParameters(tgContext, std::move(covPtr_2), pars2, pSurface);
    BoundParameters* ataPlane_from_pars_2
        = new BoundParameters(tgContext, std::move(covPtr_3), pars3, pSurface);

    // make multi bound par
    MultipleBoundParameters multi_ataPlane_from_pars(pSurface);
    multi_ataPlane_from_pars.append(0.4, ataPlane_from_pars_0);
    multi_ataPlane_from_pars.append(0.35, ataPlane_from_pars_1);
    multi_ataPlane_from_pars.append(0.25, ataPlane_from_pars_2);

    // get trackmap
    auto& trackMap = multi_ataPlane_from_pars.getTrackList();
    using TrackParMap =
        typename std::remove_reference<decltype(trackMap)>::type;
    // record the 2nd element
    typename TrackParMap::const_iterator element_1 = ++trackMap.begin();

    /// test1 : test component combination - 2 TrackParameters -> 1
    detail::ComponentCombiner combiner;
    auto                      combinedComponent
        = combiner(tgContext, *pSurface, *trackMap.begin(), *element_1);
    auto combinedParameters = 0.4 / 0.75 * ataPlane_from_pars_0->parameters()
        + 0.35 / 0.75 * ataPlane_from_pars_1->parameters();
    auto combinedCov = *ataPlane_from_pars_0->covariance() * 0.4 / 0.75
        + *ataPlane_from_pars_1->covariance() * 0.35 / 0.75;
    std::cout << "in combine component " << combinedComponent.first << ","
              << *combinedComponent.second << std::endl;
    std::cout << "in simple calculation par cov  " << combinedParameters << " "
              << combinedCov << std::endl;
    std::cout << "component 1: " << *ataPlane_from_pars_0 << std::endl;
    std::cout << "component 2: " << *ataPlane_from_pars_1 << std::endl;
    std::cout << std::endl;
    // check weight, parameters, covariance
    CHECK_CLOSE_REL(combinedComponent.first, 0.75, 1e-6);
    CHECK_CLOSE_REL(
        combinedComponent.second->parameters()[0], combinedParameters[0], 1e-6);
    CHECK_CLOSE_REL((*combinedComponent.second->covariance())(eLOC_0, eLOC_0),
                    combinedCov(eLOC_0, eLOC_0),
                    1e-6);
    CHECK_CLOSE_REL((*combinedComponent.second->covariance())(eLOC_1, eLOC_1),
                    combinedCov(eLOC_1, eLOC_1),
                    1e-6);
    CHECK_CLOSE_REL((*combinedComponent.second->covariance())(ePHI, ePHI),
                    combinedCov(ePHI, ePHI),
                    1e-6);
    CHECK_CLOSE_REL((*combinedComponent.second->covariance())(eTHETA, eTHETA),
                    combinedCov(eTHETA, eTHETA),
                    1e-6);
    CHECK_CLOSE_REL((*combinedComponent.second->covariance())(eQOP, eQOP),
                    combinedCov(eQOP, eQOP),
                    1e-6);

    /// test2 : test component reduction 3 TrackParameters -> 2
    // this test if the 1st and 2nd component are combined
    impl::reductComponent(trackMap, 2, tgContext, *pSurface);
    std::cout << "size " << trackMap.size() << std::endl;
    CHECK_CLOSE_REL(0.75, trackMap.begin()->first, 1e-6);
    //	CHECK_CLOSE_REL(combinedParameters,
    //trackMap.begin()->second->parameters(), 1e-6);
    std::cout << combinedParameters << std::endl;
    std::cout << std::endl;
    std::cout << trackMap.begin()->second->parameters() << std::endl;
  }
}
}
