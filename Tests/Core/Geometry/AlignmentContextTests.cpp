// This file is part of the Acts project.
//
// Copyright (C) 2017-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///  Boost include(s)
#define BOOST_TEST_MODULE AlignmentContext Tests
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <array>
#include <memory>

#include "Acts/Geometry/DetectorElementBase.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/PlanarBounds.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Units.hpp"

using namespace Acts::UnitLiterals;

namespace Acts {
namespace Test {

/// @class AlignmentContext
struct AlignmentContext {
  /// We have 2 different transforms
  std::shared_ptr<const std::array<Transform3D, 2>> alignmentStore = nullptr;

  /// Context index
  unsigned int alignmentIndex{0};

  /// Default contructor
  AlignmentContext() {}

  /// Constructor with Store and context index
  AlignmentContext(std::shared_ptr<const std::array<Transform3D, 2>> aStore,
                   unsigned int aIndex = 0)
      : alignmentStore(std::move(aStore)), alignmentIndex(aIndex) {}
};

/// @class AlignableDetectorElement
///
/// This is a lightweight type of detector element,
/// it simply implements the base class.
class AlignableDetectorElement : public DetectorElementBase {
 public:
  // Deleted default constructor
  AlignableDetectorElement() = delete;

  /// Constructor for single sided detector element
  /// - bound to a Plane Surface
  ///
  /// @param transform is the transform that element the layer in 3D frame
  /// @param pBounds is the planar bounds for the planar detector element
  /// @param thickness is the module thickness
  AlignableDetectorElement(std::shared_ptr<const Transform3D> transform,
                           std::shared_ptr<const PlanarBounds> pBounds,
                           double thickness)
      : DetectorElementBase(),
        m_elementTransform(std::move(transform)),
        m_elementThickness(thickness) {
    auto mutableSurface = Surface::makeShared<PlaneSurface>(pBounds, *this);
    m_elementSurface = mutableSurface;
  }

  ///  Destructor
  ~AlignableDetectorElement() override { /*nop */
  }

  /// Return local to global transform associated with this identifier
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  ///
  /// @note this is called from the surface().transform() in the PROXY mode
  const Transform3D& transform(const GeometryContext& gctx) const override;

  /// Return surface associated with this detector element
  const Surface& surface() const override;

  /// The maximal thickness of the detector element wrt normal axis
  double thickness() const override;

 private:
  /// the transform for positioning in 3D space
  std::shared_ptr<const Transform3D> m_elementTransform;
  /// the surface represented by it
  std::shared_ptr<const Surface> m_elementSurface{nullptr};
  /// the element thickness
  double m_elementThickness{0.};
};

inline const Transform3D& AlignableDetectorElement::transform(
    const GeometryContext& gctx) const {
  auto alignContext = std::any_cast<AlignmentContext>(gctx);
  if (alignContext.alignmentStore != nullptr and
      alignContext.alignmentIndex < 2) {
    return (*(alignContext.alignmentStore))[alignContext.alignmentIndex];
  }
  return (*m_elementTransform);
}

inline const Surface& AlignableDetectorElement::surface() const {
  return *m_elementSurface;
}

inline double AlignableDetectorElement::thickness() const {
  return m_elementThickness;
}

/// Unit test for creating compliant/non-compliant Surface object
BOOST_AUTO_TEST_CASE(AlignmentContextTests) {
  // The nominal and alingments
  Vector3D nominalCenter(0., 0., 0.);
  Vector3D negativeCenter(0., 0., -1.);
  Vector3D positiveCenter(0., 0., 1.);

  // Checkpoints
  Vector3D onNominal(3., 3., 0.);
  Vector3D onNegative(3., 3., -1.);
  Vector3D onPositive(3., 3., 1.);

  // Local position
  Vector2D localPosition(3., 3.);

  // A position place holder and dymmy momentum
  Vector3D globalPosition(100_cm, 100_cm, 100_cm);
  Vector3D dummyMomentum(4., 4., 4.);

  Transform3D negativeTransform = Transform3D::Identity();
  negativeTransform.translation() = negativeCenter;

  Transform3D positiveTransform = Transform3D::Identity();
  positiveTransform.translation() = positiveCenter;

  std::array<Transform3D, 2> alignmentArray = {negativeTransform,
                                               positiveTransform};

  std::shared_ptr<const std::array<Transform3D, 2>> alignmentStore =
      std::make_shared<const std::array<Transform3D, 2>>(alignmentArray);

  // The detector element at nominal position
  AlignableDetectorElement alignedElement(
      std::make_shared<const Transform3D>(Transform3D::Identity()),
      std::make_shared<const RectangleBounds>(100_cm, 100_cm), 1_mm);

  const auto& alignedSurface = alignedElement.surface();

  // The alignment centexts
  AlignmentContext defaultContext{};
  AlignmentContext negativeContext(alignmentStore, 0);
  AlignmentContext positiveContext(alignmentStore, 1);

  // Test the transforms
  BOOST_CHECK(alignedSurface.transform(defaultContext)
                  .isApprox(Transform3D::Identity()));
  BOOST_CHECK(
      alignedSurface.transform(negativeContext).isApprox(negativeTransform));
  BOOST_CHECK(
      alignedSurface.transform(positiveContext).isApprox(positiveTransform));

  // Test the centers
  BOOST_CHECK_EQUAL(alignedSurface.center(defaultContext), nominalCenter);
  BOOST_CHECK_EQUAL(alignedSurface.center(negativeContext), negativeCenter);
  BOOST_CHECK_EQUAL(alignedSurface.center(positiveContext), positiveCenter);

  // Test OnSurface
  BOOST_CHECK(
      alignedSurface.isOnSurface(defaultContext, onNominal, dummyMomentum));
  BOOST_CHECK(
      alignedSurface.isOnSurface(negativeContext, onNegative, dummyMomentum));
  BOOST_CHECK(
      alignedSurface.isOnSurface(positiveContext, onPositive, dummyMomentum));

  // Test local to Global and vice versa
  alignedSurface.localToGlobal(defaultContext, localPosition, dummyMomentum,
                               globalPosition);
  BOOST_CHECK_EQUAL(globalPosition, onNominal);
  alignedSurface.globalToLocal(defaultContext, onNominal, dummyMomentum,
                               localPosition);
  BOOST_CHECK_EQUAL(localPosition, Vector2D(3., 3.));

  alignedSurface.localToGlobal(negativeContext, localPosition, dummyMomentum,
                               globalPosition);
  BOOST_CHECK_EQUAL(globalPosition, onNegative);
  alignedSurface.globalToLocal(negativeContext, onNegative, dummyMomentum,
                               localPosition);
  BOOST_CHECK_EQUAL(localPosition, Vector2D(3., 3.));

  alignedSurface.localToGlobal(positiveContext, localPosition, dummyMomentum,
                               globalPosition);
  BOOST_CHECK_EQUAL(globalPosition, onPositive);
  alignedSurface.globalToLocal(positiveContext, onPositive, dummyMomentum,
                               localPosition);
  BOOST_CHECK_EQUAL(localPosition, Vector2D(3., 3.));
}

}  // namespace Test
}  // namespace Acts
