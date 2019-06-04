// This file is part of the Acts project.
//
// Copyright (C) 2016-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// Intersection.hpp, Acts project
///////////////////////////////////////////////////////////////////

#pragma once
#include <limits>
#include "Definitions.hpp"

namespace Acts {

/// A function typedef for the intersection correction
using CorrFnc = std::function<bool(Vector3D&, Vector3D&, double&)>;

/// @brief This enum describes the intersection status when attempting to
/// hit a surface. It still allows to define all Intersections as valid
/// except from missed and unreachable with a simple boolean cast,
/// but can be used to more accurately describe the sort of intersection.
///
/// The exact interpretation of the enum can only be done within the context
/// of the surface.intersectionEstimate(...) call, where position, direction,
/// overstep tolerance & boundary check description is given
enum class IntersectionStatus : int {
  overstepped = -1,
  missed = 0,
  unreachable = 0,
  reachable = 1,
  onSurface = 2
};

///  @struct Intersection
///
///  intersection struct used for position
struct Intersection {
  Vector3D position{0., 0., 0.};  ///< position of the intersection
  double pathLength{
      std::numeric_limits<double>::infinity()};  ///< path length to the
                                                 ///< intersection (if valid)
  double distance{
      std::numeric_limits<double>::infinity()};  ///< distance to boundary if
                                                 ///< not
  IntersectionStatus status{
      IntersectionStatus::missed};  ///< intersection status description

  /// Constructor with arguments
  ///
  /// @param sinter is the position of the intersection
  /// @param slength is the path length to the intersection
  /// @param sstatus describes the type of intersection, given the input
  /// @param dist is the distance to the closes surface boundary
  Intersection(const Vector3D& sinter, double slength,
               IntersectionStatus sstatus, double dist = 0.)
      : position(sinter),
        pathLength(slength),
        distance(dist),
        status(sstatus) {}

  /// Default constructor
  Intersection()
      : position(Vector3D(0., 0., 0.)),
        pathLength(std::numeric_limits<double>::infinity()) {}

  /// Bool() operator for validity checking
  explicit operator bool() const { return (bool)status; }

  /// Smaller operator for sorting,
  /// - it respects the validity of the intersection
  /// @param si is the intersection for testing
  bool operator<(const Intersection& si) const {
    if (status == IntersectionStatus::missed) {
      return false;
    }
    // now check the pathLength
    if (si) {
      return (pathLength < si.pathLength);
    }
    // the current path length wins
    return true;
  }

  /// Greater operator for sorting,
  /// - it respects the validity of the intersection
  /// @param si is the intersection for testing
  bool operator>(const Intersection& si) const {
    if (status == IntersectionStatus::missed) {
      return false;
    }
    // now check the pathLength
    if (si) {
      return (pathLength > si.pathLength);
    }
    // the current path length wins
    return true;
  }
};

/// class extensions to return also the object
template <typename object_t, typename representation_t = object_t>
class ObjectIntersection {
 public:
  Intersection intersection{};      ///< the intersection iself
  const object_t* object{nullptr};  ///< the object that was intersected
  const representation_t* representation{
      nullptr};  ///< the representation of the object

  /// Default constructor
  ObjectIntersection() = default;

  /// Object intersection
  ///
  /// @param sInter is the intersection
  /// @param sObject is the object to be instersected
  /// @param dir is the direction of the intersection
  template <typename O = object_t, typename R = representation_t,
            typename = std::enable_if_t<std::is_same_v<O, R>>>
  ObjectIntersection(const Intersection& sInter, const object_t* sObject)
      : intersection(sInter), object(sObject), representation(sObject) {}

  /// Object intersection constructor with different representation type
  ///
  /// @param sInter is the intersection struct
  /// @param sObject is the intersected object
  /// @param sRepresentation is the surface representation of the object
  /// @param dir is the direction
  ///
  ObjectIntersection(const Intersection& sInter, const object_t* sObject,
                     const representation_t* sRepresentation)
      : intersection(sInter),
        object(sObject),
        representation(sRepresentation) {}

  /// Bool() operator for validity checking
  explicit operator bool() const { return (bool)intersection; }

  /// @brief smaller operator for ordering & sorting
  ///
  /// @param oi is the source intersection for comparison
  bool operator<(
      const ObjectIntersection<object_t, representation_t>& oi) const {
    return (intersection < oi.intersection);
  }

  /// @brief greater operator for ordering & sorting
  ///
  /// @param oi is the source intersection for comparison
  bool operator>(
      const ObjectIntersection<object_t, representation_t>& oi) const {
    return (intersection > oi.intersection);
  }
};

/// @brief Void Direction corrector
///
/// This is used to evaluate a modified
/// intersection (e.g. curvature updated)
struct VoidIntersectionCorrector {
  // Void Corrector default constructor
  VoidIntersectionCorrector() = default;

  // Void Corrector parameter constructor
  VoidIntersectionCorrector(const Vector3D& /*unused*/,
                            const Vector3D& /*unused*/, double /*unused*/) {}

  /// Boolean() operator - returns false for void modifier
  explicit operator bool() const { return false; }

  /// empty correction interface
  bool operator()(Vector3D& /*unused*/, Vector3D& /*unused*/,
                  double /*unused*/) {
    return false;
  }

  /// Step modification call
  ///
  /// @stay put and don't do antyhing
  template <typename step_t>
  bool operator()(step_t& /*unused*/) const {
    return false;
  }
};
}  // namespace Acts
