// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// DiscSurface.h, Acts project
///////////////////////////////////////////////////////////////////

#pragma once

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryStatics.hpp"
#include "Acts/Surfaces/DiscBounds.hpp"
#include "Acts/Surfaces/InfiniteBounds.hpp"
#include "Acts/Surfaces/PolyhedronRepresentation.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Definitions.hpp"

namespace Acts {

class DetectorElementBase;

/// @class DiscSurface
///
/// Class for a DiscSurface in the, it inherits from Surface.
///
/// The DiscSurface has a polar local coordinate system, with
/// (r,phi) describing the coordinates.
///
/// The surface transform positions the disc such, that the origin
/// is at r=0, independent of the provided DiscBounds. The z-axis
/// The normal vector of the Disc, being perpendicular to the
/// radial direction.
///
/// The disc surface The only surface type for which the
/// covariance matrix is NOT given in the reference frame.
/// A conversion from polar to cartesian coordinates needs
/// to happen to transfer the local coordinates onto the
/// cartesian reference frame coordinates.
///
/// @image html DiscSurface.png
///
class DiscSurface : public Surface {
  friend Surface;

  DiscSurface() = delete;

 protected:
  /// Default Constructor is deleted

  /// Constructor for Discs from Transform3D, \f$ r_{min}, r_{max} \f$
  ///
  /// @param htrans is transform that places the disc in the global 3D space
  /// (can be nullptr)
  /// @param rmin The inner radius of the disc surface
  /// @param rmax The outer radius of the disc surface
  /// @param hphisec The opening angle of the disc surface and is optional
  ///        the default is a full disc
  DiscSurface(std::shared_ptr<const Transform3D> htrans, double rmin,
              double rmax, double hphisec = M_PI);

  /// Constructor for Discs from Transform3D, \f$ r_{min}, r_{max}, hx_{min},
  /// hx_{max} \f$
  /// This is n this case you have DiscTrapezoidalBounds
  ///
  /// @param htrans is transform that places the disc in the global 3D space
  /// (can be nullptr)
  /// @param minhalfx The half length in x at minimal r
  /// @param maxhalfx The half length in x at maximal r
  /// @param maxR The inner radius of the disc surface
  /// @param minR The outer radius of the disc surface
  /// @param avephi The position in phi (default is 0.)
  /// @param stereo The optional stereo angle
  DiscSurface(std::shared_ptr<const Transform3D> htrans, double minhalfx,
              double maxhalfx, double maxR, double minR, double avephi = 0.,
              double stereo = 0.);

  /// Constructor for Discs from Transform3D and shared DiscBounds
  ///
  /// @param htrans The transform that positions the disc in global 3D
  /// @param dbounds The disc bounds describing the surface coverage
  DiscSurface(std::shared_ptr<const Transform3D> htrans,
              std::shared_ptr<const DiscBounds> dbounds = nullptr);

  /// Constructor from DetectorElementBase : Element proxy
  ///
  /// @param dbounds The disc bounds describing the surface coverage
  /// @param detelement The detector element represented by this surface
  DiscSurface(const std::shared_ptr<const DiscBounds>& dbounds,
              const DetectorElementBase& detelement);

  /// Copy Constructor
  ///
  /// @param other The source surface for the copy
  DiscSurface(const DiscSurface& other);

  /// Copy constructor - with shift
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param other is the source cone surface
  /// @param transf is the additional transfrom applied after copying
  DiscSurface(const GeometryContext& gctx, const DiscSurface& other,
              const Transform3D& transf);

 public:
  /// Destructor - defaulted
  ~DiscSurface() override = default;

  /// Assignement operator
  ///
  /// @param other The source sourface for the assignment
  DiscSurface& operator=(const DiscSurface& other);

  /// Clone method into a concrete type of DiscSurface with shift
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param shift applied to the surface
  std::shared_ptr<DiscSurface> clone(const GeometryContext& gctx,
                                     const Transform3D& shift) const;

  /// Return the surface type
  SurfaceType type() const override;

  /// Normal vector return
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param lpos The local position is ignored
  ///
  /// @return a Vector3D by value
  const Vector3D normal(const GeometryContext& gctx,
                        const Vector2D& lpos) const final;

  /// Normal vector return without argument
  using Surface::normal;

  /// The binning position The position calcualted
  /// for a certain binning type
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param bValue The binning type to be used
  ///
  /// @return position that can beused for this binning
  const Vector3D binningPosition(const GeometryContext& gctx,
                                 BinningValue bValue) const final;

  /// This method returns the bounds by reference
  const SurfaceBounds& bounds() const final;

  /// Local to global transformation
  /// For planar surfaces the momentum is ignroed in the local to global
  /// transformation
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param lpos local 2D posittion in specialized surface frame
  /// @param mom global 3D momentum representation (optionally ignored)
  /// @param gpos global 3D position to be filled (given by reference for method
  /// symmetry)
  ///
  /// @note the momentum is ignored for Disc surfaces in this calculateion
  void localToGlobal(const GeometryContext& gctx, const Vector2D& lpos,
                     const Vector3D& mom, Vector3D& gpos) const final;

  /// Global to local transformation
  /// @note the momentum is ignored for Disc surfaces in this calculateion
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param gpos global 3D position - considered to be on surface but not
  /// inside bounds (check is done)
  /// @param mom global 3D momentum representation (optionally ignored)
  /// @param lpos local 2D position to be filled (given by reference for method
  /// symmetry)
  ///
  /// @return boolean indication if operation was successful (fail means global
  /// position was not on surface)
  bool globalToLocal(const GeometryContext& gctx, const Vector3D& gpos,
                     const Vector3D& mom, Vector2D& lpos) const final;

  /// Special method for DiscSurface : local<->local transformations polar <->
  /// cartesian
  ///
  /// @param lpolar is a local position in polar coordinates
  ///
  /// @return values is local 2D position in cartesian coordinates  @todo check
  const Vector2D localPolarToCartesian(const Vector2D& lpolar) const;

  /// Special method for Disc surface : local<->local transformations polar <->
  /// cartesian
  ///
  /// @param lcart is local 2D position in cartesian coordinates
  ///
  /// @return value is a local position in polar coordinates
  const Vector2D localCartesianToPolar(const Vector2D& lcart) const;

  /// Special method for DiscSurface : local<->local transformations polar <->
  /// cartesian
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param locpol is a local position in polar coordinates
  ///
  /// @return values is local 2D position in cartesian coordinates
  const Vector2D localPolarToLocalCartesian(const Vector2D& locpol) const;

  /// Special method for DiscSurface :  local<->global transformation when
  /// provided cartesian coordinates
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param lpos is local 2D position in cartesian coordinates
  ///
  /// @return value is a global cartesian 3D position
  const Vector3D localCartesianToGlobal(const GeometryContext& gctx,
                                        const Vector2D& lpos) const;

  /// Special method for DiscSurface : global<->local from cartesian coordinates
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param gpos is a global cartesian 3D position
  /// @param tol The absoltue tolerance parameter
  ///
  /// @return value is a local polar
  const Vector2D globalToLocalCartesian(const GeometryContext& gctx,
                                        const Vector3D& gpos,
                                        double tol = 0.) const;

  /// Initialize the jacobian from local to global
  /// the surface knows best, hence the calculation is done here.
  /// The jacobian is assumed to be initialised, so only the
  /// relevant entries are filled
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param jacobian The jacobian to be initialized
  /// @param gpos The global position of the parameters
  /// @param dir The direction at of the parameters
  ///
  /// @param pars The paranmeters vector
  void initJacobianToGlobal(const GeometryContext& gctx,
                            BoundToFreeMatrix& jacobian, const Vector3D& gpos,
                            const Vector3D& dir,
                            const BoundVector& pars) const final;

  /// Initialize the jacobian from global to local
  /// the surface knows best, hence the calculation is done here.
  /// The jacobian is assumed to be initialised, so only the
  /// relevant entries are filled
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param jacobian The jacobian to be initialized
  /// @param gpos The global position of the parameters
  /// @param dir The direction at of the parameters
  ///
  /// @return the transposed reference frame (avoids recalculation)
  const RotationMatrix3D initJacobianToLocal(const GeometryContext& gctx,
                                             FreeToBoundMatrix& jacobian,
                                             const Vector3D& gpos,
                                             const Vector3D& dir) const final;

  /// Path correction due to incident of the track
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param pos The global position as a starting point
  /// @param mom The global momentum at the starting point
  /// @return The correction factor due to incident
  double pathCorrection(const GeometryContext& gctx, const Vector3D& pos,
                        const Vector3D& mom) const final;

  /// @brief Straight line intersection schema
  ///
  /// navDir=anyDirection is to provide the closest solution
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param gpos The global position as a starting point
  /// @param gdir The global direction at the starting point
  ///        @note expected to be normalized (no checking)
  /// @param navDir is a navigation direction
  /// @param bcheck The boundary check prescription
  /// @param correct is a corrector function (e.g. for curvature correction)
  ///
  ///  <b>mathematical motivation:</b>
  ///
  /// the equation of the plane is given by: <br>
  /// @f$ \vec n \cdot \vec x = \vec n \cdot \vec p,@f$ <br>
  /// where @f$ \vec n = (n_{x}, n_{y}, n_{z})@f$ denotes the normal vector of
  /// the plane, @f$ \vec p = (p_{x}, p_{y}, p_{z})@f$ one specific point on
  /// the plane and @f$ \vec x = (x,y,z) @f$ all possible points
  /// on the plane.<br>
  /// Given a line with:<br>
  /// @f$ \vec l(u) = \vec l_{1} + u \cdot \vec v @f$, <br>
  /// the solution for @f$ u @f$ can be written:
  /// @f$ u = \frac{\vec n (\vec p - \vec l_{1})}{\vec n \vec v}@f$ <br>
  /// If the denominator is 0 then the line lies:
  /// - either in the plane
  /// - perpendicular to the normal of the plane
  ///
  /// @return is the surface intersection object
  Intersection intersectionEstimate(const GeometryContext& gctx,
                                    const Vector3D& gpos, const Vector3D& gdir,
                                    NavigationDirection navDir = forward,
                                    const BoundaryCheck& bcheck = false,
                                    CorrFnc correct = nullptr) const final;

  /// Return properly formatted class name for screen output
  std::string name() const override;

  /// Return a PolyhedronRepresentation for this object
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param l0div Number of divisions along l0 (phi)
  /// @param l1div Number of divisions along l1 (r)
  virtual PolyhedronRepresentation polyhedronRepresentation(
      const GeometryContext& gctx, size_t l0div = 10, size_t l1div = 1) const;

 protected:
  std::shared_ptr<const DiscBounds> m_bounds;  ///< bounds (shared)

 private:
  /// Clone method implementation
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param shift applied to the surface
  DiscSurface* clone_impl(const GeometryContext& gctx,
                          const Transform3D& shift) const override;
};

#include "Acts/Surfaces/detail/DiscSurface.ipp"

}  // end of namespace Acts
