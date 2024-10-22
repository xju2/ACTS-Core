// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// TrackingVolume.h, Acts project
///////////////////////////////////////////////////////////////////

#pragma once
#include <functional>
#include <map>
#include <string>
#include "Acts/Geometry/BoundarySurfaceT.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryID.hpp"
#include "Acts/Geometry/GeometrySignature.hpp"
#include "Acts/Geometry/Layer.hpp"
#include "Acts/Geometry/Volume.hpp"
#include "Acts/Geometry/detail/BoundaryIntersectionSorter.hpp"
#include "Acts/Material/IVolumeMaterial.hpp"
#include "Acts/Surfaces/BoundaryCheck.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/BinnedArray.hpp"
#include "Acts/Utilities/BoundingBox.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Frustum.hpp"
#include "Acts/Utilities/Ray.hpp"

namespace Acts {

class GlueVolumesDescriptor;
class VolumeBounds;

// master typedefs
using TrackingVolumePtr = std::shared_ptr<const TrackingVolume>;
using MutableTrackingVolumePtr = std::shared_ptr<TrackingVolume>;

using TrackingVolumeBoundaryPtr =
    std::shared_ptr<const BoundarySurfaceT<TrackingVolume>>;

// possible contained
using TrackingVolumeArray = BinnedArray<TrackingVolumePtr>;
using TrackingVolumeVector = std::vector<TrackingVolumePtr>;
using LayerArray = BinnedArray<LayerPtr>;
using LayerVector = std::vector<LayerPtr>;

// full intersection with Layer
using LayerIntersection = FullIntersection<Layer, Surface>;

// full intersection with surface
using BoundaryIntersection =
    FullIntersection<BoundarySurfaceT<TrackingVolume>, Surface>;

/// @class TrackingVolume
///
/// Full Volume description used in Tracking,
/// it inherits from Volume to get the geometrical structure.
///
///     A TrackingVolume at navigation level can provide the (layer) material
/// information / internal navigation with in
///     5 different ways:
///
///         --- a) Static confinement of Layers
///         --- b) detached sub volumes
///         --- b) unordered (arbitrarily oriented) layers
///         --- d) unordered sub volumes
///         --- e) unordered layers AND unordered subvolumes
///
///    The TrackingVolume can also be a simple container of other
///    TrackingVolumes
///
/// In addition it is capable of holding a subarray of Layers and
/// TrackingVolumes.
///
class TrackingVolume : public Volume {
  friend class TrackingGeometry;

 public:
  /// Destructor
  ~TrackingVolume() override;

  /// Factory constructor for a container TrackingVolume
  /// - by definition a Vacuum volume
  ///
  /// @param htrans is the global 3D transform to position the volume in space
  /// @param volumeBounds is the description of the volume boundaries
  /// @param containedVolumes are the static volumes that fill this volume
  /// @param volumeName is a string identifier
  ///
  /// @return shared pointer to a new TrackingVolume
  static MutableTrackingVolumePtr create(
      std::shared_ptr<const Transform3D> htrans, VolumeBoundsPtr volumeBounds,
      const std::shared_ptr<const TrackingVolumeArray>& containedVolumes =
          nullptr,
      const std::string& volumeName = "undefined") {
    return MutableTrackingVolumePtr(
        new TrackingVolume(std::move(htrans), std::move(volumeBounds),
                           containedVolumes, volumeName));
  }

  /// Factory constructor for Tracking Volume with a bounding volume hierarchy
  ///
  /// @param htrans is the global 3D transform to position the volume in space
  /// @param volBounds is the description of the volume boundaries
  /// @param boxStore Vector owning the contained bounding boxes
  /// @param descendants Vector owning the child volumes
  /// @param top The top of the hierarchy (top node)
  /// @param matprop is are materials of the tracking volume
  /// @param volumeName is a string identifier
  ///
  /// @return shared pointer to a new TrackingVolume
  static MutableTrackingVolumePtr create(
      std::shared_ptr<const Transform3D> htrans, VolumeBoundsPtr volbounds,
      std::vector<std::unique_ptr<Volume::BoundingBox>> boxStore,
      std::vector<std::unique_ptr<const Volume>> descendants,
      const Volume::BoundingBox* top,
      std::shared_ptr<const IVolumeMaterial> volumeMaterial,
      const std::string& volumeName = "undefined") {
    return MutableTrackingVolumePtr(new TrackingVolume(
        std::move(htrans), std::move(volbounds), std::move(boxStore),
        std::move(descendants), top, std::move(volumeMaterial), volumeName));
  }

  /// Factory constructor for Tracking Volumes with content
  /// - can not be a container volume
  ///
  /// @param htrans is the global 3D transform to position the volume in space
  /// @param volumeBounds is the description of the volume boundaries
  /// @param matprop is are materials of the tracking volume
  /// @param containedLayers is the confined layer array (optional)
  /// @param containedVolumes is the confined volume array (optional)
  /// @param volumeName is a string identifier
  ///
  /// @return shared pointer to a new TrackingVolume
  static MutableTrackingVolumePtr create(
      std::shared_ptr<const Transform3D> htrans, VolumeBoundsPtr volumeBounds,
      std::shared_ptr<const IVolumeMaterial> volumeMaterial,
      std::unique_ptr<const LayerArray> containedLayers = nullptr,
      std::shared_ptr<const TrackingVolumeArray> containedVolumes = nullptr,
      const std::string& volumeName = "undefined") {
    return MutableTrackingVolumePtr(new TrackingVolume(
        std::move(htrans), std::move(volumeBounds), std::move(volumeMaterial),
        std::move(containedLayers), std::move(containedVolumes), volumeName));
  }

  /// Return the associated Layer to the global position
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param gp is the associated global position
  ///
  /// @return plain pointer to layer object
  const Layer* associatedLayer(const GeometryContext& gctx,
                               const Vector3D& gp) const;

  /// @brief Resolves the volume into (compatible) Layers
  ///
  /// This is the method for the propagator/extrapolator
  /// @tparam options_t Type of navigation options object for decomposition
  /// @tparam corrector_t Type of (optional) corrector for surface intersection
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param position Position for the search
  /// @param direction Direction for the search
  /// @param options The templated navigation options
  /// @param corrfnc is the corrector struct / function
  ///
  /// @return vector of compatible intersections with layers
  template <typename options_t,
            typename corrector_t = VoidIntersectionCorrector>
  std::vector<LayerIntersection> compatibleLayers(
      const GeometryContext& gctx, const Vector3D& position,
      const Vector3D& direction, const options_t& options,
      const corrector_t& corrfnc = corrector_t()) const;

  /// @brief Resolves the volume into (compatible) Layers
  ///
  /// This is the method for the propagator/extrapolator
  /// @tparam parameters_t Type of parameters used for the decomposition
  /// @tparam options_t Type of navigation options object for decomposition
  /// @tparam corrector_t Type of (optional) corrector for surface intersection
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param parameters The templated parameters for searching
  /// @param options The templated navigation options
  /// @param corrfnc is the corrector struct / function
  ///
  /// @return vector of compatible intersections with layers
  template <typename parameters_t, typename options_t,
            typename corrector_t = VoidIntersectionCorrector>
  std::vector<LayerIntersection> compatibleLayers(
      const GeometryContext& gctx, const parameters_t& parameters,
      const options_t& options,
      const corrector_t& corrfnc = corrector_t()) const;

  /// @brief Returns all boundary surfaces sorted by the user.
  ///
  /// @tparam options_t Type of navigation options object for decomposition
  /// @tparam corrector_t Type of (optional) corrector for surface intersection
  /// @tparam sorter_t Type of the boundary surface sorter
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param position The position for searching
  /// @param direction The direction for searching
  /// @param options The templated navigation options
  /// @param corrfnc is the corrector struct / function
  /// @param sorter Sorter of the boundary surfaces
  ///
  /// @return is the templated boundary intersection
  template <typename options_t,
            typename corrector_t = VoidIntersectionCorrector,
            typename sorter_t = DefaultBoundaryIntersectionSorter>
  std::vector<BoundaryIntersection> compatibleBoundaries(
      const GeometryContext& gctx, const Vector3D& position,
      const Vector3D& direction, const options_t& options,
      const corrector_t& corrfnc = corrector_t(),
      const sorter_t& sorter = sorter_t()) const;

  /// @brief Returns all boundary surfaces sorted by the user.
  ///
  /// @tparam parameters_t Type of parameters used for the decomposition
  /// @tparam options_t Type of navigation options object for decomposition
  /// @tparam corrector_t Type of (optional) corrector for surface intersection
  /// @tparam sorter_t Type of the boundary surface sorter
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param parameters The templated parameters for searching
  /// @param options The templated navigation options
  /// @param corrfnc is the corrector struct / function
  /// @param sorter Sorter of the boundary surfaces
  ///
  /// @return is the templated boundary intersection
  template <typename parameters_t, typename options_t,
            typename corrector_t = VoidIntersectionCorrector,
            typename sorter_t = DefaultBoundaryIntersectionSorter>
  std::vector<BoundaryIntersection> compatibleBoundaries(
      const GeometryContext& gctx, const parameters_t& parameters,
      const options_t& options, const corrector_t& corrfnc = corrector_t(),
      const sorter_t& sorter = sorter_t()) const;

  /// @brief Return surfaces in given direction from bounding volume hierarchy
  /// @tparam options_t Type of navigation options object for decomposition
  /// @tparam corrector_t Type of (optional) corrector for surface intersection
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param position The position to start from
  /// @param direction The direction towards which to test
  /// @param angle The opening angle
  /// @param options The templated navigation options
  /// @param corrfnc is the corrector struct / function
  ///
  /// @return Vector of surface candidates
  template <typename options_t,
            typename corrector_t = VoidIntersectionCorrector>
  std::vector<SurfaceIntersection> compatibleSurfacesFromHierarchy(
      const GeometryContext& gctx, const Vector3D& position,
      const Vector3D& direction, double angle, const options_t& options,
      const corrector_t& corrfnc = corrector_t()) const;

  /// Return the associated sub Volume, returns THIS if no subVolume exists
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param gp is the global position associated with that search
  ///
  /// @return plain pointer to associated with the position
  const TrackingVolume* lowestTrackingVolume(const GeometryContext& gctx,
                                             const Vector3D& gp) const;

  /// Return the confined static layer array - if it exists
  /// @return the BinnedArray of static layers if exists
  const LayerArray* confinedLayers() const;

  /// Return the confined volumes of this container array - if it exists
  std::shared_ptr<const TrackingVolumeArray> confinedVolumes() const;

  /// @brief Visit all sensitive surfaces
  ///
  /// @param visitor The callable. Will be called for each sensitive surface
  /// that is found
  ///
  /// If a context is needed for the vist, the vistitor has to provide this
  /// e.g. as a private member
  void visitSurfaces(
      const std::function<void(const Acts::Surface*)>& visitor) const;

  /// Returns the VolumeName - for debug reason, might be depreciated later
  const std::string& volumeName() const;

  /// Method to return the BoundarySurfaces
  const std::vector<TrackingVolumeBoundaryPtr>& boundarySurfaces() const;

  /// Return the material of the volume
  const IVolumeMaterial* volumeMaterial() const;

  /// Return the material of the volume as shared pointer
  const std::shared_ptr<const IVolumeMaterial>& volumeMaterialSharedPtr() const;

  /// Set the volume material description
  ///
  /// The material is usually derived in a complicated way and loaded from
  /// a framework given source. As various volumes could potentially share the
  /// the same material description, it is provided as a shared object
  ///
  /// @param material Material description of this volume
  void assignVolumeMaterial(std::shared_ptr<const IVolumeMaterial> material);

  /// Glue another tracking volume to this one
  ///  - if common face is set the glued volumes are sharing the boundary, down
  /// to the last navigation volume
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param bsfMine is the boundary face indicater where to glue
  /// @param neighbor is the TrackingVolume to be glued
  /// @param bsfNeighbor is the boudnary surface of the neighbor
  void glueTrackingVolume(const GeometryContext& gctx,
                          BoundarySurfaceFace bsfMine,
                          const MutableTrackingVolumePtr& neighbor,
                          BoundarySurfaceFace bsfNeighbor);

  /// Glue another tracking volume to this one
  ///  - if common face is set the glued volumes are sharing the boundary, down
  /// to the last navigation volume
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param bsfMine is the boundary face indicater where to glue
  /// @param neighbors are the TrackingVolumes to be glued
  /// @param bsfNeighbor are the boudnary surface of the neighbors
  void glueTrackingVolumes(
      const GeometryContext& gctx, BoundarySurfaceFace bsfMine,
      const std::shared_ptr<TrackingVolumeArray>& neighbors,
      BoundarySurfaceFace bsfNeighbor);

  /// Provide a new BoundarySurface from the glueing
  ///
  ///
  /// @param bsf is the boundary face indicater where to glue
  /// @param bs is the new boudnary surface
  void updateBoundarySurface(
      BoundarySurfaceFace bsf,
      std::shared_ptr<const BoundarySurfaceT<TrackingVolume>> bs);

  /// Register the outside glue volumes -
  /// ordering is in the TrackingVolume Frame:
  ///  - negativeFaceXY
  ///  - (faces YZ, ZY, radial faces)
  ///  - positiveFaceXY
  ///
  /// @param gvd register a new GlueVolumeDescriptor
  /// @todo update to shared/unique ptr
  void registerGlueVolumeDescriptor(GlueVolumesDescriptor* gvd);

  /// Register the outside glue volumes -
  /// ordering is in the TrackingVolume Frame:
  ///  - negativeFaceXY
  ///  - (faces YZ, ZY, radial faces)
  ///  - positiveFaceXY
  GlueVolumesDescriptor& glueVolumesDescriptor();

  /// Return whether this TrackingVolume has a BoundingVolumeHierarchy
  /// associated
  /// @return If it has a BVH or not.
  bool hasBoundingVolumeHierarchy() const;

  /// Sign the volume - the geometry builder has to do that
  ///
  /// @param geosign is the volume signature
  /// @param geotype is the volume navigation type
  void sign(GeometrySignature geosign, GeometryType geotype = Static);

  /// return the Signature
  GeometrySignature geometrySignature() const;

  /// return the Signature
  GeometryType geometryType() const;

  /// Register the color code
  ///
  /// @param icolor is a color number
  void registerColorCode(unsigned int icolor);

  /// Get the color code
  unsigned int colorCode() const;

  /// Return the MotherVolume - if it exists
  const TrackingVolume* motherVolume() const;

  /// Set the MotherVolume
  ///
  /// @param mvol is the mother volume
  void setMotherVolume(const TrackingVolume* mvol);

 protected:
  /// Default constructor
  TrackingVolume();

  /// Constructor for a container Volume
  /// - vacuum filled volume either as a for other tracking volumes
  ///
  /// @param htrans is the global 3D transform to position the volume in space
  /// @param volbounds is the description of the volume boundaries
  /// @param containedVolumeArray are the static volumes that fill this volume
  /// @param volumeName is a string identifier
  TrackingVolume(std::shared_ptr<const Transform3D> htrans,
                 VolumeBoundsPtr volbounds,
                 const std::shared_ptr<const TrackingVolumeArray>&
                     containedVolumeArray = nullptr,
                 const std::string& volumeName = "undefined");

  TrackingVolume(std::shared_ptr<const Transform3D> htrans,
                 VolumeBoundsPtr volbounds,
                 std::vector<std::unique_ptr<Volume::BoundingBox>> boxStore,
                 std::vector<std::unique_ptr<const Volume>> descendants,
                 const Volume::BoundingBox* top,
                 std::shared_ptr<const IVolumeMaterial> volumeMaterial,
                 const std::string& volumeName = "undefined");

  /// Constructor for a full equipped Tracking Volume
  ///
  /// @param htrans is the global 3D transform to position the volume in space
  /// @param volumeBounds is the description of the volume boundaries
  /// @param volumeMaterial is are materials of the tracking volume
  /// @param staticLayerArray is the confined layer array (optional)
  /// @param containedVolumeArray is the confined volume array
  /// @param volumeName is a string identifier
  TrackingVolume(
      std::shared_ptr<const Transform3D> htrans, VolumeBoundsPtr volumeBounds,
      std::shared_ptr<const IVolumeMaterial> volumeMaterial,
      std::unique_ptr<const LayerArray> staticLayerArray = nullptr,
      std::shared_ptr<const TrackingVolumeArray> containedVolumeArray = nullptr,
      const std::string& volumeName = "undefined");

 private:
  /// Create Boundary Surface
  void createBoundarySurfaces();

  /// method to synchronize the layers with potentially updated volume bounds:
  /// - adapts the layer dimensions to the new volumebounds + envelope
  ///
  /// @param envelope is the clearance between volume boundary and layer
  void synchronizeLayers(double envelope = 1.) const;

  /// close the Geometry, i.e. set the GeometryID and assign material
  ///
  /// @param materialDecorator is a dedicated decorator for the
  ///        material to be assigned (surface, volume based)
  /// @param volumeMap is a map to find the a volume
  ///        by a given name
  /// @param vol is the geometry id of the volume
  ///        as calculated by the TrackingGeometry
  ///
  void closeGeometry(const IMaterialDecorator* materialDecorator,
                     std::map<std::string, const TrackingVolume*>& volumeMap,
                     size_t& vol);

  /// interlink the layers in this TrackingVolume
  void interlinkLayers();

  template <typename T>
  static std::vector<const Volume*> intersectSearchHierarchy(
      const T obj, const Volume::BoundingBox* lnode);

  /// Forbidden copy constructor - deleted
  TrackingVolume(const TrackingVolume&) = delete;

  /// Forbidden assignment - deleted
  TrackingVolume& operator=(const TrackingVolume&) = delete;

  /// The volume based material the TrackingVolume consists of
  std::shared_ptr<const IVolumeMaterial> m_volumeMaterial{nullptr};

  /// Remember the mother volume
  const TrackingVolume* m_motherVolume{nullptr};

  // the boundary surfaces
  std::vector<TrackingVolumeBoundaryPtr> m_boundarySurfaces;

  ///(a) static configuration ordered by Binned arrays
  /// static layers
  std::unique_ptr<const LayerArray> m_confinedLayers = nullptr;

  /// Array of Volumes inside the Volume when actin as container
  std::shared_ptr<const TrackingVolumeArray> m_confinedVolumes = nullptr;

  /// Volumes to glue Volumes from the outside
  GlueVolumesDescriptor* m_glueVolumeDescriptor{nullptr};

  /// The Signature done by the GeometryBuilder
  GeometrySignature m_geometrySignature{Unsigned};

  /// The gometry type for the navigation schema
  GeometryType m_geometryType{NumberOfGeometryTypes};

  //// Volume name for debug reasons & screen output
  std::string m_name;

  /// color code for displaying
  unsigned int m_colorCode{20};

  /// Bounding Volume Hierarchy (BVH)
  std::vector<std::unique_ptr<const Volume::BoundingBox>> m_boundingBoxes;
  std::vector<std::unique_ptr<const Volume>> m_descendantVolumes;
  const Volume::BoundingBox* m_bvhTop{nullptr};
};

inline const std::string& TrackingVolume::volumeName() const {
  return m_name;
}

inline const IVolumeMaterial* TrackingVolume::volumeMaterial() const {
  return m_volumeMaterial.get();
}

inline const std::shared_ptr<const IVolumeMaterial>&
TrackingVolume::volumeMaterialSharedPtr() const {
  return m_volumeMaterial;
}

inline void TrackingVolume::assignVolumeMaterial(
    std::shared_ptr<const IVolumeMaterial> material) {
  m_volumeMaterial = std::move(material);
}

inline const LayerArray* TrackingVolume::confinedLayers() const {
  return m_confinedLayers.get();
}

inline std::shared_ptr<const TrackingVolumeArray>
TrackingVolume::confinedVolumes() const {
  return m_confinedVolumes;
}

inline GeometrySignature TrackingVolume::geometrySignature() const {
  return m_geometrySignature;
}

inline GeometryType TrackingVolume::geometryType() const {
  return m_geometryType;
}

inline void TrackingVolume::registerColorCode(unsigned int icolor) {
  m_colorCode = icolor;
}

inline unsigned int TrackingVolume::colorCode() const {
  return m_colorCode;
}

inline const TrackingVolume* TrackingVolume::motherVolume() const {
  return m_motherVolume;
}

inline void TrackingVolume::setMotherVolume(const TrackingVolume* mvol) {
  m_motherVolume = mvol;
}

inline bool TrackingVolume::hasBoundingVolumeHierarchy() const {
  return m_bvhTop != nullptr;
}

#include "detail/TrackingVolume.ipp"

}  // namespace Acts
