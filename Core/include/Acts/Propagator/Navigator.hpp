// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <boost/algorithm/string.hpp>
#include <iomanip>
#include <iterator>
#include <sstream>
#include <string>
#include "Acts/Geometry/BoundarySurfaceT.hpp"
#include "Acts/Geometry/Layer.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Geometry/detail/BoundaryIntersectionSorter.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/detail/ConstrainedStep.hpp"
#include "Acts/Surfaces/Surface.hpp"

namespace Acts {

using Cstep = detail::ConstrainedStep;

/// Navigator class
///
/// This class feeds the propagator with surfaces in order to handle
/// material integration and external surface support.
///
/// The propagator will call at initialization and each step in sequence
///
/// navigator.status()
/// actionList()
/// navigator.target()
///
/// As the actionList can modify the state itself, even if status can
/// already set a target, a retargetting has to be done, in the target()
/// call.
class Navigator {
 public:
  /// Parameter wrapper struct for the stepper and state, this ensures
  /// that the necessary navigation calls respect the stepper internal
  /// representation, e.g. for multi component stepping
  ///
  /// @tparam the stepper for it's functionality
  /// @tparam the state of the stepper
  template <typename stepper_t>
  struct Parameters {
    const stepper_t& stepper;
    const typename stepper_t::State& state;

    /// There is no default parameters
    Parameters() = delete;

    /// Constructor with references
    Parameters(const stepper_t& istepper,
               const typename stepper_t::State& istate)
        : stepper(istepper), state(istate) {}

    /// Forward the position calculation to the stepper
    const Vector3D position() const { return stepper.position(state); }

    /// Forward the direction calculation to the stepper
    const Vector3D direction() const { return stepper.direction(state); }
  };

  /// @brief struct for the Navigation options that are forwarded to
  ///        the geometry
  ///
  /// @tparam propagator_state_t Type of the object for navigation state
  /// @tparam object_t Type of the object for navigation to check against
  template <typename object_t>
  struct Options {
    /// The navigation direction
    NavigationDirection navDir = forward;

    /// The boundary check directive
    BoundaryCheck boundaryCheck = true;

    // How to resolve the geometry
    /// Always look for sensitive
    bool resolveSensitive = true;
    /// Always look for material
    bool resolveMaterial = true;
    /// always look for passive
    bool resolvePassive = false;

    /// object to check against: at start
    const object_t* startObject = nullptr;
    /// object to check against: at end
    const object_t* endObject = nullptr;

    /// Target surface to exclude
    const Surface* targetSurface = nullptr;

    /// Maximum limit for this navigaiton check
    double pathLimit = std::numeric_limits<double>::max();

    /// Overstepping limit for navigation actions
    double overstepLimit = s_onSurfaceTolerance;

    /// Constructor
    ///
    /// @param nDir Navigation direction prescription
    /// @param bcheck Boundary check for the navigation action
    /// @param resolves directive to resolve sensitive surfaces
    /// @param resolve, directive to resolve material surfaces
    /// @param resolvep directive to resolve passive surfaces
    /// @param sobject Start object to check against
    /// @param eobject End object to check against
    Options(NavigationDirection ndir, BoundaryCheck bcheck,
            bool resolves = true, bool resolvem = true, bool resolvep = false,
            const object_t* sobject = nullptr,
            const object_t* eobject = nullptr)
        : navDir(ndir),
          boundaryCheck(std::move(bcheck)),
          resolveSensitive(resolves),
          resolveMaterial(resolvem),
          resolvePassive(resolvep),
          startObject(sobject),
          endObject(eobject) {}
  };

  using NavigationSurfaces = std::vector<SurfaceIntersection>;
  using NavigationSurfaceIter = NavigationSurfaces::iterator;

  using NavigationLayers = std::vector<LayerIntersection>;
  using NavigationLayerIter = NavigationLayers::iterator;

  using NavigationBoundaries = std::vector<BoundaryIntersection>;
  using NavigationBoundaryIter = NavigationBoundaries::iterator;

  /// Constructor with shared tracking geometry
  ///
  /// @param tGeometry The tracking geometry for the navigator
  Navigator(std::shared_ptr<const TrackingGeometry> tGeometry = nullptr)
      : trackingGeometry(std::move(tGeometry)) {}

  /// Tracking Geometry for this Navigator
  std::shared_ptr<const TrackingGeometry> trackingGeometry;

  /// The tolerance used to defined "reached"
  double tolerance = s_onSurfaceTolerance;

  /// Configuration for this Navigator
  /// stop at every sensitive surface (whether it has material or not)
  bool resolveSensitive = true;
  /// stop at every material surface (whether it is passive or not)
  bool resolveMaterial = true;
  /// stop at every surface regardless what it is
  bool resolvePassive = false;

  /// Nested State struct
  ///
  /// It acts as an internal state which is
  /// created for every propagation/extrapolation step
  /// and keep thread-local navigation information
  struct State {
    // External surfaces, set and given by the user
    /// the vector of external surfaces to work through
    NavigationSurfaces extSurfaces = {};
    /// the current surface iterator of the navigation state
    NavigationSurfaceIter extSurfaceIter = extSurfaces.end();

    // Navigation on surface level
    /// the vector of navigation surfaces to work through
    NavigationSurfaces navSurfaces = {};
    /// the current surface iterator of the navigation state
    NavigationSurfaceIter navSurfaceIter = navSurfaces.end();

    // Navigation on layer level
    /// the vector of navigation layers to work through
    NavigationLayers navLayers = {};
    /// the current layer iterator of the navigation state
    NavigationLayerIter navLayerIter = navLayers.end();

    // Navigation on volume level
    /// the vector of boundary surfaces to work through
    NavigationBoundaries navBoundaries = {};
    /// the current boundary iterator of the navigation state
    NavigationBoundaryIter navBoundaryIter = navBoundaries.end();

    /// Navigation state: the start volume
    const TrackingVolume* startVolume = nullptr;
    /// Navigation state: the start layer
    const Layer* startLayer = nullptr;
    /// Navigation state: the start surface
    const Surface* startSurface = nullptr;
    /// Navigation state - external state: the current surface
    const Surface* currentSurface = nullptr;
    /// Navigation state: the current volume
    const TrackingVolume* currentVolume = nullptr;
    /// Navigation state: the target surface
    const Surface* targetSurface = nullptr;

    /// Indicator for start layer treatment
    bool startLayerResolved = false;
    /// Indicator if the target is reached
    bool targetReached = false;
    /// Navigation state : a break has been detected
    bool navigationBreak = false;

    // The next target
    struct {
      /// The current target intersection
      const Surface* surface = nullptr;
      /// The pathlength to next
      double distance = std::numeric_limits<double>::max();
      /// The target boundary check for updating
      BoundaryCheck boundaryCheck = true;

    } nextTarget;
  };

  /// @brief Navigator status call, will be called in two modes
  /// --------------------------------------------------------------------
  ///
  /// (a) It initializes the Navigation stream if start volume is
  ///     not yet defined:
  ///  - initialize the volume
  ///  - establish the start layer and start volume
  ///  - set the current surface to the start surface
  ///
  /// (b) It establishes the currentSurface status during
  ///     the propagation flow, currentSurface can be
  ///  - surfaces still to be handled within a layer
  ///  - layers still to be handled within a volume
  ///  - boundaries still to be handled to exit a volume
  ///
  /// @tparam propagator_state_t is the type of Propagatgor state
  /// @tparam stepper_t is the used type of the Stepper by the Propagator
  ///
  /// @param [in,out] state is the mutable propagator state object
  /// @param [in] stepper Stepper in use
  /// @param [in] targetLost, triggered from target() call
  template <typename propagator_state_t, typename stepper_t>
  void status(propagator_state_t& state, const stepper_t& stepper, bool targetLost=false) const {
    
    // Call the navigation helper prior to actual navigation
    debugLog(state, [&] { return std::string("Status | entering navigator."); });

    // The status interseciotn to be further used 
    SurfaceIntersection sIntersection = SurfaceIntersection();

    // Check the status of the target Surface 
    if (state.navigation.targetSurface){
      // Check the status to the target surface
      sIntersection 
        = checkStatus(state, stepper, 
                      *state.navigation.targetSurface,
                      state.options.targetBoundaryCheck,
                      Cstep::aborter, 
                      "Status | target");
      // Stop only if the target is hit 
      if (sIntersection.intersection.status == IntersectionStatus::onSurface){
        // Call the navigation helper prior to actual navigation
        debugLog(state, 
          [&] { return std::string("Status | target surface hit."); });
        return;
      }
    }
    
    // Status call initializes current and next target information
    state.navigation.currentSurface = nullptr;
    state.navigation.nextTarget.surface = nullptr;
    state.navigation.nextTarget.distance = std::numeric_limits<double>::max();
    state.navigation.nextTarget.boundaryCheck = true;
  
    // (a) Pre-stepping or re-initialze call from propagator
    if (not state.navigation.startVolume or targetLost) {
      // Initialize and return
      initialize(state, stepper, targetLost);
      return;
    }
    
    // (b) Status call within propagation loop
    // Check the external surfaces first (if applicable)
    if (!state.navigation.extSurfaces.empty() and 
        state.navigation.extSurfaceIter != state.navigation.extSurfaces.end()){
      
      debugLog(state, [&] {
        return std::string("Status | testing external surface.");
      });
      // Check the status of the current external surface
      sIntersection 
        = checkStatus(state, stepper, 
                      *state.navigation.extSurfaceIter->representation,
                      state.options.externalBoundaryCheck,
                      Cstep::actor, 
                      "Status | external");
      // If the external surface is hit, return
      if (sIntersection.intersection.status == IntersectionStatus::onSurface) {
        debugLog(state,
                 [&] { return std::string("Status | external surface hit."); });
        // Switch to the next external surface and return
        ++state.navigation.extSurfaceIter;
      } else if (sIntersection) {
        // External surface still reachable, leave the iterator there
        debugLog(state, [&] {
          return std::string("Status | external surface reachable.");
        });
        // Set the next target
        updateTarget(state, stepper, sIntersection,
                     state.options.externalBoundaryCheck);
      } else {
        // External surface is not reachable anymore, skip to next
        debugLog(state, [&] {
          return std::string(
              "Status | external surface not reachable anymore, skipping it.");
        });
        // switch to the next external surface
        ++state.navigation.extSurfaceIter;
      }
    } else {
      // Navigator without external surfaces
      debugLog(state, [&] {
        return std::string("Status | no external surfaces to handle.");
      });
    }
    // Now check the navigation surfaces (if applicable)
    if (!state.navigation.navSurfaces.empty() and 
        state.navigation.navSurfaceIter != state.navigation.navSurfaces.end()) {
      debugLog(state, [&] {
        return std::string("Status | testing navigation surface.");
      });
      // Check on layer surface handling
      sIntersection 
        = checkStatus(state, stepper,
                      *state.navigation.navSurfaceIter->representation,
                      true,
                      Cstep::actor,
                      "Status | navigation");
      if (sIntersection.intersection.status == IntersectionStatus::onSurface) {
        // we are on surface
        debugLog(state, [&] {
          return std::string("Status | navigation surface hit.");
        });
        // Update the surface, either to the next surface or layer
        if (updateSurfaces(state, stepper)) { return; }
      } else if (sIntersection) {
        // Still on path to the next navigation surface
        debugLog(state, [&] {
          return std::string("Status | navigation surface still reachable.");
        });
        // Keep the same surface
        updateTarget(state, stepper, sIntersection, true);
        return;
      } else {
        debugLog(state, [&] {
          return std::string(
              "Status | navigation surface not reachable anymore, skipping it.");
        });
        // Update the surface, either to the next surface or layer
        if (updateSurfaces(state, stepper)) { return; }        
      }
    } else {
      debugLog(state, [&] {
        return std::string("Status | no navigation surfaces to handle.");
      });
    }
    // Now check the navigation layers (if applicable)
    if (!state.navigation.navLayers.empty() and 
        state.navigation.navLayerIter != state.navigation.navLayers.end()) {
      debugLog(state,
               [&] { return std::string("Status | testing layer surface."); });
      // Surface handling has not triggered a return, try layer status
      sIntersection 
        = checkStatus(state, stepper,
                      *state.navigation.navLayerIter->representation,
                      true,
                      Cstep::actor,
                      "Status | layer");
      if (sIntersection.intersection.status == IntersectionStatus::onSurface) {
        debugLog(state, [&] {
          return std::string("Status | layer surface hit.");
        });
        // Resolve the surfaces on this layer
        if (resolveSurfaces(state, stepper)) {
          return;
        }
        // This layer is done, switch to next one
        if (updateLayers(state, stepper, state.navigation.navLayerIter->object)){
          return;
        }
      } else if (sIntersection) {
        // Still on path to the next layer surface
        debugLog(state,
                 [&] { return std::string("Status | layer still reachable."); });
        updateTarget(state, stepper, sIntersection, true);
        return;
      } else {
        debugLog(state, [&] {
          return std::string(
              "Status | layer surface not reachable anymore, skipping it.");
        });
      }
    } else {
      debugLog(state, [&] {
        return std::string("Status | no layers to handle.");
      });
    }

    // Layer handling has not triggered a return, try boundary resolving
    if (!state.navigation.navBoundaries.empty() and 
         state.navigation.navBoundaryIter != state.navigation.navBoundaries.end()){
      sIntersection 
        = checkStatus(state, stepper, 
                      *state.navigation.navBoundaryIter->representation,
                      true,
                      Cstep::actor,
                      "Status | boundary");
      if (sIntersection.intersection.status == IntersectionStatus::onSurface) {
        // Set the navigation to next boundary
           debugLog(state, [&] { return std::string("Status | boundary surface hit."); });   
        if (updateVolume(state,stepper)) return;
      } else if (sIntersection) {
        // Still on path to the next boundary surface
        debugLog(state, [&] {
          return std::string("Status | boundary surface still reachable.");
        });
        updateTarget(state, stepper, sIntersection, true);
        return;
      } else {
        // switch the volume information
        return;
      }
    } else if (updateBoundaries(state, stepper)){
      // Update the boundary information
      debugLog(state, [&] {
        return std::string("Status | new boundary surface set.");
      });
      return;
    }
    // Nothing left to do, check target progressing
    debugLog(state, [&] {
      return std::string("Status | nothing to do, check target.");
    });

    return;
  }

  /// @brief Navigator target call
  /// --------------------------------------------------------------------
  ///
  /// The target call checks if a previously set target is still reachable,
  /// and updates the target information in that case, if the target it lost
  /// it has to re-establish the status()
  ///
  /// @tparam propagator_state_t is the type of Propagatgor state
  /// @tparam stepper_t is the used type of the Stepper by the Propagator
  ///
  /// @param [in,out] state is the mutable propagator state object
  /// @param [in] stepper Stepper in use
  template <typename propagator_state_t, typename stepper_t>
  void target(propagator_state_t& state, const stepper_t& stepper) const {
    
    // Call the navigation helper prior to actual navigation
    debugLog(state, [&] { return std::string("Entering navigator::target."); });

    // If you are on target already just return
    if (state.navigation.currentSurface 
      and state.navigation.currentSurface == state.navigation.targetSurface){
      return;
    }
    // Check the current target surface
    if (state.navigation.nextTarget.surface and
          checkStatus(state, stepper, 
                      *state.navigation.nextTarget.surface,
                      state.navigation.nextTarget.boundaryCheck,
                      Cstep::actor,
                      "Target | next target")){
        return;
    }
    
    debugLog(state, [&] { return std::string("Next target lost, re-initalize."); });    

    // Target was lost, re-initialize, and re-establish the status
    status(state, stepper, true);
    // Return to the propagator
    return;
  }

 private:
   
  /// Initialize call - start of propagation
  /// --------------------------------------------------------------------
  ///
  /// @tparam propagator_state_t The state type of the propagagor
  /// @tparam stepper_t The type of stepper used for the propagation
  ///
  /// @param [in,out] state is the propagation state object
  /// @param [in] stepper Stepper in use
  /// @param [in] targetLost is the re-inialisation 
  ///
  /// @return boolean return triggers exit to stepper
  template <typename propagator_state_t, typename stepper_t>
  void initialize(propagator_state_t& state, 
                 const stepper_t& stepper, bool targetLost = false) const {
    
    // Call the navigation helper prior to actual navigation
    const std::string cstring = targetLost ? "Target |" : "Status |";               
    if (!targetLost){
      debugLog(state, [&] { return cstring+std::string(" Initalize."); });
    } else {
      debugLog(state, 
      [&] { return cstring+std::string(" Initialize after target lost."); });
    }              

    // Status call initializes current and next target information
    state.navigation.currentSurface = nullptr;
    state.navigation.nextTarget.surface = nullptr;
    state.navigation.nextTarget.distance = std::numeric_limits<double>::max();
    state.navigation.nextTarget.boundaryCheck = true;
    
    // We set the current surface to the start surface
    // for eventual post-update action, e.g. material integration
    // or collection when leaving a surface at the start of
    // an extrapolation process
    state.navigation.currentSurface = targetLost ?
        nullptr : state.navigation.startSurface;
    if (state.navigation.currentSurface) {
      debugLog(state, [&] {
        std::stringstream dstream;
        dstream << cstring;
        dstream << " Current surface set to start surface ";
        dstream << state.navigation.currentSurface->geoID().toString();
        return dstream.str();
      });
    }
    // Fast Navigation initialization for start condition:
    // - short-cut through object association, saves navigation in the
    // - geometry and volume tree search for the lowest volume
    if (state.navigation.startSurface &&
        state.navigation.startSurface->associatedLayer()) {
      debugLog(state, [&] {
        return cstring+std::string(" Fast start initialization through association.");
      });
      // Assign the current layer and volume by association
      state.navigation.startLayer =
          state.navigation.startSurface->associatedLayer();
      state.navigation.startVolume =
          state.navigation.startLayer->trackingVolume();
      // Set the start volume as current volume
      state.navigation.currentVolume = state.navigation.startVolume;
    } else {
      debugLog(state, [&] {
        return cstring+std::string(" Slow start initialization through search.");
      });
      // Current volume and layer search through global search
      debugLog(state, [&] {
        std::stringstream dstream;
        dstream << cstring;
        dstream << " Starting position "
                << toString(stepper.position(state.stepping));
        return dstream.str();
      });
      debugLog(state, [&] {
        std::stringstream dstream;
        dstream << cstring;
        dstream << " Starting direction "
                << toString(stepper.direction(state.stepping));
        return dstream.str();
      });
      state.navigation.startVolume = trackingGeometry->lowestTrackingVolume(
          state.geoContext, stepper.position(state.stepping));
      state.navigation.startLayer =
          state.navigation.startVolume
              ? state.navigation.startVolume->associatedLayer(
                    state.geoContext, stepper.position(state.stepping))
              : nullptr;
      // Set the start volume as current volume
      state.navigation.currentVolume = state.navigation.startVolume;
      if (state.navigation.startVolume) {
        debugLog(state, [&] { return cstring+std::string(" Start volume resolved."); });
      }      
      // Update the Layer information
      updateLayers(state, stepper, nullptr);
    }
    return;
  }
  
  /// Status call for test surfaces (surfaces, layers, boundaries)
  /// -------------------------------------------------
  ///
  /// If there are surfaces to be handled, check if the current
  /// state is on the surface. It calls the stepper.intersectSurface
  /// method on the given object (which will update the step size)
  /// and returns the intersection result:
  /// - onSurface: currentSurface is set and usually the caller reacts
  /// - reachable/overstepped: return back to propagator
  /// - unreachable: skip this one and go to the next
  ///
  /// @tparam propagator_state_t The state type of the propagagor
  /// @tparam stepper_t The type of stepper used for the propagation
  /// @tparam navigation_surfaces_t Type of the propagagor
  /// @tparam navigation_iter_t Type of the navigation iterator
  ///
  /// @param [in,out] state is the propagation state object
  /// @param [in] stepper Stepper in use
  /// @param [in] surface is the surface to be checked
  /// @param [in] bCheck the BoundaryCheck for this status check
  ///
  /// @returns a Surface intersection 
  template <typename propagator_state_t, typename stepper_t>
  SurfaceIntersection checkStatus(propagator_state_t& state,
        const stepper_t& stepper, const Surface& surface, 
        const BoundaryCheck& bcheck, Cstep::Type stype, 
        const std::string& cstring = "Undefined") const {
      // Take the stepper and intersect  
      auto sIntersection =
          stepper.intersectSurface(state.stepping, surface, bcheck, stype);
      if (sIntersection.intersection.status == IntersectionStatus::onSurface) {
        // Set in navigation state, so actors and aborters can access it
        state.navigation.currentSurface = &surface;
        debugLog(state, [&] {
          std::stringstream dstream;
          dstream << cstring;
          dstream << " surface successfully hit, setting ";
          dstream << state.navigation.currentSurface->geoID().toString();
          return dstream.str();
        });
        // switch the iterator to the next in line & releae the actor
        stepper.releaseStepSize(state.stepping, stype);
      } else if (sIntersection) {
        debugLog(state, [&] {
          std::stringstream dstream;
          dstream << cstring;
          dstream << " surface is still reachable";
          if (sIntersection.intersection.status ==
              IntersectionStatus::overstepped) {
            dstream << " (although overstepped)";
          }
          dstream << " - continuing towards.";
          return dstream.str();
        });
      } else {
        debugLog(state, [&] {
          std::stringstream dstream;
          dstream << cstring;
          dstream << " surface unreachable.";
          return dstream.str();
        });
        // Release the according step size
        stepper.releaseStepSize(state.stepping,stype);
      }
      // Return a positive status: either on it, or on the way
      return sIntersection;
  }

  /// Update the navigation surface 
  /// -------------------------------------------------
  ///
  /// This helper method increments the navigation surface on
  /// a layer, and if it is the last one, cleans up and switches
  /// to the next layer
  ///
  /// @tparam propagator_state_t The state type of the propagagor
  /// @tparam stepper_t The type of stepper used for the propagation
  ///
  /// @param [in,out] state is the propagation state object
  /// @param [in] stepper Stepper in use
  ///
  /// @return boolean 
  template <typename propagator_state_t, typename stepper_t>
  bool updateSurfaces(propagator_state_t& state,
                     const stepper_t& stepper) const {
    // incemenet
    ++state.navigation.navSurfaceIter;
    // This is the new navigation surface
    if (state.navigation.navSurfaceIter == state.navigation.navSurfaces.end()) {
      debugLog(state, [&] {
        std::stringstream dstream;
        dstream << "Last surface of this layer reached, update layers.";
        return dstream.str();
      });
      // reset the navigation surfaces
      state.navigation.navSurfaces.clear();
      state.navigation.navSurfaceIter = state.navigation.navSurfaces.end();
      // If you are in the layer to layer loop, update the reachable ones
      if (!state.navigation.navLayers.empty()) {
        // Update the reachable layers, but exclude the current one
        auto cLayer = state.navigation.navLayerIter->object;
        return updateLayers(state, stepper, cLayer);
      }
    } 
    // we switched to the next surface, update the distance towards it
    updateTarget(state, stepper, 
      *state.navigation.navSurfaceIter, true, true);
    return true;
  }

  /// Update Target if it is closer than the current target
  template <typename propagator_state_t, typename stepper_t,
            typename target_intersection_t>
  void updateTarget(propagator_state_t& state, const stepper_t& stepper,
                    const target_intersection_t& testTarget,
                    const BoundaryCheck& bcheck,
                    bool reevaluate = false) const {
                                            
    // The current target distance
    double distance = state.navigation.nextTarget.distance;
    // The target and the distance of the test
    auto testSurface = testTarget.representation;
    double testDistance = testTarget.intersection.pathLength;
        
    // This re-evaluates the absolute distance to the target
    if (reevaluate) {
      testDistance = state.stepping.navDir *
          (stepper.position(state.stepping) - testTarget.intersection.position)
              .norm();
            
    }    
    if (state.navigation.nextTarget.surface == nullptr or
        testDistance * testDistance < distance * distance) {    
      // update the target
      state.navigation.nextTarget.surface = testSurface;
      state.navigation.nextTarget.distance = testDistance;
      state.navigation.nextTarget.boundaryCheck = bcheck;
            
      // set the step size
      stepper.updateStepSize(state.stepping,testDistance,Cstep::actor,true);
      
      debugLog(state, [&] {
        std::stringstream dstream;
        dstream << "Target | update (actor) step size to ";
        dstream << testDistance;
        return dstream.str();
      });
    }
    return;
  }
  
  /// Update the layer targets 
  /// -------------------------------------------------
  ///
  /// This helper method increments the navigation surface on
  /// a layer, and if it is the last one, cleans up and switches
  /// to the next layer
  ///
  /// @tparam propagator_state_t The state type of the propagagor
  /// @tparam stepper_t The type of stepper used for the propagation
  ///
  /// @param [in,out] state is the propagation state object
  /// @param [in] stepper Stepper in use
  ///
  /// @return boolean 
  template <typename propagator_state_t, typename stepper_t>
  bool updateLayers(propagator_state_t& state, const stepper_t& stepper,
                    const Layer* excluded) const {
    auto navCorr = stepper.corrector(state.stepping);
    // - and get the compatible layers, start layer will be excluded
    Parameters<stepper_t> navPars(stepper, state.stepping);
    Options<Layer> navOpts(state.stepping.navDir, true, resolveSensitive,
                           resolveMaterial, resolvePassive, excluded, nullptr);
    // Set also the target surface
    navOpts.targetSurface = state.navigation.targetSurface;
    navOpts.pathLimit = stepper.stepSize(state.stepping, Cstep::aborter);
    // Request the compatible layers
    state.navigation.navLayers =
        state.navigation.currentVolume->compatibleLayers(
            state.geoContext, navPars, navOpts, navCorr);

    // Layer candidates have been found
    if (!state.navigation.navLayers.empty()) {
      // Screen output where they are
      debugLog(state, [&] {
        std::stringstream dstream;
        dstream << state.navigation.navLayers.size();
        dstream << " layer candidates found at path(s): ";
        for (auto& lc : state.navigation.navLayers) {
          dstream << lc.intersection.pathLength << "  ";
        }
        return dstream.str();
      });
      // Set to the first one
      state.navigation.navLayerIter = state.navigation.navLayers.begin();
      // Now update the target
      updateTarget(state, stepper, *state.navigation.navLayerIter, true);
      return true;
    }
    return false;
  }

  /// Update the layer targets 
  /// -------------------------------------------------
  ///
  /// This helper method updates the information to the boundsary
  ///
  /// @tparam propagator_state_t The state type of the propagagor
  /// @tparam stepper_t The type of stepper used for the propagation
  ///
  /// @param [in,out] state is the propagation state object
  /// @param [in] stepper Stepper in use
  ///
  /// @return boolean 
  template <typename propagator_state_t, typename stepper_t>
  bool updateBoundaries(propagator_state_t& state,
                        const stepper_t& stepper) const {
    // The navigation corrector  
    auto navCorr = stepper.corrector(state.stepping);
    // If we have no current volume, there's little to be done
    if (!state.navigation.currentVolume) {
      debugLog(state, [&] {
        return std::string(
            "No sufficient information to resolve boundary, "
            "stopping navigation.");
      });
      // release the step size
      stepper.releaseStepSize(state.stepping, Cstep::actor);
      return false;
    }

    // Boundary update invalidates layers & surfaces
    state.navigation.navSurfaces.clear();
    state.navigation.navSurfaceIter = state.navigation.navSurfaces.end();
    state.navigation.navLayers.clear();
    state.navigation.navLayerIter = state.navigation.navLayers.end();

    // The navigation parameter and options
    Parameters<stepper_t> navPars(stepper, state.stepping);
    Options<Surface> navOpts(state.stepping.navDir, true);
    navOpts.pathLimit = state.stepping.stepSize.value(Cstep::aborter);

    // Evaluate the boundary surfaces
    state.navigation.navBoundaries =
        state.navigation.currentVolume->compatibleBoundaries(
            state.geoContext, navPars, navOpts, navCorr);
    if (!state.navigation.navBoundaries.empty()) {
      // The number of boundary candidates
      debugLog(state, [&] {
        std::stringstream dstream;
        dstream << state.navigation.navBoundaries.size();
        dstream << " boundary candidates found at path(s): ";
        for (auto& bc : state.navigation.navBoundaries) {
          dstream << bc.intersection.pathLength << "  ";
        }
        return dstream.str();
      });
      // Set the begin iterator
      state.navigation.navBoundaryIter = state.navigation.navBoundaries.begin();
      // Update the target to the Boundary
      updateTarget(state, stepper, *state.navigation.navBoundaryIter,true);
      return true;
    }
    // Nothing to do, boundary could not be resolved
    return false;
  }

  /// Update the volume information
  /// -------------------------------------------------
  ///
  /// This helper method updates the volume and is called when a 
  /// boundary surface is reached.
  ///
  /// @tparam propagator_state_t The state type of the propagagor
  /// @tparam stepper_t The type of stepper used for the propagation
  ///
  /// @param [in,out] state is the propagation state object
  /// @param [in] stepper Stepper in use
  ///
  /// @return boolean 
  template <typename propagator_state_t, typename stepper_t>
  bool updateVolume(propagator_state_t& state,
                    const stepper_t& stepper) const {
                      
    // Update volume information
    // get the attached volume information
    auto boundary = state.navigation.navBoundaryIter->object;
    state.navigation.currentVolume = boundary->attachedVolume(
        state.geoContext, stepper.position(state.stepping),
        stepper.direction(state.stepping), state.stepping.navDir);
    // No volume anymore : end of known world
    if (!state.navigation.currentVolume) {
      debugLog(state, [&] {
        return std::string(
            "No more volume to progress to, stopping navigation.");
      });
      // Navigation break & release navigation stepping
      state.navigation.navigationBreak = true;
      state.stepping.stepSize.release(Cstep::actor);
      return true;
    } 
    debugLog(state, [&] { return std::string("Volume updated."); });
    // Forget the bounday information
    state.navigation.navBoundaries.clear();
    state.navigation.navBoundaryIter =
        state.navigation.navBoundaries.end();
    return updateLayers(state,stepper,nullptr);
 }
  
  /// @brief Resolve the surfaces of this layer
  ///
  /// @tparam propagator_state_t The state type of the propagagor
  /// @tparam stepper_t The type of stepper used for the propagation
  /// @tparam is a step corrector (can be void corrector)
  ///
  /// @param [in,out] state is the propagation state object
  /// @param [in] stepper Stepper in use
  /// @param [in] navCorr is the navigation path corrector
  /// @param [in] cLayer is the already assigned current layer, e.g. start layer
  ///
  /// boolean return triggers exit to stepper
  template <typename propagator_state_t, typename stepper_t>
  bool resolveSurfaces(propagator_state_t& state, const stepper_t& stepper,
                       const Layer* cLayer = nullptr) const {
    auto navCorr = stepper.corrector(state.stepping);
    // get the layer and layer surface
    auto layerSurface = cLayer ? state.navigation.startSurface
                               : state.navigation.navLayerIter->representation;
    auto navLayer = cLayer ? cLayer : state.navigation.navLayerIter->object;
    // are we on the start layer
    Parameters<stepper_t> navPars(stepper, state.stepping);
    bool onStart = (navLayer == state.navigation.startLayer);
    auto startSurface = onStart ? state.navigation.startSurface : layerSurface;
    Options<Surface> navOpts(state.stepping.navDir, true, resolveSensitive,
                             resolveMaterial, resolvePassive, startSurface,
                             state.navigation.targetSurface);
    // Check the limit
    navOpts.pathLimit = state.stepping.stepSize.value(Cstep::aborter);
    // get the surfaces
    state.navigation.navSurfaces = navLayer->compatibleSurfaces(
        state.geoContext, navPars, navOpts, navCorr);
    // the number of layer candidates
    if (!state.navigation.navSurfaces.empty()) {
      debugLog(state, [&] {
        std::stringstream dstream;
        dstream << state.navigation.navSurfaces.size();
        dstream << " surface candidates found at path(s): ";
        for (auto& sfc : state.navigation.navSurfaces) {
          dstream << sfc.intersection.pathLength << "  ";
        }
        return dstream.str();
      });
      // Set the iterator
      state.navigation.navSurfaceIter = state.navigation.navSurfaces.begin();
      // this is the new target stream
      updateTarget(state, stepper, *state.navigation.navSurfaceIter, true);
      return true;
    }
    state.navigation.navSurfaceIter = state.navigation.navSurfaces.end();
    debugLog(state,
             [&] { return std::string("No surface candidates found."); });
    return false;
  }

  /// The private navigation debug logging
  ///
  /// It needs to be fed by a lambda function that returns a string,
  /// that guarantees that the lambda is only called in the
  /// state.options.debug == true
  /// case in order not to spend time when not needed.
  ///
  /// @tparam propagator_state_t Type of the propagator state
  ///
  /// @param[in,out] state the propagator state for the debug flag,
  ///      prefix and length
  /// @param logAction is a callable function that returns a stremable object
  template <typename propagator_state_t>
  void debugLog(propagator_state_t& state,
                const std::function<std::string()>& logAction) const {
    if (state.options.debug) {
      std::string vName = "No Volume";
      if (state.navigation.currentVolume) {
        vName = state.navigation.currentVolume->volumeName();
      }
      std::vector<std::string> lines;
      std::string input = logAction();
      boost::split(lines, input, boost::is_any_of("\n"));
      for (const auto& line : lines) {
        std::stringstream dstream;
        dstream << ">>>" << std::setw(state.options.debugPfxWidth) << vName
                << " | ";
        //dstream << std::setw(state.options.debugMsgWidth) << line << '\n';
        dstream << line << '\n';        
        state.options.debugString += dstream.str();
      }
    }
  }
};

}  // namespace Acts