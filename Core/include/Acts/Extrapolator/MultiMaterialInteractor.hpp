// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <cmath>
#include <sstream>
#include <utility>
#include "Acts/Extrapolator/MaterialInteractor.hpp"
#include "Acts/Extrapolator/detail/EmptyEffect.hpp"
#include "Acts/Extrapolator/detail/InteractionFormulas.hpp"
#include "Acts/Material/ISurfaceMaterial.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Material/MaterialProperties.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Helpers.hpp"

namespace Acts {

struct InteractionPoint
{
  /// The material surface
  const Surface* surface = nullptr;
  /// The position information of the material hit
  Vector3D position = Vector3D(0., 0., 0);
  /// The direction information of the material hit
  Vector3D direction = Vector3D(0., 0., 0);
  /// The calculated path & applied path correction factor

  bool
  operator==(const InteractionPoint& others) const
  {
    if (fabs((this->direction - others.direction).norm()) > 1e-10) {
      return false;
    }
    if (fabs((this->position - others.position).norm()) > 1e-10) {
      return false;
    }
    if (this->surface != others.surface) {
      return false;
    }
    return true;
  }
};

using InteractionPointVec = std::vector<Acts::InteractionPoint>;

/// The Material interactor struct
///
/// This is a plugin to the Propagator that
/// performs material interaction on the currentSurface
/// of the Propagagor state
struct MultiMaterialInteractor
{
  // Configuration for this MultiMaterialInteractor

  /// multiple scattering switch on/off
  bool multipleScattering = true;
  /// The scattering formula struct
  detail::HighlandScattering scattering;

  /// Bethe-Hitler struct
  /// currently empty
  detail::EmptyEffect emptyEffect;

  /// Record material in detail
  bool recordInteractions = true;

  /// Simple result struct to be returned
  /// It mainly acts as an internal state which is
  /// created for every propagation/extrapolation step
  struct this_result
  {
    /// This one is only filled when recordInteractions is switched on
    std::map<const Surface*, InteractionPointVec> multiMaterialInteractions;
    /// the number of components
    int numComponents = 0;
  };

  using result_type = this_result;

  /// @brief Interaction with detector material for the ActionList
  /// of the Propagator
  ///
  /// only contains a split function currently
  ///
  /// It checks if the state has a current surface, in which case
  /// the action is performed: the covariance is transported to the position,
  /// multiple scattering and energy loss is applied  according to the
  /// configuration.
  ///
  /// @tparam propagator_state_t is the type of Propagagor state
  /// @tparam stepper_t Type of the stepper of the propagation
  ///
  /// @param state is the mutable propagator state object
  /// @param stepper The stepper in use
  /// @param result is the mutable result state object
  //
  template <typename propagator_state_t, typename stepper_t>
  void
  operator()(propagator_state_t& state,
             const stepper_t&    stepper,
             result_type&        result) const
  {
    debugLog(state, [&] { return std::string("in MultiMaterialInteractor."); });

    // If we are on target, everything should have been done
    if (state.navigation.targetReached) {
      return;
    }

    // If switched off, then return - alows run-time configuration
    if (!multipleScattering && !recordInteractions) {
      return;
    }

    // A current surface has been already assigned by the navigator
    // check for material
    if (state.navigation.currentSurface
        && state.navigation.currentSurface->surfaceMaterial()) {
      // Let's set the pre/full/post update stage
      MaterialUpdateStage mStage = fullUpdate;
      // We are at the start surface
      if (state.navigation.startSurface == state.navigation.currentSurface) {
        debugLog(state, [&] {
          return std::string("Update on start surface: post-update mode.");
        });
        mStage = postUpdate;
        // Or is it the target surface ?
      } else if (state.navigation.targetSurface
                 == state.navigation.currentSurface) {
        debugLog(state, [&] {
          return std::string("Update on target surface: pre-update mode");
        });
        mStage = preUpdate;
      } else {
        debugLog(state, [&] {
          return std::string("Update while pass through: full mode.");
        });
      }

      /// Get the surface material & properties from them and continue if you
      /// found some
      const ISurfaceMaterial* sMaterial
          = state.navigation.currentSurface->surfaceMaterial();

      const double        m = state.options.mass;
      InteractionPointVec materialInteractionVec;

      /// typename of std::list<tuple<state,weight,status>>
      using stateColType = decltype(state.stepping.stateCol);

      /// loop the single components in the list
      typename stateColType::iterator it = state.stepping.stateCol.begin();
      while (it != state.stepping.stateCol.end()) {
        ///  apply energy-lose component split & gaussian multiple scattering
        auto& singlestate = std::get<0>((*it));

        MaterialProperties mProperties = sMaterial->materialProperties(
            singlestate.pos, singlestate.navDir, mStage);
        // Material properties (non-zero) have been found for this configuration
        if (mProperties) {
          /// more debugging output to the screen
          debugLog(state, [&] {
            return std::string("Material properties found for this surface.");
          });
          // Calculate the path correction
          double pCorrection = state.navigation.currentSurface->pathCorrection(
              state.geoContext,
              stepper.position(singlestate),
              stepper.direction(singlestate));

          // Scale the material properties
          mProperties *= pCorrection;

          // Create the material interaction class, in case we record afterwards
          // Record the material interaction if configured to do so
          Acts::InteractionPoint mInteraction;
          if (recordInteractions) {
            mInteraction.surface   = state.navigation.currentSurface;
            mInteraction.position  = stepper.position(singlestate);
            mInteraction.direction = stepper.direction(singlestate);
            materialInteractionVec.push_back(std::move(mInteraction));
          }  // end of record

          // To integrate process noise, we need to transport
          // the covariance to the current position in space
          // the 'true' indicates re-initializaiton of the further transport
          if (singlestate.covTransport) {
            stepper.covarianceTransport(singlestate, true);
          }

          // @brief if meets material surface , split the current single state
          // to get newlist
          // this should be according to the Bethe-Hitler pdf
          // in the function, return the list of multicomponents, each component
          // carries (weight,mean,variance)
          // @note current not use the Bethe-Hitler pdf, just make a list of
          // copied component, to see if they act equally in the multi-stepper,
          // the split is 2
          debugLog(state, [&] { return std::string("in Split method."); });
          stateColType newList = makeNewComponetList((*it), mProperties, m);

          // insert the newlist to the multi-component colume
          state.stepping.stateCol.insert(it, newList.begin(), newList.end());

          // delete the current component
          it = state.stepping.stateCol.erase(it);
        }
      }
      // record the material interaction
      if (recordInteractions) {
        result.multiMaterialInteractions.insert(
            std::pair<const Surface*, InteractionPointVec>(
                state.navigation.currentSurface,
                std::move(materialInteractionVec)));
      }
      // record the number of components at the last step
      result.numComponents = state.stepping.stateCol.size();
    }
  }

  /// Pure observer interface
  /// This does not apply to the surface collector
  template <typename propagator_state_t>
  void
  operator()(propagator_state_t& /*state*/) const
  {
  }

  /// @brief  simple spliting method, copy 1 component to N components with
  /// (dMean,dVariance = 0) and  equal weight
  /// @param [in] single component state
  /// @param [in] properties - the property of mater, determine the Bethe-Hitler
  /// distribution
  /// @param [in] m - mass
  /// @param [out] the list of split components
  template <typename tuplestate_t>
  std::list<tuplestate_t>
  makeNewComponetList(const tuplestate_t&       tuple_state,
                      const MaterialProperties& properties,
                      const double              m) const
  {
    auto&                   singleState = std::get<0>(tuple_state);
    double                  weight      = std::get<1>(tuple_state);
    const auto              status      = std::get<2>(tuple_state);
    const double            p           = singleState.p;
    const double            E           = std::sqrt(p * p + m * m);
    const double            lbeta       = p / E;
    const double            tInX0       = properties.thicknessInX0();
    auto                    mixture     = emptyEffect.getMixture(tInX0, p);
    std::list<tuplestate_t> splitList;
    unsigned int            iComponent = 0;
    while (iComponent < mixture.size()) {
      // energy loss
      const double dE = mixture[iComponent].mean;
      if (E + dE > m) {
        // p
        const double newP = std::sqrt((E + dE) * (E + dE) - m * m);
        // copy
        auto state = singleState;
        state.p    = std::copysign(newP, singleState.p);
        // cov for each new component
        if (state.covTransport) {
          // current variance varity : 0
          const double sigmaQoverP
              = mixture[iComponent].variance / (lbeta * newP * newP);
          // good in any case for positive direction
          if (state.navDir == forward) {
            state.cov(eQOP, eQOP) += state.navDir * sigmaQoverP * sigmaQoverP;
          } else {
            // check that covariance entry doesn't become negtive
            double sEqop = state.cov(eQOP, eQOP);
            if (sEqop > sigmaQoverP * sigmaQoverP) {
              state.cov(eQOP, eQOP) += state.navDir * sigmaQoverP * sigmaQoverP;
            }
          }
        }

        // multiple scattering
        if (multipleScattering && state.covTransport) {
          double sigmaScat = scattering(p, lbeta, tInX0);
          double sinTheta  = std::sin(VectorHelpers::theta(state.dir));
          double sigmaDeltaPhiSq
              = sigmaScat * sigmaScat / (sinTheta * sinTheta);
          double sigmaDeltaThetaSq = sigmaScat * sigmaScat;
          // Good in any case for positive direction
          if (state.navDir == forward) {
            // Just add the multiple scattering component
            state.cov(ePHI, ePHI) += state.navDir * sigmaDeltaPhiSq;
            state.cov(eTHETA, eTHETA) += state.navDir * sigmaDeltaThetaSq;
          } else {
            // We check if the covariance stays positive
            double sEphi   = state.cov(ePHI, ePHI);
            double sEtheta = state.cov(eTHETA, eTHETA);
            if (sEphi > sigmaDeltaPhiSq && sEtheta > sigmaDeltaThetaSq) {
              // Noise removal is not applied if covariance would fall below 0
              state.cov(ePHI, ePHI) -= sigmaDeltaPhiSq;
              state.cov(eTHETA, eTHETA) -= sigmaDeltaThetaSq;
            }
          }
        }

        // get the weight of the new component
        double pdfWeight = mixture[iComponent].weight;
        double newWeight = pdfWeight * weight;
        splitList.push_back(
            std::make_tuple(std::move(state), newWeight, status));
        iComponent++;
      }
    }
    return std::move(splitList);
  }

private:
  /// The private propagation debug logging
  ///
  /// It needs to be fed by a lambda function that returns a string,
  /// that guarantees that the lambda is only called in the state.debug == true
  /// case in order not to spend time when not needed.
  ///
  /// @tparam propagator_state_t Type of the propagator state
  ///
  /// @param state the propagator state for the debug flag, prefix and
  /// length
  /// @param logAction is a callable function that returns a stremable object
  template <typename propagator_state_t>
  void
  debugLog(propagator_state_t&                 state,
           const std::function<std::string()>& logAction) const
  {
    if (state.options.debug) {
      std::stringstream dstream;
      dstream << "   " << std::setw(state.options.debugPfxWidth);
      dstream << "material interaction"
              << " | ";
      dstream << std::setw(state.options.debugMsgWidth) << logAction() << '\n';
      state.options.debugString += dstream.str();
    }
  }
};

}  // end of namespace Acts
