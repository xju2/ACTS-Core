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
#include "Acts/Extrapolator/detail/InteractionFormulas.hpp"
#include "Acts/Material/ISurfaceMaterial.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Material/MaterialProperties.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Helpers.hpp"

namespace Acts {

/// @brief The Multi Material interaction struct
///
struct MultiMaterialInteraction
{
  /// The material surface
  const Surface* surface = nullptr;

  /// The position information of the material hit
  Vector3D position = Vector3D(0., 0., 0);
  /// The direction information of the material hit
  Vector3D direction = Vector3D(0., 0., 0);
  /// The calculated path & applied path correction factor
  double pathCorrection = 1.;
  /// The (passsed) material properties
  /// it is the material and the actual (corrected) path length
  MaterialProperties materialProperties = MaterialProperties();
};

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

  /// Energy loss switch on/off
  bool energyLoss = true;
  /// The energy loss formula struct
  detail::IonisationLoss ionisationloss;

  /// Record material in detail
  bool recordInteractions = true;

  /// Simple result struct to be returned
  /// It mainly acts as an internal state which is
  /// created for every propagation/extrapolation step
  struct this_result
  {
    // The accumulated materialInX0
    double materialInX0 = 0.;
    /// The accumulated materialInL0
    double materialInL0 = 0.;
    /// This one is only filled when recordInteractions is switched on
    std::vector<MultiMaterialInteraction> materialInteractions;
	///
	int numComponents = 0;
  };

  using result_type = this_result;

  /// @brief Interaction with detector material for the ActionList
  /// of the Propagator
  ///
  /// only contains a split function now
  /// @to do apply the multiple scattering cov transform
  //
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
    if (!multipleScattering && !energyLoss && !recordInteractions) {
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

      // typename of the tuple<state,weight,status>
      using stateColType   = decltype(state.stepping.stateCol);
      using tupleStateType = typename stateColType::value_type;

      // loop the single components in the list
      typename stateColType::iterator it = state.stepping.stateCol.begin();
      // N to be tuning number
      const int              N = 2;
      const double           m = state.options.mass;
      const ISurfaceMaterial* sMaterial
          = state.navigation.currentSurface->surfaceMaterial();
	  MultiMaterialInteraction mInteraction;
	  mInteraction.surface = state.navigation.currentSurface;  //*
      while (it != state.stepping.stateCol.end()) {
        auto& singlestate = std::get<0>((*it));
        // Get the surface material & properties from them and continue if you
        // found some
        MaterialProperties mProperties = sMaterial->materialProperties(
            singlestate.pos, singlestate.navDir, mStage);
        // Material properties (non-zero) have been found for this configuration
        if (mProperties) {
          /// more debugging output to the screen
          debugLog(state, [&] {
            return std::string("Material properties found for this surface.");
          });

          // Create the material interaction class, in case we record afterwards

          // To integrate process noise, we need to transport
          // the covariance to the current position in space
          // the 'true' indicates re-initializaiton of the further transport
          if (singlestate.covTransport) {
            stepper.covarianceTransport(singlestate, true);
          }

          // Calculate the path correction
          double pCorrection = state.navigation.currentSurface->pathCorrection(
            state.geoContext,
              stepper.position(singlestate), stepper.direction(singlestate));

          // Scale the material properties
          mProperties *= pCorrection;

          // the corrected thickness
          // const double cThickness = mProperties.thickness();

          // if meets material surface , split the current single state to get
          // newlist
          // here just copy
          debugLog(state, [&] { return std::string("in Split method."); });
          stateColType newList
              = split<tupleStateType>((*it), N, mProperties, m);

          // insert the newlist
          state.stepping.stateCol.insert(it, newList.begin(), newList.end());

          // delete the old state
          it = state.stepping.stateCol.erase(it);

          // This doesn't cost anything - do it regardless
          result.materialInX0 += mProperties.thicknessInX0();
          result.materialInL0 += mProperties.thicknessInL0();

          // Record the material interaction if configured to do so
          if (recordInteractions) {
            mInteraction.position           = stepper.position(singlestate);
            mInteraction.direction          = stepper.direction(singlestate);
            mInteraction.materialProperties = mProperties;
            mInteraction.pathCorrection     = pCorrection;
            result.materialInteractions.push_back(std::move(mInteraction));
          }  // end of record
        }
      }  // end of loop
	  result.numComponents  =  state.stepping.stateCol.size();
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
  /// (mean,variance) and  weight
  /// @param [in] single component state
  /// @param [in] N - the split components
  /// @param [in] properties - the property of mater, determine the Bethe-Hitler
  /// distribution
  /// @param [in] m - mass
  /// @param [out] the list of split components
  template <typename tuplestate_t>
  std::list<tuplestate_t>
  split(const tuplestate_t& tuple_state,
        const int           N,
        // const double     cThickness,
        const MaterialProperties& properties,
        const double              m) const
  {
    auto&        singleState = std::get<0>(tuple_state);
    double       weight      = std::get<1>(tuple_state);
    auto&        status      = std::get<2>(tuple_state);
    const double p           = singleState.p;
    const double E           = std::sqrt(p * p + m * m);
    const double lbeta       = p / E;
    // const double cThickness  = properties.thickness();
    const double tInX0 = properties.thicknessInX0();
    std::tuple<std::vector<double>, std::vector<double>, std::vector<double>>
                            eLossMap = Brem(N);
    std::list<tuplestate_t> splitList;
    int                     s = 0;
    while (s < N) {
      // energy loss
      const double dE = std::get<0>(eLossMap).at(s);
      if (E + dE > m) {
        // p
        const double newP  = std::sqrt((E + dE) * (E + dE) - m * m);
        auto         state = singleState;
        state.p            = std::copysign(newP, singleState.p);
        // cov for energy loss
        if (state.covTransport) {
          const double sigmaQoverP
              = std::get<1>(eLossMap).at(s) / (lbeta * newP * newP);
          // good in any case for positive direction
          if (state.navDir == forward) {
            state.cov(eQOP, eQOP) += state.navDir * sigmaQoverP * sigmaQoverP;
          } else {
            // check that covariance entry doesn't become negtive
            double sEqop = state.cov(eQOP, eQOP);
            if (sEqop > sigmaQoverP * sigmaQoverP) {
              state.cov(eQOP, eQOP) += state.navDir * sigmaQoverP * sigmaQoverP;
            }
          }  // end of forward or not
        }    // end of cov

        // multiple scattering
        if (multipleScattering) {
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
        }  // end of multiple scattering
        double pdfWeight = std::get<2>(eLossMap).at(s);
        double reweight  = pdfWeight * weight;
        splitList.push_back(std::make_tuple(state, reweight, status));
        s++;
      }  // end of (E+dE>m)
    }    // end of while
    return splitList;
  }  // end of function split

  // @brief function of Bremth - return mean,variance,weight
  std::tuple<std::vector<double>, std::vector<double>, std::vector<double>>
  Brem(const int N) const
  {
    std::vector<double> Emean(N, 0);  // here means energy loss
    std::vector<double> Vmean(N, 0);
    std::vector<double> Wmean(N, 0);
    for (int i = 0; i < N; i++) {
      Emean[i] = 0.;
      Vmean[i] = 0.;
      Wmean[i] = 1. / N;
    }
    // here should input some parameter
    // then transform, i.e.do the tranform ln(uk)-ln(1-uk)
    // weight should be unity
    return std::make_tuple(Emean, Vmean, Wmean);
  }

  // method of number of components fixed to a compositary number of (N)
  void
  reduction_by_close_components() const /*unused*/
  {
    // compare each two of the state list
    // each time combine the closest two
    // till the Maximum number
  }
  void
  reduction_by_largest_weight() const /*unused*/
  {
  }
  // method of the distance of two components
  // with the K-L distance method
  // Dkl(p1,p2) = tr[(V1-V2)(W1-W2)]+(u1-u2)^T(W1+W2)(u1-u2)
  void
  distance() const /*unused*/
  {
  }
  // method of combine two closed components
  void
  combine() const /*unused*/
  {
  }
  // method of correct the 1st 2nd order motion of components
  void
  corrected() const /*unused*/
  {
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
