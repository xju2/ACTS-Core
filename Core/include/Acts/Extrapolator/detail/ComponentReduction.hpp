// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#include "Acts/Extrapolator/detail/ComponentCombiner.hpp"
#include "Acts/Extrapolator/detail/ComponentDistance.hpp"
#include "component_reduction_impl.hpp"

namespace Acts {

namespace detail {
  /// The component reduction struct
  /// This is plugin to the Propagator that
  /// performs to reduce the number component
  struct ComponentReduction
  {
    detail::KullbackLeiblerComponentDistance klDist;
    detail::ComponentCombiner                combiner;

    /// the number of component constrained in propagate
    static const int constraintNum = 6;

    template <typename propagator_state_t, typename stepper_t>
    void
    operator()(propagator_state_t& state, const stepper_t& stepper) const
    {
      debugLog(state, [&] { return std::string("In component reduction."); });

      // If we are on target, everything should have been done
      if (state.navigation.targetReached) {
        return;
      }

      // A current surface has been already assigned by the navigator
      if (state.navigation.currentSurface) {

        // Get the multi bound parameter
        auto bs = stepper.boundState(
            state.stepping, *state.navigation.currentSurface, true);
        auto& multipleBoundPar = std::get<MultipleBoundParameters>(bs);
        auto& trackMap         = multipleBoundPar.getTrackList();
        stepper.outPut(state.stepping);

        impl::reductComponent(trackMap,
                              constraintNum,
                              state.stepping.geoContext,
                              *state.navigation.currentSurface);
        stepper.update(state.stepping, multipleBoundPar);
      }
    }

  private:
    /// The private propagation debug logging
    ///
    /// It needs to be fed by a lambda function that returns a string,
    /// that guarantees that the lambda is only called in the state.debug ==
    /// true
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
        dstream << "component reduction "
                << " | ";
        dstream << std::setw(state.options.debugMsgWidth) << logAction()
                << '\n';
        state.options.debugString += dstream.str();
      }
    }
  };

}  // end of namespace detail

}  // end of namespace Acts
