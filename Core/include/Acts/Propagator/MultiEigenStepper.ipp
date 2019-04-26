// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

template <typename B, typename C, typename E, typename A>
Acts::MultiEigenStepper<B, C, E, A>::MultiEigenStepper(B bField)
  : m_bField(std::move(bField))
{
}

template <typename B, typename C, typename E, typename A>
auto
Acts::MultiEigenStepper<B, C, E, A>::boundState(MultiState& state,
                                           const Surface& surface,
                                           bool           reinitialize) const
    -> BoundState
{
  // Transport the covariance to here
  std::unique_ptr<const Covariance> covPtr = nullptr;
  // Create the bound parameters
  BoundParameters parameters(state.geoContext,
                             std::move(covPtr),
                             state.pos,
                             state.p * state.dir,
                             state.q,
                             surface.getSharedPtr());
  // Create the bound state
  BoundState bState{
      std::move(parameters), state.jacobian, state.pathAccumulated};
  // Reset the jacobian to identity
  if (reinitialize) {
    state.jacobian = Jacobian::Identity();
  }
  /// Return the State
  return bState;
}

template <typename B, typename C, typename E, typename A>
auto
Acts::MultiEigenStepper<B, C, E, A>::curvilinearState(MultiState& state,
                                                 bool reinitialize) const
    -> CurvilinearState
{
  // Transport the covariance to here
  std::unique_ptr<const Covariance> covPtr = nullptr;
  // Create the curvilinear parameters
  CurvilinearParameters parameters(
      std::move(covPtr), state.pos, state.p * state.dir, state.q);
  // Create the bound state
  CurvilinearState curvState{
      std::move(parameters), state.jacobian, state.pathAccumulated};
  // Reset the jacobian to identity
  if (reinitialize) {
    state.jacobian = Jacobian::Identity();
  }
  /// Return the State
  return curvState;
}

template <typename B, typename C, typename E, typename A>
void
Acts::MultiEigenStepper<B, C, E, A>::update(SingleStateType& state,
                                       const BoundParameters& pars) const
{
	EigenStepper<B>::update(state, pars);
}

template <typename B, typename C, typename E, typename A>
void
Acts::MultiEigenStepper<B, C, E, A>::update(SingleStateType& state,
                                       const Vector3D& uposition,
                                       const Vector3D& udirection,
                                       double          up) const
{
	EigenStepper<B>::update(state,uposition,udirection,up);
}

template <typename B, typename C, typename E, typename A>
template <typename propagator_state_t>
Acts::Result<double>
Acts::MultiEigenStepper<B, C, E, A>::step(propagator_state_t& state) const
  {
    // a variant to record the combination of the pathlength of the compact
    double combine_h = 0;
    // loop all the components that are free
    for (auto& tuple_state : state.stepping.stateCol) {
      /// if the status is locked or dead, ignore it
      if (std::get<2>(tuple_state) != StateStatus::FREE) {
        continue;
      }
      auto& singlestate = std::get<0>(tuple_state);

      // Runge-Kutta integrator state
      auto& sd = singlestate.stepData;

      double h2, half_h;
      double error_estimate;

      // First Runge-Kutta point (at current position)
      sd.B_first = getField(singlestate, singlestate.pos);
      if (!singlestate.extension.validExtensionForStep(
              state, *this, singlestate)
          || !singlestate.extension.k1(
                 state, *this, singlestate, sd.k1, sd.B_first)) {
        return 0.;
      }

      // The following functor starts to perform a Runge-Kutta step of a certain
      // size, going up to the point where it can return an estimate of the
      // local
      // integration error. The results are stated in the local variables above,
      // allowing integration to continue once the error is deemed satisfactory
      const auto tryRungeKuttaStep = [&](const double h) -> bool {

        // State the square and half of the step size
        h2     = h * h;
        half_h = h * 0.5;

        // Second Runge-Kutta point
        const Vector3D pos1
            = singlestate.pos + half_h * singlestate.dir + h2 * 0.125 * sd.k1;
        sd.B_middle = getField(singlestate, pos1);
        if (!singlestate.extension.k2(
                state, *this, singlestate, sd.k2, sd.B_middle, half_h, sd.k1)) {
          return false;
        }

        // Third Runge-Kutta point
        if (!singlestate.extension.k3(
                state, *this, singlestate, sd.k3, sd.B_middle, half_h, sd.k2)) {
          return false;
        }

        // Last Runge-Kutta point
        const Vector3D pos2
            = singlestate.pos + h * singlestate.dir + h2 * 0.5 * sd.k3;
        sd.B_last = getField(singlestate, pos2);
        if (!singlestate.extension.k4(
                state, *this, singlestate, sd.k4, sd.B_last, h, sd.k3)) {
          return false;
        }

        // Return an estimate of the local integration error
        error_estimate = std::max(
            h2 * (sd.k1 - sd.k2 - sd.k3 + sd.k4).template lpNorm<1>(), 1e-20);
		return (error_estimate <= state.options.tolerance);
      };

      double stepSizeScaling;

      // Select and adjust the appropriate Runge-Kutta step size as given
      // ATL-SOFT-PUB-2009-001
      while (!tryRungeKuttaStep(singlestate.stepSize)) {
        stepSizeScaling
            = std::min(std::max(0.25,
                                std::pow((state.options.tolerance
                                          / std::abs(error_estimate)),
                                         0.25)),
                       4.);
        if (stepSizeScaling == 1.) {
          break;
        }
        singlestate.stepSize = singlestate.stepSize * stepSizeScaling;

        // If step size becomes too small the particle remains at the initial
        // place
        if (singlestate.stepSize * singlestate.stepSize < state.options.stepSizeCutOff * state.options.stepSizeCutOff ) {
		  // Not moving due to too low momentum needs an aborter
		  return EigenStepperError::StepSizeStalled;
        }
      }

      // use the adjusted step size
      const double h = singlestate.stepSize;

      // When doing error propagation, update the associated Jacobian matrix
      if (singlestate.covTransport) {
        // The step transport matrix in global coordinates
        ActsMatrixD<7, 7> D;
        if (!singlestate.extension.finalize(state, *this, singlestate, h, D)) {
		  return EigenStepperError::StepInvalid;
        }

        // for moment, only update the transport part
        singlestate.jacTransport = D * singlestate.jacTransport;
      } else {
        if (!singlestate.extension.finalize(state, *this, singlestate, h)) {
		  return EigenStepperError::StepInvalid;
        }
      }

      // Update the track parameters according to the equations of motion
      singlestate.pos
          += h * singlestate.dir + h2 / 6. * (sd.k1 + sd.k2 + sd.k3);
      singlestate.dir += h / 6. * (sd.k1 + 2. * (sd.k2 + sd.k3) + sd.k4);
      singlestate.dir /= singlestate.dir.norm();
      singlestate.derivative.template head<3>()     = singlestate.dir;
      singlestate.derivative.template segment<3>(3) = sd.k4;
      singlestate.pathAccumulated += h;

      // for the multistepper simply record the pathlength with combination of
      // components
      //
      state.stepping.pathAccumulated += h * std::get<1>(tuple_state);
      combine_h += h * std::get<1>(tuple_state);
    }
    // the pos/dir/mom after RKN is a combination of the all
    state.stepping.pos = position(state.stepping);
    state.stepping.dir = direction(state.stepping);
    state.stepping.p   = momentum(state.stepping);
    return combine_h;
  }
