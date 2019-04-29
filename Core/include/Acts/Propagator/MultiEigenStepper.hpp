// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Propagator/EigenStepper.hpp"

namespace Acts {

enum StateStatus : int { FREE = 0, LOCKED = 1, DEAD = 2 };

/// @brief multicomponent stepper(mcs) is based on the Single Stepper
/// implementation developed for Gaussian Sum Filter (GSF)
///
/// the state of MC contains a list of single state
/// each component contains its own position, direction and stepSize
/// in the step() method, loop all the single components and do caculating as in
/// the single stepper
/// in the surfaceReached() method, determine if all components
/// are on the aimed surface
/// with targetSurface() method ,collect the candidate surfaces inNavigator with
/// the
/// combination of components (pos,dir),
/// and update the stepSize of each components
///
/// in mcs, each single component owns a status:
/// free - not on surface
/// lock - on surface
/// dead - can not target the surface
///
/// @to do the dead surface should not to kill in the mcs,
/// should be determined in the Gaussian Sum Filter
/// @to do compliant with the MultiParameter class

template <typename BField,
          typename corrector_t     = VoidIntersectionCorrector,
          typename extensionlist_t = StepperExtensionList<DefaultExtension>,
          typename auctioneer_t    = detail::VoidAuctioneer>
class MultiEigenStepper : public EigenStepper<BField>
{

public:
  using cstep     = detail::ConstrainedStep;
  using Corrector = corrector_t;

  /// Jacobian, Covariance and State defintions
  using Jacobian   = ActsMatrixD<5, 5>;
  using Covariance = ActsSymMatrixD<5>;

  /// @note Currently the BoundState/CurvilinearState is defined for Single
  /// component which is a combination of all components, the Jacobian for this
  /// meaning is nonsense
  /// This should be replaced by
  /// std::tuple<MultiBoundParameters,list<Jacobian>,list<double> or other wise
  /// structure with Jacobians
  using BoundState       = std::tuple<BoundParameters, Jacobian, double>;
  using CurvilinearState = std::tuple<CurvilinearParameters, Jacobian, double>;
  using SingleStateType  = typename EigenStepper<BField>::State;

  /// @brief State for track parameter propagation
  ///
  /// by the propagator
  /// the MultiState behave like the SingleState in Navigator, while contains
  /// information of each single components
  struct MultiState
  {

    /// Constructor from the initial track parameters
    /// construct the multi components from one single component
    ///
    /// @param [in] gctx is the context object for the geometry
    /// @param [in] mctx is the context object for the magnetic field
    /// @param [in] par The track parameters at start
    /// @param [in] ndir The navigation direction w.r.t momentum
    /// @param [in] ssize is the maximum step size
    ///
    /// @note the coviance matrix is not in use now
    template <typename parameters_t>
    explicit MultiState(std::reference_wrapper<const GeometryContext>      gctx,
                        std::reference_wrapper<const MagneticFieldContext> mctx,
                        const parameters_t&                                par,
                        NavigationDirection ndir = forward,
                        double ssize = std::numeric_limits<double>::max())
      : pos(par.position())
      , dir(par.momentum().normalized())
      , p(par.momentum().norm())
      , q(par.charge())
      , navDir(ndir)
      , stepSize(ndir * std::abs(ssize))
      , fieldCache(mctx)
      , geoContext(gctx)
    {
      /// initialize the MC state with one component
      stateCol.push_back(
          std::make_tuple(SingleStateType(gctx, mctx, par, ndir, ssize),
                          1.,
                          StateStatus::FREE));
      /// remember the start parameters
      startPos = pos;
      startDir = dir;
    }
    /// the list of <singleState & weight && status>
    std::list<std::tuple<SingleStateType, double, StateStatus>> stateCol;

    /// Global start particle position
    Vector3D startPos = Vector3D(0., 0., 0.);

    /// Momentum start direction (normalized)
    Vector3D startDir = Vector3D(1., 0., 0.);

    /// Global particle position
    Vector3D pos = Vector3D(0., 0., 0.);

    /// Momentum direction (normalized)
    Vector3D dir = Vector3D(1., 0., 0.);

    /// Momentum
    double p = 0.;

    /// The charge
    /// Currently suppose every component has the same charge
    double q = 1.;

    /// Navigation direction, this is needed for searching
    NavigationDirection navDir;

    /// Currently not in used
    bool covTransport = false;

    /// accummulated path length state
    /// the combination of the alive components
    double pathAccumulated = 0.;

    /// the MultiStepper stepSize not used now
    cstep stepSize{std::numeric_limits<double>::max()};

    /// This caches the current magnetic field cell and stays
    /// (and interpolates) within it as long as this is valid.
    /// See step() code for details.
    typename BField::Cache fieldCache;

    /// The geometry context
    std::reference_wrapper<const GeometryContext> geoContext;

    /// List of algorithmic extensions
    extensionlist_t extension;

    /// Auctioneer for choosing the extension
    auctioneer_t auctioneer;

    /// @brief Storage of magnetic field and the sub steps during a RKN4 step
    struct
    {
      /// Magnetic fields
      Vector3D B_first, B_middle, B_last;
      /// k_i of the RKN4 algorithm
      Vector3D k1, k2, k3, k4;
    } stepData;
  };

  /// @brief Global particle position accessor
  /// get the combination of the position of all FREE components
  Vector3D
  position(const MultiState& state) const;

  /// @brief Global direction accessor
  /// get the combination of the direction of all FREE components
  Vector3D
  direction(const MultiState& state) const;

  /// @brief momentum accessor
  /// @brief Global get the combination of the momentum of all FREE components
  double
  momentum(const MultiState& state) const;

  /// Global particle position accessor
  /// get the position of each component
  Vector3D
  position(const SingleStateType& singlestate) const
  {
    return EigenStepper<BField>::position(singlestate);
  }
  /// Global direction accessor
  /// get the direction of each component
  Vector3D
  direction(const SingleStateType& singlestate) const
  {
    return EigenStepper<BField>::direction(singlestate);
  }
  /// Global momentum accessor
  /// get the momentum of each component
  double
  momentum(const SingleStateType& singlestate) const
  {
    return EigenStepper<BField>::momentum(singlestate);
  }

  /// Charge access
  double
  charge(const MultiState& state) const
  {
    return state.q;
  }

  /// always use the state_type in the Propagator
  using state_type = MultiState;

  /// Constructor requires knowledge of the detector's magnetic field
  MultiEigenStepper(BField bField = BField());

  /// @copy of the Single stepper
  /// @brief getField from single-component or multi-component
  /// @param [in,out] state is the propagation state associated with the track
  ///                 the magnetic field cell is used (and potentially updated)
  /// @param [in] pos is the field position

  template <typename stateType>
  Vector3D
  getField(stateType& state, const Vector3D& pos) const
  {
    // get the field from the cell
    return m_bField.getField(pos, state.fieldCache);
  }

  /// @brief Tests if all the single states reach a surface
  /// if all single states reach(except the dead ones) successfully,
  ///  return true; otherwise return false;
  ///
  /// the single state that successfully reach is set to locked
  /// then if all single states locked, free all of them
  ///
  /// @param [in] state State is the mcs state
  /// @param [in] surface Surface that is tested
  ///
  /// @return Boolean statement if surface is reached by state
  bool
  surfaceReached(MultiState& state, const Surface* surface) const;

  /// output template method to check
  void
  outPut(const MultiState& state) const;

  /// @brief this caculates the stepSize of all the components
  /// to the candidate surfaces/layers/boundaries in the Navigator
  /// @param [in] state State is the mcs
  /// @param [in] surface Surface that is tested
  /// @param [in] navigator options
  /// @param [in] navigator corrections
  ///
  /// @return pair of <Boolean,Double>
  /// Bool: at least one of the single components target successfully
  /// Double : the smallest distance of the component pathlengh in the compact
  //
  /// @to do: the dead component should be taken into consideration in the
  /// Fitter stage
  template <typename options_t>
  std::pair<bool, double>
  targetSurface(MultiState&      state,
                const Surface*   surface,
                const options_t& navOpts,
                const Corrector& navCorr) const;

  /// reweight the free components
  /// the free and locked components are reweighted
  void
  normalizeComponents(MultiState& state) const;

  /// @brief The method s to delete the died components
  /// acturally, this should be done at Fitting step
  void
  deleteComponents(MultiState& state) const;

  /// @brief get a sinlge parameter of combination of multi component on a
  /// surface, the jocobian is nonsence here
  /// the pathAccumulated is a combination calculated in the step()
  /// @to do parameter should be MultiBoundParameters
  BoundState
  boundState(MultiState&    state,
             const Surface& surface,
             bool           reinitialize = true) const;

  /// @brief get a sinlge parameter of combination of multi component on a
  /// surface, the jocobian is nonsence here
  /// the pathAccumulated is a combination calculated in the step()
  /// @to do parameter should be MultiCurvilinearParameters
  CurvilinearState
  curvilinearState(MultiState& state, bool reinitialize = true) const;

  /// Return a corrector
  Corrector
  corrector(MultiState& state) const
  {
    return Corrector(state.startPos, state.startDir, state.pathAccumulated);
  }

  /// updateStep method
  /// only call at navigator, use aborter value to udpate
  /// use a udpate each stepsize with the combination value ?
  void
  updateStep(MultiState&      state,
             const Corrector& navCorr,
             double           navigationStep,
             bool             release = false) const;

  void
  releaseStep(MultiState& state, cstep::Type type = cstep::actor) const;

  /// this called in StandardAborter, set a pathlimit of the combination
  /// component
  void
  updateStep(MultiState& state,
             double      abortStep,
             cstep::Type type = cstep::aborter) const;

  /// @brief update for the single state, update singlestate to some parameters
  void
  update(SingleStateType& singlestate, const BoundParameters& pars) const;

  /// @brief update for the single state, update singlestate direction and p
  void
  update(SingleStateType& singlestate,
         const Vector3D&  uposition,
         const Vector3D&  udirection,
         double           up) const;

  /// Perform a Runge-Kutta track parameter propagation step
  ///
  /// @param [in,out] state is the propagation state associated with the track
  /// parameters that are being propagated.
  ///
  ///                      the state contains the desired step size.
  ///                      It can be negative during backwards track
  ///                      propagation,
  ///                      and since we're using an adaptive algorithm, it can
  ///                      be modified by the stepper class during propagation.
  template <typename propagator_state_t>
  Result<double>
  step(propagator_state_t& state) const;

private:
  /// Magnetic field inside of the detector
  BField m_bField;
};
}  // namespace Acts

#include "Acts/Propagator/MultiEigenStepper.ipp"
