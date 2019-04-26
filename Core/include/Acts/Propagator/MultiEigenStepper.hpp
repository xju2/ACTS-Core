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

/// @brief multicomponent(MC) stepper is based on the Single Stepper
/// implementation
/// developed fot Gaussian Sum Filter (GSF)
/// the state of MC contains a list of single state
/// each component contains its own position, direction and stepSize
/// in the step() method, loop all the single components and do caculating as in
/// the single object
/// in the navigator, with surfaceReach method to determine if all components
/// are on the aimed surface
/// with targetSurface() method ,collect the candidate surfaces with the
/// combination of components,
/// and update the stepSize of each components
///
/// in MC, each single component owns a status:
/// free - not on surface
/// lock - on surface
/// dead - can not target surface
///
/// @to do the dead surface should not to kill in the MC,
/// should be determined in the Gaussian Sum Filter
/// @to do the par&cov in MultiTrackState

template <typename BField,
          typename corrector_t     = VoidIntersectionCorrector,
          typename extensionlist_t = StepperExtensionList<DefaultExtension>,
          typename auctioneer_t    = detail::VoidAuctioneer>
class MultiEigenStepper : public EigenStepper<BField>
{

private:
  // This struct is a meta-function which normally maps to BoundParameters...
  template <typename T, typename S>
  struct s
  {
    using type = BoundParameters;
  };

  // ...unless type S is int, in which case it maps to Curvilinear parameters
  template <typename T>
  struct s<T, int>
  {
    using type = CurvilinearParameters;
  };

public:
  using cstep = detail::ConstrainedStep;

  /// Jacobian, Covariance and State defintions
  using Jacobian         = ActsMatrixD<5, 5>;
  using Covariance       = ActsSymMatrixD<5>;
  using BoundState       = std::tuple<BoundParameters, Jacobian, double>;
  using CurvilinearState = std::tuple<CurvilinearParameters, Jacobian, double>;
  using SingleStateType  = typename EigenStepper<BField>::State;

  /// @brief State for track parameter propagation
  ///
  /// by the propagator
  struct MultiState
  {

    /// Constructor from the initial track parameters
    /// @param [in] par The track parameters at start
    /// @param [in] ndir The navigation direction w.r.t momentum
    /// @param [in] ssize is the maximum step size
    ///
    /// @note the covariance matrix is copied when needed
    /// the coviance matrix is not in use now
    template <typename parameters_t>
	explicit MultiState(std::reference_wrapper<const GeometryContext>      gctx,
		std::reference_wrapper<const MagneticFieldContext> mctx,
				   const parameters_t&            par,
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
      stateCol.push_back(std::make_tuple(
          SingleStateType(gctx, mctx, par, ndir, ssize), 1., StateStatus::FREE));
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
    double q = 1.;

    /// Navigation direction, this is needed for searching
    NavigationDirection navDir;

    /// no use and meaningless, should delete
    ActsMatrixD<5, 5> jacobian = ActsMatrixD<5, 5>::Identity();
    bool       covTransport = false;
    Covariance cov          = Covariance::Zero();

    /// accummulated path length state
    double pathAccumulated = 0.;

    /// adaptive step size of the runge-kutta integration
	/// should reserve it as the combination ?
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
  position(const MultiState& state) const
  {
    Vector3D pos = Vector3D(0, 0, 0);
    for (const auto& tuple_state : state.stateCol) {
      if (std::get<2>(tuple_state) != StateStatus::FREE) continue;
      Vector3D component_pos
          = EigenStepper<BField>::position(std::get<0>(tuple_state));
      pos += component_pos * std::get<1>(tuple_state);
    }
    return pos;
  }

  /// Global direction accessor
  /// get the combination of the Direction of all FREE components
  Vector3D
  direction(const MultiState& state) const
  {
    Vector3D dir = Vector3D(0, 0, 0);
    for (const auto& tuple_state : state.stateCol) {
      if (std::get<2>(tuple_state) != StateStatus::FREE) continue;
      Vector3D component_dir
          = EigenStepper<BField>::direction(std::get<0>(tuple_state));
      dir += component_dir * std::get<1>(tuple_state);
    }
    return dir;  
  }

  /// Global get the combination of the Momentum of all FREE components
  /// Actual momentum accessor
  double
  momentum(const MultiState& state) const
  {
    double mom = 0.;
    for (const auto& tuple_state : state.stateCol) {
      if (std::get<2>(tuple_state) != StateStatus::FREE) continue;
      double component_mom
          = EigenStepper<BField>::momentum(std::get<0>(tuple_state));
      mom += component_mom * std::get<1>(tuple_state);
    }
    return mom;
  }

  /// Global particle position accessor
  /// get the position of each component
  Vector3D
  position(const SingleStateType& singlestate) const
  {
    Vector3D pos = EigenStepper<BField>::position(singlestate);
    return pos;
  }
  /// Global direction accessor
  /// get the direction of each component
  Vector3D
  direction(const SingleStateType& singlestate) const
  {
    Vector3D dir = EigenStepper<BField>::direction(singlestate);
    return dir;
  }
  /// Global Momentum accessor
  /// get the Momentum of each component
  double
  momentum(const SingleStateType& singlestate) const
  {
    double mom = EigenStepper<BField>::momentum(singlestate);
    return mom;
  }

  /// Charge access
  double
  charge(const MultiState& state) const
  {
    return state.q;
  }

  /// always use the state_type in the Propagator 
  using state_type = MultiState;

  /// Return parameter types depend on the propagation mode:
  /// - when propagating to a surface we usually return BoundParameters
  /// - otherwise CurvilinearParameters
  template <typename T, typename S = int>
  using return_parameter_type = typename s<T, S>::type;

  /// Constructor requires knowledge of the detector's magnetic field
  MultiEigenStepper(BField bField = BField());

  /// Get the field for the stepping, it checks first if the access is still
  /// within the Cell, and updates the cell if necessary.
  ///
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

  /// Tests if all the single states reached a surface
  /// if all reached (except the dead ones), return true, otherwise return
  /// false;
  /// the single status that returned is set to locked
  /// if all single states locked, free all of them
  ///
  /// @param [in] state State is the MC state
  /// @param [in] surface Surface that is tested
  ///
  /// @return Boolean statement if surface is reached by state
  bool
  surfaceReached(MultiState& state, const Surface* surface) const
  {
    bool status = true;
    for (auto& tuple_state : state.stateCol) {
      auto& singlestate = std::get<0>(tuple_state);
      /// if locked, test if there is wrong
      if (std::get<2>(tuple_state) != StateStatus::FREE) {

        if (std::get<2>(tuple_state) == StateStatus::DEAD) {
            //std::cout<<" the component is dead "<<std::endl;
        } else {
            //std::cout<<" the component is locked "<<std::endl;
        }
      } else {
        /// not on surface, states set false
        if (!surface->isOnSurface(state.geoContext, EigenStepper<BField>::position(singlestate),
                                  EigenStepper<BField>::direction(singlestate),
                                  true)) {
          status = false;
        }
        // on surface: locked
        else {
          std::get<2>(tuple_state) = StateStatus::LOCKED;
        }
      }
    }
    if (status == false) {
      return false;
    } else {
      // if all the components are on surface(except the dead ones, set them
      // free and set the currentSurface in Navigator.
      for (auto& tuple_state : state.stateCol) {
        if (std::get<2>(tuple_state) == StateStatus::LOCKED) {
          std::get<2>(tuple_state) = StateStatus::FREE;
        }
      }
      return true;
    }
  }

  /// output template method to check
	void outPut(const MultiState& state) const
	{
	  for( const auto& tuple_state : state.stateCol )
	  {
		const auto& singlestate = std::get<0>(tuple_state);
		const auto& weight 		= std::get<1>(tuple_state);
		const auto& stat 		= std::get<2>(tuple_state);
		std::cout<<"the single component stat: "<<stat<<" weight "<<weight<<std::endl;
		std::cout<<"pos: "<<singlestate.pos<<std::endl;
		std::cout<<"dir: "<<singlestate.dir<<std::endl;
	  }
	}

  /// this caculates all the components the stepSize
  /// to the candidate surfaces/layers/boundaries in the Navigator
  /// @param [in] state State is the MC state
  /// @param [in] surface Surface that is tested
  /// @param [in] navigator options
  /// @param [in] navigator corrections
  ///
  /// @return pair of <Boolean,Double>
  /// Bool: at least one of the single components target successfully
  /// Double : the smallest distance of the component pathlengh in the compact
  //
  /// @to do: the dead component should be take into more consideration in the
  /// Fitter
  template <typename options_t>
  std::pair<bool, double>
  targetSurface(MultiState&        state,
                const Surface*     surface,
                const options_t&   navOpts,
                const corrector_t& navCorr) const
  {
    bool   active  = false;
    double minDist = std::numeric_limits<double>::max();
    for (auto& tuple_state : state.stateCol) {
      /// only target the free component
      if (std::get<2>(tuple_state) != StateStatus::FREE) continue;
      auto& singlestate = std::get<0>(tuple_state);
      auto  target      = EigenStepper<BField>::targetSurface(
          singlestate, surface, navOpts, navCorr);
      minDist = minDist < target.second ? minDist : target.second;
      /// a protection to avoid the flip components, which is abnormal in
      if (direction(state).dot(direction(singlestate)) < 0) {
        target.first = false;
      }
      /// this should be considered in the Fitter, currently we simply set
      /// it dead here
      if (target.first == false) {
        std::get<2>(tuple_state) = StateStatus::DEAD;
      }
      /// as long as one component is alive, return true
      else {
        active = true;
      }
    }
    deleteDied(state);
    normalize(state);
    return std::make_pair(active, minDist);
  }

  /// reweight the free components
  /// the free and locked components are reweighted
  void
  normalize(MultiState& state) const
  {
    double weight_sum = 0;
    for (auto& tuple_state : state.stateCol) {
      if (std::get<2>(tuple_state) == StateStatus::DEAD) continue;
      weight_sum += std::get<1>(tuple_state);
    }
    for (auto& tuple_state : state.stateCol) {
      if (std::get<2>(tuple_state) == StateStatus::DEAD) continue;
      std::get<1>(tuple_state) = std::get<1>(tuple_state) / weight_sum;
    }
  }
  /// @brief The method s to delete the died components
  /// acturally, this should be done at Fitting step
  void
  deleteDied(MultiState& state) const
  {
    auto& col = state.stateCol;
    typename std::list<std::tuple<SingleStateType, double, StateStatus>>::
        iterator it
        = col.begin();
    while (it != col.end()) {
      if (std::get<2>(*it) == StateStatus::DEAD) {
        it = col.erase(it);
      } else {
        it++;
      }
    }
  }

  /// @brief get a sinlge parameter of combination of multi component on a surface, the jocobian is nonsence here 
  /// the pathAccumulated is a combination calculated in the step()
  /// @to do parameter should be MultiBoundParameters
  BoundState
  boundState(MultiState&    state,
             const Surface& surface,
             bool           reinitialize = true) const;

  /// @brief get a sinlge parameter of combination of multi component on a surface, the jocobian is nonsence here 
  /// the pathAccumulated is a combination calculated in the step()
  /// @to do parameter should be MultiCurvilinearParameters
  CurvilinearState
  curvilinearState(MultiState& state, bool reinitialize = true) const;

  /// Return a corrector
  corrector_t
  corrector(MultiState& state) const
  {
    return corrector_t(state.startPos, state.startDir, state.pathAccumulated);
  }

  /// updateStep method
  /// only call at navigator, use aborter value to udpate
  /// use a udpate each stepsize with the combination value ?
  template <typename state_type>
  void
  updateStep(state_type&        state,
             const corrector_t& navCorr,
             double             navigationStep,
             bool               release = false) const
  {
    for (auto& tuple_state : state.stateCol) {
      // only deal with the free
      if (std::get<2>(tuple_state) != StateStatus::FREE) continue;
      auto& singlestate = std::get<0>(tuple_state);
      // singlestate.stepSize.update(navigationStep, cstep::actor, release);
      EigenStepper<BField>::updateStep(
          singlestate, navCorr, navigationStep, release);
      if (singlestate.pathAccumulated == 0. and navCorr(singlestate.stepSize)) {
        /*dummy*/
      }
    }
  }

  /// this called in StandardAborter, set a pathlimit of the combination
  /// component
  template <typename state_type>
  void
  updateStep(state_type& state,
             double      abortStep,
             cstep::Type type = cstep::aborter) const
  {
    for (auto& tuple_state : state.stateCol) {
      // only deal with the free
      if (std::get<2>(tuple_state) != StateStatus::FREE) continue;
      auto& singlestate = std::get<0>(tuple_state);
      singlestate.stepSize.update(abortStep, type);
    }
  }


  /// @brief update for the single state, update singlestate to some parameters
  void
  update(SingleStateType&          singlestate,
	  const BoundParameters& pars) const;

  /// @brief update for the single state, update singlestate direction and p
  void
  update(SingleStateType&          singlestate,
         const Vector3D& uposition,
         const Vector3D& udirection,
         double          up) const;

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
