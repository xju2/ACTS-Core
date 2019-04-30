// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#include <memory>
#include "Acts/EventData/MultiTrackParameters.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"

namespace Acts {

template <typename ChargePolicy>
class MultiBoundTrackParameters
    : public MultiTrackParameters<ChargePolicy>
{
public:
  /// type of covariance matrix
  using CovPtr_t = typename MultiTrackParameters<ChargePolicy>::CovPtr_t;

  template <typename T = ChargePolicy,
            std::enable_if_t<std::is_same<T, ChargedPolicy>::value, int> = 0>
  MultiBoundTrackParameters(double weight,
	  TrackParametersBase* pTrackBase,
	  std::shared_ptr<const Surface> surface )
    : MultiTrackParameters<ChargePolicy>(weight,pTrackBase),
     m_pSurface(std::move(surface))
  {
	assert(m_pSurface);
  }

  template <typename T = ChargePolicy,
            std::enable_if_t<std::is_same<T, NeutralPolicy>::value, int> = 0>
  MultiBoundTrackParameters(double weight,
	  TrackParametersBase* pTrackBase,
	  std::shared_ptr<const Surface> surface)
    : MultiTrackParameters<NeutralPolicy>(weight,pTrackBase),
	m_pSurface(std::move(surface))
  {
	assert(m_pSurface);
  }


  /// @brief copy constructor - charged/neutral
  /// @param[in] copy The source parameters
  MultiBoundTrackParameters(
      const MultiBoundTrackParameters<ChargePolicy>& copy) 
	: MultiTrackParameters<ChargePolicy>(copy),
	m_pSurface(copy.m_pSurface)
  {
  }

  /// @brief move constructor - charged/neutral
  /// @param[in] other The source parameters
  MultiBoundTrackParameters(
      MultiBoundTrackParameters<ChargePolicy>&& other)
    : MultiTrackParameters<ChargePolicy>(std::move(other))
    , m_pSurface(std::move(other.m_pSurface))
  {
  }

  /// @brief desctructor
  ~MultiBoundTrackParameters() override = default;


  /// @brief move assignment operator - charged/netural
  /// virtual constructor for type creation without casting
  MultiBoundTrackParameters<ChargePolicy>&
  operator=(MultiBoundTrackParameters<ChargePolicy>&& rhs)
  {
    // check for self-assignment
    if (this != &rhs) {
      MultiTrackParameters<ChargePolicy>::operator=(std::move(rhs));
      m_pSurface                                  = std::move(rhs.m_pSurface);
    }
    return *this;
  }

  /// @brief clone - charged/netural
  MultiTrackParameters<ChargePolicy>*
  clone() const override  /*unused*/
  {
    return new MultiBoundTrackParameters<ChargePolicy>(*this);
  }

  
  /// @brief update of the track parameterisation
  /// only possible on non-const objects, enable for local parameters
  ///
  /// @tparam ParID_t The parameter type
  ///
  /// @param newValue The new updaed value
  ///
  /// For curvilinear parameters the local parameters are forced to be
  /// (0,0), hence an update is an effective shift of the reference
  template <ParID_t par,
            std::enable_if_t<std::is_same<typename par_type<par>::type,
                                          local_parameter>::value,
                             int> = 0>
  void
  set(const GeometryContext& gctx, ParValue_t newValue, unsigned int order)
  {
    // set the parameter & update the new global position
    this->getParameterSet(order).template setParameter<par>(newValue);
    this->updateGlobalCoordinates(gctx, typename par_type<par>::type(),order);
    // recreate the surface
    m_pSurface = Surface::makeShared<PlaneSurface>(
        this->position(order), this->momentum(order).normalized());
    // reset to (0,0)
    this->getParameterSet(order).template setParameter<par>(0.);
  }

  /// @brief update of the track parameterisation
  /// only possible on non-const objects
  /// enable for parameters that are not local parameters
  /// @tparam ParID_t The parameter type
  ///
  /// @param newValue The new updaed value
  ///
  /// For curvilinear parameters the directional change of parameters
  /// causes a recalculation of the surface
  template <ParID_t par,
            std::enable_if_t<not std::is_same<typename par_type<par>::type,
                                              local_parameter>::value,
                             int> = 0>
  void
  set(const GeometryContext& gctx, ParValue_t newValue,unsigned int order)
  {
    this->getParameterSet(order).template setParameter<par>(newValue);
    this->updateGlobalCoordinates(gctx, typename par_type<par>::type(),order);
    // recreate the surface
    m_pSurface = Surface::makeShared<PlaneSurface>(
        this->position(order), this->momentum(order).normalized());
  }
  

  /// @brief access to the reference surface
  const Surface&
  referenceSurface() const final
  {
    return *m_pSurface;
  }

  /// @brief access to the measurement frame, i.e. the rotation matrix with
  /// respect to the global coordinate system, in which the local error
  /// is described.
  ///
  /// For a curvilinear track parameterisation this is identical to the
  /// rotation matrix of the intrinsic planar surface.
  RotationMatrix3D
  referenceFrame(const GeometryContext& gctx) const final
  {
    return std::move(
        m_pSurface->referenceFrame(gctx, this->position(), this->momentum()));
  }

private:
  std::shared_ptr<const Surface> m_pSurface;
};
}  // namespace Acts
