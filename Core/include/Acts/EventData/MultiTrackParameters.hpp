// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#include <type_traits>
#include "Acts/EventData/TrackParametersBase.hpp"
#include "Acts/EventData/detail/coordinate_transformations.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/GeometryContext.hpp"
#include <set>

namespace Acts {

/// @class MultiTrackParameters
///
/// @brief base class for a multi set of track parameters
///
/// This class implements the interface for charged/neutral track parameters for
/// the case that it represents a multi set of track parameters (used in Gsf or MultiTrack Fitters)

  /// pair of <weight,point_to_track>, when insert a trackParameter, that will allow to sort the track by weight 
  using WeightedTrackPars = std::pair< double, std::unique_ptr<TrackParametersBase> >;
  using TrackParSet       = std::multiset< WeightedTrackPars, std::greater<WeightedTrackPars> >;
///
/// The track parameters and their uncertainty are defined in local reference
/// frame which depends on the associated surface of the track parameters.
///
/// @tparam ChargePolicy type for distinguishing charged and neutral
/// tracks/particles
///         (must be either ChargedPolicy or NeutralPolicy)
template <class ChargePolicy>
class MultiTrackParameters : public TrackParametersBase
{
  static_assert(std::is_same<ChargePolicy, ChargedPolicy>::value
                    or std::is_same<ChargePolicy, NeutralPolicy>::value,
                "ChargePolicy must either be 'Acts::ChargedPolicy' or "
                "'Acts::NeutralPolicy");

public:
  // public typedef's

  /// vector type for stored track parameters
  using ParVector_t = typename TrackParametersBase::ParVector_t;

  /// type of covariance matrix
  using CovMatrix_t = typename TrackParametersBase::CovMatrix_t;

  /// type for unique pointer to covariance matrix
  using CovPtr_t = std::unique_ptr<const CovMatrix_t>;

  /// @brief default virtual destructor
  ~MultiTrackParameters() override = default;

  /// @brief get weighted combination of position
  ActsVectorD<3>
  position() const final
  {
	Vector3D pos = Vector3D(0, 0, 0);
	for(const auto& trackWeightedPair :  m_TrackList) {
		pos += trackWeightedPair.first * trackWeightedPair.second->position();
	}
	return pos;
  }

  /// @brief get the position of a single parameter
  ActsVectorD<3> 
  position(unsigned int order) const 
  {
//	std::assert( order < m_tracksPars.size() );
	unsigned int i = 0 ;
	for(const auto& trackWeightedPair :  m_TrackList) {
	  if( i == order ) return trackWeightedPair.second->position();
	  ++i;
	}
  }

  /// @copydoc TrackParametersBase::momentum
  ActsVectorD<3>
  momentum() const final
  {
	Vector3D mom = Vector3D(0, 0, 0);
	for(const auto& trackWeightedPair :  m_TrackList) {
		mom += trackWeightedPair.first * trackWeightedPair.second->momentum();
	}
	return mom;
  }

  ///  get the single mom of parameter
  ActsVectorD<3> 
  momentum(unsigned int order) const 
  {
//	std::assert( order < m_tracksPars.size() );
	unsigned int i = 0 ;
	for(const auto& trackWeightedPair :  m_TrackList) {
	  if( i == order ) return trackWeightedPair.second->momentum();
	  ++i;
	}
  }

  /// @brief equality operator
  ///
  /// @return @c true if all single track parameters are same
  bool
  operator==(const TrackParametersBase& rhs) const override
  {
    auto casted = dynamic_cast<decltype(this)>(&rhs);
    if (!casted) {
      return false;
    }

    return ( m_TrackList == casted->m_TrackList);
  }

  /// @copydoc TrackParametersBase::charge
  double
  charge() const final
  {
    return m_TrackList.begin()->second->charge();
  }

  /// @copydoc TrackParametersBase::getParameterSet
  /// currently get first component parameters
  /// @to do : return combination
  /// const
  const FullParameterSet&
  getParameterSet() const final
  {
    return m_TrackList.begin()->second->getParameterSet();
  }
  /// copy base class
  /// currently get first component parameters
  /// @to do : return combination
  /// writable  
  virtual FullParameterSet&
  getParameterSet() final
  {
    return m_TrackList.begin()->second->getParameterSet();
  }

  /// @brief get single component
  const FullParameterSet& 
  getParameterSet(unsigned int order) const 
  {
//	std::assert( order < m_tracksPars.size() );
	unsigned int i = 0 ;
	for( const auto& weightedTrackPars :  m_TrackList) {
	  if( i == order ) return weightedTrackPars.second->getParameterSet();
	  ++i;
	}
  }

  /// @brief get single component
  virtual FullParameterSet& 
  getParameterSet(unsigned int order) 
  {
//	std::assert( order < m_tracksPars.size() );
	unsigned int i = 0 ;
	for( auto& weightedTrackPars :  m_TrackList) {
	  if( i == order ) return weightedTrackPars.second->getParameterSet();
	  ++i;
	}
  }


  /// @brief get single component
  ParVector_t
  parameters(unsigned int order) const
  {
    return getParameterSet(order).getParameters();
  }

  /// @brief return weight
  double weight(unsigned int order) const
  {
//	std::assert( order < m_tracksPars.size() );
	unsigned int i = 0 ;
	for( auto& weightedTrackPars :  m_TrackList) {
	  if( i == order ) return weightedTrackPars.first;
	  ++i;
	}
  }

  /// @brief update mom pos 
  void 
  updateMom( ActsVectorD<3>& /*unused*/){} 

  void 
  updatePos( ActsVectorD<3>& /*unused*/ ){} 

  /// @brief append a component to the track list
  void append(double weight, TrackParametersBase* pTrackBase)
  {
	m_TrackList.insert( WeightedTrackPars( weight,std::move(std::unique_ptr<TrackParametersBase>(pTrackBase) ) ) );
  }

protected:
  /// @brief standard constructor for track parameters of charged particles
  ///
  /// @param cov unique pointer to covariance matrix (nullptr is accepted)
  /// @param parValues vector with parameter values
  /// @param position 3D vector with global position
  /// @param momentum 3D vector with global momentum
  template <typename T = ChargePolicy,
            std::enable_if_t<std::is_same<T, ChargedPolicy>::value, int> = 0>
  MultiTrackParameters( double weight, TrackParametersBase* pTrackBase )
    : TrackParametersBase()
  {
	m_TrackList.insert( WeightedTrackPars( weight,std::move(std::unique_ptr<TrackParametersBase>(pTrackBase) ) ) );
  }

  /// @brief standard constructor for track parameters of neutral particles
  ///
  /// @param cov unique pointer to covariance matrix (nullptr is accepted)
  /// @param parValues vector with parameter values
  /// @param position 3D vector with global position
  /// @param momentum 3D vector with global momentum
  template <typename T = ChargePolicy,
            std::enable_if_t<std::is_same<T, NeutralPolicy>::value, int> = 0>
  MultiTrackParameters( double weight, TrackParametersBase* pTrackBase )
    : TrackParametersBase()
  {
	m_TrackList.insert( WeightedTrackPars( weight,std::move(std::unique_ptr<TrackParametersBase>(pTrackBase) ) ) );
  }

  /// @brief default copy constructor
  /*unused*/
  MultiTrackParameters(const MultiTrackParameters<ChargePolicy>& copy)
      = default;

  /// @brief default move constructor
  MultiTrackParameters(MultiTrackParameters<ChargePolicy>&& copy) = default;

  /// @brief copy assignment operator
  ///
  /// @param rhs object to be copied
  /*unused*/
  MultiTrackParameters<ChargePolicy>&
  operator=(const MultiTrackParameters<ChargePolicy>& rhs)  {}

  /// @brief move assignment operator
  ///
  /// @param rhs object to be movied into `*this`
  MultiTrackParameters<ChargePolicy>&
  operator=(MultiTrackParameters<ChargePolicy>&& rhs)
  {
    // check for self-assignment
    if (this != &rhs) {
      m_TrackList = std::move(rhs.m_trackPars);
    }

    return *this;
  }

  /// @brief update global momentum from current parameter values
  ///
  ///
  /// @param[in] gctx is the Context object that is forwarded to the surface
  ///            for local to global coordinate transformation
  ///
  /// @note This function is triggered when called with an argument of a type
  ///       different from Acts::local_parameter
  template <typename T>
  void
  updateGlobalCoordinates(const GeometryContext& /*gctx*/, const T& /*unused*/, unsigned int order)
  {
    auto vMomentum = detail::coordinate_transformation::parameters2globalMomentum(
        getParameterSet(order).getParameters());
	unsigned int i = 0 ;
	for( auto& weightedTrackPars :  m_TrackList) {
	  if( i == order ) weightedTrackPars.second->updateMom( vMomentum );
	  ++i;
	}
  }

  /// @brief update global position from current parameter values
  ///
  /// @note This function is triggered when called with an argument of a type
  /// Acts::local_parameter
  void
  updateGlobalCoordinates(const GeometryContext& gctx,
                          const local_parameter& /*unused*/,
						  unsigned int order)
  {
    auto vPosition = detail::coordinate_transformation::parameters2globalPosition(
        gctx, getParameterSet(order).getParameters(), this->referenceSurface());
	unsigned int i = 0 ;
	for( auto& weightedTrackPars :  m_TrackList) {
	  if( i == order ) weightedTrackPars.second->updatePos( vPosition);
	  ++i;
	}
  }

  TrackParSet 	m_TrackList;
};
}  // namespace Acts
