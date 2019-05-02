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
class MultiCurvilinearTrackParameters
    : public MultiTrackParameters<ChargePolicy>
{
public:
  /// pair of <weight,point_to_track>, when insert a trackParameter, that will allow to sort the track by weight 
  using WeightedTrackPars = std::pair< double, std::unique_ptr<TrackParametersBase> >;
  using TrackParMap       = std::map<WeightedTrackPars ,const unsigned int, std::greater<WeightedTrackPars> >;
  using TrackParMapConstIter   = TrackParMap::const_iterator;
  using TrackParMapIter   = TrackParMap::iterator;
  struct TrackIndexFinder
  {
	TrackIndexFinder(const unsigned int id) : m_index(id){}
	bool operator() (const TrackParMap::value_type& trackPar)
	{
	  return trackPar.second == m_index;
	}
	const unsigned int m_index;
  };

  template <typename T = ChargePolicy,
            std::enable_if_t<std::is_same<T, ChargedPolicy>::value, int> = 0>
  MultiCurvilinearTrackParameters(double weight,TrackParametersBase* pTrackBase)
    : MultiTrackParameters<ChargePolicy>(weight,pTrackBase){}
//     m_upSurface(Surface::makeShared<PlaneSurface>(pTrackBase->position(), pTrackBase->momentum() )) {}

  template <typename T = ChargePolicy,
            std::enable_if_t<std::is_same<T, NeutralPolicy>::value, int> = 0>
  MultiCurvilinearTrackParameters(double weight,TrackParametersBase* pTrackBase)
    : MultiTrackParameters<ChargePolicy>(weight,pTrackBase){}
//     m_upSurface(Surface::makeShared<PlaneSurface>(pTrackBase->position(), pTrackBase->momentum() )) {}

  MultiCurvilinearTrackParameters() = default;

  /// @brief copy constructor - charged/neutral
  /// @param[in] copy The source parameters
  /*unused*/
  MultiCurvilinearTrackParameters(
      const MultiCurvilinearTrackParameters<ChargePolicy>& copy) 
	: MultiTrackParameters<ChargePolicy>(copy),
	m_upSurface(copy.m_upSurface)
  {
  }

  /// @brief move constructor - charged/neutral
  /// @param[in] other The source parameters
  /*unused*/
  MultiCurvilinearTrackParameters(
      MultiCurvilinearTrackParameters<ChargePolicy>&& other)
    : MultiTrackParameters<ChargePolicy>(std::move(other))
    , m_upSurface(std::move(other.m_upSurface))
  {
  }

  /// @brief desctructor
  ~MultiCurvilinearTrackParameters() override = default;

  /// @brief move assignment operator - charged/netural
  /// virtual constructor for type creation without casting
  /*unused*/
  MultiCurvilinearTrackParameters<ChargePolicy>&
  operator=(MultiCurvilinearTrackParameters<ChargePolicy>&& rhs)
  {
    // check for self-assignment
    if (this != &rhs) {
      MultiTrackParameters<ChargePolicy>::operator=(std::move(rhs));
      m_upSurface                                  = std::move(rhs.m_upSurface);
    }
    return *this;
  }

  // make reference surface with the combination pos/mom
  void makeReferenceSurface()
  {
     m_upSurface = Surface::makeShared<PlaneSurface>( this->position(), this->momentum() ) ;
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
  set(const GeometryContext& gctx, ParValue_t newValue, unsigned int id)
  {
    // set the parameter & update the new global position
    this->getParameterSet(id).template setParameter<par>(newValue);
    //this->updateGlobalCoordinates(gctx, typename par_type<par>::type(), id, this->referenceSurface(id));
    this->updateGlobalCoordinates(gctx, typename par_type<par>::type(), id);
	std::cout<<"in set pos dir "<<this->position(id)<<" "<<this->momentum(id).normalized()<<std::endl;
//	auto pos = this->position(id);
//	auto dir = this->momentum(id).normalized();
	// update the referece surface with the aimed parameter 
	this->updateReferenceSurface( this->position(id), this->momentum(id).normalized(), id);
    // recreate the surface
//    m_upSurface = Surface::makeShared<PlaneSurface>(
//        this->position(), this->momentum().normalized());
    // reset to (0,0)
    this->getParameterSet(id).template setParameter<par>(0.);
	// update the referece surface with the combination
	this->updateReferenceSurface( this->position(), this->momentum().normalized());
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
  set(const GeometryContext& gctx,ParValue_t newValue, unsigned int id)
  {
    this->getParameterSet(id).template setParameter<par>(newValue);
    this->updateGlobalCoordinates(gctx, typename par_type<par>::type(), id );
	// update the referece surface with the aimed parameter 
	this->updateReferenceSurface( this->position(id), this->momentum(id).normalized(), id);
    // recreate the surface
//    m_upSurface = Surface::makeShared<PlaneSurface>(
//        this->position(id), this->momentum(id).normalized());
	// update the referece surface with the combination
	this->updateReferenceSurface( this->position(), this->momentum().normalized());
  }
  
  /// @brief access to the reference surface
  const Surface&
  referenceSurface() const final
  {
    return *m_upSurface;
  }

  /// @brief access to the reference surface
  const Surface&
  referenceSurface(unsigned int id) const  final
  {
	TrackParMapConstIter it = std::find_if( this->m_TrackList.begin(), this->m_TrackList.end(), TrackIndexFinder(id));
	assert( it != this->m_TrackList.end() );
	return (*it).first.second->referenceSurface();
  }

  /// @brief update the reference surface
  void 
  updateReferenceSurface( const ActsVectorD<3>& pos, const ActsVectorD<3>& dir) final
  {
	m_upSurface = Surface::makeShared<PlaneSurface>(pos,dir);
  }
  void 
  updateReferenceSurface( const ActsVectorD<3>& pos, const ActsVectorD<3>& dir, unsigned int id) 
  {
	TrackParMapConstIter it = std::find_if( this->m_TrackList.begin(), this->m_TrackList.end(), TrackIndexFinder(id));
	assert( it != this->m_TrackList.end() );
	(*it).first.second->updateReferenceSurface( pos, dir );
  }

  /// @brief access to the measurement frame, i.e. the rotation matrix with
  /// respect to the global coordinate system, in which the local error
  /// is described.
  ///
  /// For a curvilinear track parameterisation this is identical to the
  /// rotation matrix of the intrinsic planar surface.
  RotationMatrix3D
  referenceFrame(const GeometryContext& gctx) const 
  {
    return m_upSurface->transform(gctx).linear();
  }
/*
  RotationMatrix3D
  referenceFrame(const GeometryContext& gctx, const unsigned int id) const 
  {
	TrackParMapIter it = std::find_if( this->m_TrackList.begin(), this->m_TrackList.end(), TrackIndexFinder(id));
	assert( it != this->m_TrackList.end() );
	return (*it).first.second->referenceSurface().transform(gctx).linear();
  }
  */

private:
  std::shared_ptr<PlaneSurface> m_upSurface;
};
}  // namespace Acts
