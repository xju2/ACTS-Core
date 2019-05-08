// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#include "Acts/Extrapolator/detail/ComponentDistance.hpp"

namespace Acts {

  namespace detail{
  /// The component reduction struct
  /// This is plugin to the Propagator that 
  /// performs to reduct the number component
struct ComponentReduction
{
  struct this_result
  {
	std::vector<double> distance;
  };
  detail::KullbackLeiblerComponentDistance klDist;

  using result_type = this_result;
  /// Jacobian, Covariance and State defintions
  using Jacobian   = ActsMatrixD<5, 5>;
  using Covariance = ActsSymMatrixD<5>;

  /// the number of component constrained in propagate
  int constraint = 12;

  template <typename propagator_state_t, typename stepper_t>
  void
  operator()(propagator_state_t& state,
             const stepper_t&    stepper,
             result_type&        result) const
  {

    // If we are on target, everything should have been done
    if (state.navigation.targetReached) {
      return;
	}

    // A current surface has been already assigned by the navigator
    if (state.navigation.currentSurface
        && state.navigation.currentSurface->surfaceMaterial()) {
		  
	  // Get the multi bound parameter
    auto  bs = stepper.boundState(state.stepping, *state.navigation.currentSurface, true);
	const auto& trackMap = std::get<MultipleBoundParameters>(bs).getTrackList();
//	result.distance.push_back ( klDist(trackMap.begin()->first,(++trackMap.begin())->first) );

	using TrackParMap = decltype(trackMap);
	TrackMap mergedMap;
	TrackMap unmergedMap;
//	unmergedMap = trackMap.clone();

	}
  } 

  template<typename TrackParMap>
  pairWithMinimumDistance(TrackParMap& unmergedMap)
  {
	std::cout<<"in pairWithMinimumDistance "<<std::endl;
	std::cout<<"size of unmergedComponentsMap "<<unmergedComponentsMap.size()<<std::endl;
	double minimumDistance = 10e10;
	bool minimumDistanceSet = false;
	typename MultiComponentStateMap::iterator minimumDistanceMarker(0);
	typename MultiComponentStateMap::iterator component = unmergedComponentsMap.begin();

	for ( ; component != unmergedComponentsMap.end(); ++component){

	  // First component of the unmerged components map is the reference and distances of all other components are determined with respect to this component
	  if ( component != unmergedComponentsMap.begin() ){
		double distance = component_distance(unmergedComponentsMap.begin()->second, component->second);
		std::cout<<"kull dis "<< distance<<std::endl;
		if (distance < minimumDistance){
		  // msg(MSG::VERBOSE) << "New minimum distance set" << endmsg;
		  minimumDistance = distance;
		  minimumDistanceMarker = component;
		  minimumDistanceSet = true;
		}            
	  }
	} 

	if ( !minimumDistanceSet ){
	  std::cout<<"minimumDistance not set!!! "<<std::endl;
	  return 0;
	}   
	const component_par_t minimumDistanceComponent = minimumDistanceMarker->second;
	return std::make_shared< std::pair<const component_par_t, typename MultiComponentStateMap::iterator> > ( std::make_pair(minimumDistanceComponent, minimumDistanceMarker) );
  } 

};

}  //end of namespace detail

}  // end of namespace Acts
