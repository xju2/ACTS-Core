// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#include "Acts/Extrapolator/detail/ComponentDistance.hpp"
#include "Acts/Extrapolator/detail/ComponentCombiner.hpp"

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
  detail::ComponentCombiner combiner;

  using result_type = this_result;
  /// Jacobian, Covariance and State defintions
  using Jacobian   = ActsMatrixD<5, 5>;
  using Covariance = ActsSymMatrixD<5>;

  /// the number of component constrained in propagate
  int constraintNum = 12;

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
	auto& trackMap = std::get<MultipleBoundParameters>(bs).getTrackList();
//	result.distance.push_back ( klDist(trackMap.begin()->first,(++trackMap.begin())->first) );

	using TrackParMap = typename std::remove_reference<decltype(trackMap)>::type;
	// the aim merging map
	TrackParMap& unmergedMap = trackMap;
	// the empty map to merge
	TrackParMap mergedMap;
	size_t numberOfComponents = unmergedMap.size();
	while ( numberOfComponents > constraintNum ){
	  if( !unmergedMap.empty() ) unmergedMap.clear();
	  while ( numberOfComponents > constraintNum && !unmergedMap.empty() ) {
		if( unmergedMap.size() > 1 ){
		  auto mergedComponentIter = pairWithMinimumDistance(trackMap);
//		  if (mergedComponentIter){
			auto& combinedComponent = *combiner(trackMap.begin()->first, mergedComponentIter->first);
			unmergedMap.erase(mergedComponentIter);
			unmergedMap.erase(trackMap.begin());
			mergedMap.insert( std::make_pair(std::move(combinedComponent),unmergedMap.size()) );
//		  }
//		  else {
//			unmergedMap.erase( trackMap.begin() );
//		  }
		  --numberOfComponents;
		}
		else{
		  mergedMap.insert( std::move( *trackMap.begin() ));
		}
	  }
	  if( unmergedMap.empty() && numberOfComponents > constraintNum ) unmergedMap = mergedMap; //move
	}

	/*
	if( unmergedMap.size()+mergedMap.size() > constraintNum ) {
	  if( unmergedMap.empty() ) unmergedMap.clear();
		if( unmergedMap.size()+mergedMap.size() > constraintNum && !unmergedMap.empty() ){
		  auto mergedComponentIter = pairWithMinimumDistance(trackMap);
		  auto combinedComponent = combiner(trackMap.begin()->first, mergedComponentIter->first);
		  unmergedMap.erase(mergedComponentIter);
		  unmergedMap.erase(trackMap.begin());
		  mergedMap.insert(std::make_pair(std::move(combinedComponent),unmergedMap.size()) );
		}
		if( unmergedMap.size()+mergedMap.size() > constraintNum) unmergedMap = mergedMap;
	}
	*/

	}
  }

  template<typename TrackParMap>
	auto 
  pairWithMinimumDistance(const TrackParMap& unmergedMap) const 
  		->typename TrackParMap::const_iterator
  {
	double minmumDist = 10e10;
	bool minmumDistSet = false;
	typename TrackParMap::const_iterator it = unmergedMap.begin();
	typename TrackParMap::const_iterator it_record ;

	for( ; it != unmergedMap.end(); it++ ){
	  if ( it != unmergedMap.begin() ){
		double distance = klDist((*it).first,unmergedMap.begin()->first);
		  if( distance < minmumDist ){
			minmumDist = distance;
			it_record = it;
			minmumDistSet = true;
		  }
	  }
	  return it;

	}
  } 

};

}  //end of namespace detail

}  // end of namespace Acts
