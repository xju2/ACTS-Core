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
  static const int constraintNum = 6;

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
    if (state.navigation.currentSurface){
    
		  
	//stepper.outPut(state.stepping);
	
	// Get the multi bound parameter
    auto  bs = stepper.boundState(state.stepping, *state.navigation.currentSurface, true);
	auto& multipleBoundPar =  std::get<MultipleBoundParameters>(bs);
	auto& trackMap = multipleBoundPar.getTrackList();
//	auto mergedComponentIter = pairWithMinimumDistance(trackMap);
//	result.distance.push_back ( klDist(trackMap.begin()->first,(++trackMap.begin())->first) );
//	auto combinedParameter = combiner( state.stepping.geoContext, *state.navigation.currentSurface, trackMap.begin()->first,mergedComponentIter->first);
//	MultipleBoundParameters newMultiBoundPar((*state.navigation.currentSurface).getSharedPtr());

	// lazy append
//	newMultiBoundPar.append( combinedParameter.first, combinedParameter.second );

	// update the stateCol
//	stepper.update(state.stepping, newMultiBoundPar);

	//stepper.outPut(state.stepping);

	
	using TrackParMap = typename std::remove_reference<decltype(trackMap)>::type;
	using WeightedTrackPars = std::pair<double, std::unique_ptr<TrackParametersBase>>;
	// the aim merging map
	TrackParMap& unmergedMap = trackMap;
	// the empty map to merge
	TrackParMap mergedMap;
	size_t numberOfComponents = unmergedMap.size();
	while ( numberOfComponents > constraintNum ){
	  if( mergedMap.empty() ) mergedMap.clear();
	  while ( numberOfComponents > constraintNum && !unmergedMap.empty() ) {
		if( unmergedMap.size() > 1 ){
		  auto mergedComponentIter = pairWithMinimumDistance(trackMap);
		  //outPutMap(trackMap);
			auto combinedComponent = combiner( state.stepping.geoContext, *state.navigation.currentSurface, *trackMap.begin(),*mergedComponentIter);
			unmergedMap.erase(mergedComponentIter);
			unmergedMap.erase(trackMap.begin());
			mergedMap.insert(std::make_pair(combinedComponent.first,std::unique_ptr<TrackParametersBase>(combinedComponent.second)));
			--numberOfComponents;
		}
		else{
		  auto& lastComponent = *unmergedMap.begin();
		  auto& lastPar    = lastComponent.second;
		  auto& lastWeight = lastComponent.first;
		  mergedMap.insert( std::make_pair(lastWeight, std::move(lastPar)) );
		  unmergedMap.erase( unmergedMap.begin());
		}
	  }
	  if( unmergedMap.empty() && numberOfComponents > constraintNum ) unmergedMap = std::move(mergedMap); //move
	}

	// merge two maps
	auto it = mergedMap.begin();
	for( ; it != mergedMap.end(); it++ ){
	  auto& par    = (*it).second;
	  auto& weight = (*it).first;
	  unmergedMap.insert( std::make_pair(weight, std::move(par)) );
	}
	std::cout<<"in reductComponent "<<trackMap.size()<<std::endl;
	stepper.update(state.stepping, multipleBoundPar);
	}
  }

  template<typename TrackParMap>
	auto 
	pairWithMinimumDistance(const TrackParMap& unmergedMap) const 
	->typename TrackParMap::const_iterator
	{
	  double minmumDist = 10e10;
	  typename TrackParMap::const_iterator it = unmergedMap.begin();
	  typename TrackParMap::const_iterator it_record ;

	  for( ; it != unmergedMap.end(); it++ ){
		if ( it == unmergedMap.begin() ) continue;
		double distance = klDist((*it),*unmergedMap.begin());
		if( distance < minmumDist ){
		  minmumDist = distance;
		  it_record = it;
		}
		return it_record;

	  }
	} 

  template<typename TrackParMap>
	void
	outPutMap(const TrackParMap& unmergedMap) const 
	{
	  typename TrackParMap::const_iterator it = unmergedMap.begin();
	  for( ; it != unmergedMap.end(); it++ ){
		std::cout<<"parameter "<<(*it).second->parameters()<<std::endl;
	  }
	}

};

}  //end of namespace detail

}  // end of namespace Acts
