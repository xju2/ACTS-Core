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
  static const int constraintNum = 12;

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
		  
	//stepper.outPut(state.stepping);
	
	// Get the multi bound parameter
    auto  bs = stepper.boundState(state.stepping, *state.navigation.currentSurface, true);
	auto& trackMap = std::get<MultipleBoundParameters>(bs).getTrackList();
//	auto mergedComponentIter = pairWithMinimumDistance(trackMap);
//	result.distance.push_back ( klDist(trackMap.begin()->first,(++trackMap.begin())->first) );
//	auto combinedParameter = combiner( state.stepping.geoContext, *state.navigation.currentSurface, trackMap.begin()->first,mergedComponentIter->first);
//	MultipleBoundParameters newMultiBoundPar((*state.navigation.currentSurface).getSharedPtr());

	// lazy append
//	newMultiBoundPar.append( combinedParameter.first, combinedParameter.second );

	// update the stateCol
//	stepper.update(state.stepping, newMultiBoundPar);

	//stepper.outPut(state.stepping);

/*	
	std::cout<<"print parameter "<<std::endl;
	std::cout<<"parameters 1 "<<std::endl;
	std::cout<<"par weight "<<trackMap.begin()->first.first<<std::endl;
	std::cout<<"par        "<<trackMap.begin()->first.second->parameters()<<std::endl;
	std::cout<<"cov 	  "<<*trackMap.begin()->first.second->covariance()<<std::endl;
	std::cout<<"parameters 2 "<<std::endl;
	std::cout<<"par weight "<<mergedComponentIter->first.first<<std::endl;
	std::cout<<"par        "<<mergedComponentIter->first.second->parameters()<<std::endl;
	std::cout<<"cov 	  "<<*mergedComponentIter->first.second->covariance()<<std::endl;
	std::cout<<"after combine "<<std::endl;
	std::cout<<"par weight "<<combinedParameter.first<<std::endl;
	std::cout<<"par "<<combinedParameter.second->parameters()<<std::endl;
	std::cout<<"cov "<<*combinedParameter.second->covariance()<<std::endl;
	std::cout<<std::endl;
*/
	


	
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
//			auto combinedComponent = combiner(state.stepping.geoContext, trackMap.begin()->first, mergedComponentIter->first);
			auto combinedComponent = combiner( state.stepping.geoContext, *state.navigation.currentSurface, trackMap.begin()->first,mergedComponentIter->first);
			unmergedMap.erase(mergedComponentIter);
			unmergedMap.erase(trackMap.begin());
			using WeightedTrackPars = std::pair<double, std::unique_ptr<TrackParametersBase>>;
			WeightedTrackPars weight_parameter( combinedComponent.first, std::unique_ptr<TrackParametersBase>(combinedComponent.second) );
			mergedMap.insert(std::make_pair(std::move(weight_parameter),0 ));
//			mergedMap.insert( std::make_pair(std::move(combinedComponent),unmergedMap.size()) );
//		  }
//		  else {
//			unmergedMap.erase( trackMap.begin() );
//		  }
		  --numberOfComponents;
		}
		else{
//		  mergedMap.insert( std::move( *trackMap.begin() ));
		}
	  }
//	  if( unmergedMap.empty() && numberOfComponents > constraintNum ) unmergedMap = mergedMap; //move
	}
	

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
	  if ( it == unmergedMap.begin() ) continue;
	  double distance = klDist((*it).first,unmergedMap.begin()->first);
	  std::cout<<"dist "<<distance<<std::endl;
	  if( distance < minmumDist ){
		minmumDist = distance;
		it_record = it;
		minmumDistSet = true;
	  }
	  return it_record;

	}
  } 

};

}  //end of namespace detail

}  // end of namespace Acts
