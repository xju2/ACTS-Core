// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "Acts/Surfaces/Surface.hpp"

namespace Acts {
  namespace impl{

	template<typename T>
	  typename T::const_iterator
	  pairWithMinimumDistance(const T& unmergedMap) ;

	template<typename T>
	  void reductComponent( T& trackMap, size_t constraintNum, const GeometryContext& geoContext, const Surface& surface)
	  {
		detail::ComponentCombiner combiner;
		using TrackParMap = typename std::remove_reference<decltype(trackMap)>::type;
		// the aim merging map
		TrackParMap& unmergedMap = trackMap;
		// the empty map to merge
		TrackParMap mergedMap;
		size_t numberOfComponents = unmergedMap.size();
		while ( numberOfComponents > constraintNum ){
		  if( !mergedMap.empty() ) mergedMap.clear();
		  while ( numberOfComponents > constraintNum && !unmergedMap.empty() ) {
			if( unmergedMap.size() > 1 ){
			  auto mergedComponentIter = pairWithMinimumDistance(trackMap);
			  auto combinedComponent = combiner( geoContext, surface, *trackMap.begin(),*mergedComponentIter);
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
		  if( unmergedMap.empty() && numberOfComponents > constraintNum ) {
			unmergedMap = std::move(mergedMap); //move
		  }
		}
		// merge two maps
		auto it = mergedMap.begin();
		for( ; it != mergedMap.end(); it++ ){
		  auto& par    = (*it).second;
		  auto& weight = (*it).first;
		  unmergedMap.insert( std::make_pair(weight, std::move(par)) );
		}
	  }

	template<typename T>
	  auto 
	  pairWithMinimumDistance(const T& unmergedMap) 
	  ->typename T::const_iterator
	  {
		detail::KullbackLeiblerComponentDistance klDist;
		double minmumDist = 10e10;
		typename T::const_iterator it = unmergedMap.begin();
		typename T::const_iterator it_record ;

		for( ; it != unmergedMap.end(); it++ ){
		  if ( it == unmergedMap.begin() ) continue;
		  double distance = klDist((*it),*unmergedMap.begin());
		  if( distance < minmumDist ){
			minmumDist = distance;
			it_record = it;
		  }
		}
		return it_record;
	  }
  }
}
