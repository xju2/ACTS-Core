// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <map>
#include <utility>
#include <list>
#include <memory>
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Utilities/GeometryContext.hpp"

namespace Acts {

  namespace detail {

	struct ComponentCombiner{

	  template< typename weight_parameter_t>
		std::pair<double, TrackParametersBase* >	
		operator()(const GeometryContext& gctx,
			const Surface& surface,
			const weight_parameter_t& uncombined1, 
			const weight_parameter_t& uncombined2) const
		{
		  using Covariance       = ActsSymMatrixD<5>;
		  double weight1 = uncombined1.first;
		  double weight2 = uncombined2.first;

		  double weightCombined  = weight1 + weight2;
		  auto parameterCombined = uncombined1.second->parameters() * weight1 + uncombined2.second->parameters() * weight2;

		  // a temporary version
		  auto covPart1 = *uncombined1.second->covariance() * weight1 + *uncombined2.second->covariance() * weight2;
		  std::unique_ptr<const Covariance> covPtr = nullptr;
		  covPtr = std::make_unique<const Covariance>(covPart1);

		  //std::unique_ptr<BoundParameters> combinedPar( new BoundParameters(gctx, covPtr, parameterCombined, surface.getSharedPtr()) );
		  TrackParametersBase* ptr = new BoundParameters(gctx, std::move(covPtr), parameterCombined, surface.getSharedPtr()) ;

//		  std::pair<double,std::unique_ptr<TrackParametersBase> > parPair( weightCombined, std::move(uptr) );
//		  return std::move(parPair);
		  return  std::make_pair( weightCombined, ptr ) ;

//		  auto trackPair = std::make_pair<double,std::unique_ptr<TrackParametersBase> >( weightCombined, std::move(uptr) );
		  //return   std::move(trackPair);
		} //end of operator()
	};

  }  //detail
}  //Acts
