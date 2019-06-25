// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/TrackParameters.hpp"

namespace Acts {
  namespace detail {
	struct ComponentCombiner{
	  template<typename weight_parameter_t>
		std::pair<double, TrackParametersBase* >	
		operator()(const GeometryContext& gctx,
			const Surface& surface,
			const weight_parameter_t& weight_parameters_1, 
			const weight_parameter_t& weight_parameters_2) const
		{
		  double weight_1 = weight_parameters_1.first;
		  double weight_2 = weight_parameters_2.first;

		  /// the simplist combination
		  double weightCombined  = weight_1 + weight_2;
		  auto parametersCombined = (weight_parameters_1.second->parameters() * weight_1 + weight_parameters_2.second->parameters() * weight_2)/weightCombined;
		  auto covCombined = (*weight_parameters_1.second->covariance() * weight_1 + *weight_parameters_2.second->covariance() * weight_2)/weightCombined;
		  std::unique_ptr<const TrackParametersBase::CovMatrix_t> covPtr = nullptr;
		  covPtr = std::make_unique<const TrackParametersBase::CovMatrix_t>(covCombined);
		  TrackParametersBase* ptr = new BoundParameters(gctx, std::move(covPtr), parametersCombined, surface.getSharedPtr()) ;
		  std::cout<<"Combine "<<weight_1<<"+"<<weight_2<<" = " <<weightCombined<<"  pararameters:  "<<weight_parameters_1.second->parameters()[4]<<","<<weight_parameters_2.second->parameters()[4]<<std::endl;
		  return  std::make_pair(weightCombined, ptr);
		} 
	};
  }  //detail
}  //Acts
