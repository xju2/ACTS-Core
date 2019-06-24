// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

namespace Acts {

  namespace detail {

	/// @brief Use the Kullback-Leibler distance 
	/// For two arbitrary probability density function, p1, p2 (Gaussian pdf), the KL distance should be 
	/// Dkl = tr[ (V1-V2)(W1-W2)]+(u1-u2).Transpose()(W1+W2)(u1-u2)
	/// e.g. track model, u1 represents meaning of parameter of track1, V1 represents meaning of covaiance of track1.
	///
	/// Currently, the KL distance is calculated in 1 dimension(eQOP) 

	/// calculate distances between two components
	/// weight_parameter_t : std::pair<weight, unique_ptr<boundpar> >
	struct KullbackLeiblerComponentDistance{

	  bool do1D = true;

	  template< typename weight_parameter_t>
		double	
		operator()(const weight_parameter_t& weight_parameters_1, const weight_parameter_t& weight_parameters_2) const
		{
		  const auto& parameters_1 = weight_parameters_1.second->parameters();
		  const auto& parameters_2 = weight_parameters_2.second->parameters();

		  const auto& cov_1 = weight_parameters_1.second->covariance();
		  const auto& cov_2 = weight_parameters_2.second->covariance();
		 
		  if ( do1D ){
			//only 1 dimension
			double firstCovTrk  = (*cov_1)(eQOP,eQOP);
			double secondCovTrk = (*cov_2)(eQOP,eQOP);
			double G1 =  firstCovTrk > 0  ? 1./firstCovTrk : 1e10;
			double G2 =  secondCovTrk > 0 ? 1./secondCovTrk : 1e10;
			double parametersDifference = parameters_1[eQOP] - parameters_2[eQOP];
			double covarianceDifference = firstCovTrk - secondCovTrk;
			double G_difference = G2 - G1;
			double G_sum        = G1 + G2;      

			double distance = covarianceDifference * G_difference + parametersDifference * G_sum * parametersDifference;
			return distance;    
		  }
		  else{
			const auto& G1 =  cov_1->inverse();
			const auto& G2 =  cov_2->inverse();
			const auto& parametersDifference = parameters_1 - parameters_2;
			const auto& covarianceDifference = *cov_1 - *cov_2;
			const auto& G_difference = G2 - G1;
			const auto& G_sum        = G1 + G2;
			double distance = (covarianceDifference * G_difference).trace() + (parametersDifference.transpose() * G_sum * parametersDifference);
			return distance;    
		  }
		} 

	};

  }  //detail
}  //Acts
