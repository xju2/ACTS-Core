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

	  bool do1D = false;

	  template< typename weight_parameter_t>
		double	
		operator()(const weight_parameter_t& firstComponent, const weight_parameter_t& secondComponent) const
		{

		  const auto& firstComponentCov  = firstComponent.second->covariance();
		  const auto& secondComponentCov = secondComponent.second->covariance();

		  const auto& firstComponentParameters  = firstComponent.second->parameters();
		  const auto& secondComponentParameters = secondComponent.second->parameters();
		  std::cout<<"cov1 cov2 "<<*firstComponentCov<<" "<<std::endl<<*secondComponentCov<<std::endl;

		  //std::cout<<"1st "<<firstComponentParameters<<std::endl;
		  //std::cout<<"2nd "<<secondComponentParameters<<std::endl;

		  if ( do1D ){
			//only 1 dimension
			double firstPars  = firstComponentParameters[eQOP];
			double secondPars = secondComponentParameters[eQOP];

			double firstCovTrk  = (*firstComponentCov)(eQOP,eQOP);
			double secondCovTrk = (*secondComponentCov)(eQOP,eQOP);

			double G1 =  firstCovTrk > 0  ? 1./firstCovTrk : 1e10;
			double G2 =  secondCovTrk > 0 ? 1./secondCovTrk : 1e10;

			double parametersDifference = firstPars - secondPars;
			double covarianceDifference = firstCovTrk - secondCovTrk;
			double G_difference = G2 - G1;
			double G_sum        = G1 + G2;      

			double distance = covarianceDifference * G_difference + parametersDifference * G_sum * parametersDifference;

			std::cout<<"in 1D kl distance "<<std::endl;
			std::cout<<"par1 par2 "<<firstPars<<" "<<std::endl<<secondPars<<std::endl;
			std::cout<<"cov1 cov2 "<<*firstCovTrk<<" "<<std::endl<<*secondCovTrk<<std::endl;
			std::cout<<"distance "<<distance<<std::endl;
			return distance;    
		  }
		  else{
			const auto& G1 =  firstComponentCov->inverse();
			const auto& G2 =  secondComponentCov->inverse();
			const auto& parametersDifference = firstComponentParameters - secondComponentParameters;
			const auto& covarianceDifference = *firstComponentCov - *secondComponentCov;
			const auto& G_difference = G2 - G1;
			const auto& G_sum        = G1 + G2;
			double distance = (covarianceDifference * G_difference).trace() + (parametersDifference.transpose() * G_sum * parametersDifference);
			std::cout<<"in ND kl distance "<<std::endl;
			std::cout<<"par1 par2 "<<firstComponentParameters<<" "<<std::endl<<secondComponentParameters<<std::endl;
			std::cout<<"cov1 cov2 "<<firstComponentCov<<" "<<std::endl<<secondComponentCov<<std::endl;
			std::cout<<"distance "<<distance<<std::endl;
			return distance;
		  }
		} //end of operator()

	};

  }  //detail
}  //Acts
