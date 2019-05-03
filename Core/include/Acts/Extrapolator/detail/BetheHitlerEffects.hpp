// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <cmath>
#include <fstream>
#include <iostream>
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Units.hpp"

namespace Acts {

namespace detail {
}
}

/*
   struct BetheHitler{

   struct Polynomial{

// Default constructor
Polynomial () {};

Polynomial (const std::vector<double>& coefficients)
:
m_coefficients(coefficients)
{};

// Evaluation of the polynomial for given material thickness (t)
double operator () (const double& t) const
{
double sum(0.);
std::vector<double>::const_iterator coefficient = m_coefficients.begin();

for ( ; coefficient != m_coefficients.end(); ++coefficient)
sum = t * sum + (*coefficient);

return sum;
}

private:
std::vector<double> m_coefficients;
};

struct ComponentValues {
ComponentValues(const double w,
const double m,
const double c)
: weight(w)
,mean(m)
,variance(c){}
double weight;
double mean;
double variance;
};


std::vector<ComponentValues>
getMixture(double pathLX0,
double momentum,
NavigationDirection navDir)  const
{
std::vector<ComponentValues> mixture;
// no Bethe-Hitler effect applied
if( pathLX0<0.0001 ){
std::cout<<"none "<<std::endl;
mixture.push_back( ComponentValues(1,0,0) );
}//end if pathLX0<0.0001

else if ( pathLX0 <0.002 ){
std::cout<<"in Single Gauss "<<std::endl;
double meanZ = exp( -1. * pathLX0);
double deltaP(0.);
double varQoverP(0.);
double sign = ( navDir == backward ) ? 1. : -1.;
if ( navDir == forward )
{
deltaP = sign * momentum * (1. - meanZ);
}
else
{
deltaP = sign * momentum * (1. / meanZ - 1.);
}
double varZ  = exp( -1. * pathLX0 * log(3.) / log(2.) ) - exp(-2. * pathLX0);
if ( navDir == forward )
{
varQoverP = 1. / (meanZ * meanZ * momentum * momentum) * varZ;
}
else
{
  varQoverP = varZ / (momentum * momentum);
}
mixture.push_back( ComponentValues(1,deltaP,varQoverP) );
}  //end if pathLX0 < 0.002

else
{
  std::cout<<"in Bethe - Hitler "<<std::endl;
  if( pathLX0 > 0.2 ) { pathLX0 = 0.2; }
  //LT
  // this should be initialize in the begining of alg
  std::string readFile =
"/afs/cern.ch/work/j/jinz/acts-core/Core/include/Acts/Extrapolator/detail/GeantSim_LT01_cdf_nC6_O5.par";
  const char* file = readFile.c_str();
  std::ifstream fin( file );
  if( !fin ) std::cout<<" read error ! "<<std::endl;
  double N;
  double order;
  double trans;
  std::vector<Polynomial> polynomial_weight;
  std::vector<Polynomial> polynomial_mean;
  std::vector<Polynomial> polynomial_val;
  fin >> N;
  fin >> order;
  fin >> trans;
  std::cout<<"Bethe - Hitler N "<<N<<std::endl;
  for( unsigned int componentIndex = 0; componentIndex < N; ++componentIndex ){
  polynomial_weight.push_back( readPolynomial(fin,order) );
  polynomial_mean.push_back( readPolynomial(fin,order) );
  polynomial_val.push_back( readPolynomial(fin,order) );
  }
  //store par into mixture
  for( unsigned int componentIndex = 0; componentIndex < N; ++componentIndex ){
  double updatedWeight = polynomial_weight[componentIndex](pathLX0);
  double updatedMean   = polynomial_mean[componentIndex](pathLX0);
  double updatedCov    = polynomial_val[componentIndex](pathLX0);
  if( trans ) {
    updatedWeight = 1./(1.+exp( -updatedWeight));
    updatedMean   = 1./(1.+exp( -updatedMean));
    updatedCov    = exp(updatedCov);
  }
  else {
    updatedCov = updatedCov * updatedCov;
  }

  double deltaP(0.);
  double deltaVar(0.);

  if ( navDir == forward )
  {
    deltaP = momentum * ( updatedMean -1 );
    deltaVar = 1./((momentum*updatedMean)*(momentum*updatedMean)*updatedCov);
  }
  else {
    deltaP = momentum * ( 1./updatedMean -1 );
    deltaVar = updatedCov/(momentum*momentum);
  }

  mixture.push_back( ComponentValues(updatedWeight,deltaP,deltaVar) );
  }

}//end if pathLX0 > 0.002
return mixture;
}

Polynomial readPolynomial (std::ifstream& fin, const int order)  const
{
  std::vector<double> coefficients(order + 1);
  int orderIndex = 0;
  for ( ; orderIndex < (order + 1); ++orderIndex ) {
  fin >> coefficients[orderIndex];
  }
  return Polynomial(coefficients);
}
};
std::ostream& operator<< (std::ostream& os, const BetheHitler::ComponentValues&
obj) {
  os<<"weight mean cov "<<obj.weight<<","<<obj.mean<<","<<obj.variance<<" ";
  return os;
}
}
}
*/
