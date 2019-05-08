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
  struct TunningEnergyEffect 
  {
    const static int N = 2;
    struct ComponentValues
    {
      double weight; 
      double mean;
      double variance; 
    };

    std::vector<ComponentValues>
    operator()(double /*noused*/, double /*noused*/) const
    {
      // make non effect, just means copy the components in material effect
	  ComponentValues comp1{0.7,-0.01,0.01};
	  ComponentValues comp2{0.3,-0.02,0.01};
	  std::vector<ComponentValues> compVec;
	  compVec.push_back(comp1);
	  compVec.push_back(comp2);
	  return std::move(compVec);

    }
  };
}
}
