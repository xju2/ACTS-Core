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
    const static int N = 6;
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
      ComponentValues              comp1{0.3, -0.0001, 0.001};
      ComponentValues              comp2{0.25, -0.0004, 0.001};
      ComponentValues              comp3{0.25, -0.0008, 0.001};
      ComponentValues              comp4{0.1, -0.002, 0.001};
      ComponentValues              comp5{0.05, -0.003, 0.001};
      ComponentValues              comp6{0.05, -0.004, 0.001};
      std::vector<ComponentValues> compVec;
      compVec.push_back(comp1);
      compVec.push_back(comp2);
      compVec.push_back(comp3);
      compVec.push_back(comp4);
      compVec.push_back(comp5);
      compVec.push_back(comp6);
      return std::move(compVec);
    }
  };
}
}
