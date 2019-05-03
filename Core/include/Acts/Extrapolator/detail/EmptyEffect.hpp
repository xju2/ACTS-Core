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
  struct EmptyEffect
  {
    const static int N = 2;
    struct ComponentValues
    {
      double weight   = 1. / N;
      double mean     = 0.;
      double variance = 0.;
    };
    std::vector<ComponentValues>
    getMixture(double /*unused*/, double /*unused*/) const
    {
      // make non effect, just means copy the components in material effect
      return std::move(std::vector<ComponentValues>(N));
    }
  };
}
}

