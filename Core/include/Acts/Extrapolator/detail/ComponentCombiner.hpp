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

	struct ComponentCombiner{

	  template< typename weight_parameter_t>
		weight_parameter_t&	
		operator()(const weight_parameter_t& uncombined1, const weight_parameter_t& uncombined2) const
		{
			return uncombined1;
		} //end of operator()

	};

  }  //detail
}  //Acts
