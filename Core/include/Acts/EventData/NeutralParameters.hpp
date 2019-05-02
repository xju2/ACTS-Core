// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#include "Acts/EventData/ChargePolicy.hpp"
#include "Acts/EventData/SingleBoundTrackParameters.hpp"
#include "Acts/EventData/SingleCurvilinearTrackParameters.hpp"
#include "Acts/EventData/SingleTrackParameters.hpp"

#include "Acts/EventData/MultiBoundTrackParameters.hpp"
#include "Acts/EventData/MultiCurvilinearTrackParameters.hpp"
#include "Acts/EventData/MultiTrackParameters.hpp"

namespace Acts {
using NeutralParameters = SingleTrackParameters<NeutralPolicy>;
using NeutralCurvilinearParameters
    = SingleCurvilinearTrackParameters<NeutralPolicy>;
using NeutralBoundParameters = SingleBoundTrackParameters<NeutralPolicy>;

using MultipleNeutralParameters = MultiTrackParameters<NeutralPolicy>;
using MultipleNeutralCurvilinearParameters
      = MultiCurvilinearTrackParameters<NeutralPolicy>;
using MultipleNeutralBoundParameters
      = MultiBoundTrackParameters<NeutralPolicy>;
}  // namespace Acts
