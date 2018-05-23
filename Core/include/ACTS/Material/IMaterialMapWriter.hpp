// This file is part of the ACTS project.
//
// Copyright (C) 2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <map>

namespace Acts {

class SurfaceMaterial;
class GeometryID;

class IMaterialMapWriter
{

public:
  virtual void
  write(const std::map<GeometryID, SurfaceMaterial*>& surfaceMaterialMap) const
      = 0;

  //virtual void
  //write(const GeometryID& geo_id, const SurfaceMaterial* surfaceMaterial)
      //= 0;
};
}
