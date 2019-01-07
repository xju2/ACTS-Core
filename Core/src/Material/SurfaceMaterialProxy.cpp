// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// ProtoSurfaceMaterial.cpp, Acts project
///////////////////////////////////////////////////////////////////

#include "Acts/Material/ProtoSurfaceMaterial.hpp"

Acts::ProtoSurfaceMaterial::ProtoSurfaceMaterial(const BinUtility& binUtility)
  : SurfaceMaterial(), m_binUtility(binUtility)
{
}

Acts::ProtoSurfaceMaterial&
Acts::ProtoSurfaceMaterial::operator*=(double /*scale*/)
{
  return (*this);
}

std::ostream&
Acts::SurfaceMaterialProxy::toStream(std::ostream& sl) const
{
  sl << "Acts::ProtoSurfaceMaterial : " << std::endl;
  if (m_binUtility.bins(0) * m_binUtility.bins(1) > 1) {
    sl << "   - Number of Material bins [0,1] : " << m_binUtility.bins(0)
       << " / " << m_binUtility.bins(1) << std::endl;
  } else {
    sl << "   - Homogeneous Material" << std::endl;
  }
  return sl;
}
