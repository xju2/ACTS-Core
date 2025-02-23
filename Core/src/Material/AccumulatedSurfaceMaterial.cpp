// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// AccumulatedSurfaceMaterial.cpp, Acts project
///////////////////////////////////////////////////////////////////

#include "Acts/Material/AccumulatedSurfaceMaterial.hpp"
#include "Acts/Material/BinnedSurfaceMaterial.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"

// Default Constructor - for homogeneous material
Acts::AccumulatedSurfaceMaterial::AccumulatedSurfaceMaterial(double splitFactor)
    : m_splitFactor(splitFactor) {
  AccumulatedVector accMat = {{AccumulatedMaterialProperties()}};
  m_accumulatedMaterial = {{accMat}};
}

// Binned Material constructor with split factor
Acts::AccumulatedSurfaceMaterial::AccumulatedSurfaceMaterial(
    const BinUtility& binUtility, double splitFactor)
    : m_binUtility(binUtility), m_splitFactor(splitFactor) {
  size_t bins0 = m_binUtility.bins(0);
  size_t bins1 = m_binUtility.bins(1);
  AccumulatedVector accVec(bins0, AccumulatedMaterialProperties());
  m_accumulatedMaterial = AccumulatedMatrix(bins1, accVec);
}

// Assign a material properites object
std::array<size_t, 3> Acts::AccumulatedSurfaceMaterial::accumulate(
    const Vector2D& lp, const MaterialProperties& mp, double pathCorrection) {
  if (m_binUtility.dimensions() == 0) {
    m_accumulatedMaterial[0][0].accumulate(mp, pathCorrection);
    return {0, 0, 0};
  }
  size_t bin0 = m_binUtility.bin(lp, 0);
  size_t bin1 = m_binUtility.bin(lp, 1);
  m_accumulatedMaterial[bin1][bin0].accumulate(mp, pathCorrection);
  return {bin0, bin1, 0};
}

// Assign a material properites object
std::array<size_t, 3> Acts::AccumulatedSurfaceMaterial::accumulate(
    const Vector3D& gp, const MaterialProperties& mp, double pathCorrection) {
  if (m_binUtility.dimensions() == 0) {
    m_accumulatedMaterial[0][0].accumulate(mp, pathCorrection);
    return {0, 0, 0};
  }
  std::array<size_t, 3> bTriple = m_binUtility.binTriple(gp);
  m_accumulatedMaterial[bTriple[1]][bTriple[0]].accumulate(mp, pathCorrection);
  return bTriple;
}

// Average the information accumulated during one event
void Acts::AccumulatedSurfaceMaterial::trackAverage(
    const std::vector<std::array<size_t, 3>>& trackBins) {
  // The touched bins are known, so you can access them directly
  if (not trackBins.empty()) {
    for (auto bin : trackBins) {
      m_accumulatedMaterial[bin[1]][bin[0]].trackAverage();
    }
  } else {
    // Run over all bins
    for (auto& matVec : m_accumulatedMaterial) {
      for (auto& mat : matVec) {
        mat.trackAverage();
      }
    }
  }
}

/// Total average creates SurfaceMaterial
std::unique_ptr<const Acts::ISurfaceMaterial>
Acts::AccumulatedSurfaceMaterial::totalAverage() {
  if (m_binUtility.bins() == 1) {
    // Return HomogeneousSurfaceMaterial
    return std::make_unique<HomogeneousSurfaceMaterial>(
        m_accumulatedMaterial[0][0].totalAverage().first, m_splitFactor);
  }
  // Create the properties matrix
  MaterialPropertiesMatrix mpMatrix(
      m_binUtility.bins(1),
      MaterialPropertiesVector(m_binUtility.bins(0), MaterialProperties()));
  // Loop over and fill
  for (size_t ib1 = 0; ib1 < m_binUtility.bins(1); ++ib1) {
    for (size_t ib0 = 0; ib0 < m_binUtility.bins(0); ++ib0) {
      mpMatrix[ib1][ib0] = m_accumulatedMaterial[ib1][ib0].totalAverage().first;
    }
  }
  // Now return the BinnedSurfaceMaterial
  return std::make_unique<const BinnedSurfaceMaterial>(
      m_binUtility, std::move(mpMatrix), m_splitFactor);
}
