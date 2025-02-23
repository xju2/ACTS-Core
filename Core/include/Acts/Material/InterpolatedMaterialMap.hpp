// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <functional>
#include <optional>
#include "Acts/Material/Material.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Interpolation.hpp"

namespace Acts {

/// @brief Struct for mapping global 3D positions to material values
///
/// Global 3D positions are transformed into a @c DIM_POS Dimensional vector
/// which is used to look up the material classification value in the
/// underlying material map.
template <typename G>
struct MaterialMapper {
 public:
  using Grid_t = G;
  static constexpr size_t DIM_POS = Grid_t::DIM;

  /// @brief Struct representing smallest grid unit in material grid
  ///
  /// This type encapsulate all required information to perform linear
  /// interpolation of material classification values within a 3D volume.
  struct MaterialCell {
    /// Number of corner points defining the confining hyper-box
    static constexpr unsigned int N = 1 << DIM_POS;

    /// @brief Default constructor
    ///
    /// @param [in] transformPos   Mapping of global 3D coordinates onto grid
    /// space
    /// @param [in] lowerLeft   Generalized lower-left corner of hyper box
    ///                         (containing the minima of the hyper box along
    ///                         each Dimension)
    /// @param [in] upperRight  Generalized upper-right corner of hyper box
    ///                         (containing the maxima of the hyper box along
    ///                         each Dimension)
    /// @param [in] materialValues Material classification values at the hyper
    /// box corners sorted in the canonical order defined in Acts::interpolate
    MaterialCell(
        std::function<ActsVectorD<DIM_POS>(const Vector3D&)> transformPos,
        std::array<double, DIM_POS> lowerLeft,
        std::array<double, DIM_POS> upperRight,
        std::array<ActsVectorF<5>, N> materialValues)
        : m_transformPos(std::move(transformPos)),
          m_lowerLeft(std::move(lowerLeft)),
          m_upperRight(std::move(upperRight)),
          m_materialValues(std::move(materialValues)) {}

    /// @brief Retrieve material at given position
    ///
    /// @param [in] position Global 3D position
    /// @return Material at the given position
    ///
    /// @pre The given @c position must lie within the current cell.
    Material getMaterial(const Vector3D& position) const {
      // defined in Interpolation.hpp
      return Material(interpolate(m_transformPos(position), m_lowerLeft,
                                  m_upperRight, m_materialValues));
    }

    /// @brief Check whether given 3D position is inside this cell
    ///
    /// @param [in] position Global 3D position
    /// @return @c true if position is inside the current cell,
    ///         otherwise @c false
    bool isInside(const Vector3D& position) const {
      const auto& gridCoordinates = m_transformPos(position);
      for (unsigned int i = 0; i < DIM_POS; ++i) {
        if (gridCoordinates[i] < m_lowerLeft.at(i) ||
            gridCoordinates[i] >= m_upperRight.at(i)) {
          return false;
        }
      }
      return true;
    }

   private:
    /// Geometric transformation applied to global 3D positions
    std::function<ActsVectorD<DIM_POS>(const Vector3D&)> m_transformPos;

    /// Generalized lower-left corner of the confining hyper-box
    std::array<double, DIM_POS> m_lowerLeft;

    /// Generalized upper-right corner of the confining hyper-box
    std::array<double, DIM_POS> m_upperRight;

    /// @brief Material component vectors at the hyper-box corners
    ///
    /// @note These values must be order according to the prescription detailed
    ///       in Acts::interpolate.
    std::array<ActsVectorF<5>, N> m_materialValues;
  };

  /// @brief Default constructor
  ///
  /// @param [in] transformPos Mapping of global 3D coordinates (cartesian)
  /// onto grid space
  /// @param [in] grid Grid storing material classification values
  MaterialMapper(
      std::function<ActsVectorD<DIM_POS>(const Vector3D&)> transformPos,
      Grid_t grid)
      : m_transformPos(std::move(transformPos)), m_grid(std::move(grid)) {}

  /// @brief Retrieve material at given position
  ///
  /// @param [in] position Global 3D position
  /// @return Material at the given position
  ///
  /// @pre The given @c position must lie within the range of the underlying
  /// map.
  Material getMaterial(const Vector3D& position) const {
    return Material(m_grid.interpolate(m_transformPos(position)));
  }

  /// @brief Retrieve material cell for given position
  ///
  /// @param [in] position Global 3D position
  /// @return material cell containing the given global position
  ///
  /// @pre The given @c position must lie within the range of the underlying
  /// map.
  MaterialCell getMaterialCell(const Vector3D& position) const {
    const auto& gridPosition = m_transformPos(position);
    size_t bin = m_grid.globalBinFromPosition(gridPosition);
    const auto& indices = m_grid.localBinsFromPosition(bin);
    const auto& lowerLeft = m_grid.lowerLeftBinEdge(indices);
    const auto& upperRight = m_grid.upperRightBinEdge(indices);

    // Loop through all corner points
    constexpr size_t nCorners = 1 << DIM_POS;
    std::array<ActsVectorF<5>, nCorners> neighbors;
    const auto& cornerIndices = m_grid.closestPointsIndices(gridPosition);

    size_t i = 0;
    for (size_t index : cornerIndices) {
      neighbors.at(i++) = m_grid.at(index);
    }

    return MaterialCell(m_transformPos, lowerLeft, upperRight,
                        std::move(neighbors));
  }

  /// @brief Get the number of bins for all axes of the map
  ///
  /// @return Vector returning number of bins for all map axes
  std::vector<size_t> getNBins() const {
    auto nBinsArray = m_grid.numLocalBins();
    return std::vector<size_t>(nBinsArray.begin(), nBinsArray.end());
  }

  /// @brief Get the minimum value of all axes of the map
  ///
  /// @return Vector returning the minima of all map axes
  std::vector<double> getMin() const {
    auto minArray = m_grid.minPosition();
    return std::vector<double>(minArray.begin(), minArray.end());
  }

  /// @brief Get the maximum value of all axes of the map
  ///
  /// @return Vector returning the maxima of all map axes
  std::vector<double> getMax() const {
    auto maxArray = m_grid.maxPosition();
    return std::vector<double>(maxArray.begin(), maxArray.end());
  }

  /// @brief Check whether given 3D position is inside look-up domain
  ///
  /// @param [in] position Global 3D position
  /// @return @c true if position is inside the defined look-up grid,
  ///         otherwise @c false
  bool isInside(const Vector3D& position) const {
    return m_grid.isInside(m_transformPos(position));
  }

  /// @brief Get a const reference on the underlying grid structure
  ///
  /// @return Grid reference
  const Grid_t& getGrid() const { return m_grid; }

 private:
  /// Geometric transformation applied to global 3D positions
  std::function<ActsVectorD<DIM_POS>(const Vector3D&)> m_transformPos;
  /// Grid storing material values
  Grid_t m_grid;
};

/// @ingroup Material
/// @brief Interpolate material classification values from material values on a
/// given grid
///
/// This class implements a material service which is initialized by a
/// material map defined by:
/// - a list of material values on a regular grid in some n-Dimensional space,
/// - a transformation of global 3D coordinates onto this n-Dimensional
/// space.
/// - a transformation of local n-Dimensional material coordinates into
/// global (cartesian) 3D coordinates
///
/// The material value for a given global position is then determined by:
/// - mapping the position onto the grid,
/// - looking up the material classification values on the closest grid points,
/// - doing a linear interpolation of these values.
/// @warning Each classification number of the material is interpolated
/// independently and thus does not consider any correlations that exists
/// between these values. This might work out since the used material is already
/// a mean of the materials in a certain bin and can therewith be treated as a
/// collection of numbers.
/// @tparam G Type of the grid
template <typename Mapper_t>
class InterpolatedMaterialMap final {
 public:
  /// @brief Temporary storage of a certain cell to improve material access
  struct Cache {
    /// Stored material cell
    std::optional<typename Mapper_t::MaterialCell> matCell;
    /// Boolean statement if the cell is initialized
    bool initialized = false;
  };

  /// @brief Create interpolated map
  ///
  /// @param [in] mapper Material map
  InterpolatedMaterialMap(Mapper_t& mapper) : m_mapper(std::move(mapper)) {}

  /// @brief Retrieve material
  ///
  /// @param [in] position Global 3D position
  ///
  /// @return Material at given position
  Material getMaterial(const Vector3D& position) const {
    return m_mapper.getMaterial(position);
  }

  /// @brief Retrieve material
  ///
  /// @param [in] position Global 3D position
  /// @param [in,out] cache Cache object. Contains material cell used for
  /// interpolation
  ///
  /// @return material at given position
  Material getMaterial(const Vector3D& position, Cache& cache) const {
    if (!cache.initialized || !(*cache.matCell).isInside(position)) {
      cache.matCell = getMaterialCell(position);
      cache.initialized = true;
    }
    return (*cache.matCell).getMaterial(position);
  }

  /// @brief Retrieve material value & its "gradient"
  ///
  /// @param [in]  position   Global 3D position
  /// @param [out] derivative "Gradient" of material as (5x5) matrix
  /// @return Material
  ///
  /// @note Currently the derivative is not calculated
  /// @todo return derivative
  Material getMaterialGradient(const Vector3D& position,
                               ActsMatrixD<5, 5>& /*derivative*/) const {
    return m_mapper.getMaterial(position);
  }

  /// @brief Retrieve material value & its "gradient"
  ///
  /// @param [in]  position   Global 3D position
  /// @param [out] derivative "Gradient" of material as (5x5) matrix
  /// @param [in,out] cache Cache object. Contains cell used for
  /// @return Material
  ///
  /// @note Currently the derivative is not calculated
  /// @note Cache is not used currently
  /// @todo return derivative
  Material getMaterialGradient(const Vector3D& position,
                               ActsMatrixD<5, 5>& /*derivative*/,
                               Cache& /*cache*/) const {
    return m_mapper.getMaterial(position);
  }

  /// @brief Convenience method to access underlying material mapper
  ///
  /// @return The material mapper
  Mapper_t getMapper() const { return m_mapper; }

  /// @brief Check whether given 3D position is inside look-up domain
  ///
  /// @param [in] position Global 3D position
  /// @return @c true if position is inside the defined map, otherwise @c false
  bool isInside(const Vector3D& position) const {
    return m_mapper.isInside(position);
  }

 private:
  /// @brief Retrieve cell for given position
  ///
  /// @param [in] position Global 3D position
  /// @return Material cell containing the given global position
  ///
  /// @pre The given @c position must lie within the range of the underlying
  /// map.
  typename Mapper_t::MaterialCell getMaterialCell(
      const Vector3D& position) const {
    return m_mapper.getMaterialCell(position);
  }

  /// @brief object for global coordinate transformation and interpolation
  ///
  /// This object performs the mapping of the global 3D coordinates onto the
  /// material grid and the interpolation of the material component values on
  /// close-by grid points.
  Mapper_t m_mapper;
};
}  // namespace Acts
