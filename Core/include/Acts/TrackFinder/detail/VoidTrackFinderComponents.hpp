// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#include <limits>

#include "Acts/EventData/SourceLinkConcept.hpp"
#include "Acts/Geometry/GeometryID.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Utilities/TypeTraits.hpp"

namespace Acts {

/// @brief This struct selects the source link most compatible with the
/// track parameter according to the chi2
///
/// The chi2 with the selected source link should satisfy provided criteria
///
struct VoidSourceLinkSelector {
  /// @brief nested config struct
  ///
  struct Config {
    // Allowed maximum chi2
    double maxChi2 = std::numeric_limits<double>::max();
  };

  /// @brief Default constructor
  VoidSourceLinkSelector() = default;

  /// @brief Constructor with config and (non-owning) logger
  ///
  /// @param config a config instance
  /// @param logger a logger instance
  VoidSourceLinkSelector(
      Config cfg,
      std::shared_ptr<const Logger> logger = std::shared_ptr<const Logger>(
          getDefaultLogger("VoidSourceLinkSelector", Logging::WARNING)
              .release()))
      : m_config(std::move(cfg)), m_logger(std::move(logger)) {}

  /// @brief Operater that select the source links compatible with
  /// the given track parameter on a surface
  ///
  /// @tparam calibrator_t The type of calibrator
  /// @tparam source_link_t The type of source link
  ///
  /// @param calibrator The measurement calibrator
  /// @param predictedParams The predicted track parameter on a surface
  /// @param sourcelinks The pool of source links
  ///
  /// @return the compatible source links
  template <typename calibrator_t, typename source_link_t>
  std::vector<source_link_t> operator()(
      const calibrator_t& calibrator, const BoundParameters& predictedParams,
      const std::vector<source_link_t>& sourcelinks) const {
    ACTS_VERBOSE("Invoked VoidSourceLinkSelector");

    using CovMatrix_t = typename BoundParameters::CovMatrix_t;

    std::vector<source_link_t> sourcelinkCandidates;
    // There is either one found source link or not
    sourcelinkCandidates.reserve(1);

    double minChi2 = std::numeric_limits<double>::max();
    for (const auto& sourcelink : sourcelinks) {
      std::visit(
          [&](const auto& calibrated) {
            // type of measurement
            using meas_t =
                typename std::remove_const<typename std::remove_reference<
                    decltype(calibrated)>::type>::type;
            // measurement (local) parameter vector
            using meas_par_t = typename meas_t::ParVector_t;
            // type of projection
            using projection_t = typename meas_t::Projection_t;

            // Take the projector (measurement mapping function)
            const projection_t& H = calibrated.projector();
            // Take the parameter covariance
            const CovMatrix_t& predicted_covariance =
                *predictedParams.covariance();
            // Get the residual
            meas_par_t residual = calibrated.residual(predictedParams);
            // Get the chi2
            double chi2 = (residual.transpose() *
                           ((calibrated.covariance() +
                             H * predicted_covariance * H.transpose()))
                               .inverse() *
                           residual)
                              .eval()(0, 0);

            // Find the source link with the min chi2
            if (chi2 < minChi2 and chi2 < m_config.maxChi2) {
              minChi2 = chi2;
              if (sourcelinkCandidates.empty()) {
                sourcelinkCandidates.push_back(sourcelink);
              } else {
                sourcelinkCandidates[0] = sourcelink;
              }
            }
          },
          calibrator(sourcelink, predictedParams));
    }

    // Check if the chi2 satisfies requirement in the config
    if (not sourcelinkCandidates.empty()) {
      ACTS_VERBOSE("Minimum Chi2: " << minChi2 << " is within chi2 criteria: "
                                    << m_config.maxChi2);
    } else {
      ACTS_DEBUG("No source link candidate");
    }

    return std::move(sourcelinkCandidates);
  }

  /// The config
  Config m_config;

  /// Pointer to a logger that is owned by the parent, track finder
  std::shared_ptr<const Logger> m_logger{nullptr};

  /// Getter for the logger, to support logging macros
  const Logger& logger() const {
    assert(m_logger);
    return *m_logger;
  }
};

}  // namespace Acts
