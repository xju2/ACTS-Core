// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#include <memory>
#include "Acts/EventData/MultiTrackParameters.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"

namespace Acts {

template <typename ChargePolicy>
class MultiBoundTrackParameters : public MultiTrackParameters<ChargePolicy>
{
public:
  /// type of covariance matrix
  using CovPtr_t = typename MultiTrackParameters<ChargePolicy>::CovPtr_t;

  /// pair of <weight,point_to_track>, when insert a trackParameter, that will
  /// allow to sort the track by weight
  using WeightedTrackPars
      = std::pair<double, std::unique_ptr<TrackParametersBase>>;
  using TrackParMap = std::multimap<double,
                                    std::unique_ptr<TrackParametersBase>,
                                    std::greater<double>>;
  using TrackParMapConstIter = TrackParMap::const_iterator;
  using TrackParMapIter      = TrackParMap::iterator;

  template <typename T = ChargePolicy,
            std::enable_if_t<std::is_same<T, ChargedPolicy>::value, int> = 0>
  MultiBoundTrackParameters(double                         weight,
                            TrackParametersBase*           pTrackBase,
                            std::shared_ptr<const Surface> surface)
    : MultiTrackParameters<ChargePolicy>(weight, pTrackBase)
    , m_pSurface(std::move(surface))
  {
    assert(m_pSurface);
  }

  template <typename T = ChargePolicy,
            std::enable_if_t<std::is_same<T, NeutralPolicy>::value, int> = 0>
  MultiBoundTrackParameters(double                         weight,
                            TrackParametersBase*           pTrackBase,
                            std::shared_ptr<const Surface> surface)
    : MultiTrackParameters<ChargePolicy>(weight, pTrackBase)
    , m_pSurface(std::move(surface))
  {
    assert(m_pSurface);
  }

  /// @brief default parameters
  template <typename T = ChargePolicy,
            std::enable_if_t<std::is_same<T, ChargedPolicy>::value, int> = 0>
  MultiBoundTrackParameters(std::shared_ptr<const Surface> surface)
    : MultiTrackParameters<ChargePolicy>(), m_pSurface(std::move(surface))
  {
    assert(m_pSurface);
  }

  /// @brief default parameters
  template <typename T = ChargePolicy,
            std::enable_if_t<std::is_same<T, NeutralPolicy>::value, int> = 0>
  MultiBoundTrackParameters(std::shared_ptr<const Surface> surface)
    : MultiTrackParameters<ChargePolicy>(), m_pSurface(std::move(surface))
  {
    assert(m_pSurface);
  }

  /// @brief copy constructor - charged/neutral
  /// @param[in] copy The source parameters
  /*unused*/
  MultiBoundTrackParameters(const MultiBoundTrackParameters<ChargePolicy>& copy)
    : MultiTrackParameters<ChargePolicy>(copy), m_pSurface(copy.m_pSurface)
  {
  }

  /// @brief move constructor - charged/neutral
  /// @param[in] other The source parameters
  /*unused*/
  MultiBoundTrackParameters(MultiBoundTrackParameters<ChargePolicy>&& other)
    : MultiTrackParameters<ChargePolicy>(std::move(other))
    , m_pSurface(std::move(other.m_pSurface))
  {
  }

  /// @brief desctructor
  ~MultiBoundTrackParameters() override = default;

  /// @brief move assignment operator - charged/netural
  /// virtual constructor for type creation without casting
  /*unused*/
  MultiBoundTrackParameters<ChargePolicy>&
  operator=(MultiBoundTrackParameters<ChargePolicy>&& rhs)
  {
    // check for self-assignment
    if (this != &rhs) {
      MultiTrackParameters<ChargePolicy>::operator=(std::move(rhs));
      m_pSurface                                  = std::move(rhs.m_pSurface);
    }
    return *this;
  }

  /// @brief set method for parameter updates
  /// obviously only allowed on non-const objects
  //
  /// @param[in] gctx is the Context object that is forwarded to the surface
  ///            for local to global coordinate transformation
  //  template <ParID_t par>
  //  void
  //  set(const GeometryContext& gctx, ParValue_t newValue, unsigned int order)
  //  {
  //    // set the parameter & update the new global position
  //    this->getParameterSet(order).template setParameter<par>(newValue);
  //    this->updateGlobalCoordinates(gctx, typename par_type<par>::type(),
  //    order);
  //  }

  /// @brief access to the reference surface
  const Surface&
  referenceSurface() const final
  {
    return *m_pSurface;
  }

  /// @brief access to the reference surface
  //  const Surface&
  //  referenceSurface(unsigned int id) const final
  //  {
  //    TrackParMapConstIter it = std::find_if(this->m_TrackList.begin(),
  //                                           this->m_TrackList.end(),
  //                                           TrackIndexFinder(id));
  //    assert(it != this->m_TrackList.end());
  //    return (*it).first.second->referenceSurface();
  //  }

  /// @brief update the reference surface
  void
  updateReferenceSurface(const ActsVectorD<3>& /*unused*/,
                         const ActsVectorD<3>& /*unused*/) final
  {
  }

  /// @brief access to the measurement frame, i.e. the rotation matrix with
  /// respect to the global coordinate system, in which the local error
  /// is described.
  ///
  /// For a curvilinear track parameterisation this is identical to the
  /// rotation matrix of the intrinsic planar surface.
  RotationMatrix3D
  referenceFrame(const GeometryContext& gctx) const final
  {
    return std::move(
        m_pSurface->referenceFrame(gctx, this->position(), this->momentum()));
  }

private:
  std::shared_ptr<const Surface> m_pSurface;
};
}  // namespace Acts
