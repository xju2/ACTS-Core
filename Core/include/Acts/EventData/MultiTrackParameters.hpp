// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#include <map>
#include <type_traits>
#include "Acts/EventData/TrackParametersBase.hpp"
#include "Acts/EventData/detail/coordinate_transformations.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/GeometryContext.hpp"

namespace Acts {

/// @class MultiTrackParameters
///
/// @brief base class for a multi set of track parameters
///
/// this class represents a multi set of track parameters (used in Gsf)

///
/// @tparam ChargePolicy type for distinguishing charged and neutral
/// tracks/particles (must be either ChargedPolicy or NeutralPolicy)
template <class ChargePolicy>
class MultiTrackParameters : public TrackParametersBase
{
  static_assert(std::is_same<ChargePolicy, ChargedPolicy>::value
                    or std::is_same<ChargePolicy, NeutralPolicy>::value,
                "ChargePolicy must either be 'Acts::ChargedPolicy' or "
                "'Acts::NeutralPolicy");

public:
  // public typedef's

  /// vector type for stored track parameters
  using ParVector_t = typename TrackParametersBase::ParVector_t;

  /// type of covariance matrix
  using CovMatrix_t = typename TrackParametersBase::CovMatrix_t;

  /// type for unique pointer to covariance matrix
  using CovPtr_t = std::unique_ptr<const CovMatrix_t>;

  /// pair of (weight,pointer_to_track), when insert a trackParameter, sort the track by weight
  using WeightedTrackPars
      = std::pair<double, std::unique_ptr<TrackParametersBase>>;
  using TrackParMap = std::multimap<double, std::unique_ptr<TrackParametersBase>, std::greater<double> >;
  using TrackParMapConstIter = TrackParMap::const_iterator;
  using TrackParMapIter      = TrackParMap::iterator;

  /// @brief default virtual destructor
  ~MultiTrackParameters() override = default;

  /// the position()/direction()/momentum() method
  /// @brief get weighted combination of position all contains two format,
  /// one is to get the combine of the parameters
  ActsVectorD<3>
  position() const final
  {
    Vector3D pos = Vector3D(0, 0, 0);
    for (const auto& trackMapIndex : m_TrackList) {
      pos += trackMapIndex.first * trackMapIndex.second->position();
    }
    return pos;
  }

  /// @copydoc TrackParametersBase::momentum
  ActsVectorD<3>
  momentum() const final
  {
    Vector3D mom = Vector3D(0, 0, 0);
    for (const auto& trackMapIndex : m_TrackList) {
      mom += trackMapIndex.first * trackMapIndex.second->momentum();
    }
    return mom;
  }

  /// @brief equality operator
  ///
  /// @return @c true if all single track parameters are same
  bool
  operator==(const TrackParametersBase& /*unused*/) const override
  {
    return true;
  }

  /// @copydoc TrackParametersBase::charge
  /*temporary*/
  double
  charge() const final
  {
    return m_TrackList.begin()->second->charge();
  }

  /// @copydoc TrackParametersBase::getParameterSet
  /// currently get first component parameters
  /// @to do : return combination
  /// const
  /*temporary*/
  const FullParameterSet&
  getParameterSet() const final
  {
    return m_TrackList.begin()->second->getParameterSet();
  }
  /// copy base class
  /// currently get first component parameters
  /// @to do : return combination
  /// writable
  /*temporary*/
  FullParameterSet&
  getParameterSet() final
  {
    return m_TrackList.begin()->second->getParameterSet();
  }

  /// @brief get single component
  /// temp
  ParVector_t
  parameters() const
  {
    return getParameterSet().getParameters();
  }

  /// @brief update mom pos
  void updateMom(ActsVectorD<3>& /*unused*/) final {}

  void updatePos(ActsVectorD<3>& /*unused*/) final {}

  /// @brief return the size of track list
  size_t
  size() const
  {
    return m_TrackList.size();
  }

  /// @brief append a component to the track list
  virtual void
  append(double weight, TrackParametersBase* pTrackBase)
  {
//    WeightedTrackPars weight_parameter(
 //       weight, std::unique_ptr<TrackParametersBase>(pTrackBase) );
    m_TrackList.insert(std::make_pair(weight, std::unique_ptr<TrackParametersBase>(pTrackBase) ));
  }

  /*
  virtual void
  append(std::pair<>double weight, TrackParametersBase* pTrackBase)
  {
    WeightedTrackPars weight_parameter(
        weight, std::unique_ptr<TrackParametersBase>(pTrackBase) );
    m_TrackList.insert(std::make_pair(std::move(weight_parameter), size()));
  }
  */

  /// @brief get the trackList
  /// const
  const TrackParMap&
  getTrackList() const
  {
    return m_TrackList;
  }
  /// @brief get the trackList
  /// writable
  TrackParMap&
  getTrackList()
  {
    return m_TrackList;
  }

protected:
  /// @brief standard constructor for track parameters of charged particles
  ///
  /// @param cov unique pointer to covariance matrix (nullptr is accepted)
  /// @param parValues vector with parameter values
  /// @param position 3D vector with global position
  /// @param momentum 3D vector with global momentum
  template <typename T = ChargePolicy,
            std::enable_if_t<std::is_same<T, ChargedPolicy>::value, int> = 0>
  MultiTrackParameters(double weight, TrackParametersBase* pTrackBase)
    : TrackParametersBase()
  {
    m_TrackList.insert(std::make_pair( weight,std::unique_ptr<TrackParametersBase>(pTrackBase)));
  }

  /// @brief standard constructor for track parameters of neutral particles
  ///
  /// @param cov unique pointer to covariance matrix (nullptr is accepted)
  /// @param parValues vector with parameter values
  /// @param position 3D vector with global position
  /// @param momentum 3D vector with global momentum
  template <typename T = ChargePolicy,
            std::enable_if_t<std::is_same<T, NeutralPolicy>::value, int> = 0>
  MultiTrackParameters(double weight, TrackParametersBase* pTrackBase)
    : TrackParametersBase()
  {
    m_TrackList.insert(std::make_pair( weight,std::unique_ptr<TrackParametersBase>(pTrackBase)));
    //WeightedTrackPars weight_parameter(
    //    weight, std::unique_ptr<TrackParametersBase>(pTrackBase));
    //m_TrackList.insert(std::make_pair(std::move(weight_parameter), 0));
  }

  /// @brief default constructor
  template <typename T = ChargePolicy,
            std::enable_if_t<std::is_same<T, ChargedPolicy>::value, int> = 0>
  MultiTrackParameters() : TrackParametersBase()
  {
  }

  /// @brief default constructor
  template <typename T = ChargePolicy,
            std::enable_if_t<std::is_same<T, NeutralPolicy>::value, int> = 0>
  MultiTrackParameters() : TrackParametersBase()
  {
  }

  /// @brief default copy constructor
  /*unused*/
  MultiTrackParameters(const MultiTrackParameters<ChargePolicy>& copy)
      = default;

  /// @brief default move constructor
  MultiTrackParameters(MultiTrackParameters<ChargePolicy>&& copy) = default;

  /// @brief copy assignment operator
  ///
  /// @param rhs object to be copied
  /*unused*/
  MultiTrackParameters<ChargePolicy>&
  operator=(const MultiTrackParameters<ChargePolicy>& /*unused*/)
  {
  }

  /// @brief move assignment operator
  ///
  /// @param rhs object to be movied into `*this`
  MultiTrackParameters<ChargePolicy>&
  operator=(MultiTrackParameters<ChargePolicy>&& /*unused*/)
  {
  }

  /// @brief update global momentum from current parameter values
  ///
  /// @param[in] gctx is the Context object that is forwarded to the surface
  ///            for local to global coordinate transformation
  ///
  /// @note This function is triggered when called with an argument of a type
  ///       different from Acts::local_parameter
//  template <typename T>
//  void
//  updateGlobalCoordinates(const GeometryContext& /*gctx*/,
//                          const T& /*unused*/,
//                          unsigned int id)
//  {
//    auto vMomentum
//        = detail::coordinate_transformation::parameters2globalMomentum(
//            getParameterSet(id).getParameters());
//    TrackParMapIter it = std::find_if(
//        m_TrackList.begin(), m_TrackList.end(), TrackIndexFinder(id));
//    assert(it != m_TrackList.end());
//    (*it).first.second->updateMom(vMomentum);
//  }
//
//  /// @brief update global position from current parameter values
//  ///
//  /// @note This function is triggered when called with an argument of a type
//  /// Acts::local_parameter
//  void
//  updateGlobalCoordinates(const GeometryContext& gctx,
//                          const local_parameter& /*unused*/,
//                          unsigned int id)
//  {
//    auto vPosition
//        = detail::coordinate_transformation::parameters2globalPosition(
//            gctx,
//            getParameterSet(id).getParameters(),
//            this->referenceSurface(id));
//    TrackParMapIter it = std::find_if(
//        m_TrackList.begin(), m_TrackList.end(), TrackIndexFinder(id));
//    assert(it != m_TrackList.end());
//    (*it).first.second->updatePos(vPosition);
//  }

  virtual const Surface&
  referenceSurface() const = 0;

  TrackParMap m_TrackList;
};
}  // namespace Acts
