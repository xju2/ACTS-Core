// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Propagator/detail/condition_list_implementation.hpp"
#include "Acts/Propagator/detail/condition_signature_check.hpp"
#include "Acts/Utilities/detail/Extendable.hpp"
#include "Acts/Utilities/detail/MPL/has_duplicates.hpp"
#include "Acts/Utilities/detail/MPL/type_collector.hpp"

#include <boost/hana/type.hpp>
#include <boost/hana/unpack.hpp>

namespace hana = boost::hana;

namespace Acts {

/// @brief ConditionList object to be used in the propagation
///
/// The condition list is a list of structs or classes that
/// is during the propagation step and can trigger a chained
/// conditional return to the caller, e.g. that could be an aborter
///
/// It can (optionally) depend on a result of an actor from
/// the actor list.
template <typename... conditions_t>
struct ConditionList : public detail::Extendable<conditions_t...> {
 private:
  static_assert(not detail::has_duplicates_v<conditions_t...>,
                "same aborter type specified several times");

  using detail::Extendable<conditions_t...>::tuple;

 public:
  // This uses the type collector
  using result_type = typename decltype(hana::unpack(
      detail::type_collector_t<detail::action_type_extractor, conditions_t...>,
      hana::template_<ConditionList>))::type;

  using detail::Extendable<conditions_t...>::get;

  /// Default constructor
  ConditionList() = default;

  /// Default copy constructor
  ///
  /// @param conditions The condition list at source
  ConditionList(const ConditionList<conditions_t...>& conditions) = default;

  /// Default move constructor
  ///
  /// @param conditions The condition list at source
  ConditionList(ConditionList<conditions_t...>&& conditions) = default;

  /// Default move assignment operator
  ///
  /// @param conditions The condition list at source
  ConditionList<conditions_t...>& operator=(
      const ConditionList<conditions_t...>& conditions) = default;

  /// Default move assignment operator
  ///
  /// @param conditions The condition list at source
  ConditionList<conditions_t...>& operator=(
      ConditionList<conditions_t...>&& conditions) = default;

  /// Constructor from tuple
  ///
  /// @param extensions Source extensions tuple
  ConditionList(const std::tuple<conditions_t...>& conditions)
      : detail::Extendable<conditions_t...>(conditions) {}

  /// Constructor from tuple move
  ///
  /// @param extensions Source extensions tuple
  ConditionList(std::tuple<conditions_t...>&& conditions)
      : detail::Extendable<conditions_t...>(std::move(conditions)) {}

  /// Append new entries and return a new condition
  template <typename... appendices_t>
  ConditionList<conditions_t..., appendices_t...> append(
      appendices_t... aps) const {
    auto catTuple =
        std::tuple_cat(tuple(), std::tuple<appendices_t...>(aps...));
    return ConditionList<conditions_t..., appendices_t...>(std::move(catTuple));
  }

  /// This is the call signature for the condition list, it broadcasts the call
  /// to the tuple() memembers of the list
  ///
  /// @tparam result_t is the result type from a certain action
  /// @tparam propagator_state_t is the state type of the propagator
  /// @tparam stepper_t Type of the stepper
  ///
  /// @param [in] result is the result object from a certain action
  /// @param [in,out] state is the state object from the propagator
  /// @param [in] stepper Stepper used for the propagation
  template <typename result_t, typename propagator_state_t, typename stepper_t>
  bool operator()(const result_t& result, propagator_state_t& state,
                  const stepper_t& stepper) const {
    // clang-format off
    static_assert(detail::all_of_v<detail::condition_signature_check_v<
                        conditions_t, 
                        propagator_state_t, stepper_t>...>,
                  "not all conditions support the specified input");
    // clang-format on

    return detail::condition_list_impl<conditions_t...>::check(tuple(), result,
                                                               state, stepper);
  }
};

}  // namespace Acts