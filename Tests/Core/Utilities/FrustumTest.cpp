// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#define BOOST_TEST_MODULE Frustum Test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <boost/test/output_test_stream.hpp>

#include "Acts/Utilities/Frustum.hpp"
#include "Acts/Utilities/PlyHelper.hpp"

#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"

using boost::test_tools::output_test_stream;

namespace Acts {
namespace Test {

BOOST_AUTO_TEST_SUITE(Utilities)
BOOST_AUTO_TEST_CASE(frustum_construction) {
  output_test_stream output;

  using Frustum2f2 = Frustum<float, 2, 2>;
  Frustum2f2 fr({1, 0}, {0, 2}, M_PI / 2.);

  BOOST_CHECK_EQUAL(fr.origin(), Vector2F(1, 0));
  CHECK_CLOSE_ABS(fr.dir(), Vector2F(0, 1), 1e-6);

  const auto& normals = fr.normals();
  BOOST_CHECK_EQUAL(normals.size(), 3);

  fr.svg(output, 200, 200);
  BOOST_CHECK(!output.is_empty(true));

  using Frustum3f3 = Frustum<float, 3, 3>;
  Frustum3f3 fr33({1, 0, 0}, {0, 2, 1}, M_PI / 2.);

  BOOST_CHECK_EQUAL(fr33.origin(), Vector3F(1, 0, 0));
  CHECK_CLOSE_ABS(fr33.dir(), Vector3F(0, 2, 1).normalized(), 1e-6);

  const auto& normals33 = fr33.normals();
  BOOST_CHECK_EQUAL(normals33.size(), 4);

  PlyHelper<float> hlp;
  // compile call to draw, does not actually test anything
  fr33.draw(hlp);

  using Frustum3f4 = Frustum<float, 3, 4>;
  Frustum3f4 fr34({1, 0, 0}, {0, 2, 1}, M_PI / 2.);

  BOOST_CHECK_EQUAL(fr34.origin(), Vector3F(1, 0, 0));
  CHECK_CLOSE_ABS(fr34.dir(), Vector3F(0, 2, 1).normalized(), 1e-6);

  const auto& normals34 = fr34.normals();
  BOOST_CHECK_EQUAL(normals34.size(), 5);

  fr34.draw(hlp);
}
BOOST_AUTO_TEST_SUITE_END()

}  // namespace Test
}  // namespace Acts
