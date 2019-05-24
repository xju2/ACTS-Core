// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// ConeSurface.ipp, Acts project
///////////////////////////////////////////////////////////////////

inline Intersection ConeSurface::intersectionEstimate(
    const GeometryContext& gctx, const Vector3D& gpos, const Vector3D& gdir,
    const BoundaryCheck& bcheck, double bwdTolerance, CorrFnc correct) const {
  // Transform into the local frame of the surface, context resolved
  const Transform3D& ctxtrans = transform(gctx);
  const Transform3D invTrans = ctxtrans.inverse();
  Vector3D point1 = invTrans * gpos;
  Vector3D direction = invTrans.linear() * gdir;

  // what you need at the and
  Vector3D solution(0, 0, 0);
  double path = 0.;

  double tan2Alpha = bounds().tanAlpha() * bounds().tanAlpha();

  // break condition for the loop
  // bool correctionDone = false;

  // lemma : the solver ->  catches & modifies current values
  auto solve = [&solution, &path, &point1, &direction, &tan2Alpha,
                &bwdTolerance]() -> IntersectionStatus {
    auto istatus = IntersectionStatus::unreachable;

    // see the header for the formula derivation
    double A = direction.x() * direction.x() + direction.y() * direction.y() -
               tan2Alpha * direction.z() * direction.z(),
           B = 2 * (direction.x() * point1.x() + direction.y() * point1.y() -
                    tan2Alpha * direction.z() * point1.z()),
           C = point1.x() * point1.x() + point1.y() * point1.y() -
               tan2Alpha * point1.z() * point1.z();
    if (A == 0.) {
      A += 1e-16;  // avoid division by zero
    }

    // quadratic equation solver
    detail::RealQuadraticEquation solns(A, B, C);
    // Only continue if you have at least one solution
    if (solns.solutions != 0) {
      // take t1 first
      double t1 = solns.first;
      Vector3D soln1Loc(point1 + t1 * direction);

      // there's only one solution for this
      if (solns.solutions == 1) {
        // set the solution
        solution = soln1Loc;
        path = t1;
        // possible statii: reachable, overstepped, unreachable
        // (on surface will be determined generally later)
        istatus = t1 > 0. ? IntersectionStatus::reachable
                          : (t1 * t1 < bwdTolerance * bwdTolerance)
                                ? IntersectionStatus::overstepped
                                : istatus;
      } else if (solns.solutions == 2) {
        // get the second solution
        double t2 = solns.second;
        Vector3D soln2Loc(point1 + t2 * direction);

        // decide between t1 and t2 - both same direction
        if (t1 * t2 > 0.) {
          if (t1 * t1 < t2 * t2) {
            path = t1;
            solution = soln1Loc;
          } else {
            path = t2;
            solution = soln2Loc;
          }
          // status is this case is: reachable, overstepped, unreachable
          // (on surface will be determined generally later)
          istatus = t1 > 0. ? IntersectionStatus::reachable
                            : (t1 * t1 < bwdTolerance * bwdTolerance)
                                  ? IntersectionStatus::overstepped
                                  : istatus;
        } else {
          // we should be reachable in any case
          istatus = IntersectionStatus::reachable;
          // different directions, in principle chose the forward one
          path = t1 > 0. ? t1 : t2;
          double otherPath = t1 < 0. ? t1 : t2;
          double otherPath2 = otherPath * otherPath;
          // otherPath is allowed to overwrite path if it's within Bwd tolerance
          // and in absolute terms smaller then the path
          if (otherPath2 < bwdTolerance * bwdTolerance &&
              otherPath2 < path * path) {
            path = otherPath;
            istatus = IntersectionStatus::overstepped;
          }
        }
      }
      // overwrite with on surface status
      istatus = (path * path < s_onSurfaceTolerance * s_onSurfaceTolerance)
                    ? IntersectionStatus::onSurface
                    : istatus;
    }
    return istatus;
  };

  auto istatus = solve();
  // if configured, correct and solve again
  if (correct && correct(point1, direction, path)) {
    istatus = solve();
  }

  // transform back into the surface frame
  solution = ctxtrans * solution;

  // Evaluate boundaries if necessary, since the solution was done in global
  // frame a global to local is necessary for the boundary check
  if (bcheck && istatus != IntersectionStatus::unreachable) {
    istatus = isOnSurface(gctx, solution, gdir, bcheck)
                  ? istatus
                  : IntersectionStatus::missed;
  }

  // set the result navigation direction
  return Intersection(solution, path, istatus);
}
