// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// CylinderSurface.ipp, Acts project
///////////////////////////////////////////////////////////////////

inline const Vector3D CylinderSurface::rotSymmetryAxis(
    const GeometryContext& gctx) const {
  // fast access via tranform matrix (and not rotation())
  return transform(gctx).matrix().block<3, 1>(0, 2);
}

inline Intersection CylinderSurface::intersectionEstimate(
    const GeometryContext& gctx, const Vector3D& gpos, const Vector3D& gdir,
    const BoundaryCheck& bcheck, double bwdTolerance, CorrFnc correct) const {
  // create line parameters
  Vector3D lpos = gpos;
  Vector3D ldir = gdir;
  // minimize the call to transform()
  const auto& tMatrix = transform(gctx).matrix();
  Vector3D caxis = tMatrix.block<3, 1>(0, 2).transpose();
  Vector3D ccenter = tMatrix.block<3, 1>(0, 3).transpose();
  // what you need at the and
  Vector3D solution(0, 0, 0);
  double path = 0.;

  // lemma : the solver ->  catches & modifies current values
  auto solve = [&solution, &path, &lpos, &ldir, &ccenter, &caxis,
                &bwdTolerance](double R) -> IntersectionStatus {
    // check documentation for explanation
    Vector3D pc = lpos - ccenter;
    Vector3D pcXcd = pc.cross(caxis);
    Vector3D ldXcd = ldir.cross(caxis);
    double a = ldXcd.dot(ldXcd);
    double b = 2. * (ldXcd.dot(pcXcd));
    double c = pcXcd.dot(pcXcd) - (R * R);
    // and solve the qaudratic equation - todo, validity check
    detail::RealQuadraticEquation qe(a, b, c);
    // check how many solution you have
    auto istatus = IntersectionStatus::unreachable;
    if (qe.solutions != 0) {
      // now try to understand the solutions
      // both are solvable into the same direction
      if (qe.first * qe.second > 0.) {
        path =
            qe.first * qe.first < qe.second * qe.second ? qe.first : qe.second;
        istatus = path > 0. ? IntersectionStatus::reachable
                            : (path * path < bwdTolerance * bwdTolerance and
                                       bwdTolerance > 0.
                                   ? IntersectionStatus::overstepped
                                   : IntersectionStatus::unreachable);
      } else {
        istatus = IntersectionStatus::reachable;
        // different directions, in principle chose the forward one
        path = qe.first > 0. ? qe.first : qe.second;
        double otherPath = qe.first < 0. ? qe.first : qe.second;
        double otherPath2 = otherPath * otherPath;
        // otherPath is allowed to overwrite path if it's within Bwd tolerance
        // and in absolute terms smaller then the path
        if (otherPath2 < bwdTolerance * bwdTolerance and
            otherPath2 < path * path and bwdTolerance > 0.) {
          path = otherPath;
          istatus = IntersectionStatus::overstepped;
        }
      }
      istatus = (path * path < s_onSurfaceTolerance * s_onSurfaceTolerance)
                    ? IntersectionStatus::onSurface
                    : istatus;

      // return the solution
      solution = lpos + path * ldir;
      // is valid if it goes into the right direction
    }
    return istatus;
  };

  // solve for radius R
  double R = bounds().r();
  auto istatus = solve(R);
  // if configured, correct and solve again
  if (correct && correct(lpos, ldir, path)) {
    istatus = solve(R);
  }
  // Evaluate boundaries if necessary, since the solution was done in global
  // frame a global to local is necessary for the boundary check
  if (bcheck && istatus != IntersectionStatus::unreachable) {
    istatus = isOnSurface(gctx, solution, gdir, bcheck)
                  ? istatus
                  : IntersectionStatus::missed;
  }
  // now return
  return Intersection(solution, path, istatus);
}
