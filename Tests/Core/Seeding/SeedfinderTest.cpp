// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Seeding/BinFinder.hpp"
#include "Acts/Seeding/InternalSeed.hpp"
#include "Acts/Seeding/InternalSpacePoint.hpp"
#include "Acts/Seeding/SeedFilter.hpp"
#include "Acts/Seeding/Seedfinder.hpp"

#include "ATLASCuts.hpp"
#include "SpacePoint.hpp"

#ifdef _OPENMP
#include "omp.h"
#else
	#define omp_get_num_threads() 0
#endif

//#include <boost/type_erasure/any_cast.hpp>

#include <chrono>
#include <fstream>
#include <iostream>
#include <sstream>
#include <utility>
#include <string>


std::vector<const SpacePoint*> readFile(std::string filename) {
  std::string line;
  int layer;
  std::vector<const SpacePoint*> readSP;

  std::ifstream spFile(filename);
  if (spFile.is_open()) {
    while (!spFile.eof()) {
      std::getline(spFile, line);
      std::stringstream ss(line);
      std::string linetype;
      ss >> linetype;
      float x, y, z, r, covr, covz;
      if (linetype == "lxyz") {
        ss >> layer >> x >> y >> z >> covr >> covz;
        r = std::sqrt(x * x + y * y);
        float f22 = covr;
        float wid = covz;
        float cov = wid * wid * .08333;
        if (cov < f22)
          cov = f22;
        if (std::abs(z) > 450.) {
          covz = 9. * cov;
          covr = .06;
        } else {
          covr = 9. * cov;
          covz = .06;
        }
        SpacePoint* sp = new SpacePoint{x, y, z, r, layer, covr, covz};
        //     if(r < 200.){
        //       sp->setClusterList(1,0);
        //     }
        readSP.push_back(sp);
      }
    }
  }
  return readSP;
}

int main(int argc, char** argv) {
	if (argc < 2) {
		std::cout << argv[0] << " file_name " << " nthreads " << std::endl;
		exit(1);
	}
	std::string input_name(argv[1]);
	int n_threads = atoi(argv[2]);
	std::vector<const SpacePoint*> spVec = readFile(input_name.c_str());
	std::cout << "input file: " << input_name << std::endl;
	std::cout << "N threads: " << n_threads << std::endl;
	std::cout << "size of read SP: " << spVec.size() << std::endl;

  Acts::SeedfinderConfig<SpacePoint> config;
  // silicon detector max
  config.rMax = 160.;
  config.deltaRMin = 5.;
  config.deltaRMax = 160.;
  config.collisionRegionMin = -250.;
  config.collisionRegionMax = 250.;
  config.zMin = -2800.;
  config.zMax = 2800.;
  config.maxSeedsPerSpM = 5;
  // 2.7 eta
  config.cotThetaMax = 7.40627;
  config.sigmaScattering = 1.00000;

  config.minPt = 500.;
  config.bFieldInZ = 0.00199724;

  config.beamPos = {-.5, -.5};
  config.impactMax = 10.;

  auto bottomBinFinder = std::make_shared<Acts::BinFinder<SpacePoint>>(
      Acts::BinFinder<SpacePoint>());
  auto topBinFinder = std::make_shared<Acts::BinFinder<SpacePoint>>(
      Acts::BinFinder<SpacePoint>());
  Acts::SeedFilterConfig sfconf;
  Acts::ATLASCuts<SpacePoint> atlasCuts = Acts::ATLASCuts<SpacePoint>();
  config.seedFilter = std::make_unique<Acts::SeedFilter<SpacePoint>>(
      Acts::SeedFilter<SpacePoint>(sfconf, &atlasCuts));
  Acts::Seedfinder<SpacePoint> seed_finder(config);

  // covariance tool, sets covariances per spacepoint as required
  auto ct = [=](const SpacePoint& sp, float, float, float) -> Acts::Vector2D {
    return {sp.covr, sp.covz};
  };

  Acts::SeedfinderState<SpacePoint> state = seed_finder.initState(
      spVec.begin(), spVec.end(), ct, bottomBinFinder, topBinFinder);
  auto start = std::chrono::system_clock::now();
///////
  std::vector<Acts::SeedfinderStateIterator<SpacePoint> > its ;
  for(Acts::SeedfinderStateIterator<SpacePoint> it = state.begin();
		  !(it == state.end()); ++it) {
		its.push_back(it);
  }
  std::cout << "number of states: " << its.size() << std::endl;

//  #pragma omp parallel num_threads(n_threads) shared(seed_finder, state, its, std::cout) default(none) 
//  {
//	#pragma omp master
//  	std::cout << "number of threads: " << omp_get_num_threads() << std::endl;
//
//	#pragma omp for
//	for(int i = 0; i < (int) its.size(); i++){
//		// Acts::SeedfinderStateIterator<SpacePoint>& it = its.at(i);
//		seed_finder.createSeedsForRegion(its[i], state);
//	}
//  }

  for (Acts::SeedfinderStateIterator<SpacePoint> it = state.begin();
        !(it == state.end()); ++it) {
     seed_finder.createSeedsForRegion(it, state);
   }

  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end - start;
  std::cout << "time to create seeds: " << elapsed_seconds.count() << std::endl;
  std::cout << "Number of regions: " << state.outputVec.size() << std::endl;
  int numSeeds = 0;
  for (auto& outVec : state.outputVec) {
    numSeeds += outVec.size();
  }
  std::cout << "Number of seeds generated: " << numSeeds << std::endl;
  for (auto& regionVec : state.outputVec) {
    for (size_t i = 0; i < regionVec.size(); i++) {
      const Acts::Seed<SpacePoint>* seed = regionVec[i].get();
      const SpacePoint* sp = seed->sp()[0];
      std::cout << sp->surface << " (" << sp->x() << ", " << sp->y() << ", " << sp->z()
                << ") ";
      sp = seed->sp()[1];
      std::cout << sp->surface << " (" << sp->x() << ", " << sp->y() << ", "
                << sp->z() << ") ";
      sp = seed->sp()[2];
      std::cout << sp->surface << " (" << sp->x() << ", " << sp->y() << ", "
                << sp->z() << ") ";
      std::cout << std::endl;
    }
  }
  return 0;
}
