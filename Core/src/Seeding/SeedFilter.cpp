#include "ACTS/Seeding/SeedFilter.hpp"
#include <utility>

namespace Acts{
namespace Seeding{
  //constructor
  SeedFilter::SeedFilter(SeedFilterConfig config,
                         std::shared_ptr<IQualityTool> qualityTool)
                         :m_cfg (config),
                          m_qualityTool (qualityTool){}

  //destructor
  SeedFilter::~SeedFilter(){}

  // function to filter seeds based on all seeds with same bottom- and middle-spacepoint.
  // return vector must contain quality of each seed 
  std::vector<std::pair<float, std::shared_ptr<InternalSeed> > >
  SeedFilter::filterSeeds_2SpFixed(std::shared_ptr<SPForSeed> bottomSP,
                                   std::shared_ptr<SPForSeed> middleSP,
                                   std::vector<std::shared_ptr<SPForSeed >>& topSpVec,
                                   std::vector<float>& invHelixRadiusVec,
                                   std::vector<float>& impactParametersVec,
                                   float zOrigin){
  
    std::vector<std::pair<float, std::shared_ptr<InternalSeed> > > selectedSeeds;
  
    // if two compatible seeds with high distance in r are found, compatible seeds span 5 layers
    // -> very good seed
    std::vector<float> compatibleSeedRadii;

    for(int i = 0; i < topSpVec.size(); i++){

      float invHelixRadius = invHelixRadiusVec.at(i);
      float lowerLimitCurv = invHelixRadius - m_cfg.deltaInvHelixRadius; 
      float upperLimitCurv = invHelixRadius + m_cfg.deltaInvHelixRadius; 
      float currentTop_r   = topSpVec.at(i)->radius();
      float impact         = impactParametersVec.at(i);

      float quality = -(impact * m_cfg.impactQualityFactor);
      quality += m_qualityTool->seedQuality(bottomSP, middleSP, topSpVec.at(i));
      for (int j = 0; j < topSpVec.size(); j++){
        if (i == j) continue;
        // compared top SP should have at least deltaRMin distance
        float otherTop_r = topSpVec.at(j)->radius();
        float deltaR = currentTop_r - otherTop_r;
        if (std::abs(deltaR) < m_cfg.deltaRMin)   continue;
        // curvature difference within limits?
        // TODO: how much slower than sorting all vectors by curvature
        // and breaking out of loop? i.e. is vector size large (e.g. in jets?)
        if (invHelixRadiusVec.at(j) < lowerLimitCurv) continue;
        if (invHelixRadiusVec.at(j) > upperLimitCurv) continue;
        bool newCompSeed = true;
        for(float previousRadius : compatibleSeedRadii){
          // original ATLAS code uses higher min distance for 2nd found compatible seed (20mm instead of 5mm)
          if(std::abs(previousRadius - otherTop_r) < m_cfg.deltaRMin) {newCompSeed = false; break;}
        }
        if(newCompSeed)
        {
          compatibleSeedRadii.push_back(otherTop_r);
          quality+= m_cfg.compatSeedQuality;
        }
        if(compatibleSeedRadii.size() >= m_cfg.compatSeedLimit) break;
      }
      // discard low quality seeds
      if (!m_qualityTool->passesQualityCut(quality, bottomSP, middleSP, topSpVec.at(i))) continue;
      selectedSeeds.push_back(std::make_pair(quality, std::make_shared<InternalSeed>(bottomSP,middleSP,topSpVec.at(i),zOrigin)));
      }
    return selectedSeeds;
  }



  // after creating all seeds with a common middle space point, filter again
  std::vector<std::shared_ptr<InternalSeed> >
  SeedFilter::filterSeeds_1SpFixed(std::vector<std::pair<float,std::shared_ptr<InternalSeed > > >& seedsPerSpM){

    //sort by quality and iterate only up to configured max number of seeds per middle SP
    std::sort(seedsPerSpM.begin(),seedsPerSpM.end(),comQuality());
    int maxSeeds = seedsPerSpM.size();
    if(maxSeeds > m_cfg.maxSeedsPerSpM){
      maxSeeds = m_cfg.maxSeedsPerSpM + 1;
    }
    auto itBegin = seedsPerSpM.begin();
    auto it = seedsPerSpM.begin();
  
    std::vector<std::shared_ptr<InternalSeed> > filteredSeeds;
    // default filter removes the last seeds if maximum amount exceeded
    // ordering by quality by filterSeeds_2SpFixed means these are the lowest quality seeds
    for(; it<itBegin+maxSeeds; ++it) {
      std::shared_ptr<InternalSeed> internalSeed       = (*it).second;
      filteredSeeds.push_back(internalSeed);
    }
    return filteredSeeds;
  }



  std::vector<std::shared_ptr<Seed> >
  SeedFilter::filterSeeds_byRegion(std::vector<std::shared_ptr<InternalSeed> >& regionSeeds){

    std::vector<std::shared_ptr<Seed> > outputSeeds;
    for(auto s : regionSeeds){
      Seed outSeed = Seed(s->spacepoint0()->spacepoint,
          s->spacepoint1()->spacepoint,
          s->spacepoint2()->spacepoint,
          s->z());
      outputSeeds.push_back(std::make_shared<Seed>(outSeed ) );
    }
    return outputSeeds;
  }
}//namespace Seeding
}//namespace Acts