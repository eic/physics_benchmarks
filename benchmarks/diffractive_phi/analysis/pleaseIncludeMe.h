#include "common_bench/benchmark.h"
#include "common_bench/mt.h"
#include "common_bench/util.h"
#include "common_bench/plot.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <utility>

#include "ROOT/RDataFrame.hxx"
#include <TH1D.h>
#include <TFitResult.h>
#include <TRandom3.h>
#include <TCanvas.h>

#include "TFile.h"
#include "TLorentzVector.h"

#include "fmt/color.h"
#include "fmt/core.h"

#include "nlohmann/json.hpp"
#include "eicd/InclusiveKinematicsData.h"
#include "eicd/ReconstructedParticleData.h"

#define PI            3.1415926
#define MASS_PION     0.13957
#define MASS_KAON     0.493667
#define MASS_AU197    183.45406466643374

//particles properties
auto momenta_from_reconstruction_plus(const std::vector<eic::ReconstructedParticleData>& parts) {
  std::vector<ROOT::Math::PxPyPzMVector> momenta{parts.size()};
  std::transform(parts.begin(), parts.end(), momenta.begin(), [](const auto& part) {
    if(part.charge>0){
      return ROOT::Math::PxPyPzMVector{part.p.x, part.p.y, part.p.z, MASS_KAON};
    }
    else{
      return ROOT::Math::PxPyPzMVector{-1e10, -1e10, -1e10, -1e10};
    }
  });
  return momenta;
}
auto momenta_from_reconstruction_minus(const std::vector<eic::ReconstructedParticleData>& parts) {
  std::vector<ROOT::Math::PxPyPzMVector> momenta{parts.size()};
  std::transform(parts.begin(), parts.end(), momenta.begin(), [](const auto& part) {
     if(part.charge<0){
      return ROOT::Math::PxPyPzMVector{part.p.x, part.p.y, part.p.z, MASS_KAON};
    }
    else{
      return ROOT::Math::PxPyPzMVector{-1e10, -1e10, -1e10, -1e10};
    }
  });
  return momenta;
}

auto vector_sum = [](std::vector<ROOT::Math::PxPyPzEVector> p1, 
  std::vector<ROOT::Math::PxPyPzEVector> p2 ){
  std::vector<ROOT::Math::PxPyPzEVector> vm;
  for(auto& i1: p1){
    if(i1.Px()<-1e9) continue;
    for(auto& i2: p2){
      if(i2.Px()<-1e9) continue;
      //pt cut
      // if(i1.Pt()<0.15||i2.Pt()<0.15) continue;
      //eta cut
      //...
      vm.push_back(i1+i2);
    }
  }
  return vm;
};

auto getPt2(const std::vector<ROOT::Math::PxPyPzEVector>& mom) {
  std::vector<double> PtVec(mom.size() );
  std::transform(mom.begin(), mom.end(), PtVec.begin(), [](const auto& part) {
    return part.Pt()*part.Pt();
  });
  return PtVec;
}
auto getMass(const std::vector<ROOT::Math::PxPyPzEVector>& mom) {
  std::vector<double> massVec(mom.size() );
  std::transform(mom.begin(), mom.end(), massVec.begin(), [](const auto& part) {
    return part.M();
  });
  return massVec;
}


const Double_t dxbin = (0.17 - 0.13) / 40; // Bin-width