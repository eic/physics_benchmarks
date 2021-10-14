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

//particles properties
auto momenta_from_reconstruction_plus(const std::vector<eic::ReconstructedParticleData>& parts) {
  std::vector<ROOT::Math::PxPyPzEVector> momenta{parts.size()};
  std::transform(parts.begin(), parts.end(), momenta.begin(), [](const auto& part) {
    if(part.charge<0) continue;
    return ROOT::Math::PxPyPzEVector{part.p.x, part.p.y, part.p.z, part.energy};
  });
  return momenta;
}
auto momenta_from_reconstruction_minus(const std::vector<eic::ReconstructedParticleData>& parts) {
  std::vector<ROOT::Math::PxPyPzEVector> momenta{parts.size()};
  std::transform(parts.begin(), parts.end(), momenta.begin(), [](const auto& part) {
    if(part.charge<0) {return ROOT::Math::PxPyPzEVector{part.p.x, part.p.y, part.p.z, part.energy};}
  });
  return momenta;
}

auto vector_sum = [](std::vector<ROOT::Math::PxPyPzEVector> p1, 
  std::vector<ROOT::Math::PxPyPzEVector> p2 ){
  std::vector<ROOT::Math::PxPyPzEVector> vm;
  for(auto& i1: p1){
    for(auto& i2: p2){
      //pt cut
      if(i1.Pt()<0.15||i2.Pt()<0.15) continue;
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