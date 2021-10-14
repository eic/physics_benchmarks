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

vector<TLorentzVector > daug_cand_1,daug_cand_2;

//particles properties
auto momenta_from_reconstruction(const std::vector<eic::ReconstructedParticleData>& parts) {
  std::vector<ROOT::Math::PxPyPzEVector> momenta{parts.size()};
  std::transform(parts.begin(), parts.end(), momenta.begin(), [](const auto& part) {
    return ROOT::Math::PxPyPzEVector{part.p.x, part.p.y, part.p.z, part.energy};
    TLorentzVector cand(part.p.x, part.p.y, part.p.z, part.energy);
    if(parts.charge()>0) daug_cand_1.push_back(cand);
    else if(parts.charge()<0) daug_cand_2.push_back(cand);
    else{
      //do something for neutrals.
    }
  });
  return momenta;
}
auto getPt(const std::vector<ROOT::Math::PxPyPzEVector>& mom) {
  std::vector<double> PtVec(mom.size() );
  ROOT::Math::PxPyPzEVector beamMom = {0, 0, -18, 18};
  std::transform(mom.begin(), mom.end(), PtVec.begin(), [beamMom](const auto& part) {
    return part.Pt();
  });
  return PtVec;
}