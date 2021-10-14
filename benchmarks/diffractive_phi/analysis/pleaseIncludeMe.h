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
#define MASS_ELECTRON 0.00051
#define MASS_PION     0.13957
#define MASS_KAON     0.493667
#define MASS_AU197    183.45406466643374


int which_vm = 1;
double vm_pid[3]={113,333,443};
double vm_mass[3]={0.77545,1.019,3.0969};
double vm_mass_width[3]={0.15,0.02,0.03};
double vm_daug_pid[3]={211,321,11};
double vm_daug_mass[3]={MASS_PION,MASS_KAON,MASS_ELECTRON};

//resolution.
auto combinatorial_diff_ratio = [] (
    const ROOT::VecOps::RVec<float>& v1,
    const ROOT::VecOps::RVec<float>& v2
) {
  std::vector<float> v;
  for (auto& i1: v1) {
    for (auto& i2: v2) {
      if (i1 != 0) {
        v.push_back((i1-i2)/i1);
      }
    }
  }
  return v;
};

auto giveme_resolution = [] (
    const std::vector<double> v1,
    const std::vector<double> v2
) {
  std::vector<float> v;
  for (auto& i1: v1) {
    if (v2.size()>0) {
      v.push_back((i1-v2[0])/i1);
    }
    else{
      v.push_back(-99);
    }
  }
  return v;
};

auto matchVectKine(ROOT::Math::PxPyPzMVector v1, ROOT::Math::PxPyPzMVector v2){
  TLorentzVector v1_L(v1.Px(),v1.Py(),v1.Pz(),v1.E());
  TLorentzVector v2_L(v2.Px(),v2.Py(),v2.Pz(),v2.E());

  if(v1_L.DeltaR(v2_L)>0.3e-1) return false;
  else return true;
}

auto scatID_cand_value = [](const ROOT::VecOps::RVec<int>& x){
  std::vector<int> value;
  for(auto& i1 : x) {value.push_back( i1 );}
  return value;
};

auto momenta_from_reconstruction_plus(const std::vector<eic::ReconstructedParticleData>& parts) {
  std::vector<ROOT::Math::PxPyPzMVector> momenta{parts.size()};
  std::transform(parts.begin(), parts.end(), momenta.begin(), [](const auto& part) {
    return ROOT::Math::PxPyPzEVector{part.p.x, part.p.y, part.p.z, part.energy};
    TLorentzVector cand(part.p.x, part.p.y, part.p.z, part.energy);
    if(part.charge>0) daug_cand_1.push_back(cand);
    else if(part.charge<0) daug_cand_2.push_back(cand);
    else{
      return ROOT::Math::PxPyPzMVector{-1e10, -1e10, -1e10, -1e10};
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
