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
#include "TVector2.h"

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

bool sort_mom_bool(ROOT::Math::PxPyPzMVector &mom1, ROOT::Math::PxPyPzMVector &mom2) {
  return  mom1.energy() > mom2.energy(); 
}

auto sort_momenta(const std::vector<ROOT::Math::PxPyPzMVector>& mom) {
  std::vector <ROOT::Math::PxPyPzMVector> sort_mom = mom;
  sort(sort_mom.begin(), sort_mom.end(), sort_mom_bool);
  return sort_mom;
}

auto findScatElec(const std::vector<eic::ReconstructedParticleData>& parts) {
  std::vector<ROOT::Math::PxPyPzMVector> momenta{parts.size()};
  std::transform(parts.begin(), parts.end(), momenta.begin(), [](const auto& part) {
    //how to find scatElec...goodElect Cut here:
    if(part.mass<0.1||part.pid==11) return ROOT::Math::PxPyPzMVector{part.p.x, part.p.y, part.p.z, MASS_ELECTRON};
    else return ROOT::Math::PxPyPzMVector{-1e10, -1e10, -1e10, -1e10};
  });
  return momenta;
}

auto tmp_findScat(const std::vector<eic::ReconstructedParticleData>& parts, 
    std::vector<int> scat_id,
  std::vector<int> scat_source) 
{
  std::vector<ROOT::Math::PxPyPzMVector> momenta{parts.size()};
  for(auto& i1 : parts){
    if(scat_id.size()>0 
        && scat_source.size()>0
          &&i1.ID.value==scat_id[0]
            &&i1.ID.source==scat_source[0])
    {
      std::cout << "scatID = " << scat_id[0] << std::endl;
      std::cout << "scatID.source = " << scat_source[0] << std::endl;
      auto scat = ROOT::Math::PxPyPzMVector{i1.p.x, i1.p.y, i1.p.z, MASS_ELECTRON};
      momenta.push_back(scat);
      std::cout << "Eta = " << scat.Eta() << std::endl;
    }
  
  }

  // std::transform(parts.begin(), parts.end(), momenta.begin(), [](const auto& part, auto ) {
  //   if(tmp==part.ID.value) return ROOT::Math::PxPyPzMVector{part.p.x, part.p.y, part.p.z, MASS_ELECTRON};
  //   else return ROOT::Math::PxPyPzMVector{-1e10, -1e10, -1e10, -1e10};
  // });
  return momenta;
}

auto vector_sum = [](std::vector<ROOT::Math::PxPyPzMVector> p1, 
  std::vector<ROOT::Math::PxPyPzMVector> p2 ){
  std::vector<ROOT::Math::PxPyPzMVector> vm;
  for(auto& i1: p1){
    if(i1.Px()<-1e9) continue;
    for(auto& i2: p2){
      if(i2.Px()<-1e9) continue;
      //pt cut
      if(i1.Pt()<0.05||i2.Pt()<0.05) continue;
      //eta cut
      if(fabs(i1.Eta())>4.0||fabs(i2.Eta())>4.0) continue;
      vm.push_back(i1+i2);
    }
  }
  return vm;
};

//cut on phi mass region and rapidity phase space
auto getPt2OfPhi(const std::vector<ROOT::Math::PxPyPzMVector>& mom) {
  std::vector<double> PtVec(mom.size() );
  std::transform(mom.begin(), mom.end(), PtVec.begin(), [](const auto& part) {
    if(fabs(part.M()-1.019)>0.02||fabs(part.Rapidity())>3.5) return -99.;
    else return part.Pt()*part.Pt();
  });
  return PtVec;
}

auto getMass(const std::vector<ROOT::Math::PxPyPzMVector>& mom) {
  std::vector<double> massVec(mom.size() );
  std::transform(mom.begin(), mom.end(), massVec.begin(), [](const auto& part) {
    return part.M();
  });
  return massVec;
}

auto getEta(const std::vector<ROOT::Math::PxPyPzMVector>& mom) {
  std::vector<double> etaVec(mom.size() );
  std::transform(mom.begin(), mom.end(), etaVec.begin(), [](const auto& part) {
    if(part.Px()<-1e9) return -10.;
    else return part.Eta();
  });
  return etaVec;
}

auto giveme_t = [](std::vector<ROOT::Math::PxPyPzMVector> vm, 
   std::vector<ROOT::Math::PxPyPzMVector> scatElec){
  std::vector<double > t_vec;
  if(scatElec[0].Px()<-1e9||vm.size()<1) {
    t_vec.push_back(-99.);
    return t_vec;
  }
  for(auto& i1: vm){
    if(fabs(i1.Rapidity())>4.0||fabs(i1.M()-1.019)>0.02) continue;
    TVector2 sum_pt(i1.Px()+scatElec[0].Px(), i1.Py()+scatElec[0].Py());
    t_vec.push_back( sum_pt.Mod2() );
  }
  return t_vec;
};
