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

auto findScatElec(const std::vector<eic::ReconstructedParticleData>& parts, 
      std::vector<int> scat_id,
    std::vector<int> scat_source) 
{
  std::vector<ROOT::Math::PxPyPzMVector> momenta;
  for(auto& i1 : parts){
    if(scat_id.size()>0 
        && scat_source.size()>0
          &&i1.ID.value==scat_id[0]
            &&i1.ID.source==scat_source[0])
    {
      auto scat = ROOT::Math::PxPyPzMVector{i1.p.x, i1.p.y, i1.p.z, MASS_ELECTRON};
      momenta.push_back(scat);
    }
    else{
      momenta.push_back(ROOT::Math::PxPyPzMVector{-1e10, -1e10, -1e10, -1e10});
    }
  
  }
  return momenta;
}

auto findScatElecMC(const std::vector<dd4pod::Geant4ParticleData>& parts)
{
  std::vector<ROOT::Math::PxPyPzMVector> momenta;
  for(auto& i1 : parts){
    if(i1.genStatus==1&&i1.pdgID==11) {
      momenta.push_back(ROOT::Math::PxPyPzMVector{i1.ps.x,i1.ps.y,i1.ps.z,i1.mass});
    }
    else {momenta.push_back(ROOT::Math::PxPyPzMVector{-1e10, -1e10, -1e10, -1e10});}
  }
  return momenta;
}

auto findGamma(const std::vector<eic::ReconstructedParticleData>& parts){
  std::vector<ROOT::Math::PxPyPzMVector> momenta;
  for(auto& i1 : parts){
    if(i1.charge==0){
      momenta.push_back(ROOT::Math::PxPyPzMVector{i1.p.x,i1.p.y,i1.p.z,0});
    }
  }
  return momenta;
}

auto findGammaMC(const std::vector<dd4pod::Geant4ParticleData>& parts) {
  std::vector<ROOT::Math::PxPyPzMVector> momenta;
  for(auto& i1 : parts){
    if(i1.genStatus==1&&i1.pdgID==22) {
      momenta.push_back(ROOT::Math::PxPyPzMVector{i1.ps.x,i1.ps.y,i1.ps.z,i1.mass});
    }
    else {momenta.push_back(ROOT::Math::PxPyPzMVector{-1e10, -1e10, -1e10, -1e10});}
  }
  return momenta;
}

auto findScatProton(const std::vector<eic::ReconstructedParticleData>& FF){
  std::vector<ROOT::Math::PxPyPzMVector> momenta;
  for(auto& i1 : FF){
    if(i1.charge==1){
      momenta.push_back(ROOT::Math::PxPyPzMVector{i1.p.x,i1.p.y,i1.p.z,i1.mass});
    }
  }
  return momenta;
}

auto resolution_MC_match_REC(const std::vector<ROOT::Math::PxPyPzMVector> MC, 
  const std::vector<ROOT::Math::PxPyPzMVector> REC)
{
  std::vector<double > resolution;
  for(auto& i1:MC){
    double res = -1.e-10;
    for(auto& i2:REC){
      if(matchVectKine(i1,i2)&&fabs(i2.M()-i1.M())<vm_mass_width[which_vm]){
        res = (i1.Pt()-i2.Pt())/i2.Pt();
      } 
    }
    resolution.push_back(res);
  }
  return resolution;
}

auto getMass(const std::vector<ROOT::Math::PxPyPzMVector>& mom) {
  std::vector<double> massVec(mom.size() );
  std::transform(mom.begin(), mom.end(), massVec.begin(), [](const auto& part) {
    return part.M();
  });
  return massVec;
}

auto getPt(const std::vector<ROOT::Math::PxPyPzMVector>& mom) {
  std::vector<double> ptVec(mom.size() );
  std::transform(mom.begin(), mom.end(), ptVec.begin(), [](const auto& part) {
    if(part.Px()<-1e9) return -99.;
    else return part.Pt();
  });
  return ptVec;
}


auto getEta(const std::vector<ROOT::Math::PxPyPzMVector>& mom) {
  std::vector<double> etaVec;
  for(auto& i1:mom){
    double eta = i1.Eta();
    if(i1.Px()<-1e9){eta=-10.;}
    etaVec.push_back(eta);
  }
  return etaVec;
}

auto getPhi(const std::vector<ROOT::Math::PxPyPzMVector>& mom) {
  std::vector<double> phiVec;
  for(auto& i1:mom){
    double phi = i1.Phi();
    if(i1.Px()<-1e9){phi=-10.;}
    phiVec.push_back(phi);
  }
  return phiVec;
}


auto giveme_t_MC(const std::vector<dd4pod::Geant4ParticleData>& parts){
  std::vector<double > t_vec;
  TLorentzVector pIn,pOut;
  for(auto& i1 : parts){
    if(i1.genStatus==4&&i1.pdgID==2212) {
      TVector3 pIn_v3(i1.ps.x,i1.ps.y,i1.ps.z);
      pIn.SetVectM(pIn_v3,i1.mass);
    }
    if(i1.genStatus==1&&i1.pdgID==2212) {
      TVector3 pOut_v3(i1.ps.x,i1.ps.y,i1.ps.z);
      pOut.SetVectM(pOut_v3,i1.mass);
    }
  }
  t_vec.push_back( -(pOut-pIn).Mag2() );
  
  return t_vec;
}

auto giveme_t_REC(const std::vector<ROOT::Math::PxPyPzMVector>& mom,
  const std::vector<dd4pod::Geant4ParticleData>& parts){
  
  std::vector<double> t_vec;
  TLorentzVector pIn,pOut;
  for(auto&i1 : parts){
    if(i1.genStatus==4&&i1.pdgID==2212) {
      TVector3 pIn_v3(i1.ps.x,i1.ps.y,i1.ps.z);
      pIn.SetVectM(pIn_v3,i1.mass);
    }
  }
  for(auto&i2: mom){
    pOut.SetPxPyPzE(i2.Px(),i2.Py(),i2.Pz(),i2.E());
  }
  t_vec.push_back( -(pOut-pIn).Mag2() );
  
  return t_vec;
}








