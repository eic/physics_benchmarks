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
  if(fabs(v1_L.M()-v2_L.M())>1e-3) return false;
  if(v1_L.DeltaR(v2_L)>0.3) return false;
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

auto findScatProtonMC(const std::vector<dd4pod::Geant4ParticleData>& parts){
  std::vector<ROOT::Math::PxPyPzMVector> momenta;
  for(auto& i1 : parts){
    if(i1.genStatus==1&&i1.pdgID==2212){
      //Double_t px = i1.ps.x;
      //Double_t py = i1.ps.y;
      //Double_t pz = i1.ps.z;
      //Double_t s = std::sin(0.025);
      //Double_t c = std::cos(0.025);
      //Double_t zz = pz;
      //Double_t xx = px;
      //pz = c*zz - s*xx;
      //px = s*zz + c*xx;
      //momenta.push_back(ROOT::Math::PxPyPzMVector{px, py, pz, i1.mass});
      momenta.push_back(ROOT::Math::PxPyPzMVector{i1.ps.x, i1.ps.y, i1.ps.z, i1.mass});
    }
  }
  return momenta;
}

auto findScatProton(const std::vector<eic::ReconstructedParticleData>& FF){
  std::vector<ROOT::Math::PxPyPzMVector> momenta;
  for(auto& i1 : FF){
    if(i1.charge==1){
      //Double_t px = i1.p.x;
      //Double_t py = i1.p.y;
      //Double_t pz = i1.p.z;
      //Double_t s = std::sin(0.025);
      //Double_t c = std::cos(0.025);
      //Double_t zz = pz;
      //Double_t xx = px;
      //pz = c*zz - s*xx;
      //px = s*zz + c*xx;
      //momenta.push_back(ROOT::Math::PxPyPzMVector{px,py,pz,i1.mass});
      momenta.push_back(ROOT::Math::PxPyPzMVector{i1.p.x,i1.p.y,i1.p.z,i1.mass});
    }
  }
  return momenta;
}

auto findPhot_MC_match_REC(const std::vector<ROOT::Math::PxPyPzMVector> MC, 
  const std::vector<ROOT::Math::PxPyPzMVector> REC)
{
  std::vector<ROOT::Math::PxPyPzMVector> ph_match;
  for(auto& i1:MC){
    auto v = ROOT::Math::PxPyPzMVector{-1e10, -1e10, -1e10, -1e10};
    if(i1.Px()<-1e9) continue;
    for(auto& i2:REC){
      if(matchVectKine(i1,i2)) v=i1;
    }
    ph_match.push_back(v);
  }
  return ph_match;
}

auto findPhot_REC_not_match_MC(const std::vector<ROOT::Math::PxPyPzMVector> REC, 
  const std::vector<ROOT::Math::PxPyPzMVector> MC)
{
  std::vector<ROOT::Math::PxPyPzMVector> ph_not_match;
  
  for(auto& i1:REC){
    auto v = ROOT::Math::PxPyPzMVector{-1e10, -1e10, -1e10, -1e10};
    if(i1.Px()<-1e9) continue;
    for(auto& i2:MC){
      if(i2.Px()<-1e9) continue;
      if(!matchVectKine(i1,i2)) v=i1;
    }
    ph_not_match.push_back(v);
  }
  return ph_not_match;
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

auto getTheta(const std::vector<ROOT::Math::PxPyPzMVector>& mom) {
  std::vector<double> thetaVec;
  for(auto& i1:mom){
    double theta = i1.Theta();
    if(i1.Px()<-1e9){theta=-10.;}
    thetaVec.push_back(theta);
  }
  return thetaVec;
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

auto getE(const std::vector<ROOT::Math::PxPyPzMVector>& mom) {
  std::vector<double> energyVec;
  for(auto& i1:mom){
    double energy = i1.E();
    if(i1.Px()<-1e9){energy=-10.;}
    energyVec.push_back(energy);
  }
  return energyVec;
}

auto getAngleDiff(const std::vector<ROOT::Math::PxPyPzMVector> ph_gen,
  const std::vector<ROOT::Math::PxPyPzMVector> ph_rec)
{
  std::vector<double> angleVec;
  for(auto& i1:ph_gen){
    if(i1.Px()<-1e9) continue;
    TLorentzVector ph_gen_4v(i1.Px(),i1.Py(),i1.Pz(),i1.E());
    for(auto& i2:ph_rec){
      if(i2.Px()<-1e9) continue;
      TLorentzVector ph_rec_4v(i2.Px(),i2.Py(),i2.Pz(),i2.E());
      double angle = ph_gen_4v.Angle(ph_rec_4v.Vect());
      angleVec.push_back(angle);
    }
  }
  return angleVec;
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
    TVector3 pOut_v3(i2.Px(),i2.Py(),i2.Pz());
    pOut.SetVectM(pOut_v3,0.93827);
  }
  t_vec.push_back( -(pOut-pIn).Mag2() );
  
  return t_vec;
}

auto giveme_t = [](std::vector<ROOT::Math::PxPyPzMVector> vm, 
   std::vector<ROOT::Math::PxPyPzMVector> scatElec){
  std::vector<double > t_vec;
  for(auto& i2: scatElec){
    for(auto& i1: vm){
      if(fabs(i1.Rapidity())>3.0) continue;
      if(i2.Px()<-1e9) continue;
      TVector2 sum_pt(i1.Px()+i2.Px(), i1.Py()+i2.Py());
      t_vec.push_back( sum_pt.Mod2() );
    }
  }
  
  return t_vec;
};






