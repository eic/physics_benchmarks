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
#include "TLorentzRotation.h"
#include "TVector2.h"

#include "fmt/color.h"
#include "fmt/core.h"

#include "nlohmann/json.hpp"
#include "eicd/InclusiveKinematicsData.h"
#include "eicd/ReconstructedParticleData.h"

#define PI            3.1415926
#define MASS_ELECTRON 0.00051
#define MASS_PROTON   0.93827
#define MASS_PION     0.13957
#define MASS_KAON     0.493667
#define MASS_AU197    183.45406466643374

int nEvent=1;
int which_vm = 1;
double vm_pid[3]={113,333,443};
double vm_mass[3]={0.77545,1.019,3.0969};
double vm_mass_width[3]={0.15,0.02,0.03};
double vm_daug_pid[3]={211,321,11};
double vm_daug_mass[3]={MASS_PION,MASS_KAON,MASS_ELECTRON};

//0 beagle, 1 sartre, 2 estarlight
int which_mc = 0;
int genStatus_scatElec[3]={21,1,1};
int genStatus_VM[3]={2,2,2};
//estarlight does not have status 2 phi, so dummy code for genStatus_VM[2].

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
    if(part.charge>0){
      return ROOT::Math::PxPyPzMVector{part.p.x, part.p.y, part.p.z, vm_daug_mass[which_vm]};
    }
    else{
      return ROOT::Math::PxPyPzMVector{-1e10, -1e10, -1e10, -1e10};
    }
  });
  return momenta;
}

auto momenta_from_reconstruction_minus(const std::vector<eic::ReconstructedParticleData>& parts,
  std::vector<int> scat_id,
    std::vector<int> scat_source) 

{
  std::vector<ROOT::Math::PxPyPzMVector> momenta;
  for(auto& i1 : parts){
   if(i1.charge<0){
      momenta.push_back(ROOT::Math::PxPyPzMVector{i1.p.x, i1.p.y, i1.p.z, vm_daug_mass[which_vm]});
    }
    else{
      momenta.push_back(ROOT::Math::PxPyPzMVector{-1e10, -1e10, -1e10, -1e10});
    }
  }
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

// auto findScatElec(const std::vector<eic::ReconstructedParticleData>& parts, 
//       std::vector<int> scat_id,
//     std::vector<int> scat_source) 
// {
//   std::vector<ROOT::Math::PxPyPzMVector> momenta;
//   for(auto& i1 : parts){
//     if(scat_id.size()>0 
//         && scat_source.size()>0
//           &&i1.ID.value==scat_id[0]
//             &&i1.ID.source==scat_source[0])
//     {
//       auto scat = ROOT::Math::PxPyPzMVector{i1.p.x, i1.p.y, i1.p.z, MASS_ELECTRON};
//       momenta.push_back(scat);
//     }
//     else{
//       momenta.push_back(ROOT::Math::PxPyPzMVector{-1e10, -1e10, -1e10, -1e10});
//     }
  
//   }
//   return momenta;
// }

auto findScatElec(const std::vector<eic::ReconstructedParticleData>& parts, 
      std::vector<int> scat_id,
    std::vector<int> scat_source) 
{
  std::vector<ROOT::Math::PxPyPzMVector> momenta;
  for(auto& i1 : parts){
    if(i1.ID.value==0)
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
    if(i1.genStatus==genStatus_scatElec[which_mc]&&i1.pdgID==11) {
      momenta.push_back(ROOT::Math::PxPyPzMVector{i1.ps.x,i1.ps.y,i1.ps.z,i1.mass});
    }
    else {momenta.push_back(ROOT::Math::PxPyPzMVector{-1e10, -1e10, -1e10, -1e10});}
  }
  return momenta;
}

auto findVMMC(const std::vector<dd4pod::Geant4ParticleData>& parts) {
  std::vector<ROOT::Math::PxPyPzMVector> momenta;
  for(auto& i1 : parts){
    if(i1.genStatus==genStatus_VM[which_mc]&&i1.pdgID==vm_pid[which_vm]) {
      momenta.push_back(ROOT::Math::PxPyPzMVector{i1.ps.x,i1.ps.y,i1.ps.z,i1.mass});
    }
    else {momenta.push_back(ROOT::Math::PxPyPzMVector{-1e10, -1e10, -1e10, -1e10});}
  }
  return momenta;
}

auto findVM_DaugPlus_MC(const std::vector<dd4pod::Geant4ParticleData>& parts) {
  std::vector<ROOT::Math::PxPyPzMVector> momenta;
  for(auto& i1 : parts){
    if(i1.genStatus==1&&i1.pdgID==vm_daug_pid[which_vm]) {
      momenta.push_back(ROOT::Math::PxPyPzMVector{i1.ps.x,i1.ps.y,i1.ps.z,i1.mass});
    }
    else {momenta.push_back(ROOT::Math::PxPyPzMVector{-1e10, -1e10, -1e10, -1e10});}
  }
  return momenta;
}

auto findVM_DaugMinus_MC(const std::vector<dd4pod::Geant4ParticleData>& parts) {
  std::vector<ROOT::Math::PxPyPzMVector> momenta;
  for(auto& i1 : parts){
    if(i1.genStatus==1&&i1.pdgID==-vm_daug_pid[which_vm]) {
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
      momenta.push_back(ROOT::Math::PxPyPzMVector{i1.ps.x, i1.ps.y, i1.ps.z, i1.mass});
    }
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

auto vector_sum = [](std::vector<ROOT::Math::PxPyPzMVector> p1, 
  std::vector<ROOT::Math::PxPyPzMVector> p2 ){
  std::vector<ROOT::Math::PxPyPzMVector> vm;
  for(auto& i1: p1){
    if(i1.Px()<-1e9) continue;
    for(auto& i2: p2){
      if(i2.Px()<-1e9) continue;
      //no detector cuts.
      vm.push_back(i1+i2);
    }
  }
  return vm;
};

auto findVM_MC_match_REC(const std::vector<ROOT::Math::PxPyPzMVector> MC, 
  const std::vector<ROOT::Math::PxPyPzMVector> REC)
{
  std::vector<ROOT::Math::PxPyPzMVector> vm_match;
  for(auto& i1:MC){
    auto v = ROOT::Math::PxPyPzMVector{-1e10, -1e10, -1e10, -1e10};
    if(i1.Px()<-1e9) continue;
    for(auto& i2:REC){
      if(matchVectKine(i1,i2)&&fabs(i2.M()-i1.M())<vm_mass_width[which_vm]) v=i1;
    }
    vm_match.push_back(v);
  }
  return vm_match;
}

auto findVM_REC_NOT_match_MC(const std::vector<ROOT::Math::PxPyPzMVector> REC, 
  const std::vector<ROOT::Math::PxPyPzMVector> MC)
{
  std::vector<ROOT::Math::PxPyPzMVector> vm_not_match;
  for(auto& i1:REC){
    auto v = ROOT::Math::PxPyPzMVector{-1e10, -1e10, -1e10, -1e10};
    //here not cutting on rec mass yet.
    bool VM_interested = true; 
    if(fabs(i1.M()-vm_mass[which_vm])>vm_mass_width[which_vm]) VM_interested=false;
    for(auto& i2:MC){
      if(i2.Px()<-1e-9) continue;
      if(!matchVectKine(i1,i2)&&VM_interested) v=i1;
    }
    vm_not_match.push_back(v);
  }
  return vm_not_match;
}

auto resolution_MC_match_REC(const std::vector<ROOT::Math::PxPyPzMVector> MC, 
  const std::vector<ROOT::Math::PxPyPzMVector> REC)
{
  std::vector<double > resolution;
  for(auto& i1:MC){
    double res = -99;
    for(auto& i2:REC){
      if(i1.Px()<-1e9||i2.Px()<-1e9) continue;
        if(matchVectKine(i1,i2)&&fabs(i2.M()-i1.M())<vm_mass_width[which_vm]){
          if(fabs(i2.Eta())<1.0){
            res = (i1.Pt()-i2.Pt())/i1.Pt();
          }
        } 
    }
    resolution.push_back(res);
  }
  return resolution;
}

auto resolution_MC_match_REC_electron(const std::vector<ROOT::Math::PxPyPzMVector> MC, 
  const std::vector<ROOT::Math::PxPyPzMVector> REC)
{
  std::vector<double > resolution;
  for(auto& i1:MC){
    double res = -99;
    for(auto& i2:REC){
      if(i1.Px()<-1e9||i2.Px()<-1e9) continue;
        if(matchVectKine(i1,i2)&&fabs(i2.M()-i1.M())<vm_mass_width[which_vm]){
            res = (i1.Pt()-i2.Pt())/i1.Pt();
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

auto getP(const std::vector<ROOT::Math::PxPyPzMVector>& mom)
{
  std::vector<double> pVec;
  for(auto& i1:mom){
    double momentum = i1.P();
    if(i1.Px()<-1e9){momentum=-10.;}
    pVec.push_back(momentum);
  }
  return pVec;

}

auto getPtVM(const std::vector<ROOT::Math::PxPyPzMVector>& mom) {
  std::vector<double> ptVec(mom.size() );
  std::transform(mom.begin(), mom.end(), ptVec.begin(), [](const auto& part) {
    if(part.Px()<-1e9||fabs(part.M()-vm_mass[which_vm])>vm_mass_width[which_vm]) return -99.;
    else return part.Pt();
  });
  return ptVec;
}

auto getPtVM_match = [](std::vector<double> pt_MC, std::vector<double> pt_REC) {
  std::vector<double> ptVec;
  bool hasREC_=false;
  for(auto& i2 : pt_REC){
    if(i2>0.) hasREC_=true;
  }
  for(auto& i1 : pt_MC){
    if(i1!=-99.&&hasREC_) ptVec.push_back(i1);
    else ptVec.push_back(-99);
  }
  return ptVec;
};

auto getEta(const std::vector<ROOT::Math::PxPyPzMVector>& mom) {
  std::vector<double> etaVec;
  for(auto& i1:mom){
    double eta = i1.Eta();
    if(i1.Px()<-1e9){eta=-10.;}
    etaVec.push_back(eta);
  }
  return etaVec;
}

auto getEtaVM(const std::vector<ROOT::Math::PxPyPzMVector>& mom) {
  std::vector<double> etaVec;
  for(auto& i1:mom){
    double eta = i1.Eta();
    if(i1.Px()<-1e9||fabs(i1.M()-vm_mass[which_vm])>vm_mass_width[which_vm]){eta=-10.;}
    etaVec.push_back(eta);
  }
  return etaVec;
}

auto getNtrk(const std::vector<eic::ReconstructedParticleData>& parts) 
{
  std::vector<int> mult;
  int n=0;
  for(auto& i1 : parts){
    if(i1.charge!=0) n++;
  }
  mult.push_back( n );
  return mult;
}

auto giveme_t_E = [](std::vector<ROOT::Math::PxPyPzMVector> vm, 
   std::vector<ROOT::Math::PxPyPzMVector> scatElec,
  const std::vector<dd4pod::Geant4ParticleData>& mc){

  TLorentzVector eIn(0,0,-18,18);
  std::vector<double > t_vec;
  for(auto& i2: scatElec){
    for(auto& i1: vm){
      if(fabs(i1.Rapidity())>3.0||fabs(i1.M()-vm_mass[which_vm])>vm_mass_width[which_vm]) continue;
      if(i2.Px()<-1e9) continue;
      TLorentzVector eOut;eOut.SetPxPyPzE(i2.Px(),i2.Py(),i2.Pz(),i2.E());
      TLorentzVector vmOut;vmOut.SetPxPyPzE(i1.Px(),i1.Py(),i1.Pz(),i1.E());
      //exact method need no boost.
      double method_E = (eIn-eOut-vmOut).Mag2();
      t_vec.push_back( -method_E );
    }
  }
  return t_vec;
};

auto giveme_t_MC_E = [](std::vector<ROOT::Math::PxPyPzMVector> vm, 
   std::vector<ROOT::Math::PxPyPzMVector> scatElec,
  const std::vector<dd4pod::Geant4ParticleData>& mc){

  TLorentzVector eIn(0,0,-18,18);
  TLorentzVector eOut(0,0,0,0.0);
  TLorentzVector vmOut(0,0,0,0.0);
  for(auto& i3 : mc){
    if(i3.genStatus==4&&i3.pdgID==11){
      TVector3 eInv3(i3.ps.x,i3.ps.y,i3.ps.z);
      eIn.SetVectM(eInv3,i3.mass);
    } 
    if(i3.genStatus==1&&i3.pdgID==11){
      TVector3 eOutv3(i3.ps.x,i3.ps.y,i3.ps.z);
      eOut.SetVectM(eOutv3,i3.mass);
    }
    if(i3.genStatus==2&&i3.pdgID==vm_pid[which_vm]){
      TVector3 vmOutv3(i3.ps.x,i3.ps.y,i3.ps.z);
      vmOut.SetVectM(vmOutv3,i3.mass);
    }
  }
  std::vector<double > t_vec;
  double method_E = (eIn-eOut-vmOut).Mag2();
  t_vec.push_back( -method_E );
  
  return t_vec;
};

auto giveme_t_A = [](std::vector<ROOT::Math::PxPyPzMVector> vm, 
   std::vector<ROOT::Math::PxPyPzMVector> scatElec,
  const std::vector<dd4pod::Geant4ParticleData>& mc){

  TLorentzVector vmOut_MC, eOut_MC;
  for(auto& i3 : mc){
    if(i3.genStatus==1&&i3.pdgID==11){
      TVector3 eOutv3(i3.ps.x,i3.ps.y,i3.ps.z);
      eOut_MC.SetVectM(eOutv3,i3.mass);
    }
    if(i3.genStatus==2&&i3.pdgID==vm_pid[which_vm]){
      TVector3 vmOutv3(i3.ps.x,i3.ps.y,i3.ps.z);
      vmOut_MC.SetVectM(vmOutv3,i3.mass);
    }
  }
  std::vector<double > t_vec;
  for(auto& i2: scatElec){
    for(auto& i1: vm){
      if(i1.Px()<-1e9) continue;
      if(fabs(i1.Rapidity())>4.0||fabs(i1.M()-vm_mass[which_vm])>vm_mass_width[which_vm]) continue;
      if(i2.Px()<-1e9) continue;

      TLorentzVector eOut;eOut.SetPxPyPzE(i2.Px(),i2.Py(),i2.Pz(),i2.E());
      TLorentzVector vmOut; vmOut.SetPxPyPzE(i1.Px(),i1.Py(),i1.Pz(),i1.E());
      
      // vmOut = vmOut_MC;
      // eOut = eOut_MC;
      
      TVector2 sum_pt(eOut.Px()+vmOut.Px(), eOut.Py()+vmOut.Py());
      t_vec.push_back( sum_pt.Mod2() );
    }
  }

  return t_vec;
};

auto giveme_t_L = [](std::vector<ROOT::Math::PxPyPzMVector> vm, 
   std::vector<ROOT::Math::PxPyPzMVector> scatElec,
  const std::vector<dd4pod::Geant4ParticleData>& mc){

  TLorentzVector eIn(0,0,-18,18);
  TLorentzVector eInTrue(0,0,-18,18);
  TLorentzVector pInTrue(0,0,109.996,110.000);
  TLorentzVector pIn(0,0,109.996,110.000);
  for(auto& i3 : mc){
    if(i3.genStatus==4&&i3.pdgID==11){
      TVector3 eInv3(i3.ps.x,i3.ps.y,i3.ps.z);
      eIn.SetVectM(eInv3,MASS_ELECTRON);
    } 
    if(i3.genStatus==4&&i3.pdgID==2212){
      TVector3 pInv3(i3.ps.x,i3.ps.y,i3.ps.z);
      pIn.SetVectM(pInv3,MASS_PROTON);
    }
  }
  std::vector<double > t_vec;
  for(auto& i2: scatElec){
    for(auto& i1: vm){
      if(fabs(i1.Rapidity())>3.0||fabs(i1.M()-vm_mass[which_vm])>vm_mass_width[which_vm]) continue;
      if(i2.Px()<-1e9) continue;
      TLorentzVector eOut;eOut.SetPxPyPzE(i2.Px(),i2.Py(),i2.Pz(),i2.E());
      TLorentzVector vmOut;vmOut.SetPxPyPzE(i1.Px(),i1.Py(),i1.Pz(),i1.E());
      TLorentzVector aInVec(pIn.Px()*197,pIn.Py()*197,pIn.Pz()*197,sqrt(pIn.Px()*197*pIn.Px()*197 + pIn.Py()*197*pIn.Py()*197 + pIn.Pz()*197*pIn.Pz()*197 + MASS_AU197*MASS_AU197) );

      double method_L = -99.;
      TLorentzVector a_beam_scattered = aInVec-(vmOut+eOut-eIn);
      double p_Aplus = a_beam_scattered.E()+a_beam_scattered.Pz();
      double p_TAsquared = TMath::Power(a_beam_scattered.Pt(),2);
      double p_Aminus = (MASS_AU197*MASS_AU197 + p_TAsquared) / p_Aplus;
      TLorentzVector a_beam_scattered_corr; 
      a_beam_scattered_corr.SetPxPyPzE(a_beam_scattered.Px(),a_beam_scattered.Py(),(p_Aplus-p_Aminus)/2., (p_Aplus+p_Aminus)/2. );
      method_L = (a_beam_scattered_corr-aInVec).Mag2();

      t_vec.push_back( -method_L );
    }
  }
  return t_vec;
};
