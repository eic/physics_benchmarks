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
double vm_mass_width[3]={0.04,0.02,0.03};
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
    // for (auto& i2: v2) {
      if (i1 != 0) {
        v.push_back((i1-v2[0])/i1);
      }
      else{
        v.push_back(-99);
      }
    // }
  }
  return v;
};

auto genMatch = [](std::vector<ROOT::Math::PxPyPzMVector> scatElec_REC,
                  std::vector<ROOT::Math::PxPyPzMVector> scatElec_MC,
                std::vector<ROOT::Math::PxPyPzMVector> vm_REC,
              std::vector<ROOT::Math::PxPyPzMVector> vm_MC )
{

  bool matchElectron=false;
  for(auto e1: scatElec_REC){
    for(auto e2: scatElec_MC){
      TLorentzVector e1_L;
       e1_L.SetPxPyPzE(e1.Px(),e1.Py(),e1.Pz(),e1.E());
      TLorentzVector e2_L;
       e2_L.SetPxPyPzE(e2.Px(),e2.Py(),e2.Pz(),e2.E());
      if( e1_L.DeltaR(e2_L)<1e-1 ) matchElectron=true ;
    }
  }

  bool matchVM=false;
  for(auto v1: vm_REC){
    for(auto v2: vm_MC){
      TLorentzVector v1_L;
       v1_L.SetPxPyPzE(v1.Px(),v1.Py(),v1.Pz(),v1.E());
      TLorentzVector v2_L;
       v2_L.SetPxPyPzE(v2.Px(),v2.Py(),v2.Pz(),v2.E());
      if( v1_L.DeltaR(v2_L)<1e-1 ) matchVM=true ;
    }
  }

  if( matchVM && matchElectron ) return 1;
  else return 0;

};

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

auto momenta_from_reconstruction_minus(const std::vector<eic::ReconstructedParticleData>& parts) {
  std::vector<ROOT::Math::PxPyPzMVector> momenta{parts.size()};
  std::transform(parts.begin(), parts.end(), momenta.begin(), [](const auto& part) {
     if(part.charge<0){
      return ROOT::Math::PxPyPzMVector{part.p.x, part.p.y, part.p.z, vm_daug_mass[which_vm]};
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
  //vector size mostly 4 sometimes 3 and 5, need to find out.  
  // std::cout << "size size = " << momenta.size() << std::endl;
  return momenta;
}

auto findScatElecMC(const std::vector<dd4pod::Geant4ParticleData>& parts)
{
  std::vector<ROOT::Math::PxPyPzMVector> momenta;
  for(auto& i1 : parts){
    if(i1.genStatus==21&&i1.pdgID==11) {
      momenta.push_back(ROOT::Math::PxPyPzMVector{i1.ps.x,i1.ps.y,i1.ps.z,i1.mass});
    }
    else {momenta.push_back(ROOT::Math::PxPyPzMVector{-1e10, -1e10, -1e10, -1e10});}
  }
  return momenta;
}

auto findVMMC(const std::vector<dd4pod::Geant4ParticleData>& parts) {
  std::vector<ROOT::Math::PxPyPzMVector> momenta;
  for(auto& i1 : parts){
    if(i1.genStatus==2&&i1.pdgID==vm_pid[which_vm]) {
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

auto getPt(const std::vector<ROOT::Math::PxPyPzMVector>& mom) {
  std::vector<double> ptVec(mom.size() );
  std::transform(mom.begin(), mom.end(), ptVec.begin(), [](const auto& part) {
    if(part.Px()<-1e9) return -99.;
    else return part.Pt();
  });
  return ptVec;
}

auto vector_sum = [](std::vector<ROOT::Math::PxPyPzMVector> p1, 
  std::vector<ROOT::Math::PxPyPzMVector> p2 ){
  std::vector<ROOT::Math::PxPyPzMVector> vm;
  for(auto& i1: p1){
    if(i1.Px()<-1e9) continue;
    for(auto& i2: p2){
      if(i2.Px()<-1e9) continue;
      //eta cut
      if(fabs(i1.Eta())>3.5||fabs(i2.Eta())>3.5) continue;
      vm.push_back(i1+i2);
    }
  }
  return vm;
};

//cut on phi mass region and rapidity phase space
auto getPt2OfVM(const std::vector<ROOT::Math::PxPyPzMVector>& mom) {
  std::vector<double> PtVec(mom.size() );
  std::transform(mom.begin(), mom.end(), PtVec.begin(), [](const auto& part) {
    if(fabs(part.M()-vm_mass[which_vm])>vm_mass_width[which_vm]||fabs(part.Rapidity())>3.5) return -99.;
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
  std::vector<double> etaVec;
  for(auto& i1:mom){
    double eta = i1.Eta();
    if(i1.Px()<-1e9){eta=-10.;}
    etaVec.push_back(eta);
  }
  return etaVec;
}

auto giveme_t = [](std::vector<ROOT::Math::PxPyPzMVector> vm, 
   std::vector<ROOT::Math::PxPyPzMVector> scatElec){
  std::vector<double > t_vec;
  if(vm.size()<1||scatElec.size()<1) {
    t_vec.push_back(-99.);
    return t_vec;
  }
  for(auto& i2: scatElec){
    for(auto& i1: vm){
      if(fabs(i1.Rapidity())>4.0||fabs(i1.M()-vm_mass[which_vm])>vm_mass_width[which_vm]) continue;
      if(i2.Px()<-1e9) continue;
      // TLorentzVector eIn(0,0,-18,18);
      // TLorentzVector eOut;eOut.SetPxPyPzE(i2.Px(),i2.Py(),i2.Pz(),i2.E());
      // TLorentzVector vmOut;vmOut.SetPxPyPzE(i1.Px(),i1.Py(),i1.Pz(),i1.E());
      // double method_E = (eIn-eOut-vmOut).Mag2();
      // t_vec.push_back( -method_E );
      TVector2 sum_pt(i1.Px()+i2.Px(), i1.Py()+i2.Py());
      t_vec.push_back( sum_pt.Mod2() );
    }
  }
  
  return t_vec;
};
