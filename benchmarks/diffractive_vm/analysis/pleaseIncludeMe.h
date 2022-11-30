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
#include "edm4eic/InclusiveKinematicsData.h"
#include "edm4eic/ReconstructedParticleData.h"
#include "edm4eic/ClusterData.h"
#include "edm4eic/MCRecoParticleAssociationData.h"
#include "edm4eic/MCRecoClusterParticleAssociationData.h"
#include "edm4hep/MCParticleData.h"

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
double vm_daug_pid[3]={211,321,-11};
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

auto momenta_from_reconstruction_plus(const std::vector<edm4eic::ReconstructedParticleData>& parts) {
  std::vector<ROOT::Math::PxPyPzMVector> momenta{parts.size()};
  std::transform(parts.begin(), parts.end(), momenta.begin(), [](const auto& part) {
    if(part.charge>0){
      return ROOT::Math::PxPyPzMVector{part.momentum.x, part.momentum.y, part.momentum.z, vm_daug_mass[which_vm]};
    }
    else{
      return ROOT::Math::PxPyPzMVector{-1e10, -1e10, -1e10, -1e10};
    }
  });
  return momenta;
}

auto momenta_from_reconstruction_minus(const std::vector<edm4eic::ReconstructedParticleData>& parts) 
{
  std::vector<ROOT::Math::PxPyPzMVector> momenta;
  for(auto& i1 : parts){
   if(i1.charge<0){
      momenta.push_back(ROOT::Math::PxPyPzMVector{i1.momentum.x, i1.momentum.y, i1.momentum.z, vm_daug_mass[which_vm]});
    }
    else{
      momenta.push_back(ROOT::Math::PxPyPzMVector{-1e10, -1e10, -1e10, -1e10});
    }
  }
  return momenta;
}

auto findScatElec(const std::vector<edm4eic::ReconstructedParticleData>& recs, 
                    const std::vector<edm4hep::MCParticleData>& mcs) 
{
  std::vector<ROOT::Math::PxPyPzMVector> momenta;
  //finding mc scat e'
  TVector3 trkMC(0,0,0);
  for(auto& i2 : mcs){
    if(i2.charge<0 && 
        i2.generatorStatus==genStatus_scatElec[which_mc]
          &&i2.PDG==11){ trkMC.SetXYZ(i2.momentum.x,i2.momentum.y,i2.momentum.z); } 
  }
  //kinematic match 
  //need to change to association and cluster matching.
  double minR=99;
  TVector3 matchRECTrk(-1e10,-1e10,-1e10);
  for(auto& i1 : recs){
    TVector3 trkREC(i1.momentum.x,i1.momentum.y,i1.momentum.z);
    if(i1.charge<0 )
    {
      if(trkREC.DeltaR(trkMC)<minR){
        minR=trkREC.DeltaR(trkMC);
        matchRECTrk=trkREC;
      }
    }
  }
  auto scat = ROOT::Math::PxPyPzMVector{matchRECTrk.Px(), matchRECTrk.Py(), matchRECTrk.Pz(), MASS_ELECTRON};
  momenta.push_back(scat);

  return momenta;
}
auto findScatElecTest(const std::vector<edm4eic::ReconstructedParticleData>& parts,
                        const std::vector<edm4eic::ClusterData>& clusters,
                            const std::vector<edm4eic::MCRecoParticleAssociationData>& assocs,
                              const std::vector<edm4eic::MCRecoClusterParticleAssociationData>& cluster_assocs) 
{
  /*
  Comment - Things to think about:
    - We need to have some matching/projection
    code to match between tracks and clusters.
    Now, everything below is just for now. 
  */
  std::vector<ROOT::Math::PxPyPzMVector> momenta;
  TLorentzVector escat(-1E10, -1E10, -1E10, -1E10);
  //EEMC
  double maxEnergy=0;
  int index=-1;
  int cluster_rec_leading_index=-1;
  for(auto& i1 : clusters){
    index++;
    auto energy=i1.energy;
    if(energy>maxEnergy){
      maxEnergy=energy;
      cluster_rec_leading_index=index;
    }
  }
  //Find sim id in cluster
  int cluster_sim_leading=-1;
  for(auto& i2 : cluster_assocs){
    int rec_clus_id=i2.recID;
    int sim_clus_id=i2.simID;

    if(rec_clus_id==cluster_rec_leading_index){
      cluster_sim_leading=sim_clus_id;
    }
  }

  //rec finding leading momentum as scat' e
  double maxMom=0.;
  TVector3 maxtrk(-1E10,-1E10,-1E10);
  int elec_index=-1;
  index=-1;
  for(auto& i2 : parts){
    index++;
    TVector3 trk(i2.momentum.x,i2.momentum.y,i2.momentum.z);
    if(i2.charge>0) continue;
    if(trk.Mag()>maxMom){
      maxMom=trk.Mag();
      maxtrk=trk;
      elec_index=index;
    }
  }

  //finding track assoc.
  int mc_elect_index=-1;
  for(auto& i3 : assocs){
    int rec_id = i3.recID;
    int sim_id = i3.simID;
    if (rec_id == elec_index) mc_elect_index=sim_id;
  }
  
  //3-second calibration.
  // maxEnergy+=0.9;
  //electron hypothesis;
  double p = sqrt(maxEnergy*maxEnergy- MASS_ELECTRON*MASS_ELECTRON );
  double eta=maxtrk.Eta();
  double phi=maxtrk.Phi();
  double pt = TMath::Sin(maxtrk.Theta())*p;
  escat.SetPtEtaPhiM(pt,eta,phi,MASS_ELECTRON);
  
  if( cluster_sim_leading == mc_elect_index && mc_elect_index != -1 ) {
    momenta.push_back(ROOT::Math::PxPyPzMVector{escat.Px(),escat.Py(),escat.Pz(),MASS_ELECTRON});
  }
  return momenta;
}

auto findScatElecMC(const std::vector<edm4hep::MCParticleData>& parts)
{
  std::vector<ROOT::Math::PxPyPzMVector> momenta;
  for(auto& i1 : parts){
    if(i1.generatorStatus==genStatus_scatElec[which_mc]&&i1.PDG==11) {
      momenta.push_back(ROOT::Math::PxPyPzMVector{i1.momentum.x,i1.momentum.y,i1.momentum.z,i1.mass});
    }
    else {momenta.push_back(ROOT::Math::PxPyPzMVector{-1e10, -1e10, -1e10, -1e10});}
  }
  return momenta;
}

auto findVMMC(const std::vector<edm4hep::MCParticleData>& parts) {
  std::vector<ROOT::Math::PxPyPzMVector> momenta;
  for(auto& i1 : parts){
    if(i1.generatorStatus==genStatus_VM[which_mc]&&i1.PDG==vm_pid[which_vm]) {
      momenta.push_back(ROOT::Math::PxPyPzMVector{i1.momentum.x,i1.momentum.y,i1.momentum.z,i1.mass});
    }
    else {momenta.push_back(ROOT::Math::PxPyPzMVector{-1e10, -1e10, -1e10, -1e10});}
  }
  return momenta;
}

auto findVM_DaugPlus_MC(const std::vector<edm4hep::MCParticleData>& parts) {
  std::vector<ROOT::Math::PxPyPzMVector> momenta;
  for(auto& i1 : parts){
    if(i1.generatorStatus==1&&i1.PDG==vm_daug_pid[which_vm]) {
      momenta.push_back(ROOT::Math::PxPyPzMVector{i1.momentum.x,i1.momentum.y,i1.momentum.z,i1.mass});
    }
    else {momenta.push_back(ROOT::Math::PxPyPzMVector{-1e10, -1e10, -1e10, -1e10});}
  }
  return momenta;
}

auto findVM_DaugMinus_MC(const std::vector<edm4hep::MCParticleData>& parts) {
  std::vector<ROOT::Math::PxPyPzMVector> momenta;
  for(auto& i1 : parts){
    if(i1.generatorStatus==1&&i1.PDG==-vm_daug_pid[which_vm]) {
      momenta.push_back(ROOT::Math::PxPyPzMVector{i1.momentum.x,i1.momentum.y,i1.momentum.z,i1.mass});
    }
    else {momenta.push_back(ROOT::Math::PxPyPzMVector{-1e10, -1e10, -1e10, -1e10});}
  }
  return momenta;
}

auto findScatProtonMC(const std::vector<edm4hep::MCParticleData>& parts){
  std::vector<ROOT::Math::PxPyPzMVector> momenta;
  for(auto& i1 : parts){
    if(i1.generatorStatus==1&&i1.PDG==2212){
      momenta.push_back(ROOT::Math::PxPyPzMVector{i1.momentum.x, i1.momentum.y, i1.momentum.z, i1.mass});
    }
  }
  return momenta;
}

auto findScatProton(const std::vector<edm4eic::ReconstructedParticleData>& FF){
  std::vector<ROOT::Math::PxPyPzMVector> momenta;
  for(auto& i1 : FF){
    if(i1.charge==1){
      momenta.push_back(ROOT::Math::PxPyPzMVector{i1.momentum.x,i1.momentum.y,i1.momentum.z,i1.mass});
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
            res = (i1.Pt()-i2.Pt())/i1.Pt();
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
        if(matchVectKine(i1,i2)
          &&fabs(i2.M()-i1.M())<vm_mass_width[which_vm]){
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

auto getRapVM(const std::vector<ROOT::Math::PxPyPzMVector>& mom) {
  std::vector<double> etaVec;
  for(auto& i1:mom){
    double eta = i1.Rapidity();
    if(i1.Px()<-1e9||fabs(i1.M()-vm_mass[which_vm])>vm_mass_width[which_vm]){eta=-10.;}
    etaVec.push_back(eta);
  }
  return etaVec;
}

auto getNtrk(const std::vector<edm4eic::ReconstructedParticleData>& parts) 
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
  const std::vector<edm4hep::MCParticleData>& mc){

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
  const std::vector<edm4hep::MCParticleData>& mc){

  TLorentzVector photIn(0.,0.,0.,0.);
  TLorentzVector vmOut(0.,0.,0.,0.);
  for(auto& i3 : mc){
    if(i3.generatorStatus==21&&i3.PDG==22){
      TVector3 photInv3(i3.momentum.x,i3.momentum.y,i3.momentum.z);
      photIn.SetVectM(photInv3,i3.mass);
    }
    if(i3.generatorStatus==2&&i3.PDG==vm_pid[which_vm]){
      TVector3 vmOutv3(i3.momentum.x,i3.momentum.y,i3.momentum.z);
      vmOut.SetVectM(vmOutv3,i3.mass);
    }
  }
  std::vector<double > t_vec;
  double method_E = (photIn-vmOut).Mag2();
  t_vec.push_back( -method_E );
  
  return t_vec;
};

//scatElec will have ambuity with J/psi decay to ee.
auto giveme_t_A = [](std::vector<ROOT::Math::PxPyPzMVector> vm, 
   std::vector<ROOT::Math::PxPyPzMVector> scatElec,
  const std::vector<edm4hep::MCParticleData>& mc){

  TLorentzVector vmOut_MC, eOut_MC;
  for(auto& i3 : mc){
    if(i3.generatorStatus==1&&i3.PDG==11){
      TVector3 eOutv3(i3.momentum.x,i3.momentum.y,i3.momentum.z);
      eOut_MC.SetVectM(eOutv3,i3.mass);
    }
    if(i3.generatorStatus==2&&i3.PDG==vm_pid[which_vm]){
      TVector3 vmOutv3(i3.momentum.x,i3.momentum.y,i3.momentum.z);
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
      
      vmOut = vmOut_MC;
      eOut = eOut_MC;

      double e_pt_res=gRandom->Gaus(0.0,0.014);//gaus fit~0.0068, set this number by full simulations for Q2>15
      double e_pt = eOut.Pt()*(1.+e_pt_res);
      eOut.SetPerp(e_pt);

      double vm_pt_res=gRandom->Gaus(0.0,0.011);//gaus fit~0.0065, set this number by full simulations for Q2>15
      double vm_pt = vmOut.Pt()*(1.+vm_pt_res);
      vmOut.SetPerp(vm_pt);
      
      TVector2 sum_pt(eOut.Px()+vmOut.Px(), eOut.Py()+vmOut.Py());
      t_vec.push_back( sum_pt.Mod2() );
    }
  }

  return t_vec;
};

auto giveme_t_L = [](std::vector<ROOT::Math::PxPyPzMVector> vm, 
   std::vector<ROOT::Math::PxPyPzMVector> scatElec,
  const std::vector<edm4hep::MCParticleData>& mc){

  TLorentzVector eIn(0,0,-18,18);
  // TLorentzVector eIn(0,0,-5,5);
  TLorentzVector pIn(-2.749,0,109.996,110.034);
  // TLorentzVector pIn(0,0,109.996,110.00);//no Crossing angle
  TLorentzVector vmOut_MC, eOut_MC;
  for(auto& i3 : mc){
     if(i3.generatorStatus==1&&i3.PDG==11){
      TVector3 eOutv3(i3.momentum.x,i3.momentum.y,i3.momentum.z);
      eOut_MC.SetVectM(eOutv3,i3.mass);
    }
    if(i3.generatorStatus==2&&i3.PDG==vm_pid[which_vm]){
      TVector3 vmOutv3(i3.momentum.x,i3.momentum.y,i3.momentum.z);
      vmOut_MC.SetVectM(vmOutv3,i3.mass);
    }
  }
  std::vector<double > t_vec;
  for(auto& i2: scatElec){
    for(auto& i1: vm){
      if(fabs(i1.Rapidity())>4.0||fabs(i1.M()-vm_mass[which_vm])>vm_mass_width[which_vm]) continue;
      if(i2.Px()<-1e9) continue;
      TLorentzVector eOut;eOut.SetPxPyPzE(i2.Px(),i2.Py(),i2.Pz(),i2.E());
      TLorentzVector vmOut;vmOut.SetPxPyPzE(i1.Px(),i1.Py(),i1.Pz(),i1.E());
      TLorentzVector aInVec(pIn.Px()*197,pIn.Py()*197,pIn.Pz()*197,sqrt(pIn.Px()*197*pIn.Px()*197 + pIn.Py()*197*pIn.Py()*197 + pIn.Pz()*197*pIn.Pz()*197 + MASS_AU197*MASS_AU197) );
      
      // vmOut = vmOut_MC;
      // eOut = eOut_MC;

      // double e_pt_res=gRandom->Gaus(0.0,0.012);//gaus fit~0.0068, set this number by full simulations for Q2>15
      // double e_pt = eOut.Pt()*(1.+e_pt_res);
      // eOut.SetPerp(e_pt);

      // double vm_pt_res=gRandom->Gaus(0.0,0.015);//gaus fit~0.0065, set this number by full simulations for Q2>15
      // double vm_pt = vmOut.Pt()*(1.+vm_pt_res);
      // vmOut.SetPerp(vm_pt);

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