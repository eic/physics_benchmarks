#pragma once
#include "ReactionKinematics.h"


namespace rad{
  namespace electro{
  

    // template<typename Tp, typename Tm>
    // PxPyPzMVector PhotonVector(const config::RVecIndexMap& react,const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m){
    //   auto phot =  beams::InitialFourVector(react[names::ElectroEleIdx()][0],px,py,pz,m);
    //   SubtractFourVector(phot,react[names::ScatEleIdx()],px,py,pz,m);
    //   return phot;
    //  }

    template<typename Tp, typename Tm>
    double Q2(const config::RVecIndexMap& react,const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m){
      auto phot = PhotoFourVector(react,px,py,pz,m);
      return -phot.M2();
    
    }
    
    template<typename Tp, typename Tm>
      double PolGammaStar(const config::RVecIndexMap& react,const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m){
      
      auto phot = PhotoFourVector(react,px,py,pz,m);
      auto scatele = FourVector(react[names::ScatEleIdx()],px,py,pz,m);
      
      auto q2 = phot.M2();
      auto GammaE = phot.E();
      auto ElScatTh = scatele.Theta();
      
      auto pol = 1./(1.+2.*(1.+GammaE*GammaE/q2)*TMath::Tan(ElScatTh/2.)*TMath::Tan(ElScatTh/2.));
      return pol;
    }
    
    //would like to use a struct here, but define cannot take structs or classes
    //must be basic types or RVecs of basic types
    //HAve added method to deal with structs by defining seperate columns for
    //each member. However this seems slower than doing the calc twice!

    struct ElCMDecay_t{
      double CosTheta=0.;
      double Phi=0.;
    };

 
    template<typename Tp, typename Tm>
    XYZVector ElectroCMDecay(const config::RVecIndexMap& react,const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m){
      //  ElCMDecay_t ElectroCMDecay(const config::RVecIndexMap& react,const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m){
      //CM frame defined by e-scattering
      auto cm = reactkine::CMVector(react,px,py,pz,m);
      auto cmBoost = cm.BoostToCM();
      auto beam = beams::InitialFourVector(react[names::ElectroEleIdx()][0],px,py,pz,m);
      auto mes = FourVector(react[names::MesonsIdx()],px,py,pz,m);
      auto photon = PhotoFourVector(react,px,py,pz,m);
      
 
      PxPyPzMVector CMBeam=boost(beam,cmBoost);
      PxPyPzMVector CMMes=boost(mes,cmBoost);
      PxPyPzMVector CMGamma=boost(photon,cmBoost);
  
      XYZVector zV=CMGamma.Vect().Unit();
      XYZVector yV=CMGamma.Vect().Cross(CMBeam.Vect()).Unit();
      XYZVector xV=yV.Cross(zV).Unit();
      
      XYZVector angles(CMMes.Vect().Dot(xV),CMMes.Vect().Dot(yV),CMMes.Vect().Dot(zV));
      //cout << names::BaryonsIdx() << " " << names::MesonsIdx() << " " << names::VirtGammaIdx() << endl;
      //if(cos(angles.Theta())==1)
      /* cout << "CMBeam: " << CMBeam << endl; */
      /* cout << "CMMeson: " << CMMes << endl; */
      /* cout << "CMGamma: " << CMGamma << endl; */
      /* cout << "zV: " << zV << endl;  */
      /* cout << "Cos(theta): " << cos(angles.Theta()) << endl; */
      /* cout << endl; */
      // ElCMDecay_t result;
      // result.CosTheta=(TMath::Cos(angles.Theta()));
      // result.Phi=angles.Phi();
      // return result;
      return angles;
    }
  
    template<typename Tp, typename Tm>
    double CosThetaCM(const config::RVecIndexMap& react,const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m){
      auto angles = ElectroCMDecay(react,px,py,pz,m);
      return TMath::Cos(angles.Theta());
    }
  
    template<typename Tp, typename Tm>
    double PhiCM(const config::RVecIndexMap& react,const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m){
      auto angles = ElectroCMDecay(react,px,py,pz,m);
      return angles.Phi();
    }

    template<typename Tp, typename Tm>
    XYZVector ElectroProtonRestDecay(const config::RVecIndexMap& react,const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m){
      //  ElCMDecay_t ElectroCMDecay(const config::RVecIndexMap& react,const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m){
      //CM frame defined by e-scattering
      auto  pr = beams::InitialFourVector(react[names::ElectroIonIdx()][0],px,py,pz,m);
      auto prBoost = pr.BoostToCM();
      auto beam = beams::InitialFourVector(react[names::ElectroEleIdx()][0],px,py,pz,m);
      auto mes = FourVector(react[names::MesonsIdx()],px,py,pz,m);
      auto photon = PhotoFourVector(react,px,py,pz,m);
      
 
      PxPyPzMVector prBeam=boost(beam,prBoost);
      PxPyPzMVector prMes=boost(mes,prBoost);
      PxPyPzMVector prGamma=boost(photon,prBoost);
  
      XYZVector zV=-prGamma.Vect().Unit();
      XYZVector yV=prGamma.Vect().Cross(prBeam.Vect()).Unit();
      XYZVector xV=yV.Cross(zV).Unit();
  
      XYZVector angles(prMes.Vect().Dot(xV),prMes.Vect().Dot(yV),prMes.Vect().Dot(zV));
      // ElCMDecay_t result;
      // result.CosTheta=(TMath::Cos(angles.Theta()));
      // result.Phi=angles.Phi();
      // return result;
      return angles;
    }

    template<typename Tp, typename Tm>
    double CosThetaProtonRest(const config::RVecIndexMap& react,const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m){
      auto angles = ElectroProtonRestDecay(react,px,py,pz,m);
      return TMath::Cos(angles.Theta());
    }
    
    template<typename Tp, typename Tm>
    double PhiProtonRest(const config::RVecIndexMap& react,const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m){
      auto angles = ElectroProtonRestDecay(react,px,py,pz,m);
      return angles.Phi();
    }



  }
}
