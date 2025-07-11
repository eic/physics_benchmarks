#pragma once

//#include "Beams.h"
#include "BasicKinematics.h"
#include "ConfigReaction.h"


namespace rad{
  namespace gn2s0s0s12{
    
   struct Angles_t{
      double CosTheta=0.;
      double Phi=0.;
    };

   ///\brief functions to compute standard reaction kinematics

    /**
    * calculate CM kinematics from beam
    */
   template<typename Tp, typename Tm>
    PxPyPzMVector PhotoCMVector(const config::RVecIndexMap& react, const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m){

     //Note PhotoFourVector must be defined in the reaction config file
     //e.g. in PhotoIonReaction.h or ElectroIonReaction.h
     return beams::InitialFourVector(react[names::InitialBotIdx()][0],px,py,pz,m) +
       PhotoFourVector(react,px,py,pz,m);
 
     }

   /**
    * calculate CM decay angles
    */
    template<typename Tp, typename Tm>
    Angles_t PhotoCMDecay(const config::RVecIndexMap& react,const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m){

      auto cm = PhotoCMVector(react,px,py,pz,m);
      auto cmBoost = cm.BoostToCM();
      auto mes = FourVector(react[names::MesonsIdx()],px,py,pz,m);
      PxPyPzMVector cm_mes=boost(mes,cmBoost);

      Angles_t result;
      result.CosTheta=(TMath::Cos(cm_mes.Theta()));
      result.Phi=cm_mes.Phi();
      return result;
     }
    /**
    * calculate Helicity decay angles
    * z-axis along -baryon in meson rest frame
    */
    template<typename Tp, typename Tm>
    Angles_t PhotoHelicityDecay(const config::RVecIndexMap& react,const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m){

      auto baryon = reactkine::BaryonFourVector(react,px,py,pz,m);
      auto meson = reactkine::MesonFourVector(react,px,py,pz,m);

      auto decBoost = meson.BoostToCM();
      //vectors in rest/decay frame of meson
      auto decBar=boost(baryon,decBoost);
      auto decGamma=boost(PhotoFourVector(react,px,py,pz,m),decBoost);
      
      XYZVector  zV=-decBar.Vect().Unit();
      XYZVector  yV=decBar.Vect().Cross(decGamma.Vect()).Unit();
      XYZVector  xV=yV.Cross(zV).Unit();

      //four vector of first [0] decay product
      auto child1 = FourVector(react[names::MesonsIdx()][0],px,py,pz,m);
      auto decChild1=boost(child1,decBoost);

      //calculate decay angles
      XYZVector angles(decChild1.Vect().Dot(xV),decChild1.Vect().Dot(yV),decChild1.Vect().Dot(zV));
      //store in angles struct
      Angles_t result;
      result.CosTheta=TMath::Cos(angles.Theta());
      result.Phi=angles.Phi();
      return result;
    }
    template<typename Tp, typename Tm>
    double CosThetaHel(const config::RVecIndexMap& react,const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m){
      auto angles = PhotoHelicityDecay(react,px,py,pz,m);
      return angles.CosTheta;
    }
     template<typename Tp, typename Tm>
    double PhiHel(const config::RVecIndexMap& react,const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m){
      auto angles = PhotoHelicityDecay(react,px,py,pz,m);
      return angles.Phi;
    }
     /**
    * calculate GJ decay angles
    * z-axis along gamma direction in meson rest frame
    */
    template<typename Tp, typename Tm>
    Angles_t PhotoGJDecay(const config::RVecIndexMap& react,const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m){

      auto baryon = reactkine::BaryonFourVector(react,px,py,pz,m);
      auto meson = reactkine::MesonFourVector(react,px,py,pz,m);

      auto decBoost = meson.BoostToCM();
      //vectors in rest/decay frame of meson
      auto decBar=boost(baryon,decBoost);
      auto decGamma=boost(PhotoFourVector(react,px,py,pz,m),decBoost);
      
      XYZVector  zV=decGamma.Vect().Unit();
      XYZVector  yV=decBar.Vect().Cross(decGamma.Vect()).Unit();
      XYZVector  xV=yV.Cross(zV).Unit();

      //four vector of first [0] decay product
      auto child1 = FourVector(react[names::MesonsIdx()][0],px,py,pz,m);
      auto decChild1=boost(child1,decBoost);

      //calculate decay angles
      XYZVector angles(decChild1.Vect().Dot(xV),decChild1.Vect().Dot(yV),decChild1.Vect().Dot(zV));
      //store in angles struct
      Angles_t result;
      result.CosTheta=TMath::Cos(angles.Theta());
      result.Phi=angles.Phi();
      return result;
    }
    template<typename Tp, typename Tm>
    double CosThetaGJ(const config::RVecIndexMap& react,const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m){
      auto angles = PhotoGJDecay(react,px,py,pz,m);
      return angles.CosTheta;
    }
     template<typename Tp, typename Tm>
    double PhiGJ(const config::RVecIndexMap& react,const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m){
      auto angles = PhotoGJDecay(react,px,py,pz,m);
      return angles.Phi;
    }
  }
}
