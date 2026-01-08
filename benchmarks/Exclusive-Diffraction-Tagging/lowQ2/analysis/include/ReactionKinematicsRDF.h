#pragma once
#include "ReactionKinematics.h"

void ReactionKinematicsRDF(){}

namespace rad{
  namespace rdf{

    void Particle_CreateBySum(config::ConfigReaction& cr,const string& name, const string& sum){
      //cr.Define(name,[](){return 6;},{}); //dummy function. particle index will be set in the following
      cr.DefineForAllTypes(name, Form("rad::reactkine::ParticleCreateBySum(%s,components_p4)",sum.data()));
      cr.AliasToPrimaryType(name);
    }
    
    void Particle_CreateByMiss(config::ConfigReaction& cr,const string& name, const string& neg){
      cr.DefineForAllTypes(name, Form("rad::reactkine::ParticleCreateByMiss(%s,%s,components_p4)",names::ReactionMap().data(),neg.data()));
      cr.AliasToPrimaryType(name);
    }
    
    void MissP(config::ConfigReaction& cr,const string& name, const string_view neg){
      cr.DefineForAllTypes(name, Form("rad::reactkine::FourVectorMissPCalc(%s,%s,components_p4)",names::ReactionMap().data(),neg.data()));
    }
    
    void MissPt(config::ConfigReaction& cr,const string& name, const string_view neg){
      cr.DefineForAllTypes(name, Form("rad::reactkine::FourVectorMissPtCalc(%s,%s,components_p4)",names::ReactionMap().data(),neg.data()));
    }
    
    void MissPz(config::ConfigReaction& cr,const string& name, const string_view neg){
      cr.DefineForAllTypes(name, Form("rad::reactkine::FourVectorMissPzCalc(%s,%s,components_p4)",names::ReactionMap().data(),neg.data()));
    }
    
    void MissTheta(config::ConfigReaction& cr,const string& name, const string_view neg){
      cr.DefineForAllTypes(name, Form("rad::reactkine::FourVectorMissThetaCalc(%s,%s,components_p4)",names::ReactionMap().data(),neg.data()));
    }
    
    void MissPhi(config::ConfigReaction& cr,const string& name, const string_view neg){
      cr.DefineForAllTypes(name, Form("rad::reactkine::FourVectorMissPhiCalc(%s,%s,components_p4)",names::ReactionMap().data(),neg.data()));
    }
    
    void MissE(config::ConfigReaction& cr,const string& name, const string_view neg){
      cr.DefineForAllTypes(name, Form("rad::reactkine::FourVectorMissECalc(%s,%s,components_p4)",names::ReactionMap().data(),neg.data()));
    }
    
    void MissMass(config::ConfigReaction& cr,const string& name, const string_view neg){
      cr.DefineForAllTypes(name, Form("rad::reactkine::FourVectorMissMassCalc(%s,%s,components_p4)",names::ReactionMap().data(),neg.data()));
    }
    
    void MissMass2(config::ConfigReaction& cr,const string& name, const string_view neg){
      cr.DefineForAllTypes(name, Form("rad::reactkine::FourVectorMissMass2Calc(%s,%s,components_p4)",names::ReactionMap().data(),neg.data()));
    }
    
    void TBot(config::ConfigReaction& cr,const string& name){
      cr.DefineForAllTypes(name, Form("rad::reactkine::TBot(%s,components_p4)",names::ReactionMap().data()));
    }
    void TTop(config::ConfigReaction& cr,const string& name){
      cr.DefineForAllTypes(name, Form("rad::reactkine::TTop(%s,components_p4)",names::ReactionMap().data()));
    }
    void TPrimeTop(config::ConfigReaction& cr,const string& name){
      cr.DefineForAllTypes(name, Form("rad::reactkine::TPrimeTop(%s,components_p4)",names::ReactionMap().data()));
    }
    void TPrimeBot(config::ConfigReaction& cr,const string& name){
      cr.DefineForAllTypes(name, Form("rad::reactkine::TPrimeBot(%s,components_p4)",names::ReactionMap().data()));
    }
    
    //instead of neg use all meson, baryon and scat_elec by default for this function
    //work in progress
    /* void ConservationCheck(config::ConfigReaction& cr){ */
    /*   cr.DefineForAllTypes(name, Form("rad::reactkine::FourVectorMissMassCalc(%s,%s,components_p4)",names::ReactionMap().data(),neg.data())); */
    
    /* } */


  }//rdf
}//rad
