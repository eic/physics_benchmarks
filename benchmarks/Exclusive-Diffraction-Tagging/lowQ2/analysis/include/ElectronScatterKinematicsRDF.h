#pragma once
#include "ElectronScatterKinematics.h"

void ElectronScatterKinematicsRDF(){}

namespace rad{
  namespace rdf{
  
    void EGammaStar(config::ConfigReaction& cr,const string& name){
      cr.DefineForAllTypes(name, Form("rad::electroion::PhotoFourVector(%s,components_p4).E()",names::ReactionMap().data()));
    }
    void Q2(config::ConfigReaction& cr,const string& name){
      cr.DefineForAllTypes(name, Form("rad::electro::Q2(%s,components_p4)",names::ReactionMap().data()));
    }
    void CosThetaCM(config::ConfigReaction& cr,const string& name){
      cr.DefineForAllTypes(name, Form("rad::electro::CosThetaCM(%s,components_p4)",names::ReactionMap().data()));
    }
    void PhiCM(config::ConfigReaction& cr,const string& name){
      cr.DefineForAllTypes(name, Form("rad::electro::PhiCM(%s,components_p4)",names::ReactionMap().data()));
    }
    void CMAngles(config::ConfigReaction& cr,const string& name){
      CosThetaCM(cr,name+"_CosTheta");
      PhiCM(cr,name+"_Phi");
    }

    void CosThetaProtonRest(config::ConfigReaction& cr,const string& name){
      cr.DefineForAllTypes(name, Form("rad::electro::CosThetaProtonRest(%s,components_p4)",names::ReactionMap().data()));
    }
    void PhiProtonRest(config::ConfigReaction& cr,const string& name){
      cr.DefineForAllTypes(name, Form("rad::electro::PhiProtonRest(%s,components_p4)",names::ReactionMap().data()));
    }
    void PRAngles(config::ConfigReaction& cr,const string& name){
      CosThetaProtonRest(cr,name+"_CosTheta");
      PhiProtonRest(cr,name+"_Phi");
    }
    
    void PolGammaStar(config::ConfigReaction& cr,const string& name){
      cr.DefineForAllTypes(name, Form("rad::electro::PolGammaStar(%s,components_p4)",names::ReactionMap().data()));
    }
    

  }//rdf
}//rad
