#pragma once
#include "Beams.h"
#include "BasicKinematics.h"
#include "ConfigReaction.h"

void ReactionKinematics(){}

namespace rad{
  namespace reactkine{
  
    ///\brief create a new particle and add it to the momentum vectors
    ///return the index to be used to access the components
    ///this also allows it to be used with Define which requires a return
    template<typename Tp, typename Tm>
    int ParticleCreateBySum(const RVecI& isum, RVec<Tp> &px, RVec<Tp> &py, RVec<Tp> &pz, RVec<Tm> &m,const RVecI& iafter){
      //sum the 4-vectors
      // std::cout<<"ParticleCreateBySum "<<isum<<m<<pz<<" "<<m.size()<<std::endl;
      PxPyPzMVector p4;
      SumFourVector(p4,isum,px,py,pz,m);
      //make particle id = last entry
      auto idx = px.size();
      //add new components
      px.push_back(p4.X());
      py.push_back(p4.Y());
      pz.push_back(p4.Z());
      m.push_back(p4.M());
      return idx;
    }
    ///\brief create a new particle and add it to the momentum vectors
    ///return the index to be used to access the components
    ///this also allows it to be used with Define which requires a return
    template<typename Tp, typename Tm>
    int  ParticleCreateByMiss(const int ibot,const int itop,const RVecI& ineg, RVec<Tp> &px, RVec<Tp> &py, RVec<Tp> &pz, RVec<Tm> &m,const RVecI& iafter){
      //sum the 4-vectors
      auto p4 = beams::InitialFourVector(itop,px,py,pz,m);
      p4+=beams::InitialFourVector(ibot,px,py,pz,m);
       
      SubtractFourVector(p4,ineg,px,py,pz,m);
 
      //make particle id = last entry
      auto idx = px.size();

      //add new components
      px.push_back(p4.X());
      py.push_back(p4.Y());
      pz.push_back(p4.Z());
      m.push_back(p4.M());
      return idx;
    }
 
    
   ///\brief missing mass fo reaction = top+bot -  neg[i]
    //const config::RVecIndexMap react must be copied for thread safety.
    template<typename Tp, typename Tm>
    double FourVectorMissMassCalc(const config::RVecIndexMap react,const RVecI ineg,const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m)
    {
      auto psum = beams::InitialFourVector(react[names::InitialTopIdx()][0],px,py,pz,m);
      psum+=beams::InitialFourVector(react[names::InitialBotIdx()][0],px,py,pz,m);
      SubtractFourVector(psum,ineg,px,py,pz,m);
      return psum.M();
    }

   ///\brief missing mass sqaured of reaction = top+bot -  neg[i]
    //const config::RVecIndexMap react must be copied for thread safety.
    template<typename Tp, typename Tm>
    double FourVectorMissMass2Calc(const config::RVecIndexMap react,const RVecI ineg,const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m)
    {
      auto psum = beams::InitialFourVector(react[names::InitialTopIdx()][0],px,py,pz,m);
      psum+=beams::InitialFourVector(react[names::InitialBotIdx()][0],px,py,pz,m);
      SubtractFourVector(psum,ineg,px,py,pz,m);
      return psum.M2();
    }
    
    ///\brief missing transverse momentum of reaction = top+bot -  neg[i]
    //const config::RVecIndexMap react must be copied for thread safety.
   template<typename Tp, typename Tm>
    double FourVectorMissPtCalc(const config::RVecIndexMap react,const RVecI ineg,const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m)
    {
      auto psum = beams::InitialFourVector(react[names::InitialTopIdx()][0],px,py,pz,m);
      psum+=beams::InitialFourVector(react[names::InitialBotIdx()][0],px,py,pz,m);
      SubtractFourVector(psum,ineg,px,py,pz,m);
      return psum.Pt();
    }

   ///\brief missing transverse momentum of reaction = top+bot -  neg[i]
    //const config::RVecIndexMap react must be copied for thread safety.
   template<typename Tp, typename Tm>
    double FourVectorMissPzCalc(const config::RVecIndexMap react,const RVecI ineg,const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m)
    {
      auto psum = beams::InitialFourVector(react[names::InitialTopIdx()][0],px,py,pz,m);
      psum+=beams::InitialFourVector(react[names::InitialBotIdx()][0],px,py,pz,m);
      SubtractFourVector(psum,ineg,px,py,pz,m);
      return psum.Pz();
    }

  ///\brief missing transverse momentum of reaction = top+bot -  neg[i]
    //const config::RVecIndexMap react must be copied for thread safety.
    template<typename Tp, typename Tm>
    double FourVectorMissThetaCalc(const config::RVecIndexMap react,const RVecI ineg,const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m)
    {
      auto psum = beams::InitialFourVector(react[names::InitialTopIdx()][0],px,py,pz,m);
      psum+=beams::InitialFourVector(react[names::InitialBotIdx()][0],px,py,pz,m);
      SubtractFourVector(psum,ineg,px,py,pz,m);
      return psum.Theta();
    }
    
    ///\brief missing transverse momentum of reaction = top+bot -  neg[i]
    //const config::RVecIndexMap react must be copied for thread safety.
   template<typename Tp, typename Tm>
    double FourVectorMissPhiCalc(const config::RVecIndexMap react,const RVecI ineg,const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m)
    {
      auto psum = beams::InitialFourVector(react[names::InitialTopIdx()][0],px,py,pz,m);
      psum+=beams::InitialFourVector(react[names::InitialBotIdx()][0],px,py,pz,m);
      SubtractFourVector(psum,ineg,px,py,pz,m);
      return psum.Phi();
    }


    ///\brief missing transverse momentum of reaction = top+bot -  neg[i]
    //const config::RVecIndexMap react must be copied for thread safety.
   template<typename Tp, typename Tm>
    double FourVectorMissPCalc(const config::RVecIndexMap react,const RVecI ineg,const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m)
    {
      auto psum = beams::InitialFourVector(react[names::InitialTopIdx()][0],px,py,pz,m);
      psum+=beams::InitialFourVector(react[names::InitialBotIdx()][0],px,py,pz,m);
      SubtractFourVector(psum,ineg,px,py,pz,m);
      return psum.P();
    }

    ///\brief missing transverse momentum of reaction = top+bot -  neg[i]
    //const config::RVecIndexMap react must be copied for thread safety.
   template<typename Tp, typename Tm>
    double FourVectorMissECalc(const config::RVecIndexMap react,const RVecI ineg,const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m)
    {
      auto psum = beams::InitialFourVector(react[names::InitialTopIdx()][0],px,py,pz,m);
      psum+=beams::InitialFourVector(react[names::InitialBotIdx()][0],px,py,pz,m);
      SubtractFourVector(psum,ineg,px,py,pz,m);
      return psum.E();
    }

   ///\brief functions to compute standard reaction kinematics
    template<typename Tp, typename Tm>
    PxPyPzMVector CMVector(const config::RVecIndexMap& react, const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m){
      RVecI icm=react[names::BaryonsIdx()];
      icm.insert(icm.end(),react[names::MesonsIdx()].begin(),react[names::MesonsIdx()].end());
    
      return FourVector(icm,px,py,pz,m);
    }

    /**
     * reaction baryon 4-vector
     */
    template<typename Tp, typename Tm>
    PxPyPzMVector BaryonFourVector(const config::RVecIndexMap& react, const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m){
      return FourVector(react[names::BaryonsIdx()],px,py,pz,m);
    }
    /**
     * reaction meson 4-vector
     */
    template<typename Tp, typename Tm>
    PxPyPzMVector MesonFourVector(const config::RVecIndexMap& react, const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m){
      return FourVector(react[names::MesonsIdx()],px,py,pz,m);
    }

   
    template<typename Tp, typename Tm>
    double T0(const config::RVecIndexMap& react,
	  const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m){

      auto tar = beams::InitialFourVector(react[names::InitialBotIdx()][0],px,py,pz,m);
      auto bar = FourVector(react[names::BaryonsIdx()],px,py,pz,m);
      //generate CM from sum of final state meson and baryon particles
      auto cm = CMVector(react,px,py,pz,m);
      auto cmBoost = cm.BoostToCM();
      PxPyPzMVector  CMTar=boost(tar,cmBoost);
      PxPyPzMVector  CMBar=boost(bar,cmBoost);
    
      //return  M1*M1 + M3*M3  - 2 * ( E1*E3 -p1*p3*costh );
      Double_t t0 = CMBar.M2() + CMTar.M2() - 2*(CMBar.E()*CMTar.E() - CMBar.P()*CMTar.P() ) ;
      /* cout << "Inside T0 Func" << endl; */
      /* cout << "CMBar: " << CMBar << endl; */
      /* cout << "CMTar: " << CMTar << endl; */
      return t0;
    }
  
   
    ///\brief return 4 momentum transfer squared of "in particles" - "out particles" on bottom vertex
    template<typename Tp, typename Tm>
    double TBot(const config::RVecIndexMap& react,const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m){
      auto psum = beams::InitialFourVector(react[names::InitialBotIdx()][0],px,py,pz,m);
      /* cout << "Inside TBot Func" << endl; */
      /* cout << "pbeam: " << psum << endl; */
      SubtractFourVector(psum,react[names::BaryonsIdx()],px,py,pz,m); 
      /* cout << "pbeam_sub_pprime: " << psum << endl; */

      return - (psum.M2());
    }
    ///\brief return 4 momentum transfer squared of "in particles" - "out particles" on top vertex
    template<typename Tp, typename Tm>
    double TTop(const config::RVecIndexMap& react,
	    const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m){
      //Get photon 4-vector
      auto phot=PhotoFourVector(react,px,py,pz,m);
      //Get meson (Top) 4-vector
      auto meso=FourVector(react[names::MesonsIdx()],px,py,pz,m);
      //subtract
      /* cout << "Inside TTop Func" << endl; */
      /* cout << "phot: " << phot << endl; */
      /* cout << "meso: " << meso << endl; */
      auto psum = phot-meso;
      //return t
      return - (psum.M2());
    }

      //For some reason this calculation only works when boosted into
      //proton beam rest frame......
      //auto pbeam = beams::InitialFourVector(react[names::InitialBotIdx()][0],px,py,pz,m);
      //auto boov=pbeam.BoostToCM();
      // auto phot=boost(PhotoFourVector(react,px,py,pz,m),boov);
      //auto meso=boost(FourVector(react[names::MesonsIdx()],px,py,pz,m),boov);

  
    ///\brief return 4 momentum transfer squared, t, minus t0 (or tmin) of "in particles" - "out particles"
    template<typename Tp, typename Tm>
    double TPrimeBot(const config::RVecIndexMap& react,const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m){
      auto tbot = TBot(react,px,py,pz,m);
      auto t0 = T0(react,px,py,pz,m);
      /* cout << "TBot: " << tbot << endl; */
      /* cout << "T0: " << t0 << endl; */
      /* cout << "TpBot; " << tbot+t0 << endl; */
      /* cout << endl; */
      return tbot + t0;
    //return TBot(react,px,py,pz,m) + T0(react,px,py,pz,m);
    }
    ///\brief return 4 momentum transfer squared, t, minus t0 (or tmin) of "in particles" - "out particles"
    template<typename Tp, typename Tm>
    double TPrimeTop(const config::RVecIndexMap& react,const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m){
      auto ttop = TTop(react,px,py,pz,m);
      auto t0 = T0(react,px,py,pz,m);
      /* cout << "TTop: " << ttop << endl; */
      /* cout << "T0: " << t0 << endl; */
      /* cout << "TpTot; " << ttop+t0 << endl; */
      /* cout << endl; */
      return ttop + t0;
      //return TTop(react,px,py,pz,m) + T0(react,px,py,pz,m);
    }
  


  }

}
