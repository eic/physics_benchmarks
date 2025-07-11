#pragma once

#include "ConfigReaction.h"

namespace rad{
  namespace config{

    class ParticleCreator {

    public:
      
      ParticleCreator() = default;
    ParticleCreator(ConfigReaction& cr):_reaction{&cr}{};
      
      
      //////////////////////////////////////////////////////////////////
      void Beam(const string& name,const string& p4name) const{
	//empty parts string as not dependent on others
	DefineParticle(name,std::vector<std::string>(),Form("rad::beams::ParticleFixedBeam(%s",p4name.data() ));
      }
      //////////////////////////////////////////////////////////////////
      // void Beam(const string& name,const string& p4name,const std::vector<std::string>& parts){
      // 	//empty parts string as not dependent on others
      // 	DefineParticle(name,parts,Form("rad::beam::ParticleEventBeam("));
      // }
      //////////////////////////////////////////////////////////////////
      void Sum(const string& name,const std::vector<std::string>& parts) const{
	//adding ( at end of string allows extra arguments to be added, like for Miss
	DefineParticle(name,parts,Form("rad::reactkine::ParticleCreateBySum(%s",VectorToString(parts).data()) );
      }
     //////////////////////////////////////////////////////////////////
      void Diff(const string& name,const std::vector<std::string>& pos,const std::vector<std::string>& neg) const{
	//adding ( at end of string allows extra arguments to be added, like for Miss
	auto parts = pos;
	parts.insert(parts.end(), neg.begin(), neg.end());
	DefineParticle(name,parts,Form("rad::config::ParticleCreateByDiff(%s,%s",VectorToString(pos).data(),VectorToString(neg).data()) );
      }
      //////////////////////////////////////////////////////////////////
      void Miss(const string& name,const std::vector<std::string>& parts) const{
	DefineParticle(name,parts,Form("rad::reactkine::ParticleCreateByMiss(%s,%s",BeamIndices().data(),VectorToString(parts).data()));
      }
      //
      
      //////////////////////////////////////////////////////////////////
      void DefineParticle(const string& name,const std::vector<std::string> parts,const string& funcExpr) const{
	
	//store all defined particle names
	std::vector<std::string> names;

	//loop over ConfigReaction types and define this particle for each
	auto types = _reaction->GetTypes();
	for(auto &atype:types){
	  //Make sure any created particle uses type_idx
	  //This ensures its Create function is called prior to this one
	  string sum ="{";
	  for(std::string p:parts){
	    if( std::find(_created.begin(),_created.end(),p)!=_created.end() ){
	      std::string type_p = atype.first + p;
	      p=type_p;
	    }
	    sum=(sum+p+",");
	  }
	  if(parts.empty()==false) sum.pop_back(); //remove last ,
	  sum+='}';

	  //Note we give all created particles as argument to ensure creation order
	  auto type_created=_created;
	  for(auto& col: type_created){
	    col=atype.first+col;
	  }
	  auto after_cols  = VectorToString(type_created);
	  
	  //format args "func(idxs,components,after_idxs")
	  //  TString type_expr = Form("%s%s,%s,%s)",funcExpr.data(),sum.data(),atype.second["components_p4"].data(),after_cols.data());
 	  TString type_expr = Form("%s,%s,%s)",funcExpr.data(),atype.second["components_p4"].data(),after_cols.data());
	  names.push_back(atype.first + name.data());
	  _reaction->Define(atype.first + name.data(),type_expr.Data());
	}
	_created.push_back(name);
	
	auto snames = VectorToString(names);
	//define name as the first type entry in names
	//this function ensures all type create particles are called at same time
	//Also we can use just name rather than type_name which should ahve same value
	//for all types
	//std::cout<<"Sum names "<<snames<<" "<<Form("ROOT::RVecU%s[0]",snames.data())<<std::endl;
	_reaction->Define(name.data(),Form("ROOT::RVecI%s[0]",snames.data()));
	Reaction()->AddParticleName(name);
	
      }
    
      //////////////////////////////////////////////////////////////////
      std::string VectorToString(const std::vector<std::string>& parts) const{
	if(parts.empty()==true) return "{}";
      
	string toString ="{";
	for(const auto& p:parts){
	  toString=(toString+p+",");
	}
	toString.pop_back(); //remove last ,
	toString+='}';
	return toString;
      }

      void SetReaction(ConfigReaction* reaction){_reaction=reaction;}
      ConfigReaction* Reaction() const { return _reaction;}
      
    private:
      
      mutable ConfigReaction* _reaction=nullptr;
      mutable std::vector<string> _created;
      
    };

    ///\brief create a new particle and add it to the momentum vectors
    ///return the index to be used to access the components
    ///this also allows it to be used with Define which requires a return

    ///Create a particle as the diffence between ipos particles and ineg
    template<typename Tp, typename Tm>
    int  ParticleCreateByDiff(const RVecI &ipos, const RVecI &ineg, RVec<Tp> &px, RVec<Tp> &py, RVec<Tp> &pz, RVec<Tm> &m,const RVecI& iafter){
      //sum the 4-vectors
      auto p4 = FourVector(ipos,px,py,pz,m);
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
    
  /* template<typename Tp, typename Tm> */
  /*   int  ParticleCreate( Tp &tpx, Tp &tpy, Tp &tpz, Tp &tm, RVec<Tp> &px, RVec<Tp> &py, RVec<Tp> &pz, RVec<Tm> &m, const RVecI& iafter){ */
  /*   auto idx = px.size(); */
    
  /*   px.push_back(tpx); */
  /*   py.push_back(tpy); */
  /*   pz.push_back(tpz); */
  /*   m.push_back(tm); */
  /*   return idx; */
  /* } */
  
  
  }//end config
}
