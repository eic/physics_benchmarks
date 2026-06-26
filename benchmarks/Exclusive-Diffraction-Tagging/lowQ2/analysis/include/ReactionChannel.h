#pragma once

//!  Base configuration class for specific reaction channel
/*!
     Used to configure reaction (ConfigReaction object) 
     to a specific final state. Based on PDG values for
     the top (meson) and bottom (baryon) vertices
 
*/
#include "ConfigReaction.h"
#include "ParticleCreator.h"
#include "Indicing.h"

//due to using PhotoFourVector I am not sure I can include both of these
//#include "PhotoIonReaction.h"
//#include "ElectroIonReaction.h"
//so instead declare rather than define the classes
//the definition will be picked up at runtime with
//the correct #include in the run script used to create
//the actual reaction

class PhotoIonReaction;
class ElectroIonReaction;

namespace rad{
  namespace config{
    using rad::names::data_type::Rec;
 
    //need to add functions for other reaction types here if needed
    void DoElectroReaction(rad::config::ElectroIonReaction* reaction,std::string pid);

    //need to add functions for other data formats here if needed
    void UsePythiaOccurances(std::map<int,int>& occur );
    void UseEpicOccurances(std::map<int,int>& occur );
    
    class ReactionChannel{

  
    public:

      template<typename Lambda>
	ReactionChannel(rad::config::ElectroIonReaction &rea,std::vector<int> mes,std::vector<int> bar,std::string pid,Lambda&& func):_reaction{&rea},_mesons{mes},_baryons{bar},_pid{pid}{

	//UseEpicOccurances(_occurances);
	//UsePythiaOccurances(_occurances);
	func(_reaction,_pid);
	ConfigureParticles();
      }
  
       std::vector<std::string> InsertScatteredElectron(std::vector<std::string> particles) const{
	auto el = dynamic_cast<rad::config::ElectroIonReaction*>(_reaction);
	if(el==nullptr) return particles; 

	//	particles.insert(particles.find('}'), string(",")+names::ScatEle().data()+"}");
	//cout<<"InsertScatteredElectron particles "<<particles<<endl;
	particles.push_back(names::ScatEle().data());
	return particles;
      }
      
      void ConfigureParticles(){

	//Assign an index to each particle and create list of names	
	SetParticleIndex();
    
	//Group particles into top and bottom vertices
	//aka Meson and Baryon components
	//this is required for calcualting reaction kinematics
	//e.g. t distributions
    
	//create intermediate meson and baryon particles
	//	rad::config::ParticleCreator particles{*_reaction};
	auto& particles = _reaction->Particles();
	
	if(_mesons.empty()){
	  //create missing particle from beams and baryon
	  particles.Miss("idxMeson",InsertScatteredElectron(_name_baryons));
	}
	else{
	  //create particle from sum of meson particles
	  particles.Sum("idxMeson",_name_mesons);
	}
	if(_baryons.empty()){
	  //create missing particle from beams and meson
	  particles.Miss("idxBaryon",InsertScatteredElectron(_name_mesons));
	}
	else{
	  //create particle from sum of baryons particles
	  particles.Sum("idxBaryon",_name_baryons);
	}
	  
	_cutString.pop_back();
    
	_reaction->setBaryonParticles({"idxBaryon"});
        _reaction->setMesonParticles({"idxMeson"});
 	// _reaction->setBaryonParticles({"Baryon0"});
        // _reaction->setMesonParticles({"Meson0","Meson1"});
   
	//must call this after all particles are configured
	_reaction->makeParticleMap();

	//apply a filter on the requested particles
	auto df0 = _reaction->CurrFrame();
	df0=df0.Define("reaction_topo",_cutString).Filter("reaction_topo==1");
	_reaction->setCurrFrame(df0);
      }

      void SetParticleIndex(){
    
	_name_mesons=SetIndexByOccurance(_mesons,"Meson");
	_name_baryons=SetIndexByOccurance(_baryons,"Baryon");
    
	for(auto& oc:_occurances){
	  auto cut = Form("(rad::helpers::Count(%s,%d)>=%d)*",_pid.data(),oc.first,oc.second-1);
	  _cutString+=cut;
	}

      };
      string CutParticleCondition(const string& var=Rec()+"pmag",const string& condition =">0."){
	string cut;
	for(auto& name:_name_mesons){
	  cut+=Form("%s[%s]%s&&",var.data(),name.data(),condition.data());
	}
	for(auto& name:_name_baryons){
	  cut+=Form("%s[%s]%s&&",var.data(),name.data(),condition.data());
	}
	cut.pop_back();cut.pop_back();//remove last &&
	return cut;
      }
      std::vector<string> SetIndexByOccurance(const std::vector<int>& parts,const string& type){
	
	std::vector<std::string> names;
	for(uint im=0;im<parts.size();im++){
	  if(_occurances.find(parts[im]) == _occurances.end()){
	    //start at 1, as 0th occurance does not exist!
	    _occurances[parts[im]]=1;
	  }
	  auto n = _occurances[parts[im]];
	  auto name = Form("%s%d",type.data(),im);
	  std::cout<<" SetIndexByOccurance( "<<name<<" "<<_pid<<std::endl;
	  names.push_back(name);
	  _reaction->setParticleIndex(name,rad::indice::useNthOccurance(n,parts[im]),{_pid},parts[im]);
	  _occurances[parts[im]]++;
	}
	return names;
	
      }
  
    private:

      // rad::config::ConfigReaction* _reaction={nullptr};
      rad::config::ElectroIonReaction* _reaction={nullptr};
      const std::vector<int> _mesons;
      const std::vector<int> _baryons;
      std::string _cutString;
      std::vector<std::string> _name_mesons;
      std::vector<std::string> _name_baryons;
      std::map<int,int> _occurances;
      std::string _pid;// column name
    };

    //specific Configuration for electro reactions
    void DoElectroReaction(rad::config::ElectroIonReaction* reaction,std::string pid){
      
      auto el = dynamic_cast<rad::config::ElectroIonReaction*>(reaction);
      if(el==nullptr){
	std::cerr<<"Error : DoElectroReaction was not give an ElectroIonReaction"<<std::endl;
	exit(0);
      }
      //Assign particles names and indices
      //indicing comes from ordering in hepmc file
      //el->setBeamIonIndex(rad::beams::InitBotFix());
      // el->setBeamElectronIndex(rad::beams::InitTopFix());
      //el->setBeamElectron(rad::indice::useNthOccurance(0,11),{"tru_pid"});
      //el->setBeamIon(rad::indice::useNthOccurance(0,2212),{"tru_pid"});
      el->setScatElectronIndex(rad::indice::useNthOccurance(1,11),{pid});


    }//DoElectroReaction()

    //specifc configuration for pythia hepmc files
    void UsePythiaOccurances(std::map<int,int>& occur ){

      //we remove beam particles from lists 
      //pythia electro + proton beam
      occur[2212]=3;
      occur[11]=3;
      //+ virtual photon
      occur[22]=2;
    }
     //specifc configuration for pythia hepmc files
    void UseEpicOccurances(std::map<int,int>& occur ){

      //we remove beam particles from lists 
     
    }
    
  }
}
