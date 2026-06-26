#pragma once

//!  Derived class to configure ePIC root files

/*!
  Use ConfigReaction classes to setup data analysis and calculations
  for particular hadronic final states.
  This derived class is configured for Electron scattering reactions
  Derive from this if your experiment is electron scattering
*/
#include "ConfigReaction.h"
#include "RVecHelpers.h"
#include "Beams.h"
#include "BasicKinematics.h"

namespace rad{
  
  namespace config{

 

    //! Class definition

    class PhotoIonReaction : public ConfigReaction {


    public:

      PhotoIonReaction(const std::string_view treeName, const std::string_view fileNameGlob, const ROOT::RDF::ColumnNames_t&  columns ={} ) : ConfigReaction{treeName,fileNameGlob,columns} {

      }
     PhotoIonReaction(const std::string_view treeName, const std::vector<std::string> &filenames, const ROOT::RDF::ColumnNames_t&  columns ={} ) : ConfigReaction{treeName,filenames,columns} {

      }

      /**
       * Make map that links particle names to indices in user functions
       * in C++ functions you can use the RVecIndexMap object indexed by 
       * name of the reaction component you need
       */
      void makeParticleMap() override{
	//note, ordering in arguments, map and names must be maintained
  	Define(names::ReactionMap().data(),
	       [](const int& tar_ion, const int& beam_gam,
		  const RVecI& baryons,const RVecI& mesons){
		 return RVecIndexMap{{tar_ion},{beam_gam},baryons,mesons};},
	       {names::TargetIon().data(),names::BeamGamma().data(),
		names::Baryons().data(),names::Mesons().data()});
      }
      
      /** 
       * Set constant index for beam gamma, and target ion
       * This assumes constant position in collection (e.g in some HepMC3 files)
       * and update the current frame to the aliased one
       */
      void setBeamGammaIndex(const int idx){
	setParticleIndex(names::BeamGamma().data(),idx);
      }
      template<typename Lambda>
      void setBeamGammaIndex(Lambda&& func,const ROOT::RDF::ColumnNames_t & columns = {} ){
	setCurrFrame(CurrFrame().Define(names::BeamGamma().data(),func,columns));
      }
      void setTargetIonIndex(const int idx){
	setParticleIndex(names::TargetIon().data(),idx);
      }

    };//PhotoIonReaction

  }

   namespace photoion{
    /**
     * beam photon 4-vector
     */
    template<typename Tp, typename Tm>
    PxPyPzMVector PhotoFourVector(const config::RVecIndexMap& react, const RVec<Tp> &px, const RVec<Tp> &py, const RVec<Tp> &pz, const RVec<Tm> &m){
      
      return beams::InitialFourVector(react[names::InitialTopIdx()][0],px,py,pz,m);
      
    }
     const std::string_view  BeamIndices() {return Form("%s,%s",names::TargetIon().data(),names::BeamGamma().data())};
  }
}

//Declare we are using this PhotoFourVector in kinematics
using rad::photoion::PhotoFourVector;
using rad::photoion::BeamIndices;
