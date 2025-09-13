#pragma once


//!  Helper functions for names of reaction parts

namespace rad{
  namespace names{

    /**
     * Names used to identify reaction components
     */
    constexpr const std::string_view  ReactionMap()  {return "reaction_map"; }
    constexpr const std::string_view  Mesons()  {return "meson"; }
    constexpr const std::string_view  Baryons() {return "baryon";}

    ///Comments not using string_view as need to return a string
    // returning static const seems to be slower, so stick with
    // returning string which should actually be slower 
    //constexpr const std::string_view  ScatEle() {return "scat_ele";}
    //Lilely due to some small string optimisations
    const std::string  ScatEle() {return "scat_ele";}
    // inline const std::string& ScatEle() {
    //   static const std::string val = "scat_ele";
    //   return val;
    // }

    const std::string  BeamEle() {return "beam_ele";}
    const std::string  BeamIon() {return "beam_ion";}
    constexpr const std::string_view  TargetIon() {return "tar_ion";}
    constexpr const std::string_view  BeamGamma() {return "beam_gamma";}
    const std::string  P4BeamEle() {return "p4beam_ele";}
    const std::string  P4BeamIon() {return "p4beam_ion";}
    const std::string  P4TargetIon() {return "p4tar_ion";}
    const std::string  P4BeamGamma() {return "p4beam_gamma";}
    //constexpr const std::string_view  Beams() {return "{beam_ele,beam_ion}";}

    /**
     * Links to reaction component links for users c++ functions
     */
    enum class InitGroup{Bot,Top}; //ordering below must match this
    enum class ElectroGroup{ BeamIon,BeamEle,ScatEle=4,VirtGam}; //4=>order in particleMap (after baryons and mesons)
    enum class PhotoGroup{ TarIon,BeamGam}; //ordering below must match this
    enum class FinalGroup{ Baryons=2,Mesons}; //ordering below must match this

    
    constexpr uint  InitialTopIdx() {return static_cast<uint>(InitGroup::Top);}
    constexpr uint  InitialBotIdx() {return static_cast<uint>(InitGroup::Bot);}

    constexpr uint  ElectroIonIdx() {return InitialBotIdx();}
    constexpr uint  ElectroEleIdx() {return InitialTopIdx();}
    
    constexpr uint  ScatEleIdx() {return static_cast<uint>(ElectroGroup::ScatEle);}
    constexpr uint  VirtGammaIdx() {return static_cast<uint>(ElectroGroup::VirtGam);}
    
    constexpr uint  BaryonsIdx() {return static_cast<uint>(FinalGroup::Baryons);}
    constexpr uint  MesonsIdx() {return static_cast<uint>(FinalGroup::Mesons);}
    
    constexpr uint  PhotoGammaIdx() {return InitialTopIdx();}
    constexpr uint  PhotoIonIdx() {return InitialBotIdx();}

    /**
     * Types of data
     */
    namespace data_type{
      const std::string  Rec() {return "rec_";}
      const std::string  Truth() {return "tru_";}
      const std::string  MC() {return "mc_";}
    }
    
  }//names
}//rad
