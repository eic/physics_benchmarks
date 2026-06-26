#include "ePICReaction.h"
#include "ReactionChannel.h"
#include "ReactionBenchmarks.h"
#include "BasicKinematicsRDF.h"
#include "ReactionKinematicsRDF.h"
#include "ElectronScatterKinematicsRDF.h"
#include <TBenchmark.h>
#include "nlohmann/json.hpp"

int ePICExcTagg2Pi(const std::string& config_name){

 // read our configuration
  std::ifstream  config_file{config_name};
  nlohmann::json config;
  config_file >> config;

  const std::string campaign = config["campaign"];
  const std::string data_files = config["rec_file"];
  const std::string ebeam = config["ebeam"];
  const std::string pbeam = config["pbeam"];
  const std::string out_dir = std::string(config["output_dir"]) + '/';

  std::cout<<"ePICPythia analysing files : "<<data_files<<std::endl;
  gROOT->ProcessLine(Form(".!ls %s",data_files.data()));
  
  ROOT::EnableImplicitMT(4);
  //increase resolution of .png files
  gStyle->SetImageScaling(4.0);
  
  gStyle->SetOptStat(0);      // Hide statistics box by default
  gStyle->SetOptTitle(0);     // Hide individual histogram titles by default
  gStyle->SetHistLineWidth(2);
  gStyle->SetFuncWidth(2);
  gStyle->SetLabelSize(0.04, "x");
  gStyle->SetLabelSize(0.05, "y");
  gStyle->SetTitleSize(0.05, "x"); // Axis titles
  gStyle->SetTitleSize(0.06, "y");
  gStyle->SetTitleOffset(0.8, "y");
  gStyle->SetPadBottomMargin(0.12); // Default for sub-pads, can be overridden per sub-pad if needed
  gStyle->SetPadTopMargin(0.02); // Default for sub-pads, can be overridden per sub-pad if needed
  gStyle->SetPadLeftMargin(0.1);
  gStyle->SetPadRightMargin(0.02); // Smaller right margin for cleaner look if no Y-axis on right
  
  // Ensure text aligns nicely, using Helvetica for clarity
  gStyle->SetTextFont(42);
  gStyle->SetTextAlign(11); // Left align  gBenchmark->Start("df");
  
  std::vector<std::vector<int>> mesons={{211,-211}};//pi+ pi-
  std::vector<std::vector<int>> baryons={{}};//and anything else, baryon will be calculated as missing particle after pi+,pi-

 
  if(mesons.size()!=baryons.size()){ cout<<"mesons "<<mesons.size()<<" "<<baryons.size()<<endl;exit(0);}
  
  std::vector<ReactionBenchmarks> allExcHists;
 
  cout<<"ePICPythia running on data "<<data_files<<endl;
  
  rad::config::ePICReaction epic{"events",data_files};
  
   //Take the beam energy and angle from the MCParticles branch
   //Here we actually take the mean over all events,
   //as this is likely to be what we have in the experiment
   epic.SetBeamsFromMC();
   auto electronP = epic.P4BeamEle().P();
   auto protonP = epic.P4BeamIon().P();
   auto WMax = (epic.P4BeamEle()+epic.P4BeamIon()).M();
   
   //Define the core columns we are interested in and match
   //Reconstructed particles with their Generated 4-momentum
   //This creates rec_ and tru_ columns which are synched to each other
   epic.AliasColumnsAndMatchWithMC();

  //plot arrangements
  std::vector<std::pair<uint,uint>> arrange;


  for(UInt_t ipy = 0; ipy<mesons.size() ; ++ipy){
    //copy from base dataframe to allow lazy execution
    //on all final states
    auto rad = epic;

    //create a reaction channel which is defined in terms of the
    //meson and baryon decay products
    //truth pid is used for filtering reactions.
    rad::config::ReactionChannel channel{rad,mesons[ipy],baryons[ipy],"tru_pid",rad::config::DoElectroReaction};  

    //Define all kinematics we are interested in
    rad::rdf::MissMass(rad,"W","{scat_ele}");
    rad::rdf::MissMass(rad,"MissMassMeson","{scat_ele,idxMeson}");
    rad::rdf::MissMass2(rad,"MissMassMeson2","{scat_ele,idxMeson}");
    rad::rdf::MissMass2(rad,"MissMassBaryon2","{scat_ele,idxBaryon}");
    rad::rdf::MissMass2(rad,"MissMass2","{scat_ele,idxBaryon,idxMeson}");
    rad::rdf::Mass(rad,"Whad","{idxMeson,idxBaryon}");
    rad::rdf::Mass(rad,"MesonMass","{idxMeson}");
    rad::rdf::Mass(rad,"BaryonMass","{idxBaryon}");

    //missing momenta
    rad::rdf::MissP(rad,"MissPMeson","{scat_ele,idxMeson}");
    rad::rdf::MissPt(rad,"MissPtMeson","{scat_ele,idxMeson}");
    rad::rdf::MissPz(rad,"MissPzMeson","{scat_ele,idxMeson}");
    rad::rdf::MissTheta(rad,"MissThetaMeson","{scat_ele,idxMeson}");

    //t distribution, column name
    rad::rdf::TBot(rad,"tb");
    rad::rdf::TPrimeBot(rad,"tbp_pn");
    rad::rdf::TTop(rad,"tt");
    rad::rdf::TPrimeTop(rad,"ttp_pn");

    //CM production angles
    rad::rdf::CMAngles(rad,"CM");
    rad::rdf::Q2(rad,"Q2");


    ///////////////////////////////////////////////////////////
    //Define benchmarks
    ///////////////////////////////////////////////////////////

    ReactionBenchmarks excTaggHists{"ExcTagger",mesons[ipy],baryons[ipy]};
    excTaggHists.SetOutDir(out_dir);
    
    //define canvas formatting as 3 by 4 plots
    arrange.push_back({4,3});
 
    //Define each variable we wish to benchmark with histogram model
    excTaggHists.AddVar("Q2",{"Q2","Q2",200,0,0.05},{200,-0.0005,0.0005});
    excTaggHists.AddVar("W",{"W","W",200,0,1.1*WMax});
    excTaggHists.AddVar("Whad",{"Whad","Whad",200,0,1.1*WMax});
    excTaggHists.AddVar("MissMass2",{"MissMass2","Missing mass squared exclusive",200,-1.1*WMax,1.1*WMax});
    excTaggHists.AddVar("MissMassMeson2",{"MissMassMeson2","Missing mass squared from meson",1000,-WMax*WMax,WMax*WMax});
    excTaggHists.AddVar("MissMassBaryon2",{"MissMassBaryon2","Missing mass squared from baryon",1000,-0.1*WMax,0.1*WMax});
    excTaggHists.AddVar("MesonMass",{"MesonMass","M(meson) [GeV]",200,0.,4});
    excTaggHists.AddVar("BaryonMass",{"BaryonMass","M(baryon) [GeV]",200,-50,50},{100,-100,100});
    excTaggHists.AddVar("tt",{"tt","-t(top) [GeV^{2}]",100,-1,5});
    excTaggHists.AddVar("tb",{"tb","-t(bottom) [GeV^{2}]",100,-1,5});
    excTaggHists.AddVar("CM_CosTheta",{"cthCM","cos(#theta_{CM})",200,0.99,1},{100,-0.0001,0.0001});
    excTaggHists.AddVar("CM_Phi",{"phCM","#phi_{CM}",100,-1.1*TMath::Pi(),1.1*TMath::Pi()},{100,-1,1});
    //some exclusivity variables
    arrange.push_back({2,2});
    excTaggHists.AddVar("MissPMeson",{"MissPMeson","Missing momentum from meson",200,0.8*protonP,protonP*1.2});
    excTaggHists.AddVar("MissPtMeson",{"MissPtMeson","Missing transvers momentum from meson",200,0,2});
    excTaggHists.AddVar("MissPzMeson",{"MissPMeson","Missing longitudinal momentum from meson",200,0.8*protonP,protonP*1.2});
    excTaggHists.AddVar("MissThetaMeson",{"MissThetaMeson","Missing theta from meson",400,0,0.02});
 
    //some individual particles : scattered electron momentum
    arrange.push_back({3,1});
    excTaggHists.AddVarElement("scat_ele","pmag",{"epmag","Scattered Electron momentum",200,-1,1.2*electronP},{100,-0.1,0.1});
    excTaggHists.AddVarElement("scat_ele","theta",{"etheta","Scattered Electron #theta",300,3.1,TMath::Pi()+0.01});
    excTaggHists.AddVarElement("scat_ele","phi",{"ephi","Scattered Electron #phi",200,-2*TMath::Pi(),2*TMath::Pi()});
    
    //add exclusive tagger events
    auto rad_tagger=rad;
    //Now filter on exclusive truth final state, i.e. recoil baryon == nucleon
    rad_tagger.Filter("TMath::Abs(tru_MissMassMeson-0.94)<0.05","exc_cut");
    excTaggHists.Declare(rad_tagger,channel.CutParticleCondition("rec_pmag",">0.05")+"&&rec_pmag[scat_ele]>0.1&&rec_theta[scat_ele]>3.1");//tagger >3.1

    allExcHists.push_back(std::move(excTaggHists));
    
 
  }
  
  ///////////////////////////////////////////////////////////
  //Draw histograms and calculate benchmarks
  ///////////////////////////////////////////////////////////
  gBenchmark->Start("processing");


  for(auto& hists:allExcHists)
    hists.Finalise(arrange);
   
 
  gBenchmark->Stop("processing");
  gBenchmark->Print("processing");
  gBenchmark->Stop("df");
  gBenchmark->Print("df");
  
  return 0;
  
}
