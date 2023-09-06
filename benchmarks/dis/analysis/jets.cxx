#include "common_bench/benchmark.h"
#include "common_bench/mt.h"
#include "common_bench/util.h"
#include "common_bench/plot.h"

#include "ROOT/RDataFrame.hxx"
#include <TCanvas.h>
#include <TH1D.h>
#include <TVector3.h>

#include "fmt/color.h"
#include "fmt/core.h"

#include "nlohmann/json.hpp"

// Jet Benchmarks

#define ETHRESH 5.0

using ints = ROOT::VecOps::RVec<int>;
using floats = ROOT::VecOps::RVec<float>;
using vecs = ROOT::VecOps::RVec<TVector3>;

int jets(const std::string& config_name)
{
  // read our configuration
  std::ifstream  config_file{config_name};
  nlohmann::json config;
  config_file >> config;

  const std::string rec_file     = config["rec_file"];
  const std::string detector     = config["detector"];
  const std::string results_path = config["results_path"];
  const std::string plot_tag     = config["plot_tag"];
  const std::string test_tag     = config["test_tag"];
  const int ebeam                = config["ebeam"];
  const int pbeam                = config["pbeam"];

  fmt::print(fmt::emphasis::bold | fg(fmt::color::forest_green),
             "Running jet analysis...\n");
  fmt::print(" - Detector package: {}\n", detector);
  fmt::print(" - input file: {}\n", rec_file);
  fmt::print(" - test tag: {}\n", test_tag);
  fmt::print(" - ebeam: {}\n", ebeam);
  fmt::print(" - pbeam: {}\n", pbeam);

  using namespace ROOT;

  ROOT::EnableImplicitMT(kNumThreads);
  ROOT::RDataFrame d("events", rec_file);

  TFile output(fmt::format("{}/jets/{}.hist.root", results_path, plot_tag).c_str(), "RECREATE");

  // Define Lambdas for Selection
  auto getNJetsThresh = [](ints &type, floats &E){
    ints mult;
    int n=0;
    for(auto i=0U; i<type.size(); i++)
      {
	if(type[i] == 0 && E[i] >= ETHRESH) n++;
      }
    mult.push_back(n);
    return mult;
  }; // Number of Jets Above Threshold

  auto getNJets = [](ints &type, floats &E){
    ints mult;
    int n=0;
    for(auto i=0U; i<type.size(); i++)
      {
	if(type[i] == 0) n++;
      }
    mult.push_back(n);
    return mult;
  };

  auto getJetEThresh = [](ints &type, floats &E){
    floats nrg;
    for(auto i=0U; i<E.size(); i++)
      {
	if(type[i] == 0 && E[i] >= ETHRESH) nrg.push_back(E[i]);
      }
    return nrg;
  }; // Jet Energy

  auto getJetE = [](ints &type, floats &E){
    floats nrg;
    for(auto i=0U; i<E.size(); i++)
      {
	if(type[i] == 0) nrg.push_back(E[i]);
      }
    return nrg;
  }; // Jet Energy

  auto getJetVecThresh = [](ints &type, floats &E, floats &x, floats &y, floats &z){
    vecs v;
    for(auto i=0U; i<x.size(); i++)
      {
	TVector3 tmp(x[i],y[i],z[i]);
	if(type[i] == 0 && E[i] >= ETHRESH) v.push_back(tmp);
      }
    return v;
  }; // Jet 3-Momentum

  auto getJetVec = [](ints &type, floats &E, floats &x, floats &y, floats &z){
    vecs v;
    for(auto i=0U; i<x.size(); i++)
      {
	TVector3 tmp(x[i],y[i],z[i]);
	if(type[i] == 0) 
	  {
	    //cout << i << " " << type[i] << " " << x[i] << " " << y[i] << " " << z[i] << " " << tmp.PseudoRapidity() << " " << tmp.Phi() << endl;
	    v.push_back(tmp);
	  }
      }
    return v;
  }; // Jet 3-Momentum

  auto getJetTrackDRThresh = [](ints &type, ints &pid, floats &E, floats &x, floats &y, floats &z){
    floats dr;
    for(auto i=0U; i<x.size(); i++)
      {
	if(type[i] == 0 && E[i] >= ETHRESH) // Look at jet with energy above thresh
	  {
	    TVector3 jetMom(x[i],y[i],z[i]);
	    int jetID = pid[i];

	    for(auto j=0U; j<x.size(); j++) // First loop over tracks
	      {
		if(type[j] == 1 && (pid[j] == jetID)) // Look at Tracks associated with jet
		  {
		    TVector3 trackA(x[j],y[j],z[j]);

		    if(j < x.size() - 1) // Overflow protection
		      {
			for(auto k=j+1; k<x.size(); k++) // Second loop over tracks
			  {
			    if(type[k] == 1 && (pid[k] == jetID)) // Look at tracks associated with jet
			      {
				TVector3 trackB(x[k],y[k],z[k]);

				float dEta = trackA.PseudoRapidity() - trackB.PseudoRapidity();
				float dPhi = TVector2::Phi_mpi_pi(trackA.Phi() - trackB.Phi());
				float dR = std::sqrt(dEta*dEta + dPhi*dPhi);

				dr.push_back(dR);
			      } // Track B
			  } // Second Loop
		      } // Overflow
		  } // Track A
	      } // First Loop
	  } // Jet Select
      } // Jet Loop
    return dr;
  };

  auto getJetTrackDR = [](ints &type, ints &pid, floats &x, floats &y, floats &z){
    floats dr;
    for(auto i=0U; i<x.size(); i++)
      {
	if(type[i] == 0) // Look at jet
	  {
	    TVector3 jetMom(x[i],y[i],z[i]);
	    int jetID = pid[i];

	    for(auto j=0U; j<x.size(); j++) // First loop over tracks
	      {
		if(type[j] == 1 && (pid[j] == jetID)) // Look at Tracks associated with jet
		  {
		    TVector3 trackA(x[j],y[j],z[j]);

		    if(j < x.size() - 1) // Overflow protection
		      {
			for(auto k=j+1; k<x.size(); k++) // Second loop over tracks
			  {
			    if(type[k] == 1 && (pid[k] == jetID)) // Look at tracks associated with jet
			      {
				TVector3 trackB(x[k],y[k],z[k]);

				float dEta = trackA.PseudoRapidity() - trackB.PseudoRapidity();
				float dPhi = TVector2::Phi_mpi_pi(trackA.Phi() - trackB.Phi());
				float dR = std::sqrt(dEta*dEta + dPhi*dPhi);

				dr.push_back(dR);
			      } // Track B
			  } // Second Loop
		      } // Overflow
		  } // Track A
	      } // First Loop
	  } // Jet Select
      } // Jet Loop
    return dr;
  }; // Look at distance between particles in jet

  auto getJetDupTrackThresh = [](ints &type, ints &pid, floats &E, floats &x, floats &y, floats &z){
    ints dup;
    for(auto i=0U; i<x.size(); i++)
      {
	if(type[i] == 0 && E[i] >= ETHRESH) // Look at jet with energy above thresh
	  {
	    int hasDupTrack = 0;

	    TVector3 jetMom(x[i],y[i],z[i]);
	    int jetID = pid[i];

	    for(auto j=0U; j<x.size(); j++) // First loop over tracks
	      {
		if(type[j] == 1 && (pid[j] == jetID)) // Look at Tracks associated with jet
		  {
		    TVector3 trackA(x[j],y[j],z[j]);

		    if(j < x.size() - 1) // Overflow protection
		      {
			for(auto k=j+1; k<x.size(); k++) // Second loop over tracks
			  {
			    if(type[k] == 1 && (pid[k] == jetID)) // Look at tracks associated with jet
			      {
				TVector3 trackB(x[k],y[k],z[k]);

				float dEta = trackA.PseudoRapidity() - trackB.PseudoRapidity();
				float dPhi = TVector2::Phi_mpi_pi(trackA.Phi() - trackB.Phi());
				float dR = std::sqrt(dEta*dEta + dPhi*dPhi);

				if(dR < 0.05) hasDupTrack = 1;
			      } // Track B
			  } // Second Loop
		      } // Overflow
		  } // Track A
	      } // First Loop

	    dup.push_back(hasDupTrack);
	  } // Jet Select
      } // Jet Loop
    return dup;
  }; // Identify jets containing duplicate tracks

  auto getPt = [](vecs &x){
    floats pt;
    for(auto i=0U; i<x.size(); i++)
      {
	pt.push_back(x[i].Perp());
      }
    return pt;
  }; // Jet Transverse Momentum

  auto getEta = [](vecs &x){
    floats eta;
    for(auto i=0U; i<x.size(); i++)
      {
	eta.push_back(x[i].PseudoRapidity());
      }
    return eta;
  }; // Jet Psudorapidity

  auto getPhi = [](vecs &x){
    floats phi;
    for(auto i=0U; i<x.size(); i++)
      {
	phi.push_back(x[i].Phi());
      }
    return phi;
  }; // Jet Phi

  auto getMatchDeltaR = [](vecs &vecR, vecs &vecG){
    floats deltaR;
    for(auto i=0U; i<vecR.size(); i++)
      {
	float minDeltaR = 10000.0;
	for(auto j=0U; j<vecG.size(); j++)
	  {
	    double dEta = vecR[i].PseudoRapidity() - vecG[j].PseudoRapidity();
	    double dPhi = TVector2::Phi_mpi_pi(vecR[i].Phi() - vecG[j].Phi());
	    double dR = TMath::Sqrt(dEta*dEta + dPhi*dPhi);
	    
	    if(dR < minDeltaR) minDeltaR = dR;
	  }  
	deltaR.push_back(minDeltaR);
      } 
    return deltaR;
  };

  auto getMatchE = [](vecs &vecR, vecs &vecG, floats &E){
    floats nrg;
    for(auto i=0U; i<vecR.size(); i++)
      {
	float minDeltaR = 10000.0;
	float minE = -1.0;
	for(auto j=0U; j<vecG.size(); j++)
	  {
	    double dEta = vecR[i].PseudoRapidity() - vecG[j].PseudoRapidity();
	    double dPhi = TVector2::Phi_mpi_pi(vecR[i].Phi() - vecG[j].Phi());
	    double dR = TMath::Sqrt(dEta*dEta + dPhi*dPhi);
	    
	    if(dR < minDeltaR) 
	      {
		minDeltaR = dR;
		minE = E[j];
	      }
	  }  
	nrg.push_back(minE);
      }
    return nrg;
  };
  
  auto getMatchEta = [](vecs &vecR, vecs &vecG){
    floats eta;
    for(auto i=0U; i<vecR.size(); i++)
      {
	float minDeltaR = 10000.0;
	float minEta = -100.0;
	for(auto j=0U; j<vecG.size(); j++)
	  {
	    double dEta = vecR[i].PseudoRapidity() - vecG[j].PseudoRapidity();
	    double dPhi = TVector2::Phi_mpi_pi(vecR[i].Phi() - vecG[j].Phi());
	    double dR = TMath::Sqrt(dEta*dEta + dPhi*dPhi);
	    
	    if(dR < minDeltaR) 
	      {
		minDeltaR = dR;
		minEta = vecG[j].PseudoRapidity();
	      }
	  }  
	eta.push_back(minEta);
      }
    return eta;
  };

  auto getMatchPhi = [](vecs &vecR, vecs &vecG){
    floats phi;
    for(auto i=0U; i<vecR.size(); i++)
      {
	float minDeltaR = 10000.0;
	float minPhi = -100.0;
	for(auto j=0U; j<vecG.size(); j++)
	  {
	    double dEta = vecR[i].PseudoRapidity() - vecG[j].PseudoRapidity();
	    double dPhi = TVector2::Phi_mpi_pi(vecR[i].Phi() - vecG[j].Phi());
	    double dR = TMath::Sqrt(dEta*dEta + dPhi*dPhi);
	    
	    if(dR < minDeltaR) 
	      {
		minDeltaR = dR;
		minPhi = vecG[j].Phi();
	      }
	  }  
	phi.push_back(minPhi);
      }
    return phi;
  };

  auto calcRes = [](floats &x, floats &y){
    floats res;
    for(auto i=0U; i<x.size(); i++)
      {
	res.push_back((x[i]-y[i])/y[i]);
      }
    return res;
  };

  auto filterDR = [](floats &x, floats &dr){
    floats e;
    for(auto i=0U; i<x.size(); i++)
      {
	if(dr[i] < 0.25) e.push_back(x[i]);
      }
    return e;
  };

  auto filterDRI = [](ints &x, floats &dr){
    ints e;
    for(auto i=0U; i<x.size(); i++)
      {
	if(dr[i] < 0.25) e.push_back(x[i]);
      }
    return e;
  };

  auto filterDup = [](floats &x, ints &dup){
    floats e;
    for(auto i=0U; i<x.size(); i++)
      {
	if(dup[i] == 0) e.push_back(x[i]);
      }
    return e;
  };

  // Define RDataFrame Columns 
  auto d1 = d.Define("jetMultThresh",getNJetsThresh,{"ReconstructedJets.type","ReconstructedJets.energy"})
             .Define("jetEThresh",getJetEThresh,{"ReconstructedJets.type","ReconstructedJets.energy"})
             .Define("jetTrackDRThresh",getJetTrackDRThresh,{"ReconstructedJets.type","ReconstructedJets.PDG","ReconstructedJets.energy","ReconstructedJets.momentum.x","ReconstructedJets.momentum.y","ReconstructedJets.momentum.z"})
             .Define("jetDupTrackThresh",getJetDupTrackThresh,{"ReconstructedJets.type","ReconstructedJets.PDG","ReconstructedJets.energy","ReconstructedJets.momentum.x","ReconstructedJets.momentum.y","ReconstructedJets.momentum.z"})
             .Define("jetVecThresh",getJetVecThresh,{"ReconstructedJets.type","ReconstructedJets.energy","ReconstructedJets.momentum.x","ReconstructedJets.momentum.y","ReconstructedJets.momentum.z"})
             .Define("jetPtThresh",getPt,{"jetVecThresh"})
             .Define("jetEtaThresh",getEta,{"jetVecThresh"})
             .Define("jetPhiThresh",getPhi,{"jetVecThresh"})
             .Define("genJetMult",getNJets,{"GeneratedJets.type","GeneratedJets.energy"})
             .Define("genJetE",getJetE,{"GeneratedJets.type","GeneratedJets.energy"})
             .Define("genJetTrackDR",getJetTrackDR,{"GeneratedJets.type","GeneratedJets.PDG","GeneratedJets.momentum.x","GeneratedJets.momentum.y","GeneratedJets.momentum.z"})
             .Define("genJetVec",getJetVec,{"GeneratedJets.type","GeneratedJets.energy","GeneratedJets.momentum.x","GeneratedJets.momentum.y","GeneratedJets.momentum.z"})
             .Define("genJetPt",getPt,{"genJetVec"})
             .Define("genJetEta",getEta,{"genJetVec"})
             .Define("genJetPhi",getPhi,{"genJetVec"})
             .Define("matchJetDeltaR",getMatchDeltaR,{"jetVecThresh","genJetVec"})
             .Define("matchJetE",getMatchE,{"jetVecThresh","genJetVec","genJetE"})
             .Define("matchJetEta",getMatchEta,{"jetVecThresh","genJetVec"})
             .Define("matchJetPhi",getMatchPhi,{"jetVecThresh","genJetVec"})
             .Define("matchJetResE",calcRes,{"jetEThresh","matchJetE"})
             .Define("jetEThreshDR",filterDR,{"jetEThresh","matchJetDeltaR"})
             .Define("matchJetEDR",filterDR,{"matchJetE","matchJetDeltaR"}) 
             .Define("jetEtaThreshDR",filterDR,{"jetEtaThresh","matchJetDeltaR"})
             .Define("jetPhiThreshDR",filterDR,{"jetPhiThresh","matchJetDeltaR"})
             .Define("matchJetEtaDR",filterDR,{"matchJetEta","matchJetDeltaR"})
             .Define("matchJetPhiDR",filterDR,{"matchJetPhi","matchJetDeltaR"})
             .Define("matchJetResEDR",calcRes,{"jetEThreshDR","matchJetEDR"})
             .Define("jetDupTrackThreshDR",filterDRI,{"jetDupTrackThresh","matchJetDeltaR"})
             .Define("jetEThreshDRDup",filterDup,{"jetEThreshDR","jetDupTrackThreshDR"})
             .Define("matchJetEDRDup",filterDup,{"matchJetEDR","jetDupTrackThreshDR"})
             .Define("jetEtaThreshDRDup",filterDup,{"jetEtaThreshDR","jetDupTrackThreshDR"})
             .Define("jetPhiThreshDRDup",filterDup,{"jetPhiThreshDR","jetDupTrackThreshDR"})
             .Define("matchJetEtaDRDup",filterDup,{"matchJetEtaDR","jetDupTrackThreshDR"})
             .Define("matchJetPhiDRDup",filterDup,{"matchJetPhiDR","jetDupTrackThreshDR"})
             .Define("matchJetResEDRDup",calcRes,{"jetEThreshDRDup","matchJetEDRDup"});

  // Book Histograms
  auto h_type = d.Histo1D({"h_type",";Type",3, 0.,3.}, "ReconstructedJets.type");
  auto h_nJetsThresh = d1.Histo1D({"h_nJetsThresh","Number of Reco Jets (E >= 5 GeV);Num Jets",20,0.,20.}, "jetMultThresh");
  auto h_jetEnergyThresh = d1.Histo1D({"h_jetEnergyThresh","Reco Jet Energy (E >= 5 GeV);Energy [GeV]",300,0.,300.}, "jetEThresh");
  auto h_jetTrackDRThresh = d1.Histo1D({"h_jetTrackDRThresh","Reco Jet: Pairwise Distance Between Tracks (E >= 5 GeV);Delta R",1000,0.,10.}, "jetTrackDRThresh");
  auto h_jetDupTrackThresh = d1.Histo1D({"h_jetDupTrackThresh","Reco Jets with Duplicate Tracks (E >= 5 GeV);0 = No Duplicate 1 = Duplicate",5,0.,5.}, "jetDupTrackThresh");
  auto h_jetPtThresh = d1.Histo1D({"h_jetPtThresh","Reco Jet Pt (E >= 5 GeV);Jet Pt [GeV]",100,0.,50.}, "jetPtThresh");
  auto h_jetEtaThresh = d1.Histo1D({"h_jetEtaThresh","Reco Jet Eta (E >= 5 GeV);Eta",100,-5.,5.}, "jetEtaThresh");
  auto h_jetPhiThresh = d1.Histo1D({"h_jetPhiThresh","Reco Jet Phi (E >= 5 GeV);Phi [Rad]",100,-TMath::Pi(),TMath::Pi()}, "jetPhiThresh");
  auto h_jetPhiVsEtaThresh = d1.Histo2D({"h_jetPhiVsEtaThresh","Reco Jet Phi Vs Eta (E >= 5 GeV);Jet Eta;Jet Phi",100,-5.,5.,100,-TMath::Pi(),TMath::Pi()}, "jetEtaThresh", "jetPhiThresh");
  auto h_jetEVsEtaThresh = d1.Histo2D({"h_jetEVsEtaThresh","Reco Jet E Vs Eta (E >= 5 GeV);Jet Eta;Jet Energy [GeV]",100,-5.,5.,300,0.,300.}, "jetEtaThresh", "jetEThresh");

  auto h_nGenJets = d1.Histo1D({"h_nGenJets","Number of Truth Jets;Num Jets",20,0.,20.}, "genJetMult");
  auto h_genJetEnergy = d1.Histo1D({"h_genJetEnergy","Truth Jet Energy;Energy [GeV]",300,0.,300.}, "genJetE");
  auto h_genJetTrackDR = d1.Histo1D({"h_genJetTrackDR","Truth Jet: Pairwise Distance Between Charged Particles;Delta R",1000,0.,10.}, "genJetTrackDR");
  auto h_genJetPt = d1.Histo1D({"h_genJetPt","Truth Jet Pt;Jet Pt [GeV]",100,0.,50.}, {"genJetPt"});
  auto h_genJetEta = d1.Histo1D({"h_genJetEta","Truth Jet Eta;Eta",100,-5.,5.}, "genJetEta");
  auto h_genJetPhi = d1.Histo1D({"h_genJetPhi","Truth Jet Phi;Phi [Rad]",100,-TMath::Pi(),TMath::Pi()}, {"genJetPhi"});
  auto h_genJetPhiVsEta = d1.Histo2D({"h_genJetPhiVsEta","Truth Jet Phi Vs Eta;Jet Eta;Jet Phi",100,-5.,5.,100,-TMath::Pi(),TMath::Pi()}, "genJetEta", "genJetPhi");
  auto h_genJetEVsEta = d1.Histo2D({"h_genJetEVsEta","Truth Jet E Vs Eta;Jet Eta;Jet Energy [GeV]",100,-5.,5.,300,0.,300.}, "genJetEta", "genJetE");

  auto h_jetMatchDeltaR = d1.Histo1D({"h_jetMatchDeltaR","Distance Between Closest Reco-Truth Jet;Delta R",1000,0.,5.}, "matchJetDeltaR");
  auto h_jetMatchTruthVsRecoE = d1.Histo2D({"h_jetMatchTruthVsRecoE","Truth Vs Reco Jet Energy;Reco Energy [GeV];Truth Energy [GeV]",300,0.,300.,300,0.,300.}, "jetEThresh","matchJetE");
  auto h_jetMatchTruthVsRecoEta = d1.Histo2D({"h_jetMatchTruthVsRecoEta","Truth Vs Reco Jet Eta;Reco Eta;Truth Eta",100,-5.,5.,100,-5.,5.}, "jetEtaThresh","matchJetEta");
  auto h_jetMatchTruthVsRecoPhi = d1.Histo2D({"h_jetMatchTruthVsRecoPhi","Truth Vs Reco Jet Phi;Reco Phi;Truth Phi",100,-TMath::Pi(),TMath::Pi(),100,-TMath::Pi(),TMath::Pi()}, "jetPhiThresh","matchJetPhi");
  auto h_jetMatchERes = d1.Histo1D({"h_jetMatchERes","(Reco - Truth)/Truth Jet Energy;Energy Diff",200,-10.,10.}, "matchJetResE");
  auto h_jetMatchEResVsE = d1.Histo2D({"h_jetMatchEResVsE","(Reco - Truth)/Truth Jet Energy Vs Reco Energy;Energy [GeV];Energy Diff",300,0.,300,2000,-10.,10.}, "jetEThresh","matchJetResE");
  auto h_jetMatchEResVsEta = d1.Histo2D({"h_jetMatchEResVsEta","(Reco - Truth)/Truth Jet Energy Vs Reco Eta;Reco Jet Eta;Energy Diff",100,-5.,5.,2000,-10.,10.}, "jetEtaThresh", "matchJetResE");

  auto h_jetMatchTruthVsRecoEDR = d1.Histo2D({"h_jetMatchTruthVsRecoEDR","Truth Vs Reco Jet Energy (Delta R < 0.25);Reco Energy [GeV];Truth Energy [GeV]",300,0.,300.,300,0.,300.}, "jetEThreshDR","matchJetEDR");
  auto h_jetMatchTruthVsRecoEtaDR = d1.Histo2D({"h_jetMatchTruthVsRecoEtaDR","Truth Vs Reco Jet Eta (Delta R < 0.25);Reco Eta;Truth Eta",100,-5.,5.,100,-5.,5.}, "jetEtaThreshDR","matchJetEtaDR");
  auto h_jetMatchTruthVsRecoPhiDR = d1.Histo2D({"h_jetMatchTruthVsRecoPhiDR","Truth Vs Reco Jet Phi (Delta R < 0.25);Reco Phi;Truth Phi",100,-TMath::Pi(),TMath::Pi(),100,-TMath::Pi(),TMath::Pi()}, "jetPhiThreshDR","matchJetPhiDR");
  auto h_jetMatchEResDR = d1.Histo1D({"h_jetMatchEResDR","(Reco - Truth)/Truth Jet Energy (Delta R < 0.25);Energy Diff",200,-10.,10.}, "matchJetResEDR");
  auto h_jetMatchEResVsEDR = d1.Histo2D({"h_jetMatchEResVsEDR","(Reco - Truth)/Truth Jet Energy Vs Reco Energy (Delta R < 0.25);Energy [GeV];Energy Diff",300,0.,300,2000,-10.,10.}, "jetEThreshDR","matchJetResEDR");
  auto h_jetMatchEResVsEtaDR = d1.Histo2D({"h_jetMatchEResVsEtaDR","(Reco - Truth)/Truth Jet Energy Vs Reco Eta (Delta R < 0.25);Reco Jet Eta;Energy Diff",100,-5.,5.,2000,-10.,10.}, "jetEtaThreshDR", "matchJetResEDR");

  auto h_jetMatchTruthVsRecoEDRDup = d1.Histo2D({"h_jetMatchTruthVsRecoEDRDup","Truth Vs Reco Jet Energy (Delta R < 0.25 No Duplicate);Reco Energy [GeV];Truth Energy [GeV]",300,0.,300.,300,0.,300.}, "jetEThreshDRDup","matchJetEDRDup");
  auto h_jetMatchTruthVsRecoEtaDRDup = d1.Histo2D({"h_jetMatchTruthVsRecoEtaDRDup","Truth Vs Reco Jet Eta (Delta R < 0.25 No Duplicate);Reco Eta;Truth Eta",100,-5.,5.,100,-5.,5.}, "jetEtaThreshDRDup","matchJetEtaDRDup");
  auto h_jetMatchTruthVsRecoPhiDRDup = d1.Histo2D({"h_jetMatchTruthVsRecoPhiDRDup","Truth Vs Reco Jet Phi (Delta R < 0.25 No Duplicate);Reco Phi;Truth Phi",100,-TMath::Pi(),TMath::Pi(),100,-TMath::Pi(),TMath::Pi()}, "jetPhiThreshDRDup","matchJetPhiDRDup");
  auto h_jetMatchEResDRDup = d1.Histo1D({"h_jetMatchEResDRDup","(Reco - Truth)/Truth Jet Energy (Delta R < 0.25 No Duplicate);Energy Diff",200,-10.,10.}, "matchJetResEDRDup");
  auto h_jetMatchEResVsEDRDup = d1.Histo2D({"h_jetMatchEResVsEDRDup","(Reco - Truth)/Truth Jet Energy Vs Reco Energy (Delta R < 0.25 No Duplicate);Energy [GeV];Energy Diff",300,0.,300,2000,-10.,10.}, "jetEThreshDRDup","matchJetResEDRDup");
  auto h_jetMatchEResVsEtaDRDup = d1.Histo2D({"h_jetMatchEResVsEtaDRDup","(Reco - Truth)/Truth Jet Energy Vs Reco Eta (Delta R < 0.25 No Duplicate);Reco Jet Eta;Energy Diff",100,-5.,5.,2000,-10.,10.}, "jetEtaThreshDRDup", "matchJetResEDRDup");

  // Render histograms
  { TCanvas c; h_type->Draw(); c.Print(fmt::format("{}/jets/{}_type.png", results_path, plot_tag).c_str()); }
  { TCanvas c; h_nJetsThresh->Draw(); c.Print(fmt::format("{}/jets/{}_nJetsThresh.png", results_path, plot_tag).c_str()); }
  { TCanvas c; h_jetEnergyThresh->Draw(); c.Print(fmt::format("{}/jets/{}_jetEnergyThresh.png", results_path, plot_tag).c_str()); }
  { TCanvas c; h_jetTrackDRThresh->Draw(); c.Print(fmt::format("{}/jets/{}_jetTrackDRThresh.png", results_path, plot_tag).c_str()); }
  { TCanvas c; h_jetDupTrackThresh->Draw(); c.Print(fmt::format("{}/jets/{}_jetDupTrackThresh.png", results_path, plot_tag).c_str()); }
  { TCanvas c; h_jetPtThresh->Draw(); c.Print(fmt::format("{}/jets/{}_jetPtThresh.png", results_path, plot_tag).c_str()); }
  { TCanvas c; h_jetEtaThresh->Draw(); c.Print(fmt::format("{}/jets/{}_jetEtaThresh.png", results_path, plot_tag).c_str()); }
  { TCanvas c; h_jetPhiThresh->Draw(); c.Print(fmt::format("{}/jets/{}_jetPhiThresh.png", results_path, plot_tag).c_str()); }
  { TCanvas c; h_jetPhiVsEtaThresh->Draw(); c.Print(fmt::format("{}/jets/{}_jetPhiVsEtaThresh.png", results_path, plot_tag).c_str()); }
  { TCanvas c; h_jetEVsEtaThresh->Draw(); c.Print(fmt::format("{}/jets/{}_jetEVsEtaThresh.png", results_path, plot_tag).c_str()); }

  { TCanvas c; h_nGenJets->Draw(); c.Print(fmt::format("{}/jets/{}_nGenJets.png", results_path, plot_tag).c_str()); }
  { TCanvas c; h_genJetEnergy->Draw(); c.Print(fmt::format("{}/jets/{}_genJetEnergy.png", results_path, plot_tag).c_str()); }
  { TCanvas c; h_genJetTrackDR->Draw(); c.Print(fmt::format("{}/jets/{}_genJetTrackDR.png", results_path, plot_tag).c_str()); }
  { TCanvas c; h_genJetPt->Draw(); c.Print(fmt::format("{}/jets/{}_genJetPt.png", results_path, plot_tag).c_str()); }
  { TCanvas c; h_genJetEta->Draw(); c.Print(fmt::format("{}/jets/{}_genJetEta.png", results_path, plot_tag).c_str()); }
  { TCanvas c; h_genJetPhi->Draw(); c.Print(fmt::format("{}/jets/{}_genJetPhi.png", results_path, plot_tag).c_str()); }
  { TCanvas c; h_genJetPhiVsEta->Draw(); c.Print(fmt::format("{}/jets/{}_genJetPhiVsEta.png", results_path, plot_tag).c_str()); }
  { TCanvas c; h_genJetEVsEta->Draw(); c.Print(fmt::format("{}/jets/{}_genJetEVsEta.png", results_path, plot_tag).c_str()); }

  { TCanvas c; h_jetMatchDeltaR->Draw(); c.Print(fmt::format("{}/jets/{}_jetMatchDeltaR.png", results_path, plot_tag).c_str()); }
  { TCanvas c; h_jetMatchTruthVsRecoE->Draw(); c.Print(fmt::format("{}/jets/{}_jetMatchTruthVsRecoE.png", results_path, plot_tag).c_str()); }
  { TCanvas c; h_jetMatchTruthVsRecoEta->Draw(); c.Print(fmt::format("{}/jets/{}_jetMatchTruthVsRecoEta.png", results_path, plot_tag).c_str()); }
  { TCanvas c; h_jetMatchTruthVsRecoPhi->Draw(); c.Print(fmt::format("{}/jets/{}_jetMatchTruthVsRecoPhi.png", results_path, plot_tag).c_str()); }
  { TCanvas c; h_jetMatchERes->Draw(); c.Print(fmt::format("{}/jets/{}_jetMatchERes.png", results_path, plot_tag).c_str()); }
  { TCanvas c; h_jetMatchEResVsE->Draw(); c.Print(fmt::format("{}/jets/{}_jetMatchEResVsE.png", results_path, plot_tag).c_str()); }
  { TCanvas c; h_jetMatchEResVsEta->Draw(); c.Print(fmt::format("{}/jets/{}_jetMatchEResVsEta.png", results_path, plot_tag).c_str()); }

  { TCanvas c; h_jetMatchTruthVsRecoEDR->Draw(); c.Print(fmt::format("{}/jets/{}_jetMatchTruthVsRecoEDR.png", results_path, plot_tag).c_str()); }
  { TCanvas c; h_jetMatchTruthVsRecoEtaDR->Draw(); c.Print(fmt::format("{}/jets/{}_jetMatchTruthVsRecoEtaDR.png", results_path, plot_tag).c_str()); }
  { TCanvas c; h_jetMatchTruthVsRecoPhiDR->Draw(); c.Print(fmt::format("{}/jets/{}_jetMatchTruthVsRecoPhiDR.png", results_path, plot_tag).c_str()); }
  { TCanvas c; h_jetMatchEResDR->Draw(); c.Print(fmt::format("{}/jets/{}_jetMatchEResDR.png", results_path, plot_tag).c_str()); }
  { TCanvas c; h_jetMatchEResVsEDR->Draw(); c.Print(fmt::format("{}/jets/{}_jetMatchEResVsEDR.png", results_path, plot_tag).c_str()); }
  { TCanvas c; h_jetMatchEResVsEtaDR->Draw(); c.Print(fmt::format("{}/jets/{}_jetMatchEResVsEtaDR.png", results_path, plot_tag).c_str()); }

  { TCanvas c; h_jetMatchTruthVsRecoEDRDup->Draw(); c.Print(fmt::format("{}/jets/{}_jetMatchTruthVsRecoEDRDup.png", results_path, plot_tag).c_str()); }
  { TCanvas c; h_jetMatchTruthVsRecoEtaDRDup->Draw(); c.Print(fmt::format("{}/jets/{}_jetMatchTruthVsRecoEtaDRDup.png", results_path, plot_tag).c_str()); }
  { TCanvas c; h_jetMatchTruthVsRecoPhiDRDup->Draw(); c.Print(fmt::format("{}/jets/{}_jetMatchTruthVsRecoPhiDRDup.png", results_path, plot_tag).c_str()); }
  { TCanvas c; h_jetMatchEResDRDup->Draw(); c.Print(fmt::format("{}/jets/{}_jetMatchEResDRDup.png", results_path, plot_tag).c_str()); }
  { TCanvas c; h_jetMatchEResVsEDRDup->Draw(); c.Print(fmt::format("{}/jets/{}_jetMatchEResVsEDRDup.png", results_path, plot_tag).c_str()); }
  { TCanvas c; h_jetMatchEResVsEtaDRDup->Draw(); c.Print(fmt::format("{}/jets/{}_jetMatchEResVsEtaDRDup.png", results_path, plot_tag).c_str()); }

  // Write Histograms
  h_type->Write();
  h_nJetsThresh->Write();
  h_jetEnergyThresh->Write();
  h_jetTrackDRThresh->Write();
  h_jetDupTrackThresh->Write();
  h_jetPtThresh->Write();
  h_jetEtaThresh->Write();
  h_jetPhiThresh->Write();
  h_jetPhiVsEtaThresh->Write();
  h_jetEVsEtaThresh->Write();

  h_nGenJets->Write();
  h_genJetEnergy->Write();
  h_genJetTrackDR->Write();
  h_genJetPt->Write();
  h_genJetEta->Write();
  h_genJetPhi->Write();
  h_genJetPhiVsEta->Write();
  h_genJetEVsEta->Write();

  h_jetMatchDeltaR->Write();
  h_jetMatchTruthVsRecoE->Write();
  h_jetMatchTruthVsRecoEta->Write();
  h_jetMatchTruthVsRecoPhi->Write();
  h_jetMatchERes->Write();
  h_jetMatchEResVsE->Write();
  h_jetMatchEResVsEta->Write();

  h_jetMatchTruthVsRecoEDR->Write();
  h_jetMatchTruthVsRecoEtaDR->Write();
  h_jetMatchTruthVsRecoPhiDR->Write();
  h_jetMatchEResDR->Write();
  h_jetMatchEResVsEDR->Write();
  h_jetMatchEResVsEtaDR->Write();

  h_jetMatchTruthVsRecoEDRDup->Write();
  h_jetMatchTruthVsRecoEtaDRDup->Write();
  h_jetMatchTruthVsRecoPhiDRDup->Write();
  h_jetMatchEResDRDup->Write();
  h_jetMatchEResVsEDRDup->Write();
  h_jetMatchEResVsEtaDRDup->Write();

  // Write output and close file
  output.Write();
  output.Close();

  return 0;
}
