#include "common_bench/benchmark.h"
#include "common_bench/mt.h"
#include "common_bench/util.h"
#include "common_bench/plot.h"

#include <TCanvas.h>
#include <TChain.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TLegend.h>
#include <TVector3.h>

#include "fmt/color.h"
#include "fmt/core.h"

#include "nlohmann/json.hpp"

double StudentT(double *x, double *par);

// Vertexing Benchmarks
// Author: K. Singla

// To run: snakemake -c1 results/epic_craterlake/dis/10on100/minQ2=1/vertexing

int vertexing(const std::string& config_name)
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
            "Running vertexing analysis...\n");
  fmt::print(" - Detector package: {}\n", detector);
  fmt::print(" - input file: {}\n", rec_file);
  fmt::print(" - test tag: {}\n", test_tag);
  fmt::print(" - ebeam: {}\n", ebeam);
  fmt::print(" - pbeam: {}\n", pbeam);

  const bool PRINT = true;

  // Input
  auto *tree = new TChain("events");
  tree->Add(rec_file.c_str());

  const int seabornRed = TColor::GetColor(213, 94, 0);

// Seaborn Green: #009E73 -> (0, 158, 115)
const int seabornGreen = TColor::GetColor(0, 158, 115);

// Seaborn Blue: #56B4E9 -> (86, 180, 233)
const int seabornBlue = TColor::GetColor(100, 149, 237);

  // Output
  TFile *ofile = TFile::Open("test.root","RECREATE");

  // TTreeReader
  TTreeReader tree_reader(tree);
  
  // Reco Vertex
  TTreeReaderArray<int> recoVtxType = {tree_reader, "CentralTrackVertices.type"};
  TTreeReaderArray<float> recoVtxX = {tree_reader, "CentralTrackVertices.position.x"};
  TTreeReaderArray<float> recoVtxY = {tree_reader, "CentralTrackVertices.position.y"};
  TTreeReaderArray<float> recoVtxZ = {tree_reader, "CentralTrackVertices.position.z"};
  
  TTreeReaderArray<unsigned int> assoPartBegin = {tree_reader, "CentralTrackVertices.associatedParticles_begin"};
  TTreeReaderArray<unsigned int> assoPartEnd = {tree_reader, "CentralTrackVertices.associatedParticles_end"};
  TTreeReaderArray<int> assoPartIndex = {tree_reader, "_CentralTrackVertices_associatedParticles.index"};

  // MC
  TTreeReaderArray<int> mcGenStat = {tree_reader, "MCParticles.generatorStatus"};
  TTreeReaderArray<int> mcPDG = {tree_reader, "MCParticles.PDG"};
  TTreeReaderArray<float> mcCharge = {tree_reader, "MCParticles.charge"};
  TTreeReaderArray<float> mcMomX = {tree_reader, "MCParticles.momentum.x"};
  TTreeReaderArray<float> mcMomY = {tree_reader, "MCParticles.momentum.y"};
  TTreeReaderArray<float> mcMomZ = {tree_reader, "MCParticles.momentum.z"};

  TTreeReaderArray<double> mcVtxX = {tree_reader, "MCParticles.vertex.x"};
  TTreeReaderArray<double> mcVtxY = {tree_reader, "MCParticles.vertex.y"};
  TTreeReaderArray<double> mcVtxZ = {tree_reader, "MCParticles.vertex.z"};

  TTreeReaderArray<unsigned int> mcParentBegin = {tree_reader, "MCParticles.parents_begin"};
  TTreeReaderArray<unsigned int> mcParentEnd = {tree_reader, "MCParticles.parents_end"};
  TTreeReaderArray<int> mcParentIndex = {tree_reader, "_MCParticles_parents.index"};

  // Reco
  TTreeReaderArray<int> recoType = {tree_reader, "ReconstructedChargedParticles.type"};

	
  // Define Histograms
  TH1D *counter = new TH1D("counter","",10,0.,10.);
  
  // Reco
  TH1D *numRecoTracksHist = new TH1D("numRecoTracks","",31,-0.5,30.5);
  TH2D *recoVsMCTracksHist = new TH2D("recoVsMCTracks","",31,-0.5,30.5,31,-0.5,30.5);
  TH1D *recoVtxEffHist = new TH1D("recoVtxEff","",4,-0.5,3.5);
  TH2D *recoVtxYvsXHist = new TH2D("recoVtxYvsX","",200,-1.,1.,200,-1,1);
  TH2D *recoVtxRvsZHist = new TH2D("recoVtxRvsZ","",200,-100.,100.,100,0.,0.8);

  // Gen
  TH1D *numGenTracksHist = new TH1D("numGenTracks","",31,-0.5,30.5);
  TH2D *genVtxYvsXHist = new TH2D("genVtxYvsXHist","",200,-1.,1.,200,-1,1);
  TH2D *genVtxRvsZHist = new TH2D("genVtxRvsZHist","",200,-100.,100.,100,0,0.5);
 
  // Resolution
  TH2D *vtxResXvsGenTrkHist = new TH2D("vtxResXvsGenTrk","",31,-0.5,30.5,200,-1,1);
  TH2D *vtxResYvsGenTrkHist = new TH2D("vtxResYvsGenTrk","",31,-0.5,30.5,200,-1,1);
  TH2D *vtxResZvsGenTrkHist = new TH2D("vtxResZvsGenTrk","",31,-0.5,30.5,200,-1,1);
  TH2D *vtxResXvsRecoTrkHist = new TH2D("vtxResXvsRecoTrk","",31,-0.5,30.5,200,-1,1);
  TH2D *vtxResYvsRecoTrkHist = new TH2D("vtxResYvsRecoTrk","",31,-0.5,30.5,200,-1,1);
  TH2D *vtxResZvsRecoTrkHist = new TH2D("vtxResZvsRecoTrk","",31,-0.5,30.5,200,-1,1);
  
  TH1D *numGenTrkswithVtxHist = new TH1D("numGenTrkswithVtx","",31,-0.5,30.5);
  TH1D *vtxEffVsGenTrkHist = new TH1D("vtxEffVsGenTrk","",31,-0.5,30.5);
  TH1D *numRecoTrkswithVtxHist = new TH1D("numRecoTrkswithVtx","",31,-0.5,30.5);
  TH1D *vtxEffVsRecoTrkHist = new TH1D("vtxEffVsRecoTrk","",31,-0.5,30.5);

  // Loop Through Events
  int NEVENTS = 0;
  while(tree_reader.Next()) {

    if(NEVENTS%10000 == 0) 
    cout << "\nEvents Processed: " << NEVENTS << endl;

    counter->Fill(0);

    //////////////////////////////////////////////////////////////////////////
    ////////////////////////  Analyze MC Tracks  /////////////////////////
    //////////////////////////////////////////////////////////////////////////
    
    //Finding MC vertex using scattered electron
    TVector3 mcEvtVtx(-999., -999., -999.);    
    for(unsigned int i=0; i<mcGenStat.GetSize(); i++)
    {
	if(mcGenStat[i] != 1) continue;
	if(mcPDG[i] != 11) continue;

	bool scatEfound = false;
	// mcParentBegin and mcParentEnd specify the entries from _MCParticles_parents.index 
	// _MCParticles_parents.index stores the MCParticle index
	for(unsigned int j=mcParentBegin[i]; j<mcParentEnd[i]; j++)
	{
          int parentPDG = mcPDG[mcParentIndex[j]];
          if(parentPDG == 11) scatEfound = true;
        }
       
	if(scatEfound == false) continue;
	//Scattered electron found
	double vtx_mc_x =  mcVtxX[i];
	double vtx_mc_y =  mcVtxY[i];
	double vtx_mc_z =  mcVtxZ[i];
	mcEvtVtx = TVector3(vtx_mc_x, vtx_mc_y, vtx_mc_z);
    }
    genVtxYvsXHist->Fill(mcEvtVtx.x(), mcEvtVtx.y());
    TVector3 mcRadius(mcEvtVtx.x(), mcEvtVtx.y(), 0);
    genVtxRvsZHist->Fill(mcEvtVtx.z(), mcRadius.Mag());
    
    //Filtering MC Tracks
    int numMCTracks=0;
    for(unsigned int i=0; i<mcGenStat.GetSize(); i++)
    {
	if(mcGenStat[i] != 1) continue;
	if(mcCharge[i] == 0) continue;
	
	TVector3 mcPartVtx(mcVtxX[i], mcVtxY[i], mcVtxZ[i]);
	TVector3 vtx_diff = mcPartVtx - mcEvtVtx;
	if(vtx_diff.Mag() > 1e-4) continue;
	
	TVector3 mcPartMom(mcMomX[i], mcMomY[i], mcMomZ[i]);
	if(fabs(mcPartMom.Eta()) > 3.5) continue;
	
	numMCTracks++;
    }
    numGenTracksHist->Fill(numMCTracks);

    //////////////////////////////////////////////////////////////////////////
    //////////////////////  Analyze Reconstructed Tracks  //////////////////////
    //////////////////////////////////////////////////////////////////////////

    numRecoTracksHist->Fill(recoType.GetSize());
	  
    //Finding Reconstructed Vertex and Vertexing Efficiency
    int nVtx=0;
    float diff=999.;
    int nAssoPart=0;
    TVector3 recoEvtVtx(-999., -999., -999.);    
    for(unsigned int i=0; i<recoVtxType.GetSize(); i++)
    {
	nVtx++;
	
	TVector3 recoVtx(recoVtxX[i], recoVtxY[i], recoVtxZ[i]);
	
	//Finding the reconstructed vertex closer to the MC vertex
	TVector3 vtx_diff = recoVtx - mcEvtVtx;
	if(vtx_diff.Mag() < diff)
	{
	    diff = vtx_diff.Mag();
	    recoEvtVtx = recoVtx;
	    
	    for(unsigned int j=assoPartBegin[i]; j<assoPartEnd[i]; j++)
	    {
              nAssoPart = j;
            }
	}
    }
    
    recoVtxEffHist->Fill(nVtx);
    recoVtxYvsXHist->Fill(recoEvtVtx.x(), recoEvtVtx.y());
    TVector3 recoRadius(recoEvtVtx.x(), recoEvtVtx.y(), 0);
    recoVtxRvsZHist->Fill(recoEvtVtx.z(), recoRadius.Mag());
    
    vtxResXvsGenTrkHist->Fill(numMCTracks, recoEvtVtx.x() - mcEvtVtx.x());
    vtxResYvsGenTrkHist->Fill(numMCTracks, recoEvtVtx.y() - mcEvtVtx.y());
    vtxResZvsGenTrkHist->Fill(numMCTracks, recoEvtVtx.z() - mcEvtVtx.z());
    
    vtxResXvsRecoTrkHist->Fill(nAssoPart, recoEvtVtx.x() - mcEvtVtx.x());
    vtxResYvsRecoTrkHist->Fill(nAssoPart, recoEvtVtx.y() - mcEvtVtx.y());
    vtxResZvsRecoTrkHist->Fill(nAssoPart, recoEvtVtx.z() - mcEvtVtx.z());
    
    if(nVtx !=0) {
    numGenTrkswithVtxHist->Fill(numMCTracks);
    numRecoTrkswithVtxHist->Fill(recoType.GetSize());
    
    recoVsMCTracksHist->Fill(numMCTracks, nAssoPart);} 
	  
    NEVENTS++;
  }
  
  gStyle->SetOptStat(0);
  ////////////////////////  MC Plots  ////////////////////////
  
  // MC vs RC associated particles
  TCanvas *c1 = new TCanvas("c1","MC vs Reco Tracks",800,600);
  c1->Clear();
  c1->Divide(1,1);

  auto func = new TF1("func","x",-1,40);
  c1->cd(1);
  recoVsMCTracksHist->Draw("COLZ");
  func->SetLineStyle(2);
  func->SetLineColor(1);
  func->Draw("SAMEL");
  recoVsMCTracksHist->SetTitle("Number of associated particles with vertex");
  recoVsMCTracksHist->GetXaxis()->SetTitle("N_{MC}");
  recoVsMCTracksHist->GetYaxis()->SetTitle("N_{RC}");


  if(PRINT) c1->Print((results_path+"/vertexing/numberAssoPart.png").c_str()); // MC versus RC Tracks
  delete c1;
 
 // MC vx versus vy Vertex
  TCanvas *c2 = new TCanvas("c2","MC Vertices",800,600);
  c2->Clear();
  c2->Divide(1,1);

  c2->cd(1);
  genVtxYvsXHist->Draw("COLZ");
  genVtxYvsXHist->SetTitle("MC Vertex: v_{x} versus v_{y}");
  genVtxYvsXHist->GetXaxis()->SetTitle("x-coordinate (in mm)");
  genVtxYvsXHist->GetYaxis()->SetTitle("y-coordinate (in mm)");
  
  if(PRINT) c2->Print((results_path+"/vertexing/MCvertexYvsX.png").c_str()); // MC Vertex vx versus vy
  delete c2;
  
  // MC vr versus vz Vertex
  TCanvas *c3 = new TCanvas("c3","MC Vertices",800,600);
  c3->Clear();
  c3->Divide(1,1);
  
  c3->cd(1);
  genVtxRvsZHist->Draw("COLZ");
  genVtxRvsZHist->SetTitle("MC Vertex: v_{r} versus v_{z}");
  genVtxRvsZHist->GetXaxis()->SetTitle("z-coordinate (in mm)");
  genVtxRvsZHist->GetYaxis()->SetTitle("#sqrt{x^{2} + y^{2}} (in mm)");


  if(PRINT) c3->Print((results_path+"/vertexing/MCvertexRvsZ.png").c_str()); // MC Vertex vr versus vz
  delete c3;
  
  ////////////////////////  Reconstructed Plots  ////////////////////////
  
  //Vertexing Efficiency
  TCanvas *c4 = new TCanvas("c4","Vertexing Efficiency",800,600);
  c4->Clear();
  c4->Divide(1,1);

  c4->cd(1);
  recoVtxEffHist->Scale(100./NEVENTS);
  recoVtxEffHist->Draw("P");
  recoVtxEffHist->SetLineColor(seabornRed);
  recoVtxEffHist->SetMarkerColor(seabornRed);
  recoVtxEffHist->SetMarkerStyle(8);
  recoVtxEffHist->SetMarkerSize(1.2);
  recoVtxEffHist->SetTitle("Vertexing Efficiency");
  recoVtxEffHist->GetXaxis()->SetTitle("# of Vertex");
  recoVtxEffHist->GetYaxis()->SetTitle("nEvents/total_events %");

  if(PRINT) c4->Print((results_path+"/vertexing/vtxEfficiency.png").c_str()); // Vertexing Efficiency
  delete c4;
  
  //Vertexing Efficiency vs MC Tracks
  TCanvas *c5 = new TCanvas("c5","Vertexing Efficiency vs MC Tracks",800,600);
  c5->Clear();
  c5->Divide(1,1);

  c5->cd(1);
  for(int i=0; i<=numGenTracksHist->GetNbinsX(); i++)
	{
		float neventsMC = numGenTracksHist->GetBinContent(i);
		float nvtxevtsMC = numGenTrkswithVtxHist->GetBinContent(i);
		
		if(neventsMC != 0)
		{
			vtxEffVsGenTrkHist->SetBinContent(i, nvtxevtsMC/neventsMC);
			vtxEffVsGenTrkHist->SetBinError(i, sqrt((nvtxevtsMC+1)/(neventsMC+2)*((nvtxevtsMC+2)/(neventsMC+3)-(nvtxevtsMC+1)/(neventsMC+2))));
		}
	}

  vtxEffVsGenTrkHist->Draw("P");
  vtxEffVsGenTrkHist->SetMarkerColor(seabornRed);
  vtxEffVsGenTrkHist->SetMarkerStyle(8);
  vtxEffVsGenTrkHist->SetMarkerSize(1.2);
  vtxEffVsGenTrkHist->SetTitle("Vertexing Efficiency vs MC Tracks");
  vtxEffVsGenTrkHist->GetXaxis()->SetTitle("# of MC Tracks");
  vtxEffVsGenTrkHist->GetYaxis()->SetTitle("Vertexing Efficiency");
  vtxEffVsGenTrkHist->GetYaxis()->SetRangeUser(0, 1.2);

  if(PRINT) c5->Print((results_path+"/vertexing/vtxEffVsMCTrks.png").c_str()); // Vertexing Efficiency versus MC Tracks
  delete c5;

  //Vertexing Efficiency vs RC Tracks
  TCanvas *c6 = new TCanvas("c6","Vertexing Efficiency vs RC Tracks",800,600);
  c6->Clear();
  c6->Divide(1,1);

  c6->cd(1);
  for(int i=0; i<=numRecoTracksHist->GetNbinsX(); i++)
	{
		float neventsRC = numRecoTracksHist->GetBinContent(i);
		float nvtxevtsRC = numRecoTrkswithVtxHist->GetBinContent(i);
		
		if(neventsRC != 0)
		{
			vtxEffVsRecoTrkHist->SetBinContent(i, nvtxevtsRC/neventsRC);
			vtxEffVsRecoTrkHist->SetBinError(i, sqrt((nvtxevtsRC+1)/(neventsRC+2)*((nvtxevtsRC+2)/(neventsRC+3)-(nvtxevtsRC+1)/(neventsRC+2))));
		}
	}

  vtxEffVsRecoTrkHist->Draw("P");
  vtxEffVsRecoTrkHist->SetMarkerColor(seabornRed);
  vtxEffVsRecoTrkHist->SetMarkerStyle(8);
  vtxEffVsRecoTrkHist->SetMarkerSize(1.2);
  vtxEffVsRecoTrkHist->SetTitle("Vertexing Efficiency vs RC Tracks");
  vtxEffVsRecoTrkHist->GetXaxis()->SetTitle("# of RC Tracks");
  vtxEffVsRecoTrkHist->GetYaxis()->SetTitle("Vertexing Efficiency");
  vtxEffVsRecoTrkHist->GetYaxis()->SetRangeUser(0, 1.2);

  if(PRINT) c6->Print((results_path+"/vertexing/vtxEffVsRCTrks.png").c_str()); // Vertexing Efficiency versus RC Tracks
  delete c6;
	
 // Reconstructed Vertex vx versus vy
  TCanvas *c7 = new TCanvas("c7","Reconstructed Vertices",800,600);
  c7->Clear();
  c7->Divide(1,1);

  c7->cd(1);
  recoVtxYvsXHist->Draw("COLZ");
  recoVtxYvsXHist->SetTitle("Reconstructed Vertex: v_{x} versus v_{y}");
  recoVtxYvsXHist->GetXaxis()->SetTitle("x-coordinate (in mm)");
  recoVtxYvsXHist->GetYaxis()->SetTitle("y-coordinate (in mm)");
  
  if(PRINT) c7->Print((results_path+"/vertexing/recoVertexYvsX.png").c_str()); // RC Vertex vx versus vy
  delete c7;
  
  // Reconstructed Vertex vr versus vz
  TCanvas *c8 = new TCanvas("c8","Reconstructed Vertices",800,600);
  c8->Clear();
  c8->Divide(1,1);
  
  c8->cd(1);
  recoVtxRvsZHist->Draw("COLZ");
  recoVtxRvsZHist->SetTitle("Reconstrcuted Vertex: v_{r} versus v_{z}");
  recoVtxRvsZHist->GetXaxis()->SetTitle("z-coordinate (in mm)");
  recoVtxRvsZHist->GetYaxis()->SetTitle("#sqrt{x^{2} + y^{2}} (in mm)");


  if(PRINT) c8->Print((results_path+"/vertexing/recoVertexRvsZ.png").c_str()); // RC Vertex vr versus vz
  delete c8;
  
  ////////////////////////  Resolution Plots  ////////////////////////
  
  //Vertex Resolution vs MC Tracks
  TCanvas *c9 = new TCanvas("c9","VtxResX vs MCTrks",800,600);
  c9->Clear();
  c9->Divide(1,1);

  c9->cd(1);
  vtxResXvsGenTrkHist->Draw("COLZ");
  vtxResXvsGenTrkHist->SetTitle("Vertex Resolution X vs MC Tracks");
  vtxResXvsGenTrkHist->GetXaxis()->SetTitle("# of MC Tracks");
  vtxResXvsGenTrkHist->GetYaxis()->SetTitle("recVtx_x - mcVtx_x (in mm)");

  TCanvas *c10 = new TCanvas("c10","VtxResY vs MCTrks",800,600);
  c10->Clear();
  c10->Divide(1,1);

  c10->cd(1);
  vtxResYvsGenTrkHist->Draw("COLZ");
  vtxResYvsGenTrkHist->SetTitle("Vertex Resolution Y vs MC Tracks");
  vtxResYvsGenTrkHist->GetXaxis()->SetTitle("# of MC Tracks");
  vtxResYvsGenTrkHist->GetYaxis()->SetTitle("recVtx_y - mcVtx_y (in mm)");
  
  TCanvas *c11 = new TCanvas("c11","VtxResZ vs MCTrks",800,600);
  c11->Clear();
  c11->Divide(1,1);

  c11->cd(1);
  vtxResZvsGenTrkHist->Draw("COLZ");
  vtxResZvsGenTrkHist->SetTitle("Vertex Resolution Z vs MC Tracks");
  vtxResZvsGenTrkHist->GetXaxis()->SetTitle("# of MC Tracks");
  vtxResZvsGenTrkHist->GetYaxis()->SetTitle("recVtx_z - mcVtx_z (in mm)");
  
  
  if(PRINT) c9->Print((results_path+"/vertexing/vtxResXvsMCTrks.png").c_str()); // Vertexing Resolution versus MC Tracks
  delete c9;
  if(PRINT) c10->Print((results_path+"/vertexing/vtxResYvsMCTrks.png").c_str()); // Vertexing Resolution versus MC Tracks
  delete c10;
  if(PRINT) c11->Print((results_path+"/vertexing/vtxResZvsMCTrks.png").c_str()); // Vertexing Resolution versus MC Tracks
  delete c11;
  
  //Resolution versus Reconstructed Particles
  TCanvas *c12 = new TCanvas("c12","VtxResX vs recoTrks",800,600);
  c12->Clear();
  c12->Divide(1,1);

  c12->cd(1);
  vtxResXvsRecoTrkHist->Draw("COLZ");
  vtxResXvsRecoTrkHist->SetTitle("Vertex Resolution X vs Reconstructed Tracks");
  vtxResXvsRecoTrkHist->GetXaxis()->SetTitle("# of RC Tracks");
  vtxResXvsRecoTrkHist->GetYaxis()->SetTitle("recVtx_x - mcVtx_x (in mm)");

  TCanvas *c13 = new TCanvas("c13","VtxResY vs recoTrks",800,600);
  c13->Clear();
  c13->Divide(1,1);

  c13->cd(1);
  vtxResYvsRecoTrkHist->Draw("COLZ");
  vtxResYvsRecoTrkHist->SetTitle("Vertex Resolution Y vs Reconstructed Tracks");
  vtxResYvsRecoTrkHist->GetXaxis()->SetTitle("# of RC Tracks");
  vtxResYvsRecoTrkHist->GetYaxis()->SetTitle("recVtx_y - mcVtx_y (in mm)");
  
  TCanvas *c14 = new TCanvas("c14","VtxResZ vs recoTrks",800,600);
  c14->Clear();
  c14->Divide(1,1);

  c14->cd(1);
  vtxResZvsRecoTrkHist->Draw("COLZ");
  vtxResZvsRecoTrkHist->SetTitle("Vertex Resolution Z vs Reconstructed Tracks");
  vtxResZvsRecoTrkHist->GetXaxis()->SetTitle("# of RC Tracks");
  vtxResZvsRecoTrkHist->GetYaxis()->SetTitle("recVtx_z - mcVtx_z (in mm)");
  
  
  if(PRINT) c12->Print((results_path+"/vertexing/vtxResXvsRCTrks.png").c_str()); // Vertexing Resolution versus MC Tracks
  delete c12;
  if(PRINT) c13->Print((results_path+"/vertexing/vtxResYvsRCTrks.png").c_str()); // Vertexing Resolution versus MC Tracks
  delete c13;
  if(PRINT) c14->Print((results_path+"/vertexing/vtxResZvsRCTrks.png").c_str()); // Vertexing Resolution versus MC Tracks
  delete c14;
  
  // Fitted plots v/s MC tracks
  TCanvas *c15 = new TCanvas("c15","Vertex Resolution vs MC Tracks",800,600);
  c15->Clear();
  c15->Divide(1,1);

  c15->cd(1);
  
  TF1 *myfunction = new TF1("fit", StudentT, -2, 2, 4);
 
  myfunction->SetParameters(vtxResXvsGenTrkHist->GetMaximum(), 0, 0.05, 1);
  vtxResXvsGenTrkHist->FitSlicesY(myfunction, 0, -1, 10);
  
  myfunction->SetParameters(vtxResYvsGenTrkHist->GetMaximum(), 0, 0.05, 1);
  vtxResYvsGenTrkHist->FitSlicesY(myfunction, 0, -1, 10);
  
  myfunction->SetParameters(vtxResZvsGenTrkHist->GetMaximum(), 0, 0.05, 1);
  vtxResZvsGenTrkHist->FitSlicesY(myfunction, 0, -1, 10);
  
  TH1D *resXsigma = (TH1D*)gDirectory->Get("vtxResXvsGenTrk_2");
  TH1D *resYsigma = (TH1D*)gDirectory->Get("vtxResYvsGenTrk_2");
  TH1D *resZsigma = (TH1D*)gDirectory->Get("vtxResZvsGenTrk_2");
  
  resXsigma->Draw("P");
  resYsigma->Draw("PSAME");
  resZsigma->Draw("PSAME");
  
  resXsigma->SetMarkerStyle(20);
  resYsigma->SetMarkerStyle(21);
  resZsigma->SetMarkerStyle(22);
  resXsigma->SetMarkerSize(1.2);
  resYsigma->SetMarkerSize(1.2);
  resZsigma->SetMarkerSize(1.2);
  resXsigma->SetMarkerColor(seabornRed);
  resYsigma->SetMarkerColor(seabornGreen);
  resZsigma->SetMarkerColor(seabornBlue);
  
  resXsigma->SetTitle("Resolution Sigma vs MC Tracks");
  resYsigma->SetTitle("Vertex Resolution vs MC Tracks");
  resZsigma->SetTitle("Vertex Resolution vs MC Tracks");
  resZsigma->GetXaxis()->SetTitle("# of MC Tracks");
  resXsigma->GetYaxis()->SetTitle("#sigma (mm)");
  resYsigma->GetYaxis()->SetTitle("#sigma (mm)");
  resZsigma->GetYaxis()->SetTitle("#sigma (mm)");
  resXsigma->GetYaxis()->SetRangeUser(0, 1);
  resYsigma->GetYaxis()->SetRangeUser(0, 1);
  resZsigma->GetYaxis()->SetRangeUser(0, 1);

  TLegend *legend15 = new TLegend(0.7, 0.7, 0.9, 0.9); // Adjust the coordinates as needed
  legend15->AddEntry(resXsigma, "v_{x}", "lep");
  legend15->AddEntry(resYsigma, "v_{y}", "lep");
  legend15->AddEntry(resZsigma, "v_{z}", "lep");
  legend15->Draw();
  
  if(PRINT) c15->Print((results_path+"/vertexing/vtxResSigmaVsMCTrks.png").c_str()); // Vertex Resolution Sigma versus MC Tracks
  delete c15;
 
  TCanvas *c16 = new TCanvas("c16","Vertex Resolution Mean vs MC Tracks",800,600);
  c16->Clear();
  c16->Divide(1,1);

  c16->cd(1);
  
  TH1D *resXmean = (TH1D*)gDirectory->Get("vtxResXvsGenTrk_1");
  TH1D *resYmean = (TH1D*)gDirectory->Get("vtxResYvsGenTrk_1");
  TH1D *resZmean = (TH1D*)gDirectory->Get("vtxResZvsGenTrk_1");
  
  resXmean->Draw("P");
  resYmean->Draw("PSAME");
  resZmean->Draw("PSAME");
  
  resXmean->SetMarkerStyle(20);
  resYmean->SetMarkerStyle(21);
  resZmean->SetMarkerStyle(22);
  resXmean->SetMarkerSize(1.2);
  resYmean->SetMarkerSize(1.2);
  resZmean->SetMarkerSize(1.2);
  resXmean->SetMarkerColor(seabornRed);
  resYmean->SetMarkerColor(seabornGreen);
  resZmean->SetMarkerColor(seabornBlue);
  
  resXmean->SetTitle("Resolution Mean vs MC Tracks");
  resYmean->SetTitle("Vertex Resolution vs MC Tracks");
  resZmean->SetTitle("Vertex Resolution vs MC Tracks");
  resZmean->GetXaxis()->SetTitle("# of MC Tracks");
  resXmean->GetYaxis()->SetTitle("#mu (mm)");
  resYmean->GetYaxis()->SetTitle("#mu (mm)");
  resZmean->GetYaxis()->SetTitle("#mu (mm)");
  resXmean->GetYaxis()->SetRangeUser(-1, 1);
  resYmean->GetYaxis()->SetRangeUser(-1, 1);
  resZmean->GetYaxis()->SetRangeUser(-1, 1);

  TLegend *legend16 = new TLegend(0.7, 0.7, 0.9, 0.9); // Adjust the coordinates as needed
  legend16->AddEntry(resXmean, "v_{x}", "lep");
  legend16->AddEntry(resYmean, "v_{y}", "lep");
  legend16->AddEntry(resZmean, "v_{z}", "lep");
  legend16->Draw();
  
  if(PRINT) c16->Print((results_path+"/vertexing/vtxResMeanVsMCTrks.png").c_str()); // Vertex Resolution Mean versus MC Tracks
  delete c16;
  
  //Fitted Plots v/s RC tracks
  TCanvas *c17 = new TCanvas("c17","Vertex Resolution vs RC Tracks",800,600);
  c17->Clear();
  c17->Divide(1,1);

  c17->cd(1);
  
  myfunction->SetParameters(vtxResXvsRecoTrkHist->GetMaximum(), 0, 0.05, 1);
  vtxResXvsRecoTrkHist->FitSlicesY(myfunction, 0, -1, 10);
  
  myfunction->SetParameters(vtxResYvsRecoTrkHist->GetMaximum(), 0, 0.05, 1);
  vtxResYvsRecoTrkHist->FitSlicesY(myfunction, 0, -1, 10);
  
  myfunction->SetParameters(vtxResZvsRecoTrkHist->GetMaximum(), 0, 0.05, 1);
  vtxResZvsRecoTrkHist->FitSlicesY(myfunction, 0, -1, 10);
  
  TH1D *resXsigma2 = (TH1D*)gDirectory->Get("vtxResXvsRecoTrk_2");
  TH1D *resYsigma2 = (TH1D*)gDirectory->Get("vtxResYvsRecoTrk_2");
  TH1D *resZsigma2 = (TH1D*)gDirectory->Get("vtxResZvsRecoTrk_2");
  
  resXsigma2->Draw("P");
  resYsigma2->Draw("PSAME");
  resZsigma2->Draw("PSAME");
  
  resXsigma2->SetMarkerStyle(20);
  resYsigma2->SetMarkerStyle(21);
  resZsigma2->SetMarkerStyle(22);
  resXsigma2->SetMarkerSize(1.2);
  resYsigma2->SetMarkerSize(1.2);
  resZsigma2->SetMarkerSize(1.2);
  resXsigma2->SetMarkerColor(seabornRed);
  resYsigma2->SetMarkerColor(seabornGreen);
  resZsigma2->SetMarkerColor(seabornBlue);
  
  resXsigma2->SetTitle("Resolution Sigma vs RC Tracks");
  resYsigma2->SetTitle("Vertex Resolution vs RC Tracks");
  resZsigma2->SetTitle("Vertex Resolution vs RC Tracks");
  resZsigma2->GetXaxis()->SetTitle("# of RC Tracks");
  resXsigma2->GetYaxis()->SetTitle("#sigma (mm)");
  resYsigma2->GetYaxis()->SetTitle("#sigma (mm)");
  resZsigma2->GetYaxis()->SetTitle("#sigma (mm)");
  resXsigma2->GetYaxis()->SetRangeUser(0, 1);
  resYsigma2->GetYaxis()->SetRangeUser(0, 1);
  resZsigma2->GetYaxis()->SetRangeUser(0, 1);

  TLegend *legend17 = new TLegend(0.7, 0.7, 0.9, 0.9); // Adjust the coordinates as needed
  legend17->AddEntry(resXsigma2, "v_{x}", "lep");
  legend17->AddEntry(resYsigma2, "v_{y}", "lep");
  legend17->AddEntry(resZsigma2, "v_{z}", "lep");
  legend17->Draw();
  
  if(PRINT) c17->Print((results_path+"/vertexing/vtxResSigmaVsRCTrks.png").c_str()); // Vertex Resolution Sigma versus RC Tracks
  delete c17;
 
  TCanvas *c18 = new TCanvas("c18","Vertex Resolution Mean vs RC Tracks",800,600);
  c18->Clear();
  c18->Divide(1,1);

  c18->cd(1);
  
  TH1D *resXmean2 = (TH1D*)gDirectory->Get("vtxResXvsRecoTrk_1");
  TH1D *resYmean2 = (TH1D*)gDirectory->Get("vtxResYvsRecoTrk_1");
  TH1D *resZmean2 = (TH1D*)gDirectory->Get("vtxResZvsRecoTrk_1");
  
  resXmean2->Draw("P");
  resYmean2->Draw("PSAME");
  resZmean2->Draw("PSAME");
  
  resXmean2->SetMarkerStyle(20);
  resYmean2->SetMarkerStyle(21);
  resZmean2->SetMarkerStyle(22);
  resXmean2->SetMarkerSize(1.2);
  resYmean2->SetMarkerSize(1.2);
  resZmean2->SetMarkerSize(1.2);
  resXmean2->SetMarkerColor(seabornRed);
  resYmean2->SetMarkerColor(seabornGreen);
  resZmean2->SetMarkerColor(seabornBlue);
 
  resXmean2->SetTitle("Resolution Mean vs RC Tracks");
  resYmean2->SetTitle("Vertex Resolution vs RC Tracks");
  resZmean2->SetTitle("Vertex Resolution vs RC Tracks");
  resZmean2->GetXaxis()->SetTitle("# of MC Tracks");
  resXmean2->GetYaxis()->SetTitle("#mu (mm)");
  resYmean2->GetYaxis()->SetTitle("#mu (mm)");
  resZmean2->GetYaxis()->SetTitle("#mu (mm)");
  resXmean2->GetYaxis()->SetRangeUser(-0.4, 0.4);
  resYmean2->GetYaxis()->SetRangeUser(-0.4, 0.4);
  resZmean2->GetYaxis()->SetRangeUser(-0.4, 0.4);

  TLegend *legend18 = new TLegend(0.7, 0.7, 0.9, 0.9); // Adjust the coordinates as needed
  legend18->AddEntry(resXmean, "v_{x}", "lep");
  legend18->AddEntry(resYmean, "v_{y}", "lep");
  legend18->AddEntry(resZmean, "v_{z}", "lep");
  legend18->Draw();
  
  if(PRINT) c18->Print((results_path+"/vertexing/vtxResMeanVsRCTrks.png").c_str()); // Vertex Resolution Mean versus RC Tracks
  delete c18;
  delete tree;
  return 0;
}
  
double StudentT(double *x, double *par){
	double norm = par[0];
  	double mean = par[1];
  	double sigma = par[2];
  	double nu = par[3];

	double pi = 3.14;
  	double st = norm * (TMath::Gamma((nu+1.0)/2.0)/(TMath::Gamma(nu/2.0)*TMath::Sqrt(pi*nu)*sigma)) * TMath::Power( (1.0+TMath::Power((x[0]-mean)/sigma,2.0)/nu), (-(nu+1.0)/2.0) );
  return st;
}
