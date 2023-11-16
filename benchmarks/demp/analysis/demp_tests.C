#include "ZDC_neutron_recon.h"

void draw_title(const char* txt){
  cout << txt <<endl;
  TLatex *title = new TLatex(0.5, 0.96,txt);
  title->SetTextColor(kBlack);
  title->SetNDC();
  title->SetTextAlign(22);
  //title->SetTextSize(0.05); // Adjust the size as needed
  title->Draw();
  
  gPad->Update();
}
//------------------
void demp_tests(const char* fname = "rec_demp.root"){

    //Define Style
    gStyle->SetOptStat(0);
    //titles in ROOT don't do what you want them to do.
    //use TLatex instead
    gStyle->SetOptTitle(0);
    gStyle->SetPadTopMargin(0.1);
    gStyle->SetOptFit(1);
    //gStyle->SetTitleSize(1,"T");
    //gStyle->SetTitleOffset(1, "T");
    //gStyle->SetTitleSize(14, "T");
    /*gStyle->SetPadBorderMode(0);
    gStyle->SetFrameBorderMode(0);
    gStyle->SetFrameLineWidth(2);
    gStyle->SetLabelSize(0.035,"X");
    gStyle->SetLabelSize(0.035,"Y");
    //gStyle->SetLabelOffset(0.01,"X");
    //gStyle->SetLabelOffset(0.01,"Y");
    
    gStyle->SetTitleXSize(0.04);
    gStyle->SetTitleXOffset(0.9);
    gStyle->SetTitleYSize(0.04);
    gStyle->SetTitleYOffset(0.9);//*/

    //Define histograms
    TH2 *h1_elec = new TH2D("h1_elec","Scattered electron true momentum vs polar angle",100,110,170,100,4.6,6.8);
    h1_elec->GetXaxis()->SetTitle("#theta [deg]"); h1_elec->GetXaxis()->CenterTitle();
    h1_elec->GetYaxis()->SetTitle("p [GeV/c]"); h1_elec->GetYaxis()->CenterTitle();
   
    TH2 *h1_pion = new TH2D("h1_pion","#pi^{+} true momentum vs. polar angle",100,0,50,100,0,50);
    h1_pion->GetXaxis()->SetTitle("#theta [deg]"); h1_pion->GetXaxis()->CenterTitle();
    h1_pion->GetYaxis()->SetTitle("p [GeV/c]"); h1_pion->GetYaxis()->CenterTitle();
 
    TH2 *h1_neut = new TH2D("h1_neut","Neutron true momentum vs. polar angle",100,0.4,2.4,100,40,100);
    h1_neut->GetXaxis()->SetTitle("#theta [deg]"); h1_neut->GetXaxis()->CenterTitle();
    h1_neut->GetYaxis()->SetTitle("p [GeV/c]"); h1_neut->GetYaxis()->CenterTitle();

    TH2 *h2_elec = new TH2D("h2_elec","Scattered electron reconstructed momentum vs. polar angle",100,110,170,100,4.6,6.8);
    h2_elec->GetXaxis()->SetTitle("#theta [deg]"); h2_elec->GetXaxis()->CenterTitle();
    h2_elec->GetYaxis()->SetTitle("p [GeV/c]"); h2_elec->GetYaxis()->CenterTitle();

    TH2 *h2_pion = new TH2D("h2_pion","#pi^{+} reconstructed momentum vs. polar angle",100,0,50,100,0,50);
    h2_pion->GetXaxis()->SetTitle("#theta [deg]"); h2_pion->GetXaxis()->CenterTitle();
    h2_pion->GetYaxis()->SetTitle("p [GeV/c]"); h2_pion->GetYaxis()->CenterTitle();

    TH2 *h2_neut = new TH2D("h2_neut","Neutron reconstructed momentum vs. polar angle",100,0.4,2.4,100,40,100);
    h2_neut->GetXaxis()->SetTitle("#theta [deg]"); h2_neut->GetXaxis()->CenterTitle();
    h2_neut->GetYaxis()->SetTitle("p [GeV/c]"); h2_neut->GetYaxis()->CenterTitle();

    TH1 *h3_neut = new TH1D("h3_neut","Neutron theta resolution",100,-1,1);
    h3_neut->GetXaxis()->SetTitle("#Delta#theta [mrad]"); h3_neut->GetXaxis()->CenterTitle();
    h3_neut->GetYaxis()->SetTitle("events"); h3_neut->GetYaxis()->CenterTitle();

    TH1 *h4_neut = new TH1D("h4_neut","Neutron energy resolution",100,-30,10);
    h4_neut->GetXaxis()->SetTitle("#Delta E/E [%]"); h4_neut->GetXaxis()->CenterTitle();
    h4_neut->GetYaxis()->SetTitle("events"); h4_neut->GetYaxis()->CenterTitle();

    TH2 *h5_neut = new TH2D("h5_neut","Neutron truth azimuth vs. polar angle",100,0.4,2.4,100,160,200);
    h5_neut->GetXaxis()->SetTitle("#theta [deg]"); h5_neut->GetXaxis()->CenterTitle();
    h5_neut->GetYaxis()->SetTitle("#phi [deg]"); h5_neut->GetYaxis()->CenterTitle();
    
    TH2 *h6_neut = new TH2D("h6_neut","Neutron reconstructed azimuth vs. polar angle",100,0.4,2.4,100,160,200);
    h6_neut->GetXaxis()->SetTitle("#theta [deg]"); h6_neut->GetXaxis()->CenterTitle();
    h6_neut->GetYaxis()->SetTitle("#phi [deg]"); h6_neut->GetYaxis()->CenterTitle();

    TH2 *h7_neut = new TH2D("h7_neut","Neutron truth azimuth vs. polar angle around p axis",100,0,0.5,100,-180,180);
    h7_neut->GetXaxis()->SetTitle("#theta^{*} [deg]"); h7_neut->GetXaxis()->CenterTitle();
    h7_neut->GetYaxis()->SetTitle("#phi^{*} [deg]"); h7_neut->GetYaxis()->CenterTitle();

    TH2 *h8_neut = new TH2D("h8_neut","Neutron reconstructed azimuth vs. polar angle around p axis",100,0,0.5,100,-180,180);
    h8_neut->GetXaxis()->SetTitle("#theta^{*} [deg]"); h8_neut->GetXaxis()->CenterTitle();
    h8_neut->GetYaxis()->SetTitle("#phi^{*} [deg]"); h8_neut->GetYaxis()->CenterTitle();

    
    TFile *f = new TFile(fname);
    TTree *tree = (TTree*) f->Get("events");

    //Create Array Reader
    TTreeReader tr(tree);

    TTreeReaderArray<int>   gen_status(tr, "MCParticles.generatorStatus");
    TTreeReaderArray<int>   gen_pid(tr, "MCParticles.PDG");
    TTreeReaderArray<float> gen_px(tr, "MCParticles.momentum.x");
    TTreeReaderArray<float> gen_py(tr, "MCParticles.momentum.y");
    TTreeReaderArray<float> gen_pz(tr, "MCParticles.momentum.z");
    TTreeReaderArray<double> gen_mass(tr, "MCParticles.mass");
    TTreeReaderArray<float> gen_charge(tr, "MCParticles.charge");
    TTreeReaderArray<double> gen_vx(tr, "MCParticles.vertex.x");
    TTreeReaderArray<double> gen_vy(tr, "MCParticles.vertex.y");
    TTreeReaderArray<double> gen_vz(tr, "MCParticles.vertex.z");
  
    TTreeReaderArray<int> rec_pid(tr, "ReconstructedChargedParticles.PDG");
    TTreeReaderArray<float> rec_px(tr, "ReconstructedChargedParticles.momentum.x");
    TTreeReaderArray<float> rec_py(tr, "ReconstructedChargedParticles.momentum.y");
    TTreeReaderArray<float> rec_pz(tr, "ReconstructedChargedParticles.momentum.z");
    TTreeReaderArray<float> rec_mass(tr, "ReconstructedChargedParticles.mass");

    //local positions and energies in the ZDC
    TTreeReaderArray<float> zdc_hit_energies(tr, "ZDCRecHits.energy");
    TTreeReaderArray<float> zdc_hit_times(tr, "ZDCRecHits.time");
    TTreeReaderArray<float> zdc_hit_local_x(tr, "ZDCRecHits.local.x");
    TTreeReaderArray<float> zdc_hit_local_y(tr, "ZDCRecHits.local.y");
    TTreeReaderArray<float> zdc_hit_local_z(tr, "ZDCRecHits.local.z");
    
    //Other variables
    TLorentzVector gen_vec;
    TVector3 gen_vertex;
    TLorentzVector rec_vec;
    int counter(0);

    //Loop over events
    while (tr.Next()) {
	
	if(counter%100==0) cout<<"Analyzing event "<<counter<<endl;
	counter++;
	double theta_neut_truth;
	double E_neut_truth;
        //Loop over generated particles, select primary electron, pi-plus, neutron
        for(int igen=0;igen<gen_status.GetSize();igen++){
        	if(gen_status[igen]==1){
                	gen_vec.SetXYZM(gen_px[igen],gen_py[igen],gen_pz[igen],gen_mass[igen]);
                	gen_vertex.SetXYZ(gen_vx[igen],gen_vy[igen],gen_vz[igen]);
			
			if(gen_pid[igen]==11)  h1_elec->Fill(gen_vec.Theta()*TMath::RadToDeg(),gen_vec.P());
			if(gen_pid[igen]==211) h1_pion->Fill(gen_vec.Theta()*TMath::RadToDeg(),gen_vec.P());
			if(gen_pid[igen]==2112) {
			  h1_neut->Fill(gen_vec.Theta()*TMath::RadToDeg(),gen_vec.P());
			  h5_neut->Fill(gen_vec.Theta()*TMath::RadToDeg(),(gen_vec.Phi()+(gen_vec.Phi()<0)*2*M_PI)*TMath::RadToDeg());
			  theta_neut_truth=gen_vec.Theta();
			  E_neut_truth=gen_vec.E();

			  TLorentzVector rotated = gen_vec;
			  rotated.RotateY(0.025);
			  h7_neut->Fill(rotated.Theta()*TMath::RadToDeg(),rotated.Phi()*TMath::RadToDeg());
			}
			  
                }
        } //End loop over generated particles
     
	//Loop over reconstructed tracks for electron and pion
	for(int irec=0;irec<rec_pid.GetSize();irec++){
		rec_vec.SetXYZM(rec_px[irec],rec_py[irec],rec_pz[irec],rec_mass[irec]);

		if(rec_pid[irec]==11)  h2_elec->Fill(rec_vec.Theta()*TMath::RadToDeg(),rec_vec.P());
                if(rec_pid[irec]==211) h2_pion->Fill(rec_vec.Theta()*TMath::RadToDeg(),rec_vec.P());

	}//End loop over reconstructed particles

	TLorentzVector neutron = ZDC_neutron_recon(zdc_hit_energies, zdc_hit_times, zdc_hit_local_x, zdc_hit_local_y,
						   zdc_hit_local_z);
	//cout << "neutron theta " << neutron.Theta() << endl;
        h2_neut->Fill(neutron.Theta()*TMath::RadToDeg(),neutron.P());
	h3_neut->Fill((neutron.Theta()-theta_neut_truth)*1000);
        h4_neut->Fill((neutron.E()/E_neut_truth-1)*100);
	h6_neut->Fill(neutron.Theta()*TMath::RadToDeg(),(neutron.Phi()+(neutron.Phi()<0)*2*M_PI)*TMath::RadToDeg());

	TLorentzVector rotated = neutron;
	rotated.RotateY(0.025);
	h8_neut->Fill(rotated.Theta()*TMath::RadToDeg(),rotated.Phi()*TMath::RadToDeg());
	
     } //End loop over events

    //Make plots
    int w=800, h=600;
    TCanvas *c1 = new TCanvas("c1", "", w,h);
    h1_elec->Draw("colz");
    draw_title(h1_elec->GetTitle());
    
    TCanvas *c2 = new TCanvas("c2", "", w,h);
    h1_pion->Draw("colz");
    draw_title(h1_pion->GetTitle());
    
    TCanvas *c3 = new TCanvas("c3", "", w,h);
    h1_neut->Draw("colz");
    draw_title(h1_neut->GetTitle());
    cout << "generated neutron count " << h1_neut->Integral();
    
    TCanvas *c4 = new TCanvas("c4", "", w,h);
    h2_elec->Draw("colz");
    draw_title(h2_elec->GetTitle());

    TCanvas *c5 = new TCanvas("c5", "", w,h);
    h2_pion->Draw("colz");
    draw_title(h2_pion->GetTitle());

    TCanvas *c6 = new TCanvas("c6", "", w,h);
    h2_neut->Draw("colz");
    draw_title(h2_neut->GetTitle());
    cout << "reconstructed neutron count " << h2_neut->Integral();


    TCanvas *c7 = new TCanvas("c7", "", w,h);
    TF1* fnc = new TF1("f1", "gaus(0)", -1,1);
    h3_neut->Draw("colz");
    h3_neut->Fit(fnc);
    draw_title(h3_neut->GetTitle());
    
    TCanvas *c8 = new TCanvas("c8", "", w,h);
    h4_neut->Draw("colz");
    h4_neut->Fit(fnc);    
    draw_title(h4_neut->GetTitle());

    TCanvas *c9 = new TCanvas("c9", "", w,h);
    h5_neut->Draw("colz");
    draw_title(h5_neut->GetTitle());

    TCanvas *c10 = new TCanvas("c10", "", w,h);
    h6_neut->Draw("colz");
    draw_title(h6_neut->GetTitle());

    TCanvas *c11 = new TCanvas("c11", "", w,h);
    h7_neut->Draw("colz");
    draw_title(h7_neut->GetTitle());

    TCanvas *c12 = new TCanvas("c12", "", w,h);
    h8_neut->Draw("colz");
    draw_title(h8_neut->GetTitle());
    
    //Print plots to file
    c1->Print("results/demp/DEMP_reco.pdf[");
    c1->Print("results/demp/DEMP_reco.pdf");
    c2->Print("results/demp/DEMP_reco.pdf");
    c3->Print("results/demp/DEMP_reco.pdf");
    c4->Print("results/demp/DEMP_reco.pdf");
    c5->Print("results/demp/DEMP_reco.pdf");  
    c6->Print("results/demp/DEMP_reco.pdf");
    c7->Print("results/demp/DEMP_reco.pdf");
    c8->Print("results/demp/DEMP_reco.pdf");
    c9->Print("results/demp/DEMP_reco.pdf");
    c10->Print("results/demp/DEMP_reco.pdf");
    c11->Print("results/demp/DEMP_reco.pdf");
    c12->Print("results/demp/DEMP_reco.pdf");
    c12->Print("results/demp/DEMP_reco.pdf]");
}
