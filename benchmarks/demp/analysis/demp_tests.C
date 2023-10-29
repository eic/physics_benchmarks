//------------------
void DEMP_reco(const char* fname = "rec_demp.root"){

    //Define Style
    gStyle->SetOptStat(0);
    gStyle->SetPadBorderMode(0);
    gStyle->SetFrameBorderMode(0);
    gStyle->SetFrameLineWidth(2);
    gStyle->SetLabelSize(0.035,"X");
    gStyle->SetLabelSize(0.035,"Y");
    //gStyle->SetLabelOffset(0.01,"X");
    //gStyle->SetLabelOffset(0.01,"Y");
    gStyle->SetTitleXSize(0.04);
    gStyle->SetTitleXOffset(0.9);
    gStyle->SetTitleYSize(0.04);
    gStyle->SetTitleYOffset(0.9);

    //Define histograms
    TH2 *h1_elec = new TH2D("h1_elec","Scattered electron true momentum vs. polar angle",100,110,170,100,4.6,6.8);
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

    //Other variables
    TLorentzVector gen_vec;
    TVector3 gen_vertex;
    TLorentzVector rec_vec;
    int counter(0);

    //Loop over events
    while (tr.Next()) {
	
	if(counter%100==0) cout<<"Analyzing event "<<counter<<endl;
	counter++;

        //Loop over generated particles, select primary electron, pi-plus, neutron
        for(int igen=0;igen<gen_status.GetSize();igen++){
        	if(gen_status[igen]==1){
                	gen_vec.SetXYZM(gen_px[igen],gen_py[igen],gen_pz[igen],gen_mass[igen]);
                	gen_vertex.SetXYZ(gen_vx[igen],gen_vy[igen],gen_vz[igen]);
			
			if(gen_pid[igen]==11)  h1_elec->Fill(gen_vec.Theta()*TMath::RadToDeg(),gen_vec.P());
			if(gen_pid[igen]==211) h1_pion->Fill(gen_vec.Theta()*TMath::RadToDeg(),gen_vec.P());
			if(gen_pid[igen]==2112) h1_neut->Fill(gen_vec.Theta()*TMath::RadToDeg(),gen_vec.P());
                }
        } //End loop over generated particles
     
	//Loop over reconstructed tracks for electron and pion
	for(int irec=0;irec<rec_pid.GetSize();irec++){
		rec_vec.SetXYZM(rec_px[irec],rec_py[irec],rec_pz[irec],rec_mass[irec]);

		if(rec_pid[irec]==11)  h2_elec->Fill(rec_vec.Theta()*TMath::RadToDeg(),rec_vec.P());
                if(rec_pid[irec]==211) h2_pion->Fill(rec_vec.Theta()*TMath::RadToDeg(),rec_vec.P());

	}//End loop over reconstructed particles

     } //End loop over events

    //Make plots
    TCanvas *c1 = new TCanvas("c1");
    h1_elec->Draw("colz");

    TCanvas *c2 = new TCanvas("c2");
    h1_pion->Draw("colz");

    TCanvas *c3 = new TCanvas("c3");
    h1_neut->Draw("colz");

    TCanvas *c4 = new TCanvas("c4");
    h2_elec->Draw("colz");

    TCanvas *c5 = new TCanvas("c5");
    h2_pion->Draw("colz");

    //Print plots to file
    c1->Print("results/demp/DEMP_reco.pdf[");
    c1->Print("results/demp/DEMP_reco.pdf");
    c2->Print("results/demp/DEMP_reco.pdf");
    c3->Print("results/demp/DEMP_reco.pdf");
    c4->Print("results/demp/DEMP_reco.pdf");
    c5->Print("results/demp/DEMP_reco.pdf");  
    c5->Print("results/demp/DEMP_reco.pdf]");

}
