#include "ZDC_neutron_recon.h" 

//------------------
void DEMP_reco(){

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
 
    TH2 *h1_neut = new TH2D("h1_neut","Neutron true energy vs. polar angle",100,0.4,2.4,100,40,100);
    h1_neut->GetXaxis()->SetTitle("#theta [deg]"); h1_neut->GetXaxis()->CenterTitle();
    h1_neut->GetYaxis()->SetTitle("E [GeV]"); h1_neut->GetYaxis()->CenterTitle();

    TH2 *h2_elec = new TH2D("h2_elec","Scattered electron reconstructed momentum vs. polar angle",100,110,170,100,4.6,6.8);
    h2_elec->GetXaxis()->SetTitle("#theta [deg]"); h2_elec->GetXaxis()->CenterTitle();
    h2_elec->GetYaxis()->SetTitle("p [GeV/c]"); h2_elec->GetYaxis()->CenterTitle();

    TH2 *h2_pion = new TH2D("h2_pion","#pi^{+} reconstructed momentum vs. polar angle",100,0,50,100,0,50);
    h2_pion->GetXaxis()->SetTitle("#theta [deg]"); h2_pion->GetXaxis()->CenterTitle();
    h2_pion->GetYaxis()->SetTitle("p [GeV/c]"); h2_pion->GetYaxis()->CenterTitle();

    TH2 *h2_neut = new TH2D("h2_neut","Neutron reconstructed energy vs. polar angle",100,0.4,2.4,100,40,100);
    h2_neut->GetXaxis()->SetTitle("#theta [deg]"); h2_neut->GetXaxis()->CenterTitle();
    h2_neut->GetYaxis()->SetTitle("E [GeV]"); h2_neut->GetYaxis()->CenterTitle();

    TH2 *h3_neut = new TH2D("h3_neut","Neutron truth azimuth vs. polar angle around p axis",100,0,0.6,100,-200,200);
    h3_neut->GetXaxis()->SetTitle("#theta^{*} [deg]"); h3_neut->GetXaxis()->CenterTitle();
    h3_neut->GetYaxis()->SetTitle("#phi^{*} [deg]"); h3_neut->GetYaxis()->CenterTitle();

    TH2 *h4_neut = new TH2D("h4_neut","Neutron reconstructed azimuth vs. polar angle around p axis",100,0,0.6,100,-200,200);
    h4_neut->GetXaxis()->SetTitle("#theta^{*} [deg]"); h4_neut->GetXaxis()->CenterTitle();
    h4_neut->GetYaxis()->SetTitle("#phi^{*} [deg]"); h4_neut->GetYaxis()->CenterTitle();

    TH1 *ht_true = new TH1D("ht_true","",100,-0.1,2); //True t
    ht_true->GetXaxis()->SetTitle("-t [GeV^{2}]"); ht_true->GetXaxis()->CenterTitle();
    ht_true->GetYaxis()->SetTitle("Yield"); ht_true->GetYaxis()->CenterTitle();
    ht_true->SetLineWidth(2);ht_true->SetLineColor(kBlue);

    TH1 *ht_rec1 = new TH1D("ht_rec1","",100,-0.1,2); //Reconstructed t -- 4 vector reconstruction using only electron and pion
    ht_rec1->GetXaxis()->SetTitle("-t [GeV^{2}]"); ht_rec1->GetXaxis()->CenterTitle();
    ht_rec1->GetYaxis()->SetTitle("Yield"); ht_rec1->GetYaxis()->CenterTitle();
    ht_rec1->SetLineWidth(2);ht_rec1->SetLineColor(kRed);

    TH1 *ht_rec2 = new TH1D("ht_rec2","",100,-0.1,2); //Reconstructed t -- ECCE method using electron and pion + neutron angle
    ht_rec2->GetXaxis()->SetTitle("-t [GeV^{2}]"); ht_rec2->GetXaxis()->CenterTitle();
    ht_rec2->GetYaxis()->SetTitle("Yield"); ht_rec2->GetYaxis()->CenterTitle();
    ht_rec2->SetLineWidth(2);ht_rec2->SetLineColor(kGreen);

    TH1 *ht_rec3 = new TH1D("ht_rec3","",100,-0.1,2); //Reconstructed t -- Pt-based reconstruction using only electron and pion
    ht_rec3->GetXaxis()->SetTitle("-t [GeV^{2}]"); ht_rec3->GetXaxis()->CenterTitle();
    ht_rec3->GetYaxis()->SetTitle("Yield"); ht_rec3->GetYaxis()->CenterTitle();
    ht_rec3->SetLineWidth(2);ht_rec3->SetLineColor(kOrange);

    TH1 *ht_rec4 = new TH1D("ht_rec4","",100,-0.1,2); //Reconstructed t -- 4 vector reconstruction using only neutron
    ht_rec4->GetXaxis()->SetTitle("-t [GeV^{2}]"); ht_rec4->GetXaxis()->CenterTitle();
    ht_rec4->GetYaxis()->SetTitle("Yield"); ht_rec4->GetYaxis()->CenterTitle();
    ht_rec4->SetLineWidth(2);ht_rec4->SetLineColor(kMagenta);

    TH1 *htw_true = new TH1D("htw_true","",100,-0.1,2); //True t weighted
    htw_true->GetXaxis()->SetTitle("-t [GeV^{2}]"); htw_true->GetXaxis()->CenterTitle();
    htw_true->GetYaxis()->SetTitle("Yield [arb.]"); htw_true->GetYaxis()->CenterTitle();
    htw_true->SetLineWidth(2);htw_true->SetLineColor(kBlue);

    TH1 *htw_rec1 = new TH1D("htw_rec1","",100,-0.1,2); //Reconstructed t weighted
    htw_rec1->GetXaxis()->SetTitle("-t [GeV^{2}]"); htw_rec1->GetXaxis()->CenterTitle();
    htw_rec1->GetYaxis()->SetTitle("Yield [arb.]"); htw_rec1->GetYaxis()->CenterTitle();
    htw_rec1->SetLineWidth(2);htw_rec1->SetLineColor(kRed);

    TH1 *htw_rec2 = new TH1D("htw_rec2","",100,-0.1,2); //Reconstructed t weighted
    htw_rec2->GetXaxis()->SetTitle("-t [GeV^{2}]"); htw_rec2->GetXaxis()->CenterTitle();
    htw_rec2->GetYaxis()->SetTitle("Yield [arb.]"); htw_rec2->GetYaxis()->CenterTitle();
    htw_rec2->SetLineWidth(2);htw_rec2->SetLineColor(kGreen);

    TH1 *htw_rec3 = new TH1D("htw_rec3","",100,-0.1,2); //Reconstructed t
    htw_rec3->GetXaxis()->SetTitle("-t [GeV^{2}]"); htw_rec3->GetXaxis()->CenterTitle();
    htw_rec3->GetYaxis()->SetTitle("Yield"); htw_rec3->GetYaxis()->CenterTitle();
    htw_rec3->SetLineWidth(2);htw_rec3->SetLineColor(kOrange);

    TH1 *htw_rec4 = new TH1D("htw_rec4","",100,-0.1,2); //Reconstructed t
    htw_rec4->GetXaxis()->SetTitle("-t [GeV^{2}]"); htw_rec4->GetXaxis()->CenterTitle();
    htw_rec4->GetYaxis()->SetTitle("Yield"); htw_rec4->GetYaxis()->CenterTitle();
    htw_rec4->SetLineWidth(2);htw_rec4->SetLineColor(kMagenta);

    TFile *f = new TFile("eicrecon_out.root");
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

    TTreeReaderArray<pair<string,vector<string>>> weight_map(tr,"_stringMap");

    //Other variables
    TLorentzVector gen_vec;
    TVector3 gen_vertex;
    TLorentzVector rec_vec;
    int counter(0);

    TLorentzVector p_beam_real; //Proton beam (incorporating event-by-event fluctuations) for true t calculation
    TLorentzVector neut_true; //True neutron for true t calculation
    TLorentzVector e_beam(0.,0.,-5.,5.); //Electron beam for t reconstruction (ignoring event-by-event fluctuations)
    TLorentzVector p_beam(100.*sin(-25./1000),0,100*cos(-25./1000),100); //Proton beam for t reconstruction (ignoring event-by-event fluctuations)
    TLorentzVector e_s_rec; //Reconstructed scattered electron from tracking detector
    TLorentzVector pi_rec; //Reconstructed pi+ from tracking detector

    //Loop over events
    while (tr.Next()) {
	
	if(counter%100==0) cout<<"Analyzing event "<<counter<<endl;
	counter++;

        //Get event weight
        std::string weight_name = weight_map[0].first;
        vector<string> weight_value = weight_map[0].second;

	double weight = std::stod( *(weight_value.begin()) );
        //cout<<"Weight = "<<weight<<endl; //Check that it is working

        //Loop over generated particles, select primary electron, pi-plus, neutron
        for(int igen=0;igen<gen_status.GetSize();igen++){

                //Get Proton beam
                if(gen_status[igen]==4 && gen_pid[igen]==2212){
			p_beam_real.SetXYZM(gen_px[igen],gen_py[igen],gen_pz[igen],gen_mass[igen]);
                } 

        	if(gen_status[igen]==1){
                	gen_vec.SetXYZM(gen_px[igen],gen_py[igen],gen_pz[igen],gen_mass[igen]);
                	gen_vertex.SetXYZ(gen_vx[igen],gen_vy[igen],gen_vz[igen]);
			
			if(gen_pid[igen]==11)  h1_elec->Fill(gen_vec.Theta()*TMath::RadToDeg(),gen_vec.P());
			if(gen_pid[igen]==211) h1_pion->Fill(gen_vec.Theta()*TMath::RadToDeg(),gen_vec.P());

			if(gen_pid[igen]==2112){ 
                                h1_neut->Fill(gen_vec.Theta()*TMath::RadToDeg(),gen_vec.E());
                                neut_true = gen_vec;

				//W.r.t proton beam direction
				TLorentzVector rotated = gen_vec;
			  	rotated.RotateY(0.025);
			  	h3_neut->Fill(rotated.Theta()*TMath::RadToDeg(),rotated.Phi()*TMath::RadToDeg());

                        }

                }
        } //End loop over generated particles

        //Calculate true t
        auto t_true = (p_beam_real - neut_true).Mag2();
     
	//Loop over reconstructed tracks for electron and pion
	for(int irec=0;irec<rec_pid.GetSize();irec++){
		rec_vec.SetXYZM(rec_px[irec],rec_py[irec],rec_pz[irec],rec_mass[irec]);

		if(rec_pid[irec]==11){
			h2_elec->Fill(rec_vec.Theta()*TMath::RadToDeg(),rec_vec.P());
			e_s_rec = rec_vec;
		}

                if(rec_pid[irec]==211){
			h2_pion->Fill(rec_vec.Theta()*TMath::RadToDeg(),rec_vec.P());
			pi_rec = rec_vec;
		}

	}//End loop over reconstructed particles

	//Neutron Reconstruction using HEXSPLIT algorithm
	TLorentzVector neut_rec = ZDC_neutron_recon(zdc_hit_energies, zdc_hit_times, 
						   zdc_hit_local_x, zdc_hit_local_y,zdc_hit_local_z);

	h2_neut->Fill(neut_rec.Theta()*TMath::RadToDeg(),neut_rec.E());

        //W.r.t proton beam direction
        TLorentzVector rotated = neut_rec;
        rotated.RotateY(0.025);
        h4_neut->Fill(rotated.Theta()*TMath::RadToDeg(),rotated.Phi()*TMath::RadToDeg());

        //Calculate reconstructed t

        //Method 1
	auto t_rec1 = (e_beam - e_s_rec - pi_rec).Mag2();
	
	//Method 2
	auto p_miss = (e_beam + p_beam - e_s_rec - pi_rec);
	
	TLorentzVector neut_opt;
	auto neut_opt_px = p_miss.P()*sin(neut_rec.Theta())*cos(neut_rec.Phi());
	auto neut_opt_py = p_miss.P()*sin(neut_rec.Theta())*sin(neut_rec.Phi());
	auto neut_opt_pz = p_miss.P()*cos(neut_rec.Theta());
	double neut_mass = 0.9396;

	neut_opt.SetXYZM(neut_opt_px,neut_opt_py,neut_opt_pz,neut_mass);

        auto t_rec2 = (p_beam - neut_opt).Mag2();

        //Method 3
        auto sum_epi = e_s_rec + pi_rec;
	sum_epi.RotateY(0.025);
	auto t_rec3 = -1.*(sum_epi).Perp2(); //Make sure t is negative

        //Method 4
	auto t_rec4 = (p_beam - neut_rec).Mag2();

        //Fill additional histograms
        ht_true->Fill(-1.*t_true);
	ht_rec1->Fill(-1.*t_rec1);
	ht_rec2->Fill(-1.*t_rec2);
	ht_rec3->Fill(-1.*t_rec3);
        ht_rec4->Fill(-1.*t_rec4);

	htw_true->Fill(-1.*t_true,weight);
        htw_rec1->Fill(-1.*t_rec1,weight);
        htw_rec2->Fill(-1.*t_rec2,weight);
	htw_rec3->Fill(-1.*t_rec3,weight);
        htw_rec4->Fill(-1.*t_rec4,weight);

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

    TCanvas *c6 = new TCanvas("c6");
    h2_neut->Draw("colz");

    TCanvas *c7 = new TCanvas("c7");
    h3_neut->Draw("colz");

    TCanvas *c8 = new TCanvas("c8");
    h4_neut->Draw("colz");

    TCanvas *c9 = new TCanvas("c9");
    ht_true->Draw("hist");
    ht_rec1->Draw("hist same");
    ht_rec2->Draw("hist same");
    ht_rec3->Draw("hist same");
    ht_rec4->Draw("hist same");

    TLegend *leg1 = new TLegend(0.5,0.6,0.85,0.85);
    leg1->SetBorderSize(0);leg1->SetTextSize(0.03);
    leg1->SetFillStyle(0);
    leg1->AddEntry(ht_true,"Truth (after beam-effects)","l");
    leg1->AddEntry(ht_rec1,"Electron+Pion 4-vector","l");
    leg1->AddEntry(ht_rec4,"Neutron only","l");
    leg1->AddEntry(ht_rec3,"Electron+Pion P_{T}-based","l");
    leg1->AddEntry(ht_rec2,"ECCE paper method","l");
    leg1->Draw();

    TCanvas *c10 = new TCanvas("c10");
    htw_true->Draw("hist");
    htw_rec1->Draw("hist same");
    htw_rec2->Draw("hist same");
    htw_rec3->Draw("hist same");
    htw_rec4->Draw("hist same");
    leg1->Draw();

    //Print plots to file
    c1->Print("DEMP_reco.pdf[");
    c1->Print("DEMP_reco.pdf");
    c2->Print("DEMP_reco.pdf");
    c3->Print("DEMP_reco.pdf");
    c4->Print("DEMP_reco.pdf");
    c5->Print("DEMP_reco.pdf");
    c6->Print("DEMP_reco.pdf");
    c7->Print("DEMP_reco.pdf");
    c8->Print("DEMP_reco.pdf");
    c9->Print("DEMP_reco.pdf");
    c10->Print("DEMP_reco.pdf");
    c10->Print("DEMP_reco.pdf]");

}
