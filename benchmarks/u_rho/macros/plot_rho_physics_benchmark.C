#include "RiceStyle.h"
using namespace std;

int setbenchstatus(double eff){
	///////////// Set benchmark status!
        // create our test definition
        // test_tag
        common_bench::Test rho_reco_eff_test{
          {
            {"name", "rho_reconstruction_efficiency"},
            {"title", "rho Reconstruction Efficiency for rho -> pi+pi- in the B0"},
            {"description", "u-channel rho->pi+pi- reconstruction efficiency "
                       "when both pions should be within B0 acceptance"},
            {"quantity", "efficiency"},
            {"target", "0.9"}
          }
        };      //these 2 need to be consistent 
        double eff_target = 0.9;    //going to find a way to use the same variable

        if(eff<0 || eff>1){
          rho_reco_eff_test.error(-1);
        }else if(eff > eff_target){
          rho_reco_eff_test.pass(eff);
        }else{
          rho_reco_eff_test.fail(eff);
        }

        // write out our test data
        common_bench::write_test(rho_reco_eff_test, "rhorecoeff.json");
	return 0;
}

void plot_rho_physics_benchmark(TString filename="./benchmark_output/plot_combined.root"){
	Ssiz_t dotPosition = filename.Last('.');
	TString figure_directory = filename(0, dotPosition);
	figure_directory += "_figures";	
	
	TFile* file = new TFile(filename);
	TString vm_label="#rho^{0}";
	TString daug_label="#pi^{#plus}#pi^{#minus}";
	//t distribution
	TH1D* h_t_MC = (TH1D*) file->Get("h_t_MC");
	TH1D* h_t_REC = (TH1D*) file->Get("h_t_REC");
	TH1D* h_t_trk_REC = (TH1D*) file->Get("h_t_trk_REC");
	TH1D* h_t_combo_REC = (TH1D*) file->Get("h_t_combo_REC");
	//mass distribution
        TH1D* h_VM_mass_MC = (TH1D*) file->Get("h_VM_mass_MC");
        TH1D* h_VM_mass_REC = (TH1D*) file->Get("h_VM_mass_REC");
	TH1D* h_VM_mass_REC_justpions = (TH1D*) file->Get("h_VM_mass_REC_justpions");
	//mass distribution within B0
	TH1D* h_VM_mass_MC_etacut = (TH1D*) file->Get("h_VM_mass_MC_etacut");
        TH1D* h_VM_mass_REC_etacut = (TH1D*) file->Get("h_VM_mass_REC_etacut");
	//dN/du distribution
        TH1D* h_dNdu_MC = (TH1D*) file->Get("h_u_MC");
        TH1D* h_dNdu_REC = (TH1D*) file->Get("h_u_REC");
	TH1D* h_dNdu_REC_justpions = (TH1D*) file->Get("h_u_REC_justpions");
	//efficiencies
	TProfile2D* h_effEtaPtPi = (TProfile2D*) file->Get("h_effEtaPtPi");
        TProfile2D* h_effEtaPtPip = (TProfile2D*) file->Get("h_effEtaPtPip"); 
        TProfile2D* h_effEtaPtPim = (TProfile2D*) file->Get("h_effEtaPtPim"); 
        TProfile2D* h_effPhiEtaPi = (TProfile2D*) file->Get("h_effPhiEtaPi"); 
        TProfile2D* h_effPhiEtaPip = (TProfile2D*) file->Get("h_effPhiEtaPip"); 
        TProfile2D* h_effPhiEtaPim = (TProfile2D*) file->Get("h_effPhiEtaPim"); 	
	//reco quality
	TH2D* h_RecoMomPi = (TH2D*) file->Get("h_RecoMomPi");
        TH2D* h_RecoMomPip= (TH2D*) file->Get("h_RecoMomPip");
        TH2D* h_RecoMomPim= (TH2D*) file->Get("h_RecoMomPim");
        TH2D* h_RecoTransMomPi = (TH2D*) file->Get("h_RecoTransMomPi");
        TH2D* h_RecoTransMomPip= (TH2D*) file->Get("h_RecoTransMomPip");
        TH2D* h_RecoTransMomPim= (TH2D*) file->Get("h_RecoTransMomPim");


	TCanvas* c1 = new TCanvas("c1","c1",1,1,600,600);
	gPad->SetLogy(1);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.15);
	gPad->SetRightMargin(0.01);
	TH1D* base1 = makeHist("base1", "", "|#it{t} | (GeV^{2})", "dN/d|#it{t} | (GeV^{-2}) ", 100,0,5.0,kBlack);
	base1->GetYaxis()->SetRangeUser(8e-2, 8e5);
	base1->GetXaxis()->SetTitleColor(kBlack);
	fixedFontHist1D(base1,1.,1.2);
	base1->GetYaxis()->SetTitleSize(base1->GetYaxis()->GetTitleSize()*1.5);
	base1->GetXaxis()->SetTitleSize(base1->GetXaxis()->GetTitleSize()*1.5);
	base1->GetYaxis()->SetLabelSize(base1->GetYaxis()->GetLabelSize()*1.5);
	base1->GetXaxis()->SetLabelSize(base1->GetXaxis()->GetLabelSize()*1.5);
	base1->GetXaxis()->SetNdivisions(4,4,0);
	base1->GetYaxis()->SetNdivisions(5,5,0);
	base1->Draw();

	h_t_MC->Draw("same");

	h_t_REC->SetMarkerStyle(20);
	h_t_REC->Draw("PEsame");

	h_t_trk_REC->SetFillColorAlpha(kBlue,0.4);
    h_t_trk_REC->SetFillStyle(1001);
	h_t_trk_REC->SetMarkerStyle(24);
	h_t_trk_REC->SetMarkerColor(kBlue);
	// h_t_trk_REC->Draw("PE3same");

	h_t_combo_REC->SetFillColorAlpha(kRed,0.4);
    h_t_combo_REC->SetFillStyle(1001);
	h_t_combo_REC->SetMarkerStyle(24);
	h_t_combo_REC->SetMarkerColor(kRed);
	// h_t_combo_REC->Draw("PE3same");

	TLatex* r42 = new TLatex(0.18, 0.91, "ep 10#times100 GeV");
	r42->SetNDC();
	r42->SetTextSize(22);
	r42->SetTextFont(43);
	r42->SetTextColor(kBlack);
	r42->Draw("same");

	TLatex* r43 = new TLatex(0.9,0.91, "EPIC");
	r43->SetNDC();
	r43->SetTextSize(0.04);
	r43->Draw("same");

	TLatex* r44 = new TLatex(0.53, 0.78, "10^{-3}<Q^{2}<10 GeV^{2}, W>2 GeV");
	r44->SetNDC();
	r44->SetTextSize(20);
	r44->SetTextFont(43);
	r44->SetTextColor(kBlack);
	r44->Draw("same");

	TLatex* r44_2 = new TLatex(0.5, 0.83, ""+vm_label+" #rightarrow "+daug_label+" eSTARlight");
	r44_2->SetNDC();
	r44_2->SetTextSize(30);
	r44_2->SetTextFont(43);
	r44_2->SetTextColor(kBlack);
	r44_2->Draw("same");

	TLegend *w7 = new TLegend(0.58,0.68,0.93,0.76);
	w7->SetLineColor(kWhite);
	w7->SetFillColor(0);
	w7->SetTextSize(17);
	w7->SetTextFont(45);
	w7->AddEntry(h_t_MC, "eSTARlight "+vm_label+" MC ", "L");
	w7->AddEntry(h_t_REC, "eSTARlight "+vm_label+" RECO ", "P");
	w7->Draw("same");

	//c1->Print("./benchmark_output/figures/benchmark_rho_dsigmadt.pdf");
	TString figure1name = figure_directory+"/benchmark_rho_dsigmadt.pdf";
        c1->Print(figure1name);

	TCanvas* c2 = new TCanvas("c2","c2",1,1,600,600);
        gPad->SetTicks();
        gPad->SetLeftMargin(0.18);
        gPad->SetBottomMargin(0.18);
        gPad->SetRightMargin(0.01);
        TH1D* base2 = makeHist("base2", "", "#pi^{#plus}#pi^{#minus} inv. mass (GeV)", "counts", 100,0.05,2.05,kBlack);
        base2->GetYaxis()->SetRangeUser(0.5, 1.2*(h_VM_mass_MC->GetMaximum()));
        base2->GetXaxis()->SetTitleColor(kBlack);
        fixedFontHist1D(base2,1.,1.2);
        base2->GetYaxis()->SetTitleSize(base2->GetYaxis()->GetTitleSize()*1.5);
        base2->GetXaxis()->SetTitleSize(base2->GetXaxis()->GetTitleSize()*1.5);
        base2->GetYaxis()->SetLabelSize(base2->GetYaxis()->GetLabelSize()*1.5);
        base2->GetXaxis()->SetLabelSize(base2->GetXaxis()->GetLabelSize()*1.5);
        base2->GetXaxis()->SetNdivisions(4,4,0);
        base2->GetYaxis()->SetNdivisions(5,5,0);
	base2->GetYaxis()->SetTitleOffset(1.3);
        base2->Draw();

	TH1D* h_VM_mass_REC_justprotons = (TH1D*)h_VM_mass_REC->Clone("h_VM_mass_REC_justprotons");
	for(int ibin=1; ibin<h_VM_mass_REC_justprotons->GetNbinsX(); ibin++){
	  h_VM_mass_REC_justprotons->SetBinContent(ibin,h_VM_mass_REC_justprotons->GetBinContent(ibin) - h_VM_mass_REC_justpions->GetBinContent(ibin));
	}

	h_VM_mass_MC->SetFillColorAlpha(kBlack,0.2);
        h_VM_mass_REC->SetFillColorAlpha(kMagenta,0.2);
	//h_VM_mass_REC_justpions->SetFillColorAlpha(kViolet+10,0.2);
	h_VM_mass_MC->SetLineColor(kBlack);
        h_VM_mass_REC->SetLineColor(kMagenta);
	h_VM_mass_REC_justpions->SetLineColor(kViolet+10);
	h_VM_mass_REC_justprotons->SetLineColor(kRed);
	h_VM_mass_MC->SetLineWidth(2);
        h_VM_mass_REC->SetLineWidth(2);
	h_VM_mass_REC_justpions->SetLineWidth(2);
	h_VM_mass_REC_justprotons->SetLineWidth(2);

	h_VM_mass_REC->Scale(3.0);
	h_VM_mass_REC_justpions->Scale(3.0);
	h_VM_mass_REC_justprotons->Scale(3.0);

        h_VM_mass_MC->Draw("HIST E same");
        h_VM_mass_REC->Draw("HIST E same");
	h_VM_mass_REC_justpions->Draw("HIST same");
        h_VM_mass_REC_justprotons->Draw("HIST same");

        r42->Draw("same");
        r43->Draw("same");
        r44->Draw("same");
        r44_2->Draw("same");

        TLegend *w8 = new TLegend(0.58,0.63,0.93,0.76);
        w8->SetLineColor(kWhite);
        w8->SetFillColor(0);
        w8->SetTextSize(17);
        w8->SetTextFont(45);
        w8->AddEntry(h_VM_mass_MC, ""+vm_label+" MC ", "L");
        w8->AddEntry(h_VM_mass_REC, vm_label+" reco.#times3", "L");
        w8->AddEntry(h_VM_mass_REC_justpions, vm_label+" reco.#times3 (#pi^{#minus}#pi^{#plus})", "L");
        w8->AddEntry(h_VM_mass_REC_justprotons, vm_label+" reco.#times3 (#pi^{#minus}p)", "L");
        w8->Draw("same");

	TString figure2name = figure_directory+"/benchmark_rho_mass.pdf";
        c2->Print(figure2name);

	///////////////////// Figure 3
        TCanvas* c3 = new TCanvas("c3","c3",1,1,600,600);
        gPad->SetTicks();
	gPad->SetLogy(1);
        gPad->SetLeftMargin(0.18);
        gPad->SetBottomMargin(0.18);
        gPad->SetRightMargin(0.01);
        TH1D* base3 = makeHist("base3", "", "-#it{u} (GeV^{2})", "dN/d#it{u} (GeV^{-2} scaled)", 100,-0.25,3.05,kBlack);
        base3->GetYaxis()->SetRangeUser(0.5, 100*(h_dNdu_MC->GetMaximum()));
        base3->GetXaxis()->SetTitleColor(kBlack);
        fixedFontHist1D(base3,1.,1.2);
        base3->GetYaxis()->SetTitleSize(base3->GetYaxis()->GetTitleSize()*1.5);
        base3->GetXaxis()->SetTitleSize(base3->GetXaxis()->GetTitleSize()*1.5);
        base3->GetYaxis()->SetLabelSize(base3->GetYaxis()->GetLabelSize()*1.5);
        base3->GetXaxis()->SetLabelSize(base3->GetXaxis()->GetLabelSize()*1.5);
	base3->GetYaxis()->SetTitleOffset(1.2);
        base3->GetXaxis()->SetNdivisions(4,4,0);
        base3->GetYaxis()->SetNdivisions(5,5,0);
        base3->Draw();

        h_dNdu_MC->SetLineColor(kBlack);
        h_dNdu_REC->SetLineColor(kMagenta);
        h_dNdu_REC_justpions->SetLineColor(kBlue);
        h_dNdu_MC->SetLineWidth(2);
        h_dNdu_REC->SetLineWidth(2);
        h_dNdu_REC_justpions->SetLineWidth(2);

        h_dNdu_MC->Draw("HIST E same");
        h_dNdu_REC->Draw("HIST E same");
        h_dNdu_REC_justpions->Draw("HIST E same");

        r42->Draw("same");
        r43->Draw("same");
        r44->Draw("same");
        r44_2->Draw("same");

	TF1* fit_mc    = new TF1("fit_mc", "[0]*exp([1]*x)",0.2, 1.2);
        TF1* fit_rec   = new TF1("fit_rec", "[0]*exp([1]*x)",0.2, 1.2);
        TF1* fit_recpi = new TF1("fit_recpi", "[0]*exp([1]*x)",0.2, 1.2);
	fit_mc->SetLineColor(kBlack);
        fit_rec->SetLineColor(kMagenta);
        fit_recpi->SetLineColor(kBlue);
        fit_mc->SetLineWidth(3);
        fit_rec->SetLineWidth(3);
        fit_recpi->SetLineWidth(3);
	TFitResultPtr r1 =  h_dNdu_MC->Fit("fit_mc","RLMSIN");
        TFitResultPtr r2 =  h_dNdu_REC->Fit("fit_rec","RLMSIN");
        TFitResultPtr r3 =  h_dNdu_REC_justpions->Fit("fit_recpi","RLMSIN");
    	double alpha_mc    = r1->Parameter(1);
        double alpha_rec   = r2->Parameter(1);
        double alpha_recpi = r3->Parameter(1);
	fit_mc->Draw("SAME");
        fit_rec->Draw("SAME");
        fit_recpi->Draw("SAME");

	TH1I* hwhite = new TH1I("hwhite","hwhite",1,0,1);
	hwhite->SetLineColor(kWhite);
        TLegend *w9 = new TLegend(0.53,0.61,0.86,0.76);
        w9->SetLineColor(kWhite);
        w9->SetFillColor(0);
        w9->SetTextSize(17);
        w9->SetTextFont(45);
        w9->AddEntry(hwhite, "~exp[#alpha(-#it{u})]", "L");
        w9->AddEntry(h_dNdu_MC, vm_label+Form(" MC, #alpha=%.1fGeV^{-2}",alpha_mc), "L");
        w9->AddEntry(h_dNdu_REC, vm_label+Form(" reco. #alpha=%.1fGeV^{-2}",alpha_rec), "L");
        w9->AddEntry(h_dNdu_REC_justpions, vm_label+Form(" reco. (#pi^{#minus}#pi^{#plus}), #alpha=%.1fGeV^{-2}",alpha_recpi), "L");
        w9->Draw("same");

        TString figure3name = figure_directory+"/benchmark_rho_dNdu.pdf";
        c3->Print(figure3name);

	///////////////////// Figure 4
        TCanvas* c4 = new TCanvas("c4","c4",1,1,600,600);
        gPad->SetTicks();
        gPad->SetLeftMargin(0.18);
        gPad->SetBottomMargin(0.18);
        gPad->SetRightMargin(0.01);
        TH1D* base4 = makeHist("base4", "", "#pi^{#plus}#pi^{#minus} inv. mass (GeV)", "counts", 100,0.05,2.05,kBlack);
        base4->GetYaxis()->SetRangeUser(0.5, 1.2*(h_VM_mass_MC_etacut->GetMaximum()));
        base4->GetXaxis()->SetTitleColor(kBlack);
        fixedFontHist1D(base4,1.,1.2);
        base4->GetYaxis()->SetTitleSize(base4->GetYaxis()->GetTitleSize()*1.5);
        base4->GetXaxis()->SetTitleSize(base4->GetXaxis()->GetTitleSize()*1.5);
        base4->GetYaxis()->SetLabelSize(base4->GetYaxis()->GetLabelSize()*1.5);
        base4->GetXaxis()->SetLabelSize(base4->GetXaxis()->GetLabelSize()*1.5);
        base4->GetXaxis()->SetNdivisions(4,4,0);
        base4->GetYaxis()->SetNdivisions(5,5,0);
        base4->GetYaxis()->SetTitleOffset(1.3);
        base4->Draw();

        h_VM_mass_MC_etacut->SetFillColorAlpha(kBlack,0.2);
        h_VM_mass_REC_etacut->SetFillColorAlpha(kMagenta,0.2);
        h_VM_mass_MC_etacut->SetLineColor(kBlack);
        h_VM_mass_REC_etacut->SetLineColor(kMagenta);
        h_VM_mass_MC_etacut->SetLineWidth(2);
        h_VM_mass_REC_etacut->SetLineWidth(2);

        h_VM_mass_MC_etacut->Draw("HIST E same");
        h_VM_mass_REC_etacut->Draw("HIST E same");

	double minbineff = h_VM_mass_MC_etacut->FindBin(0.6);
        double maxbineff = h_VM_mass_MC_etacut->FindBin(1.0);
	double thiseff = 100.0*(1.0*h_VM_mass_REC_etacut->Integral(minbineff,maxbineff))/(1.0*h_VM_mass_MC_etacut->Integral(minbineff,maxbineff));

        r42->Draw("same");
        r43->Draw("same");
        r44->Draw("same");
        r44_2->Draw("same");

        TLegend *w10 = new TLegend(0.58,0.62,0.93,0.7);
        w10->SetLineColor(kWhite);
        w10->SetFillColor(0);
        w10->SetTextSize(17);
        w10->SetTextFont(45);
        w10->AddEntry(h_VM_mass_MC_etacut, vm_label+" MC ", "L");
        w10->AddEntry(h_VM_mass_REC_etacut, vm_label+" reco. (#pi^{#minus}#pi^{#plus})", "L");
        w10->Draw("same");

        TLatex* anglelabel = new TLatex(0.59, 0.73, "9<#theta_{#pi^{#pm},MC}<13 mrad");
        anglelabel->SetNDC();
        anglelabel->SetTextSize(20);
        anglelabel->SetTextFont(43);
        anglelabel->SetTextColor(kBlack);
        anglelabel->Draw("same");

        TLatex* efflabel = new TLatex(0.59, 0.55, "reco. eff (0.6<m^{2}<1 GeV^{2})");
        efflabel->SetNDC();
        efflabel->SetTextSize(20);
        efflabel->SetTextFont(43);
        efflabel->SetTextColor(kBlack);
        efflabel->Draw("same");
        TLatex* effnlabel = new TLatex(0.59, 0.51, Form("          = %.0f%%",thiseff));
        effnlabel->SetNDC();
        effnlabel->SetTextSize(20);
        effnlabel->SetTextFont(43);
        effnlabel->SetTextColor(kBlack);
        effnlabel->Draw("same");

        TString figure4name = figure_directory+"/benchmark_rho_mass_cuts.pdf";
        c4->Print(figure4name);

        ///////////////////// Figure 5
        TCanvas* c5 = new TCanvas("c5","c5",1,1,700,560);
	TPad* p5 = new TPad("p5","Pad5",0.,0.,1.,1.);
    	p5->Divide(3,2,0,0);
	p5->Draw();
	gStyle->SetPalette(kBlueRedYellow);
	gStyle->SetOptStat(0);	

	h_effEtaPtPi  ->GetXaxis()->SetLabelSize(h_effEtaPtPi  ->GetXaxis()->GetLabelSize()*1.8);
        h_effEtaPtPip ->GetXaxis()->SetLabelSize(h_effEtaPtPip ->GetXaxis()->GetLabelSize()*1.8);
        h_effEtaPtPim ->GetXaxis()->SetLabelSize(h_effEtaPtPim ->GetXaxis()->GetLabelSize()*1.8);
        h_effEtaPtPi  ->GetYaxis()->SetLabelSize(h_effEtaPtPi  ->GetYaxis()->GetLabelSize()*1.8);
        h_effEtaPtPim ->GetZaxis()->SetLabelSize(h_effEtaPtPim ->GetZaxis()->GetLabelSize()*1.8);
        h_effEtaPtPim ->GetZaxis()->SetTitleSize(h_effEtaPtPim ->GetZaxis()->GetTitleSize()*1.8);
        h_effPhiEtaPi ->GetXaxis()->SetLabelSize(h_effPhiEtaPi ->GetXaxis()->GetLabelSize()*1.8);
        h_effPhiEtaPip->GetXaxis()->SetLabelSize(h_effPhiEtaPip->GetXaxis()->GetLabelSize()*1.8);
        h_effPhiEtaPim->GetXaxis()->SetLabelSize(h_effPhiEtaPim->GetXaxis()->GetLabelSize()*1.8);
        h_effPhiEtaPi ->GetYaxis()->SetLabelSize(h_effPhiEtaPi ->GetYaxis()->GetLabelSize()*1.8);
        h_effPhiEtaPim->GetZaxis()->SetLabelSize(h_effPhiEtaPim->GetZaxis()->GetLabelSize()*1.8);
        h_effPhiEtaPim->GetZaxis()->SetTitleSize(h_effPhiEtaPim->GetZaxis()->GetTitleSize()*1.8);

	fixedFontHist1D(h_effEtaPtPi,1.,1.2);
        fixedFontHist1D(h_effEtaPtPip,1.,1.2);
        fixedFontHist1D(h_effEtaPtPim,1.,1.2);
        fixedFontHist1D(h_effPhiEtaPi,1.,1.2);
        fixedFontHist1D(h_effPhiEtaPip,1.,1.2);
        fixedFontHist1D(h_effPhiEtaPim,1.,1.2);

    	p5->cd(1);
    	TVirtualPad* p51 = p5->cd(1);
    	p51->SetTopMargin(0.08);
    	p51->SetRightMargin(0);
    	p51->SetLeftMargin(0.21);
    	p51->SetBottomMargin(0.2);
        h_effEtaPtPi->GetXaxis()->SetRangeUser(3.9,6.05);
        h_effEtaPtPi->GetYaxis()->SetRangeUser(0,1.7);
        h_effEtaPtPi->GetZaxis()->SetRangeUser(0,1);
        h_effEtaPtPi->GetXaxis()->SetNdivisions(5);
        h_effEtaPtPi->GetYaxis()->SetNdivisions(5);
        h_effEtaPtPi->SetContour(99);
	h_effEtaPtPi->Draw("COLZ");
        TLatex* pilabel = new TLatex(0.81, 0.75, "#pi^{#pm}");
        pilabel->SetNDC();
        pilabel->SetTextSize(40);
        pilabel->SetTextFont(43);
        pilabel->SetTextColor(kBlack);
        pilabel->Draw("same");
        TLatex* r44fig5c = new TLatex(0.21, 0.93, "ep 10#times100 GeV         #rho^{0}#rightarrow#pi^{#plus}#pi^{#minus}");
        r44fig5c->SetNDC();
        r44fig5c->SetTextSize(15);
        r44fig5c->SetTextFont(43);
        r44fig5c->SetTextColor(kBlack);
        r44fig5c->Draw("same");

        p5->cd(2);
        TVirtualPad* p52 = p5->cd(2);
        p52->SetTopMargin(0.08);
        p52->SetRightMargin(0);
        p52->SetLeftMargin(0);
        p52->SetBottomMargin(0.2);
        h_effEtaPtPip->GetXaxis()->SetRangeUser(4.05,6.05);
        h_effEtaPtPip->GetYaxis()->SetRangeUser(0,1.7);
        h_effEtaPtPip->GetZaxis()->SetRangeUser(0,1);
        h_effEtaPtPip->GetXaxis()->SetNdivisions(5);
        h_effEtaPtPip->GetYaxis()->SetNdivisions(5);
        h_effEtaPtPip->SetContour(99);
        h_effEtaPtPip->Draw("COLZ");
        TLatex* piplabel = new TLatex(0.81, 0.75, "#pi^{#plus}");
        piplabel->SetNDC();
        piplabel->SetTextSize(40);
        piplabel->SetTextFont(43);
        piplabel->SetTextColor(kBlack);
        piplabel->Draw("same");
        TLatex* r44fig5a = new TLatex(0.01, 0.93, "eSTARlight        10^{-3}<Q^{2}<10 GeV^{2}");
        r44fig5a->SetNDC();
        r44fig5a->SetTextSize(15);
        r44fig5a->SetTextFont(43);
        r44fig5a->SetTextColor(kBlack);
        r44fig5a->Draw("same");

        p5->cd(3);
        TVirtualPad* p53 = p5->cd(3);
        p53->SetTopMargin(0.08);
        p53->SetRightMargin(0.2);
        p53->SetLeftMargin(0);
        p53->SetBottomMargin(0.2);
	h_effEtaPtPim->SetTitle(";#eta;;efficiency");
        h_effEtaPtPim->GetXaxis()->SetRangeUser(4.05,6.05);
        h_effEtaPtPim->GetYaxis()->SetRangeUser(0,1.7);
        h_effEtaPtPim->GetZaxis()->SetRangeUser(0,1);
        h_effEtaPtPim->GetXaxis()->SetNdivisions(5);
        h_effEtaPtPim->GetYaxis()->SetNdivisions(5);
        h_effEtaPtPim->SetContour(99);
        h_effEtaPtPim->Draw("COLZ");
        TLatex* pimlabel = new TLatex(0.62, 0.75, "#pi^{#minus}");
        pimlabel->SetNDC();
        pimlabel->SetTextSize(40);
        pimlabel->SetTextFont(43);
        pimlabel->SetTextColor(kBlack);
        pimlabel->Draw("same");
        TLatex* r43fig5 = new TLatex(0.68,0.93, "EPIC");
        r43fig5->SetNDC();
        r43fig5->SetTextSize(0.07);
        r43fig5->Draw("same");
        TLatex* r44fig5b = new TLatex(0.01, 0.93, "W>2 GeV");
        r44fig5b->SetNDC();
        r44fig5b->SetTextSize(15);
        r44fig5b->SetTextFont(43);
        r44fig5b->SetTextColor(kBlack);
        r44fig5b->Draw("same");

        p5->cd(4);
        TVirtualPad* p54 = p5->cd(4);
        p54->SetTopMargin(0.05);
        p54->SetRightMargin(0);
        p54->SetLeftMargin(0.2);
        p54->SetBottomMargin(0.21);
        h_effPhiEtaPi->GetXaxis()->SetRangeUser(0,6.2);
        h_effPhiEtaPi->GetYaxis()->SetRangeUser(4,6);
        h_effPhiEtaPi->GetZaxis()->SetRangeUser(0,1);
        h_effPhiEtaPi->GetXaxis()->SetNdivisions(8);
        h_effPhiEtaPi->GetYaxis()->SetNdivisions(5);
        h_effPhiEtaPi->SetContour(99);
        h_effPhiEtaPi->Draw("COLZ");
        TLatex* pilabela = new TLatex(0.3, 0.82, "#pi^{#pm}");
        TLatex* pilabelb = new TLatex(0.5, 0.84, "(p_{T}>0.2 GeV/c)");
        pilabela->SetNDC();
        pilabelb->SetNDC();
        pilabela->SetTextSize(40);
        pilabelb->SetTextSize(15);
        pilabela->SetTextFont(43);
        pilabelb->SetTextFont(43);
        pilabela->SetTextColor(kWhite);
        pilabelb->SetTextColor(kWhite);
        pilabela->Draw("same");
        pilabelb->Draw("same");

        p5->cd(5);
        TVirtualPad* p55 = p5->cd(5);
        p55->SetTopMargin(0.05);
        p55->SetRightMargin(0);
        p55->SetLeftMargin(0);
        p55->SetBottomMargin(0.2);
        h_effPhiEtaPip->GetXaxis()->SetRangeUser(0.15,6.2);
        h_effPhiEtaPip->GetYaxis()->SetRangeUser(4,6);
        h_effPhiEtaPip->GetZaxis()->SetRangeUser(0,1);
        h_effPhiEtaPip->GetXaxis()->SetNdivisions(8);
        h_effPhiEtaPip->GetYaxis()->SetNdivisions(5);
        h_effPhiEtaPip->SetContour(99);
        h_effPhiEtaPip->Draw("COLZ");
        TLatex* piplabela = new TLatex(0.2, 0.82, "#pi^{#plus}");
        TLatex* piplabelb = new TLatex(0.4, 0.84, "(p_{T}>0.2 GeV/c)");
        piplabela->SetNDC();
        piplabelb->SetNDC();
        piplabela->SetTextSize(40);
        piplabelb->SetTextSize(15);
        piplabela->SetTextFont(43);
        piplabelb->SetTextFont(43);
        piplabela->SetTextColor(kWhite);
        piplabelb->SetTextColor(kWhite);
        piplabela->Draw("same");
        piplabelb->Draw("same");

        p5->cd(6);
        TVirtualPad* p56 = p5->cd(6);
        p56->SetTopMargin(0.05);
        p56->SetRightMargin(0.2);
        p56->SetLeftMargin(0);
        p56->SetBottomMargin(0.2);
        h_effPhiEtaPim->SetTitle(";#phi (rad);;efficiency");
        h_effPhiEtaPim->GetXaxis()->SetRangeUser(0.15,6.2);
        h_effPhiEtaPim->GetYaxis()->SetRangeUser(4,6);
        h_effPhiEtaPim->GetZaxis()->SetRangeUser(0,1);
        h_effPhiEtaPim->GetXaxis()->SetNdivisions(8);
        h_effPhiEtaPim->GetYaxis()->SetNdivisions(5);
        h_effPhiEtaPim->SetContour(99);
        h_effPhiEtaPim->Draw("COLZ");
        TLatex* pimlabela = new TLatex(0.1, 0.82, "#pi^{#minus}");
        TLatex* pimlabelb = new TLatex(0.25, 0.84, "(p_{T}>0.2 GeV/c)");
        pimlabela->SetNDC();
        pimlabelb->SetNDC();
        pimlabela->SetTextSize(40);
        pimlabelb->SetTextSize(15);
        pimlabela->SetTextFont(43);
        pimlabelb->SetTextFont(43);
        pimlabela->SetTextColor(kWhite);
        pimlabelb->SetTextColor(kWhite);
        pimlabela->Draw("same");
        pimlabelb->Draw("same");

        TString figure5name = figure_directory+"/benchmark_rho_efficiencies.pdf";
        c5->Print(figure5name);

        ///////////////////// Figure 6
        TCanvas* c6 = new TCanvas("c6","c6",1,1,700,560);
	TPad* p6 = new TPad("p6","Pad5",0.,0.,1.,1.);
    	p6->Divide(3,2,0,0);
	p6->Draw();
	gStyle->SetPalette(kBlueRedYellow);
	gStyle->SetOptStat(0);	

	h_RecoMomPi  ->GetXaxis()->SetLabelSize(h_RecoMomPi  ->GetXaxis()->GetLabelSize()*1.8);
        h_RecoMomPip ->GetXaxis()->SetLabelSize(h_RecoMomPip ->GetXaxis()->GetLabelSize()*1.8);
        h_RecoMomPim ->GetXaxis()->SetLabelSize(h_RecoMomPim ->GetXaxis()->GetLabelSize()*1.8);
        h_RecoMomPi  ->GetYaxis()->SetLabelSize(h_RecoMomPi  ->GetYaxis()->GetLabelSize()*1.8);
        h_RecoMomPim ->GetZaxis()->SetLabelSize(h_RecoMomPim ->GetZaxis()->GetLabelSize()*1.8);
        h_RecoMomPim ->GetZaxis()->SetTitleSize(h_RecoMomPim ->GetZaxis()->GetTitleSize()*1.8);
        h_RecoTransMomPi ->GetXaxis()->SetLabelSize(h_RecoTransMomPi ->GetXaxis()->GetLabelSize()*1.8);
        h_RecoTransMomPip->GetXaxis()->SetLabelSize(h_RecoTransMomPip->GetXaxis()->GetLabelSize()*1.8);
        h_RecoTransMomPim->GetXaxis()->SetLabelSize(h_RecoTransMomPim->GetXaxis()->GetLabelSize()*1.8);
        h_RecoTransMomPi ->GetYaxis()->SetLabelSize(h_RecoTransMomPi ->GetYaxis()->GetLabelSize()*1.8);
        h_RecoTransMomPim->GetZaxis()->SetLabelSize(h_RecoTransMomPim->GetZaxis()->GetLabelSize()*1.8);
        h_RecoTransMomPim->GetZaxis()->SetTitleSize(h_RecoTransMomPim->GetZaxis()->GetTitleSize()*1.8);

	fixedFontHist1D(h_RecoMomPi,1.,1.2);
        fixedFontHist1D(h_RecoMomPip,1.,1.2);
        fixedFontHist1D(h_RecoMomPim,1.,1.2);
        fixedFontHist1D(h_RecoTransMomPi,1.,1.2);
        fixedFontHist1D(h_RecoTransMomPip,1.,1.2);
        fixedFontHist1D(h_RecoTransMomPim,1.,1.2);

	double maxz = h_RecoMomPi->GetMaximum();
        double maxzt = h_RecoTransMomPi->GetMaximum();

    	p6->cd(1);
    	TVirtualPad* p61 = p6->cd(1);
	p61->SetLogz();
    	p61->SetTopMargin(0.08);
    	p61->SetRightMargin(0);
    	p61->SetLeftMargin(0.21);
    	p61->SetBottomMargin(0.2);
        h_RecoMomPi->GetXaxis()->SetRangeUser(0,99);
        h_RecoMomPi->GetYaxis()->SetRangeUser(0,99);
        h_RecoMomPi->GetZaxis()->SetRangeUser(1,maxz);
        h_RecoMomPi->GetXaxis()->SetNdivisions(5);
        h_RecoMomPi->GetYaxis()->SetNdivisions(5);
        h_RecoMomPi->SetContour(99);
	h_RecoMomPi->Draw("COLZ");
        TBox* box1 = new TBox(7,70,32,95);
	box1->SetLineColor(kBlack);
        box1->SetFillColor(kWhite);
	box1->Draw("l same");
        TLatex* pilabelz = new TLatex(0.3, 0.75, "#pi^{#pm}");
        pilabelz->SetNDC();
        pilabelz->SetTextSize(40);
        pilabelz->SetTextFont(43);
        pilabelz->SetTextColor(kBlack);
        pilabelz->Draw("same");
        r44fig5c->Draw("same");

        p6->cd(2);
        TVirtualPad* p62 = p6->cd(2);
        p62->SetLogz();
        p62->SetTopMargin(0.08);
        p62->SetRightMargin(0);
        p62->SetLeftMargin(0);
        p62->SetBottomMargin(0.2);
        h_RecoMomPip->GetXaxis()->SetRangeUser(1,99);
        h_RecoMomPip->GetYaxis()->SetRangeUser(0,99);
        h_RecoMomPip->GetZaxis()->SetRangeUser(1,maxz);
        h_RecoMomPip->GetXaxis()->SetNdivisions(5);
        h_RecoMomPip->GetYaxis()->SetNdivisions(5);
        h_RecoMomPip->SetContour(99);
        h_RecoMomPip->Draw("COLZ");
        TBox* box2 = new TBox(9,70,30,95);
        box2->SetLineColor(kBlack);
        box2->SetFillColor(kWhite);
        box2->Draw("l same");
        TLatex* piplabelz = new TLatex(0.11, 0.75, "#pi^{#plus}");
        piplabelz->SetNDC();
        piplabelz->SetTextSize(40);
        piplabelz->SetTextFont(43);
        piplabelz->SetTextColor(kBlack);
        piplabelz->Draw("same");
        r44fig5a->Draw("same");

        p6->cd(3);
        TVirtualPad* p63 = p6->cd(3);
        p63->SetLogz();
        p63->SetTopMargin(0.08);
        p63->SetRightMargin(0.2);
        p63->SetLeftMargin(0);
        p63->SetBottomMargin(0.2);
	h_RecoMomPim->SetTitle(";p (GeV/c) MC;;counts");
        h_RecoMomPim->GetXaxis()->SetRangeUser(1,99);
        h_RecoMomPim->GetYaxis()->SetRangeUser(0,99);
        h_RecoMomPim->GetZaxis()->SetRangeUser(1,maxz);
        h_RecoMomPim->GetXaxis()->SetNdivisions(5);
        h_RecoMomPim->GetYaxis()->SetNdivisions(5);
        h_RecoMomPim->SetContour(99);
        h_RecoMomPim->Draw("COLZ");
        TBox* box3 = new TBox(12,70,40,95);
        box3->SetLineColor(kBlack);
        box3->SetFillColor(kWhite);
        box3->Draw("l same");
        TLatex* pimlabelz = new TLatex(0.12, 0.75, "#pi^{#minus}");
        pimlabelz->SetNDC();
        pimlabelz->SetTextSize(40);
        pimlabelz->SetTextFont(43);
        pimlabelz->SetTextColor(kBlack);
        pimlabelz->Draw("same");
        r43fig5->Draw("same");
        r44fig5b->Draw("same");

        p6->cd(4);
        TVirtualPad* p64 = p6->cd(4);
        p64->SetLogz();
        p64->SetTopMargin(0.05);
        p64->SetRightMargin(0);
        p64->SetLeftMargin(0.21);
        p64->SetBottomMargin(0.2);
        h_RecoTransMomPi->GetXaxis()->SetRangeUser(0.01,1.49);
        h_RecoTransMomPi->GetYaxis()->SetRangeUser(0,1.49);
        h_RecoTransMomPi->GetZaxis()->SetRangeUser(1,maxzt);
        h_RecoTransMomPi->GetXaxis()->SetNdivisions(5);
        h_RecoTransMomPi->GetYaxis()->SetNdivisions(5);
        h_RecoTransMomPi->SetContour(99);
        h_RecoTransMomPi->Draw("COLZ");
        TBox* box4 = new TBox(0.15,1.15,0.5,1.45);
        box4->SetLineColor(kBlack);
        box4->SetFillColor(kWhite);
        box4->Draw("l same");
        TLatex* pilabelaa = new TLatex(0.3, 0.8, "#pi^{#pm}");
        pilabelaa->SetNDC();
        pilabelaa->SetTextSize(40);
        pilabelaa->SetTextFont(43);
        pilabelaa->SetTextColor(kBlack);
        pilabelaa->Draw("same");

        p6->cd(5);
        TVirtualPad* p65 = p6->cd(5);
        p65->SetLogz();
        p65->SetTopMargin(0.05);
        p65->SetRightMargin(0);
        p65->SetLeftMargin(0);
        p65->SetBottomMargin(0.2);
        h_RecoTransMomPip->GetXaxis()->SetRangeUser(0.02,1.49);
        h_RecoTransMomPip->GetYaxis()->SetRangeUser(0,1.49);
        h_RecoTransMomPip->GetZaxis()->SetRangeUser(1,maxzt);
        h_RecoTransMomPip->GetXaxis()->SetNdivisions(5);
        h_RecoTransMomPip->GetYaxis()->SetNdivisions(5);
        h_RecoTransMomPip->SetContour(99);
        h_RecoTransMomPip->Draw("COLZ");
        TBox* box5 = new TBox(0.1,1.15,0.4,1.45);
        box5->SetLineColor(kBlack);
        box5->SetFillColor(kWhite);
        box5->Draw("l same");
        TLatex* piplabelaz = new TLatex(0.1, 0.8, "#pi^{#plus}");
        piplabelaz->SetNDC();
        piplabelaz->SetTextSize(40);
        piplabelaz->SetTextFont(43);
        piplabelaz->SetTextColor(kBlack);
        piplabelaz->Draw("same");

        p6->cd(6);
        TVirtualPad* p66 = p6->cd(6);
        p66->SetLogz();
        p66->SetTopMargin(0.05);
        p66->SetRightMargin(0.2);
        p66->SetLeftMargin(0);
        p66->SetBottomMargin(0.2);
        h_RecoTransMomPim->SetTitle(";p_{T} (GeV/c) MC;;counts");
        h_RecoTransMomPim->GetXaxis()->SetRangeUser(0.02,1.49);
        h_RecoTransMomPim->GetYaxis()->SetRangeUser(0,1.49);
        h_RecoTransMomPim->GetZaxis()->SetRangeUser(1,maxz);
        h_RecoTransMomPim->GetXaxis()->SetNdivisions(5);
        h_RecoTransMomPim->GetYaxis()->SetNdivisions(5);
        h_RecoTransMomPim->SetContour(99);
        h_RecoTransMomPim->Draw("COLZ");
        TBox* box6 = new TBox(0.13,1.15,0.5,1.45);
        box6->SetLineColor(kBlack);
        box6->SetFillColor(kWhite);
        box6->Draw("l same");
        TLatex* pimlabelaz = new TLatex(0.1, 0.8, "#pi^{#minus}");
        pimlabelaz->SetNDC();
        pimlabelaz->SetTextSize(40);
        pimlabelaz->SetTextFont(43);
        pimlabelaz->SetTextColor(kBlack);
        pimlabelaz->Draw("same");

        TString figure6name = figure_directory+"/benchmark_rho_recoquality.pdf";
        c6->Print(figure6name);

	double rhorecoeff = thiseff/100.0;
	rhorecoeff = rhorecoeff/5.0;//Test benchmark status failure
	setbenchstatus(rhorecoeff);

}
