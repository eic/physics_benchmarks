#include "RiceStyle.h"
using namespace std;
void plot_rho_physics_benchmark(TString filename="./benchmark_output/plot_combined.root"){
	Ssiz_t dotPosition = filename.Last('.');
	TString figure_directory = filename(0, dotPosition);
	figure_directory += "_figures";	
	
	TFile* file = new TFile(filename);
	TString vm_label="#rho^{0}";
	TString daug_label="#pi^{+}#pi^{-}";
	//t distribution
	TH1D* h_t_MC = (TH1D*) file->Get("h_t_MC");
	TH1D* h_t_REC = (TH1D*) file->Get("h_t_REC");
	TH1D* h_t_trk_REC = (TH1D*) file->Get("h_t_trk_REC");
	TH1D* h_t_combo_REC = (TH1D*) file->Get("h_t_combo_REC");
	//mass distribution
        TH1D* h_VM_mass_MC = (TH1D*) file->Get("h_VM_mass_MC");
        TH1D* h_VM_mass_REC = (TH1D*) file->Get("h_VM_mass_REC");
	TH1D* h_VM_mass_REC_justpions = (TH1D*) file->Get("h_VM_mass_REC_justpions");
	//mass distribution
        TH1D* h_dNdu_MC = (TH1D*) file->Get("h_u_MC");
        TH1D* h_dNdu_REC = (TH1D*) file->Get("h_u_REC");


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

	TLatex* r44_2 = new TLatex(0.63, 0.83, ""+vm_label+" #rightarrow "+daug_label);
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
        TH1D* base2 = makeHist("base2", "", "#pi^{+}#pi^{-} inv. mass (GeV)", "counts", 100,0.05,2.05,kBlack);
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
        w8->AddEntry(h_VM_mass_MC, "eSTARlight "+vm_label+" MC ", "L");
        w8->AddEntry(h_VM_mass_REC, vm_label+" RECO#times3", "L");
        w8->AddEntry(h_VM_mass_REC_justpions, vm_label+" RECO#times3 (#pi^{-}#pi^{+})", "L");
        w8->AddEntry(h_VM_mass_REC_justprotons, vm_label+" RECO#times3 (#pi^{-}p)", "L");
        w8->Draw("same");

	TString figure2name = figure_directory+"/benchmark_rho_mass.pdf";
        c2->Print(figure2name);


        TCanvas* c3 = new TCanvas("c3","c3",1,1,600,600);
        gPad->SetTicks();
	gPad->SetLogy(1);
        gPad->SetLeftMargin(0.18);
        gPad->SetBottomMargin(0.18);
        gPad->SetRightMargin(0.01);
        TH1D* base3 = makeHist("base3", "", "-#it{u} (GeV^{2})", "dN/d#it{u} (GeV^{-2} scaled)", 100,-0.25,3.05,kBlack);
        base3->GetYaxis()->SetRangeUser(0.5, 2.5*(h_dNdu_MC->GetMaximum()));
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

        h_dNdu_MC->SetFillColorAlpha(kBlack,0);
        h_dNdu_REC->SetFillColorAlpha(kMagenta,0);
        h_dNdu_MC->SetLineColor(kBlack);
        h_dNdu_REC->SetLineColor(kMagenta);
        h_dNdu_MC->SetLineWidth(2);
        h_dNdu_REC->SetLineWidth(2);

        h_dNdu_MC->Draw("HIST E same");
        h_dNdu_REC->Draw("HIST E same");

        r42->Draw("same");
        r43->Draw("same");
        r44->Draw("same");
        r44_2->Draw("same");

        TLegend *w9 = new TLegend(0.58,0.68,0.93,0.76);
        w9->SetLineColor(kWhite);
        w9->SetFillColor(0);
        w9->SetTextSize(17);
        w9->SetTextFont(45);
        w9->AddEntry(h_dNdu_MC, "eSTARlight "+vm_label+" MC ", "L");
        w9->AddEntry(h_dNdu_REC, "eSTARlight "+vm_label+" RECO ", "L");
        w9->Draw("same");

        TString figure3name = figure_directory+"/benchmark_rho_dNdu.pdf";
        c3->Print(figure3name);


}
