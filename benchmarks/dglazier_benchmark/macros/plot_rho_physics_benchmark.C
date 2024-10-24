#include "RiceStyle.h"

using namespace std;

void plot_rho_physics_benchmark(TString filename="./sim_output/plot_combined.root"){
  Ssiz_t dotPosition = filename.Last('.');
  TString figure_directory = filename(0, dotPosition);
  figure_directory += "_figures";	
  
  TFile* file = new TFile(filename);
  TString vm_label="#rho^{0}";
  TString daug_label="#pi^{#plus}#pi^{#minus}";
  //mass distribution
  TH1D* h_VM_mass_MC = (TH1D*) file->Get("h_VM_mass_MC");
  TH1D* h_VM_mass_REC = (TH1D*) file->Get("h_VM_mass_REC");
  TH1D* h_VM_mass_REC_justpions = (TH1D*) file->Get("h_VM_mass_REC_justpions");
  //mass distribution within B0
  TH1D* h_VM_mass_MC_etacut = (TH1D*) file->Get("h_VM_mass_MC_etacut");
  TH1D* h_VM_mass_REC_etacut = (TH1D*) file->Get("h_VM_mass_REC_etacut");
  //efficiencies
  TProfile2D* h_effEtaPtPi = (TProfile2D*) file->Get("h_effEtaPtPi");
  TProfile2D* h_effEtaPtPip = (TProfile2D*) file->Get("h_effEtaPtPip"); 
  TProfile2D* h_effEtaPtPim = (TProfile2D*) file->Get("h_effEtaPtPim"); 
  TProfile2D* h_effPhiEtaPi = (TProfile2D*) file->Get("h_effPhiEtaPi"); 
  TProfile2D* h_effPhiEtaPip = (TProfile2D*) file->Get("h_effPhiEtaPip"); 
  TProfile2D* h_effPhiEtaPim = (TProfile2D*) file->Get("h_effPhiEtaPim"); 	


  TLatex* r42 = new TLatex(0.18, 0.91, "ep 10#times100 GeV");
  r42->SetNDC();
  r42->SetTextSize(22);
  r42->SetTextFont(43);
  r42->SetTextColor(kBlack);

  TLatex* r43 = new TLatex(0.9,0.91, "#bf{EPIC}");
  r43->SetNDC();
  //r43->SetTextSize(0.04);
  r43->SetTextSize(22);

  TLatex* r44 = new TLatex(0.53, 0.78, "10^{-3}<Q^{2}<10 GeV^{2}, W>2 GeV");
  r44->SetNDC();
  r44->SetTextSize(20);
  r44->SetTextFont(43);
  r44->SetTextColor(kBlack);

  TLatex* r44_2 = new TLatex(0.5, 0.83, ""+vm_label+" #rightarrow "+daug_label+" eSTARlight");
  r44_2->SetNDC();
  r44_2->SetTextSize(30);
  r44_2->SetTextFont(43);
  r44_2->SetTextColor(kBlack);

  TCanvas* c2 = new TCanvas("c2","c2",1,1,600,600);
  gPad->SetTicks();
  gPad->SetLeftMargin(0.18);
  gPad->SetBottomMargin(0.18);
  gPad->SetTopMargin(0.10);
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

  TString figure1name = figure_directory+"/benchmark_rho_mass.pdf";
  c2->Print(figure1name);


  ///////////////////// Figure 4
  TCanvas* c4 = new TCanvas("c4","c4",1,1,600,600);
  gPad->SetTicks();
  gPad->SetLeftMargin(0.18);
  gPad->SetBottomMargin(0.18);
  gPad->SetTopMargin(0.10);
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

  TString figure2name = figure_directory+"/benchmark_rho_mass_cuts.pdf";
  c4->Print(figure2name);

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
  h_effEtaPtPim ->GetZaxis()->SetLabelSize(h_effEtaPtPim ->GetZaxis()->GetLabelSize()*0.5);
  h_effEtaPtPim ->GetZaxis()->SetTitleSize(h_effEtaPtPim ->GetZaxis()->GetTitleSize()*0.5);
  h_effPhiEtaPi ->GetXaxis()->SetLabelSize(h_effPhiEtaPi ->GetXaxis()->GetLabelSize()*1.8);
  h_effPhiEtaPip->GetXaxis()->SetLabelSize(h_effPhiEtaPip->GetXaxis()->GetLabelSize()*1.8);
  h_effPhiEtaPim->GetXaxis()->SetLabelSize(h_effPhiEtaPim->GetXaxis()->GetLabelSize()*1.8);
  h_effPhiEtaPi ->GetYaxis()->SetLabelSize(h_effPhiEtaPi ->GetYaxis()->GetLabelSize()*1.8);
  h_effPhiEtaPim->GetZaxis()->SetLabelSize(h_effPhiEtaPim->GetZaxis()->GetLabelSize()*0.5);
  h_effPhiEtaPim->GetZaxis()->SetTitleSize(h_effPhiEtaPim->GetZaxis()->GetTitleSize()*0.5);

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
  TLatex* r43fig5 = new TLatex(0.66,0.93, "#bf{EPIC}");
  r43fig5->SetNDC();
  r43fig5->SetTextSize(15);
  r43fig5->SetTextFont(43);
  r43fig5->SetTextColor(kBlack);
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

  TString figure3name = figure_directory+"/benchmark_rho_efficiencies.pdf";
  c5->Print(figure3name);
}
