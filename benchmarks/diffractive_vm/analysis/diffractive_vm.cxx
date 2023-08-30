#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include <TCanvas.h>
#include <TChain.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLorentzVector.h>
#include <TRandom3.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>

#include <TLorentzRotation.h>
#include <TVector2.h>
#include <TVector3.h>

#include "fmt/color.h"
#include "fmt/core.h"

#include "nlohmann/json.hpp"

#define PI 3.1415926
#define MASS_ELECTRON 0.00051
#define MASS_PROTON 0.93827
#define MASS_PION 0.13957
#define MASS_KAON 0.493667
#define MASS_AU197 183.45406466643374

auto giveme_t_method_L(TLorentzVector eIn, TLorentzVector eOut, TLorentzVector pIn,
                       TLorentzVector vmOut)
{
  TLorentzVector aInVec(pIn.Px() * 197, pIn.Py() * 197, pIn.Pz() * 197,
                        sqrt(pIn.Px() * 197 * pIn.Px() * 197 + pIn.Py() * 197 * pIn.Py() * 197 +
                             pIn.Pz() * 197 * pIn.Pz() * 197 + MASS_AU197 * MASS_AU197));
  double         method_L         = 0;
  TLorentzVector a_beam_scattered = aInVec - (vmOut + eOut - eIn);
  double         p_Aplus          = a_beam_scattered.E() + a_beam_scattered.Pz();
  double         p_TAsquared      = TMath::Power(a_beam_scattered.Pt(), 2);
  double         p_Aminus         = (MASS_AU197 * MASS_AU197 + p_TAsquared) / p_Aplus;
  TLorentzVector a_beam_scattered_corr;
  a_beam_scattered_corr.SetPxPyPzE(a_beam_scattered.Px(), a_beam_scattered.Py(),
                                   (p_Aplus - p_Aminus) / 2., (p_Aplus + p_Aminus) / 2.);
  method_L = -(a_beam_scattered_corr - aInVec).Mag2();

  return method_L;
}

TH1D* makeHist(const char* name, const char* title, const char* xtit, const char* ytit,
               const int nBins, const double lower, const double higher, EColor color = kBlack)
{

  TH1D* temp = new TH1D(name, title, nBins, lower, higher);

  temp->SetMarkerSize(1.0);
  temp->SetMarkerStyle(20);
  temp->SetMarkerColor(color);
  temp->SetLineColor(color);
  temp->SetStats(kFALSE);

  temp->GetXaxis()->SetTitle(xtit);
  temp->GetXaxis()->SetTitleSize(0.05);
  temp->GetXaxis()->SetTitleFont(42);
  temp->GetXaxis()->SetTitleOffset(1.25);
  temp->GetXaxis()->SetLabelSize(0.05);
  temp->GetXaxis()->SetLabelOffset(0.01);
  temp->GetXaxis()->SetLabelFont(42);
  temp->GetXaxis()->SetLabelColor(kBlack);
  temp->GetXaxis()->CenterTitle();

  temp->GetYaxis()->SetTitle(ytit);
  temp->GetYaxis()->SetTitleSize(0.05);
  temp->GetYaxis()->SetTitleFont(42);
  temp->GetYaxis()->SetTitleOffset(1.4);
  temp->GetYaxis()->SetLabelSize(0.05);
  temp->GetYaxis()->SetLabelOffset(0.01);
  temp->GetYaxis()->SetLabelFont(42);
  temp->GetYaxis()->SetLabelColor(kBlack);
  temp->GetYaxis()->CenterTitle();

  return temp;
}

void fixedFontHist1D(TH1* h, Float_t xoffset = 1.5, Float_t yoffset = 2.3)
{
  h->SetLabelFont(43, "X");
  h->SetLabelFont(43, "Y");
  // h->SetLabelOffset(0.01);
  h->SetLabelSize(16);
  h->SetTitleFont(43);
  h->SetTitleSize(20);
  h->SetLabelSize(15, "Y");
  h->SetTitleFont(43, "Y");
  h->SetTitleSize(20, "Y");
  h->SetTitleOffset(xoffset, "X");
  h->SetTitleOffset(yoffset, "Y");
  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->CenterTitle();
}

int diffractive_vm(const std::string& config_name)
{
  // read our configuration
  std::ifstream  config_file{config_name};
  nlohmann::json config;
  config_file >> config;

  const std::string rec_file      = config["rec_file"];
  const std::string vm_name       = config["vm_name"];
  const std::string detector      = config["detector"];
  const std::string output_prefix = config["output_prefix"];
  const std::string test_tag      = config["test_tag"];
  const int         ebeam         = config["ebeam"];
  const int         pbeam         = config["pbeam"];

  fmt::print(fmt::emphasis::bold | fg(fmt::color::forest_green),
             "Running DIS electron analysis...\n");
  fmt::print(" - Detector package: {}\n", detector);
  fmt::print(" - input file: {}\n", rec_file);
  fmt::print(" - output prefix: {}\n", output_prefix);
  fmt::print(" - test tag: {}\n", test_tag);
  fmt::print(" - ebeam: {}\n", ebeam);
  fmt::print(" - pbeam: {}\n", pbeam);

  auto tree = new TChain("events");
  tree->Add(rec_file.c_str());
  TTreeReader tree_reader(tree); // !the tree reader

  TTreeReaderArray<int> mc_genStatus_array = {tree_reader, "MCParticles.generatorStatus"};
  // MC particle pz array for each MC particle
  TTreeReaderArray<float>  mc_px_array   = {tree_reader, "MCParticles.momentum.x"};
  TTreeReaderArray<float>  mc_py_array   = {tree_reader, "MCParticles.momentum.y"};
  TTreeReaderArray<float>  mc_pz_array   = {tree_reader, "MCParticles.momentum.z"};
  TTreeReaderArray<double> mc_mass_array = {tree_reader, "MCParticles.mass"};
  TTreeReaderArray<int>    mc_pdg_array  = {tree_reader, "MCParticles.PDG"};

  // Reconstructed EcalEndcapNClusters
  TTreeReaderArray<float> em_energy_array     = {tree_reader, "EcalEndcapNClusters.energy"};
  TTreeReaderArray<float> em_x_array          = {tree_reader, "EcalEndcapNClusters.position.x"};
  TTreeReaderArray<float> em_y_array          = {tree_reader, "EcalEndcapNClusters.position.y"};
  TTreeReaderArray<float> emhits_x_array      = {tree_reader, "EcalEndcapNRecHits.position.x"};
  TTreeReaderArray<float> emhits_y_array      = {tree_reader, "EcalEndcapNRecHits.position.y"};
  TTreeReaderArray<float> emhits_energy_array = {tree_reader, "EcalEndcapNRecHits.energy"};

  TTreeReaderArray<unsigned int> em_rec_id_array = {tree_reader,
                                                    "EcalEndcapNClusterAssociations.recID"};
  TTreeReaderArray<unsigned int> em_sim_id_array = {tree_reader,
                                                    "EcalEndcapNClusterAssociations.simID"};

  // Reconstructed particles pz array for each reconstructed particle
  TTreeReaderArray<float> reco_px_array = {tree_reader, "ReconstructedChargedParticles.momentum.x"};
  TTreeReaderArray<float> reco_py_array = {tree_reader, "ReconstructedChargedParticles.momentum.y"};
  TTreeReaderArray<float> reco_pz_array = {tree_reader, "ReconstructedChargedParticles.momentum.z"};
  TTreeReaderArray<float> reco_charge_array = {tree_reader, "ReconstructedChargedParticles.charge"};

  TTreeReaderArray<unsigned int> rec_id = {tree_reader,
                                           "ReconstructedChargedParticleAssociations.recID"};
  TTreeReaderArray<unsigned int> sim_id = {tree_reader,
                                           "ReconstructedChargedParticleAssociations.simID"};

  // events
  TH1D* h_Q2_e      = new TH1D("h_Q2_e", ";Q^{2}_{e,MC}", 100, 0, 20);
  TH1D* h_y_e       = new TH1D("h_y_e", ";y_{e,MC}", 100, 0, 1);
  TH1D* h_energy_MC = new TH1D("h_energy_MC", ";E_{MC} (GeV)", 100, 0, 20);
  TH1D* h_t_MC      = new TH1D("h_t_MC", ";t_{MC}; counts", 100, 0, 0.2);

  TH1D* h_Q2REC_e        = new TH1D("h_Q2REC_e", ";Q^{2}_{e,REC}", 100, 0, 20);
  TH1D* h_yREC_e         = new TH1D("h_yREC_e", ";y_{e,REC}", 100, 0, 1);
  TH1D* h_energy_REC     = new TH1D("h_energy_REC", ";E_{REC} (GeV)", 100, 0, 20);
  TH1D* h_trk_energy_REC = new TH1D("h_trk_energy_REC", ";E_{REC} (GeV)", 100, 0, 20);
  TH1D* h_trk_Epz_REC    = new TH1D("h_trk_Epz_REC", ";E - p_{z} (GeV)", 200, 0, 50);

  // track
  TH1D* h_eta = new TH1D("h_eta", ";#eta", 100, -5, 5);
  TH2D* h_trk_energy_res =
      new TH2D("h_trk_energy_res", ";E_{MC} (GeV); E_{MC}-E_{REC}/E_{MC} track-base ", 100, 0, 20,
               1000, -1, 1);
  TH2D* h_trk_Pt_res =
      new TH2D("h_trk_Pt_res", ";p_{T,MC} (GeV); P_{T,MC}-P_{T,REC}/P_{T,MC} track-base ", 100, 0,
               15, 1000, -1, 1);
  TH1D* h_Epz_REC = new TH1D("h_Epz_REC", ";E - p_{z} (GeV)", 200, 0, 50);

  // VM & t
  TH1D* h_VM_mass_REC = new TH1D("h_VM_mass_REC", ";mass (GeV)", 200, 0, 4);
  TH1D* h_VM_pt_REC   = new TH1D("h_VM_pt_REC", ";p_{T} (GeV/c)", 200, 0, 2);
  TH2D* h_VM_res =
      new TH2D("h_VM_res", ";p_{T,MC} (GeV); p_{T,MC}-E_{T,REC}/p_{T,MC}", 100, 0, 2, 1000, -1, 1);
  TH1D* h_t_REC     = new TH1D("h_t_REC", ";t_{REC} (GeV^{2}); counts", 100, 0, 0.2);
  TH1D* h_t_trk_REC = new TH1D("h_t_trk_REC", ";t_{REC}(GeV^{2}) track-base; counts", 100, 0, 0.2);
  TH1D* h_t_combo_REC = new TH1D("h_t_combo_REC", ";t_{combo,REC}(GeV^{2}); counts", 100, 0, 0.2);
  TH2D* h_t_res =
      new TH2D("h_t_res", ";t_{MC} (GeV^{2}); t_{MC}-t_{REC}/t_{MC}", 100, 0, 0.2, 1000, -10, 10);
  TH2D* h_trk_t_res = new TH2D("h_trk_t_res", ";t_{MC} (GeV^{2}); t_{MC}-t_{REC}/t_{MC} track-base",
                               100, 0, 0.2, 1000, -10, 10);
  TH2D* h_t_2D = new TH2D("h_t_2D", ";t_{MC} (GeV^{2}); t_{REC} (GeV^{2}) track-base", 100, 0, 0.2,
                          100, 0, 0.2);
  TH2D* h_t_REC_2D   = new TH2D("h_t_REC_2D", ";t_{trk,REC} (GeV^{2}); t_{EEMC,REC} (GeV^{2})", 100,
                                0, 0.2, 100, 0, 0.2);
  TH2D* h_t_RECMC_2D = new TH2D("h_t_RECMC_2D", ";t_{MC} (GeV^{2}); t_{trk,REC} / t_{EEMC,REC} ",
                                100, 0, 0.2, 200, -10, 10);

  // energy clus
  TH2D* h_emClus_position_REC =
      new TH2D("h_emClus_position_REC", ";x (mm);y (mm)", 80, -800, 800, 80, -800, 800);
  TH2D* h_emHits_position_REC =
      new TH2D("h_emHits_position_REC", ";x (mm);y (mm)", 80, -800, 800, 80, -800, 800);
  TH2D* h_energy_res = new TH2D("h_energy_res", ";E_{MC} (GeV); E_{MC}-E_{REC}/E_{MC} emcal", 100,
                                0, 20, 1000, -1, 1);
  TH1D* h_energy_calibration_REC = new TH1D("h_energy_calibration_REC", ";E (GeV)", 200, 0, 2);
  TH1D* h_EoverP_REC             = new TH1D("h_EoverP_REC", ";E/p", 200, 0, 2);
  TH1D* h_ClusOverHit_REC =
      new TH1D("h_ClusOverHit_REC", ";cluster energy / new cluster energy", 200, 0, 2);

  tree_reader.SetEntriesRange(0, tree->GetEntries());
  while (tree_reader.Next()) {

    /*
    Beam particles
    */
    TLorentzVector ebeam(0, 0, 0, 0);
    TLorentzVector pbeam(0, 0, 0, 0);

    TLorentzVector vmMC(0, 0, 0, 0);
    TLorentzVector kplusMC(0, 0, 0, 0);
    TLorentzVector kminusMC(0, 0, 0, 0);

    // MC level
    TLorentzVector scatMC(0, 0, 0, 0);
    int            mc_elect_index = -1;
    double         maxPt          = -99.;
    for (int imc = 0; imc < mc_px_array.GetSize(); imc++) {
      TVector3 mctrk(mc_px_array[imc], mc_py_array[imc], mc_pz_array[imc]);
      if (mc_genStatus_array[imc] == 4) { // 4 is Sartre.
        if (mc_pdg_array[imc] == 11)
          ebeam.SetVectM(mctrk, MASS_ELECTRON);
        if (mc_pdg_array[imc] == 2212)
          pbeam.SetVectM(mctrk, MASS_PROTON);
      }
      if (mc_genStatus_array[imc] != 1)
        continue;
      if (mc_pdg_array[imc] == 11 && mctrk.Perp() > maxPt) {
        maxPt          = mctrk.Perp();
        mc_elect_index = imc;
        scatMC.SetVectM(mctrk, mc_mass_array[imc]);
      }
      if (mc_pdg_array[imc] == 321 && mc_genStatus_array[imc] == 1)
        kplusMC.SetVectM(mctrk, MASS_KAON);
      if (mc_pdg_array[imc] == -321 && mc_genStatus_array[imc] == 1)
        kminusMC.SetVectM(mctrk, MASS_KAON);
    }
    vmMC = kplusMC + kminusMC;
    // protection.
    if (ebeam.E() == pbeam.E() && ebeam.E() == 0) {
      std::cout << "problem with MC incoming beams" << std::endl;
      continue;
    }
    TLorentzVector qbeam = ebeam - scatMC;
    double         Q2    = -(qbeam).Mag2();
    double         pq    = pbeam.Dot(qbeam);
    double         y     = pq / pbeam.Dot(ebeam);

    // MC level phase space cut
    if (Q2 < 1. || Q2 > 10.)
      continue;
    if (y < 0.01 || y > 0.85)
      continue;

    h_Q2_e->Fill(Q2);
    h_y_e->Fill(y);
    h_energy_MC->Fill(scatMC.E());

    double t_MC = 0.;
    if (vmMC.E() != 0 && fabs(vmMC.Rapidity()) < 3.5) {
      double method_E = -(qbeam - vmMC).Mag2();
      t_MC            = method_E;
      h_t_MC->Fill(method_E);
    }

    // rec level
    // leading cluster
    double maxEnergy = -99.;
    double xpos      = -999.;
    double ypos      = -999.;
    for (int iclus = 0; iclus < em_energy_array.GetSize(); iclus++) {
      if (em_energy_array[iclus] > maxEnergy) {
        maxEnergy = em_energy_array[iclus];
        xpos      = em_x_array[iclus];
        ypos      = em_y_array[iclus];
      }
    }
    // leading hit energy
    double maxHitEnergy = 0.01; // threshold 10 MeV
    double xhitpos      = -999.;
    double yhitpos      = -999.;
    int    hit_index    = -1;
    for (int ihit = 0; ihit < emhits_energy_array.GetSize(); ihit++) {
      if (emhits_energy_array[ihit] > maxHitEnergy) {
        maxHitEnergy = emhits_energy_array[ihit];
        xhitpos      = emhits_x_array[ihit];
        yhitpos      = emhits_y_array[ihit];
        hit_index    = ihit;
      }
    }
    // sum over all 3x3 towers around the leading tower
    double xClus = xhitpos * maxHitEnergy;
    double yClus = yhitpos * maxHitEnergy;
    for (int ihit = 0; ihit < emhits_energy_array.GetSize(); ihit++) {
      double hitenergy = emhits_energy_array[ihit];
      double x         = emhits_x_array[ihit];
      double y         = emhits_y_array[ihit];
      double d         = sqrt((x - xhitpos) * (x - xhitpos) + (y - yhitpos) * (y - yhitpos));
      if (d < 70. && ihit != hit_index && hitenergy > 0.01) {
        maxHitEnergy += hitenergy; // clustering around leading tower 3 crystal = 60mm.
        xClus += x * hitenergy;
        yClus += y * hitenergy;
      }
    }

    h_ClusOverHit_REC->Fill(maxEnergy / maxHitEnergy);
    // weighted average cluster position.
    xClus         = xClus / maxHitEnergy;
    yClus         = yClus / maxHitEnergy;
    double radius = sqrt(xClus * xClus + yClus * yClus);
    if (radius > 550.)
      continue; // geometric acceptance cut
    // 4.4% energy calibration.
    double clusEnergy = 1.044 * maxHitEnergy;

    h_energy_REC->Fill(clusEnergy);
    // ratio of reco / truth Energy
    h_energy_calibration_REC->Fill(clusEnergy / scatMC.E());
    // energy resolution
    double res = (scatMC.E() - clusEnergy) / scatMC.E();
    h_energy_res->Fill(scatMC.E(), res);
    h_emClus_position_REC->Fill(xpos, ypos);   // default clustering position
    h_emHits_position_REC->Fill(xClus, yClus); // self clustering position

    // association of rec level scat' e
    int rec_elect_index = -1;
    for (int i = 0; i < sim_id.GetSize(); i++) {
      if (sim_id[i] == mc_elect_index) {
        // find the rec_id
        rec_elect_index = rec_id[i];
      }
    }

    TLorentzVector scatMCmatchREC(0, 0, 0, 0);
    TLorentzVector scatREC(0, 0, 0, 0);
    TLorentzVector scatClusEREC(0, 0, 0, 0);
    TLorentzVector hfs(0, 0, 0, 0);
    TLorentzVector particle(0, 0, 0, 0);
    TLorentzVector kplusREC(0, 0, 0, 0);
    TLorentzVector kminusREC(0, 0, 0, 0);
    TLorentzVector vmREC(0, 0, 0, 0);

    double maxP = -1.;
    // track loop
    for (int itrk = 0; itrk < reco_pz_array.GetSize(); itrk++) {
      TVector3 trk(reco_px_array[itrk], reco_py_array[itrk], reco_pz_array[itrk]);
      if (rec_elect_index != -1 && itrk == rec_elect_index) {
        scatMCmatchREC.SetVectM(trk, MASS_ELECTRON); // Reserved to calculate t.
      }
      if (trk.Mag() > maxP) {
        // track-base 4 vector
        maxP = trk.Mag();
        scatREC.SetVectM(trk, MASS_ELECTRON);

        // use emcal energy to define 4 vector
        double p   = sqrt(clusEnergy * clusEnergy - MASS_ELECTRON * MASS_ELECTRON);
        double eta = scatREC.Eta();
        double phi = scatREC.Phi();
        double pt  = TMath::Sin(scatREC.Theta()) * p;
        scatClusEREC.SetPtEtaPhiM(pt, eta, phi, MASS_ELECTRON);
      }
    }
    // loop over track again;
    for (int itrk = 0; itrk < reco_pz_array.GetSize(); itrk++) {
      TVector3 trk(reco_px_array[itrk], reco_py_array[itrk], reco_pz_array[itrk]);
      particle.SetVectM(trk, MASS_PION); // assume pions;
      if (itrk != rec_elect_index) {
        hfs += particle; // hfs 4vector sum.
        // selecting phi->kk daughters;
        h_eta->Fill(trk.Eta());
        if (fabs(trk.Eta()) < 3.0) {
          if (reco_charge_array[itrk] > 0)
            kplusREC.SetVectM(trk, MASS_KAON);
          if (reco_charge_array[itrk] < 0)
            kminusREC.SetVectM(trk, MASS_KAON);
        }
      }
    }
    // 4vector of VM;
    if (kplusREC.E() != 0. && kminusREC.E() != 0.) {
      vmREC = kplusREC + kminusREC;
    }

    // track-base e' energy REC;
    h_trk_energy_REC->Fill(scatMCmatchREC.E());

    // track-base e' energy resolution;
    res = (scatMC.E() - scatMCmatchREC.E()) / scatMC.E();
    h_trk_energy_res->Fill(scatMC.E(), res);

    // track-base e' pt resolution;
    res = (scatMC.Pt() - scatMCmatchREC.Pt()) / scatMC.Pt();
    h_trk_Pt_res->Fill(scatMC.Pt(), res);

    // track-base Epz scat' e
    double EpzREC = (scatMCmatchREC + hfs).E() - (scatMCmatchREC + hfs).Pz();
    h_trk_Epz_REC->Fill(EpzREC);

    // EEMC cluster Epz scat' e
    EpzREC = (scatClusEREC + hfs).E() - (scatClusEREC + hfs).Pz();
    h_Epz_REC->Fill(EpzREC);

    // E over p
    double EoverP = scatClusEREC.E() / scatMCmatchREC.P();
    h_EoverP_REC->Fill(EoverP);

    // cluster-base DIS kine;
    TLorentzVector qbeamREC = ebeam - scatClusEREC;
    double         Q2REC    = -(qbeamREC).Mag2();
    double         pqREC    = pbeam.Dot(qbeamREC);
    double         yREC     = pqREC / pbeam.Dot(ebeam);
    h_Q2REC_e->Fill(Q2REC);
    h_yREC_e->Fill(yREC);

    // Event selection:
    if (EpzREC < 27 || EpzREC > 40)
      continue;
    if (EoverP < 0.8 || EoverP > 1.18)
      continue;

    // MC level phase space cut
    if (Q2REC < 1. || Q2REC > 10.)
      continue;
    if (yREC < 0.01 || yREC > 0.85)
      continue;

    // VM rec
    if (vmREC.E() == 0)
      continue;
    double phi_mass = vmREC.M();
    h_VM_mass_REC->Fill(phi_mass);
    h_VM_pt_REC->Fill(vmREC.Pt());

    // select phi mass and rapidity window
    if (fabs(phi_mass - 1.02) < 0.02 && fabs(vmREC.Rapidity()) < 3.5) {
      // 2 versions: track and energy cluster:
      double t_trk_REC = giveme_t_method_L(ebeam, scatMCmatchREC, pbeam, vmREC);
      double t_REC     = giveme_t_method_L(ebeam, scatClusEREC, pbeam, vmREC);
      h_t_trk_REC->Fill(t_trk_REC);
      h_t_REC->Fill(t_REC);
      h_t_REC_2D->Fill(t_trk_REC, t_REC);
      if ((t_trk_REC / t_REC) > 0.5 && (t_trk_REC / t_REC) < 1.5) {
        h_t_combo_REC->Fill((t_trk_REC + t_REC) / 2.); // w=1./(fabs(1.0-(t_trk_REC/t_REC)))
      }
      h_t_RECMC_2D->Fill(t_MC, t_trk_REC / t_REC);

      // t track resolution
      res = (t_MC - t_trk_REC) / t_MC;
      h_trk_t_res->Fill(t_MC, res);

      // t EEMC resolution;
      res = (t_MC - t_REC) / t_MC;
      h_t_res->Fill(t_MC, res);

      // 2D t
      h_t_2D->Fill(t_MC, t_trk_REC);

      // VM pt resolution;
      res = (vmMC.Pt() - vmREC.Pt()) / vmMC.Pt();
      h_VM_res->Fill(vmMC.Pt(), res);
    }
  }

  TString vm_label;
  TString daug_label;
  if (vm_name == "phi") {
    vm_label   = "#phi";
    daug_label = "K^{+}K^{-}";
  } else if (vm_name == "jpsi") {
    vm_label   = "J/#psi";
    daug_label = "e^{+}e^{-}";
  } else {
    throw std::runtime_error(fmt::format("Unknown vm_name = \"{}\"", vm_name));
  }

  {

    TCanvas c("c1", "c1", 1, 1, 600, 600);
    gPad->SetLogy(1);
    gPad->SetTicks();
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.01);
    TH1D* base1 = makeHist("base1", "", "|#it{t} | (GeV^{2})", "dN/d|#it{t} | (GeV^{-2}) ", 100, 0,
                           0.18, kBlack);
    base1->GetYaxis()->SetRangeUser(8e-2, 8e5);
    base1->GetXaxis()->SetTitleColor(kBlack);
    fixedFontHist1D(base1, 1., 1.2);
    base1->GetYaxis()->SetTitleSize(base1->GetYaxis()->GetTitleSize() * 1.5);
    base1->GetXaxis()->SetTitleSize(base1->GetXaxis()->GetTitleSize() * 1.5);
    base1->GetYaxis()->SetLabelSize(base1->GetYaxis()->GetLabelSize() * 1.5);
    base1->GetXaxis()->SetLabelSize(base1->GetXaxis()->GetLabelSize() * 1.5);
    base1->GetXaxis()->SetNdivisions(4, 4, 0);
    base1->GetYaxis()->SetNdivisions(5, 5, 0);
    base1->Draw();

    h_t_MC->Draw("same");

    h_t_REC->SetMarkerStyle(20);
    h_t_REC->Draw("PEsame");

    h_t_trk_REC->SetFillColorAlpha(kBlue, 0.4);
    h_t_trk_REC->SetFillStyle(1001);
    h_t_trk_REC->SetMarkerStyle(24);
    h_t_trk_REC->SetMarkerColor(kBlue);
    // h_t_trk_REC->Draw("PE3same");

    h_t_combo_REC->SetFillColorAlpha(kRed, 0.4);
    h_t_combo_REC->SetFillStyle(1001);
    h_t_combo_REC->SetMarkerStyle(24);
    h_t_combo_REC->SetMarkerColor(kRed);
    // h_t_combo_REC->Draw("PE3same");

    TLatex* r42 = new TLatex(0.18, 0.91, "eAu 18x110 GeV");
    r42->SetNDC();
    r42->SetTextSize(22);
    r42->SetTextFont(43);
    r42->SetTextColor(kBlack);
    r42->Draw("same");

    TLatex* r43 = new TLatex(0.9, 0.91, "EPIC");
    r43->SetNDC();
    r43->SetTextSize(0.04);
    r43->Draw("same");

    TLatex* r44 = new TLatex(0.18, 0.84, "1<Q^{2}<10 GeV^{2}, 0.01 < y < 0.95");
    r44->SetNDC();
    r44->SetTextSize(20);
    r44->SetTextFont(43);
    r44->SetTextColor(kBlack);
    r44->Draw("same");

    TLatex* r44_0 = new TLatex(0.18, 0.79,
                               "  |y_{" + vm_label + "}|<3.5, |M_{inv} #minus M_{" + vm_label +
                                   "}| < 0.02 GeV");
    r44_0->SetNDC();
    r44_0->SetTextSize(20);
    r44_0->SetTextFont(43);
    r44_0->SetTextColor(kBlack);
    r44_0->Draw("same");

    TLatex* r44_2 = new TLatex(0.18, 0.18, "" + vm_label + " #rightarrow " + daug_label);
    r44_2->SetNDC();
    r44_2->SetTextSize(30);
    r44_2->SetTextFont(43);
    r44_2->SetTextColor(kBlack);
    r44_2->Draw("same");

    TLegend* w7 = new TLegend(0.48, 0.68, 0.93, 0.76);
    w7->SetLineColor(kWhite);
    w7->SetFillColor(0);
    w7->SetTextSize(17);
    w7->SetTextFont(45);
    // if(filename=="MCclusterEnergy")
    //{
    //	w7->AddEntry(h_t_MC, "Sartre "+vm_label+" MC ", "L");
    //	w7->AddEntry(h_t_REC, "Sartre "+vm_label+" RECO w. true EEMC E ", "P");
    //	w7->AddEntry(h_t_trk_REC, "Sartre "+vm_label+" RECO track only ", "P");

    //}
    // else if(filename=="MCvmAndelectron")
    //{
    //	w7->AddEntry(h_t_MC, "Sartre "+vm_label+" MC ", "L");
    //	w7->AddEntry(h_t_REC, "Sartre "+vm_label+" RECO w. true e' ", "P");
    //	w7->AddEntry(h_t_trk_REC, "Sartre "+vm_label+" MC w. RECO e' ", "P");
    //}
    // else{
    w7->AddEntry(h_t_MC, "Sartre " + vm_label + " MC ", "L");
    w7->AddEntry(h_t_REC, "Sartre " + vm_label + " RECO w. EEMC ", "P");
    // w7->AddEntry(h_t_trk_REC, "Sartre "+vm_label+" RECO track only", "P");
    // w7->AddEntry(h_t_combo_REC, "Sartre "+vm_label+" RECO best", "P");

    //}

    w7->Draw("same");

    c.Print(fmt::format("{}_benchmark-{}-dsigmadt.pdf", output_prefix, vm_name).c_str());
  }

  {
    TCanvas c("c2", "c2", 1, 1, 800, 800);
    gPad->SetTicks();
    gPad->SetLogy(1);
    gPad->SetLeftMargin(0.13);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.1);
    gPad->SetBottomMargin(0.13);

    TH1D* base1 = makeHist("base1", "", "|#it{t} | (GeV^{2})", " #delta t/t (resolution) ", 100, 0,
                           0.2, kBlack);
    base1->GetYaxis()->SetRangeUser(1e-2, 1000);
    base1->GetXaxis()->SetTitleColor(kBlack);
    fixedFontHist1D(base1, 1.2, 1.6);
    base1->GetYaxis()->SetTitleSize(base1->GetYaxis()->GetTitleSize() * 1.5);
    base1->GetXaxis()->SetTitleSize(base1->GetXaxis()->GetTitleSize() * 1.5);
    base1->GetYaxis()->SetLabelSize(base1->GetYaxis()->GetLabelSize() * 1.8);
    base1->GetXaxis()->SetLabelSize(base1->GetXaxis()->GetLabelSize() * 1.7);
    base1->GetXaxis()->SetNdivisions(5, 5, 0);
    base1->GetYaxis()->SetNdivisions(5, 5, 0);
    base1->Draw();
    TH2D* h_res;
    h_res = (TH2D*)h_t_res;
    // h_res=(TH2D*) h_trk_t_res;

    TH1D* h_res_1D = new TH1D("h_res_1D", "", 100, 0, 0.2);
    for (int ibin = 0; ibin < h_res->GetNbinsX(); ibin++) {
      TH1D*  tmp        = h_res->ProjectionY("tmp", ibin + 1, ibin + 1);
      double sigma      = tmp->GetStdDev();
      double sigmaerror = tmp->GetStdDevError();
      h_res_1D->SetBinContent(ibin + 1, sigma);
      h_res_1D->SetBinError(ibin + 1, sigmaerror);
    }

    h_res_1D->SetMarkerSize(1.6);
    h_res_1D->SetMarkerColor(kBlack);
    h_res_1D->SetLineColor(kBlack);
    h_res_1D->SetMarkerStyle(20);

    h_res_1D->Fit("pol0", "RMS0", "", 0.015, 0.025); // first  dip
    h_res_1D->Fit("pol0", "RMS0", "", 0.050, 0.060); // second dip
    h_res_1D->Fit("pol0", "RMS0", "", 0.100, 0.150); // third  dip
    h_res_1D->Draw("Psame");

    TLatex* r42 = new TLatex(0.15, 0.91, "eAu 18x110 GeV");
    r42->SetNDC();
    r42->SetTextSize(25);
    r42->SetTextFont(43);
    r42->SetTextColor(kBlack);
    r42->Draw("same");

    TLatex* r43 = new TLatex(0.9, 0.91, "EPIC");
    r43->SetNDC();
    r43->SetTextSize(0.04);
    r43->Draw("same");

    TLatex* r44_2 = new TLatex(0.15, 0.18, vm_label + " #rightarrow " + daug_label);
    r44_2->SetNDC();
    r44_2->SetTextSize(30);
    r44_2->SetTextFont(43);
    r44_2->SetTextColor(kBlack);
    r44_2->Draw("same");

    TPad* drawPad = new TPad("pad_etalab_11", "pad_etalab_11", 0.16, 0.53, 0.47, 0.83);
    drawPad->SetLeftMargin(0.18);
    drawPad->SetRightMargin(0);
    drawPad->SetTopMargin(0.0);
    drawPad->SetBottomMargin(0.18);
    drawPad->Draw("same");
    drawPad->SetTicks();
    drawPad->SetLogz(1);
    drawPad->cd();
    TH1D* base2 = makeHist("base2", "", "MC", " resolution ", 100, 0, 0.2, kBlack);
    base2->GetYaxis()->SetRangeUser(-10, 2);
    base2->GetXaxis()->SetTitleColor(kBlack);
    fixedFontHist1D(base2, 3, 3);
    base2->GetYaxis()->SetTitleSize(base2->GetYaxis()->GetTitleSize() * 1);
    base2->GetXaxis()->SetTitleSize(base2->GetXaxis()->GetTitleSize() * 1);
    base2->GetYaxis()->SetLabelSize(base2->GetYaxis()->GetLabelSize() * 1);
    base2->GetXaxis()->SetLabelSize(base2->GetXaxis()->GetLabelSize() * 1);
    base2->GetXaxis()->SetNdivisions(5, 5, 0);
    base2->GetYaxis()->SetNdivisions(5, 5, 0);
    base2->Draw();

    h_res->Draw("colzsame");

    c.Print(fmt::format("{}_benchmark-{}-t-resolution.pdf", output_prefix, vm_name).c_str());
  }

  {
    TCanvas c("c1", "c1", 1, 1, 1600, 800);
    c.Divide(4, 2, 0.01, 0.01);
    c.cd(1);
    gPad->SetLogy(1);
    gPad->SetLeftMargin(0.13);
    gPad->SetBottomMargin(0.15);
    h_Q2_e->GetXaxis()->SetTitleSize(0.8 * h_Q2_e->GetXaxis()->GetTitleSize());
    h_Q2_e->GetXaxis()->SetLabelSize(0.8 * h_Q2_e->GetXaxis()->GetLabelSize());
    h_Q2_e->GetYaxis()->SetTitleSize(0.8 * h_Q2_e->GetYaxis()->GetTitleSize());
    h_Q2_e->GetYaxis()->SetLabelSize(0.8 * h_Q2_e->GetYaxis()->GetLabelSize());
    h_Q2_e->GetXaxis()->SetTitleOffset(1.6 * h_Q2_e->GetXaxis()->GetTitleOffset());
    h_Q2_e->GetYaxis()->SetTitleOffset(2.0 * h_Q2_e->GetYaxis()->GetTitleOffset());
    h_Q2_e->GetYaxis()->SetTitle("counts");
    h_Q2_e->Draw();
    h_Q2REC_e->SetMarkerStyle(24);
    h_Q2REC_e->Draw("PEsame");
    TLegend* w7 = new TLegend(0.28, 0.7, 0.53, 0.86);
    w7->SetLineColor(kWhite);
    w7->SetFillColor(0);
    w7->SetTextSize(17);
    w7->SetTextFont(45);
    w7->AddEntry(h_energy_MC, "MC ", "L");
    w7->AddEntry(h_energy_REC, "RECO", "P");
    w7->Draw("same");

    c.cd(2);
    gPad->SetLogy(1);
    gPad->SetLeftMargin(0.13);
    gPad->SetBottomMargin(0.15);
    h_y_e->GetXaxis()->SetTitleSize(0.8 * h_y_e->GetXaxis()->GetTitleSize());
    h_y_e->GetXaxis()->SetLabelSize(0.8 * h_y_e->GetXaxis()->GetLabelSize());
    h_y_e->GetYaxis()->SetTitleSize(0.8 * h_y_e->GetYaxis()->GetTitleSize());
    h_y_e->GetYaxis()->SetLabelSize(0.8 * h_y_e->GetYaxis()->GetLabelSize());
    h_y_e->GetXaxis()->SetTitleOffset(1.6 * h_y_e->GetXaxis()->GetTitleOffset());
    h_y_e->GetYaxis()->SetTitleOffset(2.0 * h_y_e->GetYaxis()->GetTitleOffset());
    h_y_e->GetYaxis()->SetTitle("counts");
    h_y_e->Draw();
    h_yREC_e->SetMarkerStyle(24);
    h_yREC_e->Draw("PEsame");
    w7->Draw("same");

    c.cd(3);
    gPad->SetLogy(1);
    gPad->SetLeftMargin(0.13);
    gPad->SetBottomMargin(0.15);
    h_energy_MC->GetXaxis()->SetTitleSize(0.8 * h_energy_MC->GetXaxis()->GetTitleSize());
    h_energy_MC->GetXaxis()->SetLabelSize(0.8 * h_energy_MC->GetXaxis()->GetLabelSize());
    h_energy_MC->GetYaxis()->SetTitleSize(0.8 * h_energy_MC->GetYaxis()->GetTitleSize());
    h_energy_MC->GetYaxis()->SetLabelSize(0.8 * h_energy_MC->GetYaxis()->GetLabelSize());
    h_energy_MC->GetXaxis()->SetTitleOffset(1.6 * h_energy_MC->GetXaxis()->GetTitleOffset());
    h_energy_MC->GetYaxis()->SetTitleOffset(2.5 * h_energy_MC->GetYaxis()->GetTitleOffset());
    h_energy_MC->GetYaxis()->SetTitle("counts");
    h_energy_MC->Draw();
    h_energy_REC->SetMarkerStyle(24);
    h_energy_REC->Draw("PEsame");
    w7->Draw("same");
    // TLegend *w7 = new TLegend(0.28,0.7,0.53,0.86);
    // w7->SetLineColor(kWhite);
    // w7->SetFillColor(0);
    // w7->SetTextSize(17);
    // w7->SetTextFont(45);
    // w7->AddEntry(h_energy_MC, "MC ", "L");
    // w7->AddEntry(h_energy_REC, "RECO", "P");
    // w7->Draw("same");

    c.cd(4);
    gPad->SetLogy(1);
    gPad->SetLeftMargin(0.13);
    gPad->SetBottomMargin(0.15);
    h_Epz_REC->GetXaxis()->SetTitleSize(0.8 * h_Epz_REC->GetXaxis()->GetTitleSize());
    h_Epz_REC->GetXaxis()->SetLabelSize(0.8 * h_Epz_REC->GetXaxis()->GetLabelSize());
    h_Epz_REC->GetYaxis()->SetTitleSize(0.8 * h_Epz_REC->GetYaxis()->GetTitleSize());
    h_Epz_REC->GetYaxis()->SetLabelSize(0.8 * h_Epz_REC->GetYaxis()->GetLabelSize());
    h_Epz_REC->GetXaxis()->SetTitleOffset(1.6 * h_Epz_REC->GetXaxis()->GetTitleOffset());
    h_Epz_REC->GetYaxis()->SetTitleOffset(2.5 * h_Epz_REC->GetYaxis()->GetTitleOffset());
    h_Epz_REC->GetYaxis()->SetTitle("counts");
    h_Epz_REC->SetMarkerStyle(24);
    h_Epz_REC->Draw("PEsame");

    c.cd(5);
    gPad->SetLogy(1);
    gPad->SetLeftMargin(0.13);
    gPad->SetBottomMargin(0.15);
    h_energy_calibration_REC->GetXaxis()->SetTitleSize(
        0.8 * h_energy_calibration_REC->GetXaxis()->GetTitleSize());
    h_energy_calibration_REC->GetXaxis()->SetLabelSize(
        0.8 * h_energy_calibration_REC->GetXaxis()->GetLabelSize());
    h_energy_calibration_REC->GetYaxis()->SetTitleSize(
        0.8 * h_energy_calibration_REC->GetYaxis()->GetTitleSize());
    h_energy_calibration_REC->GetYaxis()->SetLabelSize(
        0.8 * h_energy_calibration_REC->GetYaxis()->GetLabelSize());
    h_energy_calibration_REC->GetXaxis()->SetTitleOffset(
        1.6 * h_energy_calibration_REC->GetXaxis()->GetTitleOffset());
    h_energy_calibration_REC->GetYaxis()->SetTitleOffset(
        2.5 * h_energy_calibration_REC->GetYaxis()->GetTitleOffset());
    h_energy_calibration_REC->GetYaxis()->SetTitle("counts");
    h_energy_calibration_REC->GetXaxis()->SetTitle("E_{reco} / E_{mc}");
    h_energy_calibration_REC->SetMarkerStyle(24);
    h_energy_calibration_REC->Draw("PEsame");

    c.cd(6);
    gPad->SetLogy(1);
    gPad->SetLeftMargin(0.13);
    gPad->SetBottomMargin(0.15);
    h_EoverP_REC->GetXaxis()->SetTitleSize(0.8 * h_EoverP_REC->GetXaxis()->GetTitleSize());
    h_EoverP_REC->GetXaxis()->SetLabelSize(0.8 * h_EoverP_REC->GetXaxis()->GetLabelSize());
    h_EoverP_REC->GetYaxis()->SetTitleSize(0.8 * h_EoverP_REC->GetYaxis()->GetTitleSize());
    h_EoverP_REC->GetYaxis()->SetLabelSize(0.8 * h_EoverP_REC->GetYaxis()->GetLabelSize());
    h_EoverP_REC->GetXaxis()->SetTitleOffset(1.6 * h_EoverP_REC->GetXaxis()->GetTitleOffset());
    h_EoverP_REC->GetYaxis()->SetTitleOffset(2.5 * h_EoverP_REC->GetYaxis()->GetTitleOffset());
    h_EoverP_REC->GetYaxis()->SetTitle("counts");
    h_EoverP_REC->SetMarkerStyle(24);
    h_EoverP_REC->Draw("PEsame");

    c.cd(7);
    gPad->SetLogy(1);
    gPad->SetLeftMargin(0.13);
    gPad->SetBottomMargin(0.15);
    h_ClusOverHit_REC->GetXaxis()->SetTitleSize(0.8 *
                                                h_ClusOverHit_REC->GetXaxis()->GetTitleSize());
    h_ClusOverHit_REC->GetXaxis()->SetLabelSize(0.8 *
                                                h_ClusOverHit_REC->GetXaxis()->GetLabelSize());
    h_ClusOverHit_REC->GetYaxis()->SetTitleSize(0.8 *
                                                h_ClusOverHit_REC->GetYaxis()->GetTitleSize());
    h_ClusOverHit_REC->GetYaxis()->SetLabelSize(0.8 *
                                                h_ClusOverHit_REC->GetYaxis()->GetLabelSize());
    h_ClusOverHit_REC->GetXaxis()->SetTitleOffset(1.6 *
                                                  h_ClusOverHit_REC->GetXaxis()->GetTitleOffset());
    h_ClusOverHit_REC->GetYaxis()->SetTitleOffset(2.5 *
                                                  h_ClusOverHit_REC->GetYaxis()->GetTitleOffset());
    h_ClusOverHit_REC->GetYaxis()->SetTitle("counts");
    h_ClusOverHit_REC->SetMarkerStyle(24);
    h_ClusOverHit_REC->Draw("PEsame");

    c.cd(8);
    gPad->SetLogz(1);
    gPad->SetLeftMargin(0.13);
    gPad->SetBottomMargin(0.15);
    h_emHits_position_REC->GetXaxis()->SetTitleSize(
        0.8 * h_emHits_position_REC->GetXaxis()->GetTitleSize());
    h_emHits_position_REC->GetXaxis()->SetLabelSize(
        0.8 * h_emHits_position_REC->GetXaxis()->GetLabelSize());
    h_emHits_position_REC->GetYaxis()->SetTitleSize(
        0.8 * h_emHits_position_REC->GetYaxis()->GetTitleSize());
    h_emHits_position_REC->GetYaxis()->SetLabelSize(
        0.8 * h_emHits_position_REC->GetYaxis()->GetLabelSize());
    h_emHits_position_REC->GetXaxis()->SetTitleOffset(
        1.6 * h_emHits_position_REC->GetXaxis()->GetTitleOffset());
    h_emHits_position_REC->GetYaxis()->SetTitleOffset(
        2. * h_emHits_position_REC->GetYaxis()->GetTitleOffset());
    h_emHits_position_REC->GetYaxis()->CenterTitle();
    h_emHits_position_REC->GetYaxis()->SetTitle("y(mm)");
    h_emHits_position_REC->Draw("colzsame");

    c.Print(fmt::format("{}_benchmark-{}-DIS-kinematics.pdf", output_prefix, vm_name).c_str());
  }

  return 0;
}
