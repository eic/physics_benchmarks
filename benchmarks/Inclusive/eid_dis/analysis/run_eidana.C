// Standalone eIDana benchmark wrapper.
// No dependency on ~/eicDIS — all source files are vendored in analysis/GlobalUtil/ and analysis/ElectronID/.
//
// Args:
//   file_list : path to text file with one stem per line
//               (stems come from xrdfs listing, without .eicrecon.edm4eic.root suffix)
//   output    : output ROOT file path
//   Ee        : electron beam energy (GeV)
//   Eh        : hadron beam energy (GeV)
//   campaign  : campaign string e.g. "26.05.0"
//   beam      : beam string e.g. "9x130"
//   q2_range  : q2 range string e.g. "q2_10to100"

#include "GlobalUtil/preLoadLib.hh"

#include <fstream>
#include <string>
#include <vector>
#include <cmath>

#include "TVector3.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLine.h"
#include "TLegend.h"
#include "TLorentzVector.h"

#include "podio/Frame.h"
#include "podio/ROOTReader.h"

#include "GlobalUtil/DrawManager.cc"
#include "GlobalUtil/getBoost.h"
#include "GlobalUtil/Constants.hh"
#include "ElectronID/ElectronID.cc"

// ---------- output tree variables (mirrors eIDana.h) ----------
TFile* outFile;
TTree* outTree;

int eID_status;
int eRecon_status;
int mc_PDG;
double EminusPz;
double EoP;

double mc_xB, mc_Q2, mc_y, mc_W2, mc_nu;
double rec_xB, rec_Q2, rec_W2, rec_y, rec_nu;

PxPyPzEVector vMC_e;
PxPyPzEVector vTRACK_e;
PxPyPzEVector vCLUSTER_e;
std::vector<PxPyPzEVector> vMC_hfs;
std::vector<PxPyPzEVector> vREC_hfs;

TH1D* h_nTPts_e;    TH1D* h_nTPts_jet_e; TH1D* h_nTPts_pi;  TH1D* h_nTPts_else;
TH1D* h_EoP_e;      TH1D* h_EoP_jet_e;   TH1D* h_EoP_pi;    TH1D* h_EoP_else;
TH1D* h_isoE_e;     TH1D* h_isoE_jet_e;  TH1D* h_isoE_pi;   TH1D* h_isoE_else;
TH1D* h_TrackEminusPz; TH1D* h_CalEminusPz;
TH1D* h_n_scat_elec;
TH2D* h_n_clusters_n_tracks;
TH1D* h_cand_mul; TH1D* h_cand_mul_eHighPt; TH1D* h_cand_mul_oHighPt;
TH1D* h_n_cluster_in_cone; TH1D* h_n_cluster_in_cone_found;

// ---------- helpers (verbatim from eIDana.C) ----------

void CalculateElectronKinematics(double fEe, double fEh, TLorentzVector kf,
    double& xB, double& Q2, double& W2, double& y, double& nu)
{
    TLorentzVector ki; ki.SetXYZM(0., 0., -fEe, MASS_ELECTRON);
    TLorentzVector P; P.SetX(fEh*sin(CROSSING_ANGLE)); P.SetY(0.);
    P.SetZ(fEh*cos(CROSSING_ANGLE)); P.SetE(std::hypot(fEh, MASS_PROTON));
    TLorentzVector q = ki - kf;
    Q2 = -(q.Dot(q));
    nu = (q.Dot(P)) / MASS_PROTON;
    xB = Q2 / (2. * MASS_PROTON * nu);
    y  = (q.Dot(P)) / (ki.Dot(P));
    W2 = MASS_PROTON*MASS_PROTON + (2.*MASS_PROTON*nu) - Q2;
}

void DefineHistograms()
{
    h_nTPts_e     = new TH1D("h_nTPts_e",     "N Track Points e; N_{TP}; Counts", 14, -0.5, 13.5);
    h_nTPts_jet_e = new TH1D("h_nTPts_jet_e", "N Track Points other e; N_{TP}; Counts", 14, -0.5, 13.5);
    h_nTPts_pi    = new TH1D("h_nTPts_pi",    "N Track Points #pi; N_{TP}; Counts", 14, -0.5, 13.5);
    h_nTPts_else  = new TH1D("h_nTPts_else",  "N Track Points others; N_{TP}; Counts", 14, -0.5, 13.5);
    h_EoP_e     = new TH1D("h_EoP_e",     "EoP e; E/p; Counts", 100, 0., 2.);
    h_EoP_jet_e = new TH1D("h_EoP_jet_e", "EoP other e; E/p; Counts", 100, 0., 2.);
    h_EoP_pi    = new TH1D("h_EoP_pi",    "EoP #pi; E/p; Counts", 100, 0., 2.);
    h_EoP_else  = new TH1D("h_EoP_else",  "EoP others; E/p; Counts", 100, 0., 2.);
    h_isoE_e     = new TH1D("h_isoE_e",     "Iso E e; Iso.E; Counts", 100, 0., 2.);
    h_isoE_jet_e = new TH1D("h_isoE_jet_e", "Iso E other e; Iso.E; Counts", 100, 0., 2.);
    h_isoE_pi    = new TH1D("h_isoE_pi",    "Iso E #pi; Iso.E; Counts", 100, 0., 2.);
    h_isoE_else  = new TH1D("h_isoE_else",  "Iso E others; Iso.E; Counts", 100, 0., 2.);
    h_TrackEminusPz = new TH1D("h_TrackEminusPz", "#Sigma(E-Pz) track; Counts", 200, 0., 50.);
    h_CalEminusPz   = new TH1D("h_CalEminusPz",   "#Sigma(E-Pz) cal; Counts",   200, 0., 50.);
    h_n_scat_elec       = new TH1D("h_n_scat_elec", "N scat e; N_{e}; Counts", 10, -0.5, 9.5);
    h_n_clusters_n_tracks = new TH2D("h_n_clusters_n_tracks", "N clusters vs N tracks; N_{trk}; N_{cls}", 5, -0.5, 4.5, 5, -0.5, 4.5);
    h_cand_mul        = new TH1D("h_cand_mul",        "e candidate mul; N; Counts", 10, -0.5, 9.5);
    h_cand_mul_eHighPt = new TH1D("h_cand_mul_eHighPt","e cand mul (e high pT); N; Counts", 10, -0.5, 9.5);
    h_cand_mul_oHighPt = new TH1D("h_cand_mul_oHighPt","e cand mul (other high pT); N; Counts", 10, -0.5, 9.5);
    h_n_cluster_in_cone       = new TH1D("h_n_cluster_in_cone",       "N clusters in cone; N; Counts", 20, -0.5, 19.5);
    h_n_cluster_in_cone_found = new TH1D("h_n_cluster_in_cone_found", "N clusters in cone (found e); N; Counts", 20, -0.5, 19.5);
}

void ResetVariables()
{
    eID_status = NO_MC; eRecon_status = NO_REC; mc_PDG = -999;
    EminusPz = -999; EoP = -999;
    mc_xB = mc_Q2 = mc_W2 = mc_y = mc_nu = -999;
    rec_xB = rec_Q2 = rec_W2 = rec_y = rec_nu = -999;
    vMC_e.SetPxPyPzE(0,0,0,0);
    vTRACK_e.SetPxPyPzE(0,0,0,0);
    vCLUSTER_e.SetPxPyPzE(0,0,0,0);
    vMC_hfs.clear(); vREC_hfs.clear();
}

void DrawVerticalLine(TCanvas* &c, double x_pos, double y_max)
{
    c->cd(); c->Modified(); c->Update();
    TLine* line = new TLine(x_pos, 0, x_pos, y_max);
    line->SetLineColor(kBlack); line->SetLineStyle(7); line->Draw("SAME");
}

void DrawParComparison(TCanvas* &c, TH1D* &h1, TH1D* &h2, TH1D* &h3, TH1D* &h4, double &draw_max)
{
    c->cd();
    h4->Draw("HIST"); h4->SetLineColor(kGray+2);
    draw_max = 1.2*std::max({h1->GetMaximum(), h2->GetMaximum(), h3->GetMaximum(), h4->GetMaximum()});
    h4->SetMaximum(draw_max);
    h3->Draw("HIST SAME"); h3->SetLineColor(kBlue);
    h2->Draw("HIST SAME"); h2->SetLineColor(kViolet);
    h1->Draw("HIST SAME"); h1->SetLineWidth(2); h1->SetLineColor(kRed); h1->SetFillColor(kRed); h1->SetFillStyle(3003);
    TLegend* leg = new TLegend(0.7, 0.6, 0.95, 0.88);
    leg->SetBorderSize(0); leg->SetFillStyle(0);
    leg->AddEntry(h1, "Electrons", "L"); leg->AddEntry(h2, "Other e's", "L");
    leg->AddEntry(h3, "Pions", "L");     leg->AddEntry(h4, "Others", "L");
    leg->Draw();
}

void DrawTCComparison(TCanvas* &c, TH1D* &ht, TH1D* &hc, double &draw_max)
{
    c->cd();
    hc->SetLineColor(kGray); hc->SetFillColor(kGray); hc->SetFillStyle(3003); hc->Draw("HIST");
    draw_max = 1.2*std::max(hc->GetMaximum(), ht->GetMaximum()); hc->SetMaximum(draw_max);
    ht->SetLineColor(kBlue); ht->SetFillColor(kBlue); ht->SetFillStyle(3003); ht->Draw("HIST SAME");
    TLegend* leg = new TLegend(0.6, 0.6, 0.8, 0.88);
    leg->SetBorderSize(0); leg->SetFillStyle(0);
    leg->AddEntry(ht, "Using E_{Track}", "L"); leg->AddEntry(hc, "Using E_{Cluster}", "L"); leg->Draw();
}

// ---------- main ----------

void run_eidana(const char* file_list, const char* output,
                int Ee, int Eh,
                const char* campaign, const char* beam, const char* q2_range)
{
    std::cout << "** eIDana benchmark: " << Ee << "x" << Eh
              << "  campaign=" << campaign
              << "  beam=" << beam
              << "  q2=" << q2_range << std::endl;

    // Build xrootd URLs from stem list
    std::vector<std::string> inFiles;
    std::ifstream fin(file_list);
    if (!fin.is_open()) { std::cerr << "Cannot open file list: " << file_list << std::endl; return; }
    std::string stem;
    while (std::getline(fin, stem)) {
        if (stem.empty()) continue;
        inFiles.push_back(
            std::string("root://epicxrd1.sdcc.bnl.gov:1095//eic/EPIC/RECO/") +
            campaign + "/epic_craterlake/DIS/pythia8.316-1.0/NC/noRad/ep/" +
            beam + "/" + q2_range + "/" + stem + ".eicrecon.edm4eic.root");
    }
    std::cout << "** " << inFiles.size() << " input files" << std::endl;

    // Podio
    podio::ROOTReader* reader = new podio::ROOTReader();
    reader->openFiles(inFiles);

    // Output
    outFile = new TFile(output, "RECREATE");
    outTree = new TTree("T_eID", "T_eID");
    outTree->Branch("eID_status",    &eID_status);
    outTree->Branch("eRecon_status", &eRecon_status);
    outTree->Branch("mc_PDG",        &mc_PDG);
    outTree->Branch("EminusPz",      &EminusPz);
    outTree->Branch("EoP",           &EoP);
    outTree->Branch("mc_xB",  &mc_xB);  outTree->Branch("mc_Q2",  &mc_Q2);
    outTree->Branch("mc_W2",  &mc_W2);  outTree->Branch("mc_y",   &mc_y);
    outTree->Branch("mc_nu",  &mc_nu);
    outTree->Branch("rec_xB", &rec_xB); outTree->Branch("rec_Q2", &rec_Q2);
    outTree->Branch("rec_W2", &rec_W2); outTree->Branch("rec_y",  &rec_y);
    outTree->Branch("rec_nu", &rec_nu);
    outTree->Branch("vMC_e",      &vMC_e);
    outTree->Branch("vTRACK_e",   &vTRACK_e);
    outTree->Branch("vCLUSTER_e", &vCLUSTER_e);
    outTree->Branch("vMC_hfs",    &vMC_hfs);
    outTree->Branch("vREC_hfs",   &vREC_hfs);

    // ElectronID
    ElectronID* eFinder = new ElectronID(Ee, Eh);
    eFinder->SetMinTrackPoints(4);
    LorentzRotation boost = getBoost(Ee, Eh, MASS_ELECTRON, MASS_PROTON);
    eFinder->SetBoost(boost);

    // Plots
    DrawManager* draw_manager = new DrawManager("ep", Form("%dx%d GeV", Ee, Eh), campaign);
    draw_manager->SetEPIC();

    DefineHistograms();

    std::cout << "** Starting event loop..." << std::endl;
    int countMCe = 0, countReconE = 0;

    for (size_t ev = 0; ev < reader->getEntries("events"); ev++) {
        auto raw = reader->readEntry("events", ev);
        if (!raw) { std::cerr << "null entry at event " << ev << "\n"; break; }
        podio::Frame event(std::move(raw));
        eFinder->SetEvent(&event);

        if (ev % 100 == 0)
            std::cout << "Event " << ev << "/" << reader->getEntries("events") << std::endl;

        edm4hep::MCParticleCollection e_mc = eFinder->GetMCElectron();
        h_n_scat_elec->Fill(e_mc.size());
        if (e_mc.size() > 0) {
            eID_status = FOUND_MC;
            TLorentzVector kprime;
            kprime.SetXYZM(e_mc[0].getMomentum().x, e_mc[0].getMomentum().y, e_mc[0].getMomentum().z, MASS_ELECTRON);
            CalculateElectronKinematics(Ee, Eh, kprime, mc_xB, mc_Q2, mc_W2, mc_y, mc_nu);
            vMC_e.SetPxPyPzE(e_mc[0].getMomentum().x, e_mc[0].getMomentum().y, e_mc[0].getMomentum().z, e_mc[0].getEnergy());
            countMCe += e_mc.size();
        }

        auto e_truth = eFinder->GetTruthReconElectron();
        if (e_truth.size() > 0) {
            eID_status = FOUND_TRUTH;
            h_n_clusters_n_tracks->Fill(e_truth[0].getTracks().size(), e_truth[0].getClusters().size());
            h_n_cluster_in_cone->Fill(eFinder->rcpart_n_clusters);
            if (e_truth[0].getClusters().size() > 0)
                h_n_cluster_in_cone_found->Fill(eFinder->rcpart_n_clusters);
            if (e_truth[0].getTracks().size() > 0 && e_truth[0].getClusters().size() > 0)
                eRecon_status = FOUND_BOTH;
            else if (e_truth[0].getTracks().size() > 0)
                eRecon_status = FOUND_TRACK_ONLY;
            else if (e_truth[0].getClusters().size() > 0)
                eRecon_status = FOUND_CLUSTER_ONLY;
        }

        auto e_candidates = eFinder->FindScatteredElectron();
        edm4eic::ReconstructedParticle e_rec;
        double TrackEminusPzSum = 0, CalEminusPzSum = 0;
        eFinder->GetEminusPzSum(TrackEminusPzSum, CalEminusPzSum);
        EminusPz = TrackEminusPzSum;

        if (e_candidates.size() > 0) {
            h_TrackEminusPz->Fill(TrackEminusPzSum);
            h_CalEminusPz->Fill(CalEminusPzSum);
            e_rec = eFinder->SelectHighestPT(e_candidates);
            EoP   = eFinder->GetCalorimeterEnergy(e_rec) / edm4hep::utils::magnitude(e_rec.getMomentum());
            vTRACK_e.SetPxPyPzE(e_rec.getMomentum().x, e_rec.getMomentum().y, e_rec.getMomentum().z, e_rec.getEnergy());
            vCLUSTER_e = eFinder->GetMomentumVectorFromCluster(e_rec, MASS_ELECTRON);

            mc_PDG = eFinder->Check_eID(e_rec);
            if (mc_PDG == 0)         eID_status = FOUND_E;
            else if (mc_PDG == -211) eID_status = FOUND_PI;
            else                     eID_status = FOUND_OTHERS;
            if (eID_status == FOUND_E) countReconE++;

            auto recoMC = eFinder->GetMC(e_rec);
            if (recoMC.isAvailable()) {
                TLorentzVector rk;
                rk.SetXYZM(recoMC.getMomentum().x, recoMC.getMomentum().y, recoMC.getMomentum().z, MASS_ELECTRON);
                CalculateElectronKinematics(Ee, Eh, rk, rec_xB, rec_Q2, rec_W2, rec_y, rec_nu);
            }

            h_cand_mul->Fill(e_candidates.size());
            if (mc_PDG == 0) h_cand_mul_eHighPt->Fill(e_candidates.size());
            else             h_cand_mul_oHighPt->Fill(e_candidates.size());

            for (const auto p : eFinder->GetMCHadronicFinalState())
                vMC_hfs.push_back({p.getMomentum().x, p.getMomentum().y, p.getMomentum().z, p.getEnergy()});
            for (const auto p : eFinder->FindHadronicFinalState(e_rec.getObjectID().index))
                vREC_hfs.push_back({p.getMomentum().x, p.getMomentum().y, p.getMomentum().z, p.getEnergy()});
        }

        for (const auto& d : eFinder->det_val) {
            if      (d.parType == 0)          { h_nTPts_e->Fill(d.nTrackPoints);     h_EoP_e->Fill(d.recon_EoP);     h_isoE_e->Fill(d.recon_isoE); }
            else if (d.parType == -211)        { h_nTPts_pi->Fill(d.nTrackPoints);    h_EoP_pi->Fill(d.recon_EoP);    h_isoE_pi->Fill(d.recon_isoE); }
            else if (abs(d.parType) == 11)     { h_nTPts_jet_e->Fill(d.nTrackPoints); h_EoP_jet_e->Fill(d.recon_EoP); h_isoE_jet_e->Fill(d.recon_isoE); }
            else                               { h_nTPts_else->Fill(d.nTrackPoints);  h_EoP_else->Fill(d.recon_EoP);  h_isoE_else->Fill(d.recon_isoE); }
        }

        outTree->Fill();
        ResetVariables();
    }

    std::cout << "** Done. eID rate: "
              << (countMCe > 0 ? (double)countReconE/countMCe*100 : 0) << "%" << std::endl;

    // Canvases
    double draw_max = 0.;

    TCanvas* c_nTPts = new TCanvas("c_nTPts", "c_nTPts", 1000, 600);
    DrawParComparison(c_nTPts, h_nTPts_e, h_nTPts_jet_e, h_nTPts_pi, h_nTPts_else, draw_max);
    DrawVerticalLine(c_nTPts, eFinder->GetMinTrackPoints()-0.5, draw_max);
    draw_manager->LableAndCollect(c_nTPts);

    TCanvas* c_EoP = new TCanvas("c_EoP", "c_EoP", 1000, 600);
    c_EoP->SetLogy();
    DrawParComparison(c_EoP, h_EoP_e, h_EoP_jet_e, h_EoP_pi, h_EoP_else, draw_max);
    DrawVerticalLine(c_EoP, eFinder->get_mEoP_min(), draw_max);
    DrawVerticalLine(c_EoP, eFinder->get_mEoP_max(), draw_max);
    draw_manager->LableAndCollect(c_EoP);

    TCanvas* c_isoE = new TCanvas("c_isoE", "c_isoE", 1000, 600);
    c_isoE->SetLogy();
    DrawParComparison(c_isoE, h_isoE_e, h_isoE_jet_e, h_isoE_pi, h_isoE_else, draw_max);
    DrawVerticalLine(c_isoE, eFinder->get_mIsoE(), draw_max);
    draw_manager->LableAndCollect(c_isoE);

    TCanvas* c_EminusPz = new TCanvas("c_EminusPz", "c_EminusPz", 1000, 600);
    DrawTCComparison(c_EminusPz, h_TrackEminusPz, h_CalEminusPz, draw_max);
    DrawVerticalLine(c_EminusPz, 2*Ee, draw_max);
    draw_manager->LableAndCollect(c_EminusPz);

    TCanvas* c_reco_mul = new TCanvas("c_reco_mul", "c_reco_mul", 1000, 600);
    c_reco_mul->SetLogy();
    h_cand_mul->Draw("HIST"); h_cand_mul->SetLineColor(kGray+2);
    h_cand_mul->GetYaxis()->SetRangeUser(1, h_cand_mul->GetMaximum()*1.5);
    h_cand_mul_eHighPt->Draw("HIST SAME"); h_cand_mul_eHighPt->SetLineColor(kBlue);
    h_cand_mul_oHighPt->Draw("HIST SAME"); h_cand_mul_oHighPt->SetLineColor(kOrange+7);
    TLegend* leg_mul = new TLegend(0.6, 0.6, 0.8, 0.88);
    leg_mul->SetBorderSize(0); leg_mul->SetFillStyle(0);
    leg_mul->AddEntry(h_cand_mul,        "All candidates", "L");
    leg_mul->AddEntry(h_cand_mul_eHighPt,"Scat. e highest p_{T}", "L");
    leg_mul->AddEntry(h_cand_mul_oHighPt,"Others highest p_{T}", "L");
    leg_mul->Draw();
    draw_manager->LableAndCollect(c_reco_mul, 2);

    TCanvas* c_ntc = new TCanvas("c_n_clusters_n_tracks", "c_n_clusters_n_tracks", 1000, 600);
    h_n_clusters_n_tracks->Scale(1.0/h_n_clusters_n_tracks->GetEntries());
    h_n_clusters_n_tracks->Draw("COLZ TEXT");
    draw_manager->LableAndCollect(c_ntc, 2);

    outFile->cd();
    outTree->Write(outTree->GetName(), 2);
    draw_manager->SaveToTree(outFile);
}
