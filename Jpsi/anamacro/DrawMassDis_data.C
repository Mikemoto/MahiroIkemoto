#include <vector>
#include <cmath>
#include <TLorentzVector.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TSpline.h>
#include <TVector2.h>
#include <TSystem.h>
#include <TParameter.h>

using namespace std;
TCanvas *c;

const float electron_mass = 0.000511; // GeV/c^2

// Matched track <-> calo pair (per matched track-calo pair)
struct MatchedTrkCalo
{
    int evt; // event id

    // track info
    float pt;
    float z0;
    float pz;
    float eta;
    float phi;
    int charge;
    int nmaps;
    int nintt;
    float chi2_trk;

    // track projection at EMCal
    float phi_emc_trk;
    float z_emc_trk;

    // calo cluster info
    float phi_calo;
    float z_calo;
    float e_calo;
    float chi2_calo;

    // for calo matching before correction
    float dphi_emc_raw; // track_phi_emc - calo_phi
    float dz_emc_raw;   // z_emc_trk   - calo_z

    // for calo matching after dphi correction
    float dphi_emc_corr;

    // calo matching after z-correction (will be filled after fit)
    float dz_emc_corr;

    // vertex (only data)
    float zvertex;

    // E/p
    float eop;
};

struct TrkCaloAllComb
{
    int evt;
    int it; // track index
    int ic; // calo index

    // track
    float pt;
    float pz;
    float eta;
    float phi;
    int charge;
    int nmaps;
    int nintt;
    float chi2_trk;

    float phi_emc_trk;
    float x_emc_trk;
    float y_emc_trk;
    float z_emc_trk;

    // calo
    float phi_calo;
    float z_calo;
    float e_calo;
    float chi2_calo; // data only (else -999)

    // matching variables (BEFORE selecting min)
    float dphi_emc_raw;  // wrapped [-pi, pi]
    float dphi_emc_corr; // wrapped [-pi, pi]
    float dz_emc;        // z_emc_trk - z_calo
    float dz_corr;       // if you want to store corrected one too
};

struct PairMass
{
    int evt;
    int ia; // index in matched_this_evt
    int ib;

    float mass;

    int q1, q2;
    float pt1, pt2;
    float eta1, eta2;
    float phi1, phi2;
    float eop1, eop2;

    bool pass_all;
    bool pass_no_eop;
    bool pass_middle_eop;
};

struct PassFlags
{
    bool pass_pt;
    bool pass_z0;
    bool pass_hits;
    bool pass_chi2;
    bool pass_dz;
    bool pass_dphi;
    bool pass_eop;
    bool pass_eop_middle;
    bool pass_vtx;
    bool pass_calomatch;
    bool pass_all;
    bool pass_no_eop;
    bool pass_middle_eop;
};

// set flags of cut value for mass reconstruction
PassFlags computeFlags(const MatchedTrkCalo &m,
                       bool isData,
                       double PtCut,
                       double Dz0VtxCut,
                       int nMapsCut,
                       int nInttCut,
                       double Chi2Cut,
                       double DzCut,
                       double phi_threshold,
                       double EopMinCut,
                       double EopMaxCut,
                       double ZvtxCut)
{
    PassFlags fl;

    // track quality cuts
    fl.pass_pt = (m.pt >= PtCut);
    fl.pass_z0 = isData ? (fabs(m.z0 - m.zvertex) < Dz0VtxCut) : true;
    fl.pass_hits = (m.nmaps == nMapsCut && m.nintt >= nInttCut);
    fl.pass_chi2 = (m.chi2_trk <= Chi2Cut);

    // Calo matching cuts (after dz correction)
    fl.pass_dz = (fabs(m.dz_emc_corr) < DzCut);
    fl.pass_dphi = (fabs(m.dphi_emc_corr) < phi_threshold);

    // E/p cuts
    fl.pass_eop = (m.eop > EopMinCut && m.eop < EopMaxCut);
    fl.pass_eop_middle = (m.eop > 0.6);

    // Vertex cut (for data only)
    fl.pass_vtx = isData ? (fabs(m.zvertex) <= ZvtxCut) : true;

    // Combined flags
    fl.pass_calomatch = (fl.pass_dz && fl.pass_dphi);
    fl.pass_all = (fl.pass_pt && fl.pass_z0 && fl.pass_hits && fl.pass_chi2 && fl.pass_dz && fl.pass_dphi && fl.pass_eop && fl.pass_vtx);
    fl.pass_no_eop = (fl.pass_pt && fl.pass_z0 && fl.pass_hits && fl.pass_chi2 && fl.pass_dz && fl.pass_dphi && fl.pass_vtx);
    fl.pass_middle_eop = (fl.pass_pt && fl.pass_z0 && fl.pass_hits && fl.pass_chi2 && fl.pass_dz && fl.pass_dphi && fl.pass_eop_middle && fl.pass_vtx);

    return fl;
}

float wrapPi(float dphi)
{
    while (dphi >= M_PI)
        dphi -= 2 * M_PI;
    while (dphi < -M_PI)
        dphi += 2 * M_PI;
    return dphi;
}

int recalcCharge(
    int trkid,
    int charge_old,
    const vector<int> *Siclus_trackid,
    const vector<float> *Siclus_x,
    const vector<float> *Siclus_y,
    float &phi_2nd_from_min,
    float &phi_max_from_min,
    float &dphi_out_2nd)
{
    if (!Siclus_trackid || !Siclus_x || !Siclus_y)
        return charge_old;

    float rmin1 = 1e9, rmin2 = 1e9, rmax = -1e9;
    float x_rmin1 = 0, y_rmin1 = 0, x_rmin2 = 0, y_rmin2 = 0, x_rmax = 0, y_rmax = 0;

    int nClus = (int)Siclus_trackid->size();
    int nMatch = 0;
    for (int ic = 0; ic < nClus; ic++)
    {
        if ((*Siclus_trackid)[ic] != trkid)
            continue;

        float x = (*Siclus_x)[ic];
        float y = (*Siclus_y)[ic];
        float r = sqrt(x * x + y * y);
        nMatch++;

        if (r < rmin1)
        {
            rmin2 = rmin1;
            x_rmin2 = x_rmin1;
            y_rmin2 = y_rmin1;
            rmin1 = r;
            x_rmin1 = x;
            y_rmin1 = y;
        }
        else if (r < rmin2)
        {
            rmin2 = r;
            x_rmin2 = x;
            y_rmin2 = y;
        }

        if (r > rmax)
        {
            rmax = r;
            x_rmax = x;
            y_rmax = y;
        }
    }

    if (nMatch < 3 || rmin2 > 9e8 || rmax < 0)
        return charge_old;

    phi_2nd_from_min = atan2(y_rmin2 - y_rmin1, x_rmin2 - x_rmin1);
    phi_max_from_min = atan2(y_rmax - y_rmin1, x_rmax - x_rmin1);

    dphi_out_2nd = wrapPi(phi_max_from_min - phi_2nd_from_min); // phi out(INTT) - in (MVTX 2nd inner)

    // if charge is + -> track curves to right side (clockwise) -> decrease phi -> dphi(out-in) < 0
    int charge_new = (dphi_out_2nd < 0 ? +1 : -1);

    return charge_new;
}

void DrawMassDis_data(const std::string &filename = "/sphenix/tg/tg01/commissioning/INTT/work/mahiro/SIliconCalo/MC/ana_jpsi/merged_200k.root", int condor_on = 0)
{
    // const std::string &filename = "/sphenix/tg/tg01/commissioning/INTT/work/mahiro/SIliconCalo/run24pp/macro/makedata/condor/data/ana_53879_00000.root";
    // const std::string &filename = "/sphenix/tg/tg01/commissioning/INTT/work/mahiro/SIliconCalo/run24pp/ana/olddst/ana_53879_50000.root";
    // const std::string &filename = "/sphenix/tg/tg01/commissioning/INTT/work/mahiro/SIliconCalo/MC/ana_jpsi/merged_200k.root";
    // const std::string &filename = "/sphenix/tg/tg01/commissioning/INTT/work/mahiro/SIliconCalo/run24pp/macro/makedata/condor/data/53879/ana_53879_00000.root";
    // const std::string &filename = "/sphenix/tg/tg01/commissioning/INTT/work/mahiro/SIliconCalo/MC/ana/merged_100m.root";
    // const std::string &filename = "/sphenix/tg/tg01/commissioning/INTT/work/mahiro/SIliconCalo/run24pp/ana/ana_53879_50000evt_test_newcalodst.root";
    // const std::string &filename = "/sphenix/tg/tg01/commissioning/INTT/work/mahiro/SIliconCalo/run24pp/ana/newgeo/calocalib/ana_53879_50000evt.root";
    // const std::string &filename = "/sphenix/tg/tg01/commissioning/INTT/work/mahiro/SIliconCalo/run24pp/ana/newgeo/calocalib/merged_1m.root";
    // const std::string &filename = "/sphenix/tg/tg01/commissioning/INTT/work/mahiro/SIliconCalo/run24pp/ana/newgeo/calocalib/recalc_charge/53879_merged_kouhan.root";
    // const std::string &filename = "/sphenix/tg/tg01/commissioning/INTT/work/mahiro/SIliconCalo/run24pp/ana/newgeo/calocalib/recalc_charge/53879_merged_zenhan.root";
    // const std::string &filename = "/sphenix/tg/tg01/commissioning/INTT/work/mahiro/SIliconCalo/run24pp/ana/newgeo/calocalib/53879_merged_5m.root";

    cout << "input file : " << filename.c_str() << endl;
    bool debug = false;
    // debug = true;

    bool matching_dphi = true; // dphi
    // bool matching_dphi = false; // dz matching

    bool use_dz_correction = false;

    int nloop = 500000;

    // const double Dz_preCut = 4.0;    // [cm] single gun sim
    // const double Dphi_preCut = 0.17; // [rad] single gun sim
    // const double Dz_preCut = 8.0;    // [cm]
    // const double Dphi_preCut = 0.25; // [rad]
    const double R_emc = 93.5; // [cm]
    // const double Wphi = 3.2;   // sim
    // double Wphi = 3.31; // run 53879
    double Wphi = 2.22; // single

    // const double Dz_preCut = 5.0;
    // const double phi_threshold = 10.;
    // const double phi_threshold = 0.062; // single
    const double phi_threshold = 0.165; // sim
    // const double phi_threshold = 0.178; // data
    const double PtCut = 0.5;
    // const double PtCut = 1.;
    const double Dz0VtxCut = 1.;
    const int nMapsCut = 3;
    const int nInttCut = 1;
    // const double Chi2Cut = 4.0513; // single
    const double Chi2Cut = 7.415; // sim pythia
    // const double Chi2Cut = 22.781; // data_
    // const double DzCut = 2.6; // single
    const double DzCut = 4.812; // sim
    // const double DzCut = 5.016; // data
    const double ZvtxCut = 10.; // cm
    const double EopMinCut = 0.8;
    const double EopMaxCut = 1.2;

    if (debug)
    {
        cout << "Cut status..." << endl;
        cout << "Track pt > " << PtCut << ", nmaps > " << nMapsCut << ", nintt > " << nInttCut << ", chi2 < " << Chi2Cut << ", dz_emc < " << DzCut << ", " << EopMinCut << "< eop < " << EopMaxCut << "|zvtx| <= " << ZvtxCut << endl;
    }

    // cout << "^^bbb" << endl;

    TFile *f = TFile::Open(filename.c_str(), "READ");
    if (!f || f->IsZombie())
    {
        std::cerr << "Error: cannot open file: " << filename.c_str() << std::endl;
        return;
    }

    TTree *trackTree = (TTree *)f->Get("trackTree");
    TTree *caloTree = (TTree *)f->Get("caloTree");
    TTree *evtTree = (TTree *)f->Get("evtTree");
    TTree *SiClusTree = (TTree *)f->Get("SiClusTree");
    if (!trackTree || !caloTree)
    {
        std::cerr << "Error: missing tree(s) in file" << std::endl;
        if (!trackTree)
            std::cerr << "  - trackTree not found" << std::endl;
        if (!caloTree)
            std::cerr << "  - caloTree  not found" << std::endl;

        f->Close();
        delete f;
        return;
    }
    //  TTree *truthTree = (TTree *)f->Get("truthTree");
    bool IS_DATA = (evtTree != nullptr);

    // ================== dphi alignment spline load ==================
    TSpline3 *sp_mu = nullptr;
    TFile *fAlign = TFile::Open("result/data/alignment/dphi_phi0_alighnment.root", "READ");
    sp_mu = (TSpline3 *)fAlign->Get("spline_mu_vs_phi0");
    if (!sp_mu)
    {
        std::cerr << "[dphi alignment] spline_mu_vs_phi0 not found. Run without alignment." << std::endl;
    }
    // ===============================================================

    TString fname_out = filename.substr(filename.find_last_of("/"), filename.size());
    // fname_out.ReplaceAll("merged", "mass_reconstruction_PYTHIA"); //
    fname_out.ReplaceAll("merged", "mass_reco"); //
    if (debug)
        fname_out.ReplaceAll(".root", Form("_debug_%devt.root", nloop));

    fname_out = (IS_DATA ? "result/data" : "result/sim") + fname_out;
    if (condor_on && IS_DATA)
    {
        TString run = gSystem->BaseName(filename.c_str()); // "ana_53876_00001.root"
        run.ReplaceAll("ana_", "");
        run = run(0, run.First("_")); // "53876"
        gSystem->mkdir(Form("result/ana_condor/%s", run.Data()), true);
        fname_out.ReplaceAll("result/data", Form("result/ana_condor/%s", run.Data()));
    }

    if (condor_on && !IS_DATA)
    {
        // fname_out.ReplaceAll("ana_", "ana_00");
        fname_out.ReplaceAll("sim", "sim/ana_condor");
    }

    Wphi = IS_DATA ? 3.31 : 3.2; // 😺

    TFile *outFile = new TFile(fname_out, "RECREATE");

    TString fname_out_pdf = fname_out;
    fname_out_pdf.ReplaceAll(".root", ".pdf");

    c = new TCanvas("c", "c", 1000, 500);
    c->Divide(2, 1);
    c->cd(1);
    c->Print(fname_out_pdf + "[");

    gStyle->SetOptStat(0); // 📦
    c->SetLeftMargin(0.20);

    // output tree of all cinbination (when debug is true)
    TTree *t_allcomb = nullptr;
    TrkCaloAllComb ac;

    if (debug)
    {
        t_allcomb = new TTree("t_allcomb", "All track-calo combinations before matching");

        t_allcomb->Branch("evt", &ac.evt, "evt/I");
        t_allcomb->Branch("it", &ac.it, "it/I");
        t_allcomb->Branch("ic", &ac.ic, "ic/I");

        t_allcomb->Branch("pt", &ac.pt, "pt/F");
        t_allcomb->Branch("pz", &ac.pz, "pz/F");
        t_allcomb->Branch("eta", &ac.eta, "eta/F");
        t_allcomb->Branch("phi", &ac.phi, "phi/F");
        t_allcomb->Branch("charge", &ac.charge, "charge/I");
        t_allcomb->Branch("nmaps", &ac.nmaps, "nmaps/I");
        t_allcomb->Branch("nintt", &ac.nintt, "nintt/I");
        t_allcomb->Branch("chi2_trk", &ac.chi2_trk, "chi2_trk/F");

        t_allcomb->Branch("phi_emc_trk", &ac.phi_emc_trk, "phi_emc_trk/F");
        t_allcomb->Branch("x_emc_trk", &ac.x_emc_trk, "x_emc_trk/F");
        t_allcomb->Branch("y_emc_trk", &ac.y_emc_trk, "y_emc_trk/F");
        t_allcomb->Branch("z_emc_trk", &ac.z_emc_trk, "z_emc_trk/F");

        t_allcomb->Branch("phi_calo", &ac.phi_calo, "phi_calo/F");
        t_allcomb->Branch("z_calo", &ac.z_calo, "z_calo/F");
        t_allcomb->Branch("e_calo", &ac.e_calo, "e_calo/F");
        t_allcomb->Branch("chi2_calo", &ac.chi2_calo, "chi2_calo/F");

        t_allcomb->Branch("dphi_emc_raw", &ac.dphi_emc_raw, "dphi_emc_raw/F");
        t_allcomb->Branch("dphi_emc_corr", &ac.dphi_emc_corr, "dphi_emc_corr/F");
        t_allcomb->Branch("dz_emc", &ac.dz_emc, "dz_emc/F");
        t_allcomb->Branch("dz_corr", &ac.dz_corr, "dz_corr/F");
    }

    TTree *t_match = nullptr;
    MatchedTrkCalo br;

    if (debug && IS_DATA)
    {
        t_match = new TTree("t_match", "Matched track-calo pairs with dz-correction");

        t_match->Branch("evt", &br.evt, "evt/I");
        t_match->Branch("pt", &br.pt, "pt/F");
        t_match->Branch("pz", &br.pz, "pz/F");
        t_match->Branch("eta", &br.eta, "eta/F");
        t_match->Branch("phi", &br.phi, "phi/F");
        t_match->Branch("charge", &br.charge, "charge/I");
        t_match->Branch("nmaps", &br.nmaps, "nmaps/I");
        t_match->Branch("nintt", &br.nintt, "nintt/I");
        t_match->Branch("chi2_trk", &br.chi2_trk, "chi2_trk/F");

        t_match->Branch("phi_emc_trk", &br.phi_emc_trk, "phi_emc_trk/F");
        t_match->Branch("z_emc_trk", &br.z_emc_trk, "z_emc_trk/F");

        t_match->Branch("phi_calo", &br.phi_calo, "phi_calo/F");
        t_match->Branch("z_calo", &br.z_calo, "z_calo/F");
        t_match->Branch("e_calo", &br.e_calo, "e_calo/F");
        t_match->Branch("chi2_calo", &br.chi2_calo, "chi2_calo/F");

        t_match->Branch("dphi_emc_raw", &br.dphi_emc_raw, "dphi_emc_raw/F");
        t_match->Branch("dz_emc_raw", &br.dz_emc_raw, "dz_emc_raw/F");
        t_match->Branch("dz_emc_corr", &br.dz_emc_corr, "dz_emc_corr/F");

        t_match->Branch("eop", &br.eop, "eop/F");
        t_match->Branch("zvertex", &br.zvertex, "zvertex/F");
    }

    // --- t_pair:  ---
    TTree *t_pair = new TTree("t_pair", "e+e- pair tree");
    PairMass pm;

    t_pair->Branch("evt", &pm.evt, "evt/I");
    t_pair->Branch("ia", &pm.ia, "ia/I");
    t_pair->Branch("ib", &pm.ib, "ib/I");
    t_pair->Branch("mass", &pm.mass, "mass/F");

    t_pair->Branch("q1", &pm.q1, "q1/I");
    t_pair->Branch("q2", &pm.q2, "q2/I");
    t_pair->Branch("pt1", &pm.pt1, "pt1/F");
    t_pair->Branch("pt2", &pm.pt2, "pt2/F");
    t_pair->Branch("eta1", &pm.eta1, "eta1/F");
    t_pair->Branch("eta2", &pm.eta2, "eta2/F");
    t_pair->Branch("phi1", &pm.phi1, "phi1/F");
    t_pair->Branch("phi2", &pm.phi2, "phi2/F");
    t_pair->Branch("eop1", &pm.eop1, "eop1/F");
    t_pair->Branch("eop2", &pm.eop2, "eop2/F");

    t_pair->Branch("pass_all", &pm.pass_all, "pass_all/O");
    t_pair->Branch("pass_no_eop", &pm.pass_no_eop, "pass_no_eop/O");
    t_pair->Branch("pass_middle_eop", &pm.pass_middle_eop, "pass_middle_eop/O");

    // input branches
    int evt = 0;
    vector<int> *track_id = nullptr;
    vector<float> *track_phi = nullptr;
    vector<float> *track_pt = nullptr;
    vector<float> *track_pz = nullptr;
    vector<float> *track_eta = nullptr;
    vector<float> *track_z = nullptr;
    vector<int> *track_nmaps = nullptr;
    vector<int> *track_nintt = nullptr;
    vector<float> *track_phi_emc = nullptr;
    vector<float> *track_x_emc = nullptr;
    vector<float> *track_y_emc = nullptr;
    vector<float> *track_z_emc = nullptr;
    vector<int> *track_charge = nullptr;
    vector<float> *track_chi2ndf = nullptr;

    trackTree->SetBranchAddress("evt", &evt);
    trackTree->SetBranchAddress("track_id", &track_id);
    trackTree->SetBranchAddress("phi0", &track_phi);
    trackTree->SetBranchAddress("pt0", &track_pt);
    trackTree->SetBranchAddress("pz0", &track_pz);
    trackTree->SetBranchAddress("eta0", &track_eta);
    trackTree->SetBranchAddress("z0", &track_z);
    trackTree->SetBranchAddress("nmaps", &track_nmaps);
    trackTree->SetBranchAddress("nintt", &track_nintt);
    trackTree->SetBranchAddress("phi_proj_emc", &track_phi_emc);
    trackTree->SetBranchAddress("x_proj_emc", &track_x_emc);
    trackTree->SetBranchAddress("y_proj_emc", &track_y_emc);
    trackTree->SetBranchAddress("z_proj_emc", &track_z_emc);
    trackTree->SetBranchAddress("charge", &track_charge);
    trackTree->SetBranchAddress("chi2ndf", &track_chi2ndf);

    int calo_evt = 0;
    vector<float> *calo_phi = nullptr;
    vector<float> *calo_energy = nullptr;
    vector<float> *calo_x = nullptr;
    vector<float> *calo_y = nullptr;
    vector<float> *calo_z = nullptr;
    vector<float> *calo_chi2 = nullptr;

    caloTree->SetBranchAddress("calo_evt", &calo_evt);
    caloTree->SetBranchAddress("phi", &calo_phi);
    caloTree->SetBranchAddress("energy", &calo_energy);
    caloTree->SetBranchAddress("x", &calo_x);
    caloTree->SetBranchAddress("y", &calo_y);
    caloTree->SetBranchAddress("z", &calo_z);
    if (IS_DATA)
        caloTree->SetBranchAddress("chi2", &calo_chi2);

    float_t xvtx = 0.0;
    float_t yvtx = 0.0;
    float_t zvtx = 0.0;

    if (IS_DATA)
    {
        evtTree->SetBranchAddress("xvtx", &xvtx);
        evtTree->SetBranchAddress("yvtx", &yvtx);
        evtTree->SetBranchAddress("zvtx", &zvtx);
    }

    int clus_evt = 0;
    vector<int> *Siclus_trackid = nullptr;
    vector<int> *Siclus_layer = nullptr;
    vector<float> *Siclus_x = nullptr;
    vector<float> *Siclus_y = nullptr;
    vector<float> *Siclus_z = nullptr;
    vector<int> *Siclus_t = nullptr;

    SiClusTree->SetBranchAddress("evt", &clus_evt);
    SiClusTree->SetBranchAddress("Siclus_trackid", &Siclus_trackid);
    SiClusTree->SetBranchAddress("Siclus_layer", &Siclus_layer);
    SiClusTree->SetBranchAddress("Siclus_x", &Siclus_x);
    SiClusTree->SetBranchAddress("Siclus_y", &Siclus_y);
    SiClusTree->SetBranchAddress("Siclus_z", &Siclus_z);
    SiClusTree->SetBranchAddress("Siclus_t", &Siclus_t);

    // takusan hists
    TH1F *h_dphi = new TH1F("h_dphi", "Track - Calo #Delta#phi;#Delta#phi;Counts", 200, -0.3, 0.3);
    TH1F *h_dphi_emc = new TH1F("h_dphi (EMCal Proj)", "Track - Calo #Delta#phi;#Delta#phi;Counts", 200, -0.3, 0.3);
    TH2F *h_dphi_emc_pt = new TH2F("h_dphi_pt", "Track - Calo #Delta#phi;pT;#Delta#phi", 200, 0, 20, 200, -0.3, 0.3);
    TH1F *h_dz = new TH1F("h_dz", "Track - Calo #Delta z;#Delta z (cm);Counts", 200, -50, 50);
    TH1F *h_dz_emc = new TH1F("h_dz_emc", "Track(EMCal Proj) - Calo #Delta z;#Delta z (cm);Counts", 200, -50, 50);
    TH1F *h_eta = new TH1F("h_eta", "1/N_{evt} dN/d#eta (Calo matched);#eta;1/N_{evt} dN/d#eta", 150, -3, 3);
    TH1F *h_eta_mass = new TH1F("h_eta_mass", "1/N_{evt} dN/d#eta (used for mass dis);#eta;1/N_{evt} dN/d#eta", 150, -3, 3);

    TH1F *h_mass = new TH1F("h_mass", "Invariant mass of matched track pairs (e^{+}e^{-});Mass (GeV/c^{2});Counts", 200, 0, 5);

    // Invariant mass histograms for different J/psi pT bins
    TH1F *h_track_chi2ndf_matched = new TH1F("h_track_chi2ndf_matched", "Chi2/NDF of matched tracks;#chi^{2}/ndf;Counts", 100, 0, 10);

    TH2F *h_reco_vs_truth_pt = new TH2F("h_reco_vs_truth_pt", "Reconstructed pT vs Truth pT;Truth pT (GeV/c);Reconstructed pT (GeV/c)", 100, 0, 20, 100, 0, 20);

    TH1F *h_zv = new TH1F("h_zv", "z vertex;[cm];Counts", 100, -50, 50);
    TH2F *h_xv_yv = new TH2F("h_xv_yv", "x vertex vs y vertex;[cm];counts", 20, -10, 10, 20, -10, 10);

    // hist of mass distribution with only each cut
    TH1F *h_mass_base = new TH1F("h_mass_base", "Mass distribution;M (GeV/c^{2});Pairs", 200, 0, 5);
    TH1F *h_mass_only_pt = new TH1F("h_mass_only_pt", "Mass (only pT cut);M (GeV/c^{2});Pairs", 200, 0, 5);
    TH1F *h_mass_only_z0 = new TH1F("h_mass_only_z0", "Mass (only z0 cut);M (GeV/c^{2});Pairs", 200, 0, 5);
    TH1F *h_mass_only_hits = new TH1F("h_mass_only_hits", "Mass (only hits cut);M (GeV/c^{2});Pairs", 200, 0, 5);
    TH1F *h_mass_only_chi2 = new TH1F("h_mass_only_chi2", "Mass (only #chi^{2} cut);M (GeV/c^{2});Pairs", 200, 0, 5);
    TH1F *h_mass_only_dz = new TH1F("h_mass_only_dz", "Mass (only |#Delta z_{EMCal}| cut);M (GeV/c^{2});Pairs", 200, 0, 5);
    TH1F *h_mass_only_dphi = new TH1F("h_mass_only_dphi", "Mass (only |#Delta#phi_{EMCal}| cut);M (GeV/c^{2});Pairs", 200, 0, 5);
    TH1F *h_mass_only_zvtx = new TH1F("h_mass_only_zvtx", "Mass (only z_{vertex} cut);M (GeV/c^{2});Pairs", 200, 0, 5);
    TH1F *h_mass_only_eop = new TH1F("h_mass_only_eop", "Mass (only E/p cut);M (GeV/c^{2});Pairs", 200, 0, 5);
    TH1F *h_mass_allcuts = new TH1F("h_mass_allcuts", "Mass (ALL cuts: pT & hits & #chi^{2} & |#Delta z| & |#Delta#phi| & E/p);M (GeV/c^{2});Pairs", 100, 0, 5);
    TH1F *h_mass_no_eop = new TH1F("h_mass_no_eop", "Mass (ALL cuts: pT & hits & #chi^{2} & |#Delta z| & |#Delta#phi|);M (GeV/c^{2});Pairs", 100, 0, 5);
    TH1F *h_mass_middle_eop = new TH1F("h_mass_middle_eop", "Mass (ALL cuts: pT & hits & #chi^{2} & |#Delta z| & |#Delta#phi| & middle E/p);M (GeV/c^{2});Pairs", 100, 0, 5);

    // mass distribution of same charge particles for check BG
    TH1F *h_mass_ss_middle_eop = new TH1F("h_mass_ss_middle_eop", "Samecharge sign mass (e^{#pm}e^{#pm}) (E/p > 0.6);M (GeV/c^{2});Pairs", 100, 0, 5);
    TH1F *h_mass_ss_allcuts = new TH1F("h_mass_ss_allcuts", "Samecharge sign mass (e^{#pm}e^{#pm}) (ALL cuts);M (GeV/c^{2});Pairs", 100, 0, 5);

    // mass distribution of all charge comb for check BG
    TH1F *h_mass_allc_middle_eop = new TH1F("h_mass_allc_middle_eop", "all cahrge comb mass (e^{#pm}e^{#pm}) (E/p > 0.6);M (GeV/c^{2});Pairs", 100, 0, 5);
    TH1F *h_mass_allc_allcuts = new TH1F("h_mass_allc_allcuts", "all cahrge comb mass (e^{#pm}e^{#pm}) (ALL cuts);M (GeV/c^{2});Pairs", 100, 0, 5);

    // hists of before/after each cut
    TH1F *h_pt_all = new TH1F("h_pt_all", "Track pt;pt [GeV/c];Tracks", 100, 0, 10);
    TH1F *h_pt_cut = new TH1F("h_pt_cut", "Track p_{T} (with only pT cut);p_{T} [GeV/c];Tracks", 100, 0, 10);
    TH1F *h_pt_chi2 = new TH1F("h_pt_chi2", "Track p_{T} after #chi^{2}/ndf cut;p_{T} [GeV/c];Tracks", 100, 0, 10);
    TH1F *h_pt_chi2_eop = new TH1F("h_pt_chi2_eop", "Track p_{T} after #chi^{2}/ndf + E/p cut;p_{T} [GeV/c];Tracks", 100, 0, 10);


    TH1F *h_dz0v_all = new TH1F("h_dz0v_all", "z0 - zvertex;z0 - zvertex [cm];Tracks", 100, -5, 5);
    TH1F *h_dz0v_cut = new TH1F("h_dz0v_cut", "z0 - zvertex (with only z0 - zv cut);z0 - zvertex [cm];Tracks", 100, -5, 5);

    TH1F *h_dz_raw = new TH1F("h_dz_raw", "#Delta z_{EMCal} (all trkr with all calo combination);#Delta z [cm];Tracks", 100, -25, 25);
    TH1F *h_dz_all = new TH1F("h_dz_all", "#Delta z_{EMCal} ;#Delta z [cm];Tracks", 100, -25, 25);
    TH1F *h_dz_corr = new TH1F("h_dz_corr", "#Delta z_{EMCal} (after correlation) ;#Delta z [cm];Tracks", 100, -25, 25);
    TH1F *h_dz_cut = new TH1F("h_dz_cut", "#Delta z_{EMCal} (with only |dz| cut);#Delta z [cm];Tracks", 100, -25, 25);

    TH1F *h_dphi_raw = new TH1F("h_dphi_raw", "#Delta#phi_{EMCal} (all trkr with all calo combination) ;#Delta#phi;Tracks", 100, -1, 1);
    TH1F *h_dphi_all = new TH1F("h_dphi_all", "#Delta#phi_{EMCal} ;#Delta#phi;Tracks", 100, -1, 1);
    TH1F *h_dphi_raw_cm = new TH1F("h_dphi_raw_cm", "#Delta#phi_{EMCal} ;#Delta#phi [cm];Tracks", 150, -60, 60);
    TH1F *h_dphi_corr_cm = new TH1F("h_dphi_corr_cm", "#Delta#phi_{EMCal} (corr, all combination) ;#Delta#phi [cm];Tracks", 150, -60, 60);
    TH1F *h_dphi_cut = new TH1F("h_dphi_cut", "#Delta#phi_{EMCal} (with only |d#phi| cut);|#Delta#phi|;Tracks", 100, -1, 1);

    TH1F *h_eop_all = new TH1F("h_eop_all", "E/p;E/p;Tracks", 100, 0, 5);
    TH1F *h_eop_cut = new TH1F("h_eop_cut", Form("E/p (with only %f <E/p< %f cut);E/p;Tracks", EopMinCut, EopMaxCut), 100, 0, 5);
    TH1F *h_eop_masscut = new TH1F("h_eop_masscut", "E/P (only mass ir around 3);E/p;Tracks", 100, 0, 5);

    TH1F *h_chi2_all = new TH1F("h_chi2_all", "#chi^{2}/ndf;#chi^{2}/ndf;Tracks", 300, 0, 30);
    TH1F *h_chi2_cut = new TH1F("h_chi2_cut", "#chi^{2}/ndf (with only cut);#chi^{2}/ndf;Tracks", 300, 0, 30);

    TH1F *h_zv_all = new TH1F("h_zv_all", "z_{vertex} distribution;[cm];Counts", 100, -50, 50);
    TH1F *h_zv_cut = new TH1F("h_zv_cut", "z_{vertex} distribution (with cut);[cm];Counts", 100, -50, 50);

    TH2F *h_nhits_all = new TH2F("h_nhits_all", "nintt vs nmaps;nmaps;nintt", 6, 0, 6, 6, 0, 6);
    TH2F *h_nhits_cut = new TH2F("h_nhits_cut", "nintt vs nmaps (with only hits cut);nmaps;nintt", 6, 0, 6, 6, 0, 6);

    TH2F *h_eop_vs_mass = new TH2F("h_eop_vs_mass", "E/p vs mass;mass;E/p", 200, 0, 5, 200, 0, 5);
    TH2F *h_eop_vs_mass_pass = new TH2F("h_eop_vs_mass_pass", "E/p vs mass (with cut);mass;E/p", 200, 0, 5, 200, 0, 5);
    TH2F *h_eta_vs_phi = new TH2F("h_eta_vs_phi", "eta vs #phi;#phi [rad];eta", 200, -0.3, 0.3, 150, -3, 3);
    TH2F *h_eta_vs_pt = new TH2F("h_eta_vs_pt", "eta vs Pt;Pt [GeV];eta", 100, 0, 10, 150, -3, 3);

    TH2F *h_dphi_vs_pt = new TH2F("h_dphi_vs_pt", "dphi vs Pt;Pt [GeV];dphi", 100, 0, 10, 100, -1, 1);
    TH2F *h_dz_vs_pt = new TH2F("h_dz_vs_pt", "dz vs Pt;Pt [GeV];dz", 100, 0, 10, 200, -10, 10);
    TH2F *h_eop_vs_pt = new TH2F("h_eop_vs_pt", "E/p vs Pt;Pt [GeV];E/p", 100, 0, 10, 150, 0, 5);

    TH2F *h_dphi_vs_tkgphi = new TH2F("h_dphi_vs_tkgphi", "dphi vs phi (silicon tracking) ;phi [rad] ;dphi [rad]", 314, -3.14, 3.14, 100, -1, 1);
    TH2F *h_dphiraw_vs_tkgphi = new TH2F("h_dphiraw_vs_tkgphi", "dphi (all combination) vs phi (silicon tracking);phi [rad] ;dphi [rad]", 314, -3.14, 3.14, 100, -1, 1);
    TH2F *h_dphi_vs_calophi = new TH2F("h_dphi_vs_calophi", "dphi vs phi (Calo cluster);phi [rad] ;dphi [rad]", 314, -3.14, 3.14, 100, -1, 1);
    TH2F *h_dphiraw_vs_calophi = new TH2F("h_dphiraw_vs_calophi", "dphi (all caobination) vs phi (Calo cluster);phi [rad];dphi [rad]", 314, -3.14, 3.14, 100, -1, 1);

    TH2F *h_dphiraw_vs_phi0 = new TH2F("h_dphiraw_vs_phi0", "dphi (all caobination) vs phi0 (phi of track start point);phi0 [rad];dphi [rad]", 314, -3.14, 3.14, 100, -1, 1);
    TH2F *h_dphicorr_vs_phi0 = new TH2F("h_dphicorr_vs_phi0", "dphi (corr, all combination) vs phi0;phi0 [rad];dphi_{corr} [rad]", 314, -3.14, 3.14, 100, -1, 1);

    TH2F *h_trkrphi_vs_calophi = new TH2F("h_trkrphi_vs_calophi", "trkr phi vs calo phi ;trkr phi [rad];calo phi [rad]", 314, -3.14, 3.14, 314, -3.14, 3.14);
    TH2F *h_trkrz_vs_caloz = new TH2F("h_trkrz_vs_caloz", "trkr z vs calo z ;trkr z [cm];calo z [cm]", 130, -130, 130, 130, -130, 130);

    TH2F *h_dzraw_vs_caloz = new TH2F("h_dzraw_vs_caloz", "all dz vs calo z ;calo z [cm] ;dz [cm]", 200, -50, 50, 200, -50, 50);
    TH2F *h_dzraw_vs_trkrz = new TH2F("h_dzraw_vs_trkrz", "all dz vs trkr z ;trkr z [cm] ; dz [cm]", 200, -50, 50, 200, -50, 50);

    TH2F *h_dzraw_vs_z0 = new TH2F("h_dzraw_vs_z0", "all dz vs z0 (z of track start point) ;z0 [cm] ; dz [cm]", 200, -20, 20, 200, -50, 50);

    TH2F *h_dz_vs_trkrz = new TH2F("h_dz_vs_trkrz", "dz vs trkr z ;trkr z [cm] ; dz [cm]", 200, -50, 50, 200, -50, 50);
    TH2F *h_dz_vs_caloz = new TH2F("h_dz_vs_caloz", "dz vs calo z ;calo z [cm] ;dz [cm]", 200, -50, 50, 200, -50, 50);
    TH2F *h_dzcorr_vs_trkrz = new TH2F("h_dzcorr_vs_trkrz", "correlated dz vs trkr z ;dz [cm];trkr z [cm]", 200, -50, 50, 100, -50, 50);
    TH2F *h_dzcorr_vs_caloz = new TH2F("h_dzcorr_vs_caloz", "correlated dz vs calo z ;dz [cm];calo z [cm]", 200, -50, 50, 100, -50, 50);

    TH2F *h_dz_vs_dphi = new TH2F("h_dz_vs_dphi", "dz vs dphi;dz_{EMCal} [cm];dphi_{EMCal} [cm]", 400, -100, 100, 400, -100, 100);

    // vector to buffer all matched track-calo pairs (for all events)
    // vector<vector<MatchedTrkCalo>> matched_evt;
    // vector<vector<PassFlags>> flags_evt;

    Long64_t nentries = trackTree->GetEntries();
    // gStyle->SetOptStat("ne");
    TProfile *pf_trkrz_vs_caloz = nullptr;
    TF1 *fit_trkrz_vs_caloz = nullptr;
    double p0, p1;
    Long64_t n_evt_all = 0;     // all events
    Long64_t n_evt_hasPair = 0; // # of events which can make at least one pair of e+e-

    if (use_dz_correction)
    {
        // ================== z correction without matching ==================
        for (Long64_t i = 0; i < nentries; i++)
        {
            trackTree->GetEntry(i);
            caloTree->GetEntry(i);

            if (evt != calo_evt)
            {
                std::cerr << "Warning: evt mismatch at entry " << i
                          << ": track evt = " << evt
                          << ", calo evt = " << calo_evt << std::endl;
                continue;
            }

            int nTrk = (int)track_z_emc->size();
            int nCalo = (int)calo_z->size();

            for (int it = 0; it < nTrk; it++)
            {
                float z_emc_trk = (*track_z_emc)[it];
                for (int ic = 0; ic < nCalo; ic++)
                {
                    float z_calo = (*calo_z)[ic];
                    // trkr with all calo (without matching)
                    if (fabs(1.13 * z_emc_trk - z_calo) < 8)
                        h_trkrz_vs_caloz->Fill(z_emc_trk, z_calo);
                }
            }

            if (debug && i + 1 >= nloop)
                break;
        }

        // ========== fitting trkr z vs calo z to get slope ==========
        gStyle->SetOptFit(1111);
        pf_trkrz_vs_caloz = h_trkrz_vs_caloz->ProfileX("pf_trkrz_vs_caloz", 1, -1, "s");
        fit_trkrz_vs_caloz = new TF1("fit_trkrz_vs_caloz", "pol1", -100, 100);
        pf_trkrz_vs_caloz->Fit(fit_trkrz_vs_caloz, "RS");

        p0 = fit_trkrz_vs_caloz->GetParameter(0); // calo_z = p0 + p1 * trkr_z
        p1 = fit_trkrz_vs_caloz->GetParameter(1);

        h_trkrz_vs_caloz->GetListOfFunctions()->Add(fit_trkrz_vs_caloz);
    }

    // ================== calo matching, set cut flag ==================

    for (Long64_t i = 0; i < nentries; ++i)
    {
        vector<MatchedTrkCalo> matched_this_evt;
        vector<PassFlags> flags_this_evt;
        trackTree->GetEntry(i);
        caloTree->GetEntry(i);
        SiClusTree->GetEntry(i);

        if (IS_DATA)
            evtTree->GetEntry(i);
        // truthTree->GetEntry(i);

        if (evt != calo_evt)
        {
            std::cerr << "Warning: evt mismatch at entry " << i << ": track evt = " << evt << ", calo evt = " << calo_evt << std::endl;
            continue;
        }

        std::cout << "--- Matching evt : " << evt << "---" << std::endl;
        n_evt_all++;

        int nTrk = (int)track_phi->size();
        int nCalo = (int)calo_phi->size();

        // First pass: recalcu charge and find matched tracks and store info
        for (size_t it = 0; it < nTrk; ++it)
        {
            float trkid = (*track_id)[it];
            int charge_old = (*track_charge)[it];
            float phi_r2mdmin = 0.0;
            float phi_rmax = 0.0;
            float dphi_out_in = 0.0;

            int charge_new = IS_DATA
                                 ? recalcCharge(
                                       trkid, charge_old,
                                       Siclus_trackid, Siclus_x, Siclus_y,
                                       phi_r2mdmin, phi_rmax, dphi_out_in)
                                 : charge_old;
            // int charge_new = charge_old;

            // if (debug)
            // {
            //     cout << "trk " << trkid
            //          << "  phi_r2mdmin=" << phi_r2mdmin
            //          << "  phi_rmax=" << phi_rmax
            //          << "  dphi_clus(out-in)=" << dphi_out_in
            //          << "  charge(old,new)=(" << charge_old << "," << charge_new << ")"
            //          << endl;
            // }

            float pt = (*track_pt)[it];
            float pz = (*track_pz)[it];
            float eta = (*track_eta)[it];
            float phi_trk = (*track_phi)[it];
            float z_trk = (*track_z)[it];
            int nmaps = (*track_nmaps)[it];
            int nintt = (*track_nintt)[it];
            float chi2trk = (*track_chi2ndf)[it];
            // float phi_emc_trk = (*track_phi_emc)[it];
            float phi_emc_trk = atan2((*track_y_emc)[it], (*track_x_emc)[it]);
            float z_emc_trk = (*track_z_emc)[it];

            // Find calo cluster with minimum dphi_emc for this track
            // ---- dphi alignment shift for this track (x = phi0 = phi_trk) ----
            float phi0_for_align = wrapPi(phi_trk);
            float shift = 0.0;
            if (IS_DATA && sp_mu)
            {
                shift = sp_mu->Eval(phi0_for_align); // [rad]
            }

            float min_dr = 1e9;
            int min_ic = -1;
            float min_dphi = 0, min_dphi_emc_raw = 0, min_dphi_emc_corr = 0;
            float min_dz = 0, min_dz_emc = 0, min_dz_corr = 0;

            // calo matching
            for (size_t ic = 0; ic < nCalo; ++ic)
            {
                float phi_calo = (*calo_phi)[ic];
                float z_calo = (*calo_z)[ic];

                float dphi_emc_raw = wrapPi(phi_emc_trk - phi_calo);
                float dphi_emc_raw_cm = dphi_emc_raw * 93.5;

                float dphi_emc_corr = wrapPi(dphi_emc_raw - shift);
                float dphi_emc_corr_cm = dphi_emc_corr * 93.5;

                float dz_emc = z_emc_trk - z_calo; // without correlation
                // === dz recalculating ===
                float dz_corr;
                if (use_dz_correction)
                    dz_corr = dz_emc + (p1 - 1.0) * z_emc_trk; // = p1*(trkr z ) - caloz
                else
                {
                    dz_corr = dz_emc;
                    h_trkrz_vs_caloz->Fill(z_emc_trk, z_calo);
                }
                h_trkrphi_vs_calophi->Fill(phi_emc_trk, phi_calo);
                h_dz_raw->Fill(dz_emc);
                h_dphi_raw->Fill(dphi_emc_raw);
                h_dphi_raw_cm->Fill(dphi_emc_raw_cm);
                h_dphiraw_vs_tkgphi->Fill(phi_emc_trk, dphi_emc_raw);
                h_dphiraw_vs_calophi->Fill(phi_calo, dphi_emc_raw);
                h_dphiraw_vs_phi0->Fill(phi_trk, dphi_emc_raw);
                h_dzraw_vs_caloz->Fill(z_calo, dz_emc);
                h_dzraw_vs_trkrz->Fill(z_emc_trk, dz_emc);
                h_dzraw_vs_z0->Fill(z_trk, dz_emc);

                h_dphi_corr_cm->Fill(dphi_emc_corr_cm);
                h_dphicorr_vs_phi0->Fill(phi_trk, dphi_emc_corr);

                if (debug && t_allcomb)
                {
                    ac.evt = evt;
                    ac.it = (int)it;
                    ac.ic = (int)ic;

                    ac.pt = pt;
                    ac.pz = pz;
                    ac.eta = eta;
                    ac.phi = phi_trk;
                    ac.charge = charge_new;
                    ac.nmaps = nmaps;
                    ac.nintt = nintt;
                    ac.chi2_trk = chi2trk;

                    ac.phi_emc_trk = phi_emc_trk;
                    ac.x_emc_trk = (*track_x_emc)[it];
                    ac.y_emc_trk = (*track_y_emc)[it];
                    ac.z_emc_trk = z_emc_trk;

                    ac.phi_calo = phi_calo;
                    ac.z_calo = z_calo;
                    ac.e_calo = (*calo_energy)[ic];

                    ac.chi2_calo = (IS_DATA && calo_chi2) ? (*calo_chi2)[ic] : -999.0f;

                    ac.dphi_emc_raw = dphi_emc_raw;
                    ac.dphi_emc_corr = dphi_emc_corr;
                    ac.dz_emc = dz_emc;
                    ac.dz_corr = dz_corr;

                    t_allcomb->Fill();
                }

                // if (fabs(dz_corr) > Dz_preCut)
                //     continue;

                // if (phi_calo > M_PI)
                //     phi_calo -= 2 * M_PI;
                // if (phi_calo < -M_PI)
                //     phi_calo += 2 * M_PI;

                // h_trkrphi_vs_calophi->Fill(phi_emc_trk, phi_calo);
                // if (fabs((*track_z_emc)[it] - (*calo_z)[ic]) < 8)
                // if (fabs(1.13 * (*track_z_emc)[it] - (*calo_z)[ic]) < 8)
                //     h_trkrz_vs_caloz->Fill(z_emc_trk, z_calo);

                // h_dz_vs_caloz->Fill(dz_emc, (*calo_z)[ic]); // before calo matching
                // h_dz_vs_trkrz->Fill(dz_emc, (*track_z_emc)[it]);

                // --- dz dphi pre-cuts ---
                // if (fabs(dz_corr) > Dz_preCut)
                //     continue;
                // if (fabs(dphi_emc) > Dphi_preCut)
                //     continue;

                // double score = matching_dphi ? fabs(dphi_emc_raw) : fabs(dz_corr); // single
                double score = sqrt(dz_corr * dz_corr + (dphi_emc_corr_cm * dphi_emc_corr_cm) / Wphi); // data or pythia
                // double score = fabs(dphi_emc); // single

                // if (matching < min_dr)
                if (score < min_dr)
                {
                    min_dr = score;
                    // min_dr = matching;
                    min_ic = ic;
                    min_dphi = wrapPi(phi_trk - phi_calo);
                    min_dphi_emc_raw = dphi_emc_raw;
                    min_dphi_emc_corr = dphi_emc_corr;
                    min_dz = (*track_z)[it] - z_calo;
                    min_dz_emc = dz_emc;
                    min_dz_corr = dz_corr;
                }
            } // calo loop

            if (min_ic < 0)
                continue; // no match

            float x_calo_match = (*calo_x)[min_ic];
            float y_calo_match = (*calo_y)[min_ic];
            float z_calo_match = (*calo_z)[min_ic];
            float phi_calo_match = (*calo_phi)[min_ic];
            float E = (*calo_energy)[min_ic];
            // if (IS_DATA)
            float chi2_calo = (IS_DATA && calo_chi2) ? (*calo_chi2)[min_ic] : -999.0f;

            (void)x_calo_match; // just in case...
            (void)y_calo_match;

            float p = pt * std::cosh(eta);
            float eop = ((p > 0) ? E / p : -999.);

            float zvertex = zvtx;

            if (eop > 0)
                h_eop_all->Fill(eop);

            // h_trkrphi_vs_calophi->Fill(phi_emc_trk, phi_calo_match);

            h_dphi->Fill(min_dphi);
            h_dphi_emc->Fill(min_dphi_emc_corr);
            h_dphi_emc_pt->Fill(pt, min_dphi_emc_corr);

            h_dz->Fill(min_dz);
            h_dz_emc->Fill(min_dz_emc);
            h_eta->Fill(eta);
            h_track_chi2ndf_matched->Fill(chi2trk);

            if (use_dz_correction)
            {
                h_dz_corr->Fill(min_dz_corr);
                h_dzcorr_vs_trkrz->Fill(min_dz_corr, z_emc_trk);
                h_dzcorr_vs_caloz->Fill(min_dz_corr, z_calo_match);
            }

            h_pt_all->Fill(pt);
            h_dz_all->Fill(min_dz_emc);
            h_dphi_all->Fill(min_dphi_emc_corr);
            // h_dphi_raw_cm->Fill(min_dphi_emc * 93.5);
            h_chi2_all->Fill(chi2trk);
            h_nhits_all->Fill(nmaps, nintt);
            h_dphi_vs_tkgphi->Fill(phi_emc_trk, min_dphi_emc_corr);
            h_dphi_vs_calophi->Fill(phi_calo_match, min_dphi_emc_corr);

            h_dz_vs_dphi->Fill(min_dz_emc, min_dphi_emc_corr * 93.5);

            h_dz_vs_caloz->Fill(z_calo_match, min_dz_emc); // after calo matching
            h_dz_vs_trkrz->Fill(z_emc_trk, min_dz_emc);

            h_eta_vs_phi->Fill(min_dphi_emc_corr, eta);
            h_eta_vs_pt->Fill(pt, eta);
            // if (pt >= 0.3)
            h_dphi_vs_pt->Fill(pt, min_dphi_emc_corr);
            h_dz_vs_pt->Fill(pt, min_dz_emc);
            h_eop_vs_pt->Fill(pt, eop);

            if (IS_DATA)
            {
                h_zv_all->Fill(zvertex);
                h_xv_yv->Fill(xvtx, yvtx);
                h_dz0v_all->Fill(z_trk - zvertex);
            }
            // if (chi2 > 5)
            //   continue; // Skip tracks with high chi2/ndf
            // if (std::abs(min_dz_emc) > 4)
            //   continue;
            // h_reco_vs_truth_pt->Fill((*truth_pt)[0], (*track_pt)[it]);

            MatchedTrkCalo mt_temp;
            mt_temp.evt = evt;
            mt_temp.pt = pt;
            mt_temp.z0 = z_trk;
            mt_temp.pz = pz;
            mt_temp.eta = eta;
            mt_temp.phi = phi_trk;
            mt_temp.charge = charge_new;
            mt_temp.nmaps = nmaps;
            mt_temp.nintt = nintt;
            mt_temp.chi2_trk = chi2trk;
            mt_temp.phi_emc_trk = phi_emc_trk;
            mt_temp.z_emc_trk = z_emc_trk;
            mt_temp.phi_calo = phi_calo_match;
            mt_temp.z_calo = z_calo_match;
            mt_temp.e_calo = E;
            // if (IS_DATA)
            mt_temp.chi2_calo = chi2_calo;

            mt_temp.dphi_emc_raw = min_dphi_emc_raw;
            mt_temp.dphi_emc_corr = min_dphi_emc_corr;
            mt_temp.dz_emc_raw = min_dz_emc;
            if (use_dz_correction)
                mt_temp.dz_emc_corr = min_dz_corr;
            else
                mt_temp.dz_emc_corr = min_dz_emc;
            mt_temp.zvertex = zvertex;
            mt_temp.eop = eop;

            PassFlags fl = computeFlags(mt_temp,
                                        IS_DATA,
                                        PtCut,
                                        Dz0VtxCut,
                                        nMapsCut,
                                        nInttCut,
                                        Chi2Cut,
                                        DzCut,
                                        phi_threshold,
                                        EopMinCut,
                                        EopMaxCut,
                                        ZvtxCut);

            // 1D distributions after each individual cut
            if (fl.pass_pt)
                h_pt_cut->Fill(mt_temp.pt);
            if (IS_DATA && fl.pass_z0)
                h_dz0v_cut->Fill(mt_temp.z0 - mt_temp.zvertex);
            if (fl.pass_dz)
                h_dz_cut->Fill(mt_temp.dz_emc_corr); // after z-correction
            if (fl.pass_dphi)
                h_dphi_cut->Fill(mt_temp.dphi_emc_corr);
            if (fl.pass_eop && mt_temp.eop > 0)
                h_eop_cut->Fill(mt_temp.eop);
            if (fl.pass_chi2)
            {
                h_chi2_cut->Fill(mt_temp.chi2_trk);
                h_pt_chi2->Fill(mt_temp.pt);
            }
            if (fl.pass_chi2 && fl.pass_eop) 
                h_pt_chi2_eop->Fill(mt_temp.pt);
            if (fl.pass_hits)
                h_nhits_cut->Fill(mt_temp.nmaps, mt_temp.nintt);
            if (IS_DATA && fl.pass_vtx)
                h_zv_cut->Fill(mt_temp.zvertex);

            matched_this_evt.push_back(mt_temp);
            flags_this_evt.push_back(fl);
        } // track loop

        // ---- (data only) fill t_match here ----
        if (debug && IS_DATA && t_match)
        {
            MatchedTrkCalo br;
            for (size_t i = 0; i < matched_this_evt.size(); ++i)
            {
                br = matched_this_evt[i];
                t_match->Fill();
            }
        }

        // ---- mass reconstruction for THIS event----
        size_t nM = matched_this_evt.size();
        if (nM < 2)
            continue;
        Long64_t n_pair_before = t_pair->GetEntries();
        for (size_t ia = 0; ia < nM; ++ia)
        {
            const MatchedTrkCalo &mt1 = matched_this_evt[ia];
            const PassFlags &f1 = flags_this_evt[ia];

            for (size_t ib = ia + 1; ib < nM; ++ib)
            {
                const MatchedTrkCalo &mt2 = matched_this_evt[ib];
                const PassFlags &f2 = flags_this_evt[ib];

                TLorentzVector v1, v2;
                v1.SetPtEtaPhiM(mt1.pt, mt1.eta, mt1.phi, electron_mass);
                v2.SetPtEtaPhiM(mt2.pt, mt2.eta, mt2.phi, electron_mass);
                float mass = (v1 + v2).M();

                // all charge comb (-+, --, ++)
                if (f1.pass_all && f2.pass_all)
                    h_mass_allc_allcuts->Fill(mass);
                if (f1.pass_middle_eop && f2.pass_middle_eop)
                    h_mass_allc_middle_eop->Fill(mass);

                if (mt1.charge * mt2.charge >= 0)
                {
                    // same charge (--, ++)
                    if (f1.pass_all && f2.pass_all)
                        h_mass_ss_allcuts->Fill(mass);
                    if (f1.pass_middle_eop && f2.pass_middle_eop)
                        h_mass_ss_middle_eop->Fill(mass);
                    continue;
                }

                //-----↓ oppopsite charge signal (-+)-----
                pm.evt = evt;
                pm.ia = (int)ia;
                pm.ib = (int)ib;
                pm.mass = mass;

                pm.q1 = mt1.charge;
                pm.q2 = mt2.charge;
                pm.pt1 = mt1.pt;
                pm.pt2 = mt2.pt;
                pm.eta1 = mt1.eta;
                pm.eta2 = mt2.eta;
                pm.phi1 = mt1.phi;
                pm.phi2 = mt2.phi;
                pm.eop1 = mt1.eop;
                pm.eop2 = mt2.eop;

                pm.pass_all = (f1.pass_all && f2.pass_all);
                pm.pass_no_eop = (f1.pass_no_eop && f2.pass_no_eop);
                pm.pass_middle_eop = (f1.pass_middle_eop && f2.pass_middle_eop);

                t_pair->Fill();

                // ---- below is your original pair-loop filling ----
                h_mass_base->Fill(mass);

                if (mt1.eop > 0)
                    h_eop_vs_mass->Fill(mass, mt1.eop, 0.5);
                if (mt2.eop > 0)
                    h_eop_vs_mass->Fill(mass, mt2.eop, 0.5);

                if (mass > 2.6 && mass < 3.6)
                {
                    if (mt1.eop > 0)
                        h_eop_masscut->Fill(mt1.eop);
                    if (mt2.eop > 0)
                        h_eop_masscut->Fill(mt2.eop);
                }

                if (f1.pass_all && f2.pass_all)
                {
                    h_mass_allcuts->Fill(mass);
                    h_eta_mass->Fill(mt1.eta);
                    h_eta_mass->Fill(mt2.eta);
                    if (mt1.eop > 0)
                        h_eop_vs_mass_pass->Fill(mass, mt1.eop, 0.5);
                    if (mt2.eop > 0)
                        h_eop_vs_mass_pass->Fill(mass, mt2.eop, 0.5);
                }

                if (f1.pass_no_eop && f2.pass_no_eop)
                    h_mass_no_eop->Fill(mass);
                if (f1.pass_middle_eop && f2.pass_middle_eop)
                    h_mass_middle_eop->Fill(mass);

                if (f1.pass_pt && f2.pass_pt)
                    h_mass_only_pt->Fill(mass);
                if (f1.pass_z0 && f2.pass_z0)
                    h_mass_only_z0->Fill(mass);
                if (f1.pass_hits && f2.pass_hits)
                    h_mass_only_hits->Fill(mass);
                if (f1.pass_chi2 && f2.pass_chi2)
                    h_mass_only_chi2->Fill(mass);
                if (f1.pass_dz && f2.pass_dz)
                    h_mass_only_dz->Fill(mass);
                if (f1.pass_dphi && f2.pass_dphi)
                    h_mass_only_dphi->Fill(mass);
                if (f1.pass_vtx && f2.pass_vtx)
                    h_mass_only_zvtx->Fill(mass);
                if (f1.pass_eop && f2.pass_eop)
                    h_mass_only_eop->Fill(mass);
            }
        }
        if (t_pair->GetEntries() > n_pair_before)
            n_evt_hasPair++;
        if (debug && i + 1 >= nloop)
            break;
    }

    // // ========== mass reconstruction ==========
    // for (size_t ie = 0; ie < matched_evt.size(); ++ie)
    // {
    //     const std::vector<MatchedTrkCalo> &v = matched_evt[ie];
    //     const std::vector<PassFlags> &flv = flags_evt[ie];

    //     size_t nM = v.size();
    //     if (nM < 2)
    //         continue;
    //     for (size_t ia = 0; ia < nM; ++ia)
    //     {
    //         const MatchedTrkCalo &mt1 = v[ia];
    //         const PassFlags &f1 = flv[ia];

    //         for (size_t ib = ia + 1; ib < nM; ++ib)
    //         {
    //             const MatchedTrkCalo &mt2 = v[ib];
    //             const PassFlags &f2 = flv[ib];

    //             // Opposite-sign pairs only
    //             if (mt1.charge * mt2.charge >= 0)
    //                 continue;

    //             TLorentzVector v1, v2;
    //             v1.SetPtEtaPhiM(mt1.pt, mt1.eta, mt1.phi, electron_mass);
    //             v2.SetPtEtaPhiM(mt2.pt, mt2.eta, mt2.phi, electron_mass);
    //             float mass = (v1 + v2).M();

    //             h_mass_base->Fill(mass);

    //             if (mt1.eop > 0)
    //                 h_eop_vs_mass->Fill(mass, mt1.eop, 0.5);
    //             if (mt2.eop > 0)
    //                 h_eop_vs_mass->Fill(mass, mt2.eop, 0.5);

    //             if (mass > 2.6 && mass < 3.6)
    //             {
    //                 if (mt1.eop > 0)
    //                     h_eop_masscut->Fill(mt1.eop);
    //                 if (mt2.eop > 0)
    //                     h_eop_masscut->Fill(mt2.eop);
    //             }

    //             // Combined selections

    //             if (f1.pass_all && f2.pass_all)
    //             {
    //                 h_mass_allcuts->Fill(mass);
    //                 h_eta_mass->Fill(mt1.eta);
    //                 h_eta_mass->Fill(mt2.eta);
    //                 if (mt1.eop > 0)
    //                     h_eop_vs_mass_pass->Fill(mass, mt1.eop, 0.5);
    //                 if (mt2.eop > 0)
    //                     h_eop_vs_mass_pass->Fill(mass, mt2.eop, 0.5);
    //             }

    //             if (f1.pass_no_eop && f2.pass_no_eop)
    //             {
    //                 h_mass_no_eop->Fill(mass);
    //             }

    //             if (f1.pass_middle_eop && f2.pass_middle_eop)
    //             {
    //                 h_mass_middle_eop->Fill(mass);
    //             }

    //             if (f1.pass_pt && f2.pass_pt)
    //                 h_mass_only_pt->Fill(mass);
    //             if (f1.pass_hits && f2.pass_hits)
    //                 h_mass_only_hits->Fill(mass);
    //             if (f1.pass_chi2 && f2.pass_chi2)
    //                 h_mass_only_chi2->Fill(mass);
    //             if (f1.pass_dz && f2.pass_dz)
    //                 h_mass_only_dz->Fill(mass);
    //             if (f1.pass_dphi && f2.pass_dphi)
    //                 h_mass_only_dphi->Fill(mass);
    //             if (f1.pass_vtx && f2.pass_vtx)
    //                 h_mass_only_zvtx->Fill(mass);
    //             if (f1.pass_eop && f2.pass_eop)
    //                 h_mass_only_eop->Fill(mass);
    //         }
    //     }
    // } // end of pair loop

    // ========== store info to output file ==========

    std::cout << "Saving histogram and ttree to root file" << std::endl;

    if (IS_DATA && debug && t_match)
        t_match->Write();

    if (debug && t_allcomb)
        t_allcomb->Write();

    t_pair->Write();

    h_dphi->Write();
    h_dphi_emc->Write();
    h_dphi_emc_pt->Write();
    h_dz->Write();
    h_dz_emc->Write();
    h_mass->Write();
    h_track_chi2ndf_matched->Write();
    h_reco_vs_truth_pt->Write();
    h_mass_base->Write();
    h_mass_only_pt->Write();
    h_mass_only_z0->Write();
    h_mass_only_hits->Write();
    h_mass_only_chi2->Write();
    h_mass_only_dz->Write();
    h_eop_vs_mass->Write();
    h_eop_vs_mass_pass->Write();
    h_mass_only_dphi->Write();
    h_mass_only_zvtx->Write();
    h_mass_only_eop->Write();
    h_mass_allcuts->Write();
    h_mass_no_eop->Write();
    h_mass_middle_eop->Write();
    h_mass_ss_allcuts->Write();
    h_mass_allc_middle_eop->Write();
    h_mass_allc_allcuts->Write();
    h_mass_ss_middle_eop->Write();
    h_pt_all->Write();
    h_pt_cut->Write();
    h_pt_chi2->Write();
    h_pt_chi2_eop->Write();
    h_dz0v_all->Write();
    h_dz0v_cut->Write();
    h_dz_all->Write();
    h_dz_corr->Write();
    h_dz_cut->Write();
    h_dphi_all->Write();
    h_dphi_cut->Write();
    h_dphi_raw_cm->Write();
    h_eop_all->Write();
    h_eop_cut->Write();
    h_eop_masscut->Write();
    h_chi2_all->Write();
    h_chi2_cut->Write();
    h_zv_all->Write();
    h_zv_cut->Write();
    h_nhits_all->Write();
    h_nhits_cut->Write();
    h_eta->Write();
    h_eta_mass->Write();
    h_eta_vs_phi->Write();
    h_eta_vs_pt->Write();
    h_dphi_vs_pt->Write();
    h_dz_vs_pt->Write();
    h_eop_vs_pt->Write();
    h_dphi_vs_tkgphi->Write();
    h_dphi_vs_calophi->Write();
    h_trkrphi_vs_calophi->Write();
    h_dz_vs_trkrz->Write();
    h_dz_vs_caloz->Write();
    h_dzcorr_vs_caloz->Write();
    h_dzcorr_vs_trkrz->Write();
    h_dz_vs_dphi->Write();
    h_dz_raw->Write();
    h_dphi_raw->Write();
    h_dphiraw_vs_tkgphi->Write();
    h_dphiraw_vs_calophi->Write();
    h_dphiraw_vs_phi0->Write();
    h_dzraw_vs_caloz->Write();
    h_dzraw_vs_trkrz->Write();
    h_dzraw_vs_z0->Write();
    h_dphi_corr_cm->Write();
    h_dphicorr_vs_phi0->Write();

    std::cout << "Δφ histogram saved to : " << fname_out << std::endl;

    h_mass_allcuts->Draw();

    c->cd(2);
    h_mass_base->Draw();
    // h_mass_allcuts->SetLineColor(TColor::GetColor(50, 205, 50));
    // h_mass_allcuts->Draw("same");
    c->Print(fname_out_pdf);

    c->cd(1);
    // h_pt_cut->Draw();
    h_pt_all->Draw();
    h_pt_cut->SetLineColor(TColor::GetColor(50, 205, 50));
    h_pt_cut->Draw("same");
    c->cd(2);
    gPad->SetLogy(1);
    h_mass_base->Draw();
    h_mass_only_pt->SetLineColor(TColor::GetColor(50, 205, 50));
    h_mass_only_pt->Draw("same");
    c->Print(fname_out_pdf);

    c->cd(1);
    h_dz0v_all->Draw();
    h_dz0v_cut->SetLineColor(TColor::GetColor(50, 205, 50));
    h_dz0v_cut->Draw("same");
    c->cd(2);
    gPad->SetLogy(1);
    h_mass_base->Draw();
    h_mass_only_z0->SetLineColor(TColor::GetColor(50, 205, 50));
    h_mass_only_z0->Draw("same");
    c->Print(fname_out_pdf);

    c->cd(1);
    h_nhits_cut->Draw();
    c->cd(2);
    gPad->SetLogy(1);
    h_mass_base->Draw();
    h_mass_only_hits->SetLineColor(TColor::GetColor(50, 205, 50));
    h_mass_only_hits->Draw("same");
    c->Print(fname_out_pdf);

    c->cd(1);
    // h_chi2_cut->Draw();
    h_chi2_all->Draw();
    h_chi2_cut->SetLineColor(TColor::GetColor(50, 205, 50));
    h_chi2_cut->Draw("same");
    c->cd(2);
    gPad->SetLogy(1);
    h_mass_base->Draw();
    h_mass_only_chi2->SetLineColor(TColor::GetColor(50, 205, 50));
    h_mass_only_chi2->Draw("same");
    c->Print(fname_out_pdf);

    c->cd(1);
    h_dz_all->Draw();
    h_dz_cut->SetLineColor(TColor::GetColor(50, 205, 50));
    h_dz_cut->Draw("same");

    c->cd(2);
    gPad->SetLogy(1);
    h_mass_base->Draw();
    h_mass_only_dz->SetLineColor(TColor::GetColor(50, 205, 50));
    h_mass_only_dz->Draw("same");
    c->Print(fname_out_pdf);

    c->cd(1);
    h_dphi_all->Draw();
    h_dphi_cut->SetLineColor(TColor::GetColor(50, 205, 50));
    h_dphi_cut->Draw("same");
    c->cd(2);
    h_mass_base->Draw();
    h_mass_only_dphi->SetLineColor(TColor::GetColor(50, 205, 50));
    h_mass_only_dphi->Draw("same");
    c->Print(fname_out_pdf);

    c->cd(1);
    // h_dphi_cut->Draw();
    h_zv_all->Draw();
    h_zv_cut->SetLineColor(TColor::GetColor(50, 205, 50));
    h_zv_cut->Draw("same");
    c->cd(2);
    h_mass_base->Draw();
    h_mass_only_zvtx->SetLineColor(TColor::GetColor(50, 205, 50));
    h_mass_only_zvtx->Draw("same");
    c->Print(fname_out_pdf);

    c->cd(1);
    // h_eop_cut->Draw();
    h_eop_all->Draw();
    h_eop_cut->SetLineColor(TColor::GetColor(50, 205, 50));
    h_eop_cut->Draw("same");
    c->cd(2);
    gPad->SetLogy(1);
    h_mass_base->Draw();
    h_mass_only_eop->SetLineColor(TColor::GetColor(50, 205, 50));
    h_mass_only_eop->Draw("same");
    c->Print(fname_out_pdf);

    c->cd(1);
    h_eta_vs_phi->Draw();
    c->cd(2);
    gPad->SetLogy(0);
    h_eta_vs_pt->Draw();
    c->Print(fname_out_pdf);

    c->cd(1);
    h_dphi_vs_pt->Draw();
    gPad->SetLogz();
    c->cd(2);
    gPad->SetLogy(0);
    h_eop_vs_pt->Draw();
    c->Print(fname_out_pdf);

    c->cd(1);
    gPad->SetLogz(0);
    h_dphi_vs_pt->Draw();
    c->cd(2);
    h_dz_vs_pt->Draw();
    c->Print(fname_out_pdf);

    c->cd(1);
    h_dphi_vs_tkgphi->Draw();
    c->cd(2);
    h_dphi_vs_calophi->Draw();
    c->Print(fname_out_pdf);

    c->cd(1);
    h_dphiraw_vs_tkgphi->Draw();
    c->cd(2);
    h_dphiraw_vs_calophi->Draw();
    c->Print(fname_out_pdf);

    c->cd(1);
    h_dphiraw_vs_phi0->Draw();
    c->cd(2);
    h_dzraw_vs_z0->Draw();
    c->Print(fname_out_pdf);

    c->cd(1);
    h_trkrphi_vs_calophi->Draw();
    c->cd(2);
    h_trkrz_vs_caloz->Draw();
    if (use_dz_correction)
    {
        pf_trkrz_vs_caloz->Draw("same");
        fit_trkrz_vs_caloz->Draw("same");
        gPad->Update();

        TPaveStats *st = (TPaveStats *)h_trkrz_vs_caloz->GetListOfFunctions()->FindObject("stats");
        if (st)
        {
            st->SetX1NDC(0.60);
            st->SetX2NDC(0.90);
            st->SetY1NDC(0.10);
            st->SetY2NDC(0.30);
        }

        gPad->Modified();
        gPad->Update();
    }
    c->Print(fname_out_pdf);

    c->cd(1);
    h_dzraw_vs_caloz->Draw();
    c->cd(2);
    h_dzraw_vs_trkrz->Draw();
    c->Print(fname_out_pdf);

    c->cd(1);
    h_dz_raw->Draw();
    c->cd(2);
    h_dphi_raw->Draw();
    c->Print(fname_out_pdf);

    c->cd(1);
    h_dz_vs_dphi->Draw();
    c->cd(2);
    h_dphi_raw_cm->Draw();
    c->Print(fname_out_pdf);

    c->cd(1);
    h_dphi_corr_cm->Draw();
    c->cd(2);
    h_dphicorr_vs_phi0->Draw();
    c->Print(fname_out_pdf);

    c->cd(1);
    gPad->SetLeftMargin(0.12);
    h_dz_vs_trkrz->Draw();
    c->cd(2);
    h_dz_vs_caloz->Draw();
    c->Print(fname_out_pdf);

    c->cd(1);
    h_dz_corr->Draw();
    h_dz_cut->SetLineColor(TColor::GetColor(50, 205, 50));
    h_dz_cut->Draw("same");
    c->cd(2);
    h_dz_all->Draw();
    c->Print(fname_out_pdf);

    c->cd(1);
    h_dzcorr_vs_trkrz->Draw();
    c->cd(2);
    h_dzcorr_vs_caloz->Draw();
    c->Print(fname_out_pdf);

    c->cd(1);
    h_mass_allcuts->Draw();
    c->cd(2);
    h_mass_ss_allcuts->Draw();
    c->Print(fname_out_pdf);

    c->cd(1);
    h_mass_middle_eop->Draw();
    c->cd(2);
    h_mass_ss_middle_eop->Draw();
    c->Print(fname_out_pdf);

    c->cd(1);
    h_mass_allc_middle_eop->Draw();
    c->cd(2);
    h_mass_allc_allcuts->Draw();
    c->Print(fname_out_pdf);

    c->cd(1);
    h_mass_middle_eop->Draw();
    h_mass_ss_middle_eop->SetLineColor(TColor::GetColor(50, 205, 50));
    h_mass_ss_middle_eop->Draw("SAME");
    c->cd(2);
    h_mass_allcuts->Draw();
    h_mass_ss_allcuts->SetLineColor(TColor::GetColor(50, 205, 50));
    h_mass_ss_allcuts->Draw("SAME");
    c->Print(fname_out_pdf);

    c->cd(1);
    h_mass_middle_eop->Draw();
    h_mass_allc_middle_eop->SetLineColor(TColor::GetColor(255, 105, 180));
    h_mass_allc_middle_eop->Draw("SAME");
    c->cd(2);
    h_mass_allcuts->Draw();
    h_mass_allc_allcuts->SetLineColor(TColor::GetColor(255, 105, 180));
    h_mass_allc_allcuts->Draw("SAME");
    c->Print(fname_out_pdf);

    c->cd(1);
    h_mass_no_eop->SetLineColor(TColor::GetColor(50, 205, 50));
    h_mass_no_eop->Draw();
    h_mass_allcuts->Draw("SAME");
    h_mass_middle_eop->SetLineColor(TColor::GetColor(255, 105, 180));
    h_mass_middle_eop->Draw("SAME");

    c->cd(2);
    gStyle->SetOptStat(0);
    h_mass_no_eop->SetStats(0);
    h_mass_allcuts->SetStats(0);
    h_mass_middle_eop->SetStats(0);
    gPad->SetLogy(1);
    h_mass_no_eop->GetYaxis()->SetRangeUser(1e-1, 1000);
    h_mass_no_eop->SetLineColor(TColor::GetColor(50, 205, 50));
    h_mass_no_eop->Draw();
    h_mass_allcuts->Draw("SAME");
    h_mass_middle_eop->SetLineColor(TColor::GetColor(255, 105, 180));
    h_mass_middle_eop->Draw("SAME");
    h_mass_no_eop->SetTitle("Mass distribution depending on the E/p cut");
    c->Print(fname_out_pdf);

    c->Print(fname_out_pdf + "]");

    h_trkrz_vs_caloz->Write();

    TParameter<Long64_t>("n_evt_all", n_evt_all).Write();
    TParameter<Long64_t>("n_evt_hasPair", n_evt_hasPair).Write();
    outFile->Close();

    std::cout << "histograms are drawn to : " << fname_out_pdf << std::endl;
    f->Close();
    fAlign->Close();
}