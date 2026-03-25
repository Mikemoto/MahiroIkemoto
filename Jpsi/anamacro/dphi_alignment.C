#include <TFile.h>
#include <TH2.h>
#include <TH1.h>
#include <TF1.h>
#include <TGraph.h>
#include <TSpline.h>
#include <TVector2.h>
#include <TCanvas.h>
#include <TDirectory.h>
#include <TString.h>
#include <TStyle.h>

void dphi_alignment()
{
    TFile *fin = TFile::Open("result/data/53879_mass_reco_zenhan_debug_500000evt.root");
    TH2F *h2 = (TH2F *)fin->Get("h_dphiraw_vs_phi0");

    const int nxb = h2->GetNbinsX();
    const int merge = 10;
    const double fit_lo = -0.2;
    const double fit_hi = 0.2;
    const int min_entries = 1;
    const double R_emc = 93.5; // cm
    TFile *fout = TFile::Open("result/data/alignment/dphi_phi0_alighnment.root", "RECREATE");
    TDirectory *dSlices = fout->mkdir("slices");

    TH2F *h2raw_cm = new TH2F("h_dphi_raw_cm_vs_phi0",
                              "raw d#phi vs #phi_{0};#phi_{0} [rad];d#phi [cm]",
                              h2->GetNbinsX(),
                              h2->GetXaxis()->GetXmin(),
                              h2->GetXaxis()->GetXmax(),
                              60, -60, 60);

    TH2F *h2corr = (TH2F *)h2->Clone("h_dphi_corr_vs_phi0");
    h2corr->Reset();
    TH2F *h2corr_cm = new TH2F("h_dphi_corr_cm_vs_phi0",
                               "d#phi corrected vs #phi_{0};#phi_{0} [rad];d#phi_{corr} [cm]",
                               h2->GetNbinsX(),
                               h2->GetXaxis()->GetXmin(),
                               h2->GetXaxis()->GetXmax(),
                               60, -60, 60);

    TCanvas *c = new TCanvas("c", "c", 900, 700);
    gStyle->SetOptStat(1110);
    gStyle->SetOptFit(1111);
    const TString pdfname = TString::Format("result/data/alignment/dphi_slices_by_phi0_%dbin.pdf", merge);
    c->Print(pdfname + "[");

    // TF1 *fg = new TF1("fg", "gaus", -0.2, 0.2);
    TF1 *fg = new TF1("fg", "gaus(0)+pol0(3)", -0.2, 0.2);
    TGraph *g_mu = new TGraph();
    int ip = 0;

    for (int ix = 1; ix <= nxb; ix += merge)
    {
        int ix1 = ix;
        int ix2 = std::min(ix + merge - 1, nxb);

        cout << "start bin : " << ix1 << "  end bin : " << ix2 << endl;

        TH1D *hy = h2->ProjectionY(Form("dphi_phi0bin_%04d_%04d", ix1, ix2), ix1, ix2);
        if (hy->GetEntries() < min_entries)
        {
            delete hy;
            continue;
        }

        double phi0_lo = h2->GetXaxis()->GetBinLowEdge(ix1);
        double phi0_hi = h2->GetXaxis()->GetBinUpEdge(ix2);
        double phi0_center = 0.5 * (phi0_lo + phi0_hi);

        hy->SetTitle(Form("d#phi slice: #phi_{0} in [%.3f, %.3f] rad;d#phi [rad];Entries", phi0_lo, phi0_hi));

        // 初期値
        int ib0 = hy->GetXaxis()->FindBin(0.0);
        fg->SetParameters(hy->GetBinContent(ib0), 0.0, 0.05);

        hy->Fit(fg, "QR"); // Quiet,  Use Range
        double mu = fg->GetParameter(1);

        g_mu->SetPoint(ip++, phi0_center, mu);

        c->cd();
        hy->Draw("E");
        fg->Draw("same");

        c->Print(pdfname);

        fout->cd();
        dSlices->cd();
        hy->Write();

        delete hy;
    }
    c->cd();
    h2->Draw("COLZ");
    c->Print(pdfname);

    TSpline3 *sp_mu = new TSpline3("sp_mu", g_mu);

    fout->cd();

    for (int ix = 1; ix <= nxb; ix++)
    {
        double phi0 = h2->GetXaxis()->GetBinCenter(ix);
        double shift = sp_mu->Eval(phi0);

        for (int iy = 1; iy <= h2->GetNbinsY(); iy++)
        {
            double w = h2->GetBinContent(ix, iy);
            if (w <= 0)
                continue;

            double dphi = h2->GetYaxis()->GetBinCenter(iy);
            double dphi_corr = TVector2::Phi_mpi_pi(dphi - shift);
            h2raw_cm->Fill(phi0, R_emc * dphi, w);
            h2corr->Fill(phi0, dphi_corr, w);
            h2corr_cm->Fill(phi0, R_emc * dphi_corr, w);
        }
    }
    c->cd();
    TH1D *h_dphi_corr_all = h2corr->ProjectionY(Form("h_dphi_corr_all_%d", (int)merge));
    h_dphi_corr_all->SetDirectory(0);
    h_dphi_corr_all->SetTitle("corrected d#phi (all #phi_{0});d#phi [rad];Entries");
    h_dphi_corr_all->Draw("HIST");
    c->Print(pdfname);

    c->cd();
    h2corr->Draw("COLZ");
    c->Print(pdfname);

    c->cd();
    h2raw_cm->Draw("COLZ");
    c->Print(pdfname);

    c->cd();
    h2corr_cm->Draw("COLZ");
    c->Print(pdfname);

    c->cd();
    TH1D *h_dphi_corr_all_cm = h2corr_cm->ProjectionY(Form("h_dphi_corr_all_cm_%d", (int)merge));
    h_dphi_corr_all_cm->SetDirectory(0);
    h_dphi_corr_all_cm->SetTitle("corrected d#phi (all #phi_{0}) [cm];d#phi [cm];Entries");
    h_dphi_corr_all_cm->Draw("HIST");
    c->Print(pdfname);

    fout->cd();
    h_dphi_corr_all_cm->Write();
    delete h_dphi_corr_all_cm;

    c->Print(pdfname + "]");

    fout->cd();
    h2->Write();
    h2corr->Write();
    h2corr_cm->Write();
    h_dphi_corr_all->Write();
    g_mu->Write("g_mu_vs_phi0");
    sp_mu->Write("spline_mu_vs_phi0");
    fout->Close();

    std::cout << "Done." << std::endl;
    std::cout << "  PDF : " << pdfname << std::endl;
}