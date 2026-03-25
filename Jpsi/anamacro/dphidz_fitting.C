// dphidz_fitting.C  (coarse fit page + final fit page)
// gaus(0)+pol0(3)
#include <iostream>
#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TStyle.h>
#include <TPaveText.h>

using namespace std;

void dphidz_fitting(
    const char *fin = "result/sim/pythia_mass_reco_100m_merged.root",
    // const char *fin = "result/sim/mass_reco_200k.root",
    const char *fout = "result/sim/fit_dphi_dz_pyhita.root",
    const char *pdf = "result/sim/fit_dphi_dz_pythia.pdf")
{
    gStyle->SetOptStat(1111);
    gStyle->SetOptFit(1111);
    gStyle->SetOptTitle(1);
    TFile *f = TFile::Open(fin, "READ");
    if (!f || f->IsZombie())
    {
        cout << "[ERROR] cannot open: " << fin << endl;
        return;
    }

    TH1 *h_dphi_raw = (TH1 *)f->Get("h_dphi_raw_cm");
    TH1 *h_dphi_corr = (TH1 *)f->Get("h_dphi_corr_cm");
    TH1 *h_dz_raw = (TH1 *)f->Get("h_dz_raw");
    TH1 *h_dz_corr = (TH1 *)f->Get("h_dz_corr");

    if (!h_dphi_raw)
        cout << "[WARN] missing: h_dphi_raw_cm" << endl;
    if (!h_dphi_corr)
        cout << "[WARN] missing: h_dphi_corr_cm" << endl;
    if (!h_dz_raw)
        cout << "[WARN] missing: h_dz_raw" << endl;
    if (!h_dz_corr)
        cout << "[WARN] missing: h_dz_corr" << endl;

    TFile *fo = new TFile(fout, "RECREATE");

    TCanvas *c = new TCanvas("c", "fit", 900, 800);
    c->Print(Form("%s(", pdf));

    // --- coarse fit windows
    // const double dphi_coarse_lo = -10.0;
    // const double dphi_coarse_hi = +10.0;
    // const double dz_coarse_lo = -20.0;
    // const double dz_coarse_hi = +20.0;

    // --- final fit window = mu1 ± nsig_final * sigma1
    const double nsig_final = 3.;

    // --- initial sigma guesses
    const double sig0_dphi = 1.0; // cm
    const double sig0_dz = 2.0;   // cm

    // =========================
    // helper: 右上に表示（hist名 + stage）
    // =========================
    auto draw_label = [&](TH1 *h, const char *stage)
    {
        TPaveText *p = new TPaveText(0.60, 0.90, 0.88, 0.97, "NDC");
        p->SetFillStyle(0);
        p->SetBorderSize(0);
        p->AddText(Form("%s", h->GetName()));
        p->AddText(Form("%s", stage));
        p->Draw();
    };

    // ============================================================
    // (1) dphi raw
    // ============================================================
    if (h_dphi_raw)
    {
        c->Clear();
        h_dphi_raw->SetStats(1); 
        h_dphi_raw->Draw();

        double sd_coarse = h_dphi_raw->GetStdDev();
        if (sd_coarse <= 0)
            sd_coarse = 1.0;
        TF1 *f1 = new TF1("f_dphi_raw", "gaus(0)+pol0(3)", -sd_coarse, +sd_coarse);

        int ib0 = h_dphi_raw->GetXaxis()->FindBin(0.0);
        int ibL = h_dphi_raw->GetXaxis()->FindBin(-sd_coarse);
        int ibH = h_dphi_raw->GetXaxis()->FindBin(+sd_coarse);
        double A0 = h_dphi_raw->GetBinContent(ib0);
        double BG0 = 0.5 * (h_dphi_raw->GetBinContent(ibL) + h_dphi_raw->GetBinContent(ibH));
        f1->SetParameters(A0 - BG0, 0.0, sig0_dphi, BG0);

        // --- coarse fit
        h_dphi_raw->Fit(f1, "QR");
        c->Modified();
        c->Update();
        c->Print(pdf);

        // 保存（coarse）
        fo->cd();
        TF1 *fcoarse = (TF1 *)f1->Clone("fit_dphi_raw_coarse");
        fcoarse->Write();

        // --- final fit range from coarse sigma
        double mu1 = f1->GetParameter(1);
        double sig1 = fabs(f1->GetParameter(2));
        if (sig1 < 1e-3 || sig1 > 1e2)
            sig1 = sig0_dphi;
        f1->SetRange(mu1 - nsig_final * sig1, mu1 + nsig_final * sig1);

        // --- final fit
        h_dphi_raw->Fit(f1, "QR");
        c->Modified();
        c->Update();
        c->Print(pdf);

        // 保存（final）＋ hist
        fo->cd();
        h_dphi_raw->Write(h_dphi_raw->GetName());
        TF1 *ffinal = (TF1 *)f1->Clone("fit_dphi_raw_final");
        ffinal->Write();
    }

    // ============================================================
    // (2) dphi corr
    // ============================================================
    if (h_dphi_corr)
    {
        c->Clear();
        h_dphi_corr->SetStats();
        h_dphi_corr->Draw();

        double sd_coarse = h_dphi_corr->GetStdDev();
        if (sd_coarse <= 0)
            sd_coarse = 1.0;
        TF1 *f1 = new TF1("f_dphi_corr", "gaus(0)+pol0(3)", -sd_coarse, +sd_coarse);

        int ib0 = h_dphi_corr->GetXaxis()->FindBin(0.0);
        int ibL = h_dphi_corr->GetXaxis()->FindBin(-sd_coarse);
        int ibH = h_dphi_corr->GetXaxis()->FindBin(+sd_coarse);

        double A0 = h_dphi_corr->GetBinContent(ib0);
        double BG0 = 0.5 * (h_dphi_corr->GetBinContent(ibL) + h_dphi_corr->GetBinContent(ibH));
        f1->SetParameters(A0 - BG0, 0.0, sig0_dphi, BG0);

        // coarse
        h_dphi_corr->Fit(f1, "QR");
        c->Modified();
        c->Update();
        c->Print(pdf);

        fo->cd();
        TF1 *fcoarse = (TF1 *)f1->Clone("fit_dphi_corr_coarse");
        fcoarse->Write();

        double mu1 = f1->GetParameter(1);
        double sig1 = fabs(f1->GetParameter(2));
        if (sig1 < 1e-3 || sig1 > 1e2)
            sig1 = sig0_dphi;
        f1->SetRange(mu1 - nsig_final * sig1, mu1 + nsig_final * sig1);

        // final
        h_dphi_corr->Fit(f1, "QR");
        c->Modified();
        c->Update();
        c->Print(pdf);

        fo->cd();
        h_dphi_corr->Write(h_dphi_corr->GetName());
        TF1 *ffinal = (TF1 *)f1->Clone("fit_dphi_corr_final");
        ffinal->Write();
    }

    // ============================================================
    // (3) dz raw
    // ============================================================
    if (h_dz_raw)
    {
        c->Clear();
        h_dz_raw->SetStats(1);
        h_dz_raw->Draw();

        double sd_coarse = h_dz_raw->GetStdDev();
        if (sd_coarse <= 0)
            sd_coarse = 1.0;

        TF1 *f1 = new TF1("f_dz_raw", "gaus(0)+pol0(3)", -sd_coarse, +sd_coarse);

        int ib0 = h_dz_raw->GetXaxis()->FindBin(0.0);
        int ibL = h_dz_raw->GetXaxis()->FindBin(-sd_coarse);
        int ibH = h_dz_raw->GetXaxis()->FindBin(+sd_coarse);
        double A0 = h_dz_raw->GetBinContent(ib0);
        double BG0 = 0.5 * (h_dz_raw->GetBinContent(ibL) + h_dz_raw->GetBinContent(ibH));
        f1->SetParameters(A0 - BG0, 0.0, sig0_dz, BG0);

        // coarse
        h_dz_raw->Fit(f1, "QR");
        c->Modified();
        c->Update();
        c->Print(pdf);

        fo->cd();
        TF1 *fcoarse = (TF1 *)f1->Clone("fit_dz_raw_coarse");
        fcoarse->Write();

        double mu1 = f1->GetParameter(1);
        double sig1 = fabs(f1->GetParameter(2));
        if (sig1 < 1e-3 || sig1 > 1e3)
            sig1 = sig0_dz;
        f1->SetRange(mu1 - nsig_final * sig1, mu1 + nsig_final * sig1);

        // final
        h_dz_raw->Fit(f1, "QR");
        c->Modified();
        c->Update();
        c->Print(pdf);

        fo->cd();
        h_dz_raw->Write(h_dz_raw->GetName());
        TF1 *ffinal = (TF1 *)f1->Clone("fit_dz_raw_final");
        ffinal->Write();
    }

    // ============================================================
    // (4) dz corr
    // ============================================================
    if (h_dz_corr)
    {
        c->Clear();
        h_dz_corr->SetStats(1);;
        h_dz_corr->Draw();

        double sd_coarse = h_dz_corr->GetStdDev();
        if (sd_coarse <= 0)
            sd_coarse = 1.0;

        TF1 *f1 = new TF1("f_dz_corr", "gaus(0)+pol0(3)", -sd_coarse, +sd_coarse);

        int ib0 = h_dz_corr->GetXaxis()->FindBin(0.0);
        int ibL = h_dz_corr->GetXaxis()->FindBin(-sd_coarse);
        int ibH = h_dz_corr->GetXaxis()->FindBin(+sd_coarse);
        double A0 = h_dz_corr->GetBinContent(ib0);
        double BG0 = 0.5 * (h_dz_corr->GetBinContent(ibL) + h_dz_corr->GetBinContent(ibH));
        f1->SetParameters(A0 - BG0, 0.0, sig0_dz, BG0);

        // coarse
        h_dz_corr->Fit(f1, "QR");
        c->Modified();
        c->Update();
        c->Print(pdf);

        fo->cd();
        TF1 *fcoarse = (TF1 *)f1->Clone("fit_dz_corr_coarse");
        fcoarse->Write();

        double mu1 = f1->GetParameter(1);
        double sig1 = fabs(f1->GetParameter(2));
        if (sig1 < 1e-3 || sig1 > 1e3)
            sig1 = sig0_dz;
        f1->SetRange(mu1 - nsig_final * sig1, mu1 + nsig_final * sig1);

        // final
        h_dz_corr->Fit(f1, "QR");
        c->Modified();
        c->Update();
        c->Print(pdf);

        fo->cd();
        h_dz_corr->Write(h_dz_corr->GetName());
        TF1 *ffinal = (TF1 *)f1->Clone("fit_dz_corr_final");
        ffinal->Write();
    }

    // ============================================================
    // 追加：見た目を整えた最終fitページ（2ページ）
    //   dphi: corr final
    //   dz  : raw  final
    // ============================================================
    {
        // styleを元に戻せるように退避
        const int optStat_bak = gStyle->GetOptStat();
        const int optFit_bak = gStyle->GetOptFit();
        const int optTitle_bak = gStyle->GetOptTitle();

        gPad->SetLeftMargin(0.15);

        // stats/fit boxを消す
        // gStyle->SetOptStat(0);
        // gStyle->SetOptFit(0);
        // gStyle->SetOptTitle(1);

        // --- Δφ (corr) clean page
        if (h_dphi_corr)
        {
            c->Clear();
            h_dphi_corr->SetStats(0);
            h_dphi_corr->SetTitle("#Delta#phi");
            h_dphi_corr->Draw();

            TF1 *f_final = (TF1 *)gROOT->FindObject("f_dphi_corr"); // このマクロ内で作ったTF1名
            if (f_final)
                f_final->Draw("same");

            c->Modified();
            c->Update();
            c->Print(pdf);
        }

        // --- Δz (raw) clean page
        if (h_dz_raw)
        {
            c->Clear();
            h_dz_raw->SetStats(0);
            h_dz_raw->SetTitle("#Deltaz");
            h_dz_raw->Draw();

            TF1 *f_final = (TF1 *)gROOT->FindObject("f_dz_raw"); // このマクロ内で作ったTF1名
            if (f_final)
                f_final->Draw("same");

            c->Modified();
            c->Update();
            c->Print(pdf);
        }

        // styleを復元
        gStyle->SetOptStat(optStat_bak);
        gStyle->SetOptFit(optFit_bak);
        gStyle->SetOptTitle(optTitle_bak);
    }

    c->Print(Form("%s)", pdf));

    fo->Close();
    f->Close();

    cout << "[DONE] wrote: " << fout << endl;
    cout << "[DONE] saved: " << pdf << endl;
}
