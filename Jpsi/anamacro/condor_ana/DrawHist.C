#include <TFile.h>
#include <TString.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TF1.h>
#include <TPad.h>
#include <TPaveStats.h>
#include <iostream>

void DrawHist(int run_nu = 1, bool use_dz_correction = true)
{
    TString outroot;
    if (run_nu == 1)
    {
        outroot = "/sphenix/tg/tg01/commissioning/INTT/work/mahiro/SIliconCalo/run24pp/macro/result/sim/pythia_mass_reco_100m_merged.root";
        // outroot = "/sphenix/tg/tg01/commissioning/INTT/work/mahiro/SIliconCalo/run24pp/macro/result/sim/mass_reco_200k.root";
    }
    else
    {
        outroot = Form("/sphenix/tg/tg01/commissioning/INTT/work/mahiro/SIliconCalo/run24pp/macro/result/data/%d_mass_reco_merged.root", run_nu);
    }
    // outroot = "/sphenix/tg/tg01/commissioning/INTT/work/mahiro/SIliconCalo/run24pp/macro/result/data/all/mass_reco_run24pp_merged.root";

    TFile *f = TFile::Open(outroot);
    if (!f || f->IsZombie())
    {
        std::cerr << "cannot open : " << outroot << std::endl;
        return;
    }

    TString fname_out_pdf = outroot;
    fname_out_pdf.ReplaceAll(".root", "_byebyecuthist.pdf");

    TCanvas *c = new TCanvas("c", "c", 1200, 600);
    c->Divide(2, 1);
    c->SetLeftMargin(0.6);

    // pdf start
    c->Print(fname_out_pdf + "[");

    // --- Get histograms
    TH1 *h_mass_allcuts = (TH1 *)f->Get("h_mass_allcuts");
    TH1 *h_mass_base = (TH1 *)f->Get("h_mass_base");

    TH1 *h_mass_ss_allcuts = (TH1 *)f->Get("h_mass_ss_allcuts");
    TH1 *h_mass_ss_middle_eop = (TH1 *)f->Get("h_mass_ss_middle_eop");

    TH1 *h_mass_allc_allcuts = (TH1 *)f->Get("h_mass_allc_allcuts");
    TH1 *h_mass_allc_middle_eop = (TH1 *)f->Get("h_mass_allc_middle_eop");

    TH1 *h_pt_all = (TH1 *)f->Get("h_pt_all");
    TH1 *h_pt_cut = (TH1 *)f->Get("h_pt_cut");
    TH1 *h_pt_chi2 = (TH1 *)f->Get("h_pt_chi2");
    TH1 *h_pt_chi2_eop = (TH1 *)f->Get("h_pt_chi2_eop");
    TH1 *h_mass_only_pt = (TH1 *)f->Get("h_mass_only_pt");

    // --- newly added: dz0v
    TH1 *h_dz0v_all = (TH1 *)f->Get("h_dz0v_all");
    TH1 *h_dz0v_cut = (TH1 *)f->Get("h_dz0v_cut");
    TH1 *h_mass_only_z0 = (TH1 *)f->Get("h_mass_only_z0");

    TH1 *h_nhits_cut = (TH1 *)f->Get("h_nhits_cut");
    TH1 *h_mass_only_hits = (TH1 *)f->Get("h_mass_only_hits");

    TH1 *h_chi2_all = (TH1 *)f->Get("h_chi2_all");
    TH1 *h_chi2_cut = (TH1 *)f->Get("h_chi2_cut");
    TH1 *h_mass_only_chi2 = (TH1 *)f->Get("h_mass_only_chi2");

    TH1 *h_dz_all = (TH1 *)f->Get("h_dz_all");
    TH1 *h_dz_cut = (TH1 *)f->Get("h_dz_cut");
    TH1 *h_mass_only_dz = (TH1 *)f->Get("h_mass_only_dz");

    TH1 *h_dphi_all = (TH1 *)f->Get("h_dphi_all");
    TH1 *h_dphi_cut = (TH1 *)f->Get("h_dphi_cut");
    TH1 *h_mass_only_dphi = (TH1 *)f->Get("h_mass_only_dphi");

    TH1 *h_zv_all = (TH1 *)f->Get("h_zv_all");
    TH1 *h_zv_cut = (TH1 *)f->Get("h_zv_cut");
    TH1 *h_mass_only_zvtx = (TH1 *)f->Get("h_mass_only_zvtx");

    TH1 *h_eop_all = (TH1 *)f->Get("h_eop_all");
    TH1 *h_eop_cut = (TH1 *)f->Get("h_eop_cut");
    TH1 *h_mass_only_eop = (TH1 *)f->Get("h_mass_only_eop");

    TH2 *h_eta_vs_phi = (TH2 *)f->Get("h_eta_vs_phi");
    TH2 *h_eta_vs_pt = (TH2 *)f->Get("h_eta_vs_pt");

    TH2 *h_dphi_vs_pt = (TH2 *)f->Get("h_dphi_vs_pt");
    TH2 *h_eop_vs_pt = (TH2 *)f->Get("h_eop_vs_pt");
    TH2 *h_dz_vs_pt = (TH2 *)f->Get("h_dz_vs_pt");

    TH2 *h_dphi_vs_tkgphi = (TH2 *)f->Get("h_dphi_vs_tkgphi");
    TH2 *h_dphi_vs_calophi = (TH2 *)f->Get("h_dphi_vs_calophi");
    TH2 *h_dphiraw_vs_tkgphi = (TH2 *)f->Get("h_dphiraw_vs_tkgphi");
    TH2 *h_dphiraw_vs_calophi = (TH2 *)f->Get("h_dphiraw_vs_calophi");

    TH2 *h_dphiraw_vs_phi0 = (TH2 *)f->Get("h_dphiraw_vs_phi0");
    TH2 *h_dzraw_vs_z0 = (TH2 *)f->Get("h_dzraw_vs_z0");

    TH2 *h_trkrphi_vs_calophi = (TH2 *)f->Get("h_trkrphi_vs_calophi");
    TH2 *h_trkrz_vs_caloz = (TH2 *)f->Get("h_trkrz_vs_caloz");
    TProfile *pf_trkrz_vs_caloz = (TProfile *)f->Get("pf_trkrz_vs_caloz");
    TF1 *fit_trkrz_vs_caloz = (TF1 *)f->Get("fit_trkrz_vs_caloz");

    TH2 *h_dzraw_vs_caloz = (TH2 *)f->Get("h_dzraw_vs_caloz");
    TH2 *h_dzraw_vs_trkrz = (TH2 *)f->Get("h_dzraw_vs_trkrz");

    TH1 *h_dz_raw = (TH1 *)f->Get("h_dz_raw");
    TH1 *h_dphi_raw = (TH1 *)f->Get("h_dphi_raw");

    TH2 *h_dz_vs_dphi = (TH2 *)f->Get("h_dz_vs_dphi");
    TH1 *h_dphi_raw_cm = (TH1 *)f->Get("h_dphi_raw_cm");

    // --- newly added: corrected dphi (cm) + vs phi0
    TH1 *h_dphi_corr_cm = (TH1 *)f->Get("h_dphi_corr_cm");
    TH2 *h_dphicorr_vs_phi0 = (TH2 *)f->Get("h_dphicorr_vs_phi0");

    TH2 *h_dz_vs_trkrz = (TH2 *)f->Get("h_dz_vs_trkrz");
    TH2 *h_dz_vs_caloz = (TH2 *)f->Get("h_dz_vs_caloz");

    TH1 *h_dz_corr = (TH1 *)f->Get("h_dz_corr");
    TH2 *h_dzcorr_vs_trkrz = (TH2 *)f->Get("h_dzcorr_vs_trkrz");
    TH2 *h_dzcorr_vs_caloz = (TH2 *)f->Get("h_dzcorr_vs_caloz");

    TH1 *h_mass_no_eop = (TH1 *)f->Get("h_mass_no_eop");
    TH1 *h_mass_middle_eop = (TH1 *)f->Get("h_mass_middle_eop");

    // ---------------- Draw sequence (match main macro) ----------------

    c->cd(1);
    if (h_mass_allcuts)
        h_mass_allcuts->Draw();
    c->cd(2);
    if (h_mass_base)
        h_mass_base->Draw();
    c->Print(fname_out_pdf);

    c->cd(1);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.03);
    gPad->SetLogy(1);
    if (h_pt_all)
        h_pt_all->GetXaxis()->SetRangeUser(0.0, 10.0);
    h_pt_all->Draw();
    if (h_pt_cut)
    {
        h_pt_cut->SetLineColor(TColor::GetColor(50, 205, 50));
        h_pt_cut->GetXaxis()->SetRangeUser(0.0, 10.0);
        h_pt_cut->Draw("same");
    }
    c->cd(2);
    gPad->SetLogy(1);

    if (h_mass_base)
        h_mass_base->Draw();
    if (h_mass_only_pt)
    {
        h_mass_only_pt->SetLineColor(TColor::GetColor(50, 205, 50));
        h_mass_only_pt->Draw("same");
    }
    c->Print(fname_out_pdf);

    if(h_pt_all&&h_pt_chi2)
    {
        c->cd(1);
        gPad->SetLogy(1);

        h_pt_all->SetTitle("p_{T} distributions of calo-matched tracks; p_{T} [GeV/c]; Tracks");
        h_pt_all->Draw();

        h_pt_chi2->SetLineColor(TColor::GetColor(50, 205, 50));
        h_pt_chi2->Draw("SAME");

        c->cd(2);
        gPad->SetLogy(1);
        h_pt_all->SetTitle("p_{T} distributions of calo-matched tracks; p_{T} [GeV/c]; Tracks");
        h_pt_all->Draw();

        h_pt_chi2->SetLineColor(TColor::GetColor(50, 205, 50));
        h_pt_chi2->Draw("SAME");

        h_pt_chi2_eop->SetLineColor(TColor::GetColor(255, 105, 180));
        h_pt_chi2_eop->Draw("SAME");
        c->Print(fname_out_pdf);
    }


    // --- dz0v block (new)
    c->cd(1);
    gPad->SetLogy(0);
    if (h_dz0v_all)
    {
        h_dz0v_all->GetXaxis()->SetRangeUser(-3, 3); // 表示範囲だけ変える
        h_dz0v_all->Draw();
    }
    if (h_dz0v_cut)
    {
        h_dz0v_cut->GetXaxis()->SetRangeUser(-3, 3); // 念のため同じに
        h_dz0v_cut->SetLineColor(TColor::GetColor(50, 205, 50));
        h_dz0v_cut->Draw("same");
    }

    c->cd(2);
    gPad->SetLogy(1);
    if (h_mass_base)
        h_mass_base->Draw();
    if (h_mass_only_z0)
    {
        h_mass_only_z0->SetLineColor(TColor::GetColor(50, 205, 50));
        h_mass_only_z0->Draw("same");
    }
    c->Print(fname_out_pdf);

    c->cd(1);
    if (h_nhits_cut)
        h_nhits_cut->Draw();
    c->cd(2);
    gPad->SetLogy(1);
    if (h_mass_base)
        h_mass_base->Draw();
    if (h_mass_only_hits)
    {
        h_mass_only_hits->SetLineColor(TColor::GetColor(50, 205, 50));
        h_mass_only_hits->Draw("same");
    }
    c->Print(fname_out_pdf);

    c->cd(1);
    if (h_chi2_all)
        h_chi2_all->Draw();
    // if (h_chi2_cut)
    // {
    //     h_chi2_cut->SetLineColor(TColor::GetColor(50, 205, 50));
    //     h_chi2_cut->Draw("same");
    // }
    c->cd(2);
    gPad->SetLogy(1);
    if (h_mass_base)
        h_mass_base->Draw();
    if (h_mass_only_chi2)
    {
        h_mass_only_chi2->SetLineColor(TColor::GetColor(50, 205, 50));
        h_mass_only_chi2->Draw("same");
    }
    c->Print(fname_out_pdf);

    c->cd(1);
    if (h_dz_all)
        h_dz_all->Draw();
    c->cd(2);
    gPad->SetLogy(1);
    if (h_mass_base)
        h_mass_base->Draw();
    if (h_mass_only_dz)
    {
        h_mass_only_dz->SetLineColor(TColor::GetColor(50, 205, 50));
        h_mass_only_dz->Draw("same");
    }
    c->Print(fname_out_pdf);

    c->cd(1);
    if (h_dphi_all)
        h_dphi_all->Draw();
    if (h_dphi_cut)
    {
        h_dphi_cut->SetLineColor(TColor::GetColor(50, 205, 50));
        h_dphi_cut->Draw("same");
    }
    c->cd(2);
    if (h_mass_base)
        h_mass_base->Draw();
    if (h_mass_only_dphi)
    {
        h_mass_only_dphi->SetLineColor(TColor::GetColor(50, 205, 50));
        h_mass_only_dphi->Draw("same");
    }
    c->Print(fname_out_pdf);

    c->cd(1);
    if (h_zv_all)
        h_zv_all->Draw();
    // if (h_zv_cut)
    // {
    //     h_zv_cut->SetLineColor(TColor::GetColor(50, 205, 50));
    //     h_zv_cut->Draw("same");
    // }
    c->cd(2);
    if (h_mass_base)
        h_mass_base->Draw();
    if (h_mass_only_zvtx)
    {
        h_mass_only_zvtx->SetLineColor(TColor::GetColor(50, 205, 50));
        h_mass_only_zvtx->Draw("same");
    }
    c->Print(fname_out_pdf);

    c->cd(1);

    if (h_eop_all)
    {
        h_eop_all->GetXaxis()->SetRangeUser(0.0, 2.5);
        h_eop_all->Draw();
    }

    // if (h_eop_cut)
    // {
    //     h_eop_cut->SetLineColor(TColor::GetColor(50, 205, 50));
    //     h_eop_cut->Draw("same");
    // }
    c->cd(2);
    gPad->SetLogy(1);
    if (h_mass_base)
        h_mass_base->Draw();
    if (h_mass_only_eop)
    {
        h_mass_only_eop->SetLineColor(TColor::GetColor(50, 205, 50));
        h_mass_only_eop->Draw("same");
    }
    c->Print(fname_out_pdf);

    c->cd(1);
    if (h_eta_vs_phi)
        h_eta_vs_phi->Draw();
    c->cd(2);
    gPad->SetLogy(0);
    if (h_eta_vs_pt)
        h_eta_vs_pt->Draw();
    c->Print(fname_out_pdf);

    c->cd(1);
    if (h_dphi_vs_pt)
        h_dphi_vs_pt->Draw();
    gPad->SetLogz();
    c->cd(2);
    gPad->SetLogy(0);
    if (h_eop_vs_pt)
        h_eop_vs_pt->Draw();
    c->Print(fname_out_pdf);

    c->cd(1);
    gPad->SetLogz(0);
    if (h_dphi_vs_pt)
        h_dphi_vs_pt->Draw();
    c->cd(2);
    if (h_dz_vs_pt)
        h_dz_vs_pt->Draw();
    c->Print(fname_out_pdf);

    c->cd(1);
    if (h_dphi_vs_tkgphi)
        h_dphi_vs_tkgphi->Draw();
    c->cd(2);
    if (h_dphi_vs_calophi)
        h_dphi_vs_calophi->Draw();
    c->Print(fname_out_pdf);

    c->cd(1);
    if (h_dphiraw_vs_tkgphi)
        h_dphiraw_vs_tkgphi->Draw();
    c->cd(2);
    if (h_dphiraw_vs_calophi)
        h_dphiraw_vs_calophi->Draw();
    c->Print(fname_out_pdf);

    c->cd(1);
    if (h_dphiraw_vs_phi0)
        h_dphiraw_vs_phi0->Draw();
    c->cd(2);
    if (h_dzraw_vs_z0)
        h_dzraw_vs_z0->Draw();
    c->Print(fname_out_pdf);

    c->cd(1);
    if (h_trkrphi_vs_calophi)
        h_trkrphi_vs_calophi->Draw();
    c->cd(2);
    if (h_trkrz_vs_caloz)
        h_trkrz_vs_caloz->Draw();

    if (use_dz_correction)
    {
        if (pf_trkrz_vs_caloz)
            pf_trkrz_vs_caloz->Draw("same");
        if (fit_trkrz_vs_caloz)
            fit_trkrz_vs_caloz->Draw("same");
        gPad->Update();

        if (h_trkrz_vs_caloz)
        {
            TPaveStats *st = (TPaveStats *)h_trkrz_vs_caloz->GetListOfFunctions()->FindObject("stats");
            if (st)
            {
                st->SetX1NDC(0.60);
                st->SetX2NDC(0.90);
                st->SetY1NDC(0.10);
                st->SetY2NDC(0.30);
            }
        }
        gPad->Modified();
        gPad->Update();
    }
    c->Print(fname_out_pdf);

    c->cd(1);
    if (h_dzraw_vs_caloz)
        h_dzraw_vs_caloz->Draw();
    c->cd(2);
    if (h_dzraw_vs_trkrz)
        h_dzraw_vs_trkrz->Draw();
    c->Print(fname_out_pdf);

    c->cd(1);
    if (h_dz_raw)
        h_dz_raw->Draw();
    c->cd(2);
    if (h_dphi_raw)
        h_dphi_raw->Draw();
    c->Print(fname_out_pdf);

    c->cd(1);
    if (h_dz_vs_dphi)
        h_dz_vs_dphi->Draw();
    c->cd(2);
    if (h_dphi_raw_cm)
        h_dphi_raw_cm->Draw();
    c->Print(fname_out_pdf);

    // --- corrected dphi (cm) block (new)
    c->cd(1);
    if (h_dphi_corr_cm)
        h_dphi_corr_cm->Draw();
    c->cd(2);
    gPad->SetLeftMargin(0.13);
    gPad->SetRightMargin(0.02);
    gPad->Modified();
    gPad->Update();
    if (h_dphicorr_vs_phi0)
        h_dphicorr_vs_phi0->Draw();
    c->Print(fname_out_pdf);

    c->cd(1);
    gPad->SetLeftMargin(0.12);
    if (h_dz_vs_trkrz)
        h_dz_vs_trkrz->Draw();
    c->cd(2);
    if (h_dz_vs_caloz)
        h_dz_vs_caloz->Draw();
    c->Print(fname_out_pdf);

    c->cd(1);
    if (h_dz_corr)
        h_dz_corr->Draw();
    if (h_dz_cut)
    {
        h_dz_cut->SetLineColor(TColor::GetColor(50, 205, 50));
        h_dz_cut->Draw("same");
    }
    c->cd(2);
    if (h_dz_all)
        h_dz_all->Draw();
    c->Print(fname_out_pdf);

    c->cd(1);
    if (h_dzcorr_vs_trkrz)
        h_dzcorr_vs_trkrz->Draw();
    c->cd(2);
    if (h_dzcorr_vs_caloz)
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

    // same plots with log
    c->cd(1);
    gPad->SetLogy(1);
    h_mass_middle_eop->Draw();
    h_mass_ss_middle_eop->SetLineColor(TColor::GetColor(50, 205, 50));
    h_mass_ss_middle_eop->Draw("SAME");
    c->cd(2);
    gPad->SetLogy(1);
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

    // --- E/p-dependent mass overlay (same as main macro)
    c->cd(1);
    gPad->SetLogy(0);
    if (h_mass_no_eop)
    {
        h_mass_no_eop->SetLineColor(TColor::GetColor(50, 205, 50));
        h_mass_no_eop->Draw();
    }
    if (h_mass_allcuts)
        h_mass_allcuts->Draw("SAME");
    if (h_mass_middle_eop)
    {
        h_mass_middle_eop->SetLineColor(TColor::GetColor(255, 105, 180));
        h_mass_middle_eop->Draw("SAME");
    }

    c->cd(2);
    gStyle->SetOptStat(0);
    TH1 *h_noeop_log = h_mass_no_eop ? (TH1 *)h_mass_no_eop->Clone("h_mass_no_eop_log") : nullptr;
    TH1 *h_all_log = h_mass_allcuts ? (TH1 *)h_mass_allcuts->Clone("h_mass_allcuts_log") : nullptr;
    TH1 *h_mid_log = h_mass_middle_eop ? (TH1 *)h_mass_middle_eop->Clone("h_mass_middle_log") : nullptr;

    gPad->SetLogy(1);

    h_noeop_log->SetStats(0);
    // h_noeop_log->GetYaxis()->SetRangeUser(1, 10000);
    h_noeop_log->SetLineColor(TColor::GetColor(50, 205, 50));
    h_noeop_log->SetTitle("Mass distribution depending on the E/p cut");
    h_noeop_log->Draw();
    h_all_log->SetStats(0);
    h_all_log->Draw("SAME");
    h_mid_log->SetStats(0);
    h_mid_log->SetLineColor(TColor::GetColor(255, 105, 180));
    h_mid_log->Draw("SAME");

    c->Print(fname_out_pdf);

    c->cd(1);
    gPad->SetLogy(1);
    h_mass_allcuts->Draw();
    c->Print(fname_out_pdf);

    // pdf end
    c->Print(fname_out_pdf + "]");
    std::cout << "histograms are drawn to : " << fname_out_pdf << std::endl;

    f->Close();
}
