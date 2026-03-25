#include <cmath>
#include <TLorentzVector.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TLegend.h>
#include <TStyle.h>
#include <iostream>

using namespace std;
TCanvas *c;

void DrawBG()
{
    const std::string file1 = "/sphenix/tg/tg01/commissioning/INTT/work/mahiro/SIliconCalo/run24pp/macro/result/sim/pythia_mass_reco_100m_merged.root"; // p+p (100M)
    // const std::string file2 = "/sphenix/tg/tg01/commissioning/INTT/work/mahiro/SIliconCalo/run24pp/macro/result/sim/mass_reco_200k.root";               // J/psi single inclusive (200k)
    const std::string file2 = "/sphenix/tg/tg01/commissioning/INTT/work/mahiro/SIliconCalo/run24pp/macro/result/sim/dphidez/mass_reco_200k.root";               // J/psi single inclusive (200k)

    const char *h_mass_base = "h_mass_base";
    const char *h_mass_allcuts = "h_mass_allcuts";

    // --- Normalization parameter
    const double cs_jpsi = 180.1;                       // nb,   cross section of J/ψ → e+e-
    const double rcs_pp = 42.0e6;                       // 42mb = 42,000,000nb,  reaction cross section of p+p collision
    const double Npp = 1e8;                             // # of events for PYTHIA p+p
    const double Ec_jpsi_pp = Npp * (cs_jpsi / rcs_pp); // # of expected J/ψ count in p+p

    const double Br_ee = 0.0596;            // J/psi -> e+ e- branching ratio (5.96%)
    const double Nref = Ec_jpsi_pp / Br_ee; // Reference Number of produced J/psi for PYTHIA p+p in all events
    const double Ngen_jpsi = 200000.0;      // # of events for J/psi single sim (# of J/psi meson produced)
    const double fGun = Nref / Ngen_jpsi;   // scale factor for J/psi single gun sim

    c = new TCanvas("c", "c", 1000, 500);
    c->Divide(2, 1);
    gStyle->SetOptStat("0");

    TFile *f1 = TFile::Open(file1.c_str(), "READ");
    TFile *f2 = TFile::Open(file2.c_str(), "READ");
    if (!f1 || f1->IsZombie() || !f2 || f2->IsZombie())
    {
        std::cerr << "file cannot open" << endl;
        return;
    }

    // ---- Left pad: h_mass_base ----
    // gPad->SetLogy(); // Enable if background is large and you want log scale

    TH1 *h_massbase_pp = (TH1 *)f1->Get(h_mass_base);
    TH1 *h_massbase_gun = (TH1 *)f2->Get(h_mass_base);
    int col_bg = TColor::GetColor("#E6862E");

    if (!h_massbase_pp || !h_massbase_gun)
    {
        std::cerr << "we cannot find histogram" << h_mass_base << endl;
        // If you still want to draw the right pad even when the left is missing, do not return here
    }
    else
    {
        TH1 *hpp = (TH1 *)h_massbase_pp->Clone("h_massbase_pp_copy");
        hpp->SetDirectory(0);
        TH1 *hgun = (TH1 *)h_massbase_gun->Clone("h_massbase_gun_copy");
        hgun->SetDirectory(0);

        // ★ Key change: do nothing to p+p; scale only the gun by fGun
        if (fGun > 0)
            hgun->Scale(fGun);
        cout << "scale factor" << fGun << endl;
        // Styling
        hpp->SetLineColor(col_bg);
        // hpp->SetLineWidth(2);
        // hgun->SetLineColor(kRed + 1);
        // hgun->SetLineWidth(2);

        // // Axis labels (set if empty)
        // if (hpp->GetXaxis()->GetTitle() == nullptr || strlen(hpp->GetXaxis()->GetTitle()) == 0)
        //     hpp->GetXaxis()->SetTitle("M_{ee} [GeV]");
        // if (hpp->GetYaxis()->GetTitle() == nullptr || strlen(hpp->GetYaxis()->GetTitle()) == 0)
        //     hpp->GetYaxis()->SetTitle("Counts (scaled)");

        // Draw: p+p first, then gun on top
        // hpp->Draw("HIST");
        // hgun->Draw("HIST SAME");

        // auto leg1 = new TLegend(0.58, 0.70, 0.88, 0.88);
        // leg1->SetBorderSize(0);
        // leg1->AddEntry(hpp, "p+p (100M, raw)", "l");
        // leg1->AddEntry(hgun, Form("J/psi gun (200k)  x %.4g", fGun), "l");
        // leg1->Draw();
    }

    // ---- Right pad: h_mass_allcuts ----
    c->cd(1);

    // gPad->SetLogy(); // Enable if needed

    TH1 *h_mass_pp = (TH1 *)f1->Get(h_mass_allcuts);
    TH1 *h_mass_gun = (TH1 *)f2->Get(h_mass_allcuts);

    if (!h_mass_pp || !h_mass_gun)
    {
        std::cerr << "we cannot find histogram" << h_mass_allcuts << endl;
    }
    else
    {
        TH1 *hpp2 = (TH1 *)h_mass_pp->Clone("h_mass_pp_copy");
        hpp2->SetDirectory(0);
        TH1 *hgun2 = (TH1 *)h_mass_gun->Clone("h_mass_gun_copy");
        hgun2->SetDirectory(0);
        hpp2->SetTitle("Signal vs p+p background");

        // ★ Same: scale only the gun by fGun
        if (fGun > 0)
            hgun2->Scale(fGun);

        // Styling
        hpp2->SetLineColor(col_bg);
        // hpp2->SetLineWidth(2);
        // hgun2->SetLineColor(kRed + 1);
        // hgun2->SetLineWidth(2);

        // // Axis labels (set if empty)
        // if (hpp2->GetXaxis()->GetTitle() == nullptr || strlen(hpp2->GetXaxis()->GetTitle()) == 0)
        //     hpp2->GetXaxis()->SetTitle("M_{ee} [GeV]");
        // if (hpp2->GetYaxis()->GetTitle() == nullptr || strlen(hpp2->GetYaxis()->GetTitle()) == 0)
        //     hpp2->GetYaxis()->SetTitle("Counts (scaled)");

        // Draw
        hpp2->SetStats(0);
        hgun2->SetStats(0);

        hpp2->Draw("HIST");
        hgun2->Draw("HIST SAME");
        // hgun2->Draw("HIST");

        c->cd(2);
        gPad->SetLogy(1);
        hpp2->GetYaxis()->SetRangeUser(1e-1, 1000);
        hpp2->SetStats(0);
        hgun2->SetStats(0);

        hpp2->Draw("HIST");
        hgun2->Draw("HIST SAME");

        // auto leg2 = new TLegend(0.58, 0.70, 0.88, 0.88);
        // leg2->SetBorderSize(0);
        // leg2->AddEntry(hpp2, "p+p (100M, raw)", "l");
        // leg2->AddEntry(hgun2, Form("J/psi gun (200k)  x %.4g", fGun), "l");
        // leg2->Draw();
    }

    // Save (if "result" folder does not exist, switch to the safer line below)
    c->SaveAs("/sphenix/tg/tg01/commissioning/INTT/work/mahiro/SIliconCalo/run24pp/macro/result/sim/mass_BG.pdf");
    // c->SaveAs("mass_BG.pdf"); // Safer option when no folder exists

    // Cleanup
    f1->Close();
    f2->Close();
    delete f1;
    delete f2;
}
