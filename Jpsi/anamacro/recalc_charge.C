// This macro is just for drawing the event display, not for recalculating the charge of track... sorry for your confusion.

#include <vector>
#include <iostream>
#include <cmath>
#include <TLorentzVector.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TH2F.h>
#include <TMarker.h>
#include <TColor.h>
#include <algorithm>
#include <vector>
#include <cmath>

using namespace std;
TCanvas *c;
const float electron_mass = 0.000511; // GeV/c^2

void draw_mvtx(int mode, int run_nu) // Drawing MVTX on event displays
{
    const double kLayer_radii[6] = {2.24, 2.67, 3.01, 3.46, 3.78, 4.21};
    double dx;
    if (run_nu == 1)
        dx = 0.;
    else
        dx = 0.75;
    if (mode == 0) // x-y plane
    {
        for (int i = 0; i < 6; i++)
        {
            auto circle = new TEllipse(dx, 0.0, kLayer_radii[i]);
            circle->SetLineColorAlpha(kGray, 0.5);
            circle->SetLineWidth(2);
            circle->SetFillStyle(0);
            circle->Draw("same");
        }
    }
    else if (mode == 1) // z-r plane
    {
        for (int j = 0; j < 2; j++)
        {
            for (int i = 0; i < 6; i++)
            {
                TLine *line = new TLine(-13.55, (2 * j - 1) * kLayer_radii[i], 13.55, (2 * j - 1) * kLayer_radii[i]);
                line->SetLineColorAlpha(kGray, 0.5);
                line->SetLineWidth(2);
                line->Draw("same");
            }
        }
    }
}

void draw_intt(int mode) // Drawing INTT on event displays
{
    const double kLayer_radii[4] = {7.1888, 7.800, 9.680, 10.330};
    if (mode == 0) // x-y plane
    {
        for (int i = 0; i < 4; i++)
        {
            auto circle = new TEllipse(0.0, 0.0, kLayer_radii[i]);
            circle->SetLineColorAlpha(kGray, 0.5);
            circle->SetLineWidth(2);
            circle->SetFillStyle(0);
            circle->Draw("same");
        }
    }
    else if (mode == 1) // z-r plane
    {
        for (int j = 0; j < 2; j++)
        {
            for (int i = 0; i < 4; i++)
            {
                TLine *line = new TLine(-22.8, (2 * j - 1) * kLayer_radii[i], 22.8, (2 * j - 1) * kLayer_radii[i]);
                line->SetLineColorAlpha(kGray, 0.5);
                line->SetLineWidth(2);
                line->Draw("same");
            }
        }
    }
}

void draw_emc(int mode) // Drawing emc on event displays
{
    const double kLayer_radii[2] = {0.90 * 100, 1.16 * 100};
    if (mode == 0) // x-y plane
    {
        for (int i = 0; i < 2; i++)
        {
            auto circle = new TEllipse(0.0, 0.0, kLayer_radii[i]);
            circle->SetLineColorAlpha(kGray, 0.5);
            circle->SetLineWidth(2);
            circle->SetFillStyle(0);
            circle->Draw("same");
        }
    }
    else if (mode == 1) // z-r plane
    {
        for (int j = 0; j < 2; j++)
        {
            auto box = new TBox(-130, -(2 * j - 1) * kLayer_radii[0], 130, -(2 * j - 1) * kLayer_radii[1]);
            box->SetLineColorAlpha(kGray, 0.5);
            box->SetLineWidth(2);
            box->SetFillStyle(0);
            box->Draw("same");
        }
    }
}

TH1 *draw_frame(int mode, int detect, int run_nu, int ievt) // Drawing frame on event displays
{
    TH1 *flame;
    string title;

    if (mode == 0) // x-y plane
    {
        if (run_nu == 1) // MC
            title = Form("x-y plane {MC event%d};x [cm];y [cm]", ievt);
        else // Data
            title = Form("x-y plane {run%d event%d};x [cm];y [cm]", run_nu, ievt);
        if (detect == 0) // INTT
        {
            flame = gPad->DrawFrame(-15, -15, 15, 15, title.c_str());
            draw_intt(0);
            draw_mvtx(0, run_nu);
        }
        else if (detect == 1) // emc
        {
            gPad->SetLeftMargin(0.12);
            flame = gPad->DrawFrame(-120, -120, 120, 120, title.c_str());
            draw_intt(0);
            draw_mvtx(0, run_nu);
            draw_emc(0);
        }
    }
    else if (mode == 1) // z-r plane
    {
        if (run_nu == 1) // MC
            title = Form("z-r plane {MC event%d};z [cm];r [cm]", ievt);
        else // Data
            title = Form("z-r plane {run%d event%d};z [cm];r [cm]", run_nu, ievt);
        if (detect == 0) // INTT + MVTX
        {
            flame = gPad->DrawFrame(-25, -15, 25, 15, title.c_str());
            draw_intt(1);
            draw_mvtx(1, run_nu);
        }
        else if (detect == 1) // emc
        {
            gPad->SetLeftMargin(0.12);
            flame = gPad->DrawFrame(-140, -120, 140, 120, title.c_str());
            draw_intt(1);
            draw_mvtx(1, run_nu);
            draw_emc(1);
        }
    }
    return flame;
}

void draw_clusters(int mode,
                   const vector<float> &xclus,
                   const vector<float> &yclus,
                   const vector<float> &zclus,
                   const vector<float> &rclus,
                   int color,
                   float z0,
                   float r0)
{
    if (xclus.size() == 0)
        return;
    TGraph *g_temp = new TGraph();
    g_temp->SetMarkerColorAlpha(color, 1.0);
    // g_temp->SetLineColor(color);
    g_temp->SetMarkerStyle(20);
    g_temp->SetMarkerSize(0.8);
    for (int i = 0; i < (int)xclus.size(); i++)
    {
        if (mode == 0)
        {
            g_temp->SetPoint(i, xclus[i], yclus[i]);
        }
        else if (mode == 1)
        {
            // line->Draw("same");
            g_temp->SetPoint(i, zclus[i], rclus[i]);
        }
    }
    if (mode == 1)
    {
        TMarker *mz0 = new TMarker(z0, r0, 34);
        mz0->SetMarkerColorAlpha(color, 0.9);
        mz0->SetMarkerSize(1.2);
        mz0->Draw("same");
    }
    g_temp->Draw("P same");
}

// void draw_track_polyline_zr(
//     float z0, float r0,
//     const vector<float> &zclus,
//     const vector<float> &rclus,
//     int color)
// {
//     const int n = (int)zclus.size();
//     if (n == 0)
//         return;

//     // ★ 半径 |r| が小さい順（内側→外側）に並べる
//     vector<int> idx(n);
//     for (int i = 0; i < n; i++)
//         idx[i] = i;

//     std::sort(idx.begin(), idx.end(),
//               [&](int a, int b)
//               {
//                   return std::fabs(rclus[a]) < std::fabs(rclus[b]);
//               });

//     // ★ (z0,r0) + クラスター点列で polyline
//     TGraph *g = new TGraph(n + 1);
//     g->SetPoint(0, z0, r0);

//     for (int i = 0; i < n; i++)
//     {
//         const int ic = idx[i];
//         g->SetPoint(i + 1, zclus[ic], rclus[ic]);
//     }

//     g->SetLineColor(color);
//     g->SetLineWidth(2);
//     g->Draw("L same");
// }

// 💬
// ---- 直線フィット（r = a z + b） ----
static bool FitLine_r_vs_z(const std::vector<double> &z, const std::vector<double> &r,
                           double &a, double &b)
{
    const int n = (int)z.size();
    if (n < 2)
        return false;

    double Sz = 0, Sr = 0, Szz = 0, Szr = 0;
    for (int i = 0; i < n; i++)
    {
        Sz += z[i];
        Sr += r[i];
        Szz += z[i] * z[i];
        Szr += z[i] * r[i];
    }

    const double denom = n * Szz - Sz * Sz;
    if (std::fabs(denom) < 1e-12)
        return false;

    a = (n * Szr - Sz * Sr) / denom;
    b = (Sr - a * Sz) / (double)n;
    return true;
}

// ---- 枠（zmin/zmax/rmin/rmax）との交点を作って、枠内の“見える線分”を描く ----
static bool DrawClippedLineZR(double a, double b,
                              double zMin, double zMax,
                              double rMin, double rMax,
                              int color, int width = 2)
{
    std::vector<std::pair<double, double>> pts; // (z,r)

    auto add_if_inside = [&](double z, double r)
    {
        if (z >= zMin - 1e-9 && z <= zMax + 1e-9 &&
            r >= rMin - 1e-9 && r <= rMax + 1e-9)
        {
            pts.emplace_back(z, r);
        }
    };

    // 1) z = zMin, zMax での r
    add_if_inside(zMin, a * zMin + b);
    add_if_inside(zMax, a * zMax + b);

    // 2) r = rMin, rMax での z（a=0 だと割れないので注意）
    if (std::fabs(a) > 1e-12)
    {
        add_if_inside((rMin - b) / a, rMin);
        add_if_inside((rMax - b) / a, rMax);
    }

    // 重複っぽい点を軽く間引く
    auto close_pt = [&](const std::pair<double, double> &p1,
                        const std::pair<double, double> &p2)
    {
        return (std::fabs(p1.first - p2.first) < 1e-6 &&
                std::fabs(p1.second - p2.second) < 1e-6);
    };
    std::vector<std::pair<double, double>> uniq;
    for (auto &p : pts)
    {
        bool ok = true;
        for (auto &q : uniq)
        {
            if (close_pt(p, q))
            {
                ok = false;
                break;
            }
        }
        if (ok)
            uniq.push_back(p);
    }
    pts.swap(uniq);

    if ((int)pts.size() < 2)
        return false;

    // “いちばん離れてる2点”を線分の端点として採用
    int i1 = 0, i2 = 1;
    double dmax = -1;
    for (int i = 0; i < (int)pts.size(); i++)
    {
        for (int j = i + 1; j < (int)pts.size(); j++)
        {
            const double dz = pts[i].first - pts[j].first;
            const double dr = pts[i].second - pts[j].second;
            const double d2 = dz * dz + dr * dr;
            if (d2 > dmax)
            {
                dmax = d2;
                i1 = i;
                i2 = j;
            }
        }
    }

    TLine *ln = new TLine(pts[i1].first, pts[i1].second,
                          pts[i2].first, pts[i2].second);
    ln->SetLineColor(color);
    ln->SetLineWidth(width);
    ln->Draw("same");
    return true;
}

// ---- 片方向（出発点→クラスター側）にだけ、枠まで延長して描く ----
static bool DrawClippedRayZR(
    double a, double b,
    double z0, double r0,     // ★出発点
    double zDir, double rDir, // ★方向決め用（クラスター代表点）
    double zMin, double zMax,
    double rMin, double rMax,
    int color, int width = 2)
{
    // 方向：クラスター側へ
    double vz = zDir - z0;
    double vr = rDir - r0;

    // もし代表点がほぼ出発点と同じなら、傾きから方向を作る
    if (std::fabs(vz) < 1e-9 && std::fabs(vr) < 1e-9)
    {
        vz = 1.0;
        vr = a * vz;
    }

    // 直線 r = a z + b に沿う方向に寄せる（vzの符号だけ使って整える）
    if (std::fabs(vz) < 1e-12)
    {
        // ほぼ垂直（zが動かない）なら r方向だけで決める
        vr = (vr >= 0) ? +1.0 : -1.0;
        vz = 0.0;
    }
    else
    {
        const double sgn = (vz > 0) ? +1.0 : -1.0;
        vz = sgn;
        vr = a * vz;
    }

    // t>=0 で枠に当たる交点候補を探す（最小の正のt）
    double best_t = 1e100;
    bool found = false;

    auto try_z = [&](double zB)
    {
        if (std::fabs(vz) < 1e-12)
            return;
        const double t = (zB - z0) / vz;
        if (t <= 0)
            return;
        const double r = r0 + vr * t;
        if (r < rMin - 1e-9 || r > rMax + 1e-9)
            return;
        if (t < best_t)
        {
            best_t = t;
            found = true;
        }
    };

    auto try_r = [&](double rB)
    {
        if (std::fabs(vr) < 1e-12)
            return;
        const double t = (rB - r0) / vr;
        if (t <= 0)
            return;
        const double z = z0 + vz * t;
        if (z < zMin - 1e-9 || z > zMax + 1e-9)
            return;
        if (t < best_t)
        {
            best_t = t;
            found = true;
        }
    };

    // 枠4辺との交点チェック
    try_z(zMin);
    try_z(zMax);
    try_r(rMin);
    try_r(rMax);

    if (!found)
        return false;

    const double z1 = z0 + vz * best_t;
    const double r1 = r0 + vr * best_t;

    TLine *ln = new TLine(z0, r0, z1, r1);
    ln->SetLineColor(color);
    ln->SetLineWidth(width);
    ln->Draw("same");
    return true;
}

// ---- トラック1本： (z0,r0)+cluster を直線フィットして、枠まで延長して描く ----
void draw_track_fitline_zr(
    float z0, float r0,
    const std::vector<float> &zclus,
    const std::vector<float> &rclus,
    double zMin, double zMax,
    double rMin, double rMax,
    int color)
{
    // フィット点を作る（z0も含める）
    std::vector<double> zz;
    std::vector<double> rr;
    zz.reserve(zclus.size() + 1);
    rr.reserve(rclus.size() + 1);

    zz.push_back((double)z0);
    rr.push_back((double)r0);

    for (int i = 0; i < (int)zclus.size(); i++)
    {
        zz.push_back((double)zclus[i]);
        rr.push_back((double)rclus[i]);
    }

    double a = 0, b = 0;
    if (!FitLine_r_vs_z(zz, rr, a, b))
        return;

    // ★クラスターがある方向だけに伸ばすため、代表点（平均）を作る
    double zMean = 0.0, rMean = 0.0;
    const int nclus = (int)zclus.size();
    if (nclus > 0)
    {
        for (int i = 0; i < nclus; i++)
        {
            zMean += (double)zclus[i];
            rMean += (double)rclus[i];
        }
        zMean /= (double)nclus;
        rMean /= (double)nclus;
    }
    else
    {
        // クラスターが無いなら何もしない
        return;
    }

    // ★片方向（出発点→クラスター側）にだけ、枠まで延長して描画
    DrawClippedRayZR(a, b, (double)z0, (double)r0, zMean, rMean, zMin, zMax, rMin, rMax, color, 2);
}

// =====================
// ★ Track arc drawer ★
// =====================
static TGraph *MakeArcXY(
    double x0_cm, double y0_cm,
    double px, double py,
    double q, double Bz_T,
    double phi_span,
    int npts)
{
    const double pt = std::hypot(px, py);
    if (pt <= 0 || Bz_T == 0 || q == 0)
        return nullptr;

    // 円の中心 (cm)
    const double factor = 100.0 / (q * 0.3 * Bz_T); // [cm/GeV]
    const double xc = x0_cm + (-py) * factor;
    const double yc = y0_cm + (px)*factor;

    // 半径 (cm)
    const double Rcm = std::hypot(x0_cm - xc, x0_cm - xc) == 0 ? std::hypot(x0_cm - xc, y0_cm - yc)
                                                               : std::hypot(x0_cm - xc, y0_cm - yc);

    // 開始角
    const double theta0 = std::atan2(y0_cm - yc, x0_cm - xc);

    // 回転方向（符号）
    const double sgn = (q * Bz_T > 0) ? +1.0 : -1.0;

    auto *g = new TGraph(npts);
    for (int i = 0; i < npts; ++i)
    {
        const double t = (npts == 1) ? 0.0 : (double)i / (double)(npts - 1);
        const double theta = theta0 + sgn * (t * phi_span);
        const double x = xc + Rcm * std::cos(theta);
        const double y = yc + Rcm * std::sin(theta);
        g->SetPoint(i, x, y);
    }
    return g;
}

static TGraph *MakeArcZR(
    double x0_cm, double y0_cm, double z0_cm,
    double px, double py, double pz,
    double q, double Bz_T,
    double phi_span,
    int npts,
    double zMax, // ★追加
    double rMax  // ★追加
)
{
    const double pt = std::hypot(px, py);
    if (pt <= 0 || Bz_T == 0 || q == 0)
        return nullptr;

    const double factor = 100.0 / (q * 0.3 * Bz_T); // [cm/GeV]
    const double xc = x0_cm + (-py) * factor;
    const double yc = y0_cm + (px)*factor;
    const double Rcm = std::hypot(x0_cm - xc, y0_cm - yc);

    const double theta0 = std::atan2(y0_cm - yc, x0_cm - xc);
    const double sgn = (q * Bz_T > 0) ? +1.0 : -1.0;

    auto *g = new TGraph();
    int ip = 0;

    for (int i = 0; i < npts; ++i)
    {
        const double t = (npts == 1) ? 0.0 : (double)i / (double)(npts - 1);

        const double dtheta = sgn * (t * phi_span);
        const double theta = theta0 + dtheta;

        const double x = xc + Rcm * std::cos(theta);
        const double y = yc + Rcm * std::sin(theta);

        // 進行距離 s [cm]
        const double s = Rcm * dtheta;

        const double z = z0_cm + (pz / pt) * s;

        double r = std::hypot(x, y);
        if (y < 0)
            r = -r; // あなたの r 表示に合わせる

        // ★描画枠を超えたら打ち切り
        if (std::fabs(z) > zMax)
            break;
        if (std::fabs(r) > rMax)
            break;

        g->SetPoint(ip++, z, r);
    }

    if (ip < 2)
    {
        delete g;
        return nullptr;
    }
    return g;
}

float wrapPi(float dphi)
{
    while (dphi >= M_PI)
        dphi -= 2 * M_PI;
    while (dphi < -M_PI)
        dphi += 2 * M_PI;
    return dphi;
}

void recalc_charge(int run_nu = 53879)
{
    string filename;
    // const string &filename = Form("/sphenix/tg/tg01/commissioning/INTT/work/mahiro/SIliconCalo/run24pp/macro/makedata/condor/data/53879/ana_53879_00000.root");
    // const string &filename = Form("/sphenix/tg/tg01/commissioning/INTT/work/mahiro/SIliconCalo/MC/ana/condor_merge/ana_merge_00000.root");
    if (run_nu == 1)
        filename = Form("/sphenix/tg/tg01/commissioning/INTT/work/mahiro/SIliconCalo/MC/ana_jpsi/merged_200k.root");
    // filename = Form("/sphenix/tg/tg01/commissioning/INTT/work/mahiro/SIliconCalo/MC/ana/condor_merge/ana_merge_00000.root");
    else
        filename = Form("/sphenix/tg/tg01/commissioning/INTT/work/mahiro/SIliconCalo/run24pp/macro/makedata/condor/data/53879/ana_%d_00000.root", run_nu);
    // filename = Form("/sphenix/tg/tg01/commissioning/INTT/work/mahiro/SIliconCalo/run24pp/ana/newgeo/calocalib/%d_merged_5m.root", run_nu);

    // const string &filename = Form("/sphenix/tg/tg01/commissioning/INTT/work/mahiro/SIliconCalo/run24pp/ana/newgeo/calocalib/recalc_charge/ana_%d_1000evt.root", run_nu);

    bool debug = false;
    debug = true;

    // =========================
    // ★ Track line drawing param ★
    // =========================
    const double Bz_T = -1.4;            // 磁場 [T]
    const int npts = 300;                // 円弧の点数（滑らかさ）
    const double phi_span = TMath::Pi(); // 描く角度（とりあえずπ）

    int nloop = 450;
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
    bool IS_DATA = (evtTree != nullptr);

    TString fname_out = filename.substr(filename.find_last_of("/"), filename.size());
    // fname_out.ReplaceAll("merged", "mass_reconstruction_PYTHIA"); //
    fname_out.ReplaceAll("ana", "eventdisplay"); //
    if (debug)
        fname_out.ReplaceAll(".root", "_debug.root");

    // if (!IS_DATA)
    fname_out.ReplaceAll("merged", "eventdisplay"); //

    fname_out = (IS_DATA ? "result/data/evtdis" : "result/sim/evtdis") + fname_out;
    TFile *outFile = new TFile(fname_out, "RECREATE");

    TString fname_out_pdf = fname_out;
    fname_out_pdf.ReplaceAll(".root", ".pdf");

    c = new TCanvas("c", "c", 1000, 500);
    c->Divide(2, 1);
    c->cd(1);
    c->Print(fname_out_pdf + "[");

    // input branches
    int evt = 0, calo_evt = 0;
    vector<int> *track_id = nullptr;
    vector<float> *track_phi = nullptr;
    vector<float> *track_px = nullptr;
    vector<float> *track_py = nullptr;
    vector<float> *track_pt = nullptr;
    vector<float> *track_pz = nullptr;
    vector<float> *track_eta = nullptr;
    vector<float> *track_x = nullptr;
    vector<float> *track_y = nullptr;
    vector<float> *track_z = nullptr;
    vector<int> *track_nmaps = nullptr;
    vector<int> *track_nintt = nullptr;
    vector<float> *track_phi_emc = nullptr;
    vector<float> *track_x_emc = nullptr;
    vector<float> *track_y_emc = nullptr;
    vector<float> *track_z_emc = nullptr;
    vector<int> *track_charge = nullptr;
    vector<float> *track_chi2ndf = nullptr;

    vector<float> *calo_phi = nullptr;
    vector<float> *calo_energy = nullptr;
    vector<float> *calo_x = nullptr;
    vector<float> *calo_y = nullptr;
    vector<float> *calo_z = nullptr;
    vector<float> *calo_chi2 = nullptr;

    int clus_evt = 0;
    vector<int> *Siclus_trackid = nullptr;
    vector<int> *Siclus_layer = nullptr;
    vector<float> *Siclus_x = nullptr;
    vector<float> *Siclus_y = nullptr;
    vector<float> *Siclus_z = nullptr;
    vector<int> *Siclus_t = nullptr;

    float_t xvtx = 0.0;
    float_t yvtx = 0.0;
    float_t zvtx = 0.0;
    // vector<float> *truth_pt = 0;

    trackTree->SetBranchAddress("evt", &evt);
    trackTree->SetBranchAddress("track_id", &track_id);
    trackTree->SetBranchAddress("phi0", &track_phi);
    trackTree->SetBranchAddress("px0", &track_px);
    trackTree->SetBranchAddress("py0", &track_py);
    trackTree->SetBranchAddress("pt0", &track_pt);
    trackTree->SetBranchAddress("pz0", &track_pz);
    trackTree->SetBranchAddress("eta0", &track_eta);
    trackTree->SetBranchAddress("x0", &track_x);
    trackTree->SetBranchAddress("y0", &track_y);
    trackTree->SetBranchAddress("z0", &track_z);
    trackTree->SetBranchAddress("nmaps", &track_nmaps);
    trackTree->SetBranchAddress("nintt", &track_nintt);
    trackTree->SetBranchAddress("phi_proj_emc", &track_phi_emc);
    trackTree->SetBranchAddress("x_proj_emc", &track_x_emc);
    trackTree->SetBranchAddress("y_proj_emc", &track_y_emc);
    trackTree->SetBranchAddress("z_proj_emc", &track_z_emc);
    trackTree->SetBranchAddress("charge", &track_charge);
    trackTree->SetBranchAddress("chi2ndf", &track_chi2ndf);

    // evtTree is present for data; still guard in case tree is missing
    if (evtTree)
    {
        evtTree->SetBranchAddress("xvtx", &xvtx);
        evtTree->SetBranchAddress("yvtx", &yvtx);
        evtTree->SetBranchAddress("zvtx", &zvtx);
    }

    caloTree->SetBranchAddress("calo_evt", &calo_evt);
    caloTree->SetBranchAddress("phi", &calo_phi);
    caloTree->SetBranchAddress("energy", &calo_energy);
    caloTree->SetBranchAddress("x", &calo_x);
    caloTree->SetBranchAddress("y", &calo_y);
    caloTree->SetBranchAddress("z", &calo_z);
    if (IS_DATA)
        caloTree->SetBranchAddress("chi2", &calo_chi2);

    SiClusTree->SetBranchAddress("evt", &clus_evt);
    SiClusTree->SetBranchAddress("Siclus_trackid", &Siclus_trackid);
    SiClusTree->SetBranchAddress("Siclus_layer", &Siclus_layer);
    SiClusTree->SetBranchAddress("Siclus_x", &Siclus_x);
    SiClusTree->SetBranchAddress("Siclus_y", &Siclus_y);
    SiClusTree->SetBranchAddress("Siclus_z", &Siclus_z);
    SiClusTree->SetBranchAddress("Siclus_t", &Siclus_t);

    Long64_t nentries = trackTree->GetEntries();
    for (Long64_t i = 0; i < nentries; i++)
    {
        cout << "--------- evt : " << i << "----------" << endl;

        trackTree->GetEntry(i);
        SiClusTree->GetEntry(i);
        if (evtTree)
            evtTree->GetEntry(i);
        if (debug && evtTree)
        {
            cout << "vtx (x,y,z) = (" << xvtx << ", " << yvtx << ", " << zvtx << ")" << endl;
        }

        if (evt != clus_evt)
        {
            std::cerr << "Warning: evt mismatch at entry " << i
                      << ": track evt = " << evt
                      << ", clus evt = " << clus_evt << std::endl;
            continue;
        }

        int nTrk = (track_pt ? (int)track_pt->size() : 0);
        int nClus = (Siclus_trackid ? (int)Siclus_trackid->size() : 0);
        // cout << "ntrk : " << nTrk << "  nclus : " << nClus << endl;

        vector<vector<float>> x_by_trk, y_by_trk;
        vector<vector<float>> z_by_trk, r_by_trk;
        vector<float> z0_by_trk;
        vector<float> r0_by_trk;

        vector<float> x0_used, y0_used, z0_used;
        vector<float> px_used, py_used, pz_used;
        vector<int> q_used;

        int it_used = 0;
        for (int it = 0; it < nTrk; it++)
        {
            cout << "--- track : " << it << "---" << endl;

            int trkid = (*track_id)[it];
            int charge_trk = (*track_charge)[it];
            // float track_pt = (*track_pt)[it];

            // if (charge_trk > 0 || (*track_pt)[it] > 1.0)
            //     // if ((*track_pt)[it] < 0.5)
            //     continue;
            x_by_trk.emplace_back();
            y_by_trk.emplace_back();
            z_by_trk.emplace_back();
            r_by_trk.emplace_back();

            x0_used.push_back((*track_x)[it]);
            y0_used.push_back((*track_y)[it]);
            z0_used.push_back((*track_z)[it]);

            px_used.push_back((*track_px)[it]);
            py_used.push_back((*track_py)[it]);
            pz_used.push_back((*track_pz)[it]);

            q_used.push_back(charge_trk);

            z0_by_trk.push_back((*track_z)[it]);
            float r0 = sqrt((*track_x)[it] * (*track_x)[it] + (*track_y)[it] * (*track_y)[it]);
            if ((*track_y)[it] < 0)
                r0 = -r0;
            r0_by_trk.push_back(r0);
            cout << "r0 : " << r0 << endl;
            // for recalculate charge

            float rmin1 = 1e9, rmin2 = 1e9, rmax = -1e9;
            float x_rmin1 = 0, y_rmin1 = 0;
            float x_rmin2 = 0, y_rmin2 = 0;
            float x_rmax = 0, y_rmax = 0;

            float phi_2nd_from_min = 0.0;
            float phi_max_from_min = 0.0;
            float dphi_out_2nd = 0.0;
            int nclus_trk = 0;

            for (int ic = 0; ic < nClus; ic++)
            {
                int trkid_clus = (*Siclus_trackid)[ic];
                if (trkid != trkid_clus)
                    continue;

                float clus_x = (*Siclus_x)[ic];
                float clus_y = (*Siclus_y)[ic];
                float clus_z = (*Siclus_z)[ic];
                float clus_r = sqrt(clus_x * clus_x + clus_y * clus_y);
                if (clus_y < 0)
                    clus_r = -clus_r;

                x_by_trk[it_used].push_back(clus_x);
                y_by_trk[it_used].push_back(clus_y);
                z_by_trk[it_used].push_back(clus_z);
                r_by_trk[it_used].push_back(clus_r);

                float r = sqrt(clus_x * clus_x + clus_y * clus_y);

                if (r < rmin1)
                {
                    rmin2 = rmin1;
                    x_rmin2 = x_rmin1;
                    y_rmin2 = y_rmin1;

                    rmin1 = r;
                    x_rmin1 = clus_x;
                    y_rmin1 = clus_y;
                }
                else if (r < rmin2)
                {
                    rmin2 = r;
                    x_rmin2 = clus_x;
                    y_rmin2 = clus_y;
                }

                if (r > rmax)
                {
                    rmax = r;
                    x_rmax = clus_x;
                    y_rmax = clus_y;
                }

                // if (debug)
                // {
                //     cout << "track id : " << trkid << "  cluster id : " << trkid_clus
                //          << "  clus (x,y,z) = (" << clus_x << ", " << clus_y << ", " << clus_z << ")"
                //          << "  phi of cluster from vtx = " << phi_from_vtx << " , degree = " << degree
                //          << endl;
                // }
                nclus_trk++;
            }
            // if (nclus_trk >= 3 && rmin2 < 9e8 && rmax > -9e8)
            // {
            //     phi_2nd_from_min = atan2(y_rmin2 - y_rmin1, x_rmin2 - x_rmin1);
            //     phi_max_from_min = atan2(y_rmax - y_rmin1, x_rmax - x_rmin1);
            //     dphi_out_2nd = wrapPi(phi_max_from_min - phi_2nd_from_min); // out - in
            //     int charge_recalc = (dphi_out_2nd < 0) ? +1 : -1;           // dphi<0 -> charge = +1, +track curves to right side (clockwise)
            //     if (debug)
            //     {
            //         cout << "trk " << trkid
            //              << "  phi(2ndin,out)=(" << phi_2nd_from_min << "," << phi_max_from_min << ")"
            //              << "  dphi_clus(out-in)=" << dphi_out_2nd
            //              << "  charge(old,new)=(" << charge_trk << "," << charge_recalc << ")"
            //              << endl;
            //     }
            // }
            it_used++;
        }

        c->cd(1);
        TH1 *fr_xy = draw_frame(0, 0, run_nu, i);
        TMarker *mv_xy = new TMarker(xvtx, yvtx, 29);
        if (run_nu != 1)
            mv_xy->Draw("same");

        c->cd(2);
        TH1 *fr_zr = draw_frame(1, 0, run_nu, i);
        TMarker *mv_zr = new TMarker(zvtx, sqrt(xvtx * xvtx + yvtx * yvtx), 29);
        if (run_nu != 1)
            mv_zr->Draw("same");

        c->cd(1);
        const int colPlus = TColor::GetColor("#EE7F86");  // やわらかい赤（サーモン寄り）
        const int colMinus = TColor::GetColor("#7FB0FF"); // やわらかい青（淡い水色寄り）

        double zMin = -25, zMax = 25, rMin = -15, rMax = 15;
        if (fr_zr)
        {
            zMin = fr_zr->GetXaxis()->GetXmin();
            zMax = fr_zr->GetXaxis()->GetXmax();
            rMin = fr_zr->GetYaxis()->GetXmin();
            rMax = fr_zr->GetYaxis()->GetXmax();
        }

        for (int it = 0; it < x_by_trk.size(); it++)
        {
            const int q = q_used[it]; // +1 / -1 / 0
            const int col = (q > 0) ? colPlus : (q < 0 ? colMinus : kGray + 1);

            c->cd(1);
            draw_clusters(0, x_by_trk[it], y_by_trk[it], z_by_trk[it], r_by_trk[it], col, z0_by_trk[it], r0_by_trk[it]);
            // draw_clusters(0, x_by_trk[it], y_by_trk[it], z_by_trk[it], r_by_trk[it], 2 + (it % 8), z0_by_trk[it], r0_by_trk[it]);
            c->cd(2);
            draw_clusters(1, x_by_trk[it], y_by_trk[it], z_by_trk[it], r_by_trk[it], col, z0_by_trk[it], r0_by_trk[it]);
            // draw_clusters(1, x_by_trk[it], y_by_trk[it], z_by_trk[it], r_by_trk[it], 2 + (it % 8), z0_by_trk[it], r0_by_trk[it]);
        }

        c->cd(1);
        for (int it = 0; it < (int)x0_used.size(); ++it)
        {
            const double q = (q_used[it] > 0) ? +1.0 : (q_used[it] < 0 ? -1.0 : 0.0);
            if (q == 0.0)
                continue;

            TGraph *gArc = MakeArcXY(
                x0_used[it], y0_used[it],
                px_used[it], py_used[it],
                q, Bz_T,
                phi_span, npts);
            if (!gArc)
                continue;
            const int lineCol = (q > 0) ? colPlus : colMinus;

            gArc->SetLineWidth(2);
            // gArc->SetLineColorAlpha(2 + (it % 8), 0.4); // クラスター点と色を揃える
            gArc->SetLineColor(lineCol); // クラスター点と色を揃える
            gArc->Draw("L same");
        }

        // --- ★ z-r の線（追加）
        // const double zMax = 25.0;
        // const double rMax = 15.0;

        // c->cd(2);
        // for (int it = 0; it < (int)x0_used.size(); ++it)
        // {
        //     const double q = (q_used[it] > 0) ? +1.0 : (q_used[it] < 0 ? -1.0 : 0.0);
        //     if (q == 0.0)
        //         continue;

        //     TGraph *gZR = MakeArcZR(
        //         x0_used[it], y0_used[it], z0_used[it],
        //         px_used[it], py_used[it], pz_used[it],
        //         q, Bz_T,
        //         phi_span, npts,
        //         zMax, rMax);
        //     if (!gZR)
        //         continue;

        //     gZR->SetLineWidth(2);
        //     gZR->SetLineColorAlpha(2 + (it % 8), 0.6);
        //     gZR->Draw("L same");
        // }

        c->cd(2);
        for (int it = 0; it < (int)z_by_trk.size(); it++)
        {
            const int q = q_used[it];
            const int col = (q > 0) ? colPlus : (q < 0 ? colMinus : kGray + 1);

            draw_track_fitline_zr(
                z0_by_trk[it], r0_by_trk[it],
                z_by_trk[it], r_by_trk[it],
                zMin, zMax, rMin, rMax,
                col);
        }

        c->Print(fname_out_pdf);

        if (debug && i >= nloop)
            break;
    }
    c->Print(fname_out_pdf + "]");
}
