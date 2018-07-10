#include "TFile.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TH1.h"
#include "TH2.h"
#include "TRandom.h"
#include "TClassTable.h"
#include "TSystem.h"
#include "TROOT.h"

// void PhaseSpace_tgt_plot(){
void PhaseSpace_tgt_plot(TString filename = "musr_370.root"){


    gROOT->SetStyle("Plain");
    gROOT->Reset();

    ////////// DEFINE THE VARIABLES //////////
    double muIniMomX;
    double muIniPosX;
    double muIniMomY;
    double muIniPosY;
    double muIniMomZ;
    // vector<Double_t> save_px[2];
    // vector<Double_t> save_x[2];
    // vector<Double_t> save_py[2];
    // vector<Double_t> save_y[2];
    // vector<Double_t> save_pz[2];
    double save_px[2];
    double save_x[2];
    double save_py[2];
    double save_y[2];
    double save_pz[2];
    // double save_px;
    // double save_x;
    // double save_py;
    // double save_y;
    // double save_pz;


    ////////// GET THE TREE FILE //////////
    TString dir = "/Users/YASUDA/data/muonLinac/srsim/";
    dir.Append(filename);
    TFile *myfile = new TFile(dir,"READ");
    // TFile *myfile = new TFile("/Users/YASUDA/data/muonLinac/srsim/musr_370.root","READ"); // for local directory run 370
    TTree *tree = (TTree*)myfile->Get("t1");
    // tree = (TTree*)myfile->Get("t1") ;
    tree->SetBranchAddress("save_px", &save_px);
    tree->SetBranchAddress("save_x",  &save_x );
    tree->SetBranchAddress("save_py", &save_py);
    tree->SetBranchAddress("save_y",  &save_y );
    tree->SetBranchAddress("save_pz", &save_pz);

    int nEntry = tree->GetEntries();
    // nEntry = 1000; // for test
    // nEntry = 100; // for test
    // nEntry = 10; // for test

    /////////// CHECK THE FILE CONTENTS //////////
    cout << "CHECK THE FILE CONTENTS" << endl;
    cout << "save_px"   << "\t" << "save_x"    << "\t" << "save_py"   << "\t" << "save_y"    << "\t" << "save_pz"   << endl;
    for(int iEntry = 0 ; iEntry < nEntry ; iEntry++){
        tree->GetEntry(iEntry);
        cout << "iEntry = " << iEntry << endl;
        cout << save_px[0]   << "\t" << save_x[0]    << "\t" << save_py[0]   << "\t" << save_y[0]    << "\t" << save_pz[0]   << endl; // for debug
        cout << save_px[1]   << "\t" << save_x[1]    << "\t" << save_py[1]   << "\t" << save_y[1]    << "\t" << save_pz[1]   << endl; // for debug
        // cout << save_px   << "\t" << save_x    << "\t" << save_py   << "\t" << save_y    << "\t" << save_pz   << endl; // for debug
    }
    // save_px   = 0;
    // save_x    = 0;
    // save_py   = 0;
    // save_y    = 0;
    // save_pz   = 0;

    // save_px   = 0;
    // save_x    = 0;
    // save_py   = 0;
    // save_y    = 0;
    // save_pz   = 0;

    ////////// DRAW 2D HISTOGRAM //////////

    double xp_tgt = 0;
    double yp_tgt = 0;
    double x_tgt  = 0;
    double y_tgt  = 0;
    double xm_tgt = 0;
    double ym_tgt = 0;

    double xp_tgt_sum = 0;
    double yp_tgt_sum = 0;
    double x_tgt_sum  = 0;
    double y_tgt_sum  = 0;
    double xm_tgt_sum = 0;
    double ym_tgt_sum = 0;


    ////////// PLOT THE INITIAL BEAM PHASE SPACE //////////
    // TH2D (const char *name, const char *title, Int_t nbinsx, Double_t xlow, Double_t xup, Int_t nbinsy, Double_t ylow, Double_t yup) //
    // TH2D * h2_pmX_tgt = new TH2D("h2_pmX_tgt", "h2_pmX_tgt", 500, -60, 10, 60, -0.3, 0.3); // res : 0.01
    // TH2D * h2_pmY_tgt = new TH2D("h2_pmY_tgt", "h2_pmY_tgt", 500, -0.3, 0.3, 60, -0.3, 0.3); // res : 0.01
    // TH2D * h2_psX_tgt = new TH2D("h2_psX_tgt", "h2_psX_tgt", 500, -60, 10, 60, -0.3, 0.3);
    // TH2D * h2_psY_tgt = new TH2D("h2_psY_tgt", "h2_psY_tgt", 500, -0.3, 0.3, 60, -0.3, 0.3);
    TH2D * h2_psX_tgt = new TH2D("h2_psX_tgt", "h2_psX_tgt", 60, -60, 0., 140, -0.08, 0.06);
    TH2D * h2_psY_tgt = new TH2D("h2_psY_tgt", "h2_psY_tgt", 60, -30, 30., 60, -0.03, 0.03);
    TH2D * h2_pmX_tgt = new TH2D("h2_pmX_tgt", "h2_pmX_tgt", 60, -60, 0., 140, -0.8, 0.6); // res : 0.01
    TH2D * h2_pmY_tgt = new TH2D("h2_pmY_tgt", "h2_pmY_tgt", 60, -30, 30., 60, -0.3, 0.3); // res : 0.01


    cout << "x" << "\t" << "xp" << "\t" << "y" << "\t" << "yp" << endl;  // for debug
    for(int iEntry = 0 ; iEntry < nEntry ; iEntry++){
        tree->GetEntry(iEntry);
        // xp_tgt = save_px/save_pz;
        // yp_tgt = save_py/save_pz;
        // x_tgt  = save_x;
        // y_tgt  = save_y;
        // xm_tgt = save_px;
        // ym_tgt = save_py;

        xp_tgt = save_px[0]/save_pz[0];
        yp_tgt = save_py[0]/save_pz[0];
        x_tgt  = save_x[0];
        y_tgt  = save_y[0];
        xm_tgt = save_px[0];
        ym_tgt = save_py[0];

        xp_tgt_sum += xp_tgt;
        yp_tgt_sum += yp_tgt;
        x_tgt_sum  += x_tgt;
        y_tgt_sum  += y_tgt;
        xm_tgt_sum += xm_tgt;
        ym_tgt_sum += ym_tgt;

        h2_psX_tgt->Fill(x_tgt,xp_tgt);
        h2_psY_tgt->Fill(y_tgt,yp_tgt);
        h2_pmX_tgt->Fill(x_tgt,xm_tgt);
        h2_pmY_tgt->Fill(y_tgt,ym_tgt);

        // h2_psX_tgt->Fill(save_x,save_px/save_pz);
        // h2_psY_tgt->Fill(save_y,save_py/save_pz);
        // h2_pmX_tgt->Fill(save_x,save_px);
        // h2_pmY_tgt->Fill(save_y,save_py);



    }

    cout << "PHASE SPACE AT TARGET" << endl;
    cout << "mean X "    << "\t" << "mean X' "    << "\t" << "mean Y "    << "\t" << "mean Y' "    << endl;
    cout << x_tgt_sum/nEntry << "\t" << xp_tgt_sum/nEntry << "\t" << y_tgt_sum/nEntry << "\t" << yp_tgt_sum/nEntry << endl;
    cout << "mean X "    << "\t" << "mean Px "    << "\t" << "mean Y "    << "\t" << "mean Py "    << endl;
    cout << x_tgt_sum/nEntry << "\t" << xm_tgt_sum/nEntry << "\t" << y_tgt_sum/nEntry << "\t" << ym_tgt_sum/nEntry << endl;

    TCanvas *c_tgt = new TCanvas("c_tgt","Phase Space at Target",1000,800);
    c_tgt->Divide(2,2);
    c_tgt->cd(1);
    h2_psX_tgt->Draw("colz");
    c_tgt->cd(2);
    h2_psY_tgt->Draw("colz");
    c_tgt->cd(3);
    h2_pmX_tgt->Draw("colz");
    c_tgt->cd(4);
    h2_pmY_tgt->Draw("colz");

}
