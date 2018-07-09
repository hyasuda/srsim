#include "TFile.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TH1.h"
#include "TH2.h"
#include "TRandom.h"
#include "TClassTable.h"
#include "TSystem.h"
#include "TROOT.h"

void PhaseSpace_plot(){


    gROOT->SetStyle("Plain");
    gROOT->Reset();

    ////////// DEFINE THE VARIABLES //////////
    double muIniMomX;
    double muIniPosX;
    double muIniMomY;
    double muIniPosY;
    double muIniMomZ;
    double save_px;
    double save_x;
    double save_py;
    double save_y;
    double save_pz;



    ////////// GET THE TREE FILE //////////
    TFile *myfile = new TFile("/Users/YASUDA/data/muonLinac/srsim/musr_370.root","READ"); // for local directory run 370
    TTree *tree = (TTree*)myfile->Get("t1") ;
    tree->SetBranchAddress("muIniMomX",&muIniMomX);
    tree->SetBranchAddress("muIniPosX",&muIniPosX);
    tree->SetBranchAddress("muIniMomY",&muIniMomY);
    tree->SetBranchAddress("muIniPosY",&muIniPosY);
    tree->SetBranchAddress("muIniMomZ",&muIniMomZ);

    int nEntry = tree->GetEntries();
    nEntry = 1000; // for test
    // nEntry = 100; // for test
    // nEntry = 10; // for test

    /////////// CHECK THE FILE CONTENTS //////////
    cout << "CHECK THE FILE CONTENTS" << endl;
    cout << "muIniMomX" << "\t" << "muIniPosX" << "\t" << "muIniMomY" << "\t" << "muIniPosY" << "\t" << "muIniMomZ" << endl;
    cout << "save_px"   << "\t" << "save_x"    << "\t" << "save_py"   << "\t" << "save_y"    << "\t" << "save_pz"   << endl;
    for(int iEntry = 0 ; iEntry < nEntry ; iEntry++){
        tree->GetEntry(iEntry);
        cout << "iEntry = " << iEntry << endl;
        cout << muIniMomX << "\t" << muIniPosX << "\t" << muIniMomY << "\t" << muIniPosY << "\t" << muIniMomZ << endl; // for debug
        cout << save_px   << "\t" << save_x    << "\t" << save_py   << "\t" << save_y    << "\t" << save_pz   << endl; // for debug
    }
    muIniMomX = 0;
    muIniPosX = 0;
    muIniMomY = 0;
    muIniPosY = 0;
    muIniMomZ = 0;
    save_px   = 0;
    save_x    = 0;
    save_py   = 0;
    save_y    = 0;
    save_pz   = 0;

    ////////// DRAW 2D HISTOGRAM //////////

    double xp = 0;
    double yp = 0;
    double x  = 0;
    double y  = 0;
    double xm = 0;
    double ym = 0;

    double xp_tgt = 0;
    double yp_tgt = 0;
    double x_tgt  = 0;
    double y_tgt  = 0;
    double xm_tgt = 0;
    double ym_tgt = 0;

    double xp_sum = 0;
    double yp_sum = 0;
    double x_sum  = 0;
    double y_sum  = 0;
    double xm_sum = 0;
    double ym_sum = 0;

    double xp_tgt_sum = 0;
    double yp_tgt_sum = 0;
    double x_tgt_sum  = 0;
    double y_tgt_sum  = 0;
    double xm_tgt_sum = 0;
    double ym_tgt_sum = 0;


    ////////// PLOT THE INITIAL BEAM PHASE SPACE //////////
    // TH2D (const char *name, const char *title, Int_t nbinsx, Double_t xlow, Double_t xup, Int_t nbinsy, Double_t ylow, Double_t yup) //
    TH2D * h2_pmX = new TH2D("h2_pmX", "h2_pmX", 60, -0.3, 0.3, 60, -0.3, 0.3); // res : 0.01
    TH2D * h2_pmY = new TH2D("h2_pmY", "h2_pmY", 60, -0.3, 0.3, 60, -0.3, 0.3); // res : 0.01
    TH2D * h2_psX = new TH2D("h2_psX", "h2_psX", 60, -0.3, 0.3, 60, -0.03, 0.03);
    TH2D * h2_psY = new TH2D("h2_psY", "h2_psY", 60, -0.3, 0.3, 60, -0.03, 0.03);
    TH2D * h2_pmX_tgt = new TH2D("h2_pmX_tgt", "h2_pmX_tgt", 100, -60, 10, 60, -0.3, 0.3); // res : 0.01
    TH2D * h2_pmY_tgt = new TH2D("h2_pmY_tgt", "h2_pmY_tgt", 100, -0.3, 0.3, 60, -0.3, 0.3); // res : 0.01
    TH2D * h2_psX_tgt = new TH2D("h2_psX_tgt", "h2_psX_tgt", 100, -60, 10, 60, -0.3, 0.3);
    TH2D * h2_psY_tgt = new TH2D("h2_psY_tgt", "h2_psY_tgt", 100, -0.3, 0.3, 60, -0.3, 0.3);


    // cout << "x" << "\t" << "xp" << "\t" << "y" << "\t" << "yp" << endl;  // for debug
    for(int iEntry = 0 ; iEntry < nEntry ; iEntry++){
        tree->GetEntry(iEntry);
        // cout << muIniMomX << "\t" << muIniMomZ << endl;
        xp = muIniMomX/muIniMomZ;
        yp = muIniMomY/muIniMomZ;
        x  = muIniPosX;
        y  = muIniPosY;
        xm = muIniMomX;
        ym = muIniMomY;

        xp_tgt = save_px/save_pz;
        yp_tgt = save_py/save_pz;
        x_tgt  = save_x;
        y_tgt  = save_y;
        xm_tgt = save_px;
        ym_tgt = save_py;

        xp_sum += xp;
        yp_sum += yp;
        x_sum  += x;
        y_sum  += y;
        xm_sum += xm;
        ym_sum += ym;

        xp_tgt_sum += xp_tgt;
        yp_tgt_sum += yp_tgt;
        x_tgt_sum  += x_tgt;
        y_tgt_sum  += y_tgt;
        xm_tgt_sum += xm_tgt;
        ym_tgt_sum += ym_tgt;

        h2_psX->Fill(x,xp);
        h2_psY->Fill(y,yp);
        h2_pmX->Fill(x,xm);
        h2_pmY->Fill(y,ym);

        h2_psX_tgt->Fill(x_tgt,xp_tgt);
        h2_psY_tgt->Fill(y_tgt,yp_tgt);
        h2_pmX_tgt->Fill(x_tgt,xm_tgt);
        h2_pmY_tgt->Fill(y_tgt,ym_tgt);

    }

    cout << "PHASE SPACE OF INITIAL BEAM" << endl;
    cout << "mean X "    << "\t" << "mean X' "    << "\t" << "mean Y "    << "\t" << "mean Y' "    << endl;
    cout << x_sum/nEntry << "\t" << xp_sum/nEntry << "\t" << y_sum/nEntry << "\t" << yp_sum/nEntry << endl;
    cout << "mean X "    << "\t" << "mean Px "    << "\t" << "mean Y "    << "\t" << "mean Py "    << endl;
    cout << x_sum/nEntry << "\t" << xm_sum/nEntry << "\t" << y_sum/nEntry << "\t" << ym_sum/nEntry << endl;
    cout << "PHASE SPACE AT TARGET" << endl;
    cout << "mean X "    << "\t" << "mean X' "    << "\t" << "mean Y "    << "\t" << "mean Y' "    << endl;
    cout << x_tgt_sum/nEntry << "\t" << xp_tgt_sum/nEntry << "\t" << y_tgt_sum/nEntry << "\t" << yp_tgt_sum/nEntry << endl;
    cout << "mean X "    << "\t" << "mean Px "    << "\t" << "mean Y "    << "\t" << "mean Py "    << endl;
    cout << x_tgt_sum/nEntry << "\t" << xm_tgt_sum/nEntry << "\t" << y_tgt_sum/nEntry << "\t" << ym_tgt_sum/nEntry << endl;

    TCanvas *c_init = new TCanvas("c_init","Phase Space of Initial Beam",1000,800);
    c_init->Divide(2,2);
    c_init->cd(1);
    h2_psX->Draw("colz");
    c_init->cd(2);
    h2_psY->Draw("colz");
    c_init->cd(3);
    h2_pmX->Draw("colz");
    c_init->cd(4);
    h2_pmY->Draw("colz");

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
