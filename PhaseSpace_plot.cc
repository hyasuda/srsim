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
    // TFile *myfile = new TFile("/home/had/hyasuda/srsim/run/data/musr_330.root","READ"); // for KEKCC directory
    // TFile *myfile = new TFile("/Users/YASUDA/data/muonLinac/srsim/musr_330.root","READ"); // for local directory run 330
    TFile *myfile = new TFile("/Users/YASUDA/data/muonLinac/srsim/musr_370.root","READ"); // for local directory run 370
    // TTree *tree = myfile->Get("t1") ;
    TTree *tree = (TTree*)myfile->Get("t1") ;
    tree->SetBranchAddress("muIniMomX",&muIniMomX);
    tree->SetBranchAddress("muIniPosX",&muIniPosX);
    tree->SetBranchAddress("muIniMomY",&muIniMomY);
    tree->SetBranchAddress("muIniPosY",&muIniPosY);
    tree->SetBranchAddress("muIniMomZ",&muIniMomZ);

    tree->SetBranchAddress("save_px", &save_px);
    tree->SetBranchAddress("save_x",  &save_x );
    tree->SetBranchAddress("save_py", &save_py);
    tree->SetBranchAddress("save_y",  &save_y );
    tree->SetBranchAddress("save_pz", &save_pz);

    int nEntry = tree->GetEntries();
    // nEntry = 1000; // for test
    nEntry = 100; // for test
    // nEntry = 10; // for test

    /////////// CHECK THE FILE CONTENTS //////////
    // cout << "CHECK THE FILE CONTENTS" << endl;
    cout << "muIniMomX" << "\t" << "muIniPosX" << "\t" << "muIniMomY" << "\t" << "muIniPosY" << "\t" << "muIniMomZ" << endl;
    for(int iEntry = 0 ; iEntry < nEntry ; iEntry++){
        tree->GetEntry(iEntry);
        cout << muIniMomX << "\t" << muIniPosX << "\t" << muIniMomY << "\t" << muIniPosY << "\t" << muIniMomZ << endl; // for debug
    }
    muIniMomX = 0;
    muIniPosX = 0;
    muIniMomY = 0;
    muIniPosY = 0;
    muIniMomZ = 0;

    ////////// DRAW 2D HISTOGRAM //////////

    double xp = 0;
    double yp = 0;
    double x  = 0;
    double y  = 0;
    double xm = 0;
    double ym = 0;

    double xp_sum = 0;
    double yp_sum = 0;
    double x_sum = 0;
    double y_sum = 0;
    double xm_sum = 0;
    double ym_sum = 0;

    ////////// PLOT THE INITIAL BEAM PHASE SPACE //////////
    // TH2D (const char *name, const char *title, Int_t nbinsx, Double_t xlow, Double_t xup, Int_t nbinsy, Double_t ylow, Double_t yup) //
    TH2D * h2_pmX = new TH2D("h2_pmX", "h2_pmX", 60, -0.3, 0.3, 60, -0.3, 0.3); // res : 0.01
    TH2D * h2_pmY = new TH2D("h2_pmY", "h2_pmY", 60, -0.3, 0.3, 60, -0.3, 0.3); // res : 0.01
    TH2D * h2_psX = new TH2D("h2_psX", "h2_psX", 60, -0.3, 0.3, 60, -0.03, 0.03);
    TH2D * h2_psY = new TH2D("h2_psY", "h2_psY", 60, -0.3, 0.3, 60, -0.03, 0.03);


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
        // cout << xp << "\t" << yp << endl;
        // cout << xm << "\t" << ym << endl;
        // cout << x << "\t" << xp << "\t" << y << "\t" << ym << endl;
        // cout << endl;

        xp_sum += xp;
        yp_sum += yp;
        x_sum  += x;
        y_sum  += y;
        xm_sum += xm;
        ym_sum += ym;

        h2_psX->Fill(x,xp);
        h2_psY->Fill(y,yp);

        h2_pmX->Fill(x,xm);
        h2_pmY->Fill(y,ym);

    }

    cout << "mean X "    << "\t" << "mean X' "    << "\t" << "mean Y "    << "\t" << "mean Y' "    << endl;
    cout << x_sum/nEntry << "\t" << xp_sum/nEntry << "\t" << y_sum/nEntry << "\t" << yp_sum/nEntry << endl;

    cout << "mean X "    << "\t" << "mean Px "    << "\t" << "mean Y "    << "\t" << "mean Py "    << endl;
    cout << x_sum/nEntry << "\t" << xm_sum/nEntry << "\t" << y_sum/nEntry << "\t" << ym_sum/nEntry << endl;

    // TCanvas *c1 = new TCanvas("c1","c1");
    // TCanvas *c2 = new TCanvas("c2","c2");
    // TCanvas *c3 = new TCanvas("c3","c3");
    // TCanvas *c4 = new TCanvas("c4","c4");
    TCanvas *c_init = new TCanvas("c_init","Phase Space of Initial Beam",1000,800);
    c_init->Divide(2,2);
    // c1->cd();
    c_init->cd(1);
    h2_psX->Draw("colz");
    // c2->cd();
    c_init->cd(2);
    h2_psY->Draw("colz");
    // c3->cd();
    c_init->cd(3);
    h2_pmX->Draw("colz");
    // c4->cd();
    c_init->cd(4);
    h2_pmY->Draw("colz");

    cout << "////////// BEAM AT TARGET //////////" << endl;

    xp = 0;
    yp = 0;
    x  = 0;
    y  = 0;
    xm = 0;
    ym = 0;
    xp_sum = 0;
    yp_sum = 0;
    x_sum = 0;
    y_sum = 0;
    xm_sum = 0;
    ym_sum = 0;

    // TH2D (const char *name, const char *title, Int_t nbinsx, Double_t xlow, Double_t xup, Int_t nbinsy, Double_t ylow, Double_t yup) //
    TH2D * h2_pmX_tgt = new TH2D("h2_pmX_tgt", "h2_pmX_tgt", 100, -60, 10, 60, -0.3, 0.3); // res : 0.01
    TH2D * h2_pmY_tgt = new TH2D("h2_pmY_tgt", "h2_pmY_tgt", 100, -0.3, 0.3, 60, -0.3, 0.3); // res : 0.01
    TH2D * h2_psX_tgt = new TH2D("h2_psX_tgt", "h2_psX_tgt", 100, -60, 10, 60, -0.3, 0.3);
    TH2D * h2_psY_tgt = new TH2D("h2_psY_tgt", "h2_psY_tgt", 100, -0.3, 0.3, 60, -0.3, 0.3);

    /////////// CHECK THE FILE CONTENTS //////////
    // cout << "CHECK THE FILE CONTENTS" << endl;
    cout << "save_px" << "\t" << "save_x" << "\t" << "save_py" << "\t" << "save_y" << "\t" << "save_pz" << endl;
    for(int iEntry = 0 ; iEntry < nEntry ; iEntry++){
        tree->GetEntry(iEntry);
        cout << save_px << "\t" << save_x << "\t" << save_py << "\t" << save_y << "\t" << save_pz << endl; // for debug
    }
    save_px = 0;
    save_x = 0;
    save_py = 0;
    save_y = 0;
    save_pz = 0;

    for(int iEntry = 0 ; iEntry < nEntry ; iEntry++){
        tree->GetEntry(iEntry);
        xp = save_px/save_pz;
        yp = save_py/save_pz;
        x  = save_x;
        y  = save_y;
        xm = save_px;
        ym = save_py;
        // cout << xp << "\t" << yp << endl;
        // cout << xm << "\t" << ym << endl;
        cout << xp << "\t";

        xp_sum += xp;
        yp_sum += yp;
        x_sum  += x;
        y_sum  += y;
        xm_sum += xm;
        ym_sum += ym;

        h2_psX_tgt->Fill(x,xp);
        h2_psY_tgt->Fill(y,yp);
        h2_pmX_tgt->Fill(x,xm);
        h2_pmY_tgt->Fill(y,ym);

    }

    cout << "////////// BEAM PARAMETER AT TARGET //////////" << endl;
    cout << "mean X "    << "\t" << "mean X' "    << "\t" << "mean Y "    << "\t" << "mean Y' "    << endl;
    cout << x_sum/nEntry << "\t" << xp_sum/nEntry << "\t" << y_sum/nEntry << "\t" << yp_sum/nEntry << endl;

    cout << "mean X "    << "\t" << "mean Px "    << "\t" << "mean Y "    << "\t" << "mean Py "    << endl;
    cout << x_sum/nEntry << "\t" << xm_sum/nEntry << "\t" << y_sum/nEntry << "\t" << ym_sum/nEntry << endl;

    // TCanvas *c5 = new TCanvas("c5","c5");
    // TCanvas *c6 = new TCanvas("c6","c6");
    // TCanvas *c7 = new TCanvas("c7","c7");
    // TCanvas *c8 = new TCanvas("c8","c8");
    TCanvas *c_tgt = new TCanvas("c_tgt","Phase Space at Target",1000,800);
    c_tgt->Divide(2,2);
    // c5->cd();
    c_tgt->cd(1);
    h2_psX_tgt->Draw("colz");
    // c6->cd();
    c_tgt->cd(2);
    h2_psY_tgt->Draw("colz");
    // c7->cd();
    c_tgt->cd(3);
    h2_pmX_tgt->Draw("colz");
    // c8->cd();
    c_tgt->cd(4);
    h2_pmY_tgt->Draw("colz");

}
