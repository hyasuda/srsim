#include "TFile.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TH1.h"
#include "TH2.h"
#include "TRandom.h"
#include "TClassTable.h"
#include "TSystem.h"
#include "TROOT.h"

void PhaseSpace_init_plot(TString filename = "musr_370.root"){
// void PhaseSpace_init_plot(){


    gROOT->SetStyle("Plain");
    gROOT->Reset();

    ////////// DEFINE THE VARIABLES //////////
    double muIniMomX = 0;
    double muIniPosX = 0;
    double muIniMomY = 0;
    double muIniPosY = 0;
    double muIniMomZ = 0;


    ////////// GET THE TREE FILE //////////
    TString dir = "/Users/YASUDA/data/muonLinac/srsim/";
    dir.Append(filename);
    TFile *myfile = new TFile(dir,"READ");
    // TFile *myfile = new TFile("/Users/YASUDA/data/muonLinac/srsim/musr_370.root","READ"); // for local directory run 370
    TTree *tree = (TTree*)myfile->Get("t1") ;
    tree->SetBranchAddress("muIniMomX",&muIniMomX);
    tree->SetBranchAddress("muIniPosX",&muIniPosX);
    tree->SetBranchAddress("muIniMomY",&muIniMomY);
    tree->SetBranchAddress("muIniPosY",&muIniPosY);
    tree->SetBranchAddress("muIniMomZ",&muIniMomZ);

    // tree->Branch("muIniMomX",&muIniMomX);
    // tree->Branch("muIniPosX",&muIniPosX);
    // tree->Branch("muIniMomY",&muIniMomY);
    // tree->Branch("muIniPosY",&muIniPosY);
    // tree->Branch("muIniMomZ",&muIniMomZ);

    int nEntry = tree->GetEntries();
    // nEntry = 1000; // for test
    // nEntry = 100; // for test
    nEntry = 10; // for test

    /////////// CHECK THE FILE CONTENTS //////////
    cout << "CHECK THE FILE CONTENTS" << endl;
    cout << "muIniMomX" << "\t" << "muIniPosX" << "\t" << "muIniMomY" << "\t" << "muIniPosY" << "\t" << "muIniMomZ" << endl;
    for(int iEntry = 0 ; iEntry < nEntry ; iEntry++){
        tree->GetEntry(iEntry);
        cout << "iEntry = " << iEntry << endl;
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
    double x_sum  = 0;
    double y_sum  = 0;
    double xm_sum = 0;
    double ym_sum = 0;

    ////////// PLOT THE INITIAL BEAM PHASE SPACE //////////
    // TH2D (const char *name, const char *title, Int_t nbinsx, Double_t xlow, Double_t xup, Int_t nbinsy, Double_t ylow, Double_t yup) //
    TH2D * h2_pmX = new TH2D("h2_pmX", "h2_pmX", 60, -0.3, 0.3, 60, -0.3, 0.3); // res : 0.01
    TH2D * h2_pmY = new TH2D("h2_pmY", "h2_pmY", 60, -0.3, 0.3, 60, -0.3, 0.3); // res : 0.01
    TH2D * h2_psX = new TH2D("h2_psX", "h2_psX", 60, -0.3, 0.3, 60, -0.03, 0.03); // res : 0.001
    TH2D * h2_psY = new TH2D("h2_psY", "h2_psY", 60, -0.3, 0.3, 60, -0.03, 0.03); // res : 0.001


    cout << "x" << "\t" << "xp" << "\t" << "y" << "\t" << "yp" << endl;  // for debug
    for(int iEntry = 0 ; iEntry < nEntry ; iEntry++){
        tree->GetEntry(iEntry);
        xp = muIniMomX/muIniMomZ;
        yp = muIniMomY/muIniMomZ;
        x  = muIniPosX;
        y  = muIniPosY;
        xm = muIniMomX;
        ym = muIniMomY;

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

    cout << "PHASE SPACE OF INITIAL BEAM" << endl;
    cout << "mean X "    << "\t" << "mean X' "    << "\t" << "mean Y "    << "\t" << "mean Y' "    << endl;
    cout << x_sum/nEntry << "\t" << xp_sum/nEntry << "\t" << y_sum/nEntry << "\t" << yp_sum/nEntry << endl;
    cout << "mean X "    << "\t" << "mean Px "    << "\t" << "mean Y "    << "\t" << "mean Py "    << endl;
    cout << x_sum/nEntry << "\t" << xm_sum/nEntry << "\t" << y_sum/nEntry << "\t" << ym_sum/nEntry << endl;


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

}
