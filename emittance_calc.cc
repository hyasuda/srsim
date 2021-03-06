// This macro is to calculate the beam emittance //
// This file is made by H.Yasuda on Aug. 16, 2018 //
#include "TFile.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TH1.h"
#include "TH2.h"
#include "TRandom.h"
#include "TClassTable.h"
#include "TSystem.h"
#include "TROOT.h"

// Physics constants //
const double RF = 324.;   // RF frequency [MHz]
const double m_mu = 105.66; // muon mass [MeV]
const double c = 299792458; // light speed [m/s]

// contents of tree file //
int nEntry;
TTree *tree;
TFile *myfile;
TString filename;
const int nval = 13;
double val[nval] = {};
double val_mean[nval] = {};
const int nvalemit = 9;
double valemit[2][nvalemit] = {}; // xx, x'x', xx', yy, y'y', yy', phiphi, WW, phiW
double valemit_mean[2][nvalemit] = {};
const int nxyz = 3;
double emit[2][nxyz] = {};

const double beta  = 0.0796704;
const double gam = 1.00319;

const double xmaxcut =  5.;
const double xmincut = -5.;


TH2D* h2[8];
TH2D* h2_EXY[4];

void read();
void check();
void mean();
void emit_calc();
void hist_set();
void hist_plot();
void EXY_plot();
// void cut();

// MAIN FUNCTION //
void emittance_calc(string input = "musr_561.root"){

    filename = input;

    cout << scientific;

    gROOT->SetStyle("Plain");
    gROOT->Reset();
    gStyle->SetPadGridX(1);
    gStyle->SetPadGridY(1);

    read();
    nEntry = tree->GetEntries();
    // nEntry = 2; // for test
    // check();
    mean();
    emit_calc();
    hist_set();
    hist_plot();
    // EXY_plot();
}

void read(){
    TString dir = "/Users/YASUDA/data/muonLinac/srsim/";
    dir.Append(filename);
    myfile = new TFile(dir, "READ");
    tree = new TTree("tree","tree");
    tree = (TTree*)myfile->Get("t1") ;

    tree->SetBranchAddress("muIniMomX",    &val[0]);  // MeV
    tree->SetBranchAddress("muIniPosX",    &val[1]);  // mm
    tree->SetBranchAddress("muIniMomY",    &val[2]);  // MeV
    tree->SetBranchAddress("muIniPosY",    &val[3]);  // mm
    tree->SetBranchAddress("muIniMomZ",    &val[4]);  // MeV
    tree->SetBranchAddress("save_px",      &val[5]);
    tree->SetBranchAddress("save_x",       &val[6]);
    tree->SetBranchAddress("save_py",      &val[7]);
    tree->SetBranchAddress("save_y",       &val[8]);
    tree->SetBranchAddress("save_pz",      &val[9]);  // MeV
    tree->SetBranchAddress("muIniTime",    &val[10]); // us
    tree->SetBranchAddress("muTargetTime", &val[11]); // us
    tree->SetBranchAddress("save_ke",      &val[12]); // MeV
}

/////////// CHECK THE FILE CONTENTS //////////
void check(){
    cout << "CHECK THE FILE CONTENTS" << endl;
    // cout << "muIniMomX"  << "\t" << "muIniPosX" << "\t" << "muIniMomY"  << "\t" << "muIniPosY" << "\t" << "muIniMomZ"  << endl;
    // cout << "save_px" << "\t" << "save_x" << "\t" << "save_py" << "\t" << "save_y" << "\t" << "save_pz" << endl;
    for(int iEntry = 0 ; iEntry < nEntry ; iEntry++){
        tree->GetEntry(iEntry);
        cout << "iEntry = " << iEntry << ", nval = " << nval << endl;
        cout << scientific;
        for(int i = 0 ; i < nval ; i++){
            cout << val[i] << "\t";
        }
        cout << endl;
    }
}

void mean(){
    cout << " Calculation of mean : nEntry = " << nEntry << endl;

    for(int iEntry = 0 ; iEntry < nEntry ; iEntry++){
        tree->GetEntry(iEntry);
        for(int ival = 0 ; ival < nval ; ival++){
            if(val[6] > 5. || val[6] < -5.) continue;
            val_mean[ival] += val[ival];
        }
    }

    cout << "muIniMomX \t muIniPosX \t muIniMomY \t muIniPosY \t muIniMomZ \t save_px \t save_x \t save_py \t save_y \t save_pz \t muIniTime \t muTargetTime" << endl;
    for(int ival = 0 ; ival < nval ; ival++){
        val_mean[ival] = val_mean[ival]/nEntry;
        cout << val_mean[ival] << "\t";
    }
    cout << endl;
}

void emit_calc(){

    cout << "////////// CALCULATION OF EMITTANCE VARIABLES MEAN //////////" << endl;

    double xp_mean = 0;
    double yp_mean = 0;
    double xpout_mean = 0;
    double ypout_mean = 0;

    for(int iEntry = 0 ; iEntry < nEntry ; iEntry++){
        tree->GetEntry(iEntry);
        xp_mean += val[0]/val[4];
        yp_mean += val[2]/val[4];
        xpout_mean += val[5]/val[9];
        ypout_mean += val[7]/val[9];
    }
    xp_mean = xp_mean/nEntry;
    yp_mean = yp_mean/nEntry;
    xpout_mean = xpout_mean/nEntry;
    ypout_mean = ypout_mean/nEntry;

    int emitex = 0;

    for(int iEntry = 0 ; iEntry < nEntry ; iEntry++){
        tree->GetEntry(iEntry);

        // calculation of initial beam //
        valemit[0][0] = (val[1] - val_mean[1]) * (val[1] - val_mean[1]); //xx
        valemit[0][1] = (val[0]/val[4] - xp_mean) * (val[0]/val[4] - xp_mean); //xpxp
        valemit[0][2] = (val[1] - val_mean[1]) * (val[0]/val[4] - xp_mean); //xxp
        valemit[0][3] = (val[3] - val_mean[3]) * (val[3] - val_mean[3]); //yy
        valemit[0][4] = (val[2]/val[4] - yp_mean) * (val[2]/val[4] - yp_mean);//ypyp
        valemit[0][5] = (val[3] - val_mean[3]) * (val[2]/val[4] - yp_mean);//yyp
        valemit[0][6] = 2*3.1415*RF*(val[10] - val_mean[10]); //Delta_phi [deg]
        valemit[0][7] = beta*(val[4] - val_mean[4]); // Delta_W [MeV]
        valemit[0][8] = valemit[0][6]*valemit[0][7]; // Delta_phi * Delta_W
        valemit[0][6] = valemit[0][6]*valemit[0][6]; // (Delta_phi)^2
        valemit[0][7] = valemit[0][7]*valemit[0][7]; // (Delta_W)^2

        if(val[6] > xmaxcut || val[6] < xmincut) {
            emitex++;
            continue;
        }
        // calculation of output beam //
        valemit[1][0] = (val[6] - val_mean[6]) * (val[6] - val_mean[6]); //xx
        valemit[1][1] = (val[5]/val[9] - xpout_mean) * (val[5]/val[9] - xpout_mean); //xpxp
        valemit[1][2] = (val[6] - val_mean[6]) * (val[5]/val[9] - xpout_mean); //xxp
        valemit[1][3] = (val[8] - val_mean[8]) * (val[8] - val_mean[8]); //yy
        valemit[1][4] = (val[7]/val[9] - ypout_mean) * (val[7]/val[9] - ypout_mean);//ypyp
        valemit[1][5] = (val[8] - val_mean[8]) * (val[7]/val[9] - ypout_mean);//yyp
        valemit[1][6] = 2*3.1415*RF*(val[11] - val_mean[11]); //Delta_phi [deg]
        valemit[1][7] = beta*(val[9] - val_mean[9]); // Delta_W [MeV]
        valemit[1][8] = valemit[1][6]*valemit[1][7]; // Delta_phi * Delta_W
        valemit[1][6] = valemit[1][6]*valemit[1][6]; // (Delta_phi)^2
        valemit[1][7] = valemit[1][7]*valemit[1][7]; // (Delta_W)^2

        for(int ival = 0 ; ival < nvalemit ; ival++ ){
            valemit_mean[0][ival] += valemit[0][ival];
            valemit_mean[1][ival] += valemit[1][ival];
        }
    }
    nEntry = nEntry - emitex;
    cout << "emitex = " << emitex << endl;
    cout << "nEntry - emitex = " << nEntry << endl;
    cout << "initial" << endl;
    cout << "xx      \t xpxp      \t xxp      \t yy      \t ypyp      \t yyp      \t phiphi      \t WW      \t phiW " << endl;
    for(int ival = 0 ; ival < nvalemit ; ival++){
        valemit_mean[0][ival] = valemit_mean[0][ival]/nEntry;
        cout << valemit_mean[0][ival] << "\t";
    }
    cout << endl;
    cout << "output" << endl;
    cout << "xx      \t xpxp      \t xxp      \t yy      \t ypyp      \t yyp      \t phiphi      \t WW      \t phiW " << endl;
    for(int ival = 0 ; ival < nvalemit ; ival++){
        valemit_mean[1][ival] = valemit_mean[1][ival]/nEntry;
        cout << valemit_mean[1][ival] << "\t";
    }
    cout << endl;

    cout << "/////////// EMITTANCE ///////////" << endl;
    cout << "/////        INITIAL        /////" << endl;
    emit[0][0] = sqrt(valemit_mean[0][0]*valemit_mean[0][1] - valemit_mean[0][2]*valemit_mean[0][2]);
    emit[0][1] = sqrt(valemit_mean[0][3]*valemit_mean[0][4] - valemit_mean[0][5]*valemit_mean[0][5]);
    emit[0][2] = sqrt(valemit_mean[0][6]*valemit_mean[0][7] - valemit_mean[0][8]*valemit_mean[0][8]);
    cout << "x[mm mrad] \ty[mm mrad] \tz[deg MeV]" << endl;
    cout << emit[0][0] << "\t" << emit[0][1] << "\t" << emit[0][2] << endl;
    cout << "x,n[mm mrad] \ty,n[mm mrad] \tz,n[deg MeV]" << endl;
    cout << beta*gam*emit[0][0]*1.e+3 << "\t" << beta*gam*emit[0][1]*1.e+3 << "\t" << emit[0][2] << endl;

    cout << "/////        OUTPUT         /////" << endl;
    emit[1][0] = sqrt(valemit_mean[1][0]*valemit_mean[1][1] - valemit_mean[1][2]*valemit_mean[1][2]);
    emit[1][1] = sqrt(valemit_mean[1][3]*valemit_mean[1][4] - valemit_mean[1][5]*valemit_mean[1][5]);
    emit[1][2] = sqrt(valemit_mean[1][6]*valemit_mean[1][7] - valemit_mean[1][8]*valemit_mean[1][8]);
    cout << "x[mm mrad] \ty[mm mrad] \tz[deg MeV]" << endl;
    cout << emit[1][0] << "\t" << emit[1][1] << "\t" << emit[1][2] << endl;
    cout << "x,n[mm mrad] \ty,n[mm mrad] \tz,n[deg MeV]" << endl;
    cout << beta*gam*emit[1][0]*1.e+3 << "\t" << beta*gam*emit[1][1]*1.e+3 << "\t" << emit[1][2] << endl;

    cout << "/////        GROWTH         /////" << endl;
    cout << "x \t y \t z" << endl;
    cout << emit[1][0]/emit[0][0] << "\t" << emit[1][1]/emit[0][1] << "\t" << emit[1][2]/emit[0][2] << endl;
}


void hist_set(){
    int    xnbin = 100;
    double xmin  = -0.3;
    double xmax  = 0.3;
    int    ynbin = 100;
    double ymin  = -100;
    double ymax  = 100;
    int ih=0;

    // initial
    ih=0; xmin = -15. ; xmax = 15. ; ymin = -50. ; ymax = 50.;
    h2[ih] = new TH2D(Form("h2_%d", ih),Form("h2_%d", ih), xnbin, xmin, xmax, ynbin, ymin, ymax);
    ih=1; xmin = -15. ; xmax = 15. ; ymin = -50. ; ymax = 50.;
    h2[ih] = new TH2D(Form("h2_%d", ih),Form("h2_%d", ih), xnbin, xmin, xmax, ynbin, ymin, ymax);
    ih=2; xmin = -15. ; xmax = 15. ; ymin = -10. ; ymax = 10.;
    h2[ih] = new TH2D(Form("h2_%d", ih),Form("h2_%d", ih), xnbin, xmin, xmax, ynbin, ymin, ymax);
    ih=3; xmin = -180. ; xmax = 180. ; ymin = -0.02 ; ymax = 0.02;
    h2[ih] = new TH2D(Form("h2_%d", ih),Form("h2_%d", ih), xnbin, xmin, xmax, ynbin, ymin, ymax);

    // output
    // ih=4; xmin = -15. ; xmax = 15. ; ymin = -8. ; ymax = 8.;
    ih=4; xmin = -15. ; xmax = 15. ; ymin = -50. ; ymax = 50.;
    h2[ih] = new TH2D(Form("h2_%d", ih),Form("h2_%d", ih), xnbin, xmin, xmax, ynbin, ymin, ymax);
    ih=5; xmin = -25. ; xmax = 25. ; ymin = -50. ; ymax = 50.;
    h2[ih] = new TH2D(Form("h2_%d", ih),Form("h2_%d", ih), xnbin, xmin, xmax, ynbin, ymin, ymax);
    ih=6; xmin = -30. ; xmax = 30. ; ymin = -30. ; ymax = 30.;
    h2[ih] = new TH2D(Form("h2_%d", ih),Form("h2_%d", ih), xnbin, xmin, xmax, ynbin, ymin, ymax);
    ih=7; xmin = -180. ; xmax = 180. ; ymin = -0.02 ; ymax = 0.02;
    h2[ih] = new TH2D(Form("h2_%d", ih),Form("h2_%d", ih), xnbin, xmin, xmax, ynbin, ymin, ymax);

    for(int iEntry = 0 ; iEntry < nEntry ; iEntry++){
        tree->GetEntry(iEntry);

        // initial
        h2[0]->Fill(val[1], val[0]/val[4]*1.e+03);
        h2[1]->Fill(val[3], val[2]/val[4]*1.e+03);
        h2[2]->Fill(val[1], val[3]);
        h2[3]->Fill(2*3.1415*RF*(val[10]-val_mean[10]), beta*(val[4] - val_mean[4]));

        // output
        if(val[6] > xmaxcut || val[6] < xmincut) continue;
        h2[4]->Fill(val[6], val[5]/val[9]*1.e+03);
        h2[5]->Fill(val[8], val[7]/val[9]*1.e+03);
        h2[6]->Fill(val[6], val[8]);
        h2[7]->Fill(2*3.1415*RF*(val[11]-val_mean[11]), beta*(val[9] - val_mean[9]));
    }
}

void hist_plot(){
    TCanvas * c1 = new TCanvas("c1","input",800,600);
    TCanvas * c2 = new TCanvas("c2","output",800,600);
    c1->Divide(2,2);
    c2->Divide(2,2);
    h2[3]->GetXaxis()->SetRangeUser(-180.,180.);
    for(int i = 0 ; i < 4 ; i++){
        c1->cd(i+1);
        h2[i]->Draw("colz");
        h2[i]->GetXaxis()->SetNdivisions(504);
        h2[i]->GetYaxis()->SetNdivisions(504);
    }
    for(int i = 0 ; i < 4 ; i++){
        c2->cd(i+1);
        c2->SetGrid();
        h2[i+4]->Draw("colz");
        h2[i+4]->GetXaxis()->SetNdivisions(504);
        h2[i+4]->GetYaxis()->SetNdivisions(504);
    }
}


// For the search to the energy dependence of position //
void EXY_plot(){
    TCanvas * c_EXY = new TCanvas("c_EXY", "c_EXY", 1000, 800);
    for(int i = 0 ; i < 4 ; i++){
        h2_EXY[i] = new TH2D(Form("h2_EXY_%d", i), Form("h2_EXY_%d", i), 100, -15, 15, 100, -0.015, 0.015);
    }
    for(int iEntry = 0 ; iEntry < nEntry ; iEntry++){
        tree->GetEntry(iEntry);
        h2_EXY[0]->Fill(val[1], beta*(val[4] - val_mean[4]));
        h2_EXY[1]->Fill(val[3], beta*(val[4] - val_mean[4]));
        h2_EXY[2]->Fill(val[6], beta*(val[9] - val_mean[9]));
        h2_EXY[3]->Fill(val[8], beta*(val[9] - val_mean[9]));
    }
    c_EXY->Divide(2,2);
    for(int i = 0 ; i < 4 ; i++){
        c_EXY->cd(i+1);
        h2_EXY[i]->Draw("colz");
    }

}
