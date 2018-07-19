// This macro is to calculate the beam emittance //
// This file is made by H.Yasuda on July 10, 2018 //
#include "TFile.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TH1.h"
#include "TH2.h"
#include "TRandom.h"
#include "TClassTable.h"
#include "TSystem.h"
#include "TROOT.h"
// const int nh=2;
#define X 0
#define Y 1
// TH2F* fH[nh];
// TCanvas* c1;
const int    xnbin= 100;
const double xmin = -0.3;
const double xmax = 0.3;
const int    ynbin= 100;
const double ymin = -0.1;
const double ymax = 0.1;
// void emittance_calc(){
void emittance_calc(TString filename = "musr_370.root"){

    gROOT->SetStyle("Plain");
    gROOT->Reset();

    double muIniMomX;
    double muIniPosX;
    double muIniMomY;
    double muIniPosY;
    double muIniMomZ;
    double save_px[2];
    double save_x[2];
    double save_py[2];
    double save_y[2];
    double save_pz[2];

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

    ////////// READING THE TREE FILE //////////
    TString dir = "/Users/YASUDA/data/muonLinac/srsim/";
    dir.Append(filename);
    TFile *myfile = new TFile(dir, "READ");
    // TFile *myfile = new TFile("/Users/YASUDA/data/muonLinac/srsim/musr_370.root","READ"); // for local directory run 370
    TTree *tree = new TTree("tree","tree");
    tree = (TTree*)myfile->Get("t1") ;
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
    // nEntry = 100; // for test
    // nEntry = 10; // for test

    /////////// CHECK THE FILE CONTENTS //////////
    // cout << "CHECK THE FILE CONTENTS" << endl;
    // cout << "muIniMomX"  << "\t" << "muIniPosX" << "\t" << "muIniMomY"  << "\t" << "muIniPosY" << "\t" << "muIniMomZ"  << endl;
    // cout << "save_px[0]" << "\t" << "save_x[0]" << "\t" << "save_py[0]" << "\t" << "save_y[0]" << "\t" << "save_pz[0]" << endl;
    // cout << "save_px[1]" << "\t" << "save_x[1]" << "\t" << "save_py[1]" << "\t" << "save_y[1]" << "\t" << "save_pz[1]" << endl;
    // for(int iEntry = 0 ; iEntry < nEntry ; iEntry++){
    //     tree->GetEntry(iEntry);
    //     cout << "iEntry = " << iEntry << endl;
    //     cout << muIniMomX << "\t" << muIniPosX << "\t" << muIniMomY << "\t" << muIniPosY << "\t" << muIniMomZ << endl;
    //     cout << save_px[0]   << "\t" << save_x[0]    << "\t" << save_py[0]   << "\t" << save_y[0]    << "\t" << save_pz[0]   << endl;
    //     cout << save_px[1]   << "\t" << save_x[1]    << "\t" << save_py[1]   << "\t" << save_y[1]    << "\t" << save_pz[1]   << endl;
    // }

    // for(int i=0; i<nh; i++){
    //     fH[i] = new TH2F(Form("fH%d", i), "", xnbin, xmin, xmax, ynbin, ymin, ymax);
    // }
    // c1 = new TCanvas("c1","",10,10,800,600);


    ////////// CALCULATION OF BEAM CENTER //////////
    double x_m = 0;
    double y_m = 0;
    double xp_m = 0;
    double yp_m = 0;
    double zp_m = 0;

    double x_m_tgt = 0;
    double y_m_tgt = 0;
    double xp_m_tgt = 0;
    double yp_m_tgt = 0;

    for(int iEntry = 0 ; iEntry < nEntry ; iEntry++){
        tree->GetEntry(iEntry);
        x_m += muIniPosX;
        y_m += muIniPosY;
        xp_m += muIniMomX/muIniMomZ;
        yp_m += muIniMomY/muIniMomZ;

        x_m_tgt += save_x[0];
        y_m_tgt += save_y[0];
        xp_m_tgt += save_px[0]/save_pz[0];
        yp_m_tgt += save_py[0]/save_pz[0];
    }

    x_m = x_m/nEntry;
    y_m = y_m/nEntry;
    xp_m = xp_m/nEntry;
    yp_m = yp_m/nEntry;

    x_m_tgt = x_m_tgt/nEntry;
    y_m_tgt = y_m_tgt/nEntry;
    xp_m_tgt = xp_m_tgt/nEntry;
    yp_m_tgt = yp_m_tgt/nEntry;

    cout << "///// AVERAGE /////" << endl;
    cout << "x_m"  << "\t" << "y_m"  << endl;
    cout << "xp_m" << "\t" << "yp_m" << endl;
    cout << x_m  << "\t" << y_m  << endl;
    cout << xp_m << "\t" << yp_m << endl;
    cout << "x_m_tgt"  << "\t" << "y_m_tgt"  << endl;
    cout << "xp_m_tgt" << "\t" << "yp_m_tgt" << endl;
    cout << x_m_tgt  << "\t" << y_m_tgt  << endl;
    cout << xp_m_tgt << "\t" << yp_m_tgt << endl;

    ////////// CALCULATION OF EMITTANCE //////////
    double xx_sum   = 0;
    double xpxp_sum = 0;
    double xxp_sum  = 0;
    double yy_sum   = 0;
    double ypyp_sum = 0;
    double yyp_sum  = 0;
    double zm = 0;
    double zm_sum = 0;
    // double zmzm_mean = 0 // This is not used here because of no time data.
    double xx_tgt_sum   = 0;
    double xpxp_tgt_sum = 0;
    double xxp_tgt_sum  = 0;
    double yy_tgt_sum   = 0;
    double ypyp_tgt_sum = 0;
    double yyp_tgt_sum  = 0;
    double zm_tgt = 0;
    double zm_tgt_sum = 0;
    // double zmzm_tgt_mean = 0 // This is not used here because of no time data.



    for(int iEntry = 0 ; iEntry < nEntry ; iEntry++){
        tree->GetEntry(iEntry);
        xp = muIniMomX/muIniMomZ;
        yp = muIniMomY/muIniMomZ;
        zm = muIniMomZ;
        x  = muIniPosX;
        y  = muIniPosY;
        // xm = muIniMomX;
        // ym = muIniMomY;

        // fH[X] -> Fill(x, xp);
        // fH[Y] -> Fill(y, yp);

        xp_tgt = save_px[0]/save_pz[0];
        yp_tgt = save_py[0]/save_pz[0];
        zm_tgt = save_pz[0];
        x_tgt  = save_x[0];
        y_tgt  = save_y[0];
        // xm_tgt = save_px[0];
        // ym_tgt = save_py[0];

        xx_sum   += (x  - x_m)  * (x  - x_m );
        xpxp_sum += (xp - xp_m) * (xp - xp_m);
        xxp_sum  += (x  - x_m)  * (xp - xp_m);
        yy_sum   += (y  - y_m)  * (y  - y_m );
        ypyp_sum += (yp - yp_m) * (yp - yp_m);
        yyp_sum  += (y  - y_m)  * (yp - yp_m);
        zm_sum   += zm;

        xx_tgt_sum   += (x_tgt  - x_m_tgt)  * (x_tgt  - x_m_tgt );
        xpxp_tgt_sum += (xp_tgt - x_m_tgt)  * (xp_tgt - x_m_tgt );
        xxp_tgt_sum  += (x_tgt  - x_m_tgt)  * (xp_tgt - xp_m_tgt);
        yy_tgt_sum   += (y_tgt  - y_m_tgt)  * (y_tgt  - y_m_tgt );
        ypyp_tgt_sum += (yp_tgt - yp_m_tgt) * (yp_tgt - yp_m_tgt);
        yyp_tgt_sum  += (y_tgt  - yp_m_tgt) * (yp_tgt - yp_m_tgt);
        zm_tgt_sum   += zm_tgt;
    }

    cout << "Summary of Beam Data at INITIAL " << endl;
    cout << "Mean X^2" << "\t" << "Mean X'^2" << "\t" << "Mean Y^2" << "\t" << "Mean Y'^2" << endl;
    cout << xx_sum     << "\t" << xpxp_sum    << "\t" << yy_sum     << "\t" << ypyp_sum    << endl;
    cout << "Mean X*X'"<< "\t" << "Mean Y *Y'"<< endl;
    cout << xxp_sum    << "\t" << yyp_sum     << endl;
    cout << "Mean Pz" << endl;
    cout << zm_sum/nEntry << endl;
    cout << "nEntry " << endl;
    cout << nEntry << endl;

    cout << "Summary of Beam Data at TARGET " << endl;
    cout << "Mean X^2" << "\t" << "Mean X'^2" << "\t" << "Mean Y^2" << "\t" << "Mean Y'^2" << endl;
    cout << xx_tgt_sum     << "\t" << xpxp_tgt_sum    << "\t" << yy_tgt_sum     << "\t" << ypyp_tgt_sum    << endl;
    cout << "Mean X*X'"<< "\t" << "Mean Y *Y'"<< endl;
    cout << xxp_tgt_sum    << "\t" << yyp_tgt_sum     << endl;
    cout << "Mean Pz" << endl;
    cout << zm_tgt_sum/nEntry << endl;

    double ex_rms = 0; // RMS Emittance of PhaseSpace X
    double ey_rms = 0; // RMS Emittance of PhaseSpace Y
    double ex_tgt_rms = 0; // RMS Emittance of PhaseSpace X
    double ey_tgt_rms = 0; // RMS Emittance of PhaseSpace Y


    ex_rms = sqrt(xx_sum*xpxp_sum - xxp_sum*xxp_sum)/nEntry;
    ey_rms = sqrt(yy_sum*ypyp_sum - yyp_sum*yyp_sum)/nEntry;
    // ex_rms = sqrt(xx_sum*xpxp_sum/(nEntry*nEntry) - xxp_sum*xxp_sum/(nEntry*nEntry));
    // ey_rms = sqrt(yy_sum*ypyp_sum/(nEntry*nEntry) - yyp_sum*yyp_sum/(nEntry*nEntry));
    ex_tgt_rms = sqrt(xx_tgt_sum*xpxp_tgt_sum - xxp_tgt_sum*xxp_tgt_sum)/nEntry;
    ey_tgt_rms = sqrt(yy_tgt_sum*ypyp_tgt_sum - yyp_tgt_sum*yyp_tgt_sum)/nEntry;


    cout << "RMS Emitttance of X at INITIAL: " << ex_rms << endl;
    cout << "RMS Emitttance of Y at INITIAL: " << ey_rms << endl;
    cout << "RMS Emitttance of X at TARGET : " << ex_tgt_rms << endl;
    cout << "RMS Emitttance of Y at TARGET : " << ey_tgt_rms << endl;

    double beta  = 0.0796704;
    double gamma = 1.00319;
    // cout << "Normalized RMS Emitttance of X at INITIAL: " << ex_rms*beta*gamma*3.1415*3.1415*1000 << endl;
    // cout << "Normalized RMS Emitttance of Y at INITIAL: " << ey_rms*beta*gamma*3.1415*3.1415*1000 << endl;
    // cout << "Normalized RMS Emitttance of X at TARGET : " << ex_tgt_rms*beta*gamma*3.1415*3.1415*1000 << endl;
    // cout << "Normalized RMS Emitttance of Y at TARGET : " << ey_tgt_rms*beta*gamma*3.1415*3.1415*1000 << endl;
    cout << "Normalized RMS Emitttance of X at INITIAL: " << ex_rms*beta*gamma*1000 << endl;
    cout << "Normalized RMS Emitttance of Y at INITIAL: " << ey_rms*beta*gamma*1000 << endl;
    cout << "Normalized RMS Emitttance of X at TARGET : " << ex_tgt_rms*beta*gamma*1000 << endl;
    cout << "Normalized RMS Emitttance of Y at TARGET : " << ey_tgt_rms*beta*gamma*1000 << endl;

    cout << "Growth of Emittance of X : " << ex_tgt_rms/ex_rms << endl;
    cout << "Growth of Emittance of Y : " << ey_tgt_rms/ey_rms << endl;
    // cout << "Growth of Normalized Emittance of X : " << ex_tgt_rms*beta/ex_rms << end;
    // cout << "Growth of Normalized Emittance of X : " << ex_tgt_rms/ex_rms << end;


    // c1->Divide(2,1);
    // c1->cd(1); fH[X] -> Draw("colz");
    // c1->cd(2); fH[Y] -> Draw("colz");
    // c1->cd();
}
