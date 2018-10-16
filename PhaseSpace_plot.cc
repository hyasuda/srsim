#include "TFile.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TH1.h"
#include "TH2.h"
#include "TRandom.h"
#include "TClassTable.h"
#include "TSystem.h"
#include "TROOT.h"

const double RF = 324. ; // Frequency of RF 324 MHz
const double m_mu = 105.6 ; // Muon mass 105.6 MeV

// void PhaseSpace_plot(){
void PhaseSpace_plot(TString filename = "musr_370.root"){


    gROOT->SetStyle("Plain");
    gROOT->Reset();

    ////////// DEFINE THE VARIABLES //////////
    double muIniMomX;
    double muIniPosX;
    double muIniMomY;
    double muIniPosY;
    double muIniMomZ;

    double muIniTime;
    double muTargetTime;
    double save_ke[2];

    double save_px[2];
    double save_x[2];
    double save_py[2];
    double save_y[2];
    double save_pz[2];

    double muIniPolX;
    double muIniPolY;
    double muIniPolZ;
    double muTargetPolX;
    double muTargetPolY;
    double muTargetPolZ;



    ////////// GET THE TREE FILE //////////
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
    tree->SetBranchAddress("muIniPolX", &muIniPolX);
    tree->SetBranchAddress("muIniPolY", &muIniPolY);
    tree->SetBranchAddress("muIniPolZ", &muIniPolZ);
    tree->SetBranchAddress("muTargetPolX", &muTargetPolX);
    tree->SetBranchAddress("muTargetPolY", &muTargetPolY);
    tree->SetBranchAddress("muTargetPolZ", &muTargetPolZ);
    tree->SetBranchAddress("muIniTime", &muIniTime);
    tree->SetBranchAddress("save_ke", &save_ke);
    tree->SetBranchAddress("muTargetTime", &muTargetTime);

    // int nEntry = tree->GetEntries();
    int nEntry = 10;

    /////////// CHECK THE FILE CONTENTS //////////
    // cout << "CHECK THE FILE CONTENTS" << endl;
    // cout << "muIniMomX"  << "\t" << "muIniPosX" << "\t" << "muIniMomY"  << "\t" << "muIniPosY" << "\t" << "muIniMomZ"  << endl;
    // cout << "save_px[0]" << "\t" << "save_x[0]" << "\t" << "save_py[0]" << "\t" << "save_y[0]" << "\t" << "save_pz[0]" << endl;
    // cout << "save_px[1]" << "\t" << "save_x[1]" << "\t" << "save_py[1]" << "\t" << "save_y[1]" << "\t" << "save_pz[1]" << endl;
    // for(int iEntry = 0 ; iEntry < nEntry ; iEntry++){
    //     tree->GetEntry(iEntry);
    //     cout << "iEntry = " << iEntry << endl;
    //     cout << muIniMomX << "\t" << muIniPosX << "\t" << muIniMomY << "\t" << muIniPosY << "\t" << muIniMomZ << endl; // for debug
    //     cout << save_px[0]   << "\t" << save_x[0]    << "\t" << save_py[0]   << "\t" << save_y[0]    << "\t" << save_pz[0]   << endl; // for debug
    //     cout << save_px[1]   << "\t" << save_x[1]    << "\t" << save_py[1]   << "\t" << save_y[1]    << "\t" << save_pz[1]   << endl; // for debug
    // }
    // muIniMomX = 0;
    // muIniPosX = 0;
    // muIniMomY = 0;
    // muIniPosY = 0;
    // muIniMomZ = 0;
    // save_px   = 0;
    // save_x    = 0;
    // save_py   = 0;
    // save_y    = 0;
    // save_pz   = 0;

    ////////// DRAW 2D HISTOGRAM //////////

    double xp = 0;
    double yp = 0;
    double x  = 0;
    double y  = 0;
    double xm = 0;
    double ym = 0;

    double t0 = 0;
    double t1 = 0;
    double phi0 = 0;
    double phi1 = 0;
    double ke0 = 0;
    double ke1 = 0;

    double t0_sum = 0;
    double t1_sum = 0;
    double phi0_sum = 0;
    double phi1_sum = 0;
    double ke0_sum = 0;
    double ke1_sum = 0;

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


    ////////// PLOT THE INITIAL BEAM PHASE SPACE AND BEAM POLARIZATION //////////

    TH1D * h1_polX = new TH1D("h1_polX", "h1_polX", 1, -0.00000001, 0.00000001);
    TH1D * h1_polY = new TH1D("h1_polY", "h1_polY", 1, -0.00000001, 0.00000001);
    TH1D * h1_polZ = new TH1D("h1_polZ", "h1_polZ", 500, -1, 1);
    TH1D * h1_polX_tgt = new TH1D("h1_polX_tgt", "h1_polX_tgt", 100, -0.08, 0.02);
    TH1D * h1_polY_tgt = new TH1D("h1_polY_tgt", "h1_polY_tgt", 100, -0.0001, 0.0001);
    TH1D * h1_polZ_tgt = new TH1D("h1_polZ_tgt", "h1_polZ_tgt", 100, 0.998, 1.);
    TH1D * h1_polZ_diff = new TH1D("h1_polZ_diff", "h1_polZ_diff", 500, -0.005, 0.005);

    // TH2D (const char *name, const char *title, Int_t nbinsx, Double_t xlow, Double_t xup, Int_t nbinsy, Double_t ylow, Double_t yup) //
    TH2D * h2_psX = new TH2D("h2_psX", "h2_psX", 100, -10., 10., 100, -0.0065, 0.0065);
    TH2D * h2_psY = new TH2D("h2_psY", "h2_psY", 100, -10., 10., 100, -0.006, 0.006);
    TH2D * h2_xy  = new TH2D("h2_xy",  "h2_xy" , 100, -10., 10., 100, -10. , 10. );
    // TH2D * h2_psZ = new TH2D("h2_psZ", "h2_psZ", 100, -180., 180., 100, -0.010, 0.010);
    TH2D * h2_psZ = new TH2D("h2_psZ", "h2_psZ", 100, -10., 10., 100, -0.015, 0.015);

    // ANY TIME //
    TH2D * h2_psX_tgt = new TH2D("h2_psX_tgt", "h2_psX_tgt", 100, -10., 10., 100, -0.0065, 0.0065);
    TH2D * h2_psY_tgt = new TH2D("h2_psY_tgt", "h2_psY_tgt", 100, -10., 10., 100, -0.0065, 0.0065);
    TH2D * h2_xy_tgt  = new TH2D("h2_xy_tgt",  "h2_xy_tgt" , 100, -10., 10., 100, -10. , 10. );
    // TH2D * h2_psZ_tgt = new TH2D("h2_psZ_tgt", "h2_psZ_tgt", 100, -180., 180., 100, -0.010, 0.010);
    TH2D * h2_psZ_tgt = new TH2D("h2_psZ_tgt", "h2_psZ_tgt", 100, -10., 10., 100, -0.015, 0.015);

    // TH2D * h2_psX_tgt = new TH2D();
    // TH2D * h2_psY_tgt = new TH2D();
    // TH2D * h2_xy_tgt  = new TH2D();

    if(!strcmp(filename, "musr_532.root")){
        delete h2_psX_tgt;
        delete h2_psY_tgt;
        delete h2_xy_tgt;
        h1_polY_tgt = new TH1D( "h1_polY_tgt", "h1_polY_tgt", 500,   -0.0001,  0.0001 );
        h2_psX_tgt  = new TH2D(  "h2_psX_tgt",  "h2_psX_tgt", 100,  -3.    , 3.    , 100,  -0.1 ,  0.1);
        h2_psY_tgt  = new TH2D(  "h2_psY_tgt",  "h2_psY_tgt", 100,  -3.    , 3.    , 100,  -0.1,  0.1);
        h2_xy_tgt   = new TH2D(       "h2_xy",      "h2_xy" , 100,  -3.    , 3.    , 100, -3.   , 3.);
    }

    if(!strcmp(filename, "musr_548.root")){
        delete h2_psX_tgt;
        delete h2_psY_tgt;
        delete h2_xy_tgt;
        h1_polY_tgt = new TH1D( "h1_polY_tgt", "h1_polY_tgt", 500,   -0.0001,  0.0001 );
        h2_psX_tgt  = new TH2D(  "h2_psX_tgt",  "h2_psX_tgt", 100,  -15.    , 15.    , 100,  -0.005 ,  0.005);
        h2_psY_tgt  = new TH2D(  "h2_psY_tgt",  "h2_psY_tgt", 100,  -10.    , 10.    , 100,  -0.006,  0.006);
        h2_xy_tgt   = new TH2D(       "h2_xy",      "h2_xy" , 100,  -15.    , 15.    , 100, -10.   , 10.);
    }

    if(!strcmp(filename, "musr_545.root")){
        delete h2_psX_tgt;
        delete h2_psY_tgt;
        delete h2_xy_tgt;
        h1_polY_tgt = new TH1D( "h1_polY_tgt", "h1_polY_tgt", 500,   -0.0001,  0.0001 );
        h2_psX_tgt  = new TH2D(  "h2_psX_tgt",  "h2_psX_tgt", 100,  -15.    , 15.    , 100,  -0.005 ,  0.005);
        h2_psY_tgt  = new TH2D(  "h2_psY_tgt",  "h2_psY_tgt", 100,  -10.    , 10.    , 100,  -0.006,  0.006);
        h2_xy_tgt   = new TH2D(       "h2_xy",      "h2_xy" , 100,  -15.    , 15.    , 100, -10.   , 10.);
    }

    if(!strcmp(filename, "musr_541.root")){
        delete h2_psX_tgt;
        delete h2_psY_tgt;
        delete h2_xy_tgt;
        h1_polY_tgt = new TH1D("h1_polY_tgt", "h1_polY_tgt", 500, -0.0001, 0.0001);
        h2_psX_tgt = new TH2D("h2_psX_tgt", "h2_psX_tgt", 100,   -40.,   40., 100, -0.01, 0.01);
        h2_psY_tgt = new TH2D("h2_psY_tgt", "h2_psY_tgt", 100,   -10.,  10., 100,  -0.006,  0.006);
        h2_xy_tgt  = new TH2D(     "h2_xy",     "h2_xy" , 100,   -40.,   40., 100,    -10,    10);
    }

    if(!strcmp(filename, "musr_530.root")){
        delete h2_psX_tgt;
        delete h2_psY_tgt;
        delete h2_xy_tgt;
        h1_polY_tgt = new TH1D("h1_polY_tgt", "h1_polY_tgt", 500, -0.001, 0.001);
        h2_psX_tgt = new TH2D("h2_psX_tgt", "h2_psX_tgt", 100,   -8.,   8., 100, -0.025, 0.025);
        h2_psY_tgt = new TH2D("h2_psY_tgt", "h2_psY_tgt", 100,  -30.,  30., 100,  -0.04,  0.04);
        h2_xy_tgt  = new TH2D(     "h2_xy",     "h2_xy" , 100,   -5.,   5., 100,    -30,    30);

    }

    if(!strcmp(filename, "musr_531.root")){
        delete h2_psX_tgt;
        delete h2_psY_tgt;
        delete h2_xy_tgt;
        h1_polY_tgt = new TH1D("h1_polY_tgt", "h1_polY_tgt", 500, -0.001, 0.001);
        h2_psX_tgt = new TH2D("h2_psX_tgt", "h2_psX_tgt", 100,   -8.,   8., 100, -0.025, 0.025);
        h2_psY_tgt = new TH2D("h2_psY_tgt", "h2_psY_tgt", 100,   -2.,   2., 100,  -0.04,  0.04);
        h2_xy_tgt  = new TH2D(     "h2_xy",     "h2_xy" , 100,   -5.,   5., 100,    -30,    30);

    }



    if(!strcmp(filename, "musr_511.root")){
        delete h2_psX_tgt;
        delete h2_psY_tgt;
        delete h2_xy_tgt;
        h1_polY_tgt = new TH1D("h1_polY_tgt", "h1_polY_tgt", 500, -0.001, 0.001);
        h2_psX_tgt = new TH2D("h2_psX_tgt", "h2_psX_tgt", 100,  -10.,  10., 100, -0.025, 0.025);
        h2_psY_tgt = new TH2D("h2_psY_tgt", "h2_psY_tgt", 100,  -30.,  30., 100,  -0.03,  0.03);
        h2_xy_tgt  = new TH2D(     "h2_xy",     "h2_xy" , 100,   -5.,   5., 100,    -30,    30);

    }

    if(!strcmp(filename, "musr_521.root")){
        delete h2_psX_tgt;
        delete h2_psY_tgt;
        delete h2_xy_tgt;
        h2_psX_tgt = new TH2D("h2_psX_tgt", "h2_psX_tgt", 100,  -20.,  40., 100, -0.03, 0.1);
        h2_psY_tgt = new TH2D("h2_psY_tgt", "h2_psY_tgt", 100, -100., 100., 100,  -0.3, 0.3);
        h2_xy_tgt  = new TH2D(     "h2_xy",     "h2_xy" , 100,  -30.,  60., 100,  -100, 100);

    }

    if(filename == "musr_370.root"){
        // delete h2_psX_tgt;
        // delete h2_psY_tgt;
        // delete h2_xy_tgt;
        h2_psX_tgt = new TH2D("h2_psX_tgt", "h2_psX_tgt", 60, -60, 0., 140, -0.08, 0.06);
        h2_psY_tgt = new TH2D("h2_psY_tgt", "h2_psY_tgt", 60, -30, 30., 60, -0.03, 0.03);
        h2_xy_tgt  = new TH2D("h2_xy",  "h2_xy" , 100, -60, 10, 100, -35 , 35 );
    }
    if(filename == "musr_430.root"){
        // delete h2_psX_tgt;
        // delete h2_psY_tgt;
        // delete h2_xy_tgt;
        h2_psX_tgt = new TH2D("h2_psX_tgt", "h2_psX_tgt", 100, -15., 55.,  100, -0.01, 0.08   );
        h2_psY_tgt = new TH2D("h2_psY_tgt", "h2_psY_tgt", 100, -150, 150., 100, -0.28, 0.28);
        h2_xy_tgt  = new TH2D("h2_xy_tgt",  "h2_xy_tgt" , 100, -15., 55,   100, -150., 150. );
    }
    if(filename == "musr_450.root"){
        // delete h2_psX_tgt;
        // delete h2_psY_tgt;
        // delete h2_xy_tgt;
        h2_psX_tgt = new TH2D("h2_psX_tgt", "h2_psX_tgt", 100, -5., 5.,  100, -0.03, 0.03   );
        h2_psY_tgt = new TH2D("h2_psY_tgt", "h2_psY_tgt", 100, -30, 30., 100, -0.03, 0.03);
        h2_xy_tgt  = new TH2D("h2_xy_tgt",  "h2_xy_tgt" , 100, -5., 5,   100, -30., 30. );
    }
    if(filename == "musr_500.root"){
        // delete h2_psX_tgt;
        // delete h2_psY_tgt;
        // delete h2_xy_tgt;
        h2_psX_tgt = new TH2D("h2_psX_tgt", "h2_psX_tgt", 100, -70., 0.,  100, -0.1, 0.1   );
        h2_psY_tgt = new TH2D("h2_psY_tgt", "h2_psY_tgt", 100, -30, 20., 100, -0.05, 0.05 );
        h2_xy_tgt  = new TH2D("h2_xy_tgt",  "h2_xy_tgt" , 100, -70., 0.,   100, -30., 20. );
    }

    for(int iEntry = 0 ; iEntry < nEntry ; iEntry++){
        tree->GetEntry(iEntry);
        t0 = muIniTime;
        t1 = muTargetTime;

        t0_sum += t0;
        t1_sum += t1;

        // ke0 = (muIniMomZ*muIniMomZ + muIniMomX*muIniMomX + muIniMomY*muIniMomY)/(2*m_mu);
        ke0 = muIniMomZ*muIniMomZ/(2*m_mu);
        ke1 = save_ke[0];

        ke0_sum += ke0;
        ke1_sum += ke1;
    }

    cout << "= = = AVERAGE = = =" << endl;
    cout << t0_sum/nEntry << "\t" << t1_sum/nEntry << endl;
    cout << ke0_sum/nEntry << "\t" << ke1_sum/nEntry << endl;


    // cout << "x" << "\t" << "xp" << "\t" << "y" << "\t" << "yp" << endl;  // for debug
    for(int iEntry = 0 ; iEntry < nEntry ; iEntry++){
        tree->GetEntry(iEntry);
        // cout << muIniMomX << "\t" << muIniMomZ << endl;
        xp = muIniMomX/muIniMomZ;
        yp = muIniMomY/muIniMomZ;
        x  = muIniPosX;
        y  = muIniPosY;

        xp_tgt = save_px[0]/save_pz[0];
        yp_tgt = save_py[0]/save_pz[0];
        x_tgt  = save_x[0];
        y_tgt  = save_y[0];
        xm_tgt = save_px[0];
        ym_tgt = save_py[0];

        t0 = muIniTime;
        t1 = muTargetTime;
        // phi0 = RF*(t0-t0_sum/nEntry)/180; // phase [deg]
        // phi1 = RF*(t1-t1_sum/nEntry)/180; // phase [deg]
        phi0 = RF*(t0-t0_sum/nEntry); // phase [deg]
        cout << "CHECK" << endl;
        cout << 1/RF << " [usec]" << endl;
        cout << t0-t0_sum/nEntry << " [usec]" << endl;
        phi1 = RF*(t1-t1_sum/nEntry); // phase [deg]

        // ke0 = (muIniMomZ*muIniMomZ + muIniMomX*muIniMomX+ muIniMomY*muIniMomY)/(2*m_mu);
        ke0 = muIniMomZ*muIniMomZ/(2*m_mu);
        ke1 = save_ke[0];

        h2_psZ->Fill(phi0,ke0-(ke0_sum/nEntry));
        h2_psZ_tgt->Fill(phi1,ke1-(ke1_sum/nEntry));

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

        h2_psX_tgt->Fill(x_tgt,xp_tgt);
        h2_psY_tgt->Fill(y_tgt,yp_tgt);

        h2_xy->Fill(x,y);
        h2_xy_tgt->Fill(x_tgt,y_tgt);

        h1_polX->Fill(muIniPolX/nEntry);
        h1_polY->Fill(muIniPolY/nEntry);
        h1_polZ->Fill(muIniPolZ/nEntry);
        h1_polX_tgt->Fill(muTargetPolX);
        h1_polY_tgt->Fill(muTargetPolY);
        h1_polZ_tgt->Fill(muTargetPolZ);
        h1_polZ_diff->Fill(muTargetPolZ + muIniPolZ);

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
    h2_psX->SetXTitle("x[mm]");
    h2_psX->SetYTitle("x'[rad]");
    h2_psX->Draw("colz");
    c_init->cd(2);
    h2_psY->SetXTitle("y[mm]");
    h2_psY->SetYTitle("y'[rad]");
    h2_psY->Draw("colz");
    c_init->cd(3);
    h2_xy->SetXTitle("x[mm]");
    h2_xy->SetYTitle("y[mm]");
    h2_xy->Draw("colz");
    c_init->cd(4);
    h2_psZ->Draw("colz");

    TCanvas *c_tgt = new TCanvas("c_tgt","Phase Space at Target",1000,800);
    c_tgt->Divide(2,2);
    c_tgt->cd(1);
    h2_psX_tgt->SetXTitle("x[mm]");
    h2_psX_tgt->SetYTitle("x'[rad]");
    h2_psX_tgt->Draw("colz");
    c_tgt->cd(2);
    h2_psY_tgt->SetXTitle("y[mm]");
    h2_psY_tgt->SetYTitle("y'[rad]");
    h2_psY_tgt->Draw("colz");
    c_tgt->cd(3);
    h2_xy_tgt->SetXTitle("x[mm]");
    h2_xy_tgt->SetYTitle("y[mm]");
    h2_xy_tgt->Draw("colz");
    c_tgt->cd(4);
    h2_psZ_tgt->Draw("colz");

    TCanvas * c_pol = new TCanvas("c_pol","Polarization of Beam", 1000, 800);
    c_pol->Divide(2,2);
    c_pol->cd(1);
    h1_polX_tgt->Scale(1./nEntry);
    h1_polX_tgt->SetLineColor(1);
    h1_polX_tgt->SetLineWidth(2);
    h1_polX_tgt->SetXTitle("X Polarization");
    h1_polX_tgt->SetYTitle("Event/Entry");
    h1_polX_tgt->Draw();
    h1_polX->SetLineColor(2);
    h1_polX->SetLineWidth(4);
    h1_polX->Draw("SAME");
    c_pol->cd(2);
    h1_polY_tgt->Scale(1./nEntry);
    h1_polY_tgt->SetLineColor(1);
    h1_polY_tgt->SetLineWidth(2);
    h1_polY_tgt->SetXTitle("Y Polarization");
    h1_polY_tgt->SetYTitle("Event/Entry");
    h1_polY_tgt->Draw();
    h1_polY->SetLineColor(2);
    h1_polY->SetLineWidth(4);
    h1_polY->Draw("SAME");
    c_pol->cd(3);
    h1_polZ_tgt->Scale(1./nEntry);
    h1_polZ_tgt->SetXTitle("Z Polarization");
    h1_polZ_tgt->SetYTitle("Event/Entry");
    h1_polZ_tgt->SetLineColor(1);
    h1_polZ_tgt->SetLineWidth(2);
    h1_polZ_tgt->Draw();
    // c_pol->cd(4);
    // h1_polZ_diff->Scale(1./nEntry);
    // h1_polZ_diff->SetLineColor(1);
    // h1_polZ_diff->SetLineWidth(2);
    // h1_polZ_diff->SetXTitle("Z Polarization diff. (Target - Initial)");
    // h1_polZ_diff->SetYTitle("Event/Entry");
    // h1_polZ_diff->Draw();
}
