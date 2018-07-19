#include "TFile.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TH1.h"
#include "TH2.h"
#include "TRandom.h"
#include "TClassTable.h"
#include "TSystem.h"
#include "TROOT.h"

void ana_test(){

  gROOT->Reset();

  ////////// DEFINE THE VARIABLES //////////
  double muIniMomX;
  double muIniPosX;
  double muIniMomY;
  double muIniPosY;
  double muIniMomZ;


  ////////// GET THE TREE FILE //////////
  // TFile *myfile = new TFile("/home/had/hyasuda/srsim/run/data/musr_330.root","READ"); // for KEKCC directory
  TFile *myfile = new TFile("/Users/YASUDA/data/muonLinac/srsim/musr_330.root","READ"); // for local directory
  // TTree *tree = myfile->Get("t1") ;
  TTree *tree = (TTree*)myfile->Get("t1") ;
  tree->SetBranchAddress("muIniMomX",&muIniMomX);
  tree->SetBranchAddress("muIniPosX",&muIniPosX);
  tree->SetBranchAddress("muIniMomY",&muIniMomY);
  tree->SetBranchAddress("muIniPosY",&muIniPosY);
  tree->SetBranchAddress("muIniMomZ",&muIniMomZ);

  int nEntry = tree->GetEntries();
  // nEntry = 1000; // for test
  // nEntry = 10; // for test

  /////////// CHECK THE FILE CONTENTS //////////
  // cout << "CHECK THE FILE CONTENTS" << endl;
  // cout << "muIniMomX" << "\t" << "muIniPosX" << "\t" << "muIniMomY" << "\t" << "muIniPosY" << "\t" << "muIniMomZ" << endl;
  for(int iEntry = 0 ; iEntry < nEntry ; iEntry++){
    tree->GetEntry(iEntry);
    // cout << muIniMomX << "\t" << muIniPosX << "\t" << muIniMomY << "\t" << muIniPosY << "\t" << muIniMomZ << endl; // for debug
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


  // TH2D * h2_psX = new TH2D();
  // TH2D * h2_psY = new TH2D();
  // TH2D (const char *name, const char *title, Int_t nbinsx, Double_t xlow, Double_t xup, Int_t nbinsy, Double_t ylow, Double_t yup) //
  TH2D * h2_pmX = new TH2D("h2_pmX", "h2_pmX", 60, -0.3, 0.3, 60, -0.3, 0.3); // res : 0.01
  TH2D * h2_pmY = new TH2D("h2_pmY", "h2_pmY", 60, -0.3, 0.3, 60, -0.3, 0.3); // res : 0.01
  TH2D * h2_psX = new TH2D("h2_psX", "h2_psX", 60, -0.3, 0.3, 60, -0.03, 0.03);
  TH2D * h2_psY = new TH2D("h2_psY", "h2_psY", 60, -0.3, 0.3, 60, -0.03, 0.03);


  // cout << "x" << "\t" << "xp" << "\t" << "y" << "\t" << "yp" << endl;  // for debug
  for(int iEntry = 0 ; iEntry < nEntry ; iEntry++){
    tree->GetEntry(iEntry);
    xp = muIniMomX/muIniMomZ;
    yp = muIniMomY/muIniMomZ;
    x  = muIniPosX;
    y  = muIniPosY;
    xm = muIniMomX;
    ym = muIniMomY;
    cout << xp << "\t" << yp << endl;
    cout << xm << "\t" << ym << endl;

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

  TCanvas *c1 = new TCanvas("c1","c1");
  TCanvas *c2 = new TCanvas("c2","c2");
  TCanvas *c3 = new TCanvas("c3","c3");
  TCanvas *c4 = new TCanvas("c4","c4");

  c1->cd();
  h2_psX->Draw("colz");

  c2->cd();
  h2_psY->Draw("colz");

  c3->cd();
  h2_pmX->Draw("colz");

  c4->cd();
  h2_pmY->Draw("colz");


  ////////// CALCULATION OF EMITTANCE //////////

  xp_sum = 0;
  yp_sum = 0;
  x_sum  = 0;
  y_sum  = 0;
  double xx_sum   = 0;
  double xpxp_sum = 0;
  double xxp_sum  = 0;
  double yy_sum   = 0;
  double ypyp_sum = 0;
  double yyp_sum  = 0;
  double zm = 0;
  double zm_sum = 0;
  // double zmzm_mean = 0



  for(int iEntry = 0 ; iEntry < nEntry ; iEntry++){
    tree->GetEntry(iEntry);
    xp = muIniMomX/muIniMomZ;
    yp = muIniMomY/muIniMomZ;
    zm = muIniMomZ;
    x  = muIniPosX;
    y  = muIniPosY;
    // xm = muIniMomX;
    // ym = muIniMomY;

    xx_sum += x*x;
    xpxp_sum += xp*xp;
    xxp_sum += x*xp;

    yy_sum += y*y;
    ypyp_sum += yp*yp;
    yyp_sum += y*yp;

    zm_sum += zm;
}

  cout << "Mean X^2" << "\t" << "Mean X'^2" << "\t" << "Mean Y^2" << "\t" << "Mean Y'^2" << endl;
  cout << xx_sum     << "\t" << xpxp_sum    << "\t" << yy_sum     << "\t" << ypyp_sum    << endl;
  cout << "Mean X*X'"<< "\t" << "Mean Y *Y'"<< endl;
  cout << xxp_sum    << "\t" << yyp_sum     << endl;
  cout << "Mean Pz" << endl;
  cout << zm_sum/nEntry << endl;

  double ex_rms = 0; // RMS Emittance of PhaseSpace X
  double ey_rms = 0; // RMS Emittance of PhaseSpace Y

  // ex_rms = sqrt(xx_sum*xpxp_sum/(nEntry*nEntry) - xxp_sum*xxp_sum/(nEntry*nEntry));
  // ey_rms = sqrt(yy_sum*ypyp_sum/(nEntry*nEntry) - yyp_sum*yyp_sum/(nEntry*nEntry));
  ex_rms = sqrt(xx_sum*xpxp_sum - xxp_sum*xxp_sum)/nEntry;
  ey_rms = sqrt(yy_sum*ypyp_sum - yyp_sum*yyp_sum)/nEntry;

  cout << "RMS Emitttance of X : " << ex_rms << endl;
  cout << "RMS Emitttance of Y : " << ey_rms << endl;

  double beta  = 0.0796704;
  double gamma = 1.00319;
  // cout << "Normalized RMS Emitttance of X : " << ex_rms*beta*gamma*3.1415*3.1415*1000 << endl;
  // cout << "Normalized RMS Emitttance of Y : " << ey_rms*beta*gamma*3.1415*3.1415*1000 << endl;
  cout << "Normalized RMS Emitttance of X : " << ex_rms*beta*gamma*3.1415*3.1415*1000 << endl;
  cout << "Normalized RMS Emitttance of Y : " << ey_rms*beta*gamma*3.1415*3.1415*1000 << endl;

}
