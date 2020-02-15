#include <TAxis.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <THn.h>
#include <THnSparse.h>
#include <TLine.h>
#include <TMath.h>
#include <TROOT.h>
#include <TStyle.h>
#include <iostream>
#include <string>
#include "../include/constants.hpp"

int make_plots(std::string inFileName =
                   "/home/krishna/Downloads/clas12_elastic/build/elastic_clas12_alldata_total_01_30_2020.root") {
  gStyle->SetHistMinimumZero();
  TFile *root_data = new TFile(inFileName.c_str());

  for (auto &&det : detector_name) {
    TCanvas *can = new TCanvas(Form("Phi_%s", det.second.c_str()), Form("Phi_%s", det.second.c_str()), 1920, 1080);
    can->Divide(3, 2);
    for (int sec = 1; sec < 7; sec++) {
      can->cd(sec);
      TH1D *phi = (TH1D *)root_data->Get(Form("Phi/Phie_minus_Phip_%s_%d", det.second.c_str(), sec));
      phi->GetXaxis()->SetRangeUser(2, 4);
      phi->Draw("");
      TLine *l = new TLine(phi_min_cut, 0, phi_min_cut, phi->GetMaximum() * 1.05);
      l->Draw("SAME");
      TLine *r = new TLine(phi_max_cut, 0, phi_max_cut, phi->GetMaximum() * 1.05);
      r->Draw("SAME");
    }
    can->SaveAs(Form("%s.png", can->GetName()));
  }

  {
    TCanvas *can = new TCanvas("MissingMass", "MissingMass", 1920, 1080);
    can->Divide(3, 2);
    for (int sec = 1; sec < 7; sec++) {
      can->cd(sec);
      TH1D *missingMass2 = (TH1D *)root_data->Get(Form("W_vs_Q2/MM2_hist_sec_%d", sec));
      missingMass2->GetXaxis()->SetRangeUser(-1, 1);
      missingMass2->Draw("");
      TLine *l = new TLine(-MM2_cut, 0, -MM2_cut, missingMass2->GetMaximum() * 1.05);
      l->Draw("SAME");
      TLine *r = new TLine(MM2_cut, 0, MM2_cut, missingMass2->GetMaximum() * 1.05);
      r->Draw("SAME");
    }
    can->SaveAs(Form("%s.png", can->GetName()));
  }

  for (auto &&det : detector_name) {
    TCanvas *can =
        new TCanvas(Form("W_at180_%s", det.second.c_str()), Form("W_at180_%s", det.second.c_str()), 1920, 1080);
    can->Divide(3, 2);
    for (int sec = 1; sec < 7; sec++) {
      can->cd(sec);
      TH1D *W = (TH1D *)root_data->Get(Form("at180/W_1pos_at180_%s_%d", det.second.c_str(), sec));
      W->Draw("");
    }
    can->SaveAs(Form("%s.png", can->GetName()));
  }

  for (auto &&det : detector_name) {
    TCanvas *can =
        new TCanvas(Form("W_at180_MM_%s", det.second.c_str()), Form("W_at180_MM_%s", det.second.c_str()), 1920, 1080);
    can->Divide(3, 2);
    for (int sec = 1; sec < 7; sec++) {
      can->cd(sec);
      TH1D *W = (TH1D *)root_data->Get(Form("at180_MM/W_1pos_at180_MM_%s_%d", det.second.c_str(), sec));
      W->Draw("");
    }
    can->SaveAs(Form("%s.png", can->GetName()));
  }

  for (auto &&det : detector_name) {
    TCanvas *can = new TCanvas(Form("W_at180_MMvsWo_%s", det.second.c_str()),
                               Form("W_at180_MMvsWo_%s", det.second.c_str()), 1920, 1080);
    can->Divide(3, 2);
    for (int sec = 1; sec < 7; sec++) {
      can->cd(sec);
      TH1D *W_wo = (TH1D *)root_data->Get(Form("at180/W_1pos_at180_%s_%d", det.second.c_str(), sec));
      W_wo->SetLineColor(kRed);
      W_wo->Draw("same");

      TH1D *W = (TH1D *)root_data->Get(Form("at180_MM/W_1pos_at180_MM_%s_%d", det.second.c_str(), sec));
      W->Draw("same");
    }
    can->SaveAs(Form("%s.png", can->GetName()));
  }

  {
    TCanvas *can = new TCanvas("MomVsTheta", "MomVsTheta", 1920, 1080);
    can->Divide(2, 2);
    can->cd(1);
    TH2D *MomVsTheta = (TH2D *)root_data->Get("Mom Vs Beta/MomVsTheta_pos_both_0");
    MomVsTheta->Draw("");
    can->cd(2);
    TH2D *MomVsTheta_lowW = (TH2D *)root_data->Get("Mom Vs Beta/MomVsTheta_lowW_both_0");
    MomVsTheta_lowW->Draw("");
    can->cd(3);
    TH2D *MomVsTheta_cent = (TH2D *)root_data->Get("Mom Vs Beta/MomVsTheta_lowW_in_Central_0");
    MomVsTheta_cent->Draw("");
    can->cd(4);
    TH2D *MomVsTheta_for = (TH2D *)root_data->Get("Mom Vs Beta/MomVsTheta_lowW_in_Forward_0");
    MomVsTheta_for->Draw("");
    can->SaveAs(Form("%s.png", can->GetName()));
  }
  /*
    {
      TCanvas *can = new TCanvas("W in Q^2 bins", "W in Q^2 bins", 1920, 1080);
      // can->DivideSquare(10);
      THnSparseD *nsparce = (THnSparseD *)root_data->Get("nsparce");
      can->cd(0);
      nsparce->Projection(1, 0)->Draw("COLZ");
      nsparce->Projection(1, 0)->SetShowProjectionX(true);
    }
  */
  return 0;
}
