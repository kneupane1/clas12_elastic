/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/
#include "histogram.hpp"

Histogram::Histogram(const std::string& output_file) {
  RootOutputFile = std::make_shared<TFile>(output_file.c_str(), "RECREATE");
  def = std::make_shared<TCanvas>("def");

  makeHists();
  Nsparce = std::make_shared<THnSparseD>("nsparce", "nsparce", 3, sparce_bins, sparce_xmin, sparce_xmax);
  makeHists_electron_cuts();
  makeHistSF();
}

Histogram::~Histogram() { this->Write(); }

void Histogram::Write() {
  std::cout << GREEN << "Writting" << DEF << std::endl;
  Write_SF();
  // Nsparce->Sumw2();
  Nsparce->Write();
  std::cout << BOLDBLUE << "WvsQ2()" << DEF << std::endl;
  Write_WvsQ2();

  std::cout << BOLDBLUE << "Write_MomVsBeta()" << DEF << std::endl;
  TDirectory* Write_MomVsBeta_folder = RootOutputFile->mkdir("Mom Vs Beta");
  Write_MomVsBeta_folder->cd();
  Write_MomVsBeta();
  deltaT_proton[0]->Write();
  deltaT_proton[1]->Write();

  std::cerr << BOLDBLUE << "Write_Electron_cuts()" << DEF << std::endl;
  TDirectory* Electron_Cuts = RootOutputFile->mkdir("Electron_Cuts");
  Electron_Cuts->cd();
  Write_Electron_cuts();

  std::cerr << BOLDBLUE << "Write_SF_1D()" << DEF << std::endl;
  TDirectory* SF_1D = RootOutputFile->mkdir("sampling_fraction");
  SF_1D->cd();
  write_histSf();

  std::cout << BOLDBLUE << "Done Writing!!!" << DEF << std::endl;
}

void Histogram::makeHists() {
  deltaT_proton[0] = std::make_shared<TH2D>("DeltaTProton", "DeltaTProton", bins, zero, 10, bins, -5, 5);
  deltaT_proton[1] = std::make_shared<TH2D>("DeltaTProton_cut", "DeltaTProton_cut", bins, zero, 10, bins, -5, 5);
  Int_t n = 10;

  // for (Int_t i = 0; i < n; i++) {
  //   x[i] = i * 0.1;
  //   y[i] = 10 * sin(x[i] + 0.2);
  // }

  Double_t P_e[10] = {0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 10};
  Double_t SF_upper[10] = {.25 + 0.0,       .236 + 3 * .021, .244 + 3 * .019, .249 + 3 * .016, .254 + 3 * .012,
                           .252 + 3 * .010, .252 + 3 * .009, .252 + 3 * .008, .25 + 3 * .008,  .237 + 3 * .011};
  Double_t SF_lower[10] = {.25 - 0.0,       .236 - 3 * .021, .244 - 3 * .019, .249 - 3 * .016, .254 - 3 * .012,
                           .252 - 3 * .010, .252 - 3 * .009, .252 - 3 * .008, .25 - 3 * .008,  .237 - 3 * .011};
  SF_gr_upper = std::make_shared<TGraph>(n, P_e, SF_upper);
  SF_gr_lower = std::make_shared<TGraph>(n, P_e, SF_lower);

  for (short sec = 0; sec < num_sectors; sec++) {
    MissingMass[sec] =
        std::make_shared<TH1D>(Form("MM2_hist_sec_%d", sec), Form("MM2_hist_sec_%d", sec), bins, -w_max, w_max);

    mass_pi0_hist[0][sec] =
        std::make_shared<TH1D>(Form("mass_pi0_hist_%d", sec), Form("mass_pi0_hist_%d", sec), bins, 0, 0.5);
    mass_pi0_hist[1][sec] = std::make_shared<TH1D>(Form("mass_pi0_hist_aferPcuts_%d", sec),
                                                   Form("mass_pi0_hist_aferPcuts_%d", sec), bins, 0, 0.5);

    W_hist_all_events[sec] =
        std::make_shared<TH1D>(Form("W_hist_sec_%d", sec), Form("W_hist_sec_%d", sec), bins, zero, w_max);
    W_hist_1pos[sec] =
        std::make_shared<TH1D>(Form("W_hist_1pos_%d", sec), Form("W_hist_1pos_%d", sec), bins, zero, w_max);
    W_hist_1pos_0charge[sec] =
        std::make_shared<TH1D>(Form("W_hist_1pos_MM0_%d", sec), Form("W_hist_1pos_MM0_%d", sec), bins, zero, w_max);
    W_hist_1pos_noOther[sec] = std::make_shared<TH1D>(Form("W_hist_1pos_noOther_%d", sec),
                                                      Form("W_hist_1pos_noOther_%d", sec), bins, zero, w_max);

    W_vs_q2_all_events[sec] =
        std::make_shared<TH2D>(Form("WQ2_sec_%d", sec), Form("WQ2_sec_%d", sec), bins, zero, w_max, bins, zero, q2_max);
    W_vs_q2_1pos[sec] = std::make_shared<TH2D>(Form("WQ2_1pos_%d", sec), Form("WQ2_1pos_%d", sec), bins, zero, w_max,
                                               bins, zero, q2_max);
    W_vs_q2_1pos_0charge[sec] = std::make_shared<TH2D>(Form("WQ2_1pos_MM0_%d", sec), Form("WQ2_1pos_MM0_%d", sec), bins,
                                                       zero, w_max, bins, zero, q2_max);
    W_vs_q2_1pos_noOther[sec] = std::make_shared<TH2D>(
        Form("WQ2_1pos_noOther_%d", sec), Form("WQ2_1pos_noOther_%d", sec), bins, zero, w_max, bins, zero, q2_max);
    for (auto&& det : detector_name) {
      int d = detector_fill[det.first];

      ThetaVsP[d][sec] =
          std::make_shared<TH2D>(Form("MomVsTheta_pos_%s_%d", det.second.c_str(), sec),
                                 Form("MomVsTheta_pos_%s_%d", det.second.c_str(), sec), 500, zero, 6.0, 500, 0, 90);

      ThetaVsP_lowW[d][sec] =
          std::make_shared<TH2D>(Form("MomVsTheta_lowW_%s_%d", det.second.c_str(), sec),
                                 Form("MomVsTheta_lowW_%s_%d", det.second.c_str(), sec), 500, zero, 6.0, 500, 0, 90);

      ThetaVsPCalc[d][sec] =
          std::make_shared<TH2D>(Form("MomVsTheta_Calc_%s_%d", det.second.c_str(), sec),
                                 Form("MomVsTheta_Calc_%s_%d", det.second.c_str(), sec), 500, zero, 6.0, 500, 0, 90);
      MomVsBeta[d][sec] =
          std::make_shared<TH2D>(Form("MomVsBeta_%s_%d", det.second.c_str(), sec),
                                 Form("MomVsBeta_%s_%d", det.second.c_str(), sec), 500, zero, p_max, 500, zero, 1.2);
      Phie_vs_Phip[d][sec] =
          std::make_shared<TH2D>(Form("Phie_vs_Phip_%s_%d", det.second.c_str(), sec),
                                 Form("Phie_vs_Phip_%s_%d", det.second.c_str(), sec), 500, -PI, PI, 500, -PI, PI);
      Phie_Phip_hist[d][sec] =
          std::make_shared<TH1D>(Form("Phie_minus_Phip_%s_%d", det.second.c_str(), sec),
                                 Form("Phie_minus_Phip_%s_%d", det.second.c_str(), sec), 500, zero, 2 * PI);
      W_hist_1pos_at180[d][sec] =
          std::make_shared<TH1D>(Form("W_1pos_at180_%s_%d", det.second.c_str(), sec),
                                 Form("W_1pos_at180_%s_%d", det.second.c_str(), sec), bins, zero, w_max);
      W_vs_q2_1pos_at180[d][sec] = std::make_shared<TH2D>(Form("WQ2_1pos_at180_%s_%d", det.second.c_str(), sec),
                                                          Form("WQ2_1pos_at180_%s_%d", det.second.c_str(), sec), bins,
                                                          zero, w_max, bins, zero, q2_max);
      W_hist_1pos_at180_MM[d][sec] =
          std::make_shared<TH1D>(Form("W_1pos_at180_MM_%s_%d", det.second.c_str(), sec),
                                 Form("W_1pos_at180_MM_%s_%d", det.second.c_str(), sec), bins, zero, w_max);
      W_vs_q2_1pos_at180_MM[d][sec] = std::make_shared<TH2D>(Form("WQ2_1pos_at180_MM_%s_%d", det.second.c_str(), sec),
                                                             Form("WQ2_1pos_at180_MM_%s_%d", det.second.c_str(), sec),
                                                             bins, zero, w_max, bins, zero, q2_max);
    }
  }
}

void Histogram::makeHists_electron_cuts() {
  for (short c = 0; c < num_cuts; c++) {
    EC_sampling_fraction[c] = std::make_shared<TH2D>(Form("EC_sampling_fraction%1.12s", cut_name[c].c_str()),
                                                     Form("EC_sampling_fraction%1.12s", cut_name[c].c_str()), bins,
                                                     p_min, p_max, bins, zero, 1.0);
    vz_position[c] = std::make_shared<TH1D>(Form("vz_position%1.12s", cut_name[c].c_str()),
                                            Form("vz_position%1.12s", cut_name[c].c_str()), bins, -40, 40);
    pcal_sec[c] = std::make_shared<TH2D>(Form("pcal_sec%1.12s", cut_name[c].c_str()),
                                         Form("pcal_sec%1.12s", cut_name[c].c_str()), bins, -420, 420, bins, -420, 420);
    dcr1_sec[c] = std::make_shared<TH2D>(Form("dcr1_sec%1.12s", cut_name[c].c_str()),
                                         Form("dcr1_sec%1.12s", cut_name[c].c_str()), bins, -180, 180, bins, -180, 180);
    dcr2_sec[c] = std::make_shared<TH2D>(Form("dcr2_sec%1.12s", cut_name[c].c_str()),
                                         Form("dcr2_sec%1.12s", cut_name[c].c_str()), bins, -270, 270, bins, -270, 270);
    dcr3_sec[c] = std::make_shared<TH2D>(Form("dcr3_sec%1.12s", cut_name[c].c_str()),
                                         Form("dcr3_sec%1.12s", cut_name[c].c_str()), bins, -320, 320, bins, -320, 320);
  }
}

void Histogram::makeHistSF() {
  for (size_t i = 0; i < W_BINS; i++) {
    SF_1D[i] = std::make_shared<TH1D>(Form("sf_1d_%1.12s_GeV", W_BINS_NAME[i].c_str()),
                                      Form("sf_1d%1.12s_GeV", W_BINS_NAME[i].c_str()), bins, 0.0, 0.5);
  }
}
void Histogram::Fill_SF(const std::shared_ptr<Branches12>& _d) {
  sf_hist->Fill(_d->p(0), _d->ec_tot_energy(0) / _d->p(0));

  if (_d->p(0) > 0 && _d->p(0) < 1)
    SF_1D[0]->Fill(_d->ec_tot_energy(0) / _d->p(0));
  else if (_d->p(0) > 1. && _d->p(0) < 2.)
    SF_1D[1]->Fill(_d->ec_tot_energy(0) / _d->p(0));
  else if (_d->p(0) > 2 && _d->p(0) < 3)
    SF_1D[2]->Fill(_d->ec_tot_energy(0) / _d->p(0));
  else if (_d->p(0) > 3 && _d->p(0) < 4)
    SF_1D[3]->Fill(_d->ec_tot_energy(0) / _d->p(0));
  else if (_d->p(0) > 4 && _d->p(0) < 5)
    SF_1D[4]->Fill(_d->ec_tot_energy(0) / _d->p(0));
  else if (_d->p(0) > 5 && _d->p(0) < 6)
    SF_1D[5]->Fill(_d->ec_tot_energy(0) / _d->p(0));
  else if (_d->p(0) > 6 && _d->p(0) < 7)
    SF_1D[6]->Fill(_d->ec_tot_energy(0) / _d->p(0));
  else if (_d->p(0) > 7 && _d->p(0) < 8)
    SF_1D[7]->Fill(_d->ec_tot_energy(0) / _d->p(0));
  else if (_d->p(0) > 8 && _d->p(0) < 9)
    SF_1D[8]->Fill(_d->ec_tot_energy(0) / _d->p(0));
  else if (_d->p(0) > 9 && _d->p(0) < 12)
    SF_1D[9]->Fill(_d->ec_tot_energy(0) / _d->p(0));
}
void Histogram::FillHists_electron_cuts(const std::shared_ptr<Branches12>& _d) {
  vz_position[0]->Fill(_d->vz(0));
  pcal_sec[0]->Fill(_d->ec_pcal_x(0), _d->ec_pcal_y(0));
  dcr1_sec[0]->Fill(_d->dc_r1_x(0), _d->dc_r1_y(0));
  dcr2_sec[0]->Fill(_d->dc_r2_x(0), _d->dc_r2_y(0));
  dcr3_sec[0]->Fill(_d->dc_r3_x(0), _d->dc_r3_y(0));
  EC_sampling_fraction[0]->Fill(_d->p(0), _d->ec_tot_energy(0) / _d->p(0));
}

void Histogram::FillHists_electron_with_cuts(const std::shared_ptr<Branches12>& _d) {
  vz_position[1]->Fill(_d->vz(0));
  pcal_sec[1]->Fill(_d->ec_pcal_x(0), _d->ec_pcal_y(0));
  dcr1_sec[1]->Fill(_d->dc_r1_x(0), _d->dc_r1_y(0));
  dcr2_sec[1]->Fill(_d->dc_r2_x(0), _d->dc_r2_y(0));
  dcr3_sec[1]->Fill(_d->dc_r3_x(0), _d->dc_r3_y(0));
  EC_sampling_fraction[1]->Fill(_d->p(0), _d->ec_tot_energy(0) / _d->p(0));
}
void Histogram::Write_SF() {
  sf_hist->Write();
  // gStyle->SetOptFit(1111);
  gStyle->SetOptFit(1111);
  SF_gr_upper->Fit("pol2");
  SF_gr_upper->Write();
  SF_gr_lower->Fit("pol2");
  SF_gr_lower->Write();
}
void Histogram::write_histSf() {
  for (size_t i = 0; i < W_BINS; i++) {
    SF_1D[i]->SetXTitle("sf (E/P)");
    SF_1D[i]->Fit("gaus", "QMR+", "QMR+", 0.18, 0.30);
    //  gROOT->SetStyle("Plain");
    gStyle->SetOptFit(1111);
    // if (SF_1D[i]->GetEntries())
    SF_1D[i]->Write();
  }
}
void Histogram::Write_Electron_cuts() {
  for (short c = 0; c < num_cuts; c++) {
    vz_position[c]->SetXTitle("vz (GeV)");
    if (vz_position[c]->GetEntries()) vz_position[c]->Write();
    pcal_sec[c]->SetXTitle("x/cm");
    pcal_sec[c]->SetYTitle("y/cm");
    pcal_sec[c]->SetOption("COLZ1");
    if (pcal_sec[c]->GetEntries()) pcal_sec[c]->Write();

    dcr1_sec[c]->SetXTitle("x/cm");
    dcr1_sec[c]->SetYTitle("y/cm");
    dcr1_sec[c]->SetOption("COLZ1");
    if (dcr1_sec[c]->GetEntries()) dcr1_sec[c]->Write();

    dcr2_sec[c]->SetXTitle("x/cm");
    dcr2_sec[c]->SetYTitle("y/cm");
    dcr2_sec[c]->SetOption("COLZ1");
    if (dcr2_sec[c]->GetEntries()) dcr2_sec[c]->Write();

    dcr3_sec[c]->SetXTitle("x/cm");
    dcr3_sec[c]->SetYTitle("y/cm");
    dcr3_sec[c]->SetOption("COLZ1");
    if (dcr3_sec[c]->GetEntries()) dcr3_sec[c]->Write();

    EC_sampling_fraction[c]->SetXTitle("Momentum (GeV)");
    EC_sampling_fraction[c]->SetYTitle("Sampling Fraction");
    EC_sampling_fraction[c]->SetOption("COLZ1");
    EC_sampling_fraction[c]->Write();
  }
}
void Histogram::Fill_Sparce(const std::shared_ptr<Reaction>& _e) {
  std::lock_guard<std::mutex> lk(mutex);
  double ret[NUM_DIM] = {_e->W(), _e->Q2(), static_cast<double>(_e->sec())};
  Nsparce->Fill(ret);
}
void Histogram::Fill_WvsQ2(const std::shared_ptr<Reaction>& _e) {
  short sec = _e->sec();
  short pos_det = _e->pos_det();
  if ((sec > 0 && sec < num_sectors) || pos_det != -1) {
    W_hist_all_events[0]->Fill(_e->W());
    W_vs_q2_all_events[0]->Fill(_e->W(), _e->Q2());
    W_hist_all_events[sec]->Fill(_e->W());
    W_vs_q2_all_events[sec]->Fill(_e->W(), _e->Q2());

    if (_e->onePositive()) {
      W_hist_1pos[0]->Fill(_e->W());
      W_vs_q2_1pos[0]->Fill(_e->W(), _e->Q2());
      W_hist_1pos[sec]->Fill(_e->W());
      W_vs_q2_1pos[sec]->Fill(_e->W(), _e->Q2());

      Phie_vs_Phip[0][0]->Fill(_e->phi_e(), _e->phi_p());
      Phie_Phip_hist[0][0]->Fill(_e->phi_diff());
      Phie_vs_Phip[0][sec]->Fill(_e->phi_e(), _e->phi_p());
      Phie_Phip_hist[0][sec]->Fill(_e->phi_diff());

      Phie_vs_Phip[pos_det][0]->Fill(_e->phi_e(), _e->phi_p());
      Phie_Phip_hist[pos_det][0]->Fill(_e->phi_diff());
      Phie_vs_Phip[pos_det][sec]->Fill(_e->phi_e(), _e->phi_p());
      Phie_Phip_hist[pos_det][sec]->Fill(_e->phi_diff());
    }
    if (_e->onePositive_MM0()) {
      W_hist_1pos_0charge[0]->Fill(_e->W());
      W_vs_q2_1pos_0charge[0]->Fill(_e->W(), _e->Q2());
      W_hist_1pos_0charge[sec]->Fill(_e->W());
      W_vs_q2_1pos_0charge[sec]->Fill(_e->W(), _e->Q2());
    }
    if (_e->onePositive_noOther()) {
      W_hist_1pos_noOther[0]->Fill(_e->W());
      W_vs_q2_1pos_noOther[0]->Fill(_e->W(), _e->Q2());
      W_hist_1pos_noOther[sec]->Fill(_e->W());
      W_vs_q2_1pos_noOther[sec]->Fill(_e->W(), _e->Q2());
    }
    if (_e->onePositive_at180()) {
      MissingMass[0]->Fill(_e->MM2());
      MissingMass[sec]->Fill(_e->MM2());

      W_hist_1pos_at180[0][0]->Fill(_e->W());
      W_vs_q2_1pos_at180[0][0]->Fill(_e->W(), _e->Q2());
      W_hist_1pos_at180[0][sec]->Fill(_e->W());
      W_vs_q2_1pos_at180[0][sec]->Fill(_e->W(), _e->Q2());

      W_hist_1pos_at180[pos_det][0]->Fill(_e->W());
      W_vs_q2_1pos_at180[pos_det][0]->Fill(_e->W(), _e->Q2());
      W_hist_1pos_at180[pos_det][sec]->Fill(_e->W());
      W_vs_q2_1pos_at180[pos_det][sec]->Fill(_e->W(), _e->Q2());
    }
    if (_e->onePositive_at180_MM0()) {
      W_hist_1pos_at180_MM[0][0]->Fill(_e->W());
      W_vs_q2_1pos_at180_MM[0][0]->Fill(_e->W(), _e->Q2());
      W_hist_1pos_at180_MM[0][sec]->Fill(_e->W());
      W_vs_q2_1pos_at180_MM[0][sec]->Fill(_e->W(), _e->Q2());

      W_hist_1pos_at180_MM[pos_det][0]->Fill(_e->W());
      W_vs_q2_1pos_at180_MM[pos_det][0]->Fill(_e->W(), _e->Q2());
      W_hist_1pos_at180_MM[pos_det][sec]->Fill(_e->W());
      W_vs_q2_1pos_at180_MM[pos_det][sec]->Fill(_e->W(), _e->Q2());
    }
  }
}

void Histogram::Write_WvsQ2() {
  TDirectory* phi_folder = RootOutputFile->mkdir("Phi");
  phi_folder->cd();
  for (short j = 0; j < detector_name.size(); j++) {
    for (int i = 0; i < num_sectors; i++) {
      Phie_vs_Phip[j][i]->SetXTitle("Phie");
      Phie_vs_Phip[j][i]->SetYTitle("Phip");
      Phie_vs_Phip[j][i]->SetOption("COLZ");
      Phie_vs_Phip[j][i]->Write();
    }
    for (int i = 0; i < num_sectors; i++) {
      Phie_Phip_hist[j][i]->SetXTitle("Phi");
      Phie_Phip_hist[j][i]->Write();
    }
  }

  TDirectory* at180_folder = RootOutputFile->mkdir("at180");
  at180_folder->cd();
  for (short j = 0; j < detector_name.size(); j++) {
    for (int i = 0; i < num_sectors; i++) {
      W_hist_1pos_at180[j][i]->SetXTitle("W (GeV)");
      W_hist_1pos_at180[j][i]->Write();
    }
    for (int i = 0; i < num_sectors; i++) {
      W_vs_q2_1pos_at180[j][i]->SetXTitle("W (GeV)");
      W_vs_q2_1pos_at180[j][i]->SetYTitle("Q^2 (GeV^2)");
      W_vs_q2_1pos_at180[j][i]->SetOption("COLZ");
      W_vs_q2_1pos_at180[j][i]->Write();
    }
  }

  TDirectory* at180_MM_folder = RootOutputFile->mkdir("at180_MM");
  at180_MM_folder->cd();
  for (short j = 0; j < detector_name.size(); j++) {
    for (int i = 0; i < num_sectors; i++) {
      W_hist_1pos_at180_MM[j][i]->SetXTitle("W (GeV)");
      W_hist_1pos_at180_MM[j][i]->Write();
    }
    for (int i = 0; i < num_sectors; i++) {
      W_vs_q2_1pos_at180_MM[j][i]->SetXTitle("W (GeV)");
      W_vs_q2_1pos_at180_MM[j][i]->SetYTitle("Q^2 (GeV^2)");
      W_vs_q2_1pos_at180_MM[j][i]->SetOption("COLZ");
      W_vs_q2_1pos_at180_MM[j][i]->Write();
    }
  }

  TDirectory* W_vs_Q2_folder = RootOutputFile->mkdir("W_vs_Q2");
  W_vs_Q2_folder->cd();
  for (int i = 0; i < num_sectors; i++) {
    MissingMass[i]->SetXTitle("MM^2 (GeV)");
    MissingMass[i]->Write();

    mass_pi0_hist[0][i]->SetXTitle("MM(GeV)");
    mass_pi0_hist[0][i]->Write();

    mass_pi0_hist[1][i]->SetXTitle("MM(GeV)");
    mass_pi0_hist[1][i]->Write();

    W_hist_all_events[i]->SetXTitle("W (GeV)");
    W_hist_all_events[i]->Write();
    W_hist_1pos[i]->SetXTitle("W (GeV)");
    W_hist_1pos[i]->Write();
    W_hist_1pos_0charge[i]->SetXTitle("W (GeV)");
    W_hist_1pos_0charge[i]->Write();
    W_hist_1pos_noOther[i]->SetXTitle("W (GeV)");
    W_hist_1pos_noOther[i]->Write();

    W_vs_q2_all_events[i]->SetXTitle("W (GeV)");
    W_vs_q2_all_events[i]->SetYTitle("Q^2 (GeV^2)");
    W_vs_q2_all_events[i]->SetOption("COLZ");
    W_vs_q2_all_events[i]->Write();

    W_vs_q2_1pos[i]->SetXTitle("W (GeV)");
    W_vs_q2_1pos[i]->SetYTitle("Q^2 (GeV^2)");
    W_vs_q2_1pos[i]->SetOption("COLZ");
    W_vs_q2_1pos[i]->Write();

    W_vs_q2_1pos_0charge[i]->SetXTitle("W (GeV)");
    W_vs_q2_1pos_0charge[i]->SetYTitle("Q^2 (GeV^2)");
    W_vs_q2_1pos_0charge[i]->SetOption("COLZ");
    W_vs_q2_1pos_0charge[i]->Write();

    W_vs_q2_1pos_noOther[i]->SetXTitle("W (GeV)");
    W_vs_q2_1pos_noOther[i]->SetYTitle("Q^2 (GeV^2)");
    W_vs_q2_1pos_noOther[i]->SetOption("COLZ");
    W_vs_q2_1pos_noOther[i]->Write();
  }
}

void Histogram::Fill_MomVsBeta(const std::shared_ptr<Reaction>& _e) {
  if (!_e->onePositive_at180()) return;
  if (_e->pos_det() == -1) return;

  MomVsBeta[0][0]->Fill(_e->pos_P(), _e->pos_beta());
  MomVsBeta[0][_e->sec()]->Fill(_e->pos_P(), _e->pos_beta());

  ThetaVsP[0][0]->Fill(_e->pos_P(), _e->pos_theta());
  ThetaVsP[0][_e->sec()]->Fill(_e->pos_P(), _e->pos_theta());

  ThetaVsPCalc[0][0]->Fill(_e->pos_P(), _e->pos_theta_calc());
  ThetaVsPCalc[0][_e->sec()]->Fill(_e->pos_P(), _e->pos_theta_calc());

  MomVsBeta[_e->pos_det()][0]->Fill(_e->pos_P(), _e->pos_beta());
  MomVsBeta[_e->pos_det()][_e->sec()]->Fill(_e->pos_P(), _e->pos_beta());

  ThetaVsP[_e->pos_det()][0]->Fill(_e->pos_P(), _e->pos_theta());
  ThetaVsP[_e->pos_det()][_e->sec()]->Fill(_e->pos_P(), _e->pos_theta());

  ThetaVsPCalc[_e->pos_det()][0]->Fill(_e->pos_P(), _e->pos_theta_calc());
  ThetaVsPCalc[_e->pos_det()][_e->sec()]->Fill(_e->pos_P(), _e->pos_theta_calc());

  if (_e->W() < 2.0) {
    ThetaVsP_lowW[0][0]->Fill(_e->pos_P(), _e->pos_theta());
    ThetaVsP_lowW[0][_e->sec()]->Fill(_e->pos_P(), _e->pos_theta());
    ThetaVsP_lowW[_e->pos_det()][0]->Fill(_e->pos_P(), _e->pos_theta());
    ThetaVsP_lowW[_e->pos_det()][_e->sec()]->Fill(_e->pos_P(), _e->pos_theta());
  }
}

void Histogram::Write_MomVsBeta() {
  for (short i = 0; i < detector_name.size(); i++) {
    for (short p = 0; p < num_sectors; p++) {
      MomVsBeta[i][p]->SetXTitle("Momentum (GeV)");
      MomVsBeta[i][p]->SetYTitle("#beta");
      MomVsBeta[i][p]->SetOption("COLZ1");
      MomVsBeta[i][p]->Write();
    }
    for (short p = 0; p < num_sectors; p++) {
      ThetaVsP[i][p]->SetXTitle("Momentum (GeV)");
      ThetaVsP[i][p]->SetYTitle("#theta");
      ThetaVsP[i][p]->SetOption("COLZ1");
      ThetaVsP[i][p]->Write();
    }
    for (short p = 0; p < num_sectors; p++) {
      ThetaVsP_lowW[i][p]->SetXTitle("Momentum (GeV)");
      ThetaVsP_lowW[i][p]->SetYTitle("#theta");
      ThetaVsP_lowW[i][p]->SetOption("COLZ1");
      ThetaVsP_lowW[i][p]->Write();
    }
    for (short p = 0; p < num_sectors; p++) {
      ThetaVsPCalc[i][p]->SetXTitle("Momentum (GeV)");
      ThetaVsPCalc[i][p]->SetYTitle("#theta");
      ThetaVsPCalc[i][p]->SetOption("COLZ1");
      ThetaVsPCalc[i][p]->Write();
    }
  }
}

void Histogram::Fill_Dt(const std::shared_ptr<Delta_T>& dt) {
  for (int i = 0; i < dt->gpart(); i++) {
    if (dt->charge(i) == 1) deltaT_proton[0]->Fill(dt->mom(i), dt->dt_P(i));
  }
}

void Histogram::Fill_Dt(const std::shared_ptr<Delta_T>& dt, int part) {
  deltaT_proton[1]->Fill(dt->mom(part), dt->dt_P(part));
}

void Histogram::Fill_pi0(const std::shared_ptr<Reaction>& _e) {
  if (_e->pi0_mass() < 0.0001) return;
  short sec = _e->sec();
  short pos_det = _e->pos_det();
  mass_pi0_hist[0][0]->Fill(_e->pi0_mass());
  if ((sec > 0 && sec < num_sectors) || pos_det != -1) {
    mass_pi0_hist[0][sec]->Fill(_e->pi0_mass());
    if (_e->onePositive_at180_MM0()) {
      mass_pi0_hist[1][0]->Fill(_e->pi0_mass());
      mass_pi0_hist[1][sec]->Fill(_e->pi0_mass());
    }
  }
}
