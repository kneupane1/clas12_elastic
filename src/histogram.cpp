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
  // makeHistSF();
  //  makeHists_x_mu();
}

Histogram::~Histogram() { this->Write(); }

void Histogram::Write() {
  std::cout << GREEN << "Writting" << DEF << std::endl;
  //  Write_SF();
  Nsparce->Write();
  std::cout << BOLDBLUE << "WvsQ2()" << DEF << std::endl;
  Write_WvsQ2();

  std::cout << BOLDBLUE << "Write_MomVsBeta()" << DEF << std::endl;
  TDirectory* Write_MomVsBeta_folder = RootOutputFile->mkdir("Mom Vs Beta");
  Write_MomVsBeta_folder->cd();
  Write_MomVsBeta();
  deltaT_proton[before_cut]->Write();
  deltaT_proton[after_cut]->Write();

  std::cerr << BOLDBLUE << "Write_Electron_cuts()" << DEF << std::endl;
  TDirectory* Electron_Cuts = RootOutputFile->mkdir("Electron_Cuts");
  Electron_Cuts->cd();
  Write_Electron_cuts();

  // std::cerr << BOLDBLUE << "write_histSf()" << DEF << std::endl;
  // TDirectory* histSf = RootOutputFile->mkdir("sampling_fraction");
  // histSf->cd();
  // write_histSf();
  // std::cerr << BOLDBLUE << "write_hist_x_mu()" << DEF << std::endl;
  // TDirectory* hists_x_mu = RootOutputFile->mkdir("hists_x_mu");
  // hists_x_mu->cd();
  // write_hist_x_mu();
  std::cout << BOLDBLUE << "Done Writing!!!" << DEF << std::endl;
}

void Histogram::makeHists() {
  deltaT_proton[before_cut] = std::make_shared<TH2D>("DeltaTProton", "DeltaTProton", bins, zero, 10, bins, -5, 5);
  deltaT_proton[after_cut] =
      std::make_shared<TH2D>("DeltaTProton_cut", "DeltaTProton_cut", bins, zero, 10, bins, -5, 5);

  // Double_t P_e[30] = {/*0.15, 0.45, 0.75, 1.05,*/
  //                     1.35, 1.65, 1.95, 2.25, 2.55, 2.85, 3.15, 3.45, 3.75, 4.05, 4.35, 4.65, 4.95, 5.25, 5.55,
  //                     5.85, 6.15, 6.45, 6.75, 7.05, 7.35, 7.65, 7.95, 8.25, 8.55, 8.85, 9.15, 9.45, 9.75, 10.05};
  //
  // Double_t SF_upper[30] = {/*.249 + (3 * .032), .197 + (3 * .042), .217 + (3 * .036), .230 + (3 * .035), */
  //                          .248 + (3 * .022), .241 + (3 * .022), .242 + (3 * .020), .244 + (3 * .019),
  //                          .245 + (3 * .018), .246 + (3 * .018), .246 + (3 * .017), .247 + (3 * .017),
  //                          .247 + (3 * .016), .247 + (3 * .016), .247 + (3 * .016), .247 + (3 * .015),
  //                          .247 + (3 * .015), .247 + (3 * .015), .247 + (3 * .015), .246 + (3 * .014),
  //                          .246 + (3 * .014), .245 + (3 * .014), .244 + (3 * .014), .243 + (3 * .014),
  //                          .242 + (3 * .014), .241 + (3 * .014), .240 + (3 * .014), .238 + (3 * .014),
  //                          .236 + (3 * .014), .234 + (3 * .014), .232 + (3 * .014), .230 + (3 * .014),
  //                          .217 + (3 * .015), .210 + (3 * .014) /*,.112+(3*0.027)*/};
  //
  // Double_t SF_lower[30] = {/*.249 - (3 * .032), .197 - (3 * .042), .217 - (3 * .036), .230 - (3 * .035), */
  //                          .248 - (3 * .022), .241 - (3 * .022), .242 - (3 * .020), .244 - (3 * .019),
  //                          .245 - (3 * .018), .246 - (3 * .018), .246 - (3 * .017), .247 - (3 * .017),
  //                          .247 - (3 * .016), .247 - (3 * .016), .247 - (3 * .016), .247 - (3 * .015),
  //                          .247 - (3 * .015), .247 - (3 * .015), .247 - (3 * .015), .246 - (3 * .014),
  //                          .246 - (3 * .014), .245 - (3 * .014), .244 - (3 * .014), .243 - (3 * .014),
  //                          .242 - (3 * .014), .241 - (3 * .014), .240 - (3 * .014), .238 - (3 * .014),
  //                          .236 - (3 * .014), .234 - (3 * .014), .232 - (3 * .014), .230 - (3 * .014),
  //                          .217 - (3 * .015), .210 - (3 * .014) /*,.112-(3*0.027)*/};
  // SF_gr_upper = std::make_shared<TGraph>(30, P_e, SF_upper);
  // SF_gr_lower = std::make_shared<TGraph>(30, P_e, SF_lower);

  for (short sec = 0; sec < NUM_SECTORS; sec++) {
    MissingMass[sec] =
        std::make_shared<TH1D>(Form("MM2_hist_sec_%d", sec), Form("MM2_hist_sec_%d", sec), bins, -w_max, w_max);

    mass_pi0_hist[before_cut][sec] =
        std::make_shared<TH1D>(Form("mass_pi0_hist_%d", sec), Form("mass_pi0_hist_%d", sec), bins, 0, 0.5);
    mass_pi0_hist[after_cut][sec] = std::make_shared<TH1D>(Form("mass_pi0_hist_aferPcuts_%d", sec),
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
  for (auto&& cut : before_after_cut) {
    int c = cut.first;
    auto type = cut.second.c_str();
    EC_sampling_fraction[c] =
        std::make_shared<TH2D>(Form("EC_sampling_fraction%s", type), Form("EC_sampling_fraction%s", type), bins, p_min,
                               p_max, bins, zero, 0.50);
    vz_position[c] = std::make_shared<TH1D>(Form("vz_position%s", type), Form("vz_position%s", type), bins, -40, 40);
    pcal_sec[c] =
        std::make_shared<TH2D>(Form("pcal_sec%s", type), Form("pcal_sec%s", type), bins, -420, 420, bins, -420, 420);
    dcr1_sec[c] =
        std::make_shared<TH2D>(Form("dcr1_sec%s", type), Form("dcr1_sec%s", type), bins, -180, 180, bins, -180, 180);
    dcr2_sec[c] =
        std::make_shared<TH2D>(Form("dcr2_sec%s", type), Form("dcr2_sec%s", type), bins, -270, 270, bins, -270, 270);
    dcr3_sec[c] =
        std::make_shared<TH2D>(Form("dcr3_sec%s", type), Form("dcr3_sec%s", type), bins, -320, 320, bins, -320, 320);
  }
}
//
// void Histogram::makeHistSF() {
//   for (size_t i = 0; i < W_BINS; i++) {
//     SF_1D[i] = std::make_shared<TH1D>(Form("sf_1D %1.12s (GeV)", W_BINS_NAME[i].c_str()),
//                                       Form("sf_1D %1.12s (GeV)", W_BINS_NAME[i].c_str()), bins, 0.0, 0.5);
//   }
// }
// void Histogram::populate_SF(const std::shared_ptr<Branches12>& _d, double min, double max, int index_sf) {
//   if (_d->p(0) > min && _d->p(0) < max) {
//     SF_1D[index_sf]->Fill(_d->ec_tot_energy(0) / _d->p(0));
//   }
// }
//
// void Histogram::Fill_SF(const std::shared_ptr<Branches12>& _d) {
//   EI_P_PCAL_P->Fill((_d->ec_ecin_energy(0) + _d->ec_ecout_energy(0)) / _d->p(0), (_d->ec_pcal_energy(0) / _d->p(0)));
//   sf_hist->Fill(_d->p(0), _d->ec_tot_energy(0) / _d->p(0));
//
//   int index_sf = 0;
//   // populate_SF(_d, 0., 0.5, 0);
//   // if (_d->p(0) > 1. && _d->p(0) < 1.5) SF_1D[2]->Fill(_d->ec_tot_energy(0) / _d->p(0));
//   for (float x = 0.; x <= 10.6; x = x + .3) {
//     float y = x + 0.3;
//     populate_SF(_d, x, y, index_sf);
//     index_sf++;
//   }
// }
void Histogram::FillHists_electron_cuts(const std::shared_ptr<Branches12>& _d) {
  vz_position[before_cut]->Fill(_d->vz(0));
  pcal_sec[before_cut]->Fill(_d->ec_pcal_x(0), _d->ec_pcal_y(0));
  dcr1_sec[before_cut]->Fill(_d->dc_r1_x(0), _d->dc_r1_y(0));
  dcr2_sec[before_cut]->Fill(_d->dc_r2_x(0), _d->dc_r2_y(0));
  dcr3_sec[before_cut]->Fill(_d->dc_r3_x(0), _d->dc_r3_y(0));
  EC_sampling_fraction[before_cut]->Fill(_d->p(0), _d->ec_tot_energy(0) / _d->p(0));
}

void Histogram::FillHists_electron_with_cuts(const std::shared_ptr<Branches12>& _d) {
  vz_position[after_cut]->Fill(_d->vz(0));
  pcal_sec[after_cut]->Fill(_d->ec_pcal_x(0), _d->ec_pcal_y(0));
  dcr1_sec[after_cut]->Fill(_d->dc_r1_x(0), _d->dc_r1_y(0));
  dcr2_sec[after_cut]->Fill(_d->dc_r2_x(0), _d->dc_r2_y(0));
  dcr3_sec[after_cut]->Fill(_d->dc_r3_x(0), _d->dc_r3_y(0));
  EC_sampling_fraction[after_cut]->Fill(_d->p(0), _d->ec_tot_energy(0) / _d->p(0));
}
//
// void Histogram::makeHists_x_mu() {
//   for (short i = 0; i < NUM_CONDITIONS; i++) {
//     E_x_mu_hist[i] = std::make_shared<TH1D>(Form("E_x_mu_ %12.20s ", NUM_CONDITIONS_NAME[i].c_str()),
//                                             Form("E_x_mu %12.20s ", NUM_CONDITIONS_NAME[i].c_str()), bins,
//                                             -2.0, 11.0);
//     diff_E2_P2_x_mu_hist[i] =
//         std::make_shared<TH1D>(Form("diff_E2_P2_x_mu_ %12.20s ", NUM_CONDITIONS_NAME[i].c_str()),
//                                Form("diff_E2_P2_x_mu %12.20s ", NUM_CONDITIONS_NAME[i].c_str()), bins, -1.0, 1.5);
//     diff_E_P_x_mu_hist[i] =
//         std::make_shared<TH1D>(Form("diff_E_P_x_mu_ %12.20s ", NUM_CONDITIONS_NAME[i].c_str()),
//                                Form("diff_E_P_x_mu %12.20s)", NUM_CONDITIONS_NAME[i].c_str()), bins, -1.0, 1.5);
//     mom_vs_E_x_mu_hist[i] = std::make_shared<TH2D>(Form("mom_vs_E_ %12.20s ", NUM_CONDITIONS_NAME[i].c_str()),
//                                                    Form("mom_vs_E %12.20s ", NUM_CONDITIONS_NAME[i].c_str()), bins,
//                                                    -1.0, 3.0, bins, 0.0, 3.0);
//
//     theta_elec_hist[i] =
//         std::make_shared<TH1D>(Form("theta_elec_scattering_ %12.20s ", NUM_CONDITIONS_NAME[i].c_str()),
//                                Form("theta_elec_scattering %12.20s ", NUM_CONDITIONS_NAME[i].c_str()), bins,
//                                0.0, 50.0);
//     theta_x_mu_hist[i] =
//         std::make_shared<TH1D>(Form("theta_x_mu_ %12.20s ", NUM_CONDITIONS_NAME[i].c_str()),
//                                Form("theta_x_mu %12.20s ", NUM_CONDITIONS_NAME[i].c_str()), bins, 0.0, 180.0);
//
//     diff_theta_elec_x_mu_hist[i] = std::make_shared<TH1D>(
//         Form("diff_theta_elec_x_mu_ %12.20s ", NUM_CONDITIONS_NAME[i].c_str()),
//         Form("diff_theta_elec_x_mu %12.20s ", NUM_CONDITIONS_NAME[i].c_str()), bins, -50.0, 180.0);
//   }
// }
// void Histogram::Fill_x_mu(const std::shared_ptr<Reaction>& _e) {
//   theta_elec_hist[0]->Fill(_e->theta_elec());
//   theta_x_mu_hist[0]->Fill(_e->theta_x_mu());
//   if (_e->onePositive()) {
//     diff_theta_in_x_mu->Fill(-_e->theta_beam() + _e->theta_x_mu());
//     E_x_mu_hist[0]->Fill(_e->E_x_mu());
//     diff_E2_P2_x_mu_hist[0]->Fill(_e->M2_x_mu());
//     diff_E_P_x_mu_hist[0]->Fill(_e->M_x_mu());
//     mom_vs_E_x_mu_hist[0]->Fill(_e->E_x_mu(), _e->P_x_mu());
//     diff_theta_elec_x_mu_hist[0]->Fill(_e->theta_x_mu() - _e->theta_elec());
//   }
//
//   if (_e->onePositive_noOther()) {
//     E_x_mu_hist[1]->Fill(_e->E_x_mu());
//     diff_E2_P2_x_mu_hist[1]->Fill(_e->M2_x_mu());
//     diff_E_P_x_mu_hist[1]->Fill(_e->M_x_mu());
//     mom_vs_E_x_mu_hist[1]->Fill(_e->E_x_mu(), _e->P_x_mu());
//     theta_x_mu_hist[1]->Fill(_e->theta_x_mu());
//     theta_elec_hist[1]->Fill(_e->theta_elec());
//     diff_theta_elec_x_mu_hist[1]->Fill(_e->theta_x_mu() - _e->theta_elec());
//   }
//
//   if (_e->onePositive_at180()) {
//     E_x_mu_hist[2]->Fill(_e->E_x_mu());
//     diff_E2_P2_x_mu_hist[2]->Fill(_e->M2_x_mu());
//     diff_E_P_x_mu_hist[2]->Fill(_e->M_x_mu());
//     mom_vs_E_x_mu_hist[2]->Fill(_e->E_x_mu(), _e->P_x_mu());
//     theta_x_mu_hist[2]->Fill(_e->theta_x_mu());
//     theta_elec_hist[2]->Fill(_e->theta_elec());
//     diff_theta_elec_x_mu_hist[2]->Fill(_e->theta_x_mu() - _e->theta_elec());
//   }
//   if (_e->onePositive_at180_MM0()) {
//     E_x_mu_hist[3]->Fill(_e->E_x_mu());
//     diff_E2_P2_x_mu_hist[3]->Fill(_e->M2_x_mu());
//     diff_E_P_x_mu_hist[3]->Fill(_e->M_x_mu());
//     mom_vs_E_x_mu_hist[3]->Fill(_e->E_x_mu(), _e->P_x_mu());
//     theta_x_mu_hist[3]->Fill(_e->theta_x_mu());
//     theta_elec_hist[3]->Fill(_e->theta_elec());
//     diff_theta_elec_x_mu_hist[3]->Fill(_e->theta_x_mu() - _e->theta_elec());
//
//     if (abs(_e->E_x_mu() * _e->E_x_mu()) > 0.1 && abs(_e->P_x_mu() * _e->P_x_mu()) > 0.1 &&
//         (_e->E_x_mu() * _e->E_x_mu() - _e->P_x_mu() * _e->P_x_mu()) < 0.2 && (_e->theta_x_mu() > 4.0)) {
//       E_x_mu_hist[4]->Fill(_e->E_x_mu());
//       diff_E2_P2_x_mu_hist[4]->Fill(_e->M2_x_mu());
//       diff_E_P_x_mu_hist[4]->Fill(_e->M_x_mu());
//       mom_vs_E_x_mu_hist[4]->Fill(_e->E_x_mu(), _e->P_x_mu());
//       theta_elec_hist[4]->Fill(_e->theta_elec());
//       theta_x_mu_hist[4]->Fill(_e->theta_x_mu());
//       diff_theta_elec_x_mu_hist[4]->Fill(_e->theta_x_mu() - _e->theta_elec());
//
//       if (_e->E_x_mu() > 0.2) {
//         E_x_mu_hist[5]->Fill(_e->E_x_mu());
//         diff_E2_P2_x_mu_hist[5]->Fill(_e->M2_x_mu());
//         diff_E_P_x_mu_hist[5]->Fill(_e->M_x_mu());
//         mom_vs_E_x_mu_hist[5]->Fill(_e->E_x_mu(), _e->P_x_mu());
//         theta_elec_hist[5]->Fill(_e->theta_elec());
//         theta_x_mu_hist[5]->Fill(_e->theta_x_mu());
//         diff_theta_elec_x_mu_hist[5]->Fill(_e->theta_x_mu() - _e->theta_elec());
//
//         P_x_mu->Fill(_e->P_x_mu());
//         Px_x_mu->Fill(_e->Px_x_mu());
//         Py_x_mu->Fill(_e->Py_x_mu());
//         Pz_x_mu->Fill(_e->Pz_x_mu());
//         diff_theta_ph_x_mu->Fill(-_e->theta_beam() + _e->theta_x_mu());
//       }
//     }
//   }
// }
// void Histogram::write_hist_x_mu() {
//   for (short i; i < NUM_CONDITIONS; i++) {
//     E_x_mu_hist[i]->Fit("gaus", "QMR+", "QMR+", -0.45, 0.38);
//     //  gROOT->SetStyle("Plain");
//     gStyle->SetOptFit(1111);
//     E_x_mu_hist[i]->SetXTitle("Energy comp (GeV)");
//     E_x_mu_hist[i]->Write();
//     diff_E2_P2_x_mu_hist[i]->SetXTitle("E2-P2 (GeV)");
//     diff_E2_P2_x_mu_hist[i]->Write();
//     diff_E_P_x_mu_hist[i]->SetXTitle("E-P (GeV)");
//     diff_E_P_x_mu_hist[i]->Write();
//     mom_vs_E_x_mu_hist[i]->SetXTitle("Energy (GeV)");
//     mom_vs_E_x_mu_hist[i]->SetYTitle("Mom (GeV)");
//     mom_vs_E_x_mu_hist[i]->SetOption("COLZ");
//     mom_vs_E_x_mu_hist[i]->Write();
//     theta_elec_hist[i]->SetXTitle("theta (deg)");
//     theta_elec_hist[i]->Write();
//     theta_x_mu_hist[i]->SetXTitle("theta (deg)");
//     theta_x_mu_hist[i]->Write();
//     diff_theta_elec_x_mu_hist[i]->SetXTitle("theta (deg)");
//     diff_theta_elec_x_mu_hist[i]->Write();
//   }
// }
// void Histogram::Write_SF() {
//   sf_hist->SetOption("COLZ");
//   sf_hist->Write();
//   gStyle->SetOptFit(1111);
//   // set stat
//   SF_gr_upper->Fit("pol2");
//   SF_gr_upper->SetTitle(" sf upper 3#sigma cut");
//
//   SF_gr_upper->Write();
//   SF_gr_lower->Fit("pol2");
//   SF_gr_lower->SetTitle(" sf lower 3#sigma cut");
//
//   SF_gr_lower->Write();
//   EI_P_PCAL_P->SetOption("COLZ");
//   EI_P_PCAL_P->SetXTitle("PCAL/P");
//   EI_P_PCAL_P->SetYTitle("Etot/P");
//   EI_P_PCAL_P->Write();
//
//   P_x_mu->SetXTitle("Momentum (GeV)");
//   P_x_mu->Write();
//   Px_x_mu->SetXTitle("Px (GeV)");
//   Px_x_mu->Write();
//   Py_x_mu->SetXTitle("Py (GeV)");
//   Py_x_mu->Write();
//   Pz_x_mu->SetXTitle("Pz (GeV)");
//   Pz_x_mu->Write();
//   diff_theta_ph_x_mu->SetXTitle("theta_diff (deg)");
//   diff_theta_in_x_mu->SetXTitle("theta_diff (deg)");
//   diff_theta_ph_x_mu->Write();
//   diff_theta_in_x_mu->Write();
//
//   // Dthtea_vs_Dphi->SetXTitle("#Delata#Theta (deg)");
//   // Dthtea_vs_Dphi->SetYTitle("#Delata#Phi (deg)");
//   // Dthtea_vs_Dphi->SetOption("COLZ");
//   // Dthtea_vs_Dphi->Write();
// }
// void Histogram::write_histSf() {
//   for (size_t i = 0; i < W_BINS; i++) {
//     SF_1D[i]->SetXTitle("sf (E/P)");
//     SF_1D[i]->Fit("gaus", "QMR+", "QMR+", 0.225, 0.28);
//     //  gROOT->SetStyle("Plain");
//     gStyle->SetOptFit(1111);
//     if (SF_1D[i]->GetEntries()) SF_1D[i]->Write();
//   }
// }
void Histogram::Write_Electron_cuts() {
  for (auto&& cut : before_after_cut) {
    int c = cut.first;
    vz_position[c]->SetXTitle("vz (GeV)");
    vz_position[c]->Fit("gaus", "QMR+", "QMR+", -7.089, 2.0);
    gStyle->SetOptFit(1111);
    if (vz_position[c]->GetEntries()) vz_position[c]->Write();
    pcal_sec[c]->SetXTitle("x (cm)");
    pcal_sec[c]->SetYTitle("y (cm)");
    pcal_sec[c]->SetOption("COLZ1");
    if (pcal_sec[c]->GetEntries()) pcal_sec[c]->Write();

    dcr1_sec[c]->SetXTitle("x (cm)");
    dcr1_sec[c]->SetYTitle("y (cm)");
    dcr1_sec[c]->SetOption("COLZ1");
    if (dcr1_sec[c]->GetEntries()) dcr1_sec[c]->Write();

    dcr2_sec[c]->SetXTitle("x (cm)");
    dcr2_sec[c]->SetYTitle("y (cm)");
    dcr2_sec[c]->SetOption("COLZ1");
    if (dcr2_sec[c]->GetEntries()) dcr2_sec[c]->Write();

    dcr3_sec[c]->SetXTitle("x (cm)");
    dcr3_sec[c]->SetYTitle("y (cm)");
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
  //  if ((sec > 0 && sec < NUM_SECTORS) || pos_det != -1) {
  if (_e->NPip()) {
    MM_hist_NPip->Fill(_e->MM());
    MM2_hist_NPip->Fill(_e->MM2());
    // if (MM2<)
    W_hist_all_events[all_sectors]->Fill(_e->W());
    W_vs_q2_all_events[all_sectors]->Fill(_e->W(), _e->Q2());
    W_hist_all_events[sec]->Fill(_e->W());
    W_vs_q2_all_events[sec]->Fill(_e->W(), _e->Q2());

    if (_e->onePositive()) {
      W_hist_1pos[all_sectors]->Fill(_e->W());
      W_vs_q2_1pos[all_sectors]->Fill(_e->W(), _e->Q2());
      W_hist_1pos[sec]->Fill(_e->W());
      W_vs_q2_1pos[sec]->Fill(_e->W(), _e->Q2());

      Phie_vs_Phip[both_detectors][all_sectors]->Fill(_e->phi_e(), _e->phi_p());
      Phie_Phip_hist[both_detectors][all_sectors]->Fill(_e->phi_diff());
      Phie_vs_Phip[both_detectors][sec]->Fill(_e->phi_e(), _e->phi_p());
      Phie_Phip_hist[both_detectors][sec]->Fill(_e->phi_diff());

      Phie_vs_Phip[pos_det][all_sectors]->Fill(_e->phi_e(), _e->phi_p());
      Phie_Phip_hist[pos_det][all_sectors]->Fill(_e->phi_diff());
      Phie_vs_Phip[pos_det][sec]->Fill(_e->phi_e(), _e->phi_p());
      Phie_Phip_hist[pos_det][sec]->Fill(_e->phi_diff());
    }
    if (_e->onePositive_MM0()) {
      W_hist_1pos_0charge[all_sectors]->Fill(_e->W());
      W_vs_q2_1pos_0charge[all_sectors]->Fill(_e->W(), _e->Q2());
      W_hist_1pos_0charge[sec]->Fill(_e->W());
      W_vs_q2_1pos_0charge[sec]->Fill(_e->W(), _e->Q2());
    }
    if (_e->onePositive_noOther()) {
      W_hist_1pos_noOther[all_sectors]->Fill(_e->W());
      W_vs_q2_1pos_noOther[all_sectors]->Fill(_e->W(), _e->Q2());
      W_hist_1pos_noOther[sec]->Fill(_e->W());
      W_vs_q2_1pos_noOther[sec]->Fill(_e->W(), _e->Q2());
    }
    if (_e->onePositive_at180()) {
      MissingMass[all_sectors]->Fill(_e->MM2());
      MissingMass[sec]->Fill(_e->MM2());

      W_hist_1pos_at180[both_detectors][all_sectors]->Fill(_e->W());
      W_vs_q2_1pos_at180[both_detectors][all_sectors]->Fill(_e->W(), _e->Q2());
      W_hist_1pos_at180[both_detectors][sec]->Fill(_e->W());
      W_vs_q2_1pos_at180[both_detectors][sec]->Fill(_e->W(), _e->Q2());

      W_hist_1pos_at180[pos_det][all_sectors]->Fill(_e->W());
      W_vs_q2_1pos_at180[pos_det][all_sectors]->Fill(_e->W(), _e->Q2());
      W_hist_1pos_at180[pos_det][sec]->Fill(_e->W());
      W_vs_q2_1pos_at180[pos_det][sec]->Fill(_e->W(), _e->Q2());
    }
    if (_e->onePositive_at180_MM0()) {
      W_hist_1pos_at180_MM[both_detectors][all_sectors]->Fill(_e->W());
      W_vs_q2_1pos_at180_MM[both_detectors][all_sectors]->Fill(_e->W(), _e->Q2());
      W_hist_1pos_at180_MM[both_detectors][sec]->Fill(_e->W());
      W_vs_q2_1pos_at180_MM[both_detectors][sec]->Fill(_e->W(), _e->Q2());

      W_hist_1pos_at180_MM[pos_det][all_sectors]->Fill(_e->W());
      W_vs_q2_1pos_at180_MM[pos_det][all_sectors]->Fill(_e->W(), _e->Q2());
      W_hist_1pos_at180_MM[pos_det][sec]->Fill(_e->W());
      W_vs_q2_1pos_at180_MM[pos_det][sec]->Fill(_e->W(), _e->Q2());
    }
  }
}

void Histogram::Write_WvsQ2() {
  MM_hist_NPip->SetXTitle("MM (GeV)");
  MM_hist_NPip->Write();
  MM2_hist_NPip->SetXTitle("MM2 (GeV)");
  MM2_hist_NPip->Write();
  TDirectory* phi_folder = RootOutputFile->mkdir("Phi");
  phi_folder->cd();
  for (short j = 0; j < detector_name.size(); j++) {
    for (int i = 0; i < NUM_SECTORS; i++) {
      Phie_vs_Phip[j][i]->SetXTitle("Phie");
      Phie_vs_Phip[j][i]->SetYTitle("Phip");
      Phie_vs_Phip[j][i]->SetOption("COLZ");
      Phie_vs_Phip[j][i]->Write();
    }
    for (int i = 0; i < NUM_SECTORS; i++) {
      Phie_Phip_hist[j][i]->SetXTitle("Phi");
      Phie_Phip_hist[j][i]->Write();
    }
  }
  phi_folder->Write();
  delete phi_folder;

  TDirectory* at180_folder = RootOutputFile->mkdir("at180");
  at180_folder->cd();
  for (short j = 0; j < detector_name.size(); j++) {
    for (int i = 0; i < NUM_SECTORS; i++) {
      W_hist_1pos_at180[j][i]->SetXTitle("W (GeV)");
      W_hist_1pos_at180[j][i]->Write();
    }
    for (int i = 0; i < NUM_SECTORS; i++) {
      W_vs_q2_1pos_at180[j][i]->SetXTitle("W (GeV)");
      W_vs_q2_1pos_at180[j][i]->SetYTitle("Q^2 (GeV^2)");
      W_vs_q2_1pos_at180[j][i]->SetOption("COLZ");
      W_vs_q2_1pos_at180[j][i]->Write();
    }
  }
  at180_folder->Write();
  delete at180_folder;

  TDirectory* at180_MM_folder = RootOutputFile->mkdir("at180_MM");
  at180_MM_folder->cd();
  for (short j = 0; j < detector_name.size(); j++) {
    for (int i = 0; i < NUM_SECTORS; i++) {
      W_hist_1pos_at180_MM[j][i]->SetXTitle("W (GeV)");
      W_hist_1pos_at180_MM[j][i]->Write();
    }
    for (int i = 0; i < NUM_SECTORS; i++) {
      W_vs_q2_1pos_at180_MM[j][i]->SetXTitle("W (GeV)");
      W_vs_q2_1pos_at180_MM[j][i]->SetYTitle("Q^2 (GeV^2)");
      W_vs_q2_1pos_at180_MM[j][i]->SetOption("COLZ");
      W_vs_q2_1pos_at180_MM[j][i]->Write();
    }
  }
  at180_MM_folder->Write();
  delete at180_MM_folder;

  TDirectory* W_vs_Q2_folder = RootOutputFile->mkdir("W_vs_Q2");
  W_vs_Q2_folder->cd();
  for (int i = 0; i < NUM_SECTORS; i++) {
    MissingMass[i]->SetXTitle("MM^2 (GeV)");
    MissingMass[i]->Write();

    mass_pi0_hist[before_cut][i]->SetXTitle("MM(GeV)");
    mass_pi0_hist[before_cut][i]->Write();

    mass_pi0_hist[after_cut][i]->SetXTitle("MM(GeV)");
    mass_pi0_hist[after_cut][i]->Write();
    W_hist_all_events[i]->Fit("gaus", "QMR+", "QMR+", 0.8, 1.1);
    gStyle->SetOptFit(1111);
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
  W_vs_Q2_folder->Write();
  delete W_vs_Q2_folder;
}

void Histogram::Fill_MomVsBeta(const std::shared_ptr<Reaction>& _e) {
  if (!_e->onePositive_at180()) return;
  if (_e->pos_det() == -1) return;

  MomVsBeta[both_detectors][all_sectors]->Fill(_e->pos_P(), _e->pos_beta());
  MomVsBeta[both_detectors][_e->sec()]->Fill(_e->pos_P(), _e->pos_beta());

  ThetaVsP[both_detectors][all_sectors]->Fill(_e->pos_P(), _e->pos_theta());
  ThetaVsP[both_detectors][_e->sec()]->Fill(_e->pos_P(), _e->pos_theta());

  ThetaVsPCalc[both_detectors][all_sectors]->Fill(_e->pos_P(), _e->pos_theta_calc());
  ThetaVsPCalc[both_detectors][_e->sec()]->Fill(_e->pos_P(), _e->pos_theta_calc());

  MomVsBeta[_e->pos_det()][all_sectors]->Fill(_e->pos_P(), _e->pos_beta());
  MomVsBeta[_e->pos_det()][_e->sec()]->Fill(_e->pos_P(), _e->pos_beta());

  ThetaVsP[_e->pos_det()][all_sectors]->Fill(_e->pos_P(), _e->pos_theta());
  ThetaVsP[_e->pos_det()][_e->sec()]->Fill(_e->pos_P(), _e->pos_theta());

  ThetaVsPCalc[_e->pos_det()][all_sectors]->Fill(_e->pos_P(), _e->pos_theta_calc());
  ThetaVsPCalc[_e->pos_det()][_e->sec()]->Fill(_e->pos_P(), _e->pos_theta_calc());

  if (_e->W() < 2.0) {
    ThetaVsP_lowW[both_detectors][all_sectors]->Fill(_e->pos_P(), _e->pos_theta());
    ThetaVsP_lowW[both_detectors][_e->sec()]->Fill(_e->pos_P(), _e->pos_theta());
    ThetaVsP_lowW[_e->pos_det()][all_sectors]->Fill(_e->pos_P(), _e->pos_theta());
    ThetaVsP_lowW[_e->pos_det()][_e->sec()]->Fill(_e->pos_P(), _e->pos_theta());
  }
}

void Histogram::Write_MomVsBeta() {
  for (short i = 0; i < detector_name.size(); i++) {
    for (short p = 0; p < NUM_SECTORS; p++) {
      MomVsBeta[i][p]->SetXTitle("Momentum (GeV)");
      MomVsBeta[i][p]->SetYTitle("#beta");
      MomVsBeta[i][p]->SetOption("COLZ1");
      MomVsBeta[i][p]->Write();
    }
    for (short p = 0; p < NUM_SECTORS; p++) {
      ThetaVsP[i][p]->SetXTitle("Momentum (GeV)");
      ThetaVsP[i][p]->SetYTitle("#theta");
      ThetaVsP[i][p]->SetOption("COLZ1");
      ThetaVsP[i][p]->Write();
    }
    for (short p = 0; p < NUM_SECTORS; p++) {
      ThetaVsP_lowW[i][p]->SetXTitle("Momentum (GeV)");
      ThetaVsP_lowW[i][p]->SetYTitle("#theta");
      ThetaVsP_lowW[i][p]->SetOption("COLZ1");
      ThetaVsP_lowW[i][p]->Write();
    }
    for (short p = 0; p < NUM_SECTORS; p++) {
      ThetaVsPCalc[i][p]->SetXTitle("Momentum (GeV)");
      ThetaVsPCalc[i][p]->SetYTitle("#theta");
      ThetaVsPCalc[i][p]->SetOption("COLZ1");
      ThetaVsPCalc[i][p]->Write();
    }
  }
}

void Histogram::Fill_Dt(const std::shared_ptr<Delta_T>& dt) {
  for (int i = 0; i < dt->gpart(); i++) {
    if (dt->charge(i) == 1) deltaT_proton[before_cut]->Fill(dt->mom(i), dt->dt_P(i));
  }
}

void Histogram::Fill_Dt(const std::shared_ptr<Delta_T>& dt, int part) {
  deltaT_proton[after_cut]->Fill(dt->mom(part), dt->dt_P(part));
}

void Histogram::Fill_pi0(const std::shared_ptr<Reaction>& _e) {
  if (_e->pi0_mass() < 0.0001) return;
  short sec = _e->sec();
  short pos_det = _e->pos_det();
  mass_pi0_hist[before_cut][all_sectors]->Fill(_e->pi0_mass());
  if ((sec > 0 && sec < NUM_SECTORS) || pos_det != -1) {
    mass_pi0_hist[before_cut][sec]->Fill(_e->pi0_mass());
    if (_e->onePositive_at180_MM0()) {
      mass_pi0_hist[after_cut][all_sectors]->Fill(_e->pi0_mass());
      mass_pi0_hist[after_cut][sec]->Fill(_e->pi0_mass());
    }
  }
}
