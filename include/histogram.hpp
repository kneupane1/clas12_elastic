/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

#ifndef HIST_H_GUARD
#define HIST_H_GUARD
#include <mutex>
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "THnSparse.h"
#include "TLegend.h"
#include "TLorentzVector.h"
#include "TPaveStats.h"
#include "TROOT.h"
#include "TStyle.h"
#include "colors.hpp"
#include "constants.hpp"
#include "cuts.hpp"
#include "deltat.hpp"
#include "reaction.hpp"

using TH2D_ptr = std::shared_ptr<TH2D>;
using TH1D_ptr = std::shared_ptr<TH1D>;
using TGraph_ptr = std::shared_ptr<TGraph>;

class Histogram {
 protected:
  static const short NUM_SECTORS = 7;
  static const short NUM_DET = 3;
  static const short CUTS = 2;
  // enum to easily access detector and sector information
  enum cuts { before_cut, after_cut };
  enum detector { both_detectors, forward, central };
  enum sector { all_sectors, one, two, three, four, five };

  // Mutex needed for filling some histograms
  std::mutex mutex;

  // Output file
  std::shared_ptr<TFile> RootOutputFile;
  // Default canvas
  std::shared_ptr<TCanvas> def;

  int bins = 500;
  double p_min = 0.0;
  double Dt_max = 10.0;
  double Dt_min = -Dt_max;
  double q2_max = 10.6;
  double w_max = 4.5;
  double p_max = 10.6;

  double zero = 0.0;

  static const short NUM_DIM = 3;
  //// W, Q2, sector
  int sparce_bins[NUM_DIM] = {bins, 10, 7};
  double sparce_xmin[NUM_DIM] = {zero, zero, 0};
  double sparce_xmax[NUM_DIM] = {w_max, q2_max, 6};

  static const short NUM_CUT = 2;

  TH2D_ptr sf_hist = std::make_shared<TH2D>("SF", "SF", 500, 0, 10.5, 500, 0, 1);
  TH1D_ptr vz_position[NUM_CUT];
  TH2D_ptr pcal_sec[NUM_CUT];
  TH2D_ptr dcr1_sec[NUM_CUT];
  TH2D_ptr dcr2_sec[NUM_CUT];
  TH2D_ptr dcr3_sec[NUM_CUT];
  TH2D_ptr EC_sampling_fraction[NUM_CUT];
  // Kinematics
  TH1D_ptr W_hist_all_events[NUM_SECTORS];
  TH1D_ptr W_hist_1pos[NUM_SECTORS];
  TH1D_ptr W_hist_1pos_0charge[NUM_SECTORS];
  TH1D_ptr W_hist_1pos_noOther[NUM_SECTORS];
  TH1D_ptr W_hist_1pos_at180[NUM_DET][NUM_SECTORS];
  TH1D_ptr W_hist_1pos_at180_MM[NUM_DET][NUM_SECTORS];

  TH2D_ptr W_vs_q2_all_events[NUM_SECTORS];
  TH2D_ptr W_vs_q2_1pos[NUM_SECTORS];
  TH2D_ptr W_vs_q2_1pos_0charge[NUM_SECTORS];
  TH2D_ptr W_vs_q2_1pos_noOther[NUM_SECTORS];
  TH2D_ptr W_vs_q2_1pos_at180[NUM_DET][NUM_SECTORS];
  TH2D_ptr W_vs_q2_1pos_at180_MM[NUM_DET][NUM_SECTORS];

  TH2D_ptr ThetaVsP[NUM_DET][NUM_SECTORS];
  TH2D_ptr ThetaVsPCalc[NUM_DET][NUM_SECTORS];
  TH2D_ptr ThetaVsP_lowW[NUM_DET][NUM_SECTORS];
  TH2D_ptr MomVsBeta[NUM_DET][NUM_SECTORS];

  TH2D_ptr Phie_vs_Phip[NUM_DET][NUM_SECTORS];
  TH1D_ptr Phie_Phip_hist[NUM_DET][NUM_SECTORS];

  TH1D_ptr MissingMass[NUM_SECTORS];
  TH1D_ptr mass_pi0_hist[CUTS][NUM_SECTORS];
  TH2D_ptr deltaT_proton[CUTS];

  std::shared_ptr<THnSparse> Nsparce;

 public:
  Histogram(const std::string& output_file);
  ~Histogram();

  // W and Q^2
  void makeHists();
  void Fill_WvsQ2(const std::shared_ptr<Reaction>& _e);
  void Write_WvsQ2();
  // sampling Fraction
  void makeHistSF();
  void write_histSf();
  // P and E
  void Fill_MomVsBeta(const std::shared_ptr<Reaction>& _e);
  void Write_MomVsBeta();

  void Fill_SF(const std::shared_ptr<Branches12>& _d);
  void Write_SF();
  void makeHists_electron_cuts();
  void FillHists_electron_cuts(const std::shared_ptr<Branches12>& _d);
  void FillHists_electron_with_cuts(const std::shared_ptr<Branches12>& _d);

  void Write_Electron_cuts();
  void Fill_Sparce(const std::shared_ptr<Reaction>& _e);
  void Fill_Dt(const std::shared_ptr<Delta_T>& dt);
  void Fill_Dt(const std::shared_ptr<Delta_T>& dt, int part);
  void Fill_pi0(const std::shared_ptr<Reaction>& _e);

  //
  void Write();
};

#endif
