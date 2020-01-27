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
// static const int W_BINS = 35;
// for loop le array create garna sakinchha?
static const short NUM_CONDITIONS = 6;
std::string NUM_CONDITIONS_NAME[NUM_CONDITIONS] = {"onePositive ",      "noOther_onePos",  " onePos_at180",
                                                   " at180_MM0_onePos", "E2_P2_Condition", "higher_then_200 MeV"};

// std::string W_BINS_NAME[W_BINS] = {"0.0-0.3", "0.3-0.6", "0.6=0.9", "0.9-1.2", "1.2-1.5", "1.5-1.8",  "1.8-2.1",
//                                    "2.1-2.4", "2.4-2.7", "2.7-3.0", "3.0-3.3", "3.3-3.6", "3.6-3.9",  "3.9-4.2",
//                                    "4.2-4.5", "4.5-4.8", "4.8-5.1", "5.1-5.4", "5.4-5.7", "5.7-6.0",  " 6.0-6.3",
//                                    "6.3-6.6", "6.6-6.9", "6.9-7.2", "7.2-7.5", "7.5-7.8", "7.8-8.1",  "8.1-8.4",
//                                    "8.4-8.7", "8.7-9.0", "9.0-9.3", "9.3-9.6", "9.6-9.9", "9.9-10.2", "10.2-10.5"};

TH2D_ptr sf_hist = std::make_shared<TH2D>("SF", "SF", 500, 0, 10.5, 500, 0, 0.5);
TH2D_ptr EI_P_PCAL_P = std::make_shared<TH2D>("EI/P VS PCAL/P", "EI/P VS PCAL/P", 500, 0, 0.35, 500, 0, 0.35);
TH1D_ptr P_x_mu = std::make_shared<TH1D>("mom (x_mu = e(p,p')e')", "mom (x_mu = e(p,p')e')", 500, -0.50, 11.0);
TH1D_ptr Px_x_mu = std::make_shared<TH1D>("Px (x_mu = e(p,p')e')", "Px (x_mu = e(p,p')e')", 500, -0.50, 2.0);
TH1D_ptr Py_x_mu = std::make_shared<TH1D>("Py (x_mu = e(p,p')e')", "Py (x_mu = e(p,p')e')", 500, -0.50, 2.0);
TH1D_ptr Pz_x_mu = std::make_shared<TH1D>("Pz (x_mu = e(p,p')e')", "Pz (x_mu = e(p,p')e')", 500, -1.0, 11.0);

TH1D_ptr diff_theta_in_x_mu = std::make_shared<TH1D>("diff#theta x_mu and initial electron",
                                                     "diff#theta x_mu and initial electron", 500, -5.0, 180);
TH1D_ptr diff_theta_ph_x_mu = std::make_shared<TH1D>("diff#theta x_mu and initial electron ph",
                                                     "diff#theta x_mu and initial electron ph", 500, 0.0, 180);

TH1D_ptr MM_hist_NPip_before_cut = std::make_shared<TH1D>("MM_hist_NPip_before_cut", "MM_hist_NPip_before_cut", 500, p_min, w_max);
TH1D_ptr MM2_hist_NPip_before_cut = std::make_shared<TH1D>("MM2_hist_NPip_before_cut", "MM_hist_NPip_before_cut", 500, p_min, 16.0);
TH1D_ptr MM_hist_NPip = std::make_shared<TH1D>("MM_hist_NPip", "MM_hist_NPip", 500, p_min, w_max);
TH1D_ptr MM2_hist_NPip = std::make_shared<TH1D>("MM2_hist_NPip", "MM_hist_NPip", 500, p_min, w_max);

TH1D_ptr E_x_mu_hist[NUM_CONDITIONS];
TH1D_ptr diff_E2_P2_x_mu_hist[NUM_CONDITIONS];
TH1D_ptr diff_E_P_x_mu_hist[NUM_CONDITIONS];
TH2D_ptr mom_vs_E_x_mu_hist[NUM_CONDITIONS];
TH1D_ptr theta_elec_hist[NUM_CONDITIONS];
TH1D_ptr theta_x_mu_hist[NUM_CONDITIONS];
TH1D_ptr diff_theta_elec_x_mu_hist[NUM_CONDITIONS];

TH1D_ptr vz_position[NUM_CUT];
TH2D_ptr pcal_sec[NUM_CUT];
TH2D_ptr dcr1_sec[NUM_CUT];
TH2D_ptr dcr2_sec[NUM_CUT];
TH2D_ptr dcr3_sec[NUM_CUT];
TH2D_ptr EC_sampling_fraction[NUM_CUT];
// SF 1D
// TH1D_ptr SF_1D[W_BINS];
// TGraph_ptr SF_gr_upper;
// TGraph_ptr SF_gr_lower;
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
void populate_SF(const std::shared_ptr<Branches12>& _d, double min, double max, int index_sf);
void Fill_SF(const std::shared_ptr<Branches12>& _d);
void write_histSf();
void Write_SF();
void makeHists_x_mu();
void Fill_x_mu(const std::shared_ptr<Reaction>& _e);
void write_hist_x_mu();

// P and E
void Fill_MomVsBeta(const std::shared_ptr<Reaction>& _e);
void Write_MomVsBeta();

// ecectron cuts
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
