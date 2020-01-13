/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/
#include "cuts.hpp"
#include <iostream>

Cuts::Cuts(const std::shared_ptr<Branches12>& data) : _data(data) { _dt = std::make_shared<Delta_T>(data); }
Cuts::Cuts(const std::shared_ptr<Branches12>& data, const std::shared_ptr<Delta_T>& dt) : _data(data), _dt(dt) {}
Cuts::~Cuts() {}

bool Cuts::ElectronCuts() {
  bool _elec = true;
  // Number of good particles is greater than 0
  // So that we can check at(0) without errors
  _elec &= (_data->gpart() > 0);
  if (!_elec) return false;
  //_elec &= !std::isnan(_data->cc_nphe_tot(0));

  _elec &= (_data->gpart() < 20);
  _elec &= (_data->charge(0) == NEGATIVE);
  _elec &= (_data->pid(0) == ELECTRON);
  // Why 1.0 for minimumm momentum cut?
  _elec &= (_data->p(0) > 1.0);
  _elec &= ((abs(_data->status(0)) >= 2000) && abs(_data->status(0)) < 4000);
  _elec &= (_data->vz(0) > -7.9 && _data->vz(0) < 2.0);
  // Use the chi2pid instead of straight line cuts on SF
  _elec &= (abs(_data->chi2pid(0)) < 3);

  // FiducialCuts is the slowest of the cuts because of all the calcuations
  // If it already fails a different cut we will quit before
  // calulating for the FiducialCuts to save time
  if (!_elec) return _elec;
  _elec &= FiducialCuts();

  return _elec;
}

bool Cuts::FiducialCuts() {
  bool _fid_cut = true;
  // DC sector never changes so get it once and store it to use all the time
  short dc_sec = (_data->dc_sec(0) - 1);
  // Same with these values
  float sin_dc_sec = sinf(dc_sec * ROTATE);
  float cos_dc_sec = cosf(dc_sec * ROTATE);

  float x_PCAL_rot = _data->ec_pcal_y(0) * sin_dc_sec + _data->ec_pcal_x(0) * cos_dc_sec;
  float y_PCAL_rot = _data->ec_pcal_y(0) * cos_dc_sec - _data->ec_pcal_x(0) * sin_dc_sec;

  float left_PCAL = (HEIGHT_PCAL - SLOPE * y_PCAL_rot);
  float right_PCAL = (HEIGHT_PCAL + SLOPE * y_PCAL_rot);
  float radius2_PCAL = X_SQUARE_PCAL - (y_PCAL_rot * y_PCAL_rot);  // circle radius r^2 = x^2 + y^2

  // I do this to clean up what is happening and makse sure that the cuts are not ambiguous
  _fid_cut &= (x_PCAL_rot > left_PCAL);
  _fid_cut &= (x_PCAL_rot > right_PCAL);
  _fid_cut &= (x_PCAL_rot * x_PCAL_rot > radius2_PCAL);
  _fid_cut &= (x_PCAL_rot < 372);

  // If it fails pcal cut return before calculating DC cut to save time
  if (!_fid_cut) return _fid_cut;

  float x1_rot = _data->dc_r1_y(0) * sin_dc_sec + _data->dc_r1_x(0) * cos_dc_sec;
  float y1_rot = _data->dc_r1_y(0) * cos_dc_sec - _data->dc_r1_x(0) * sin_dc_sec;
  float left_r1 = (DCR1_HEIGHT - SLOPE * y1_rot);
  float right_r1 = (DCR1_HEIGHT + SLOPE * y1_rot);
  float radius2_DCr1 = DCR1_SQUARE - (y1_rot * y1_rot);

  _fid_cut &= (x1_rot > left_r1);
  _fid_cut &= (x1_rot > right_r1);
  _fid_cut &= (x1_rot * x1_rot > radius2_DCr1);

  // If it fails cut return before calculating cut to save time
  if (!_fid_cut) return _fid_cut;

  float x2_rot = _data->dc_r2_y(0) * sin_dc_sec + _data->dc_r2_x(0) * cos_dc_sec;
  float y2_rot = _data->dc_r2_y(0) * cos_dc_sec - _data->dc_r2_x(0) * sin_dc_sec;
  float left_r2 = (DCR2_HEIGHT - SLOPE * y2_rot);
  float right_r2 = (DCR2_HEIGHT + SLOPE * y2_rot);
  float radius2_DCr2 = DCR2_SQUARE - (y2_rot * y2_rot);

  _fid_cut &= (x2_rot > left_r2);
  _fid_cut &= (x2_rot > right_r2);
  _fid_cut &= ((x2_rot * x2_rot) > radius2_DCr2);

  // If it fails cut return before calculating cut to save time
  if (!_fid_cut) return _fid_cut;

  float x3_rot = _data->dc_r3_y(0) * sin_dc_sec + _data->dc_r3_x(0) * cos_dc_sec;
  float y3_rot = _data->dc_r3_y(0) * cos_dc_sec - _data->dc_r3_x(0) * sin_dc_sec;
  float left_r3 = (DCR3_HEIGHT - SLOPE * y3_rot);
  float right_r3 = (DCR3_HEIGHT + SLOPE * y3_rot);
  float radius2_DCr3 = DCR3_SQUARE - pow(y3_rot, 2);

  _fid_cut &= (x3_rot > left_r3);
  _fid_cut &= (x3_rot > right_r3);
  _fid_cut &= ((x3_rot * x3_rot) > radius2_DCr3);

  return _fid_cut;
}

bool Cuts::IsPip(int i) {
  if (_data->gpart() <= i) return false;
  bool _pip = true;
  _pip &= (_data->charge(i) == POSITIVE);
  _pip &= !(abs(_dt->dt_P(i)) < 0.5 || abs(_dt->dt_ctof_P(i)) < 0.2);
  return _pip;
}

bool Cuts::IsProton(int i) {
  if (_data->gpart() <= i) return false;
  bool _proton = true;
  _proton &= (_data->charge(i) == POSITIVE);
  //_proton &= _data->pid(i) == PROTON;
  //_proton = (_data->p(i) > 1.0);
  if (!std::isnan(_dt->dt_P(i))) _proton &= (abs(_dt->dt_P(i)) < 0.5);
  // if (!std::isnan(_dt->dt_ctof_P(i))) _proton &= (abs(_dt->dt_ctof_P(i)) < 0.2);
  return _proton;
}

bool Cuts::IsPim(int i) {
  if (_data->gpart() <= i) return false;
  bool _pim = true;
  _pim &= (_data->charge(i) == NEGATIVE);
  _pim &= (abs(_dt->dt_Pi(i)) < 0.5 || abs(_dt->dt_ctof_Pi(i)) < 0.2);
  return _pim;
}
