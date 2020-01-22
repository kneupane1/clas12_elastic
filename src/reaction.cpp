/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/
#include "reaction.hpp"

Reaction::Reaction(const std::shared_ptr<Branches12>& data, float beam_energy) {
  _data = data;
  _beam = std::make_unique<TLorentzVector>();
  _beam_energy = beam_energy;

  _beam->SetPxPyPzE(0.0, 0.0, 10.587, 10.6);  // sqrt(_beam_energy * _beam_energy - MASS_E * MASS_E), _beam_energy);
  _gamma = std::make_unique<TLorentzVector>();
  _target = std::make_unique<TLorentzVector>(0.0, 0.0, 0.0, MASS_P);
  _elec = std::make_unique<TLorentzVector>();
  //  _x_mu = std::make_unique<TLorentzVector>();

  this->SetElec();
}

Reaction::~Reaction() {}

void Reaction::SetElec() {
  _hasE = true;
  _elec->SetXYZM(_data->px(0), _data->py(0), _data->pz(0), MASS_E);

  *_gamma += *_beam - *_elec;
  // Can calculate W and Q2 here
  _W = physics::W_calc(*_beam, *_elec);
  _Q2 = physics::Q2_calc(*_beam, *_elec);
  _beam_theta = _beam->Theta() * RAD2DEG;
  _elec_theta = _elec->Theta() * RAD2DEG;
}

void Reaction::SetPositive(int i) {
  _pos_beta.push_back(_data->beta(i));
  _numPos++;
  _hasPos = true;
  _pos_det.push_back(abs(_data->status(i) / 1000));
  _pos.push_back(std::make_unique<TLorentzVector>());
  _pos.back()->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_P);
}

bool Reaction::PosStats() {
  if (_pos.size() != 2) return false;
  if (abs(_pos.front()->Phi() - _pos.back()->Phi()) > 0.5) return false;
  if (abs(_pos.front()->Theta() - _pos.back()->Theta()) > 0.5) return false;
  if (_pos_det.front() == _pos_det.back()) return false;
  /*
  std::cout << "num_pos: " << _pos.size() << std::endl;
  std::cout << "num_phi: ";
  for (auto& _d : _pos_det) std::cout << detector_name[_d] << "\t";
  std::cout << std::endl;
  std::cout << "pos_phi: ";
  for (auto& _p : _pos) std::cout << _p->Phi() << "\t";
  std::cout << std::endl;
  std::cout << "pos_theta: ";
  for (auto& _p : _pos) std::cout << _p->Theta() << "\t";
  std::cout << std::endl;
  */

  return true;
}

void Reaction::SetOther(int i) {
  _numOther++;
  _hasOther = true;
  if (_data->charge(i) != 0) {
    _other.push_back(std::make_unique<TLorentzVector>());
    _other.back()->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), mass[_data->pid(i)]);
  }
  if (_data->pid(i) == PHOTON) {
    _photons.push_back(std::make_unique<TLorentzVector>());
    _photons.back()->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), 0);
  }
}

void Reaction::CalcMissMass() {
  if (_pos.size() > 0) {
    auto _x_mu = std::make_unique<TLorentzVector>();
    //      auto mm = std::make_unique<TLorentzVector>();
    //  *mm += (*_gamma + *_target);

    *_x_mu += (*_gamma + *_target);
    // for (auto& _p : _photons) *mm -= *_p;
    for (auto& _p : _pos) *_x_mu -= *_p;
    _MM = _x_mu->M();
    _MM2 = _x_mu->M2();

    _x_mu_E = _x_mu->E();
    _x_mu_P = _x_mu->P();
    _x_mu_Px = _x_mu->Px();
    _x_mu_Py = _x_mu->Py();
    _x_mu_Pz = _x_mu->Pz();

    _x_mu_theta = _x_mu->Theta() * RAD2DEG;

    _x_mu_m2 = _x_mu->M2();
    _x_mu_m = _x_mu->M();

    //  _x_mu->SetPxPyPzE(mm->Px(), mm->Py(), mm->Pz(), mm->E());
    //(*_gamma + *_target - *_positive);
  }
}

float Reaction::MM() {
  if (_MM != _MM) CalcMissMass();
  return _MM;
}
float Reaction::MM2() {
  if (_MM2 != _MM2) CalcMissMass();

  return _MM2;
}
float Reaction::M_x_mu() {
  if (_x_mu_m != _x_mu_m) CalcMissMass();

  return _x_mu_m;
}
float Reaction::M2_x_mu() {
  if (_x_mu_m2 != _x_mu_m2) CalcMissMass();

  return _x_mu_m2;
}
float Reaction::Px_x_mu() {
  if (_x_mu_Px != _x_mu_Px) CalcMissMass();

  return _x_mu_Px;
}
float Reaction::Py_x_mu() {
  if (_x_mu_Py != _x_mu_Py) CalcMissMass();

  return _x_mu_Py;
}
float Reaction::Pz_x_mu() {
  if (_x_mu_Pz != _x_mu_Pz) CalcMissMass();

  return _x_mu_Pz;
}
float Reaction::P_x_mu() {
  if (_x_mu_P != _x_mu_P) CalcMissMass();

  return _x_mu_P;
}
float Reaction::E_x_mu() {
  if (_x_mu_E != _x_mu_E) CalcMissMass();
  //  std::cout << "_x_mu_p  " << _x_mu->E() << '\n';
  //  if (_x_mu_E > 0)
  return _x_mu_E;
  // else
  // return NAN;
}
float Reaction::theta_x_mu() {
  if (_x_mu_theta != _x_mu_theta) CalcMissMass();

  return _x_mu_theta;
}
float Reaction::theta_beam() { return _beam_theta; }
float Reaction::theta_elec() { return _elec_theta; }

void Reaction::CalcMassPi0() {
  _pi0_mass = 0;
  if (_photons.size() == 2) {
    auto mass = std::make_unique<TLorentzVector>();
    for (auto& _p : _photons) *mass -= *_p;
    _pi0_mass = mass->M();
  }
}

float Reaction::pi0_mass() {
  if (_pi0_mass != _pi0_mass) CalcMassPi0();
  return _pi0_mass;
}
