/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/

#ifndef MAIN_H_GUARD
#define MAIN_H_GUARD

#include <iostream>
#include "TFile.h"
#include "TH1.h"
#include "branches.hpp"
#include "colors.hpp"
#include "cuts.hpp"
#include "histogram.hpp"
#include "reaction.hpp"

size_t run(std::shared_ptr<TChain> _chain, const std::shared_ptr<Histogram> &_hists, int thread_id);
size_t run_files(std::vector<std::string> inputs, const std::shared_ptr<Histogram> &hists, int thread_id);

size_t run_files(std::vector<std::string> inputs, const std::shared_ptr<Histogram> &hists, int thread_id) {
  // Called once for each thread
  // Make a new chain to process for this thread
  auto chain = std::make_shared<TChain>("clas12");
  // Add every file to the chain
  for (auto in : inputs) chain->Add(in.c_str());

  // Run the function over each thread
  return run(chain, hists, thread_id);
}

size_t run(std::shared_ptr<TChain> _chain, const std::shared_ptr<Histogram> &_hists, int thread_id) {
  // Get the number of events in this thread
  size_t num_of_events = (int)_chain->GetEntries();
  float beam_energy = NAN;
  if (getenv("BEAM_E") != NULL) beam_energy = atof(getenv("BEAM_E"));

  // Print some information for each thread
  std::cout << "=============== " << RED << "Thread " << thread_id << DEF << " =============== " << BLUE
            << num_of_events << " Events " << DEF << "===============\n";

  // Make a data object which all the branches can be accessed from
  auto data = std::make_shared<Branches12>(_chain);

  // Total number of events "Processed"
  size_t total = 0;
  // For each event
  for (size_t current_event = 0; current_event < num_of_events; current_event++) {
    // Get current event
    _chain->GetEntry(current_event);
    // If we are the 0th thread print the progress of the thread every 1000 events
    if (thread_id == 0 && current_event % 10000 == 0)
      std::cerr << "\t" << (100 * current_event / num_of_events) << " %\r" << std::flush;

    if (data->charge(0) == NEGATIVE) _hists->FillHists_electron_cuts(data);

    auto dt = std::make_shared<Delta_T>(data);
    auto cuts = std::make_shared<Cuts>(data, dt);
    if (!cuts->ElectronCuts()) continue;
    _hists->Fill_SF(data);
    _hists->FillHists_electron_with_cuts(data);

    _hists->Fill_Dt(dt);

    // Make a reaction class from the data given
    auto event = std::make_shared<Reaction>(data, beam_energy);
    // For each particle in the event
    for (int part = 1; part < data->gpart(); part++) {
      if (cuts->IsProton(part)) {
        event->SetPositive(part);
        _hists->Fill_Dt(dt, part);
      } else {
        event->SetOther(part);
      }
    }
    _hists->Fill_pi0(event);
    if (event->onePositive_at180_MM0()) {
      total++;
      _hists->Fill_Sparce(event);
    }
    _hists->Fill_WvsQ2(event);
    _hists->Fill_MomVsBeta(event);
  }
  std::cout << "Percent = " << 100.0 * total / num_of_events << std::endl;
  // Return the total number of events
  return num_of_events;
}
#endif

/*
ep -> e x+

W for all events
W for 2 particles
W for 2 Part 2nd positive
hist Phi_e - Phi_pos ~ 90
W for cut around 90
pos_mom vs pos_theta
theta_p_pos_calc_from_electron - theta_pos_measured

calc theta from magnitude of pos momentum

*/
