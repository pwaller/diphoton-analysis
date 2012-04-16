#ifndef _ANALYSIS_H_
#define _ANALYSIS_H_

#include <string>
#include <iostream>

#include <a4/processor.h>
#include <a4/atlas/ntup/photon/Event.pb.h>
#include <a4/atlas/Event.pb.h>

#include "config.h"

class Photon;

namespace ana {


class Configuration;

using a4::process::ProcessorOf;
namespace ntup = a4::atlas::ntup::photon;


class Analysis : public ProcessorOf<ntup::Event, a4::atlas::EventMetaData> {
protected:
    Configuration C;
    friend class Configuration;
    
    shared<EnergyRescaler> _rescaler;
    shared<Root::TPileupReweighting> _pileup_tool;
    
    bool _should_ptcut, _should_masscut, _simulation, _is_sm_diphoton_sample;
    double _pt_low, _pt_high, _mass_low, _mass_high;
    uint32_t _current_run, _current_sample;
    
    double _sum_mc_weights;
    uint64_t _event_count;
    
public:
    static Analysis* construct(const std::string& name, Configuration* c);

    Analysis(Configuration* c)
        : C(*c),
          _should_ptcut(false), _should_masscut(false),
          _simulation(false), _is_sm_diphoton_sample(false),
          _current_run(false), _current_sample(false),
          _sum_mc_weights(0), _event_count(0)
        
    {
        set_metadata_behavior(MANUAL_BACKWARD);
    }
    
    virtual void process(const ntup::Event& event);
    void process_end_metadata();
    
    // Called when a new mc_channel is encountered
    void new_sample(const ntup::Event& event);
    bool pass_grl(const ntup::Event& event);
    double compute_mass(const ntup::Event& event, const Photon& lead, const Photon& sublead) const;
    double resolution_plots(double mgg, double true_mgg);
    
    void get_smdiph_weight(const double mass_gev, double& w, double& err);

};

class Filter : public Analysis {
public:
    Filter(Configuration* c) : Analysis(c) {}
    void filter_photons(const ntup::Event& event, ntup::Event& new_event);
    void process(const ntup::Event& event);
};


inline Analysis* Analysis::construct(const std::string& name, Configuration* c) {
    if (name == "Analysis") return new Analysis(c);
    if (name == "Filter")   return new Filter(c);
    FATAL("Unknown processor: '", name, "'");
}


}

#endif

