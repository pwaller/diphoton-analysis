#ifndef _ANALYSIS_H_
#define _ANALYSIS_H_

#include <string>
#include <iostream>

#include <a4/processor.h>
#include <a4/atlas/ntup/photon/Event.pb.h>

#include "config.h"

class Photon;

namespace ana {


class Configuration;

using a4::process::ProcessorOf;
namespace ntup = a4::atlas::ntup::photon;


class Analysis : public ProcessorOf<ntup::Event> {
protected:
    Configuration C;
    friend class Configuration;
    
    shared<EnergyRescaler> _rescaler;
    shared<Root::TPileupReweighting> _pileup_tool;
    
    bool _should_ptcut, _should_masscut, _simulation;
    double _pt_low, _pt_high, _mass_low, _mass_high;
    uint32_t _current_run, _current_sample;
    
public:
    static Analysis* construct(const std::string& name, Configuration* c);

    Analysis(Configuration* c) : C(*c) {}
    void process(const ntup::Event& event);
    
    // Called when a new mc_channel is encountered
    void new_sample(const ntup::Event& event);
    bool pass_grl(const ntup::Event& event);
    double compute_mass(const ntup::Event& event, const Photon& lead, const Photon& sublead) const;
    double resolution_plots(double mgg, double true_mgg);

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

