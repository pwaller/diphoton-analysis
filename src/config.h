#ifndef _CONFIG_H_
#define _CONFIG_H_

#include <a4/application.h>
#include <a4/grl.h>
using a4::atlas::GRL;
using a4::atlas::FileGRL;
using a4::atlas::NoGRL;

#include "event_list.h"
#include "external.h"


namespace ana {


class Configuration : public a4::process::Configuration {
public:
    std::string _grl_name, _processor;
    
    shared<GRL> _grl;
    shared<EventList> _ee_events;
    
    std::string _pileup_mc_file, _pileup_data_file, _ee_event_file;
    bool _do_pileup_reweighting,
         _do_plot,
         _do_sf_reweighting,
         _require_mc_match,
         _write_anatree,
         _filter_reco_photons;

    virtual void add_options(po::options_description_easy_init opt) {
        opt("grl", po::value(&_grl_name), "GRL file");
        opt("processor,P", po::value(&_processor)->default_value("Analysis"), "Processor");
        opt("pileup-mc", po::value<std::string>(&_pileup_mc_file)->default_value("pileup_reweighting.root"), "File with MC distributions for Pileup reweighting");
        opt("pileup-data", po::value<std::string>(&_pileup_data_file)->default_value("pileup_data.root"), "File with data distributions for Pileup reweighting");
        opt("rw-pileup", po::bool_switch(&_do_pileup_reweighting)->default_value(false), "Reweight Pileup");
        opt("rw-scalefactor", po::bool_switch(&_do_sf_reweighting)->default_value(false), "Scale factor reweighting");
        opt("do-plot", po::bool_switch(&_do_plot)->default_value(false), "Make plots");
        opt("require-mc-match", po::bool_switch(&_require_mc_match)->default_value(false), "Only pass photons which pass the hard process match");
        opt("write-anatree", po::bool_switch(&_write_anatree)->default_value(false), "Write analysis tree with corrected photons");
        opt("ee-event-file", po::value(&_ee_event_file), "Filename of list of events to exclude for ee cut");
        opt("filter-reco-ph", po::bool_switch(&_filter_reco_photons)->default_value(false), "Filter reconstructed photons");
    }

    virtual void read_arguments(po::variables_map& arguments) {
        if (arguments.count("grl")) 
            _grl.reset(new FileGRL(arguments["grl"].as<std::string>()));
        else 
            _grl.reset(new NoGRL());
            
        if (_ee_event_file != "")
            _ee_events.reset(new EventList(_ee_event_file));
    }

    void setup_processor(a4::process::Processor&);
    a4::process::Processor* new_processor();
};


}

#endif
