#include "config.h"
#include "analysis.h"

namespace ana {

a4::process::Processor* Configuration::new_processor() {
    return Analysis::construct(_processor, this);
}


void Configuration::setup_processor(a4::process::Processor& _g) {
    auto& g = static_cast<ana::Analysis&>(_g);
    
    g._rescaler.reset(new EnergyRescaler());
    g._rescaler->useDefaultCalibConstants("2011");
    
    if (_do_pileup_reweighting)
        g._pileup_tool.reset(get_prw(_pileup_mc_file.c_str(), _pileup_data_file.c_str()));
}


}

