#include "config.h"
#include "analysis.h"

namespace ana {


Configuration::Configuration() 
    : _smdiph_reweight("_smdiph_reweight", "_smdiph_reweight", 13, 0, 1300)
{
    
    // kFactorVsMass_DiphoxNLO_mstw2008nlo_iso7GeV_16642_phD3PD000131_v_B_H4_iso_filt.root
    
    TH1D& h = _smdiph_reweight;
    h.SetBinContent(0, 0.77524548769); h.SetBinError(0, 0.193811371922);
    h.SetBinContent(1, 0.77524548769); h.SetBinError(1, 0.193811371922);
    h.SetBinContent(2, 1.61900544167); h.SetBinError(2, 0.404751360416);
    h.SetBinContent(3, 1.64430201054); h.SetBinError(3, 0.411075502634);
    h.SetBinContent(4, 1.49268329144); h.SetBinError(4, 0.373170822859);
    h.SetBinContent(5, 1.44546496868); h.SetBinError(5, 0.36136624217);
    h.SetBinContent(6, 1.31024003029); h.SetBinError(6, 0.327560007572);
    h.SetBinContent(7, 1.28618884087); h.SetBinError(7, 0.321547210217);
    h.SetBinContent(8, 1.16020214558); h.SetBinError(8, 0.290050536394);
    h.SetBinContent(9, 1.09157836437); h.SetBinError(9, 0.272894591093);
    h.SetBinContent(10, 1.06538069248); h.SetBinError(10, 0.26634517312);
    h.SetBinContent(11, 1.01491463184); h.SetBinError(11, 0.253728657961);
    h.SetBinContent(12, 0.994816064835); h.SetBinError(12, 0.248704016209);
    h.SetBinContent(13, 0.977248132229); h.SetBinError(13, 0.244312033057);
    h.SetBinContent(14, 1); h.SetBinError(14, 0.20);
}



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

