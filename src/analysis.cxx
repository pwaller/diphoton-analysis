#include "all.h"

#include <iostream>
#include <unordered_map>
#include <math.h>

#include <a4/application.h>
using a4::store::ObjectStore;

#include <a4/histogram.h>
using a4::hist::H1;
using a4::hist::H2;
using a4::hist::Cutflow;
#include <a4/axis.h>
using a4::hist::VariableAxis;

#include <a4/utility.h>
using a4::process::utility::vector_of_ptr;
using a4::process::utility::in_map;

#include "analysis.h"
#include "config.h"

#include "constants.h"
#include "event_view.h"

//using a4::atlas::ntup::photon::Event;
#include <a4/atlas/ntup/photon/Event.pb.h>

namespace ntup = a4::atlas::ntup::photon;

namespace ana {


const std::unordered_map<int, ResonanceSample> resonance_samples = {
    { 105324, {105324, "105324/", 1250, 0.05, 4.725, 9.08 } },
    { 105623, {105623, "105623/", 500, 0.01, 0.075, 82.5 } },
    { 105833, {105833, "105833/", 800, 0.03, 1.088, 54.4 } },
    { 105834, {105834, "105834/", 1000, 0.03, 1.361, 13.9 } },
    { 105835, {105835, "105835/", 700, 0.05, 2.646, 329.5 } },
    { 105836, {105836, "105836/", 1000, 0.05, 3.78, 39.0 } },
    { 105837, {105837, "105837/", 1500, 0.05, 4.725, 2.64 } },
    { 105838, {105838, "105838/", 800, 0.1, 12.09, 600.9 } },
    { 105839, {105839, "105839/", 1000, 0.1, 15.12, 152.6 } },
    { 105841, {105841, "105841/", 1250, 0.1, 18.9, 36.1 } },
    { 105842, {105842, "105842/", 1250, 0.1, 18.9, 36.1 } },
    { 106623, {106623, "106623/", 800, 0.01, 0.12, 6.0 } },
    { 106643, {106643, "106643/", 1000, 0.01, 0.151, 1.56 } },
    { 106644, {106644, "106644/", 500, 0.03, 0.68, 741.8 } },
    { 106684, {106684, "106684/", 300, 0.01, 0.045, 1052.0 } },
    { 115557, {115557, "115557/", 700, 0.01, 0.089, 12.98 } },
    { 115558, {115558, "115558/", 700, 0.03, 0.871, 116.3 } },
    { 115559, {115559, "115559/", 800, 0.05, 2.79, 150.1 } },
    { 115560, {115560, "115560/", 900, 0.03, 1.15, 26.94 } },
    { 115561, {115561, "115561/", 900, 0.05, 2.98, 74.99 } },
    { 115562, {115562, "115562/", 900, 0.07, 5.67, 142.8 } },
    { 115563, {115563, "115563/", 900, 0.1, 12.06, 293.1 } },
    { 115564, {115564, "115564/", 1100, 0.05, 3.75, 20.99 } },
    { 115565, {115565, "115565/", 1100, 0.07, 7.44, 41.15 } },
    { 115566, {115566, "115566/", 1100, 0.1, 14.93, 83.46 } },
    { 115567, {115567, "115567/", 1100, 0.2, 57.6, 326.5 } },
    { 115568, {115568, "115568/", 1500, 0.1, 8.23, 10.5 } },
    { 115569, {115569, "115569/", 1250, 0.15, 37.3, 81.07 } },
    { 115570, {115570, "115570/", 1250, 0.2, 65.2, 142.2 } },
    { 119870, {119870, "119870/", 1750, 0.1, 26.46, 3.44 } },
    { 119871, {119871, "119871/", 2000, 0.1, 30.24, 1.21 } },
    { 119872, {119872, "119872/", 2250, 0.1, 34.02, 0.46 } },
};

const ResonanceSample* get_resonance(const ntup::Event& e) {
    if (e.has_mc_channel_number() == 0) return NULL;
    const auto& i = resonance_samples.find(e.mc_channel_number());
    if (i != resonance_samples.end())
        return &i->second;
    return NULL;
}

// Why!?
// Gives 10^-12 compatibility with reducible template
const double low_mass_edge = 409.40104882461151;
VariableAxis mass_logbins       (VariableAxis::log_bins(53, low_mass_edge, 3000));
VariableAxis mass_logbins_full  (VariableAxis::log_bins(200, 10, 3000));

void Analysis::process_end_metadata() {
    auto m = metadata();
    auto& processing_step = *m.add_processing_steps();
    processing_step.set_name("pwanalysis");
    
    m.set_sum_mc_weights(_sum_mc_weights); _sum_mc_weights = 0;
    m.set_event_count(_event_count); _event_count = 0;
    
    metadata_end_block(m);
}

inline void Analysis::get_smdiph_weight(const double mass_gev,
                                         double& w, double& err) {
    const auto bin = C._smdiph_reweight.FindBin(mass_gev);
    w = C._smdiph_reweight.GetBinContent(bin);
    err = C._smdiph_reweight.GetBinError(bin);
}


inline void Analysis::make_resonance_plots(
    ObjectStore D, const ntup::Event& event,
    const Photon& lead, const Photon& sublead,
    const double mgg, const double mgg_true, 
    double mass, double km, double width, double cross_section) 
{    
    const double resolution = 1.208688 + 0.009822771*mgg_true;
    const double resol_width = hypot(width, resolution);
    
    D.T<H1>("mgg")     (100, mass - resol_width*5, mass + resol_width*5).fill(mgg / 1000);
    D.T<H1>("mgg_true")(100, mass - width*5,       mass + width*5      ).fill(mgg_true / 1000);
    
}

inline void Analysis::make_resonances_plots(
    const uint32_t number,
    const ntup::Event& event,
    const Photon& lead, const Photon& sublead,
    const double mgg, const double mgg_true) {
    
    switch (number) {
        case 105324: make_resonance_plots(S("resonances/105324/"), event, lead, sublead, mgg, mgg_true, 1250, 0.05, 4.725, 9.08); break;
        case 105623: make_resonance_plots(S("resonances/105623/"), event, lead, sublead, mgg, mgg_true, 500, 0.01, 0.075, 82.5); break;
        case 105833: make_resonance_plots(S("resonances/105833/"), event, lead, sublead, mgg, mgg_true, 800, 0.03, 1.088, 54.4); break;
        case 105834: make_resonance_plots(S("resonances/105834/"), event, lead, sublead, mgg, mgg_true, 1000, 0.03, 1.361, 13.9); break;
        case 105835: make_resonance_plots(S("resonances/105835/"), event, lead, sublead, mgg, mgg_true, 700, 0.05, 2.646, 329.5); break;
        case 105836: make_resonance_plots(S("resonances/105836/"), event, lead, sublead, mgg, mgg_true, 1000, 0.05, 3.78, 39.0); break;
        case 105837: make_resonance_plots(S("resonances/105837/"), event, lead, sublead, mgg, mgg_true, 1500, 0.05, 4.725, 2.64); break;
        case 105838: make_resonance_plots(S("resonances/105838/"), event, lead, sublead, mgg, mgg_true, 800, 0.1, 12.09, 600.9); break;
        case 105839: make_resonance_plots(S("resonances/105839/"), event, lead, sublead, mgg, mgg_true, 1000, 0.1, 15.12, 152.6); break;
        case 105841: make_resonance_plots(S("resonances/105841/"), event, lead, sublead, mgg, mgg_true, 1250, 0.1, 18.9, 36.1); break;
        case 106623: make_resonance_plots(S("resonances/106623/"), event, lead, sublead, mgg, mgg_true, 800, 0.01, 0.12, 6.0); break;
        case 106643: make_resonance_plots(S("resonances/106643/"), event, lead, sublead, mgg, mgg_true, 1000, 0.01, 0.151, 1.56); break;
        case 106644: make_resonance_plots(S("resonances/106644/"), event, lead, sublead, mgg, mgg_true, 500, 0.03, 0.68, 741.8); break;
        case 106684: make_resonance_plots(S("resonances/106684/"), event, lead, sublead, mgg, mgg_true, 300, 0.01, 0.045, 1052.0); break;
        case 115557: make_resonance_plots(S("resonances/115557/"), event, lead, sublead, mgg, mgg_true, 700, 0.01, 0.089, 12.98); break;
        case 115558: make_resonance_plots(S("resonances/115558/"), event, lead, sublead, mgg, mgg_true, 700, 0.03, 0.871, 116.3); break;
        case 115559: make_resonance_plots(S("resonances/115559/"), event, lead, sublead, mgg, mgg_true, 800, 0.05, 2.79, 150.1); break;
        case 115560: make_resonance_plots(S("resonances/115560/"), event, lead, sublead, mgg, mgg_true, 900, 0.03, 1.15, 26.94); break;
        case 115561: make_resonance_plots(S("resonances/115561/"), event, lead, sublead, mgg, mgg_true, 900, 0.05, 2.98, 74.99); break;
        case 115562: make_resonance_plots(S("resonances/115562/"), event, lead, sublead, mgg, mgg_true, 900, 0.07, 5.67, 142.8); break;
        case 115563: make_resonance_plots(S("resonances/115563/"), event, lead, sublead, mgg, mgg_true, 900, 0.1, 12.06, 293.1); break;
        case 115564: make_resonance_plots(S("resonances/115564/"), event, lead, sublead, mgg, mgg_true, 1100, 0.05, 3.75, 20.99); break;
        case 115565: make_resonance_plots(S("resonances/115565/"), event, lead, sublead, mgg, mgg_true, 1100, 0.07, 7.44, 41.15); break;
        case 115566: make_resonance_plots(S("resonances/115566/"), event, lead, sublead, mgg, mgg_true, 1100, 0.1, 14.93, 83.46); break;
        case 115567: make_resonance_plots(S("resonances/115567/"), event, lead, sublead, mgg, mgg_true, 1100, 0.2, 57.6, 326.5); break;
        case 115568: make_resonance_plots(S("resonances/115568/"), event, lead, sublead, mgg, mgg_true, 1500, 0.1, 8.23, 10.5); break;
        case 115569: make_resonance_plots(S("resonances/115569/"), event, lead, sublead, mgg, mgg_true, 1250, 0.15, 37.3, 81.07); break;
        case 115570: make_resonance_plots(S("resonances/115570/"), event, lead, sublead, mgg, mgg_true, 1250, 0.2, 65.2, 142.2); break;
        case 119870: make_resonance_plots(S("resonances/119870/"), event, lead, sublead, mgg, mgg_true, 1750, 0.1, 26.46, 3.44); break;
        case 119871: make_resonance_plots(S("resonances/119871/"), event, lead, sublead, mgg, mgg_true, 2000, 0.1, 30.24, 1.21); break;
        case 119872: make_resonance_plots(S("resonances/119872/"), event, lead, sublead, mgg, mgg_true, 2250, 0.1, 34.02, 0.46); break;
        
        case 145536:
        {
            // Template sample
            const int samples[] = {105324, 105623, 105833, 105834, 105835, 
                105836, 105837, 105838, 105839, 105841, 106623, 106643, 106644,
                106684, 115557, 115558, 115559, 115560, 115561, 115562, 115563, 
                115564, 115565, 115566, 115567, 115568, 115569, 115570, 119870, 
                119871, 119872};
            foreach (const int&s, samples) {
                // TODO(pwaller): Reweighting!
                make_resonances_plots(s, event, lead, sublead, mgg, mgg_true);
            }
            break;
        }
        
        default:
            break;
    }
}

void Analysis::process(const ntup::Event& event) {
    if (event.run_number() != _current_run || 
        event.mc_channel_number() != _current_sample)
        new_sample(event);
    
    #define PASSED(x) \
        S.T<Cutflow>("cutflow").passed(x);
        
    #define EFFPLOT(name) \
        S.T<H1>("eff/mgg_true/" name)(7000, 0, 7e6).fill(mgg_true); \
        
    #define EFFPLOT_1(name) \
    do { \
        if (is_mc) { \
            EFFPLOT(name); \
            if (phtr_1 && phtr_2) { \
                S.T<H1>("eff/true_lead/pt/" name)(3000, 0, 3e6).fill(phtr_1->pt()); \
                S.T<H1>("eff/true_lead/eta/" name)(100, -2.5, 2.5).fill(phtr_1->eta()); \
                S.T<H1>("eff/true_sublead/pt/" name)(3000, 0, 3e6).fill(phtr_2->pt()); \
                S.T<H1>("eff/true_sublead/eta/" name)(100, -2.5, 2.5).fill(phtr_2->eta()); \
            } \
        } \
    } while(false)
        
        
    bool is_mc = event.photon_truth_particles_size() != 0;
    bool data = !is_mc;
    
    S.set_weight(1);

    if (event.run_number() != _current_run ||
        event.mc_channel_number() != _current_sample)
        new_sample(event);
        
    uint32_t event_number = event.event_number();
    uint32_t run_number = event.run_number();
    
    if (is_mc)
        S.mul_weight(event.mc_event_weight());
    
    if (rerun_systematics_current == NULL) {
        _event_count += 1;
        _sum_mc_weights += S.weight();
    }
    
    // MC: Pileup reweighting and run number reassignment
    float pileup_weight = 1.0;
    if (!data && C._do_pileup_reweighting && !systematic("noprw")) {
        _pileup_tool->SetRandomSeed(event_number+1);
        pileup_weight = _pileup_tool->GetCombinedWeight(
            event.run_number(),
            event.mc_channel_number(),
            event.averageintperxing());
        //pileup_tool->SetRandomSeed(314159 + event.mc_channel_number()*2718 + event_number);
        //run_number = pileup_tool->GetRandomRunNumber(run_number);
    }
    
    S.T<H1>("pileup_weight")(100, 0, 2).fill(pileup_weight);
    
    if (!data) {
        S.mul_weight(pileup_weight);
    }
    
    //auto hard_process_photons = vector_of<ana::TruePhoton>(event.photon_truth_particles());
    auto hard_process_photons = vector_of_ptr(event.photon_truth_particles());
    REMOVE_IF(hard_process_photons, ph, !ph->ishardprocphoton());
    SORT_KEY (hard_process_photons, ph, - ph->pt());
    
    const ntup::PhotonTruthParticle *phtr_1 = NULL, *phtr_2 = NULL;
    
    ALorentzVector phtr_1_lv, phtr_2_lv;
                   
    if (is_mc && hard_process_photons.size() >= 2) {
        phtr_1 = hard_process_photons.at(0);
        phtr_2 = hard_process_photons.at(1);
        phtr_1_lv = ALorentzVector::from_ptetaphim(*phtr_1);
        phtr_2_lv = ALorentzVector::from_ptetaphim(*phtr_2);
    }
    
    auto gravitons = vector_of_ptr(event.photon_truth_particles());
    REMOVE_IF(gravitons, ph, ph->pdgid() != 5000039);
    double mgg_true = -999;
    if (gravitons.size())
        mgg_true = gravitons[0]->m();
    else if (is_mc && hard_process_photons.size() >= 2) {
        mgg_true = (phtr_1_lv + phtr_2_lv).m();
    }
    
    if (_is_sm_diphoton_sample) {
        double w = 1, err = 0;
        get_smdiph_weight(mgg_true / 1000., w, err);
        if (systematic("kfac_up"))
            w += err;
        else if (systematic("kfac_down"))
            w -= err;
        if (systematic("kfac_off"))
            w = 1;
        S.mul_weight(w);
    }
    
    
    //S.T<H1>("weight")(100, 0, 3, "weight").fill(S.weight());
    
    PASSED("Total");
    
    if (is_mc && _should_masscut)
        if (mgg_true < _mass_low || mgg_true >= _mass_high)
            return;
    
    if (is_mc && C._require_mc_match && hard_process_photons.size() < 2)
        return;
    
    PASSED("Good MC");
    
    EFFPLOT("0_total");
    
    if (!event.ef()._2g20_loose()) return;
    PASSED("2g20_loose");
    EFFPLOT("1_trigger");
    
    if (!pass_grl(event)) return;
    PASSED("GRL");
    
    bool good_vertex = false;
    foreach (auto& vertex, event.primary_vertices())
        if (vertex.ntracks() >= 3) {
            good_vertex = true;
            break;
        }
        
    if (!good_vertex) return;
    PASSED("PV");
    
    // This is here so that we can get the single-photon reco-efficiency
    EFFPLOT_1("2_pv");
    
    auto true_photon = [&](const ntup::Photon* reco_photon)
                        -> const ntup::PhotonTruthParticle* {
        const int index = reco_photon->truth_index();
        if (index < 0) return NULL;
        const auto size = event.photon_truth_particles_size();
        if (index >= size) {
            DEBUG("Requested skimmed photon :-( ", index, " ", reco_photon->pt(), " ", size);
            return NULL;
        }
        return &event.photon_truth_particles(index);
    };
    
    {
        // Only compute quantities here which are missing
        ntup::Event& mutable_event = const_cast<ntup::Event&>(event);
        foreach (auto& mutable_ph, *mutable_event.mutable_photons())
            Photon::compute_extra_quantities(mutable_ph);
    }

    auto good_photons = Photon::make_vector(event.photons());
    foreach_enumerate (i, auto& ph, good_photons) {
        auto index = ph->has_original_index() ? ph->original_index() : i;
        ph.compute_corrections(event, index, *_rescaler);
    }
    SORT_KEY (good_photons, ph, -ph->pt());
    
    if (event.photons_size() < 2) return;
    Photon ph_1, ph_2;
    
    if (is_mc) {
        foreach (auto ph, good_photons) {
            const auto* phtr = true_photon(ph);
            if (!phtr)
                continue;
            if (phtr == phtr_1) {
                // If we already have a match, take the one with the better matching pt.
                if (!ph_1 || abs(ph_1->pt() - phtr->pt()) > abs(ph->pt() - phtr->pt())) {
                    ph_1 = ph;
                }
            } else if (phtr == phtr_2) {
                if (!ph_2 || abs(ph_2->pt() - phtr->pt()) > abs(ph->pt() - phtr->pt())) {
                    ph_2 = ph;
                }
            }
        }
        if (C._require_mc_match && !(ph_1 && ph_2))
            return;
    }
    
    PASSED("Reco 2#gamma");
    EFFPLOT_1("3_reco");
    
    S.T<H1>("actualintperxing")(120, 0, 30, "#mu").fill(event.actualintperxing());
    S.T<H1>("averageintperxing")(120, 0, 30, "#mu").fill(event.averageintperxing());
    
    #define CUT(name, o, cut) \
        REMOVE_IF(good_photons, o, cut); \
        if (good_photons.size() < 2) return; \
        { \
            if (is_mc && C._require_mc_match) { \
                auto cutfunc = [&](const Photon& o) { return cut; }; \
                if (cutfunc(ph_1) || cutfunc(ph_2)) \
                    return; \
                EFFPLOT_1(name); \
            } else if (is_mc) { \
                EFFPLOT_1(name); \
            } \
        } \
        PASSED(name);
    //plot_boson(S("cut/" name "/"), phtr_1_lv, phtr_2_lv);
    
    CUT("4_pt", ph, ph.corrected().pt() <= 25e3);
    
    CUT("5_eta", ph, abs(ph->etas2()) >= 1.37 && abs(ph->etas2()) <= 1.52
                     || abs(ph->etas2()) >= 2.37);
    
    CUT("6_oq", ph, ph->oq() & OQ_BAD_BITS);
    
    CUT("7_phclean", ph, ((ph->oq() & LARBITS_PHOTON_CLEANING) != 0
                          && (ph->reta() > 0.98
                              || ph->rphi() > 1.0
                              || ((ph->oq() & LARBITS_OUTOFTIME_CLUSTER) != 0)
                            )
                         ));
    
    CUT("8_loose", ph, !ph.corrected().loose());
            
    SORT_KEY (good_photons, ph, -ph.corrected().pt());
    auto lead = good_photons[0], sublead = good_photons[1];
    
    if (is_mc && C._do_sf_reweighting) {
        S.mul_weight(lead.scale_factor());
        S.mul_weight(sublead.scale_factor());
    }
    
    PASSED("preselection");
    
    if (!lead.corrected().tight() || !sublead.corrected().tight())
        return;
        
    CUT("9_tight", ph, !ph.corrected().tight());
    
    //plot_boson(S("2_tight/"), *lead, *sublead);
    
    auto isolated = [&](const Photon& ph) {
        return ph.corrected().analysis_isolation() < 5000.;
    };
    
    if (!isolated(lead) || !isolated(sublead))
        return;
    
    CUT("10_iso", ph, !isolated(ph));
    
    const double mgg = compute_mass(event, lead, sublead);
    
    S.T<H1>("sel_reco_mgg")(7000, 0, 7e3, "m_{#gamma#gamma} [GeV]").fill(mgg / 1000);
    S.T<H1>("sel_reco_mgg_log")
        .with_axis(mass_logbins, "m_{#gamma#gamma} [GeV]")
        .fill(mgg / 1000);
    S.T<H1>("sel_reco_mgg_log_full")
        .with_axis(mass_logbins_full, "m_{#gamma#gamma} [GeV]")
        .fill(mgg / 1000);
    
    if (is_mc)
        resolution_plots(mgg, mgg_true);
        
    if (event.larerror() > 1) return;
    PASSED("LarOK");
    
    if (C._ee_events)
        if (C._ee_events->present(event.run_number(), event.event_number()))
            return;
    PASSED("!ee");
    
    if (mgg < 140e3) return;
    PASSED("mgg140");
        
    EFFPLOT_1("11_mass");
    
    make_resonances_plots(event.mc_channel_number(), event,
                          lead, sublead, mgg, mgg_true);
}

void Analysis::new_sample(const ntup::Event& event) {
    _current_run = event.run_number();
    _current_sample = event.mc_channel_number();
    _simulation = event.issimulation();
    
    _current_resonance = get_resonance(event);
    
    switch (_current_sample) {
        // Samples:
        //       JF:    PhotonJet:           DP:
        case 105802: case 108087: case 115802:
            _should_ptcut = true; _pt_low = 20000; _pt_high = 45000;
            break; // e.g. JF17
        
        case 105807: case 108081: case 115803:
            _should_ptcut = true; _pt_low = 45000; _pt_high = 85000;
            break; // e.g. JF35
        
        case 105814: case 108082: case 115804:
            _should_ptcut = true; _pt_low = 85000; _pt_high = 150000;
            break; // e.g. JF70
        
        case 105815: case 108083: case 115809:
            _should_ptcut = true; _pt_low = 150000; _pt_high = 260000;
            break; // e.g. JF140
        
        case 105812: case 108084: case 115810:
            _should_ptcut = true; _pt_low = 260000; _pt_high = 550000;
            break; // e.g. JF240
        
        case 119075: case 119079: case 119077:
            _should_ptcut = true; _pt_low = 550000; _pt_high = 9999999999999;
            break; // e.g. JF500
        default:
            _should_ptcut = false;
    }
    
    switch (_current_sample) {
        case 119584: // mc11_7TeV.119584.Pythiagamgam15_highmass
            _is_sm_diphoton_sample = _should_masscut = true;
            _mass_low = 0;      _mass_high = 800e3;
            true;
            break;
        case 145606: // mc11_7TeV.145606.Pythiagamgam15_M_gt_800
            _is_sm_diphoton_sample = _should_masscut = true;
            _mass_low = 800e3;  _mass_high = 1500e3;
            break;
        case 145607: // mc11_7TeV.145607.Pythiagamgam15_M_gt_1500
            _is_sm_diphoton_sample = _should_masscut = true;
            _mass_low = 1500e3; _mass_high = 9999999999999;
            break;
        default:
            _is_sm_diphoton_sample = _should_masscut = false;
    }
}

bool Analysis::pass_grl(const ntup::Event& event) {
    if (event.issimulation()) return true;
    return C._grl->pass(event.run_number(), event.lbn());
}
    
double Analysis::compute_mass(const ntup::Event& event, const Photon& lead, const Photon& sublead) const {
    double E1 = lead.corrected().e();
    double eta1 = lead->etas1();
    assert(lead->has_etas1());
    double phi1 = lead->phi();
    
    double E2 = sublead.corrected().e();
    double eta2 = sublead->etas1();
    assert(sublead->has_etas1());
    double phi2 = sublead->phi();
    
    double PV_ID = event.primary_vertices(0).z();
    double mgg = GetCorrectedInvMass(E1, eta1, phi1, E2, eta2, phi2, PV_ID);
    return mgg;
}

double Analysis::resolution_plots(double mgg, double mgg_true) {
    S.T<H1>("sel_mgg_true")(7000, 0, 7e3, "m_{#gamma#gamma} [GeV]").fill(mgg_true/1000.);
    
    S.T<H1>("sel_mgg_true_log")
        .with_axis(mass_logbins, "m_{#gamma#gamma} [GeV]")
        .fill(mgg_true / 1000);
    S.T<H1>("sel_mgg_true_log_full")
        .with_axis(mass_logbins_full, "m_{#gamma#gamma} [GeV]")
        .fill(mgg_true / 1000);
    
    S.T<H2>("sel_true_v_reco_mgg")
        (70, 0, 7e6, "m_{gg} (true)")
        (70, 0, 7e6, "m_{gg} (reco)")
        .fill(mgg_true, mgg);
        
    S.T<H2>("sel_true_v_reco_mgg_limited")
        (400, 0, 2e6, "m_{gg} (true)")
        (400, 0, 2e6, "m_{gg} (reco)")
        .fill(mgg_true, mgg);
        
    S.T<H2>("sel_true_v_recores_mgg")
        (70, 0, 7e6, "m_{gg} (true)")
        (200, -100e3, 100e3, "m_{gg} (reco) - m_{gg} (true)")
        .fill(mgg_true, mgg - mgg_true);
        
    S.T<H2>("sel_true_v_recoresrel_mgg")
        (70, 0, 7e6, "m_{gg} (true)")
        (400, -0.1, 0.1, "(m_{gg} (reco) - m_{gg} (true)) / m_{gg} (true)")
        .fill(mgg_true, (mgg - mgg_true) / mgg_true);
}


}


int main(int argc, const char* argv[]) {
    return a4::process::a4_main_configuration<ana::Configuration>(argc, argv);
}

