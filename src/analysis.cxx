#include "all.h"

#include <iostream>

#include <a4/histogram.h>
using a4::hist::H1;
using a4::hist::H2;
using a4::hist::Cutflow;

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


void Analysis::process(const ntup::Event& event) {

    if (event.run_number() != _current_run || 
        event.mc_channel_number() != _current_sample)
        new_sample(event);
    
    #define PASSED(x) \
        S.T<Cutflow>("cutflow").passed(x);
        
    #define EFFPLOT(name) \
        S.T<H1>("eff/true_mgg/" name)(7000, 0, 7e6).fill(true_mgg); \
        
    #define EFFPLOT_1(name) \
    do { \
        EFFPLOT(name); \
        if (is_mc && phtr_1 && phtr_2) { \
            S.T<H1>("eff/true_lead/pt/" name)(3000, 0, 3e6).fill(phtr_1->pt()); \
            S.T<H1>("eff/true_lead/eta/" name)(100, -2.5, 2.5).fill(phtr_1->eta()); \
            S.T<H1>("eff/true_sublead/pt/" name)(3000, 0, 3e6).fill(phtr_2->pt()); \
            S.T<H1>("eff/true_sublead/eta/" name)(100, -2.5, 2.5).fill(phtr_2->eta()); \
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
    
    // MC: Pileup reweighting and run number reassignment
    float pileup_weight = 1.0;
    if (!data && C._do_pileup_reweighting && systematic("prw")) {
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
    
    PASSED("Total");
    
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
    double true_mgg = -999;
    if (gravitons.size())
        true_mgg = gravitons[0]->m();
    else if (is_mc && hard_process_photons.size() >= 2) {
        true_mgg = (phtr_1_lv + phtr_2_lv).m();
    }
    
    if (is_mc && _should_masscut)
        if (true_mgg < _mass_low || true_mgg >= _mass_high)
            return;
    
    if (is_mc && C._require_mc_match && hard_process_photons.size() < 2)
        return;
            
    PASSED("Hardproc 2#gamma");
    
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
    
    S.T<H1>("sel_reco_mgg")(7000, 0, 7e6).fill(mgg);
    
    if (is_mc)
        resolution_plots(mgg, true_mgg);
        
    if (event.larerror() > 1) return;
    PASSED("LarOK");
    
    if (C._ee_events)
        if (C._ee_events->present(event.run_number(), event.event_number()))
            return;
    PASSED("!ee");
    
    if (mgg < 140e3) return;
    PASSED("mgg140");
    
    EFFPLOT_1("11_mass");
}

void Analysis::new_sample(const ntup::Event& event) {
    _current_run = event.run_number();
    _current_sample = event.mc_channel_number();
    _simulation = event.issimulation();
    
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
            _should_masscut = true; _mass_low = 0;    _mass_high = 800e3;
            break;
        case 145606: // mc11_7TeV.145606.Pythiagamgam15_M_gt_800
            _should_masscut = true; _mass_low = 800e3;  _mass_high = 1500e3;
            break;
        case 145607: // mc11_7TeV.145607.Pythiagamgam15_M_gt_1500
            _should_masscut = true; _mass_low = 1500e3; _mass_high = 9999999999999;
            break;
        default:
            _should_masscut = false;
    }
}

bool Analysis::pass_grl(const ntup::Event& event) {
    if (event.issimulation()) return true;
    return C._grl->pass(event.run_number(), event.lbn());
}
    
double Analysis::compute_mass(const ntup::Event& event, const Photon& lead, const Photon& sublead) const {
    double E1 = lead.corrected().e();
    double eta1 = lead->etas1();
    double phi1 = lead->phi();
    
    double E2 = sublead.corrected().e();
    double eta2 = sublead->etas1();
    double phi2 = sublead->phi();
    
    double PV_ID = event.primary_vertices(0).z();
    double mgg = GetCorrectedInvMass(E1, eta1, phi1, E2, eta2, phi2, PV_ID);
    return mgg;
}

double Analysis::resolution_plots(double mgg, double true_mgg) {
    S.T<H1>("sel_true_mgg")(7000, 0, 7e6).fill(true_mgg);
    
    S.T<H2>("sel_true_v_reco_mgg")
        (70, 0, 7e6, "m_{gg} (true)")
        (70, 0, 7e6, "m_{gg} (reco)")
        .fill(true_mgg, mgg);
        
    S.T<H2>("sel_true_v_reco_mgg_limited")
        (400, 0, 2e6, "m_{gg} (true)")
        (400, 0, 2e6, "m_{gg} (reco)")
        .fill(true_mgg, mgg);
        
    S.T<H2>("sel_true_v_recores_mgg")
        (70, 0, 7e6, "m_{gg} (true)")
        (200, -100e3, 100e3, "m_{gg} (reco) - m_{gg} (true)")
        .fill(true_mgg, mgg - true_mgg);
        
    S.T<H2>("sel_true_v_recoresrel_mgg")
        (70, 0, 7e6, "m_{gg} (true)")
        (400, -0.1, 0.1, "(m_{gg} (reco) - m_{gg} (true)) / m_{gg} (true)")
        .fill(true_mgg, (mgg - true_mgg) / true_mgg);
}


}


int main(int argc, const char* argv[]) {
    return a4::process::a4_main_configuration<ana::Configuration>(argc, argv);
}

