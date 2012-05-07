#include "all.h"

#include <unordered_map>

#include <a4/utility.h>
// SORT_KEY, REMOVE_IF
using a4::process::utility::vector_of_ptr;
using a4::process::utility::in_map;

#include "analysis.h"

#include <a4/atlas/ntup/photon/Event.pb.h>
namespace ntup = a4::atlas::ntup::photon;

#define KEEP(whats, what, ...) do { \
    auto new_##what = new_event.add_##whats(); \
    new_##what->CopyFrom(what); \
    __VA_ARGS__ \
    } while (0)
#define FILTER(idx, def, whats) \
    new_event.clear_##whats(); \
    foreach_enumerate (idx, def, event.whats())

typedef std::unordered_map<unsigned int, ntup::Photon*> IndexOfPhoton;
typedef std::unordered_map<const ntup::PhotonTruthParticle*, 
                           const ntup::Photon*> TruthToRecoPhotonMap;

namespace ana {


void Filter::process(const ntup::Event& event) {
    ntup::Event new_event(event);
    filter_photons(event, new_event);
    write(new_event);
}
    
void Filter::filter_photons(const ntup::Event& event, ntup::Event& new_event) {
    IndexOfPhoton photon_truth_to_keep, photon_ef_to_keep;
    
    if (C._filter_reco_photons) {
        FILTER(original_index, const ntup::Photon& ph, photons) {
         
            // Reject photons in the crack
            //if (etas2_crack(ph)) continue;
                        
            // Reject low pt photons
            if (ph.pt() < 20e3)
                continue;
            
            KEEP(photons, ph,
                new_ph->set_original_index(original_index);
                // Record true photon to keep
                if (ph.has_truth_matched() && ph.truth_matched())
                    if (ph.truth_index() != -1)
                        photon_truth_to_keep[ph.truth_index()] = new_ph;
                
                // Record EFPhoton record to keep
                if (ph.has_ef_index() && ph.ef_index() != -1)
                    if (ph.ef_index() != -1)
                        photon_ef_to_keep[ph.ef_index()] = new_ph;
            );
        }
    } else {
        foreach (auto& new_ph, *new_event.mutable_photons()) {
            if (new_ph.truth_index() != -1)
                photon_truth_to_keep[new_ph.truth_index()] = &new_ph;
            if (new_ph.ef_index() != -1)
                photon_ef_to_keep[new_ph.ef_index()] = &new_ph;
        }
    }
    
    auto interesting_true_ph = [](const ntup::PhotonTruthParticle& ph) {
        return ph.pdgid() == 5000039 || 
               ph.ishardprocphoton() ||
               ph.pt() > 20e3;
    };
    
    // Get the parent barcodes of interesting particles
    std::unordered_set<uint32_t> interesting_parents;
    foreach (auto& ph, event.photon_truth_particles())
        if (interesting_true_ph(ph))
            foreach (auto& parent, ph.parents())
                interesting_parents.insert(parent);
            
    // Keep photons from the hard process or which are matched to 
    // reconstructed photons we kept.
    unsigned int new_index = 0;
    FILTER(old_index, const ntup::PhotonTruthParticle& ph, photon_truth_particles) {
        
        bool has_kept_reconstruction = in_map(photon_truth_to_keep, old_index);
        bool interesting = interesting_true_ph(ph)
                           || interesting_parents.count(ph.barcode());
        
        if (has_kept_reconstruction || interesting) {
        
            KEEP(photon_truth_particles, ph);
            
            // Update .photons()[].truth_index()
            if (has_kept_reconstruction)
                photon_truth_to_keep[old_index]->set_truth_index(new_index);
                
            new_index++;
        }
    }
    
    // Keep relevent EF records
    new_index = 0;
    FILTER(old_index, const ntup::EFPhoton& ph, efphotons) {
        if (in_map(photon_ef_to_keep, old_index)) {
                    
            KEEP(efphotons, ph);
            
            photon_ef_to_keep[old_index]->set_ef_index(new_index);
            new_index++;
        }
    }
}


}

