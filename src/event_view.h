#ifndef EVENT_VIEW_H
#define EVENT_VIEW_H

#include <vector>

#include <a4/types.h>
#include <a4/atlas/ntup/photon/Event.pb.h>

#include <a4/alorentzvector.h>

namespace ntup = a4::atlas::ntup::photon;

template <class TransientClass, class PersistentClass>
class PersistentWrapper {
protected:
    PersistentClass const* _object;
    
public:    
    PersistentWrapper() : _object(NULL) {}
    PersistentWrapper(PersistentClass const& object) : _object(&object) {}
    PersistentWrapper(PersistentClass const* const object) : _object(object) {}
    
    PersistentClass const* operator->() const { assert(_object); return _object; }
    const PersistentClass& operator* () const { assert(_object); return *_object; }
    
    operator bool() const { return _object != NULL; }
    operator PersistentClass const*() const { assert(_object); return _object; }
    //TransientClass& operator=(PersistentClass const* object) { _object = object; return *this; }
    
    template<class Container>
    static std::vector<TransientClass> make_vector(Container const& container) { 
        std::vector<TransientClass> result;
        foreach (auto const& object, container)
            result.push_back(TransientClass(object));
        return result;
    }
    
};

//class Event;

class Photon : public PersistentWrapper<Photon, ntup::Photon> {
private:
    shared<ntup::Photon> _corrected;
    mutable shared<ALorentzVector> _lv;
    
public:
    Photon() { init(); }
    explicit Photon(ntup::Photon const& ph) : PersistentWrapper<Photon, ntup::Photon>(ph) { init(); }
    explicit Photon(ntup::Photon const* ph) : PersistentWrapper<Photon, ntup::Photon>(ph) { init(); }
    
    void init() {
    }
    
    /// Corrections applied
    /// 1. Energy scale correction (data), Smearing correction (MC)
    /// 2. Shower shape fudge factors
    /// 3. PhotonIDTool: isem, loose, tight
    void compute_corrections(const ntup::Event& event, 
        const int original_index, EnergyRescaler& rescaler) {
        
        auto& orig_ph = *_object;
        compute_extra_quantities(*const_cast<ntup::Photon*>(_object));
        
        _corrected.reset(new ntup::Photon());
        auto& ph = *_corrected;
        
        #define COPY(what) ph.set_##what(orig_ph.what())
        COPY(pt);
        COPY(etas2);
        COPY(cl_e);
        COPY(cl_eta);
        COPY(cl_phi);
        COPY(phi);
        COPY(isconv);
        COPY(etcone40);
        COPY(etcone40_ed_corrected);
        COPY(etap);
        
        COPY(ethad);
        COPY(ethad1);
        COPY(rhad);
        COPY(rhad1);
        COPY(e277);
        COPY(reta);
        COPY(rphi);
        COPY(weta2);
        COPY(f1);
        COPY(fside);
        COPY(wstot);
        COPY(ws3);
        COPY(deltae);
        COPY(eratio);
        
        COPY(loose);
        COPY(tight);
        COPY(isem);
        #undef COPY
        
        auto is_mc = event.issimulation();
        
        double factor = 1;
        double new_e = -999;
        
        if (is_mc) {
            const bool not_mc11c = false;
            rescaler.SetRandomSeed(1771561 + event.event_number() + (original_index * 10));
            factor = rescaler.getSmearingCorrectionMeV(
                ph.cl_eta(), ph.cl_e(), 0, not_mc11c, "2011");
            
            //DEBUG("  Smearing factor: ", factor);
            
            new_e = ph.cl_e() * factor;
            //DEBUG("  new energy: ", new_e);
        } else {
            new_e = rescaler.applyEnergyCorrectionMeV(
                ph.cl_eta(), ph.cl_phi(), ph.cl_e(), 
                ph.cl_e() / cosh(ph.cl_eta()),
                0, "PHOTON");
            factor = new_e / ph.cl_e();
        }
                
        // Update corrected photon after fudging
        ph.set_original_index(original_index);
        
        const double new_et = new_e / cosh(ph.etas2());
        const double original_pt = ph.pt();
        
        ph.set_rhad(ph.ethad() / new_et);
        ph.set_rhad1(ph.ethad1() / new_et);
                
        ph.set_e(new_e);
        ph.set_pt(new_et);

        if (is_mc) {
            
            // Shower fudging
        
            double rhad = ph.rhad(),
                    rhad1 = ph.rhad1(),
                    e277 = ph.e277(),
                    reta = ph.reta(),
                    rphi = ph.rphi(),
                    weta2 = ph.weta2(),
                    f1 = ph.f1(),
                    fside = ph.fside(),
                    wstot = ph.wstot(),
                    ws3 = ph.ws3(),
                    deltae = ph.deltae(),
                    eratio = ph.eratio();
                    
            FudgeMCTool fmc(new_et, ph.etas2(), ph.isconv(), 0);
            fmc.SetPreselection(8);
            
            fmc.FudgeShowers(
                new_et, ph.etas2(), rhad1, rhad, e277, reta, rphi, weta2,
                f1, fside, wstot, ws3, deltae, eratio, ph.isconv());
                
            ph.set_rhad(rhad);
            ph.set_rhad1(rhad1);
            ph.set_reta(reta);
            ph.set_rphi(rphi);
            ph.set_weta2(weta2);
            ph.set_f1(f1);
            ph.set_fside(fside);
            ph.set_wstot(wstot);
            ph.set_ws3(ws3);
            ph.set_deltae(deltae);
            ph.set_eratio(eratio);
            
            PhotonIDTool selection = PhotonIDTool(
                //ph.cl_e() * factor / cosh(ph.etas2()),
                new_et,
                ph.etas2(),
                ph.rhad1(),
                ph.rhad(), 
                e277,
                ph.reta(),
                ph.rphi(),
                weta2,
                f1,
                fside,
                wstot,
                ws3,
                deltae,
                eratio,
                ph.isconv());
                
            ph.set_isem(selection.isEM(3, 6));
            ph.set_loose(selection.PhotonCutsLoose(3));
            ph.set_tight(selection.PhotonCutsTight(6));
        }
        
        // Isolation
                            
                //new_e,
        double isolation = CaloIsoCorrection::GetPtEDCorrectedIsolation(
            ph.etcone40(),
            ph.etcone40_ed_corrected(),
            //orig_ph.cl_e(),
            new_e,
            ph.etas2(),
            ph.etap(),
            ph.cl_eta(),
            40,
            is_mc, 
            ph.etcone40(), 
            ph.isconv(),
            CaloIsoCorrection::PHOTON);
            
        ph.set_analysis_isolation(isolation);
        
    }
        
    bool isolated() {
        return corrected().analysis_isolation() < 5000;
    }
    
    double scale_factor() {
        auto& ph = **this;
        if (not ph.isconv() and 1.81 < abs(ph.etas2()) && abs(ph.etas2()) < 2.37 )
            return scaleForFFUncovertedPhoton(ph.cl_pt() / 1e3);
        return 1;
    }
    
    static void compute_extra_quantities(ntup::Photon& ph) {
        if (!ph.has_rhad()) ph.set_rhad(ph.ethad()  / (ph.cl_e() / cosh(ph.etas2())));
        if (!ph.has_rhad1()) ph.set_rhad1(ph.ethad1() / (ph.cl_e() / cosh(ph.etas2())));
        if (!ph.has_deltae()) ph.set_deltae(ph.emax2() - ph.emins1());
        
        double eratio = 0;
        if (fabs(ph.emaxs1() + ph.emax2()) > 0) 
            eratio = (ph.emaxs1() - ph.emax2()) / (ph.emaxs1() + ph.emax2()); 
        if (!ph.has_eratio()) ph.set_eratio(eratio);
    }
    
    const ntup::Photon& corrected() const {
        if (!_corrected)
            FATAL("Corrections used before compute_corrections(ntup::Event) called");
        return *_corrected;
    }
    
    const ALorentzVector& lv() const {
        if (!_lv)
            _lv.reset(new ALorentzVector(ALorentzVector::from_ptetaphie(**this)));
        return *_lv;
    }
};

/*
class Event : public PersistentWrapper<Event, ntup::Event> {
public:
    Event() {}
    Event(ntup::Event const& o) : PersistentWrapper<Event, ntup::Event>(o) {}
    Event(ntup::Event const* o) : PersistentWrapper<Event, ntup::Event>(o) {}
    
    mutable shared<std::vector<Photon>> _photons;
    
    uint32_t ntrack_vertices(int ntrack=2) {
        uint32_t result = 0;
        foreach (auto& v, (*this)->primary_vertices())
            if (v.ntracks() >= ntrack)
                result++;
        return result
    }
    
    const std::vector<Photon>& photons() const {
        if (!_photons) {
            _photons.reset(new std::vector<Photon>(Photon::make_vector((*this)->photons())));
            foreach (auto ph, *_photons)
                ph._event = this;
        }
        return *_photons; 
    }
};
*/

#endif
