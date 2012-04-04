#include <TLorentzVector.h>
#include "external.h"

/*
TLorentzVector TLV(ALorentzVector alv) {
    return TLorentzVector(alv.px, alv.py, alv.pz, alv.E);
}*/

Root::TPileupReweighting* get_prw(const char* mc_file, const char* data_file) {
    auto* pileup_tool = new Root::TPileupReweighting("pileup_reweighting");
    
    pileup_tool->AddConfigFile(mc_file);
    pileup_tool->AddLumiCalcFile(data_file);
    
    //pileup_tool->AddMCDistribution(mc_file, "MCPileupReweighting");
    //pileup_tool->AddDataDistribution(data_file, "LumiMetaData");
    pileup_tool->SetUnrepresentedDataAction(2);
    //pileup_tool->SetDefaultChannel(default_channel);
    pileup_tool->Initialize();
    return pileup_tool;
}

double ReturnRZ_1stSampling_cscopt2(double eta_1st_sampling)
{
    //adapted from CaloDepthTool.cxx, double CaloDepthTool::cscopt2_parametrized(const CaloCell_ID::CaloSample sample,
    //  const double eta, const double /*phi*/ )
    // No warranty !!!

    double millimeter=1;
    double radius = -99999;
    float aeta_1st_sampling = fabs(eta_1st_sampling);
    if (aeta_1st_sampling<1.5) { //barrel
        if (aeta_1st_sampling < 0.8)
            radius = (1558.859292 - 4.990838*aeta_1st_sampling - 21.144279*aeta_1st_sampling*aeta_1st_sampling)*millimeter;
        else
            radius = (1522.775373 + 27.970192*aeta_1st_sampling - 21.104108*aeta_1st_sampling*aeta_1st_sampling)*millimeter;
    }
    else { //endcap
        if (aeta_1st_sampling < 1.5)
            radius = (12453.297448 - 5735.787116*aeta_1st_sampling)*millimeter;
        else
            radius = 3790.671754*millimeter;
        if (eta_1st_sampling < 0.) radius = -radius;
    }
    return radius;
}

double GetCorrectedInvMass(double Elead, double etaS1lead, double philead,
    double Esublead, double etaS1sublead, double phisublead, double PV_ID)
{
    double R_photon_front;
    double Z_photon_front;
    //leading photon
    if (fabs(etaS1lead) < 1.5) {                                    // barrel
        R_photon_front=ReturnRZ_1stSampling_cscopt2(etaS1lead);
        Z_photon_front=R_photon_front*sinh(etaS1lead);
    }
    else {                                                          // endcap
        Z_photon_front=ReturnRZ_1stSampling_cscopt2(etaS1lead);
        R_photon_front=Z_photon_front/sinh(etaS1lead);
    }
    double eta_corrected_leading=asinh((Z_photon_front- PV_ID)/R_photon_front);
    //subleading photon
    if (fabs(etaS1sublead) < 1.5) {      // barrel
        R_photon_front=ReturnRZ_1stSampling_cscopt2(etaS1sublead);
        Z_photon_front=R_photon_front*sinh(etaS1sublead);
    }
    else {                    // endcap
        Z_photon_front=ReturnRZ_1stSampling_cscopt2(etaS1sublead);
        R_photon_front=Z_photon_front/sinh(etaS1sublead);
    }

    double eta_corrected_subleading=asinh((Z_photon_front- PV_ID)/R_photon_front);
    double mass_gg_corrected;
    float energy_corrected_EMscale_leading = Elead;  //apply energy corrections due to EMscale
    float ph_pt_corrected_leading= energy_corrected_EMscale_leading / cosh(eta_corrected_leading);
    TLorentzVector Vp1_corrected;
    Vp1_corrected.SetPtEtaPhiE(ph_pt_corrected_leading,eta_corrected_leading, philead ,energy_corrected_EMscale_leading);
    //subleading photon :
    float energy_corrected_EMscale_subleading= Esublead;
    float ph_pt_corrected_subleading=energy_corrected_EMscale_subleading/cosh(eta_corrected_subleading);
    TLorentzVector Vp2_corrected;
    Vp2_corrected.SetPtEtaPhiE(ph_pt_corrected_subleading,eta_corrected_subleading,phisublead,energy_corrected_EMscale_subleading);
    TLorentzVector Vdiphoton_corrected=Vp1_corrected+Vp2_corrected;
    mass_gg_corrected=Vdiphoton_corrected.M();
    return mass_gg_corrected;
}

double scaleForFFUncovertedPhoton(const double pT)
{
    // temporary obtained by comparing FFMC to electron extrapolation
    // apply this scale factor for FFed unconverted photon
    // in 1.81 < |eta| < 2.37 (use ph_eta2 in DPD)
    // pT in GeV
    if (pT<25.) return 0.86;
    else if (pT<30.) return 0.89;
    else if (pT<35.) return 0.94;
    else if (pT<40.) return 0.92;
    else if (pT<45.) return 0.98;
    else return 0.97;
}
