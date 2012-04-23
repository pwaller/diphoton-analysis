#ifndef A4PW_EXTERNAL_H_
#define A4PW_EXTERNAL_H_
#include "egammaAnalysisUtils/egammaSFclass.h"
#include "egammaAnalysisUtils/EnergyRescaler.h"
#include "egammaAnalysisUtils/checkOQ.h"
#include "egammaAnalysisUtils/CaloIsoCorrection.h"
#include "PileupReweighting/TPileupReweighting.h"

#include "FudgeMCTool.h"
#include "PhotonIDTool.h"

using eg2011::EnergyRescaler;

#include <TLorentzVector.h>

//TLorentzVector TLV(ALorentzVector alv);
Root::TPileupReweighting* get_prw(const char* mc_file, const char* data_file);
double ReturnRZ_1stSampling_cscopt2(double eta_1st_sampling);
double GetCorrectedInvMass(double Elead, double etaS1lead, double philead,
    double Esublead, double etaS1sublead, double phisublead, double PV_ID);
double scaleForFFUncovertedPhoton(const double pT);

double ComputeWeight(float truemass, float gravitonPoleMass, float Coupling);

#endif
