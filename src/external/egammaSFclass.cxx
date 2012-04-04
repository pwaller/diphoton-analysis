// Efficiency scale factors macro
// Date: 11/01/2010
// Author: Olivier Arnaez <olivier.arnaez@cern.ch>
//         Jan Kretzschmar <jan.kretzschmar@cern.ch>
//
// Usage: 
// std::pair<float,float> sf_error = objsf->scaleFactor(eta(cluster), ET(MeV), set, range, rel, etcorrection)
//
// Please note: alternative accessors scaleFactorLoose, scaleFactorMedium, ... 
//              have been disabled, as the amount of sets is expanding and this was not maintainable
//              Read the "set=X" documentation below and use the correct arguments for scaleFactor() function
//
// The first number (sf_error.first) returns the efficiency scale factor,
// the second number is its uncertainty (sf_errror.second)
// 
// The combined (W/Z/jpsi) scale factor and uncertainty vs eta and ET (in MeV) are given
// 
// different "sets" of numbers available (not all with all sets):
//    * Loose SF (set=0)
//    * Medium SF (set=1)
//    * Tight SF (set=2)
//    * e20_medium trigger SF (set=3) (use set >=8 for release 17 2011 data/MC11a)
//    * reco+trkqual SF (set=4)
//    * Loose++ SF (set=5)
//    * Medium++ SF (set=6)
//    * Tight++ SF (set=7)
//    * Forward Loose SF, detailed G4 FCAL (set=23)
//    * Forward Tight SF, detailed G4 FCAL  (set=24)
//    * Forward Loose SF, Frozen showers FCAL (set=25)
//    * Forward Tight SF, Frozen showers FCAL  (set=26)
//
//    * e20_medium trigger SF w.r.t Medium++ offline (set=8)
//    * e20_medium MC efficiency w.r.t Medium++ offline (set=9)
//    * e20_medium trigger SF w.r.t Tight++ offline (set=10)
//    * e20_medium MC efficiency w.r.t Tight++ offline (set=11)
//    * e22_medium trigger SF w.r.t Medium++ offline (set=12)
//    * e22_medium MC efficiency w.r.t Medium++ offline (set=13)
//    * e22_medium trigger SF w.r.t Tight++ offline (set=14)
//    * e22_medium MC efficiency w.r.t Tight++ offline (set=15)
//    * e22vh_medium1 trigger SF (using e22_medium1 on MC11a) w.r.t Medium++ offline (set=16)
//    * e22_medium1 MC efficiency w.r.t Medium++ offline (set=17)
//    * e22vh_medium1 trigger SF (using e22_medium1 on MC11a) w.r.t Tight++ offline (set=18)
//    * e22_medium1 MC efficiency w.r.t Tight++ offline (set=19)
//    * e20_medium MC efficiency w.r.t Loose++ offline (set=20)
//    * e22_medium1 MC efficiency w.r.t Loose++ offline (set=21)
//    * e22vh_medium1 MC efficiency w.r.t Loose++ offline (set=22)
//
//    * e20_medium SF  w.r.t Loose++ offline (set=27)  
//    * e22_medium1 SF  w.r.t Loose++ offline (set=28)
//    * e22vh_medium1 SF  w.r.t Loose++ offline (set=29)
///   


// data and MC release selection:
//    * release 15 2010 data/MC09 (rel=0)
//    * release 16 2010 data/MC10 (rel=1)
//    * release 16.6 2011 data/MC10ab (estimated from 2010 data) (rel=2)
//    * release 16.6 estimated from 2011 data "EPS recommendations" and MC10b (rel=3)
//    * release 16.6 estimated from 2011 data "EPS recommendations" including Jpsi measurements (rel=4)
//    * release 17 estimated from 2011 data/MC11a "CERN council recommendations" (rel=5)
//    * release 17 estimated from 2011 data/MC11a/b/c "Moriond recommendations" G4 FullSim MC (rel=6)
//    * release 17 estimated from 2011 data/MC11a/b/c "Moriond recommendations" AFII MC (rel=7)
// measured with probes in the ranges:
//    * 20-50 GeV range (range=0)
//    * 30-50 GeV (range=1)
// and correcting (etcorrection=1) or not (etcorrection=0) for the ET-dependence
//
// Eta binning is changing from release to release
//
// Note that for rel>=4 range should be left at 0 and etcorrection=1 always
// 
// For now separete function for Forward Electrons (|eta|>=2.5)
// and ONLY for release 16.6, 2011 data in analogy to function for central electrons
// std::pair<float,float> sf_error = objsf->scaleFactorForward(eta, set)
// where eta = electron eta = cluster eta (as no track)
// cuts are ForwardLoose (set=0) or ForwardTight (set=2)
// 

#include "egammaAnalysisUtils/egammaSFclass.h"
#include <cmath>

egammaSFclass::egammaSFclass()
{
  //Definition of the eta binning
  m_Etabins.push_back(-2.47);
  m_Etabins.push_back(-2.01); 
  m_Etabins.push_back(-1.52); 
  m_Etabins.push_back(-1.37); 
  m_Etabins.push_back(-0.8); 
  m_Etabins.push_back(0); 
  m_Etabins.push_back(0.8); 
  m_Etabins.push_back(1.37); 
  m_Etabins.push_back(1.52); 
  m_Etabins.push_back(2.01); 
  m_Etabins.push_back(2.47);
  //Definition of the fine eta binning
  m_FineEtabins.push_back(-2.47);
  m_FineEtabins.push_back(-2.37);
  m_FineEtabins.push_back(-2.01);
  m_FineEtabins.push_back(-1.81);
  m_FineEtabins.push_back(-1.52);
  m_FineEtabins.push_back(-1.37);
  m_FineEtabins.push_back(-1.15);
  m_FineEtabins.push_back(-0.8 );
  m_FineEtabins.push_back(-0.6 );
  m_FineEtabins.push_back(-0.1 );
  m_FineEtabins.push_back( 0.  );
  m_FineEtabins.push_back( 0.1 );
  m_FineEtabins.push_back( 0.6 );
  m_FineEtabins.push_back( 0.8 );
  m_FineEtabins.push_back( 1.15);
  m_FineEtabins.push_back( 1.37);
  m_FineEtabins.push_back( 1.52);
  m_FineEtabins.push_back( 1.81);
  m_FineEtabins.push_back( 2.01);
  m_FineEtabins.push_back( 2.37);
  m_FineEtabins.push_back( 2.47);
  //Definition of the eta binning with 11 bins
  m_11Etabins.push_back(-2.47);
  m_11Etabins.push_back(-2.01); 
  m_11Etabins.push_back(-1.52); 
  m_11Etabins.push_back(-1.37); 
  m_11Etabins.push_back(-0.8); 
  m_11Etabins.push_back(-0.1); 
  m_11Etabins.push_back(0.1); 
  m_11Etabins.push_back(0.8); 
  m_11Etabins.push_back(1.37); 
  m_11Etabins.push_back(1.52); 
  m_11Etabins.push_back(2.01); 
  m_11Etabins.push_back(2.47);
  //Definition of the eta binning for forward electrons
  m_FwdEtabins.push_back(2.5);
  m_FwdEtabins.push_back(2.6);
  m_FwdEtabins.push_back(2.7);
  m_FwdEtabins.push_back(2.8);
  m_FwdEtabins.push_back(2.9);
  m_FwdEtabins.push_back(3.0);
  m_FwdEtabins.push_back(3.16);
  m_FwdEtabins.push_back(3.35);
  m_FwdEtabins.push_back(3.6);
  m_FwdEtabins.push_back(4.0);
  m_FwdEtabins.push_back(4.9);

  //Definition of the ET binning
  m_ETbins.push_back(0.);
  m_ETbins.push_back(20000.); 
  m_ETbins.push_back(25000.); 
  m_ETbins.push_back(30000.); 
  m_ETbins.push_back(35000.); 
  m_ETbins.push_back(40000.); 
  m_ETbins.push_back(45000.); 
  m_ETbins.push_back(500000000.); 
  //Definition of the ET binning on the full range
  m_ETbinsFullRange.push_back(    0.);
  m_ETbinsFullRange.push_back( 7000.);
  m_ETbinsFullRange.push_back(10000.);
  m_ETbinsFullRange.push_back(15000.);
  m_ETbinsFullRange.push_back(20000.); 
  m_ETbinsFullRange.push_back(25000.); 
  m_ETbinsFullRange.push_back(30000.); 
  m_ETbinsFullRange.push_back(35000.); 
  m_ETbinsFullRange.push_back(40000.); 
  m_ETbinsFullRange.push_back(45000.); 
  m_ETbinsFullRange.push_back(500000000.); 
  //Definition of the ET binning for trigger
  m_ETbinsTrigger.push_back(21000.);
  m_ETbinsTrigger.push_back(23000.);
  m_ETbinsTrigger.push_back(25000.); 
  m_ETbinsTrigger.push_back(30000.); 
  m_ETbinsTrigger.push_back(35000.); 
  m_ETbinsTrigger.push_back(40000.); 
  m_ETbinsTrigger.push_back(500000000.); 


  //For the scale factors of the standard egamma cuts 

  //Release 15
  //Probes between 30 and 50 GeV (plateau region)
  //Loose
  efficienciesRel15Loose3050.push_back(98.1); 
  efficienciesRel15Loose3050.push_back(99.0); 
  efficienciesRel15Loose3050.push_back(0.); 
  efficienciesRel15Loose3050.push_back(98.6); 
  efficienciesRel15Loose3050.push_back(99.5); 
  efficienciesRel15Loose3050.push_back(99.1); 
  efficienciesRel15Loose3050.push_back(98.8); 
  efficienciesRel15Loose3050.push_back(0.); 
  efficienciesRel15Loose3050.push_back(99.9); 
  efficienciesRel15Loose3050.push_back(98.2);
  uncertaintiesRel15Loose3050.push_back(1.6); 
  uncertaintiesRel15Loose3050.push_back( 1.5); 
  uncertaintiesRel15Loose3050.push_back(0.); 
  uncertaintiesRel15Loose3050.push_back( 1.5); 
  uncertaintiesRel15Loose3050.push_back( 1.5); 
  uncertaintiesRel15Loose3050.push_back( 1.5); 
  uncertaintiesRel15Loose3050.push_back( 1.5); 
  uncertaintiesRel15Loose3050.push_back(0.); 
  uncertaintiesRel15Loose3050.push_back( 1.5); 
  uncertaintiesRel15Loose3050.push_back( 1.6);
  //Medium
  efficienciesRel15Medium3050.push_back(95.4); 
  efficienciesRel15Medium3050.push_back(98.7);
  efficienciesRel15Medium3050.push_back(0.); 
  efficienciesRel15Medium3050.push_back(97.9);
  efficienciesRel15Medium3050.push_back(98.1);
  efficienciesRel15Medium3050.push_back(97.7); 
  efficienciesRel15Medium3050.push_back(97.9); 
  efficienciesRel15Medium3050.push_back(0.); 
  efficienciesRel15Medium3050.push_back(99.9); 
  efficienciesRel15Medium3050.push_back(97.4);
  uncertaintiesRel15Medium3050.push_back(1.7);
  uncertaintiesRel15Medium3050.push_back( 1.6);
  uncertaintiesRel15Medium3050.push_back(0.); 
  uncertaintiesRel15Medium3050.push_back( 1.6);
  uncertaintiesRel15Medium3050.push_back( 1.5); 
  uncertaintiesRel15Medium3050.push_back( 1.5); 
  uncertaintiesRel15Medium3050.push_back( 1.5); 
  uncertaintiesRel15Medium3050.push_back(0.); 
  uncertaintiesRel15Medium3050.push_back( 1.6);
  uncertaintiesRel15Medium3050.push_back( 1.7);
  //Tight
  efficienciesRel15Tight3050.push_back(92.3); 
  efficienciesRel15Tight3050.push_back(99.2); 
  efficienciesRel15Tight3050.push_back(0.);
  efficienciesRel15Tight3050.push_back(101.5); 
  efficienciesRel15Tight3050.push_back(98.9); 
  efficienciesRel15Tight3050.push_back(99.9);
  efficienciesRel15Tight3050.push_back(104.2); 
  efficienciesRel15Tight3050.push_back(0.);
  efficienciesRel15Tight3050.push_back(102.6); 
  efficienciesRel15Tight3050.push_back(95.5);
  uncertaintiesRel15Tight3050.push_back(3.3);
  uncertaintiesRel15Tight3050.push_back( 2.3); 
  uncertaintiesRel15Tight3050.push_back(0.);
  uncertaintiesRel15Tight3050.push_back( 2.0); 
  uncertaintiesRel15Tight3050.push_back( 1.8); 
  uncertaintiesRel15Tight3050.push_back( 1.8);
  uncertaintiesRel15Tight3050.push_back( 2.5); 
  uncertaintiesRel15Tight3050.push_back(0.); 
  uncertaintiesRel15Tight3050.push_back( 5.0); 
  uncertaintiesRel15Tight3050.push_back( 3.2);

  //Probes between 20 and 50 GeV
  //Loose
  efficienciesRel15Loose2050.push_back(97.6); 
  efficienciesRel15Loose2050.push_back(99.0); 
  efficienciesRel15Loose2050.push_back(0.); 
  efficienciesRel15Loose2050.push_back(98.2); 
  efficienciesRel15Loose2050.push_back(99.1); 
  efficienciesRel15Loose2050.push_back(98.8); 
  efficienciesRel15Loose2050.push_back(98.2); 
  efficienciesRel15Loose2050.push_back(0.); 
  efficienciesRel15Loose2050.push_back(99.6); 
  efficienciesRel15Loose2050.push_back(97.4);
  uncertaintiesRel15Loose2050.push_back(1.6); 
  uncertaintiesRel15Loose2050.push_back(1.5); 
  uncertaintiesRel15Loose2050.push_back(0.); 
  uncertaintiesRel15Loose2050.push_back( 1.5); 
  uncertaintiesRel15Loose2050.push_back( 1.5); 
  uncertaintiesRel15Loose2050.push_back( 1.5); 
  uncertaintiesRel15Loose2050.push_back( 1.5); 
  uncertaintiesRel15Loose2050.push_back(0.); 
  uncertaintiesRel15Loose2050.push_back( 1.5); 
  uncertaintiesRel15Loose2050.push_back( 1.6);
  //Medium
  efficienciesRel15Medium2050.push_back(94.5); 
  efficienciesRel15Medium2050.push_back(98.8);
  efficienciesRel15Medium2050.push_back(0.); 
  efficienciesRel15Medium2050.push_back(97.2);
  efficienciesRel15Medium2050.push_back(97.4);
  efficienciesRel15Medium2050.push_back(97.2); 
  efficienciesRel15Medium2050.push_back(96.7); 
  efficienciesRel15Medium2050.push_back(0.); 
  efficienciesRel15Medium2050.push_back(99.5); 
  efficienciesRel15Medium2050.push_back(96.1);
  uncertaintiesRel15Medium2050.push_back(1.7);
  uncertaintiesRel15Medium2050.push_back( 1.6);
  uncertaintiesRel15Medium2050.push_back(0.); 
  uncertaintiesRel15Medium2050.push_back( 1.6);
  uncertaintiesRel15Medium2050.push_back( 1.5); 
  uncertaintiesRel15Medium2050.push_back( 1.5); 
  uncertaintiesRel15Medium2050.push_back( 1.5); 
  uncertaintiesRel15Medium2050.push_back(0.); 
  uncertaintiesRel15Medium2050.push_back( 2.9);
  uncertaintiesRel15Medium2050.push_back( 1.7);
  //Tight
  efficienciesRel15Tight2050.push_back(92.5); 
  efficienciesRel15Tight2050.push_back(99.5); 
  efficienciesRel15Tight2050.push_back(0.);
  efficienciesRel15Tight2050.push_back(100.6); 
  efficienciesRel15Tight2050.push_back(98.2); 
  efficienciesRel15Tight2050.push_back(98.7);
  efficienciesRel15Tight2050.push_back(103.3); 
  efficienciesRel15Tight2050.push_back(0.);
  efficienciesRel15Tight2050.push_back(102.8); 
  efficienciesRel15Tight2050.push_back(93.6);
  uncertaintiesRel15Tight2050.push_back(3.4);
  uncertaintiesRel15Tight2050.push_back( 2.4); 
  uncertaintiesRel15Tight2050.push_back(0.);
  uncertaintiesRel15Tight2050.push_back( 2.1); 
  uncertaintiesRel15Tight2050.push_back( 1.8); 
  uncertaintiesRel15Tight2050.push_back( 1.8);
  uncertaintiesRel15Tight2050.push_back( 2.5); 
  uncertaintiesRel15Tight2050.push_back(0.); 
  uncertaintiesRel15Tight2050.push_back( 4.5); 
  uncertaintiesRel15Tight2050.push_back( 3.4);


  //Release 16
  //Probes between 30 and 50 GeV (plateau region)
  //Medium
  efficienciesRel16Medium3050.push_back(98.8); 
  efficienciesRel16Medium3050.push_back(98.0);
  efficienciesRel16Medium3050.push_back(96.9); 
  efficienciesRel16Medium3050.push_back(98.0);
  efficienciesRel16Medium3050.push_back(97.4);
  efficienciesRel16Medium3050.push_back(98.1); 
  efficienciesRel16Medium3050.push_back(98.1); 
  efficienciesRel16Medium3050.push_back(98.3); 
  efficienciesRel16Medium3050.push_back(98.6); 
  efficienciesRel16Medium3050.push_back(97.5);
  uncertaintiesRel16Medium3050.push_back(0.8);
  uncertaintiesRel16Medium3050.push_back(0.9);
  uncertaintiesRel16Medium3050.push_back(2.5); 
  uncertaintiesRel16Medium3050.push_back(0.8);
  uncertaintiesRel16Medium3050.push_back(0.7); 
  uncertaintiesRel16Medium3050.push_back(0.7); 
  uncertaintiesRel16Medium3050.push_back(0.8); 
  uncertaintiesRel16Medium3050.push_back(2.6); 
  uncertaintiesRel16Medium3050.push_back(0.8);
  uncertaintiesRel16Medium3050.push_back(0.8);
  //Tight
  efficienciesRel16Tight3050.push_back(102.0); 
  efficienciesRel16Tight3050.push_back(102.7); 
  efficienciesRel16Tight3050.push_back(114.4);
  efficienciesRel16Tight3050.push_back(106.7); 
  efficienciesRel16Tight3050.push_back( 99.0); 
  efficienciesRel16Tight3050.push_back(100.1);
  efficienciesRel16Tight3050.push_back(105.7); 
  efficienciesRel16Tight3050.push_back(110.8);
  efficienciesRel16Tight3050.push_back(104.2); 
  efficienciesRel16Tight3050.push_back(102.7);
  uncertaintiesRel16Tight3050.push_back(3.0);
  uncertaintiesRel16Tight3050.push_back(1.1); 
  uncertaintiesRel16Tight3050.push_back(3.9);
  uncertaintiesRel16Tight3050.push_back(1.1); 
  uncertaintiesRel16Tight3050.push_back(0.8); 
  uncertaintiesRel16Tight3050.push_back(0.8);
  uncertaintiesRel16Tight3050.push_back(0.9); 
  uncertaintiesRel16Tight3050.push_back(4.6); 
  uncertaintiesRel16Tight3050.push_back(2.6); 
  uncertaintiesRel16Tight3050.push_back(1.2);

  //Probes between 20 and 50 GeV
  //Medium
  efficienciesRel16Medium2050.push_back(97.6); 
  efficienciesRel16Medium2050.push_back(96.8);
  efficienciesRel16Medium2050.push_back(97.7); 
  efficienciesRel16Medium2050.push_back(97.1);
  efficienciesRel16Medium2050.push_back(96.8);
  efficienciesRel16Medium2050.push_back(97.6); 
  efficienciesRel16Medium2050.push_back(97.2); 
  efficienciesRel16Medium2050.push_back(98.2); 
  efficienciesRel16Medium2050.push_back(97.9); 
  efficienciesRel16Medium2050.push_back(96.2);
  uncertaintiesRel16Medium2050.push_back(1.0);
  uncertaintiesRel16Medium2050.push_back(1.0);
  uncertaintiesRel16Medium2050.push_back(3.3); 
  uncertaintiesRel16Medium2050.push_back(1.1);
  uncertaintiesRel16Medium2050.push_back(0.8); 
  uncertaintiesRel16Medium2050.push_back(0.8); 
  uncertaintiesRel16Medium2050.push_back(0.9); 
  uncertaintiesRel16Medium2050.push_back(3.2); 
  uncertaintiesRel16Medium2050.push_back(1.0);
  uncertaintiesRel16Medium2050.push_back(2.8);
  //Tight
  efficienciesRel16Tight2050.push_back(100.2); 
  efficienciesRel16Tight2050.push_back(101.5); 
  efficienciesRel16Tight2050.push_back(117.9);
  efficienciesRel16Tight2050.push_back(105.7); 
  efficienciesRel16Tight2050.push_back( 98.1); 
  efficienciesRel16Tight2050.push_back( 99.1);
  efficienciesRel16Tight2050.push_back(105.2); 
  efficienciesRel16Tight2050.push_back(113.9);
  efficienciesRel16Tight2050.push_back(103.8); 
  efficienciesRel16Tight2050.push_back(101.2);
  uncertaintiesRel16Tight2050.push_back(1.1);
  uncertaintiesRel16Tight2050.push_back(1.2); 
  uncertaintiesRel16Tight2050.push_back(4.4);
  uncertaintiesRel16Tight2050.push_back(1.5); 
  uncertaintiesRel16Tight2050.push_back(0.9); 
  uncertaintiesRel16Tight2050.push_back(1.0);
  uncertaintiesRel16Tight2050.push_back(1.1); 
  uncertaintiesRel16Tight2050.push_back(5.2); 
  uncertaintiesRel16Tight2050.push_back(3.0); 
  uncertaintiesRel16Tight2050.push_back(1.3);


  //For the ET-corrections of the scale factors
  //Medium
  ETCorrectionsMediumRel16.push_back( 79.6);
  ETCorrectionsMediumRel16.push_back( 93.9);
  ETCorrectionsMediumRel16.push_back( 96.2);
  ETCorrectionsMediumRel16.push_back( 99.7);
  ETCorrectionsMediumRel16.push_back(100.6);
  ETCorrectionsMediumRel16.push_back(100.4);
  ETCorrectionsMediumRel16.push_back(101.00);
  uncertaintiesETCorrectionsMediumRel16.push_back( 9.4);
  uncertaintiesETCorrectionsMediumRel16.push_back( 3.6);
  uncertaintiesETCorrectionsMediumRel16.push_back( 1.4);
  uncertaintiesETCorrectionsMediumRel16.push_back( 0.7);
  uncertaintiesETCorrectionsMediumRel16.push_back( 0.5);
  uncertaintiesETCorrectionsMediumRel16.push_back( 0.7);
  uncertaintiesETCorrectionsMediumRel16.push_back( 1.7);
  //Medium
  ETCorrectionsTightRel16.push_back( 76.7);
  ETCorrectionsTightRel16.push_back( 93.6);
  ETCorrectionsTightRel16.push_back( 95.1);
  ETCorrectionsTightRel16.push_back( 99.9);
  ETCorrectionsTightRel16.push_back(100.4);
  ETCorrectionsTightRel16.push_back(100.0);
  ETCorrectionsTightRel16.push_back(100.7);
  uncertaintiesETCorrectionsTightRel16.push_back(10.0);
  uncertaintiesETCorrectionsTightRel16.push_back( 3.7);
  uncertaintiesETCorrectionsTightRel16.push_back( 1.6);
  uncertaintiesETCorrectionsTightRel16.push_back( 0.9);
  uncertaintiesETCorrectionsTightRel16.push_back( 0.7);
  uncertaintiesETCorrectionsTightRel16.push_back( 0.9);
  uncertaintiesETCorrectionsTightRel16.push_back( 1.8);



  //Release 16.6 Data 2010
  //Probes between 30 and 50 GeV (plateau region)
  //Medium
  efficienciesRel166Data2010Medium3050.push_back(98.44); 
  efficienciesRel166Data2010Medium3050.push_back(96.93);
  efficienciesRel166Data2010Medium3050.push_back(96.61); 
  efficienciesRel166Data2010Medium3050.push_back(96.87);
  efficienciesRel166Data2010Medium3050.push_back(97.06);
  efficienciesRel166Data2010Medium3050.push_back(97.49); 
  efficienciesRel166Data2010Medium3050.push_back(97.04); 
  efficienciesRel166Data2010Medium3050.push_back(97.17); 
  efficienciesRel166Data2010Medium3050.push_back(97.31); 
  efficienciesRel166Data2010Medium3050.push_back(97.51);
  uncertaintiesRel166Data2010Medium3050.push_back(2.14);
  uncertaintiesRel166Data2010Medium3050.push_back(2.20);
  uncertaintiesRel166Data2010Medium3050.push_back(2.84); 
  uncertaintiesRel166Data2010Medium3050.push_back(2.13);
  uncertaintiesRel166Data2010Medium3050.push_back(2.18); 
  uncertaintiesRel166Data2010Medium3050.push_back(2.10); 
  uncertaintiesRel166Data2010Medium3050.push_back(2.13); 
  uncertaintiesRel166Data2010Medium3050.push_back(2.89); 
  uncertaintiesRel166Data2010Medium3050.push_back(2.13);
  uncertaintiesRel166Data2010Medium3050.push_back(2.21);
  //Tight
  efficienciesRel166Data2010Tight3050.push_back(101.47); 
  efficienciesRel166Data2010Tight3050.push_back(104.02); 
  efficienciesRel166Data2010Tight3050.push_back(112.70);
  efficienciesRel166Data2010Tight3050.push_back(106.82); 
  efficienciesRel166Data2010Tight3050.push_back( 99.35); 
  efficienciesRel166Data2010Tight3050.push_back(100.13);
  efficienciesRel166Data2010Tight3050.push_back(105.94); 
  efficienciesRel166Data2010Tight3050.push_back(113.57);
  efficienciesRel166Data2010Tight3050.push_back(105.48); 
  efficienciesRel166Data2010Tight3050.push_back(101.99);
  uncertaintiesRel166Data2010Tight3050.push_back(3.46);
  uncertaintiesRel166Data2010Tight3050.push_back(2.65); 
  uncertaintiesRel166Data2010Tight3050.push_back(3.65);
  uncertaintiesRel166Data2010Tight3050.push_back(2.49); 
  uncertaintiesRel166Data2010Tight3050.push_back(2.33); 
  uncertaintiesRel166Data2010Tight3050.push_back(2.28);
  uncertaintiesRel166Data2010Tight3050.push_back(2.45); 
  uncertaintiesRel166Data2010Tight3050.push_back(3.72); 
  uncertaintiesRel166Data2010Tight3050.push_back(3.38); 
  uncertaintiesRel166Data2010Tight3050.push_back(2.70);

  //Probes between 20 and 50 GeV
  //Medium
  efficienciesRel166Data2010Medium2050.push_back(97.35); 
  efficienciesRel166Data2010Medium2050.push_back(95.86);
  efficienciesRel166Data2010Medium2050.push_back(96.25); 
  efficienciesRel166Data2010Medium2050.push_back(95.80);
  efficienciesRel166Data2010Medium2050.push_back(96.01);
  efficienciesRel166Data2010Medium2050.push_back(96.84); 
  efficienciesRel166Data2010Medium2050.push_back(96.04); 
  efficienciesRel166Data2010Medium2050.push_back(96.54); 
  efficienciesRel166Data2010Medium2050.push_back(96.59); 
  efficienciesRel166Data2010Medium2050.push_back(96.33);
  uncertaintiesRel166Data2010Medium2050.push_back(2.21);
  uncertaintiesRel166Data2010Medium2050.push_back(2.25);
  uncertaintiesRel166Data2010Medium2050.push_back(3.22); 
  uncertaintiesRel166Data2010Medium2050.push_back(2.27);
  uncertaintiesRel166Data2010Medium2050.push_back(2.23); 
  uncertaintiesRel166Data2010Medium2050.push_back(2.13); 
  uncertaintiesRel166Data2010Medium2050.push_back(2.17); 
  uncertaintiesRel166Data2010Medium2050.push_back(3.20); 
  uncertaintiesRel166Data2010Medium2050.push_back(2.24);
  uncertaintiesRel166Data2010Medium2050.push_back(2.41);
  //Tight
  efficienciesRel166Data2010Tight2050.push_back(99.90); 
  efficienciesRel166Data2010Tight2050.push_back(103.11); 
  efficienciesRel166Data2010Tight2050.push_back(116.16);
  efficienciesRel166Data2010Tight2050.push_back(105.70); 
  efficienciesRel166Data2010Tight2050.push_back( 97.98); 
  efficienciesRel166Data2010Tight2050.push_back( 99.08);
  efficienciesRel166Data2010Tight2050.push_back(105.23); 
  efficienciesRel166Data2010Tight2050.push_back(115.12);
  efficienciesRel166Data2010Tight2050.push_back(104.91); 
  efficienciesRel166Data2010Tight2050.push_back(101.99);
  uncertaintiesRel166Data2010Tight2050.push_back(2.28);
  uncertaintiesRel166Data2010Tight2050.push_back(2.89); 
  uncertaintiesRel166Data2010Tight2050.push_back(4.35);
  uncertaintiesRel166Data2010Tight2050.push_back(2.72); 
  uncertaintiesRel166Data2010Tight2050.push_back(2.40); 
  uncertaintiesRel166Data2010Tight2050.push_back(2.24);
  uncertaintiesRel166Data2010Tight2050.push_back(2.48); 
  uncertaintiesRel166Data2010Tight2050.push_back(4.17); 
  uncertaintiesRel166Data2010Tight2050.push_back(2.45); 
  uncertaintiesRel166Data2010Tight2050.push_back(3.29);
  //For the ET-corrections of the scale factors
  //Medium
  ETCorrectionsMediumRel166Data2010.push_back(80.60);
  ETCorrectionsMediumRel166Data2010.push_back(92.07);
  ETCorrectionsMediumRel166Data2010.push_back(96.34);
  ETCorrectionsMediumRel166Data2010.push_back(100.19);
  ETCorrectionsMediumRel166Data2010.push_back(101.54);
  ETCorrectionsMediumRel166Data2010.push_back(101.25);
  ETCorrectionsMediumRel166Data2010.push_back(102.29);
  uncertaintiesETCorrectionsMediumRel166Data2010.push_back( 9.60);
  uncertaintiesETCorrectionsMediumRel166Data2010.push_back( 3.27);
  uncertaintiesETCorrectionsMediumRel166Data2010.push_back( 1.40);
  uncertaintiesETCorrectionsMediumRel166Data2010.push_back( 0.70);
  uncertaintiesETCorrectionsMediumRel166Data2010.push_back( 0.53);
  uncertaintiesETCorrectionsMediumRel166Data2010.push_back( 0.74);
  uncertaintiesETCorrectionsMediumRel166Data2010.push_back( 1.59);
  //Tight
  ETCorrectionsTightRel166Data2010.push_back(77.78);
  ETCorrectionsTightRel166Data2010.push_back(91.84);
  ETCorrectionsTightRel166Data2010.push_back(95.67);
  ETCorrectionsTightRel166Data2010.push_back(100.86);
  ETCorrectionsTightRel166Data2010.push_back(101.83);
  ETCorrectionsTightRel166Data2010.push_back(101.33);
  ETCorrectionsTightRel166Data2010.push_back(102.10);
  uncertaintiesETCorrectionsTightRel166Data2010.push_back(10.29);
  uncertaintiesETCorrectionsTightRel166Data2010.push_back( 3.47);
  uncertaintiesETCorrectionsTightRel166Data2010.push_back( 1.52);
  uncertaintiesETCorrectionsTightRel166Data2010.push_back( 1.04);
  uncertaintiesETCorrectionsTightRel166Data2010.push_back( 0.66);
  uncertaintiesETCorrectionsTightRel166Data2010.push_back( 0.92);
  uncertaintiesETCorrectionsTightRel166Data2010.push_back( 1.90);


  //Release 16.6 Data 2011 EPS recommendations
  //Identification for probes between 20 and 50 GeV
  //Medium
  efficienciesRel166EPSMedium2050.push_back(95.7273);
  efficienciesRel166EPSMedium2050.push_back(95.5243);
  efficienciesRel166EPSMedium2050.push_back(96.403);
  efficienciesRel166EPSMedium2050.push_back(96.3494);
  efficienciesRel166EPSMedium2050.push_back(97.9518);
  efficienciesRel166EPSMedium2050.push_back(96.3292);
  efficienciesRel166EPSMedium2050.push_back(97.0952);
  efficienciesRel166EPSMedium2050.push_back(96.3317);
  efficienciesRel166EPSMedium2050.push_back(97.1977);
  efficienciesRel166EPSMedium2050.push_back(97.8678);
  efficienciesRel166EPSMedium2050.push_back(96.5697);
  efficienciesRel166EPSMedium2050.push_back(96.7783);
  efficienciesRel166EPSMedium2050.push_back(97.0532);
  efficienciesRel166EPSMedium2050.push_back(96.4621);
  efficienciesRel166EPSMedium2050.push_back(95.3501);
  efficienciesRel166EPSMedium2050.push_back(97.9656);
  efficienciesRel166EPSMedium2050.push_back(96.3031);
  efficienciesRel166EPSMedium2050.push_back(97.3978);
  efficienciesRel166EPSMedium2050.push_back(95.7546);
  efficienciesRel166EPSMedium2050.push_back(97.2443);
  uncertaintiesRel166EPSMedium2050.push_back(0.758538);
  uncertaintiesRel166EPSMedium2050.push_back(1.48083);
  uncertaintiesRel166EPSMedium2050.push_back(0.778086);
  uncertaintiesRel166EPSMedium2050.push_back(0.496963);
  uncertaintiesRel166EPSMedium2050.push_back(1.0011);
  uncertaintiesRel166EPSMedium2050.push_back(0.694056);
  uncertaintiesRel166EPSMedium2050.push_back(0.603261);
  uncertaintiesRel166EPSMedium2050.push_back(0.719089);
  uncertaintiesRel166EPSMedium2050.push_back(0.635625);
  uncertaintiesRel166EPSMedium2050.push_back(0.825545);
  uncertaintiesRel166EPSMedium2050.push_back(0.777055);
  uncertaintiesRel166EPSMedium2050.push_back(0.655198);
  uncertaintiesRel166EPSMedium2050.push_back(0.736623);
  uncertaintiesRel166EPSMedium2050.push_back(0.633197);
  uncertaintiesRel166EPSMedium2050.push_back(1.04172);
  uncertaintiesRel166EPSMedium2050.push_back(0.612204);
  uncertaintiesRel166EPSMedium2050.push_back(0.47725);
  uncertaintiesRel166EPSMedium2050.push_back(1.32532);
  uncertaintiesRel166EPSMedium2050.push_back(0.74313);
  uncertaintiesRel166EPSMedium2050.push_back(1.44683);
  //Tight
  efficienciesRel166EPSTight2050.push_back( 99.9569);
  efficienciesRel166EPSTight2050.push_back( 99.1664);
  efficienciesRel166EPSTight2050.push_back(103.421 );
  efficienciesRel166EPSTight2050.push_back(102.688 );
  efficienciesRel166EPSTight2050.push_back(113.028 );
  efficienciesRel166EPSTight2050.push_back(111.078 );
  efficienciesRel166EPSTight2050.push_back(103.481 );
  efficienciesRel166EPSTight2050.push_back( 99.5783);
  efficienciesRel166EPSTight2050.push_back( 98.4303);
  efficienciesRel166EPSTight2050.push_back(100.837 );
  efficienciesRel166EPSTight2050.push_back( 99.1868);
  efficienciesRel166EPSTight2050.push_back( 98.1188);
  efficienciesRel166EPSTight2050.push_back(100.492 );
  efficienciesRel166EPSTight2050.push_back(102.816 );
  efficienciesRel166EPSTight2050.push_back(109.09  );
  efficienciesRel166EPSTight2050.push_back(113.772 );
  efficienciesRel166EPSTight2050.push_back(103.355 );
  efficienciesRel166EPSTight2050.push_back(103.454 );
  efficienciesRel166EPSTight2050.push_back( 98.4376);
  efficienciesRel166EPSTight2050.push_back(102.174 );
  uncertaintiesRel166EPSTight2050.push_back(2.82899);
  uncertaintiesRel166EPSTight2050.push_back(1.47076);
  uncertaintiesRel166EPSTight2050.push_back(2.64305);
  uncertaintiesRel166EPSTight2050.push_back(0.692373);
  uncertaintiesRel166EPSTight2050.push_back(2.0146 );
  uncertaintiesRel166EPSTight2050.push_back(0.967662);
  uncertaintiesRel166EPSTight2050.push_back(0.714802);
  uncertaintiesRel166EPSTight2050.push_back(0.807023);
  uncertaintiesRel166EPSTight2050.push_back(0.686988);
  uncertaintiesRel166EPSTight2050.push_back(1.4562);
  uncertaintiesRel166EPSTight2050.push_back(0.984975);
  uncertaintiesRel166EPSTight2050.push_back(0.703155);
  uncertaintiesRel166EPSTight2050.push_back(0.80346);
  uncertaintiesRel166EPSTight2050.push_back(0.742777);
  uncertaintiesRel166EPSTight2050.push_back(1.78409);
  uncertaintiesRel166EPSTight2050.push_back(1.13598);
  uncertaintiesRel166EPSTight2050.push_back(0.716145);
  uncertaintiesRel166EPSTight2050.push_back(2.28302);
  uncertaintiesRel166EPSTight2050.push_back(1.13891);
  uncertaintiesRel166EPSTight2050.push_back(2.02877);
  //Identification for low ET probes
  //Medium
  efficienciesRel166EPSMediumLowET.push_back(91.16);
  efficienciesRel166EPSMediumLowET.push_back(99.84);
  efficienciesRel166EPSMediumLowET.push_back( 0.00);
  efficienciesRel166EPSMediumLowET.push_back(101.4);
  efficienciesRel166EPSMediumLowET.push_back(96.76);
  efficienciesRel166EPSMediumLowET.push_back(98.11);
  efficienciesRel166EPSMediumLowET.push_back(96.75);
  efficienciesRel166EPSMediumLowET.push_back( 0.00);
  efficienciesRel166EPSMediumLowET.push_back(86.38);
  efficienciesRel166EPSMediumLowET.push_back(84.37);
  uncertaintiesRel166EPSMediumLowET.push_back(11.0);
  uncertaintiesRel166EPSMediumLowET.push_back( 8.5);
  uncertaintiesRel166EPSMediumLowET.push_back( 0.0);
  uncertaintiesRel166EPSMediumLowET.push_back(10.8);
  uncertaintiesRel166EPSMediumLowET.push_back( 6.7);
  uncertaintiesRel166EPSMediumLowET.push_back( 7.0);
  uncertaintiesRel166EPSMediumLowET.push_back( 7.2);
  uncertaintiesRel166EPSMediumLowET.push_back( 0.0);
  uncertaintiesRel166EPSMediumLowET.push_back(10.1);
  uncertaintiesRel166EPSMediumLowET.push_back(10.2);
  //Tight
  efficienciesRel166EPSTightLowET.push_back(91.67);
  efficienciesRel166EPSTightLowET.push_back(100.6);
  efficienciesRel166EPSTightLowET.push_back( 0.00);
  efficienciesRel166EPSTightLowET.push_back(101.1);
  efficienciesRel166EPSTightLowET.push_back(96.88);
  efficienciesRel166EPSTightLowET.push_back(98.14);
  efficienciesRel166EPSTightLowET.push_back(98.23);
  efficienciesRel166EPSTightLowET.push_back( 0.00);
  efficienciesRel166EPSTightLowET.push_back(86.59);
  efficienciesRel166EPSTightLowET.push_back(84.39);
  uncertaintiesRel166EPSTightLowET.push_back(10.9);
  uncertaintiesRel166EPSTightLowET.push_back( 9.6);
  uncertaintiesRel166EPSTightLowET.push_back( 0.0);
  uncertaintiesRel166EPSTightLowET.push_back(10.5);
  uncertaintiesRel166EPSTightLowET.push_back( 6.1);
  uncertaintiesRel166EPSTightLowET.push_back( 6.1);
  uncertaintiesRel166EPSTightLowET.push_back( 9.5);
  uncertaintiesRel166EPSTightLowET.push_back( 0.0);
  uncertaintiesRel166EPSTightLowET.push_back(11.3);
  uncertaintiesRel166EPSTightLowET.push_back( 8.6);
  //For the ET-corrections of the identification scale factors
  //Medium
  ETCorrectionsMediumRel166EPS.push_back( 87.0781);
  ETCorrectionsMediumRel166EPS.push_back( 90.9091);
  ETCorrectionsMediumRel166EPS.push_back( 97.3568);
  ETCorrectionsMediumRel166EPS.push_back(100.453);
  ETCorrectionsMediumRel166EPS.push_back(101.55);
  ETCorrectionsMediumRel166EPS.push_back(101.365);
  ETCorrectionsMediumRel166EPS.push_back(102.087);
  uncertaintiesETCorrectionsMediumRel166EPS.push_back(6.00538);
  uncertaintiesETCorrectionsMediumRel166EPS.push_back(2.62057);
  uncertaintiesETCorrectionsMediumRel166EPS.push_back(0.93479);
  uncertaintiesETCorrectionsMediumRel166EPS.push_back(0.94788);
  uncertaintiesETCorrectionsMediumRel166EPS.push_back(0.43064);
  uncertaintiesETCorrectionsMediumRel166EPS.push_back(0.40351);
  uncertaintiesETCorrectionsMediumRel166EPS.push_back(0.53891);
  //Tight
  ETCorrectionsTightRel166EPS.push_back( 84.3469);
  ETCorrectionsTightRel166EPS.push_back( 89.3899);
  ETCorrectionsTightRel166EPS.push_back( 97.1825);
  ETCorrectionsTightRel166EPS.push_back(100.33);
  ETCorrectionsTightRel166EPS.push_back(101.319);
  ETCorrectionsTightRel166EPS.push_back(101.238);
  ETCorrectionsTightRel166EPS.push_back(101.552);
  uncertaintiesETCorrectionsTightRel166EPS.push_back(6.52625);
  uncertaintiesETCorrectionsTightRel166EPS.push_back(2.75939);
  uncertaintiesETCorrectionsTightRel166EPS.push_back(1.6303);
  uncertaintiesETCorrectionsTightRel166EPS.push_back(1.29104);
  uncertaintiesETCorrectionsTightRel166EPS.push_back(0.420933);
  uncertaintiesETCorrectionsTightRel166EPS.push_back(0.435997);
  uncertaintiesETCorrectionsTightRel166EPS.push_back(1.05739);
  //For the low ET electrons
  //Medium
  ETCorrectionsMediumRel166EPSFullRange.push_back(0.000/0.9666);
  ETCorrectionsMediumRel166EPSFullRange.push_back(97.36/0.9666);
  ETCorrectionsMediumRel166EPSFullRange.push_back(93.55/0.9666);
  uncertaintiesETCorrectionsMediumRel166EPSFullRange.push_back(7.25/0.9666);
  uncertaintiesETCorrectionsMediumRel166EPSFullRange.push_back(7.41/0.9666);
  uncertaintiesETCorrectionsMediumRel166EPSFullRange.push_back(8.57/0.9666);
  ETCorrectionsMediumRel166EPSFullRange.push_back( 87.0781);
  ETCorrectionsMediumRel166EPSFullRange.push_back( 90.9091);
  ETCorrectionsMediumRel166EPSFullRange.push_back( 97.3568);
  ETCorrectionsMediumRel166EPSFullRange.push_back(100.453);
  ETCorrectionsMediumRel166EPSFullRange.push_back(101.55);
  ETCorrectionsMediumRel166EPSFullRange.push_back(101.365);
  ETCorrectionsMediumRel166EPSFullRange.push_back(102.087);
  uncertaintiesETCorrectionsMediumRel166EPSFullRange.push_back(9.18078);
  uncertaintiesETCorrectionsMediumRel166EPSFullRange.push_back(2.62057);
  uncertaintiesETCorrectionsMediumRel166EPSFullRange.push_back(0.93479);
  uncertaintiesETCorrectionsMediumRel166EPSFullRange.push_back(0.94788);
  uncertaintiesETCorrectionsMediumRel166EPSFullRange.push_back(0.43064);
  uncertaintiesETCorrectionsMediumRel166EPSFullRange.push_back(0.40351);
  uncertaintiesETCorrectionsMediumRel166EPSFullRange.push_back(0.53891);
  //Tight
  ETCorrectionsTightRel166EPSFullRange.push_back(0.000/0.9673);
  ETCorrectionsTightRel166EPSFullRange.push_back(105.8/0.9673);
  ETCorrectionsTightRel166EPSFullRange.push_back(98.8/0.9673);
  uncertaintiesETCorrectionsTightRel166EPSFullRange.push_back(10.24/0.9673);
  uncertaintiesETCorrectionsTightRel166EPSFullRange.push_back(10.43/0.9673);  
  uncertaintiesETCorrectionsTightRel166EPSFullRange.push_back(10.50/0.9673);
  ETCorrectionsTightRel166EPSFullRange.push_back( 84.3469);
  ETCorrectionsTightRel166EPSFullRange.push_back( 89.3899);
  ETCorrectionsTightRel166EPSFullRange.push_back( 97.1825);
  ETCorrectionsTightRel166EPSFullRange.push_back(100.33);
  ETCorrectionsTightRel166EPSFullRange.push_back(101.319);
  ETCorrectionsTightRel166EPSFullRange.push_back(101.238);
  ETCorrectionsTightRel166EPSFullRange.push_back(101.552);
  uncertaintiesETCorrectionsTightRel166EPSFullRange.push_back(10.1599);
  uncertaintiesETCorrectionsTightRel166EPSFullRange.push_back(2.75939);
  uncertaintiesETCorrectionsTightRel166EPSFullRange.push_back(1.6303);
  uncertaintiesETCorrectionsTightRel166EPSFullRange.push_back(1.29104);
  uncertaintiesETCorrectionsTightRel166EPSFullRange.push_back(0.420933);
  uncertaintiesETCorrectionsTightRel166EPSFullRange.push_back(0.435997);
  uncertaintiesETCorrectionsTightRel166EPSFullRange.push_back(1.05739);
  //Trigger efficiency scale factors
  efficienciesRel166EPSTrigger.push_back(96.5517);
  efficienciesRel166EPSTrigger.push_back(97.3861);
  efficienciesRel166EPSTrigger.push_back(98.4245);
  efficienciesRel166EPSTrigger.push_back(98.6712);
  efficienciesRel166EPSTrigger.push_back(97.7936);
  efficienciesRel166EPSTrigger.push_back(99.7033);
  efficienciesRel166EPSTrigger.push_back(98.9571);
  efficienciesRel166EPSTrigger.push_back(98.4703);
  efficienciesRel166EPSTrigger.push_back(99.3016);
  efficienciesRel166EPSTrigger.push_back(99.1186);
  efficienciesRel166EPSTrigger.push_back(99.2838);
  efficienciesRel166EPSTrigger.push_back(99.2266);
  efficienciesRel166EPSTrigger.push_back(99.709);
  efficienciesRel166EPSTrigger.push_back(99.1478);
  efficienciesRel166EPSTrigger.push_back(99.5733);
  efficienciesRel166EPSTrigger.push_back(98.9866);
  efficienciesRel166EPSTrigger.push_back(99.8198);
  efficienciesRel166EPSTrigger.push_back(97.821);
  efficienciesRel166EPSTrigger.push_back(97.862);
  efficienciesRel166EPSTrigger.push_back(97.901);
  uncertaintiesRel166EPSTrigger.push_back(0.645476);
  uncertaintiesRel166EPSTrigger.push_back(0.588429);
  uncertaintiesRel166EPSTrigger.push_back(0.432384);
  uncertaintiesRel166EPSTrigger.push_back(0.43052);
  uncertaintiesRel166EPSTrigger.push_back(0.579508);
  uncertaintiesRel166EPSTrigger.push_back(0.410817);
  uncertaintiesRel166EPSTrigger.push_back(0.457);
  uncertaintiesRel166EPSTrigger.push_back(0.515013);
  uncertaintiesRel166EPSTrigger.push_back(0.402588);
  uncertaintiesRel166EPSTrigger.push_back(0.418344);
  uncertaintiesRel166EPSTrigger.push_back(0.415669);
  uncertaintiesRel166EPSTrigger.push_back(0.404291);
  uncertaintiesRel166EPSTrigger.push_back(0.407594);
  uncertaintiesRel166EPSTrigger.push_back(0.460203);
  uncertaintiesRel166EPSTrigger.push_back(0.410275);
  uncertaintiesRel166EPSTrigger.push_back(0.53542);
  uncertaintiesRel166EPSTrigger.push_back(0.425722);
  uncertaintiesRel166EPSTrigger.push_back(0.667037);
  uncertaintiesRel166EPSTrigger.push_back(0.426163);
  uncertaintiesRel166EPSTrigger.push_back(0.976323);
  //Reco+trackquality efficiencies
  efficienciesRel166EPSRecoTrkQual.push_back( 97.59);
  efficienciesRel166EPSRecoTrkQual.push_back(100.91);
  efficienciesRel166EPSRecoTrkQual.push_back(100.91);
  efficienciesRel166EPSRecoTrkQual.push_back(100.91);
  efficienciesRel166EPSRecoTrkQual.push_back(100.91);
  efficienciesRel166EPSRecoTrkQual.push_back(100.91);
  efficienciesRel166EPSRecoTrkQual.push_back(100.91);
  efficienciesRel166EPSRecoTrkQual.push_back( 99.84);
  efficienciesRel166EPSRecoTrkQual.push_back( 99.84);
  efficienciesRel166EPSRecoTrkQual.push_back( 99.84);
  efficienciesRel166EPSRecoTrkQual.push_back( 99.84);
  efficienciesRel166EPSRecoTrkQual.push_back( 99.84);
  efficienciesRel166EPSRecoTrkQual.push_back( 99.84);
  efficienciesRel166EPSRecoTrkQual.push_back(100.91);
  efficienciesRel166EPSRecoTrkQual.push_back(100.91);
  efficienciesRel166EPSRecoTrkQual.push_back(100.91);
  efficienciesRel166EPSRecoTrkQual.push_back(100.91);
  efficienciesRel166EPSRecoTrkQual.push_back(100.91);
  efficienciesRel166EPSRecoTrkQual.push_back(100.91);
  efficienciesRel166EPSRecoTrkQual.push_back( 97.59);
  uncertaintiesRel166EPSRecoTrkQual.push_back(1.84);
  uncertaintiesRel166EPSRecoTrkQual.push_back(0.70);
  uncertaintiesRel166EPSRecoTrkQual.push_back(0.70);
  uncertaintiesRel166EPSRecoTrkQual.push_back(0.70);
  uncertaintiesRel166EPSRecoTrkQual.push_back(0.70);
  uncertaintiesRel166EPSRecoTrkQual.push_back(0.70);
  uncertaintiesRel166EPSRecoTrkQual.push_back(0.70);
  uncertaintiesRel166EPSRecoTrkQual.push_back(0.66);
  uncertaintiesRel166EPSRecoTrkQual.push_back(0.66);
  uncertaintiesRel166EPSRecoTrkQual.push_back(0.66);
  uncertaintiesRel166EPSRecoTrkQual.push_back(0.66);
  uncertaintiesRel166EPSRecoTrkQual.push_back(0.66);
  uncertaintiesRel166EPSRecoTrkQual.push_back(0.66);
  uncertaintiesRel166EPSRecoTrkQual.push_back(0.70);
  uncertaintiesRel166EPSRecoTrkQual.push_back(0.70);
  uncertaintiesRel166EPSRecoTrkQual.push_back(0.70);
  uncertaintiesRel166EPSRecoTrkQual.push_back(0.70);
  uncertaintiesRel166EPSRecoTrkQual.push_back(0.70);
  uncertaintiesRel166EPSRecoTrkQual.push_back(0.70);
  uncertaintiesRel166EPSRecoTrkQual.push_back(1.84);

  // 2011 data with rel. 17 and MC11a ("CERN Council SF")
  // for technical reasons the values are first stored in float[],
  // then converted to vector<float>
  // Raw Rel17CC Reco+TQ SF
const float Sf_RecoTrkQ_Eta[]={ 102.01, 100.67, 100.97, 100.17, 99.40, 99.16, 99.25, 100.13, 100.73, 100.57, 102.30};
const float Sf_RecoTrkQ_Eta_err[]={ 0.70, 0.57, 0.70, 0.57, 1.11, 1.16, 0.99, 0.55, 0.90, 0.60, 0.71};
  // Raw Rel17CC Identification SF
const float sfLoosePP_Combined_eta[] = {0.978162, 0.989691, 0.9892, 1.00281, 0.993113, 0.994409, 0.995224, 1.00113, 0.9927, 0.990337, 0.98053};
const float errsfLoosePP_Combined_eta[] = {0.0184629, 0.015968, 0.00871837, 0.00385742, 0.00430604, 0.00414063, 0.00707358, 0.003712, 0.00843564, 0.0164764, 0.0178917};
const float sfLoosePP_Jpsi_eta[] = {0.928469, 0.876753, 0.947689, 0.940677, 0.933882, 0.932504, 0.943054, 0.924861, 1.07193, 0.909942, 0.94223};
const float errsfLoosePP_Jpsi_eta[] = {0.0442547, 0.0651155, 0.100367, 0.0459643, 0.0318983, 0.0337912, 0.0316421, 0.0362685, 0.0843151, 0.0566668, 0.0470655};
const float sfLoosePP_Combined_pt[] = {0., 1.04564, 1.02127, 0.950536, 0.956266, 0.985196, 1.00014, 1.00734, 1.00668, 1.00266};
const float errsfLoosePP_Combined_pt[] = {1., 0.0577688, 0.0532959, 0.0192058, 0.0159554, 0.0120122, 0.00643931, 0.00608316, 0.00608894, 0.00670763};
const float sfMediumPP_Combined_eta[] = {0.956077, 0.984517, 0.9933, 0.998451, 0.998374, 1.01566, 0.999115, 0.995048, 0.9972, 0.98697, 0.957895};
const float errsfMediumPP_Combined_eta[] = {0.013147, 0.0124841, 0.00889719, 0.00400233, 0.00446367, 0.00438371, 0.00441865, 0.00390813, 0.0090824, 0.0131541, 0.0154712};
const float sfMediumPP_Jpsi_eta[] = {0.913436, 0.892599, 0.981171, 0.918171, 0.939638, 0.935174, 0.934618, 0.907705, 1.09734, 0.874291, 0.903363};
const float errsfMediumPP_Jpsi_eta[] = {0.0451658, 0.0664901, 0.10942, 0.0456451, 0.0314451, 0.0334607, 0.0309696, 0.0362455, 0.0959015, 0.0564854, 0.0509632};
const float sfMediumPP_Combined_pt[] = {0., 1.06787, 1.0114, 0.949246, 0.940358, 0.974558, 0.994974, 1.0084, 1.00916, 1.0066};
const float errsfMediumPP_Combined_pt[] = {1., 0.0569981, 0.0482483, 0.0216574, 0.0173227, 0.0114571, 0.00633696, 0.00606375, 0.00609331, 0.00677809};
const float sfTightPP_Combined_eta[] = {0.970385, 1.00039, 1.0294, 1.02121, 1.00159, 1.01284, 1.00105, 1.01674, 1.0349, 1.00659, 0.971479};
const float errsfTightPP_Combined_eta[] = {0.0144101, 0.0116894, 0.00947048, 0.00424625, 0.00453523, 0.00451146, 0.00448671, 0.00412469, 0.00987978, 0.0116082, 0.0147539};
const float sfTightPP_Jpsi_eta[] = {0.961754, 0.913472, 1.00017, 0.920565, 0.940924, 0.930151, 0.934168, 0.898207, 1.19533, 0.887737, 0.949335};
const float errsfTightPP_Jpsi_eta[] = {0.0504488, 0.0706803, 0.117501, 0.0490665, 0.0352094, 0.0382616, 0.035019, 0.0403119, 0.109184, 0.0611907, 0.055913};
const float sfTightPP_Combined_pt[] = {0., 1.067, 1.0142, 0.953088, 0.94455, 0.974825, 0.995567, 1.00683, 1.00781, 1.00327};
const float errsfTightPP_Combined_pt[] = {1., 0.0635091, 0.0501458, 0.0228872, 0.0181984, 0.0118053, 0.00635714, 0.00609709, 0.00613041, 0.00679589};
  // Raw Rel17CC Trigger SF and eff
//////////////////////////////////////
/////// trigger e20_medium
//////////////////////////////////////

///// MC Efficiencies  vs Et
const float mcEff_e20_loo1_Et[6]={0.864088, 0.89944, 0.939482, 0.97356, 0.988114, 1.0045};
const float mcEff_e20_med1_Et[6]={0.932642, 0.959893, 0.98097, 0.990982, 0.995659, 1.00106};
const float mcEff_e20_tig1_Et[6]={0.940908, 0.965677, 0.984083, 0.992407, 0.996765, 1.00114};
//////  SF vs Et
const float SF_e20_med1_Et[6]={1.00245, 0.999368, 0.998914, 0.999601, 1.00048, 0.999809};
const float  SF_e20_med1_Et_toterror[6]={0.00597661,0.00298427,0.00154714,0.00105366,0.000587446,0.000577095};
const float SF_e20_tig1_Et[6]={1.0025, 0.995979, 0.997492, 0.999587, 1.0005, 1.00004}; 
const float  SF_e20_tig1_Et_toterror[6]={0.00570995,0.00346066,0.00152597,0.000993986,0.000568664,0.000466317};
//// MC Efficiencies vs eta
const float  mcEff_e20_loo1_eta[20]={0.83473,0.949204,0.944577,0.944431,0.828458,0.964727,0.970421,0.965379,0.955024,0.950054,0.940576,0.955574,0.966962,0.971468,0.963756,0.840913,0.95043,0.944207,0.951791,0.840843};
const float  mcEff_e20_med1_eta[20]={0.872355,0.969229,0.977964,0.973087,0.849374,0.985378,0.986858,0.988503,
0.979529,0.980543,0.974007,0.980666,0.98927,0.988601,0.98482,0.860735,0.978141,0.977674,0.970493,0.875797};
const float  mcEff_e20_tig1_eta[20]={0.886963,0.975401,0.983324,0.979107,0.853658,0.988534,0.990043,
0.991176,0.981667,0.982593,0.975875,0.982813,0.992012,0.99155,0.987791,0.865031,0.983154,0.982475,
0.975659,0.888129};
//// SF vs eta
const float  SF_e20_med1_eta[20]={1.01132,0.988154,0.98865,0.987197,1.03248,1.00244,0.994201,0.981812,1.00942,
0.923586,0.97601,1.00496,1.00136,0.995827,1.00157,1.02941,0.988385,0.983446,0.99408,1.01889};
const float SF_e20_med1_eta_toterror[20]={0.0138381, 0.010173, 0.0101778, 0.0105304, 0.0112517, 0.0100553, 
0.0100303, 0.0100703, 0.0100135, 0.0108708, 0.0105664, 0.0100226, 0.0100271, 0.0100309, 
0.0100894, 0.0114349, 0.0101388, 0.010223, 0.010139, 0.0129749};
const float  SF_e20_tig1_eta[20]={1.00768,0.98708,0.988363,0.986254,1.03095,1.00141,0.993596,0.980864,
1.00878,0.926747,0.975567,1.00439,0.999926,0.995449,1.00165,1.02746,0.989032,0.984948,0.994878,1.01632};
const float SF_e20_tig1_eta_toterror[20]={0.0133341, 0.0102052, 0.0101942, 0.0104491, 0.0114524, 0.0100434, 
0.0100267, 0.010074, 0.0100116, 0.0109426, 0.0105505, 0.010022, 0.0100252, 0.0100245, 0.0100841, 
0.0115263, 0.0101251, 0.0101966, 0.0101346, 0.0131077};
 /////////////////////////////////////////
////////  Trigger e22medium
/////////////////////////////////////////
///// MC efficiencies vs Et

float mcEff_e22_loo1_Et[6]={0., 0.877805, 0.933197, 0.973957, 0.990786, 1.00939}; 
float mcEff_e22_med1_Et[6]={0., 0.938168, 0.97598, 0.990546, 0.996888, 1.00369}; 
float mcEff_e22_tig1_Et[6]={0., 0.945008, 0.978626, 0.991608, 0.997479, 1.00331};
///   SF vs Et
const float SF_e22_med1_Et[6]={0., 1.00106, 0.997813, 1.00152, 1.00105, 0.999557};
const float  SF_e22_med1_Et_toterror[6]={1.,0.00788436,0.00348746,0.00196079,0.00128428,0.000730529};
const float SF_e22_tig1_Et[6]={0., 0.997016, 0.996317, 1.00213, 1.00131, 0.999715};
const float  SF_e22_tig1_Et_toterror[6]={1.,0.00795901,0.00301526,0.00186497,0.00135485,0.000586026};
///  MC efficiencies vs eta
const float  mcEff_e22_loo1_eta[20]={0.80987,0.935756,0.930274,0.933041,0.774026,0.955098,0.962545,0.958426,0.946747,0.940417,0.932461,0.946313,0.957602,0.961308,0.951713,0.796754,0.937367,0.935115,0.938189,0.820834};
const float  mcEff_e22_med1_eta[20]={0.850781,0.957769,0.965943,0.962989,0.796618,0.976691,0.98031,0.981781,0.9729,0.973739,0.966018,0.972751,0.981814,0.98023,0.973673,0.818503,0.967127,0.968771,0.959493,0.855468};
const float  mcEff_e22_tig1_eta[20]={0.866554,0.963796,0.97137,0.968614,0.800967,0.98014,0.983543,
0.98483,0.975155,0.975365,0.967583,0.974961,0.984649,0.983594,0.977511,0.821548,0.972552,
0.973661,0.965113,0.867067};
/// SF vs eta
const float  SF_e22_med1_eta[20]={1.0429,0.993361,0.990606,0.983569,1.07278,1.00356,0.99341,0.983135,
1.00858,0.922439,0.975137,1.00435,1.00485,0.998841,1.00251,1.05195,0.988802,0.974716,0.998945,1.03681};
const float SF_e22_med1_eta_toterror[20]={0.0170972, 0.0104995, 0.010535, 0.0113039, 0.0164707, 0.0101648, 
0.0101347, 0.0102433, 0.0100746, 0.01125, 0.0113536, 0.0100524, 0.010093, 0.0100712, 0.0102185, 
0.0129822, 0.0106089, 0.0108853, 0.0103578, 0.0177669};
const float  SF_e22_tig1_eta[20]={1.04354,0.9923,0.991662,0.983054,1.07381,1.00123,0.993014,0.982268,
1.00798,0.923537,0.977255,1.00411,1.00372,0.997861,1.00232,1.05236,0.988432,0.977005,0.997971,1.03317};
const float SF_e22_tig1_eta_toterror[20]={0.0167165, 0.0104359, 0.0106771, 0.0114704, 0.0167785, 0.0101834, 
0.0101823, 0.0102809, 0.0100571, 0.0113901, 0.011439, 0.0100538, 0.0100842, 0.0100766, 0.0101922, 
0.0131861, 0.0105969, 0.010921, 0.0103099, 0.0168021};
 
//////////////////////////////////////////////
///////  trigger e22vh medium1
////////////////////////////////////////////
//////////  MC efficiencies vs Et

const float mcEff_e22vh_loo1_Et[6]={0., 0.867613, 0.925484, 0.971542, 0.996169, 1.02254};
const float mcEff_e22vh_med1_Et[6]={0., 0.935306, 0.976425, 0.989257, 0.998406, 1.00712};
const float mcEff_e22vh_tig1_Et[6]={0., 0.946316, 0.984708, 0.993602, 1.00022, 1.00552};
///  SF vs Et
const float SF_e22vh_med1_Et[6]={0., 0.976255, 0.990213, 1.00065, 0.999608, 1.00088};
const float  SF_e22vh_med1_Et_toterror[6]={1.,0.00839289,0.00617629,0.00476838,0.00253987,0.00111677};
const float SF_e22vh_tig1_Et[6]={0., 0.971424, 0.986322, 0.998, 0.999025, 1.00115};
const float  SF_e22vh_tig1_Et_toterror[6]={1.,0.00840815,0.00611464,0.00418324,0.00233599,0.00100704};

/// MC efficiencies vs eta
const float  mcEff_e22vh_loo1_eta[20]={0.70588,0.850925,0.885496,0.88813,0.708504,0.898123,0.918406,0.916392,0.898269,0.815695,0.81267,0.897603,0.920762,0.919525,0.902682,0.726406,0.872302,0.884492,0.855155,0.72748};

const float  mcEff_e22vh_med1_eta[20]={0.806431,0.934885,0.950239,0.936858,0.783633,0.948002,0.9631,0.960062,
0.94823,0.909058,0.905434,0.948098,0.967413,0.963304,0.946392,0.793201,0.933363,0.945796,0.94782,0.828012};
const float  mcEff_e22vh_tig1_eta[20]={0.840512,0.949399,0.959665,0.95542,0.801649,0.969443,0.972239,
0.970461,0.959291,0.956205,0.956761,0.960181,0.977778,0.973809,0.965969,0.804513,0.953652,0.956252,
0.957446,0.847437};
// SF vs eta
const float  SF_e22vh_med1_eta[20]={0.984147,0.980365,0.970567,0.984624,0.97203,1.01202,0.999753,0.991051,
1.02403,0.990278,1.02291,1.01998,1.00556,1.00275,1.00858,1.02538,0.993383,0.965577,0.973939,0.953943};
const float SF_e22vh_med1_eta_toterror[20]={0.0305649, 0.0112576, 0.0118168, 0.0127388, 0.0190914, 0.0112009, 
0.0103002, 0.0105304, 0.0102212, 0.0128408, 0.0125235, 0.0103448, 0.0104107, 0.010408, 0.0108098, 
0.016527, 0.0125478, 0.0111028, 0.0107405, 0.0223372};
const float  SF_e22vh_tig1_eta[20]={0.964279,0.977745,0.975187,0.978157,0.961275,0.999952,0.995353,0.984084,
1.01769,0.956689,0.985402,1.01267,0.998659,0.99684,0.997916,1.0241,0.985229,0.970021,0.976287,0.947805};
const float SF_e22vh_tig1_eta_toterror[20]={0.033432, 0.0204931, 0.0208272, 0.0212713, 0.026605, 0.0203521, 
0.0201282, 0.020186, 0.0200903, 0.0213355, 0.0208187, 0.020117, 0.0201333, 0.0201497, 0.0203914, 
0.0245803, 0.0209498, 0.0205718, 0.0205268, 0.0272371};

/*
////////////////// NEW MC11a


const  float mcEff_e20_loo1_eta[]={0.8347, 0.949163, 0.944883, 0.944589, 0.828705, 0.964345, 
0.97052, 0.965388, 0.95503, 0.950373, 0.940983, 0.955517, 0.966908, 0.971548, 0.963624, 
0.84242, 0.949762, 0.944818, 0.951371, 0.844253}; 
 const float mcEff_e20_med1_eta[]={0.872217, 0.969157, 0.977825, 0.972825, 0.849464, 0.985315, 
0.987004, 0.988374, 0.979534, 0.979508, 0.972804, 0.980607, 0.989346, 0.988724, 0.984624, 
0.861959, 0.977611, 0.978006, 0.970531, 0.878625}; 
const float mcEff_e20_tig1_eta[]={0.886826, 0.975179, 0.983357, 0.97897, 0.853501, 0.98854, 
0.990268, 0.991087, 0.981716, 0.982607, 0.976196, 0.982777, 0.992041, 0.991794, 0.987693, 
0.866187, 0.982831, 0.983067, 0.975679, 0.890348}; 
const float mcEff_e20_loo1_Et[]={0.868757, 0.906014, 0.944793, 0.979182, 0.993821, 1.01042};
const float mcEff_e20_med1_Et[]={0.936213, 0.962769, 0.983821, 0.99346, 0.998437, 1.00404};
const float mcEff_e20_tig1_Et[]={0.942965, 0.967463, 0.98611, 0.994182, 0.998706, 1.0033};

const float  SF_e20_loo1_eta[20]={1.01534,0.984086,0.992276,0.968912,1.02685,0.991268,0.988064,0.977688,
1.00953,0.905909,0.967138,1.00372,0.996045,0.989519,0.988762,1.03063,0.973772,0.991522,0.98947,1.01551};
const float SF_e20_loo1_eta_toterror[]={0.0130578, 0.00161489, 0.001948, 0.00672797, 0.00534218, 0.00150788, 0.00118971, 0.00156628, 0.000765469, 0.00338748, 0.00294344, 0.000816807, 0.00124604, 0.00103267, 0.00137789, 0.00434045, 0.0050157, 0.00214738, 0.00156567, 0.0118714}; 
const  float SF_e20_med1_eta[]={1.0117, 0.987966, 0.988809, 0.986705, 1.03136, 1.00264, 0.994252, 0.981866, 1.00938, 0.905025, 0.967256, 1.00502, 1.00116, 0.995734, 1.00169, 1.02746, 0.988488, 0.982986, 0.993766, 1.01487}; 
const float SF_e20_med1_eta_toterror[]={0.0121997, 0.00137878, 0.00147188, 0.00255856, 0.00392385, 
0.0007532, 0.000626649, 0.000961292, 0.000422854, 0.00326331, 0.00263292, 0.000514347, 
0.000603586, 0.000632032, 0.000899679, 0.0038529, 0.00126918, 0.00160459, 0.00128164, 0.0113221}; 
const float SF_e20_tig1_eta[]={1.00802, 0.986775, 0.98832, 0.985623, 1.02958, 1.00147, 0.993581, 
0.98081, 1.00873, 0.926821, 0.975146, 1.00443, 0.999768, 0.995227, 1.00175, 1.02615, 0.989041, 
0.984289, 0.994168, 1.01223}; 
const float SF_e20_tig1_eta_toterror[]={0.0119843, 0.00141746, 0.00154014, 0.00238675, 0.00407454, 
0.000713527, 0.000607688, 0.000965085, 0.000406054, 0.00322731, 0.00237672, 0.00050872, 
0.000579617, 0.000582791, 0.000869127, 0.00389916, 0.00120008, 0.00154486, 0.00124813, 0.0114889}; 
 
const float SF_e20_loo1_Et[]={0.974092, 0.981163, 0.991234, 0.996109, 0.998977, 1.00161};
const float SF_e20_loo1_Et_toterror[]={0.0075647, 0.00481608, 0.00320545, 0.00218902, 0.00105023, 
0.000459886};
const float SF_e20_med1_Et[]={1.00018, 0.999127, 0.998588, 0.999919, 1.0006, 0.999954};
const float SF_e20_med1_Et_toterror[]={0.00463147, 0.00261407, 0.00110315, 0.000750933, 0.000517627, 
0.000325451};
const float SF_e20_tig1_Et[]={1.00141, 0.995718, 0.997332, 0.999669, 1.00066, 1.0001};
const float SF_e20_tig1_Et_toterror[]={0.00441438, 0.0027716, 0.00106634, 0.000718039, 0.000501996, 
0.000325547};

////  e22medium


const float  mcEff_e22_loo1_eta[20]={0.832022,0.948435,0.942234,0.944573,0.788029,0.963523,0.969334,
0.965267,0.953677,0.948645,0.940754,0.953558,0.964152,0.968087,0.960429,0.810809,0.94815,
0.945998,0.95014,0.843545};
const float mcEff_e22_med1_eta[20] = {0.85961, 0.97046, 0.97739, 0.97383, 0.80818, 0.98484, 0.98707, 
0.98873, 0.97817, 0.97512, 0.97110, 0.97991, 0.98828, 0.98567, 0.97380, 0.83092, 0.97625, 
0.97433, 0.97088, 0.88236}; 
const float mcEff_e22_tig1_eta[20] = {0.89017, 0.97628, 0.98268, 0.97931, 0.81145, 0.98831, 0.99030, 
0.99143, 0.98106, 0.97664, 0.97423, 0.98199, 0.99078, 0.98849, 0.97871, 0.83364, 0.98155, 0.98400, 
0.97624, 0.89184}; 
const float mcEff_e22_loo1_Et[]={0., 0.878219, 0.934702, 0.975375, 0.992215, 1.01084};
const float mcEff_e22_med1_Et[]={0., 0.939262, 0.976914, 0.99151, 0.99784, 1.00449};
const float mcEff_e22_tig1_Et[]={0., 0.945666, 0.979538, 0.992275, 0.99805, 1.00381};

const float  SF_e22_loo1_eta[20]={1.05578,0.98749,0.993409,0.967593,1.06721,0.98955,0.989584,0.977458,
1.01024,0.907711,0.964577,1.00326,0.999489,0.992449,0.990272,1.05787,0.976958,
0.985512,0.996182,1.03177};
const float SF_e22_loo1_eta_toterror[]={0.0153876, 0.00283897, 0.00365395, 0.00647611, 0.00888457, 0.00226249, 0.0016356, 0.00242826, 0.00124083, 0.00488543, 0.0047627, 0.00139215, 0.00187997, 0.00149441, 0.00244067, 0.00712401, 0.00438553, 0.00367472, 0.00245468, 0.0144138};
const float SF_e22_med1_eta[]={1.04513, 0.990634, 0.988736, 0.98363, 1.06736, 1.00261, 0.99458, 
0.982739, 1.00856, 0.905826, 0.962418, 1.00457, 1.0033, 0.998387, 1.00061, 1.05027, 0.989164, 
0.975638, 0.998642, 1.02912}; 
const float SF_e22_med1_eta_toterror[]={0.013873, 0.00236579, 0.00251913, 0.00398823, 0.00782065, 
0.00125015, 0.00109308, 0.00162961, 0.000724821, 0.00434826, 0.0037365, 0.000812905, 0.000981582, 
0.000955442, 0.0015643, 0.00653014, 0.002118, 0.00283932, 0.00198101, 0.0136735}; 
const float SF_e22_tig1_eta[]={1.04489, 0.988417, 0.990143, 0.983469, 1.06933, 1.00059, 0.993938, 
0.98201, 1.00796, 0.92504, 0.975935, 1.0042, 1.00241, 0.997385, 1.00008, 1.05144, 0.988913, 
0.977766, 0.99687, 1.02637}; 
const float SF_e22_tig1_eta_toterror[]={0.0138419, 0.0022802, 0.0024961, 0.00406588, 0.00817411, 
0.00127067, 0.0011713, 0.00161794, 0.000705816, 0.00420083, 0.00358342, 0.000816685, 0.000896891, 
0.00090756, 0.00145694, 0.00689821, 0.00215051, 0.00277435, 0.00198534, 0.0137613}; 

const float SF_e22_loo1_Et[]={0., 0.98665, 0.989572, 0.997938, 0.999318, 1.00108};
const float SF_e22_loo1_Et_toterror[]={1, 0.0081315, 0.00410451, 0.00235089, 0.00126392, 0.00069558};
const float SF_e22_med1_Et[]={0.,0.999837, 0.995717, 0.998179, 1.00083, 1.00035};
const float SF_e22_med1_Et_toterror[]={1.,0.0544177, 0.0333886, 0.0231181, 0.0126187, 0.00601838};
const float SF_e22_tig1_Et[]={0.,0.999476, 0.994537, 0.999034, 1.00004, 1.0005};
const float SF_e22_tig1_Et_toterror[]={1.,0.0304421, 0.0167658, 0.0165566, 0.00902732, 0.00431354};

///////   e22vh medium
const float mcEff_e22vh_loo1_eta[20] = {0.723515,0.863297,0.892337,0.897324,0.719559,0.903321,0.925592,0.921486,0.906042,0.827351,0.815623,0.904621,0.925475,0.927009,0.91249,0.733977,0.889296,0.895898,0.868753,0.756952};
const float  mcEff_e22vh_med1_eta[20]={0.827976,0.946008,0.95749,0.944501,0.793274,0.947288,0.964607,0.963974,0.952484,0.881975,0.873494,0.951732,0.970088,0.96586,0.952049,0.798397,0.944486,0.954791,0.956044,0.85463};
const float  mcEff_e22vh_tig1_eta[20]={0.859277,0.959705,0.969998,0.962972,0.809268,0.974332,0.977381,0.974905,0.96647,0.96496,0.959858,0.966156,0.981071,0.979541,0.974038,0.812969,0.965375,0.966633,0.966164,0.874814};

const float mcEff_e22vh_loo1_Et[]={0., 0.860133, 0.920322, 0.966649, 0.991088, 1.01929};
const float mcEff_e22vh_med1_Et[]={0., 0.934748, 0.975172, 0.988596, 0.997157, 1.00696};
const float mcEff_e22vh_tig1_Et[]={0., 0.947078, 0.98204, 0.991989, 0.997974, 1.00466};

const  float SF_e22vh_loo1_eta[]={0.96539, 0.951172, 0.979968, 0.966965, 0.985037, 1.01246, 0.999731, 
0.986406, 1.02924, 1.01289, 1.05758, 1.02527, 1.00649, 0.999593, 0.998486, 1.03959, 0.984195, 
0.968464, 0.952544, 0.92481}; 
const float SF_e22vh_loo1_eta_toterror[]={0.0220327, 0.00532938, 0.00587379, 0.00848225, 0.0126295, 
0.00422658, 0.00279119, 0.0034246, 0.0023838, 0.00809234, 0.00858341, 0.00242839, 0.00338899, 
0.00275655, 0.00441373, 0.0117198, 0.0061174, 0.00548657, 0.00511056, 0.0188316}; 
 const float SF_e22vh_med1_eta[]={0.975829, 0.976335, 0.971136, 0.985203, 0.97232, 1.01855, 1.00378, 
0.991932, 1.02498, 1.0054, 1.04892, 1.02104, 1.00689, 1.0043, 1.00784, 1.02998, 0.991559, 
0.963973, 0.974938, 0.938505}; 
const float SF_e22vh_med1_eta_toterror[]={0.0182564, 0.00343797, 0.00369781, 0.00480054, 0.0103036, 
0.003283, 0.00193349, 0.00242209, 0.00166597, 0.00680114, 0.00710823, 0.0017541, 0.00216051, 0.00189069, 
0.00338253, 0.00995842, 0.00385163, 0.00371609, 0.00299345, 0.0160339}; 
const float SF_e22vh_tig1_eta[]={0.959364, 0.974746, 0.972884, 0.978999, 0.965078, 1.00193, 0.997124, 
0.985397, 1.01645, 0.954836, 0.988628, 1.01212, 0.999971, 0.996856, 0.995909, 1.02474, 0.983732, 
0.967716, 0.976645, 0.933241}; 
const float SF_e22vh_tig1_eta_toterror[]={0.0182264, 0.00319429, 0.00348506, 0.00507881, 0.0105802, 
0.00239738, 0.00163221, 0.0021311, 0.0014303, 0.00407252, 0.00401897, 0.00146839, 0.00173545, 
0.00150506, 0.00242264, 0.0111007, 0.00340827, 0.00350712, 0.00294779, 0.0164273}; 

const float SF_e22vh_loo1_Et[]={0., 0.941736, 0.975829, 0.993597, 0.998106, 1.00442};
const float SF_e22vh_loo1_Et_toterror[]={1., 0.00961248, 0.00688575, 0.00435265, 0.00241682, 0.00124295};
const float SF_e22vh_med1_Et[]={0, 0.974353, 0.994114, 0.999104, 1.0007, 1.00082};
const float SF_e22vh_med1_Et_toterror[]={1., 0.0586086, 0.0200998, 0.0191683, 0.0107912, 0.00893949};
const float SF_e22vh_tig1_Et[]={0, 0.972387, 0.991234, 0.997235, 1, 1.00153};
const float SF_e22vh_tig1_Et_toterror[]={1., 0.035747, 0.0181833, 0.016924, 0.00904792, 0.00446041};


*/

 copyToVector(Sf_RecoTrkQ_Eta, 11, efficienciesRel17CCRecoTrkQual, 1.);
 copyToVector(Sf_RecoTrkQ_Eta_err, 11, uncertaintiesRel17CCRecoTrkQual, 1.);
 
 // Identification eta for probes between 15 and 50 GeV
 copyToVector(sfLoosePP_Combined_eta, 11, efficienciesRel17CCLoosePP1550);
 copyToVector(errsfLoosePP_Combined_eta, 11, uncertaintiesRel17CCLoosePP1550);
 copyToVector(sfMediumPP_Combined_eta, 11, efficienciesRel17CCMediumPP1550);
 copyToVector(errsfMediumPP_Combined_eta, 11, uncertaintiesRel17CCMediumPP1550);
 copyToVector(sfTightPP_Combined_eta, 11, efficienciesRel17CCTightPP1550);
 copyToVector(errsfTightPP_Combined_eta, 11, uncertaintiesRel17CCTightPP1550);
 //Identification eta for low ET probes
 copyToVector(sfLoosePP_Jpsi_eta, 11, efficienciesRel17CCLoosePP415);
 copyToVector(errsfLoosePP_Jpsi_eta, 11, uncertaintiesRel17CCLoosePP415);
 copyToVector(sfMediumPP_Jpsi_eta, 11, efficienciesRel17CCMediumPP415);
 copyToVector(errsfMediumPP_Jpsi_eta, 11, uncertaintiesRel17CCMediumPP415);
 copyToVector(sfTightPP_Jpsi_eta, 11, efficienciesRel17CCTightPP415);
 copyToVector(errsfTightPP_Jpsi_eta, 11, uncertaintiesRel17CCTightPP415);
 // ET correction
 copyToVector(sfLoosePP_Combined_pt, 10, ETCorrectionsRel17CCLoosePP);
 copyToVector(errsfLoosePP_Combined_pt, 10, uncertaintiesETCorrectionsRel17CCLoosePP);
 copyToVector(sfMediumPP_Combined_pt, 10, ETCorrectionsRel17CCMediumPP);
 copyToVector(errsfMediumPP_Combined_pt, 10, uncertaintiesETCorrectionsRel17CCMediumPP);
 copyToVector(sfTightPP_Combined_pt, 10, ETCorrectionsRel17CCTightPP);
 copyToVector(errsfTightPP_Combined_pt, 10, uncertaintiesETCorrectionsRel17CCTightPP);
 
 // Trigger efficiencies
 // e20_medium B-J
  copyToVector(SF_e20_med1_eta, 20, efficienciesRel17CCe20_mediumMediumPP);
  copyToVector(SF_e20_med1_eta_toterror, 20, uncertaintiesRel17CCe20_mediumMediumPP);
  copyToVector(SF_e20_med1_Et, 6, efficienciesRel17CCe20_mediumMediumPPET);
  copyToVector(SF_e20_med1_Et_toterror, 6, uncertaintiesRel17CCe20_mediumMediumPPET);

  copyToVector(mcEff_e20_med1_eta, 20, MCefficienciesRel17CCe20_mediumMediumPP);
  copyToVector(mcEff_e20_med1_Et, 6, MCefficienciesRel17CCe20_mediumMediumPPET);

  copyToVector(SF_e20_tig1_eta, 20, efficienciesRel17CCe20_mediumTightPP);
  copyToVector(SF_e20_tig1_eta_toterror, 20, uncertaintiesRel17CCe20_mediumTightPP);
  copyToVector(SF_e20_tig1_Et, 6, efficienciesRel17CCe20_mediumTightPPET);
  copyToVector(SF_e20_tig1_Et_toterror, 6, uncertaintiesRel17CCe20_mediumTightPPET);

  copyToVector(mcEff_e20_tig1_eta, 20, MCefficienciesRel17CCe20_mediumTightPP);
  copyToVector(mcEff_e20_tig1_Et, 6, MCefficienciesRel17CCe20_mediumTightPPET);

  copyToVector(mcEff_e20_loo1_eta, 20, MCefficienciesRel17CCe20_mediumLoosePP);
  copyToVector(mcEff_e20_loo1_Et, 6, MCefficienciesRel17CCe20_mediumLoosePPET);


  // e22_medium K
  copyToVector(SF_e22_med1_eta, 20, efficienciesRel17CCe22_mediumMediumPP);
  copyToVector(SF_e22_med1_eta_toterror, 20, uncertaintiesRel17CCe22_mediumMediumPP);
  copyToVector(SF_e22_med1_Et, 6, efficienciesRel17CCe22_mediumMediumPPET);
  copyToVector(SF_e22_med1_Et_toterror, 6, uncertaintiesRel17CCe22_mediumMediumPPET);

  copyToVector(mcEff_e22_med1_eta, 20, MCefficienciesRel17CCe22_mediumMediumPP);
  copyToVector(mcEff_e22_med1_Et, 6, MCefficienciesRel17CCe22_mediumMediumPPET);

  copyToVector(SF_e22_tig1_eta, 20, efficienciesRel17CCe22_mediumTightPP);
  copyToVector(SF_e22_tig1_eta_toterror, 20, uncertaintiesRel17CCe22_mediumTightPP);
  copyToVector(SF_e22_tig1_Et, 6, efficienciesRel17CCe22_mediumTightPPET);
  copyToVector(SF_e22_tig1_Et_toterror, 6, uncertaintiesRel17CCe22_mediumTightPPET);

  copyToVector(mcEff_e22_tig1_eta, 20, MCefficienciesRel17CCe22_mediumTightPP);
  copyToVector(mcEff_e22_tig1_Et, 6, MCefficienciesRel17CCe22_mediumTightPPET);

  copyToVector(mcEff_e22_loo1_eta, 20, MCefficienciesRel17CCe22_mediumLoosePP);
  copyToVector(mcEff_e22_loo1_Et, 6, MCefficienciesRel17CCe22_mediumLoosePPET);


  // e22vh_medium1 L-M
  copyToVector(SF_e22vh_med1_eta, 20, efficienciesRel17CCe22vh_medium1MediumPP);
  copyToVector(SF_e22vh_med1_eta_toterror, 20, uncertaintiesRel17CCe22vh_medium1MediumPP);
  copyToVector(SF_e22vh_med1_Et, 6, efficienciesRel17CCe22vh_medium1MediumPPET);
  copyToVector(SF_e22vh_med1_Et_toterror, 6, uncertaintiesRel17CCe22vh_medium1MediumPPET);

  copyToVector(mcEff_e22vh_med1_eta, 20, MCefficienciesRel17CCe22vh_medium1MediumPP);
  copyToVector(mcEff_e22vh_med1_Et, 6, MCefficienciesRel17CCe22vh_medium1MediumPPET);

  copyToVector(SF_e22vh_tig1_eta, 20, efficienciesRel17CCe22vh_medium1TightPP);
  copyToVector(SF_e22vh_tig1_eta_toterror, 20, uncertaintiesRel17CCe22vh_medium1TightPP);
  copyToVector(SF_e22vh_tig1_Et, 6, efficienciesRel17CCe22vh_medium1TightPPET);
  copyToVector(SF_e22vh_tig1_Et_toterror, 6, uncertaintiesRel17CCe22vh_medium1TightPPET);

  copyToVector(mcEff_e22vh_tig1_eta, 20, MCefficienciesRel17CCe22vh_medium1TightPP);
  copyToVector(mcEff_e22vh_tig1_Et, 6, MCefficienciesRel17CCe22vh_medium1TightPPET);


  copyToVector(mcEff_e22vh_loo1_eta, 20, MCefficienciesRel17CCe22vh_medium1LoosePP);
  copyToVector(mcEff_e22vh_loo1_Et, 6, MCefficienciesRel17CCe22vh_medium1LoosePPET);


  // 2011 data with rel. 17 and MC11a/b/c G4 FullSim ("Moriond SF")
const float sfLoosePP_Combined_Moriond_eta[] = {0.983727, 0.992645, 0.9892, 1.00425, 0.994238, 0.99583, 0.995955, 1.00307, 0.9907, 0.994061, 0.984558};
const float errsfLoosePP_Combined_Moriond_eta[] = {0.0190335, 0.0158532, 0.00838153, 0.00394904, 0.00432761, 0.00425582, 0.00429913, 0.00357115, 0.0079806, 0.0164504, 0.0178845};
const float sfLoosePP_Jpsi_Moriond_eta[] = {0.928469, 0.876753, 0.947689, 0.940677, 0.933882, 0.932504, 0.943054, 0.924861, 1.07193, 0.909942, 0.94223};
const float errsfLoosePP_Jpsi_Moriond_eta[] = {0.0442547, 0.0651155, 0.100367, 0.0459643, 0.0318983, 0.0337912, 0.0316421, 0.0362685, 0.0843151, 0.0566668, 0.0470655};

const float sfLoosePP_Combined_Moriond_pt[] = {0., 1.04564, 1.02127, 0.952484, 0.958019, 0.986368, 1.00122, 1.00889, 1.00725, 1.00164};
const float errsfLoosePP_Combined_Moriond_pt[] = {1., 0.0584143, 0.0539949, 0.0220578, 0.0191925, 0.011236, 0.00582933, 0.00661811, 0.00532031, 0.00642074};


const float sfMedium_Combined_Moriond_eta[] = {0.980389, 0.984739, 0.987, 0.995007, 0.992199, 0.999625, 0.992335, 0.991465, 0.9896, 0.990364, 0.982351};
const float errsfMedium_Combined_Moriond_eta[] = {0.0187827, 0.01235, 0.00807527, 0.00418403, 0.0044624, 0.00429145, 0.00446588, 0.00383729, 0.00793473, 0.0138236, 0.0179662};
const float sfMedium_Jpsi_Moriond_eta[] = {0.942597, 0.974524, 1.08641, 1.01779, 0.975308, 0.97664, 0.982064, 0.971485, 1.07755, 0.981647, 0.910087};
const float errsfMedium_Jpsi_Moriond_eta[] = {0.042, 0.085, 0.091, 0.049, 0.028, 0.03, 0.028, 0.049, 0.095, 0.06, 0.049};
const float sfMedium_Combined_Moriond_pt[] = {0., 1.02283, 0.980082, 0.954413, 0.950976, 0.982899, 0.998409, 1.00998, 1.00825, 1.00394};
const float errsfMedium_Combined_Moriond_pt[] = {1., 0.0562668, 0.0482163, 0.0225589, 0.0206991, 0.00776764, 0.00602853, 0.00685129, 0.00541447, 0.00645837};


const float sfMediumPP_Combined_Moriond_eta[] = {0.967144, 0.990423, 0.9926, 0.998398, 1.00044, 1.01609, 1.00082, 0.995716, 0.9922, 0.994387, 0.967035};
const float errsfMediumPP_Combined_Moriond_eta[] = {0.0152302, 0.0126984, 0.00866025, 0.00419149, 0.00447153, 0.00445324, 0.00447221, 0.00391371, 0.00843564, 0.0148758, 0.0164505};
const float sfMediumPP_Jpsi_Moriond_eta[] = {0.913436, 0.892599, 0.981171, 0.918171, 0.939638, 0.935174, 0.934618, 0.907705, 1.09734, 0.874291, 0.903363};
const float errsfMediumPP_Jpsi_Moriond_eta[] = {0.0451658, 0.0664901, 0.10942, 0.0456451, 0.0314451, 0.0334607, 0.0309696, 0.0362455, 0.0959015, 0.0564854, 0.0509632};
const float sfMediumPP_Combined_Moriond_pt[] = {0., 1.06787, 1.0114, 0.952377, 0.942511, 0.977914, 0.995868, 1.00973, 1.01033, 1.0053};
const float errsfMediumPP_Combined_Moriond_pt[] = {1., 0.0576523, 0.0490194, 0.0245875, 0.0217735, 0.00821047, 0.00626121, 0.00701316, 0.00551153, 0.00653204};


const float sfTightPP_Combined_Moriond_eta[] = {0.987679, 1.00704, 1.027, 1.02319, 1.00345, 1.01298, 1.00215, 1.0204, 1.0276, 1.01574, 0.985344};
const float errsfTightPP_Combined_Moriond_eta[] = {0.00770372, 0.0116938, 0.00920923, 0.00471006, 0.00464143, 0.0051928, 0.0045983, 0.00446568, 0.0090824, 0.0140841, 0.0168724};
const float sfTightPP_Jpsi_Moriond_eta[] = {0.961754, 0.913472, 1.00017, 0.920565, 0.940924, 0.930151, 0.934168, 0.898207, 1.19533, 0.887737, 0.949335};
const float errsfTightPP_Jpsi_Moriond_eta[] = {0.0504488, 0.0706803, 0.117501, 0.0490665, 0.0352094, 0.0382616, 0.035019, 0.0403119, 0.109184, 0.0611907, 0.055913};
const float sfTightPP_Combined_Moriond_pt[] = {0., 1.067, 1.0142, 0.957576, 0.949916, 0.979087, 0.996125, 1.00763, 1.00897, 1.0016};
const float errsfTightPP_Combined_Moriond_pt[] = {1., 0.0640969, 0.0508882, 0.0259462, 0.0232052, 0.0085569, 0.00700376, 0.00719543, 0.00578251, 0.00667381};



  // reco+TQ are by choice identical to the "Cern council" results!
  copyToVector(Sf_RecoTrkQ_Eta, 11, efficienciesRel17MoriondRecoTrkQual, 1.);
  copyToVector(Sf_RecoTrkQ_Eta_err, 11, uncertaintiesRel17MoriondRecoTrkQual, 1.);
 
  // Identification eta for probes between 15 and 50 GeV
  copyToVector(sfLoosePP_Combined_Moriond_eta, 11, efficienciesRel17MoriondLoosePP1550);
  copyToVector(errsfLoosePP_Combined_Moriond_eta, 11, uncertaintiesRel17MoriondLoosePP1550);
  copyToVector(sfMedium_Combined_Moriond_eta, 11, efficienciesRel17MoriondMedium1550);
  copyToVector(errsfMedium_Combined_Moriond_eta, 11, uncertaintiesRel17MoriondMedium1550);
  copyToVector(sfMediumPP_Combined_Moriond_eta, 11, efficienciesRel17MoriondMediumPP1550);
  copyToVector(errsfMediumPP_Combined_Moriond_eta, 11, uncertaintiesRel17MoriondMediumPP1550);
  copyToVector(sfTightPP_Combined_Moriond_eta, 11, efficienciesRel17MoriondTightPP1550);
  copyToVector(errsfTightPP_Combined_Moriond_eta, 11, uncertaintiesRel17MoriondTightPP1550);
  //Identification eta for low ET probes
  copyToVector(sfLoosePP_Jpsi_Moriond_eta, 11, efficienciesRel17MoriondLoosePP415);
  copyToVector(errsfLoosePP_Jpsi_Moriond_eta, 11, uncertaintiesRel17MoriondLoosePP415);
  copyToVector(sfMedium_Jpsi_Moriond_eta, 11, efficienciesRel17MoriondMedium415);
  copyToVector(errsfMedium_Jpsi_Moriond_eta, 11, uncertaintiesRel17MoriondMedium415);
  copyToVector(sfMediumPP_Jpsi_Moriond_eta, 11, efficienciesRel17MoriondMediumPP415);
  copyToVector(errsfMediumPP_Jpsi_Moriond_eta, 11, uncertaintiesRel17MoriondMediumPP415);
  copyToVector(sfTightPP_Jpsi_Moriond_eta, 11, efficienciesRel17MoriondTightPP415);
  copyToVector(errsfTightPP_Jpsi_Moriond_eta, 11, uncertaintiesRel17MoriondTightPP415);
  // ET correction
  copyToVector(sfLoosePP_Combined_Moriond_pt, 10, ETCorrectionsRel17MoriondLoosePP);
  copyToVector(errsfLoosePP_Combined_Moriond_pt, 10, uncertaintiesETCorrectionsRel17MoriondLoosePP);
  copyToVector(sfMedium_Combined_Moriond_pt, 10, ETCorrectionsRel17MoriondMedium);
  copyToVector(errsfMedium_Combined_Moriond_pt, 10, uncertaintiesETCorrectionsRel17MoriondMedium);
  copyToVector(sfMediumPP_Combined_Moriond_pt, 10, ETCorrectionsRel17MoriondMediumPP);
  copyToVector(errsfMediumPP_Combined_Moriond_pt, 10, uncertaintiesETCorrectionsRel17MoriondMediumPP);
  copyToVector(sfTightPP_Combined_Moriond_pt, 10, ETCorrectionsRel17MoriondTightPP);
  copyToVector(errsfTightPP_Combined_Moriond_pt, 10, uncertaintiesETCorrectionsRel17MoriondTightPP);
 
 // Trigger efficiencies
 // e20_medium B-J  


const float  mcEff_e20_loo1c_eta[20]={0.848389,0.956334,0.950973,0.953063,0.832028,0.970957,0.975289,
0.971215,0.961237,0.960304,0.934168, 0.946466,0.955915,0.959373,0.954392, 0.848096,0.956117,
0.950492,0.957751,0.854038};
 const float mcEff_e20_med1c_eta[]={0.885794, 0.977261, 0.983728, 0.980814, 0.852281, 0.990963, 
0.991499, 0.991901, 0.984694, 0.987225, 0.96583, 0.970429, 0.977551, 0.975742, 0.973367, 
0.869728, 0.984179, 0.985036, 0.978288, 0.887523}; 
const float mcEff_e20_tig1c_eta[]={0.902422, 0.98492, 0.98942, 0.985926, 0.855551, 0.993568, 
0.993645, 0.994026, 0.986244, 0.989781, 0.968321, 0.972046, 0.979683, 0.978163, 0.975803, 
0.871911, 0.988501, 0.990043, 0.985543, 0.906235}; 

const float mcEff_e20_loo1c_Et[]={0.876544, 0.911258, 0.946972, 0.980841, 0.994258, 1.00986};
const float mcEff_e20_med1c_Et[]={0.943346, 0.966699, 0.984657, 0.994935, 0.999148, 1.00367};
const float mcEff_e20_tig1c_Et[]={0.949817, 0.972587, 0.98651, 0.995306, 0.999683, 1.003};
//SF MC11c

const float  SF_e20_loo1c_eta[20]={0.998949,0.976708,0.985923,0.960299,1.02274,0.984517,0.983232,
0.971822,1.00301,0.896542,0.974195,1.01332,1.0075,1.00208,0.998327,1.02373,0.967299,0.985603,
0.982879,1.00387};
const float SF_e20_loo1c_eta_toterror[20]={0.0133927, 0.00183132, 0.00237557, 0.00727585, 0.00635603, 
0.00182113, 0.00122957, 0.00174555, 0.000928082, 0.00381891, 0.00305131, 0.00109265, 0.0015499, 
0.00131237, 0.0017381, 0.00527935, 0.00418972, 0.00252536, 0.00185791, 0.0125777};
 const float SF_e20_med1c_eta[]={0.996195, 0.979772, 0.982875, 0.978667, 1.02795, 0.996923, 0.989745, 
0.978374, 1.00409, 0.89795, 0.974241, 1.01556, 1.01324, 1.00898, 1.01327, 1.01829, 0.981891, 
0.97597, 0.985885, 1.00469}; 
const float SF_e20_med1c_eta_toterror[]={0.0126418, 0.00152712, 0.00169908, 0.00279707, 0.00496428, 
0.000840994, 0.000693205, 0.00104571, 0.00056198, 0.00333248, 0.00253498, 0.000735072, 
0.000894095, 0.000821252, 0.0011766, 0.004989, 0.00136793, 0.001719, 0.00147176, 0.0120504}; 
const float SF_e20_tig1c_eta[]={0.990602, 0.977016, 0.982264, 0.978669, 1.02712, 0.9964, 0.990204, 
0.97791, 1.0041, 0.920103, 0.983076, 1.01551, 1.01238, 1.0091, 1.01396, 1.01941, 0.983368, 
0.977353, 0.984219, 0.994481}; 
const float SF_e20_tig1c_eta_toterror[]={0.0126946, 0.0015388, 0.00158688, 0.00270277, 0.00521395, 
0.00081376, 0.000654352, 0.00103651, 0.000539949, 0.00334691, 0.0024402, 0.000736938, 
0.00087863, 0.00079522, 0.00112249, 0.00498023, 0.00130398, 0.00159166, 0.00135454, 0.0120652}; 
 
const float SF_e20_loo1c_Et[]={0.965482, 0.975569, 0.989, 0.994475, 0.998588, 1.00222};
const float SF_e20_loo1c_Et_toterror[]={0.0097343, 0.00707489, 0.00621119, 0.00541167, 0.00513511, 0.00503009};
const float SF_e20_med1c_Et[]={0.992814, 0.995266, 0.997939, 0.998636, 1.00009, 1.00052};
const float SF_e20_med1c_Et_toterror[]={0.00740121, 0.00586089, 0.00520539, 0.00506951, 0.00504291, 0.00501729};
const float SF_e20_tig1c_Et[]={0.994419, 0.99071, 0.997164, 0.998778, 0.999922, 1.00064};
const float SF_e20_tig1c_Et_toterror[]={0.00717398, 0.00590631, 0.00517393, 0.0050728, 0.00504031, 0.00501678};

  copyToVector(SF_e20_loo1c_eta, 20, efficienciesRel17Morionde20_mediumLoosePP);
  copyToVector(SF_e20_loo1c_eta_toterror, 20, uncertaintiesRel17Morionde20_mediumLoosePP);
  copyToVector(SF_e20_loo1c_Et, 6, efficienciesRel17Morionde20_mediumLoosePPET);
  copyToVector(SF_e20_loo1c_Et_toterror, 6, uncertaintiesRel17Morionde20_mediumLoosePPET);

  copyToVector(SF_e20_med1c_eta, 20, efficienciesRel17Morionde20_mediumMediumPP);
  copyToVector(SF_e20_med1c_eta_toterror, 20, uncertaintiesRel17Morionde20_mediumMediumPP);
  copyToVector(SF_e20_med1c_Et, 6, efficienciesRel17Morionde20_mediumMediumPPET);
  copyToVector(SF_e20_med1c_Et_toterror, 6, uncertaintiesRel17Morionde20_mediumMediumPPET);

  copyToVector(SF_e20_tig1c_eta, 20, efficienciesRel17Morionde20_mediumTightPP);
  copyToVector(SF_e20_tig1c_eta_toterror, 20, uncertaintiesRel17Morionde20_mediumTightPP);
  copyToVector(SF_e20_tig1c_Et, 6, efficienciesRel17Morionde20_mediumTightPPET);
  copyToVector(SF_e20_tig1c_Et_toterror, 6, uncertaintiesRel17Morionde20_mediumTightPPET);

  copyToVector(mcEff_e20_loo1c_eta, 20, MCefficienciesRel17Morionde20_mediumLoosePP);
  copyToVector(mcEff_e20_loo1c_Et, 6, MCefficienciesRel17Morionde20_mediumLoosePPET);
  copyToVector(mcEff_e20_med1c_eta, 20, MCefficienciesRel17Morionde20_mediumMediumPP);
  copyToVector(mcEff_e20_med1c_Et, 6, MCefficienciesRel17Morionde20_mediumMediumPPET);
  copyToVector(mcEff_e20_tig1c_eta, 20, MCefficienciesRel17Morionde20_mediumTightPP);
  copyToVector(mcEff_e20_tig1c_Et, 6, MCefficienciesRel17Morionde20_mediumTightPPET);

  //// e22_medium K 

const float  mcEff_e22_loo1c_eta[20]={0.840584,0.953721,0.950048,0.952885,0.795771,0.972006,0.974945,0.969896,0.96188,0.956744,0.901899,0.914529,0.925091,0.925385,0.922256,0.808596,0.957785,0.951094,0.956664,0.84567};
const float mcEff_e22_med1c_eta[20] = {0.88053, 0.97644, 0.98375, 0.97965, 0.81019, 0.99031, 0.99011, 0.99125, 0.98457, 0.98497, 0.93506, 0.94347, 0.95187, 0.94754, 0.94672, 0.82641, 0.98257, 0.98388, 0.97803, 0.87764}; 
const float mcEff_e22_tig1c_eta[20] = {0.90645, 0.98682, 0.98994, 0.98426, 0.81090, 0.99335, 0.99266, 0.99259, 0.98572, 0.98822, 0.93821, 0.94517, 0.95526, 0.95105, 0.95012, 0.83115, 0.98861, 0.98896, 0.98797, 0.90253}; 

const float mcEff_e22_loo1c_Et[]={0., 0.888972, 0.939452, 0.977535, 0.994082, 1.01171};
const float mcEff_e22_med1c_Et[]={0., 0.944215, 0.977993, 0.992536, 0.999043, 1.00482};
const float mcEff_e22_tig1c_Et[]={0., 0.951214, 0.980248, 0.993901, 0.999523, 1.00402};

//
 //SF 11c 
const float  SF_e22_loo1c_eta[20]={1.04505,0.982017,0.98524,0.959152,1.05681,0.980914,0.983889,0.972793,1.00162,0.900027,1.00613,1.04607,1.04169,1.03825,1.03126,1.06077,0.96713,0.980231,0.989389,1.02917};
const float SF_e22_loo1c_eta_toterror[20]={0.0175488, 0.00292899, 0.00391901, 0.00658872, 0.010585, 0.00239365, 0.00164762, 0.0024999, 0.00131486, 0.0050669, 0.011244, 0.0101584, 0.0103494, 0.0102053, 0.010486, 0.00810613, 
0.00424727, 0.00377396, 0.00260471, 0.0152136};
const  float SF_e22_med1c_eta[]={1.03024, 0.98367, 0.981344, 0.975884, 1.05909, 0.997067, 0.990773, 
0.979462, 1.00306, 0.900019, 1.00263, 1.04896, 1.04779, 1.04571, 1.04256, 1.05124, 0.981777, 
0.970867, 0.989822, 1.02413}; 
const float SF_e22_med1c_eta_toterror[]={0.0154226, 0.00239801, 0.00274037, 0.00393086, 0.00868051, 
0.00133035, 0.0011408, 0.00169364, 0.000833678, 0.00446543, 0.0108841, 0.0100862, 0.0101717, 0.0101252, 0.0102804, 0.00747541, 0.00219934, 0.0028696, 0.00209202, 0.0146176}; 
const float SF_e22_tig1c_eta[]={1.03125, 0.97921, 0.982814, 0.976566, 1.06166, 0.9957, 0.990672, 
0.979767, 1.00323, 0.919104, 1.01636, 1.04919, 1.0467, 1.04592, 1.04398, 1.04873, 0.981864, 
0.97313, 0.986599, 1.0129}; 
const float SF_e22_tig1c_eta_toterror[]={0.0154282, 0.00235617, 0.0026783, 0.00396651, 0.00881431, 
0.00133136, 0.00124223, 0.0016658, 0.000797657, 0.0043214, 0.0108864, 0.0100887, 0.0101736, 0.0101333, 
0.0102847, 0.00765196, 0.00218572, 0.00280942, 0.00204223, 0.0147474}; 
const float SF_e22_loo1c_Et[]={0., 0.976139, 0.986007, 0.997185, 0.998897, 1.00168};
const float SF_e22_loo1c_Et_toterror[]={1., 0.010067, 0.00635541, 0.00571717, 0.0051924, 0.00506772};
const float SF_e22_med1c_Et[]={0., 0.995079, 0.996567, 1.00054, 1.00046, 1.00007};
const float SF_e22_med1c_Et_toterror[]={1, 0.00803895, 0.00556503, 0.00525804, 0.00510948, 0.00505457};
const float SF_e22_tig1c_Et[]={0., 0.991279, 0.995014, 1.00048, 1.00033, 1.00026};
const float SF_e22_tig1c_Et_toterror[]={1, 0.00827327, 0.00547239, 0.00526172, 0.0051088, 0.00503921};


  copyToVector(SF_e22_loo1c_eta, 20, efficienciesRel17Morionde22_mediumLoosePP);
  copyToVector(SF_e22_loo1c_eta_toterror, 20, uncertaintiesRel17Morionde22_mediumLoosePP);
  copyToVector(SF_e22_loo1c_Et, 6, efficienciesRel17Morionde22_mediumLoosePPET);
  copyToVector(SF_e22_loo1c_Et_toterror, 6, uncertaintiesRel17Morionde22_mediumLoosePPET);

  copyToVector(SF_e22_med1c_eta, 20, efficienciesRel17Morionde22_mediumMediumPP);
  copyToVector(SF_e22_med1c_eta_toterror, 20, uncertaintiesRel17Morionde22_mediumMediumPP);
  copyToVector(SF_e22_med1c_Et, 6, efficienciesRel17Morionde22_mediumMediumPPET);
  copyToVector(SF_e22_med1c_Et_toterror, 6, uncertaintiesRel17Morionde22_mediumMediumPPET);

  copyToVector(SF_e22_tig1c_eta, 20, efficienciesRel17Morionde22_mediumTightPP);
  copyToVector(SF_e22_tig1c_eta_toterror, 20, uncertaintiesRel17Morionde22_mediumTightPP);
  copyToVector(SF_e22_tig1c_Et, 6, efficienciesRel17Morionde22_mediumTightPPET);
  copyToVector(SF_e22_tig1c_Et_toterror, 6, uncertaintiesRel17Morionde22_mediumTightPPET);

  copyToVector(mcEff_e22_loo1c_eta, 20, MCefficienciesRel17Morionde22_mediumLoosePP);
  copyToVector(mcEff_e22_loo1c_Et, 6, MCefficienciesRel17Morionde22_mediumLoosePPET);
  copyToVector(mcEff_e22_med1c_eta, 20, MCefficienciesRel17Morionde22_mediumMediumPP);
  copyToVector(mcEff_e22_med1c_Et, 6, MCefficienciesRel17Morionde22_mediumMediumPPET);
  copyToVector(mcEff_e22_tig1c_eta, 20, MCefficienciesRel17Morionde22_mediumTightPP);
  copyToVector(mcEff_e22_tig1c_Et, 6, MCefficienciesRel17Morionde22_mediumTightPPET);

  // e22vh_medium1 L-M
const float  mcEff_e22vh_loo1c_eta[20]={0.753385,0.878486,0.89769,0.906778,0.726772,0.932791,0.937339,0.92971,0.930485,0.915429,0.907636,0.930841,0.932583,0.938241,0.930442,0.742618,0.914729,0.899687,0.879789,0.757878};
const float  mcEff_e22vh_med1c_eta[20]={0.855469,0.957704,0.962118,0.957255,0.79774,0.975774,0.975337,0.971563,0.974653,0.973731,0.968879,0.976654,0.973704,0.977386,0.973054,0.811455,0.963224,0.96084,0.959001,0.858417};
const float  mcEff_e22vh_tig1c_eta[20]={0.873634,0.969314,0.971734,0.966338,0.803909,0.981871,0.980063,0.975153,0.977292,0.977283,0.972042,0.978912,0.976689,0.98102,0.977944,0.81777,0.97238,0.969662,0.969989,0.884153};

const float mcEff_e22vh_loo1c_Et[]={0., 0.841141, 0.916559, 0.969456, 0.989858, 1.01729};
const float mcEff_e22vh_med1c_Et[]={0., 0.917413, 0.970792, 0.991354, 0.997288, 1.00565};
const float mcEff_e22vh_tig1c_Et[]={0., 0.9261, 0.974482, 0.992761, 0.998153, 1.00453};

// SF 11c
const float  SF_e22vh_loo1c_eta[20]={0.927087,0.934726,0.974119,0.956884,0.975238,0.980472,0.987201,0.97768,1.00221,0.91543,0.950363,0.996395,0.998819,0.987623,0.979209,1.02748,0.956823,0.964382,0.940593,0.923664};
const float SF_e22vh_loo1c_eta_toterror[20]={0.0145067, 0.00322803, 0.00300779, 0.00775627, 0.00862159, 
0.00191749, 0.00135194, 0.00186679, 0.001231, 0.00377874, 0.00299133, 0.00114196, 0.00182424, 
0.0013885, 0.00201595, 0.00657098, 0.00553785, 0.00313311, 0.00284801, 0.014375};
 const float SF_e22vh_med1c_eta[]={0.94445, 0.96441, 0.966463, 0.972076, 0.966836, 0.98881, 0.992735, 
0.984184, 1.00167, 0.910663, 0.945655, 0.994987, 1.00315, 0.992454, 0.98608, 1.01339, 0.972269, 
0.957904, 0.971932, 0.934367}; 
const float SF_e22vh_med1c_eta_toterror[]={0.0124695, 0.00182701, 0.00207682, 0.00424349, 0.0058064, 
0.00118148, 0.000888236, 0.00123393, 0.000644093, 0.00317546, 0.00220527, 0.000653257, 0.00100326, 
0.000829406, 0.00128047, 0.00488073, 0.00246987, 0.00214208, 0.00168504, 0.012745}; 
const float SF_e22vh_tig1c_eta[]={0.943578, 0.965082, 0.971145, 0.975588, 0.971474, 0.994236, 0.994395, 
0.985146, 1.00519, 0.942796, 0.976236, 0.998936, 1.00446, 0.995353, 0.991931, 1.01868, 0.976645, 0.964693, 0.972793, 0.923387}; 
const float SF_e22vh_tig1c_eta_toterror[]={0.0124287, 0.00172295, 0.00196625, 0.00442023, 0.00585806, 
0.00107372, 0.000884212, 0.0012172, 0.00058651, 0.00276835, 0.00213339, 0.000617185, 0.000966887, 
0.000788727, 0.00119594, 0.00525708, 0.00246356, 0.00214231, 0.00163689, 0.0128776}; 

//
const float SF_e22vh_loo1c_Et[]={0., 0.961035, 0.977873, 0.988738, 0.997348, 1.00438};
const float SF_e22vh_loo1c_Et_toterror[]={1., 0.00925124, 0.00750312, 0.00589792, 0.00522709, 0.00505072};
const float SF_e22vh_med1c_Et[]={0., 0.991858, 0.993749, 0.996976, 0.99943, 1.00129};
const float SF_e22vh_med1c_Et_toterror[]={1., 0.00698902, 0.00587858, 0.00527041, 0.00506882, 0.00501655};
const float SF_e22vh_tig1c_Et[]={0., 0.990188, 0.994625, 0.997235, 0.999762, 1.00104};
const float SF_e22vh_tig1c_Et_toterror[]={1., 0.0070191, 0.00567937, 0.00522457, 0.00507731, 0.00501631};

  copyToVector(SF_e22vh_loo1c_eta, 20, efficienciesRel17Morionde22vh_medium1LoosePP);
  copyToVector(SF_e22vh_loo1c_eta_toterror, 20, uncertaintiesRel17Morionde22vh_medium1LoosePP);
  copyToVector(SF_e22vh_loo1c_Et, 6, efficienciesRel17Morionde22vh_medium1LoosePPET);
  copyToVector(SF_e22vh_loo1c_Et_toterror, 6, uncertaintiesRel17Morionde22vh_medium1LoosePPET);

  copyToVector(SF_e22vh_med1c_eta, 20, efficienciesRel17Morionde22vh_medium1MediumPP);
  copyToVector(SF_e22vh_med1c_eta_toterror, 20, uncertaintiesRel17Morionde22vh_medium1MediumPP);
  copyToVector(SF_e22vh_med1c_Et, 6, efficienciesRel17Morionde22vh_medium1MediumPPET);
  copyToVector(SF_e22vh_med1c_Et_toterror, 6, uncertaintiesRel17Morionde22vh_medium1MediumPPET);

  copyToVector(SF_e22vh_tig1c_eta, 20, efficienciesRel17Morionde22vh_medium1TightPP);
  copyToVector(SF_e22vh_tig1c_eta_toterror, 20, uncertaintiesRel17Morionde22vh_medium1TightPP);
  copyToVector(SF_e22vh_tig1c_Et, 6, efficienciesRel17Morionde22vh_medium1TightPPET);
  copyToVector(SF_e22vh_tig1c_Et_toterror, 6, uncertaintiesRel17Morionde22vh_medium1TightPPET);

  copyToVector(mcEff_e22vh_loo1c_eta, 20, MCefficienciesRel17Morionde22vh_medium1LoosePP);
  copyToVector(mcEff_e22vh_loo1c_Et, 6, MCefficienciesRel17Morionde22vh_medium1LoosePPET);
  copyToVector(mcEff_e22vh_med1c_eta, 20, MCefficienciesRel17Morionde22vh_medium1MediumPP);
  copyToVector(mcEff_e22vh_med1c_Et, 6, MCefficienciesRel17Morionde22vh_medium1MediumPPET);
  copyToVector(mcEff_e22vh_tig1c_eta, 20, MCefficienciesRel17Morionde22vh_medium1TightPP);
  copyToVector(mcEff_e22vh_tig1c_Et, 6, MCefficienciesRel17Morionde22vh_medium1TightPPET);



  efficienciesRel17MoriondForwardLoose.push_back(99.9);
  efficienciesRel17MoriondForwardLoose.push_back(100.6);
  efficienciesRel17MoriondForwardLoose.push_back(98.6);
  efficienciesRel17MoriondForwardLoose.push_back(101.3);
  efficienciesRel17MoriondForwardLoose.push_back(100.2);
  efficienciesRel17MoriondForwardLoose.push_back(100.6);
  efficienciesRel17MoriondForwardLoose.push_back(0.000);
  efficienciesRel17MoriondForwardLoose.push_back(101.3);
  efficienciesRel17MoriondForwardLoose.push_back(97.7);
  efficienciesRel17MoriondForwardLoose.push_back(85.6);
  uncertaintiesRel17MoriondForwardLoose.push_back(3.1);
  uncertaintiesRel17MoriondForwardLoose.push_back(3.0);
  uncertaintiesRel17MoriondForwardLoose.push_back(2.1);
  uncertaintiesRel17MoriondForwardLoose.push_back(2.8);
  uncertaintiesRel17MoriondForwardLoose.push_back(2.0);
  uncertaintiesRel17MoriondForwardLoose.push_back(1.8);
  uncertaintiesRel17MoriondForwardLoose.push_back(100.);
  uncertaintiesRel17MoriondForwardLoose.push_back(7.5);
  uncertaintiesRel17MoriondForwardLoose.push_back(4.5);
  uncertaintiesRel17MoriondForwardLoose.push_back(6.9);

  efficienciesRel17MoriondForwardTight.push_back(97.8);
  efficienciesRel17MoriondForwardTight.push_back(93.2);
  efficienciesRel17MoriondForwardTight.push_back(83.6);
  efficienciesRel17MoriondForwardTight.push_back(90.1);
  efficienciesRel17MoriondForwardTight.push_back(84.9);
  efficienciesRel17MoriondForwardTight.push_back(87.4);
  efficienciesRel17MoriondForwardTight.push_back(0.000);
  efficienciesRel17MoriondForwardTight.push_back(90.9);
  efficienciesRel17MoriondForwardTight.push_back(90.0);
  efficienciesRel17MoriondForwardTight.push_back(75.8);
  uncertaintiesRel17MoriondForwardTight.push_back(4.3);
  uncertaintiesRel17MoriondForwardTight.push_back(3.7);
  uncertaintiesRel17MoriondForwardTight.push_back(3.1);
  uncertaintiesRel17MoriondForwardTight.push_back(2.6);
  uncertaintiesRel17MoriondForwardTight.push_back(2.5);
  uncertaintiesRel17MoriondForwardTight.push_back(2.1);
  uncertaintiesRel17MoriondForwardTight.push_back(100.);
  uncertaintiesRel17MoriondForwardTight.push_back(5.2);
  uncertaintiesRel17MoriondForwardTight.push_back(3.7);
  uncertaintiesRel17MoriondForwardTight.push_back(5.2);

  // Frozen Shower corrections in FCAL
  efficienciesRel17MoriondFrozenShowersForwardLoose = efficienciesRel17MoriondForwardLoose;
  efficienciesRel17MoriondFrozenShowersForwardTight = efficienciesRel17MoriondForwardTight;

  efficienciesRel17MoriondFrozenShowersForwardLoose[7] *= 1.015;
  efficienciesRel17MoriondFrozenShowersForwardLoose[8] *= 1.008;
  efficienciesRel17MoriondFrozenShowersForwardLoose[9] *= 1.018;

  efficienciesRel17MoriondFrozenShowersForwardTight[7] *= 1.010;
  efficienciesRel17MoriondFrozenShowersForwardTight[8] *= 1.006;
  efficienciesRel17MoriondFrozenShowersForwardTight[9] *= 1.041;


  // 2011 data with rel. 17 and MC11a/b/c AFII ("Moriond SF")
const float sfLoosePP_Combined_Moriond_AFII_eta[] = {0.990782, 0.984431, 1.0473, 1.00123, 0.999666, 0.998526, 1.00158, 1.0011, 1.0422, 0.985262, 0.990121};
const float errsfLoosePP_Combined_Moriond_AFII_eta[] = {0.0206951, 0.0160733, 0.00920923, 0.00394199, 0.00433968, 0.00427718, 0.00434684, 0.00354063, 0.00832827, 0.0159877, 0.0185366};
const float sfLoosePP_Jpsi_Moriond_AFII_eta[] = {0.944222, 0.867498, 1.01712, 0.917386, 0.938785, 0.923568, 0.945621, 0.920033, 1.03348, 0.897483, 0.949521};
const float errsfLoosePP_Jpsi_Moriond_AFII_eta[] = {0.0528257, 0.0924515, 0.116728, 0.0549565, 0.0425954, 0.0456152, 0.0422733, 0.0611321, 0.129607, 0.0822144, 0.068279};
const float sfLoosePP_Combined_Moriond_AFII_pt[] = {0, 1.03524, 1.02731, 0.946599, 0.952656, 0.985613, 1.00222, 1.01042, 1.00727, 1.00107};
const float errsfLoosePP_Combined_Moriond_AFII_pt[] = {1., 0.0575071, 0.0554038, 0.0221221, 0.0191394, 0.0122292, 0.00584829, 0.00683113, 0.00533981, 0.0064188};

const float sfMedium_Combined_Moriond_AFII_eta[] = {0.993134, 0.973003, 1.0982, 0.997358, 1.0093, 1.01942, 1.01548, 0.998222, 1.0948, 0.97789, 0.992183};
const float errsfMedium_Combined_Moriond_AFII_eta[] = {0.0209136, 0.0138568, 0.00914549, 0.00418142, 0.00855772, 0.00459038, 0.0089923, 0.00384268, 0.00838153, 0.0146522, 0.0187527};
const float sfMedium_Jpsi_Moriond_AFII_eta[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
const float errsfMedium_Jpsi_Moriond_AFII_eta[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
const float sfMedium_Combined_Moriond_AFII_pt[] = {0., 1, 1, 0.94817, 0.938252, 0.976575, 0.997979, 1.01049, 1.01606, 1.01232};
const float errsfMedium_Combined_Moriond_AFII_pt[] = {1.00005, 1.00005, 1.00005, 0.0240558, 0.0205825, 0.00782366, 0.00611601, 0.00752741, 0.00545777, 0.00649618};

const float sfMediumPP_Combined_Moriond_AFII_eta[] = {0.985436, 0.979718, 1.1627, 1.00117, 1.01279, 1.02311, 1.01838, 1.00231, 1.154, 0.982628, 0.983485};
const float errsfMediumPP_Combined_Moriond_AFII_eta[] = {0.0182925, 0.0136327, 0.00914549, 0.00414907, 0.00754674, 0.00472048, 0.00843378, 0.00387936, 0.00866025, 0.0150049, 0.0178148};
const float sfMediumPP_Jpsi_Moriond_AFII_eta[] = {0.939877, 0.877112, 1.11866, 0.906626, 0.953615, 0.947036, 0.95061, 0.904586, 1.10336, 0.856909, 0.918436};
const float errsfMediumPP_Jpsi_Moriond_AFII_eta[] = {0.0544204, 0.0932347, 0.133983, 0.055122, 0.0422624, 0.0459001, 0.0418465, 0.0611701, 0.141489, 0.0819855, 0.0715578};
const float sfMediumPP_Combined_Moriond_AFII_pt[] = {0., 1.06326, 1.01785, 0.935893, 0.932939, 0.973115, 0.996294, 1.00996, 1.01637, 1.0127};
const float errsfMediumPP_Combined_Moriond_AFII_pt[] = {1., 0.0561516, 0.0509324, 0.0246582, 0.0217461, 0.00823226, 0.0063997, 0.00788586, 0.00555109, 0.00656721};

const float sfTightPP_Combined_Moriond_AFII_eta[] = {1.00032, 0.997738, 1.2036, 1.0263, 1.01493, 1.03012, 1.01828, 1.02795, 1.1994, 1.00119, 0.999282};
const float errsfTightPP_Combined_Moriond_AFII_eta[] = {0.0180631, 0.013335, 0.0120847, 0.00472454, 0.00465211, 0.00587229, 0.00919978, 0.00446069, 0.00974115, 0.0128199, 0.0176403};
const float sfTightPP_Jpsi_Moriond_AFII_eta[] = {0.961117, 0.878841, 1.12571, 0.908417, 0.962517, 0.947508, 0.958219, 0.897417, 1.15715, 0.859185, 0.943734};
const float errsfTightPP_Jpsi_Moriond_AFII_eta[] = {0.058636, 0.0955017, 0.14339, 0.0581671, 0.0453877, 0.0501971, 0.0461154, 0.0635471, 0.151759, 0.0846577, 0.0746562};
const float sfTightPP_Combined_Moriond_AFII_pt[] = {0., 1.05978, 1.02195, 0.944854, 0.938838, 0.974907, 0.995292, 1.0075, 1.01237, 1.01166};
const float errsfTightPP_Combined_Moriond_AFII_pt[] = {1., 0.0632592, 0.0528317, 0.0260368, 0.0231078, 0.00852419, 0.00717741, 0.00797035, 0.00573223, 0.00678183};


  // reco+TQ are by choice identical to the "Cern council" results!
  copyToVector(Sf_RecoTrkQ_Eta, 11, efficienciesRel17MoriondAFIIRecoTrkQual, 1.);
  copyToVector(Sf_RecoTrkQ_Eta_err, 11, uncertaintiesRel17MoriondAFIIRecoTrkQual, 1.);
 
 // Identification eta for probes between 15 and 50 GeV
  copyToVector(sfLoosePP_Combined_Moriond_AFII_eta, 11, efficienciesRel17MoriondAFIILoosePP1550);
  copyToVector(errsfLoosePP_Combined_Moriond_AFII_eta, 11, uncertaintiesRel17MoriondAFIILoosePP1550);
  copyToVector(sfMedium_Combined_Moriond_AFII_eta, 11, efficienciesRel17MoriondAFIIMedium1550);
  copyToVector(errsfMedium_Combined_Moriond_AFII_eta, 11, uncertaintiesRel17MoriondAFIIMedium1550);
  copyToVector(sfMediumPP_Combined_Moriond_AFII_eta, 11, efficienciesRel17MoriondAFIIMediumPP1550);
  copyToVector(errsfMediumPP_Combined_Moriond_AFII_eta, 11, uncertaintiesRel17MoriondAFIIMediumPP1550);
  copyToVector(sfTightPP_Combined_Moriond_AFII_eta, 11, efficienciesRel17MoriondAFIITightPP1550);
  copyToVector(errsfTightPP_Combined_Moriond_AFII_eta, 11, uncertaintiesRel17MoriondAFIITightPP1550);
  //Identification eta for low ET probes
  copyToVector(sfLoosePP_Jpsi_Moriond_AFII_eta, 11, efficienciesRel17MoriondAFIILoosePP415);
  copyToVector(errsfLoosePP_Jpsi_Moriond_AFII_eta, 11, uncertaintiesRel17MoriondAFIILoosePP415);
  copyToVector(sfMedium_Jpsi_Moriond_AFII_eta, 11, efficienciesRel17MoriondAFIIMedium415);
  copyToVector(errsfMedium_Jpsi_Moriond_AFII_eta, 11, uncertaintiesRel17MoriondAFIIMedium415);
  copyToVector(sfMediumPP_Jpsi_Moriond_AFII_eta, 11, efficienciesRel17MoriondAFIIMediumPP415);
  copyToVector(errsfMediumPP_Jpsi_Moriond_AFII_eta, 11, uncertaintiesRel17MoriondAFIIMediumPP415);
  copyToVector(sfTightPP_Jpsi_Moriond_AFII_eta, 11, efficienciesRel17MoriondAFIITightPP415);
  copyToVector(errsfTightPP_Jpsi_Moriond_AFII_eta, 11, uncertaintiesRel17MoriondAFIITightPP415);
  // ET correction
  copyToVector(sfLoosePP_Combined_Moriond_AFII_pt, 10, ETCorrectionsRel17MoriondAFIILoosePP);
  copyToVector(errsfLoosePP_Combined_Moriond_AFII_pt, 10, uncertaintiesETCorrectionsRel17MoriondAFIILoosePP);
  copyToVector(sfMedium_Combined_Moriond_AFII_pt, 10, ETCorrectionsRel17MoriondAFIIMedium);
  copyToVector(errsfMedium_Combined_Moriond_AFII_pt, 10, uncertaintiesETCorrectionsRel17MoriondAFIIMedium);
  copyToVector(sfMediumPP_Combined_Moriond_AFII_pt, 10, ETCorrectionsRel17MoriondAFIIMediumPP);
  copyToVector(errsfMediumPP_Combined_Moriond_AFII_pt, 10, uncertaintiesETCorrectionsRel17MoriondAFIIMediumPP);
  copyToVector(sfTightPP_Combined_Moriond_AFII_pt, 10, ETCorrectionsRel17MoriondAFIITightPP);
  copyToVector(errsfTightPP_Combined_Moriond_AFII_pt, 10, uncertaintiesETCorrectionsRel17MoriondAFIITightPP);
  
 // Trigger efficiencies 
 // e20_medium B-J

const float mcEff_e20_loo1AF_eta[]={0.837784, 0.958348, 0.954186, 0.944921, 0.834737, 0.971725, 
0.97076, 0.965811, 0.953751, 0.953137, 0.93615, 0.944741, 0.959838, 0.964714, 0.963277, 0.849676, 
0.951548, 0.954741, 0.958616, 0.846399}; 
 const float mcEff_e20_med1AF_eta[]={0.898123, 0.979152, 0.982123, 0.976245, 0.908137, 0.989558, 
0.991127, 0.993397, 0.984975, 0.988081, 0.976318, 0.979741, 0.987332, 0.986481, 0.981413, 0.916848, 
0.981398, 0.981793, 0.979101, 0.902182}; 
const float mcEff_e20_tig1AF_eta[]={0.91683, 0.986025, 0.98714, 0.981276, 0.912122, 0.992078, 
0.993612, 0.995387, 0.986526, 0.991699, 0.980298, 0.981418, 0.989378, 0.989054, 0.984199, 0.92104, 
0.986463, 0.986933, 0.985685, 0.921703}; 

const float mcEff_e20_loo1AF_Et[]={0.875913, 0.914573, 0.951531, 0.984388, 0.996399, 1.00911};
const float mcEff_e20_med1AF_Et[]={0.93754, 0.967349, 0.987042, 0.99702, 0.99998, 1.00245};
const float mcEff_e20_tig1AF_Et[]={0.945612, 0.971708, 0.988945, 0.99764, 1.00018, 1.00189};

const  float SF_e20_loo1AF_eta[]={1.01159, 0.974655, 0.982602, 0.968572, 1.01941, 0.983739, 0.987819, 
0.97726, 1.01088, 0.903296, 0.972191, 1.01517, 1.00338, 0.996529, 0.989118, 1.02178, 0.971944, 
0.981215, 0.981992, 1.01294}; 
const float SF_e20_loo1AF_eta_toterror[]={0.0127233, 0.00162322, 0.00192803, 0.00676287, 0.0112238, 
0.00184649, 0.00106451, 0.00145809, 0.00081248, 0.00363176, 0.00642653, 0.000852247, 0.00118298, 
0.00105471, 0.00147156, 0.00642835, 0.00484407, 0.00192131, 0.00152101, 0.0116767}; 
const  float SF_e20_med1AF_eta[]={0.982521, 0.97788, 0.98448, 0.983248, 0.964709, 0.998338, 0.990116, 
0.976901, 1.0038, 0.897178, 0.963795, 1.00591, 1.0032, 0.997999, 1.00497, 0.965928, 0.984674, 
0.979194, 0.985067, 0.988369}; 
const float SF_e20_med1AF_eta_toterror[]={0.0120236, 0.00133611, 0.00140082, 0.00237377, 0.0074995, 
0.000718399, 0.000621108, 0.000937908, 0.000408613, 0.00297396, 0.00374874, 0.000519895, 
0.000600476, 0.000664762, 0.000916474, 0.00516008, 0.00122018, 0.00152137, 0.00124171, 0.0112851}; 
const float SF_e20_tig1AF_eta[]={0.975036, 0.975921, 0.984533, 0.983307, 0.963394, 0.997896, 0.990237, 
0.976572, 1.00381, 0.918325, 0.971076, 1.00582, 1.00246, 0.997985, 1.00531, 0.965012, 0.985399, 
0.980433, 0.984076, 0.977792}; 
const float SF_e20_tig1AF_eta_toterror[]={0.0119353, 0.00134182, 0.00144698, 0.00244823, 0.00715443, 
0.000685729, 0.000590407, 0.000939461, 0.000393648, 0.0025467, 0.00279181, 0.000515491, 
0.000565587, 0.000603549, 0.000878349, 0.00568995, 0.00114848, 0.00145738, 0.00123273, 0.0113562}; 

const float SF_e20_loo1AF_Et[]={0.966609, 0.972459, 0.984701, 0.991333, 0.996887, 1.00341};
const float SF_e20_loo1AF_Et_toterror[]={0.0101678, 0.00776877, 0.00648891, 0.00567344, 0.00510576, 
0.00503267};
const float SF_e20_med1AF_Et[]={0.998711, 0.994343, 0.995275, 0.996295, 0.999002, 1.00149};
const float SF_e20_med1AF_Et_toterror[]={0.00721775, 0.00577872, 0.00530762, 0.00520579, 0.00503048, 0.00501212};
const float SF_e20_tig1AF_Et[]={0.998593, 0.991357, 0.994461, 0.996193, 0.99918, 1.0015};
const float SF_e20_tig1AF_Et_toterror[]={0.00736422, 0.00580596, 0.00526614, 0.00519258, 0.00503818, 0.005011};


  copyToVector(SF_e20_loo1AF_eta, 20, efficienciesRel17MoriondAFIIe20_mediumLoosePP);
  copyToVector(SF_e20_loo1AF_eta_toterror, 20, uncertaintiesRel17MoriondAFIIe20_mediumLoosePP);
  copyToVector(SF_e20_loo1AF_Et, 6, efficienciesRel17MoriondAFIIe20_mediumLoosePPET);
  copyToVector(SF_e20_loo1AF_Et_toterror, 6, uncertaintiesRel17MoriondAFIIe20_mediumLoosePPET);
  
  copyToVector(SF_e20_med1AF_eta, 20, efficienciesRel17MoriondAFIIe20_mediumMediumPP);
  copyToVector(SF_e20_med1AF_eta_toterror, 20, uncertaintiesRel17MoriondAFIIe20_mediumMediumPP);
  copyToVector(SF_e20_med1AF_Et, 6, efficienciesRel17MoriondAFIIe20_mediumMediumPPET);
  copyToVector(SF_e20_med1AF_Et_toterror, 6, uncertaintiesRel17MoriondAFIIe20_mediumMediumPPET);

  copyToVector(SF_e20_tig1AF_eta, 20, efficienciesRel17MoriondAFIIe20_mediumTightPP);
  copyToVector(SF_e20_tig1AF_eta_toterror, 20, uncertaintiesRel17MoriondAFIIe20_mediumTightPP);
  copyToVector(SF_e20_tig1AF_Et, 6, efficienciesRel17MoriondAFIIe20_mediumTightPPET);
  copyToVector(SF_e20_tig1AF_Et_toterror, 6, uncertaintiesRel17MoriondAFIIe20_mediumTightPPET);

  copyToVector(mcEff_e20_loo1AF_eta, 20, MCefficienciesRel17MoriondAFIIe20_mediumLoosePP);
  copyToVector(mcEff_e20_loo1AF_Et, 6, MCefficienciesRel17MoriondAFIIe20_mediumLoosePPET);
  copyToVector(mcEff_e20_med1AF_eta, 20, MCefficienciesRel17MoriondAFIIe20_mediumMediumPP);
  copyToVector(mcEff_e20_med1AF_Et, 6, MCefficienciesRel17MoriondAFIIe20_mediumMediumPPET);
  copyToVector(mcEff_e20_tig1AF_eta, 20, MCefficienciesRel17MoriondAFIIe20_mediumTightPP);
  copyToVector(mcEff_e20_tig1AF_Et, 6, MCefficienciesRel17MoriondAFIIe20_mediumTightPPET);


  // e22_medium K

const float mcEff_e22_loo1AF_eta[20] = {0.833747,0.956846,0.951756,0.943097,0.809434,0.969812,0.970325,0.965237,0.953058,0.951805,0.919851,0.927237,0.942461,0.946344,0.944457,0.825211,0.94888,0.953134,0.95764,0.843531};
const float mcEff_e22_med1AF_eta[20] = {0.90366, 0.97809, 0.98258, 0.97693, 0.88845, 0.98838, 0.99082, 0.99282, 0.98426, 0.98784, 0.96275, 0.96594, 0.97370, 0.97213, 0.96527, 0.89546, 0.98029, 0.98190, 0.97815, 0.89741};
const float mcEff_e22_tig1AF_eta[20] = {0.93069, 0.98667, 0.98724, 0.98215, 0.89277, 0.99071, 0.99317, 0.99507, 0.98575, 0.99166, 0.96768, 0.96734, 0.97543, 0.97595, 0.96925, 0.90150, 0.98575, 0.98697, 0.98643, 0.92243};
const float mcEff_e22_loo1AF_Et[]={0., 0.885186, 0.939948, 0.981269, 0.996204, 1.01011};
const float mcEff_e22_med1AF_Et[]={0., 0.940987, 0.979161, 0.994945, 0.999995, 1.00335};
const float mcEff_e22_tig1AF_Et[]={0., 0.946754, 0.981608, 0.995751, 1.00006, 1.00283};

const float  SF_e22_loo1AF_eta[20]={1.0536,0.97881,0.983469,0.969107,1.03897,0.983133,0.988574,0.977489,1.01089,0.904715,0.986548,1.03174,1.02249,1.01525,1.00702,1.03936,0.976207,0.978133,0.98838,1.03178};
const float SF_e22_loo1AF_eta_toterror[20]={0.0150737, 0.00273343, 0.00338387, 0.00654829, 0.0147355, 
0.00229121, 0.00153903, 0.002332, 0.00129245, 0.00545286, 0.00650272, 0.00137271, 0.00199354, 
0.0015861, 0.00240757, 0.00941302, 0.00430818, 0.00362917, 0.00235633, 0.0140102};
const float SF_e22_med1AF_eta[]={1.01241, 0.981162, 0.983998, 0.98209, 0.977993, 0.998123, 0.990079, 
0.977828, 1.003, 0.898257, 0.975397, 1.02197, 1.02169, 1.01657, 1.01914, 0.978953, 0.986046, 
0.973089, 0.989486, 1.00536}; 
const float SF_e22_med1AF_eta_toterror[]={0.0135984, 0.0022501, 0.00247966, 0.00378395, 0.010849, 
0.00121093, 0.0010623, 0.00158136, 0.000705665, 0.00462989, 0.00497647, 0.000843935, 0.00108361, 0.00103298, 
0.0016114, 0.00881784, 0.00192929, 0.00276556, 0.00189862, 0.0131987}; 
const float SF_e22_tig1AF_eta[]={1.00792, 0.978083, 0.98605, 0.982132, 0.979026, 0.997125, 0.990273, 
0.977802, 1.00294, 0.916215, 0.987263, 1.02216, 1.02163, 1.0161, 1.0188, 0.979266, 0.986212, 
0.974763, 0.98686, 0.994493}; 
const float SF_e22_tig1AF_eta_toterror[]={0.0134268, 0.00216703, 0.00238488, 0.00398101, 0.0107513, 
0.00120725, 0.00110686, 0.00157572, 0.00067908, 0.00427834, 0.00403234, 0.000845761, 0.00101562, 
0.00102069, 0.00155189, 0.00832287, 0.00198942, 0.00271315, 0.00190448, 0.0134414}; 

const float SF_e22_loo1AF_Et[]={0., 0.979921, 0.985088, 0.992992, 0.996368, 1.00287};
const float SF_e22_loo1AF_Et_toterror[]={1., 0.00943126, 0.00686203, 0.00605053, 0.00516196, 0.00506721};
const float SF_e22_med1AF_Et[]={0., 0.998053, 0.994937, 0.997683, 0.99907, 1.00108};
const float SF_e22_med1AF_Et_toterror[]={1., 0.00732952, 0.00586389, 0.00554881, 0.0050874, 0.00502854};
const float SF_e22_tig1AF_Et[]={0., 0.99551, 0.993202, 0.998189, 0.999363, 1.001};
const float SF_e22_tig1AF_Et_toterror[]={1., 0.00750782, 0.00566467, 0.00551505, 0.00508131, 0.00502837};


  copyToVector(SF_e22_loo1AF_eta, 20, efficienciesRel17MoriondAFIIe22_mediumLoosePP);
  copyToVector(SF_e22_loo1AF_eta_toterror, 20, uncertaintiesRel17MoriondAFIIe22_mediumLoosePP);
  copyToVector(SF_e22_loo1AF_Et, 6, efficienciesRel17MoriondAFIIe22_mediumLoosePPET);
  copyToVector(SF_e22_loo1AF_Et_toterror, 6, uncertaintiesRel17MoriondAFIIe22_mediumLoosePPET);

  copyToVector(SF_e22_med1AF_eta, 20, efficienciesRel17MoriondAFIIe22_mediumMediumPP);
  copyToVector(SF_e22_med1AF_eta_toterror, 20, uncertaintiesRel17MoriondAFIIe22_mediumMediumPP);
  copyToVector(SF_e22_med1AF_Et, 6, efficienciesRel17MoriondAFIIe22_mediumMediumPPET);
  copyToVector(SF_e22_med1AF_Et_toterror, 6, uncertaintiesRel17MoriondAFIIe22_mediumMediumPPET);

  copyToVector(SF_e22_tig1AF_eta, 20, efficienciesRel17MoriondAFIIe22_mediumTightPP);
  copyToVector(SF_e22_tig1AF_eta_toterror, 20, uncertaintiesRel17MoriondAFIIe22_mediumTightPP);
  copyToVector(SF_e22_tig1AF_Et, 6, efficienciesRel17MoriondAFIIe22_mediumTightPPET);
  copyToVector(SF_e22_tig1AF_Et_toterror, 6, uncertaintiesRel17MoriondAFIIe22_mediumTightPPET);

  copyToVector(mcEff_e22_loo1AF_eta, 20, MCefficienciesRel17MoriondAFIIe22_mediumLoosePP);
  copyToVector(mcEff_e22_loo1AF_Et, 6, MCefficienciesRel17MoriondAFIIe22_mediumLoosePPET);
  copyToVector(mcEff_e22_med1AF_eta, 20, MCefficienciesRel17MoriondAFIIe22_mediumMediumPP);
  copyToVector(mcEff_e22_med1AF_Et, 6, MCefficienciesRel17MoriondAFIIe22_mediumMediumPPET);
  copyToVector(mcEff_e22_tig1AF_eta, 20, MCefficienciesRel17MoriondAFIIe22_mediumTightPP);
  copyToVector(mcEff_e22_tig1AF_Et, 6, MCefficienciesRel17MoriondAFIIe22_mediumTightPPET);

  // e22vh_medium1 L-M

const float mcEff_e22vh_loo1AF_eta[20]={0.687217,0.869343,0.889377,0.894433,0.707848,0.927745,0.930883,0.926332,0.923035,0.914771,0.904406,0.921625,0.929644,0.931911,0.922862,0.717363,0.902978,0.888726,0.868617,0.710998};
const float mcEff_e22vh_med1AF_eta[20]={0.863631,0.956778,0.953849,0.948926,0.864313,0.97087,0.976042,0.977363,0.973809,0.975869,0.971452,0.976045,0.978746,0.978596,0.966552,0.872116,0.956014,0.952087,0.95684,0.868906};
const float  mcEff_e22vh_tig1AF_eta[20]={0.886873,0.967241,0.963728,0.957384,0.871022,0.976454,0.980351,0.980797,0.976152,0.979879,0.97632,0.978193,0.981469,0.982378,0.971531,0.879474,0.964959,0.962242,0.967308,0.890388};
const float mcEff_e22vh_loo1AF_Et[]={0., 0.840015, 0.913787, 0.97277, 0.994508, 1.01668};
const float mcEff_e22vh_med1AF_Et[]={0., 0.915316, 0.968868, 0.993366, 1.00034, 1.00444};
const float mcEff_e22vh_tig1AF_Et[]={0., 0.922225, 0.972056, 0.994273, 1.0006, 1.00361};

 const float SF_e22vh_loo1AF_eta[]={1.01634, 0.944556, 0.983225, 0.97009, 1.00127, 0.985803, 0.994047, 
0.981245, 1.0103, 0.916105, 0.953837, 1.00636, 1.00198, 0.994332, 0.987251, 1.06357, 0.969275, 
0.976274, 0.952691, 0.98456}; 
const float SF_e22vh_loo1AF_eta_toterror[]={0.0135111, 0.00261142, 0.00272109, 0.00769375, 0.0157165, 
0.00170863, 0.00119628, 0.00157821, 0.00106624, 0.00361085, 0.00797615, 0.00112804, 0.00152361, 
0.00120121, 0.00197525, 0.0117977, 0.00563133, 0.00274614, 0.0024229, 0.0136814}; 
const float SF_e22vh_med1AF_eta[]={0.935521, 0.965343, 0.974842, 0.980607, 0.892346, 0.993804, 0.992018, 
0.978343, 1.00254, 0.908673, 0.943179, 0.995608, 0.997979, 0.991227, 0.992714, 0.942859, 0.979602, 
0.96671, 0.974127, 0.923082}; 
const float SF_e22vh_med1AF_eta_toterror[]={0.0117289, 0.00160082, 0.00194001, 0.00325921, 0.00863573, 
0.00103898, 0.000756795, 0.00106191, 0.000555511, 0.00270058, 0.00484516, 0.000592502, 0.000825087, 
0.000741658, 0.00113719, 0.00614673, 0.00193288, 0.00193402, 0.00147208, 0.0117276}; 
const float SF_e22vh_tig1AF_eta[]={0.929489, 0.967151, 0.979212, 0.98471, 0.896601, 0.999751, 0.994103, 
0.979477, 1.00637, 0.940301, 0.971972, 0.99967, 0.999566, 0.993977, 0.998478, 0.947168, 0.984156, 
0.972131, 0.975489, 0.916916}; 
const float SF_e22vh_tig1AF_eta_toterror[]={0.0117703, 0.00153752, 0.00180921, 0.00352061, 0.008609, 
0.000942758, 0.000719669, 0.00103639, 0.000486768, 0.0022468, 0.0031436, 0.000543901, 0.000767379, 
0.000686692, 0.00103902, 0.00604, 0.0019273, 0.00192993, 0.0014464, 0.0120508}; 

const float SF_e22vh_loo1AF_Et[]={0., 0.96301, 0.981539, 0.986071, 0.993393, 1.0057};
const float SF_e22vh_loo1AF_Et_toterror[]={1., 0.00907074, 0.00832292, 0.00640083, 0.00523584, 0.00513875};
const float SF_e22vh_med1AF_Et[]={0., 0.994123, 0.995714, 0.994947, 0.996375, 1.00249};
const float SF_e22vh_med1AF_Et_toterror[]={1., 0.00670693, 0.0061998, 0.0055103, 0.00512157, 0.00502586};
const float SF_e22vh_tig1AF_Et[]={0., 0.994298, 0.997056, 0.995666, 0.997267, 1.00191};
const float SF_e22vh_tig1AF_Et_toterror[]={1., 0.00682907, 0.00591816, 0.0054974, 0.00514654, 0.00501857};


  copyToVector(SF_e22vh_loo1AF_eta, 20, efficienciesRel17MoriondAFIIe22vh_medium1LoosePP);
  copyToVector(SF_e22vh_loo1AF_eta_toterror, 20, uncertaintiesRel17MoriondAFIIe22vh_medium1LoosePP);
  copyToVector(SF_e22vh_loo1AF_Et, 6, efficienciesRel17MoriondAFIIe22vh_medium1LoosePPET);
  copyToVector(SF_e22vh_loo1AF_Et_toterror, 6, uncertaintiesRel17MoriondAFIIe22vh_medium1LoosePPET);

  copyToVector(SF_e22vh_med1AF_eta, 20, efficienciesRel17MoriondAFIIe22vh_medium1MediumPP);
  copyToVector(SF_e22vh_med1AF_eta_toterror, 20, uncertaintiesRel17MoriondAFIIe22vh_medium1MediumPP);
  copyToVector(SF_e22vh_med1AF_Et, 6, efficienciesRel17MoriondAFIIe22vh_medium1MediumPPET);
  copyToVector(SF_e22vh_med1AF_Et_toterror, 6, uncertaintiesRel17MoriondAFIIe22vh_medium1MediumPPET);

  copyToVector(SF_e22vh_tig1AF_eta, 20, efficienciesRel17MoriondAFIIe22vh_medium1TightPP);
  copyToVector(SF_e22vh_tig1AF_eta_toterror, 20, uncertaintiesRel17MoriondAFIIe22vh_medium1TightPP);
  copyToVector(SF_e22vh_tig1AF_Et, 6, efficienciesRel17MoriondAFIIe22vh_medium1TightPPET);
  copyToVector(SF_e22vh_tig1AF_Et_toterror, 6, uncertaintiesRel17MoriondAFIIe22vh_medium1TightPPET);

  copyToVector(mcEff_e22vh_loo1AF_eta, 20, MCefficienciesRel17MoriondAFIIe22vh_medium1LoosePP);
  copyToVector(mcEff_e22vh_loo1AF_Et, 6, MCefficienciesRel17MoriondAFIIe22vh_medium1LoosePPET);
  copyToVector(mcEff_e22vh_med1AF_eta, 20, MCefficienciesRel17MoriondAFIIe22vh_medium1MediumPP);
  copyToVector(mcEff_e22vh_med1AF_Et, 6, MCefficienciesRel17MoriondAFIIe22vh_medium1MediumPPET);
  copyToVector(mcEff_e22vh_tig1AF_eta, 20, MCefficienciesRel17MoriondAFIIe22vh_medium1TightPP);
  copyToVector(mcEff_e22vh_tig1AF_Et, 6, MCefficienciesRel17MoriondAFIIe22vh_medium1TightPPET);

}

std::pair<float,float> egammaSFclass::scaleFactor(float eta, float ET, int set, int range, int rel, bool etcorrection) {

   std::vector<float> * vectEff=0;
   std::vector<float> * vectUnc=0;
   std::vector<float> * vectEtaBinning=0;

   bool doAbsEta = false;

   if (rel==7) { //release 17 for 2011 data and AFII MC11a/b/c, "Moriond recommendations"
     // range is ignored here
     vectEtaBinning = &m_11Etabins;
     if (set == 0 || set == 2 || set == 3 || (set > 22 && set<27) || set>29) {
       std::cout << "egammaSFclass: ERROR : unknown correction set" << std::endl;
       return make_pair(-1.,-1.);
     }
     else if (set==4) {//Reco + track quality requirements
       // this has implicit ET dependence, so don't confuse the user
       etcorrection = false;
       vectEff = &efficienciesRel17MoriondAFIIRecoTrkQual;
       vectUnc = &uncertaintiesRel17MoriondAFIIRecoTrkQual;
       if (ET<15000.) {
	 float eff=1.,unc=0.05;
	 if (fabs(eta)<1.37) {
	   eff=1.;unc=0.02;
	 }
	 return make_pair(eff,unc);
       }
     }
     else if (set==1) {//Medium
       if (ET>=15000.) {
	 vectEff = &efficienciesRel17MoriondAFIIMedium1550;
	 vectUnc = &uncertaintiesRel17MoriondAFIIMedium1550;
       } else {
	 vectEff = &efficienciesRel17MoriondAFIIMedium415;
	 vectUnc = &uncertaintiesRel17MoriondAFIIMedium415;
       }
     }
     else if (set==5) {//Loose++
       if (ET>=15000.) {
	 vectEff = &efficienciesRel17MoriondAFIILoosePP1550;
	 vectUnc = &uncertaintiesRel17MoriondAFIILoosePP1550;
       } else {
	 vectEff = &efficienciesRel17MoriondAFIILoosePP415;
	 vectUnc = &uncertaintiesRel17MoriondAFIILoosePP415;
       }
     }
     else if (set==6) {//Medium++
       if (ET>=15000.) {
	 vectEff = &efficienciesRel17MoriondAFIIMediumPP1550;
	 vectUnc = &uncertaintiesRel17MoriondAFIIMediumPP1550;
       } else {
	 vectEff = &efficienciesRel17MoriondAFIIMediumPP415;
	 vectUnc = &uncertaintiesRel17MoriondAFIIMediumPP415;
       }
     }
     else if (set==7) {//Tight++
       if (ET>=15000.) {
	 vectEff = &efficienciesRel17MoriondAFIITightPP1550;
	 vectUnc = &uncertaintiesRel17MoriondAFIITightPP1550;
       } else {
	 vectEff = &efficienciesRel17MoriondAFIITightPP415;
	 vectUnc = &uncertaintiesRel17MoriondAFIITightPP415;
       }
     }
     else if (set==8) {//e20_medium trigger SF w.r.t Medium++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &efficienciesRel17MoriondAFIIe20_mediumMediumPP;
       vectUnc = &uncertaintiesRel17MoriondAFIIe20_mediumMediumPP;
     }
     else if (set==9) {//e20_medium MC trigger efficiency w.r.t Medium++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &MCefficienciesRel17MoriondAFIIe20_mediumMediumPP;
       vectUnc = 0; // no error
     }
     else if (set==10) {//e20_medium trigger SF w.r.t Tight++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &efficienciesRel17MoriondAFIIe20_mediumTightPP;
       vectUnc = &uncertaintiesRel17MoriondAFIIe20_mediumTightPP;
     }
     else if (set==11) {//e20_medium MC trigger efficiency w.r.t Tight++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &MCefficienciesRel17MoriondAFIIe20_mediumTightPP;
       vectUnc = 0; // no error
     }
     else if (set==12) {//e22_medium trigger SF w.r.t Medium++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &efficienciesRel17MoriondAFIIe22_mediumMediumPP;
       vectUnc = &uncertaintiesRel17MoriondAFIIe22_mediumMediumPP;
     }
     else if (set==13) {//e22_medium MC trigger efficiency w.r.t Medium++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &MCefficienciesRel17MoriondAFIIe22_mediumMediumPP;
       vectUnc = 0; // no error
     }
     else if (set==14) {//e22_medium trigger SF w.r.t Tight++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &efficienciesRel17MoriondAFIIe22_mediumTightPP;
       vectUnc = &uncertaintiesRel17MoriondAFIIe22_mediumTightPP;
     }
     else if (set==15) {//e22_medium MC trigger efficiency w.r.t Tight++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &MCefficienciesRel17MoriondAFIIe22_mediumTightPP;
       vectUnc = 0; // no error
     }
     else if (set==16) {//e22vh_medium1 trigger SF w.r.t Medium++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &efficienciesRel17MoriondAFIIe22vh_medium1MediumPP;
       vectUnc = &uncertaintiesRel17MoriondAFIIe22vh_medium1MediumPP;
     }
     else if (set==17) {//e22vh_medium1 MC trigger efficiency w.r.t Medium++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &MCefficienciesRel17MoriondAFIIe22vh_medium1MediumPP;
       vectUnc = 0; // no error
     }
     else if (set==18) {//e22vh_medium1 trigger SF w.r.t Tight++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &efficienciesRel17MoriondAFIIe22vh_medium1TightPP;
       vectUnc = &uncertaintiesRel17MoriondAFIIe22vh_medium1TightPP;
     }
     else if (set==19) {//e22vh_medium1 MC trigger efficiency w.r.t Tight++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &MCefficienciesRel17MoriondAFIIe22vh_medium1TightPP;
       vectUnc = 0; // no error
     }
     else if (set==20) {//e20_medium MC trigger efficiency w.r.t Loose++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &MCefficienciesRel17MoriondAFIIe20_mediumLoosePP;
       vectUnc = 0; // no error
     }
     else if (set==21) {//e22_medium MC trigger efficiency w.r.t Loose++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &MCefficienciesRel17MoriondAFIIe22_mediumLoosePP;
       vectUnc = 0; // no error
     }
     else if (set==22) {//e22vh_medium1 MC trigger efficiency w.r.t Loose++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &MCefficienciesRel17MoriondAFIIe22vh_medium1LoosePP;
       vectUnc = 0; // no error
     }
     else if (set==27) {//e20_medium trigger SF ATLF w.r.t Loose++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &efficienciesRel17MoriondAFIIe20_mediumLoosePP;
       vectUnc = &uncertaintiesRel17MoriondAFIIe20_mediumLoosePP;
     }
     else if (set==28) {//e22_medium trigger SF ATLF w.r.t Loose++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &efficienciesRel17MoriondAFIIe22_mediumLoosePP;
       vectUnc = &uncertaintiesRel17MoriondAFIIe22_mediumLoosePP;
     }
     else if (set==29) {//e22vh_medium1 trigger SF ATLF w.r.t Loose++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &efficienciesRel17MoriondAFIIe22vh_medium1LoosePP;
       vectUnc = &uncertaintiesRel17MoriondAFIIe22vh_medium1LoosePP;
     }

   }
   else if (rel==6) { //release 17 for 2011 data and G4 FullSim MC11a/b/c, "Moriond recommendations"
     // range is ignored here
     vectEtaBinning = &m_11Etabins;
     if (set == 0 || set == 2 || set == 3 || set > 29) {
       std::cout << "egammaSFclass: ERROR : only Reco+TrackQuality, Medium, Loose++, Medium++, Tight++, FwdLoose, FwdTight and 3 single electron triggers exist" << std::endl;
       return make_pair(-1.,-1.);
     }
     else if (set==4) {//Reco + track quality requirements
       // this has implicit ET dependence, so don't confuse the user
       etcorrection = false;
       vectEff = &efficienciesRel17MoriondRecoTrkQual;
       vectUnc = &uncertaintiesRel17MoriondRecoTrkQual;
       if (ET<15000.) {
	 float eff=1.,unc=0.05;
	 if (fabs(eta)<1.37) {
	   eff=1.;unc=0.02;
	 }
	 return make_pair(eff,unc);
       }
     }
     else if (set==1) {//Medium
       if (ET>=15000.) {
	 vectEff = &efficienciesRel17MoriondMedium1550;
	 vectUnc = &uncertaintiesRel17MoriondMedium1550;
       } else {
	 vectEff = &efficienciesRel17MoriondMedium415;
	 vectUnc = &uncertaintiesRel17MoriondMedium415;
       }
     }
     else if (set==5) {//Loose++
       if (ET>=15000.) {
	 vectEff = &efficienciesRel17MoriondLoosePP1550;
	 vectUnc = &uncertaintiesRel17MoriondLoosePP1550;
       } else {
	 vectEff = &efficienciesRel17MoriondLoosePP415;
	 vectUnc = &uncertaintiesRel17MoriondLoosePP415;
       }
     }
     else if (set==6) {//Medium++
       if (ET>=15000.) {
	 vectEff = &efficienciesRel17MoriondMediumPP1550;
	 vectUnc = &uncertaintiesRel17MoriondMediumPP1550;
       } else {
	 vectEff = &efficienciesRel17MoriondMediumPP415;
	 vectUnc = &uncertaintiesRel17MoriondMediumPP415;
       }
     }
     else if (set==7) {//Tight++
       if (ET>=15000.) {
	 vectEff = &efficienciesRel17MoriondTightPP1550;
	 vectUnc = &uncertaintiesRel17MoriondTightPP1550;
       } else {
	 vectEff = &efficienciesRel17MoriondTightPP415;
	 vectUnc = &uncertaintiesRel17MoriondTightPP415;
       }
     }
     else if (set==8) {//e20_medium trigger SF w.r.t Medium++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &efficienciesRel17Morionde20_mediumMediumPP;
       vectUnc = &uncertaintiesRel17Morionde20_mediumMediumPP;
     }
     else if (set==9) {//e20_medium MC trigger efficiency w.r.t Medium++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &MCefficienciesRel17Morionde20_mediumMediumPP;
       vectUnc = 0; // no error
     }
     else if (set==10) {//e20_medium trigger SF w.r.t Tight++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &efficienciesRel17Morionde20_mediumTightPP;
       vectUnc = &uncertaintiesRel17Morionde20_mediumTightPP;
     }
     else if (set==11) {//e20_medium MC trigger efficiency w.r.t Tight++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &MCefficienciesRel17Morionde20_mediumTightPP;
       vectUnc = 0; // no error
     }
     else if (set==12) {//e22_medium trigger SF w.r.t Medium++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &efficienciesRel17Morionde22_mediumMediumPP;
       vectUnc = &uncertaintiesRel17Morionde22_mediumMediumPP;
     }
     else if (set==13) {//e22_medium MC trigger efficiency w.r.t Medium++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &MCefficienciesRel17Morionde22_mediumMediumPP;
       vectUnc = 0; // no error
     }
     else if (set==14) {//e22_medium trigger SF w.r.t Tight++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &efficienciesRel17Morionde22_mediumTightPP;
       vectUnc = &uncertaintiesRel17Morionde22_mediumTightPP;
     }
     else if (set==15) {//e22_medium MC trigger efficiency w.r.t Tight++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &MCefficienciesRel17Morionde22_mediumTightPP;
       vectUnc = 0; // no error
     }
     else if (set==16) {//e22vh_medium1 trigger SF w.r.t Medium++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &efficienciesRel17Morionde22vh_medium1MediumPP;
       vectUnc = &uncertaintiesRel17Morionde22vh_medium1MediumPP;
     }
     else if (set==17) {//e22vh_medium1 MC trigger efficiency w.r.t Medium++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &MCefficienciesRel17Morionde22vh_medium1MediumPP;
       vectUnc = 0; // no error
     }
     else if (set==18) {//e22vh_medium1 trigger SF w.r.t Tight++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &efficienciesRel17Morionde22vh_medium1TightPP;
       vectUnc = &uncertaintiesRel17Morionde22vh_medium1TightPP;
     }
     else if (set==19) {//e22vh_medium1 MC trigger efficiency w.r.t Tight++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &MCefficienciesRel17Morionde22vh_medium1TightPP;
       vectUnc = 0; // no error
     }
     else if (set==20) {//e20_medium MC trigger efficiency w.r.t Loose++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &MCefficienciesRel17Morionde20_mediumLoosePP;
       vectUnc = 0; // no error
     }
     else if (set==21) {//e22_medium MC trigger efficiency w.r.t Loose++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &MCefficienciesRel17Morionde22_mediumLoosePP;
       vectUnc = 0; // no error
     }
     else if (set==22) {//e22vh_medium1 MC trigger efficiency w.r.t Loose++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &MCefficienciesRel17Morionde22vh_medium1LoosePP;
       vectUnc = 0; // no error
     }
     else if (set==23) {//ForwardLoose electrons
       doAbsEta = true;
       vectEtaBinning = &m_FwdEtabins;
       vectEff = &efficienciesRel17MoriondForwardLoose;
       vectUnc = &uncertaintiesRel17MoriondForwardLoose;
     }
     else if (set==24) {//ForwardTight electrons
       doAbsEta = true;
       vectEtaBinning = &m_FwdEtabins;
       vectEff = &efficienciesRel17MoriondForwardTight;
       vectUnc = &uncertaintiesRel17MoriondForwardTight;
     }
     else if (set==25) {//ForwardLoose electrons FrozenShowers
       doAbsEta = true;
       vectEtaBinning = &m_FwdEtabins;
       vectEff = &efficienciesRel17MoriondFrozenShowersForwardLoose;
       vectUnc = &uncertaintiesRel17MoriondForwardLoose;
     }
     else if (set==26) {//ForwardTight electrons FrozenShowers
       doAbsEta = true;
       vectEtaBinning = &m_FwdEtabins;
       vectEff = &efficienciesRel17MoriondFrozenShowersForwardTight;
       vectUnc = &uncertaintiesRel17MoriondForwardTight;
     }
     else if (set==27) {//e20_medium trigger SF MC11C w.r.t Loose++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &efficienciesRel17Morionde20_mediumLoosePP;
       vectUnc = &uncertaintiesRel17Morionde20_mediumLoosePP;
     }
     else if (set==28) {//e22_medium trigger SF MC11C w.r.t Loose++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &efficienciesRel17Morionde22_mediumLoosePP;
       vectUnc = &uncertaintiesRel17Morionde22_mediumLoosePP;
     }
     else if (set==29) {//e22vh_medium trigger SF MC11C w.r.t Loose++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &efficienciesRel17Morionde22vh_medium1LoosePP;
       vectUnc = &uncertaintiesRel17Morionde22vh_medium1LoosePP;
     }

   }
   else if (rel==5) { //release 17 for 2011 data and MC11a, "CERN council recommendations"
     // range is ignored here
     vectEtaBinning = &m_11Etabins;
     if (set < 4 || set > 22) {
       std::cout << "egammaSFclass: ERROR : only Reco+TrackQuality, IsEM++ menu, and 3 single electron triggers exist" << std::endl;
       return make_pair(-1.,-1.);
     }
     else if (set==4) {//Reco + track quality requirements
       // this has implicit ET dependence, so don't confuse the user
       etcorrection = false;
       vectEff = &efficienciesRel17CCRecoTrkQual;
       vectUnc = &uncertaintiesRel17CCRecoTrkQual;
       if (ET<15000.) {
	 float eff=1.,unc=0.05;
	 if (fabs(eta)<1.37) {
	   eff=1.;unc=0.02;
	 }
	 return make_pair(eff,unc);
       }
     }
     else if (set==5) {//Loose++
       if (ET>=15000.) {
	 vectEff = &efficienciesRel17CCLoosePP1550;
	 vectUnc = &uncertaintiesRel17CCLoosePP1550;
       } else {
	 vectEff = &efficienciesRel17CCLoosePP415;
	 vectUnc = &uncertaintiesRel17CCLoosePP415;
       }
     }
     else if (set==6) {//Medium++
       if (ET>=15000.) {
	 vectEff = &efficienciesRel17CCMediumPP1550;
	 vectUnc = &uncertaintiesRel17CCMediumPP1550;
       } else {
	 vectEff = &efficienciesRel17CCMediumPP415;
	 vectUnc = &uncertaintiesRel17CCMediumPP415;
       }
     }
     else if (set==7) {//Tight++
       if (ET>=15000.) {
	 vectEff = &efficienciesRel17CCTightPP1550;
	 vectUnc = &uncertaintiesRel17CCTightPP1550;
       } else {
	 vectEff = &efficienciesRel17CCTightPP415;
	 vectUnc = &uncertaintiesRel17CCTightPP415;
       }
     }
     else if (set==8) {//e20_medium trigger SF w.r.t Medium++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &efficienciesRel17CCe20_mediumMediumPP;
       vectUnc = &uncertaintiesRel17CCe20_mediumMediumPP;
     }
     else if (set==9) {//e20_medium MC trigger efficiency w.r.t Medium++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &MCefficienciesRel17CCe20_mediumMediumPP;
       vectUnc = 0; // no error
     }
     else if (set==10) {//e20_medium trigger SF w.r.t Tight++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &efficienciesRel17CCe20_mediumTightPP;
       vectUnc = &uncertaintiesRel17CCe20_mediumTightPP;
     }
     else if (set==11) {//e20_medium MC trigger efficiency w.r.t Tight++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &MCefficienciesRel17CCe20_mediumTightPP;
       vectUnc = 0; // no error
     }
     else if (set==12) {//e22_medium trigger SF w.r.t Medium++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &efficienciesRel17CCe22_mediumMediumPP;
       vectUnc = &uncertaintiesRel17CCe22_mediumMediumPP;
     }
     else if (set==13) {//e22_medium MC trigger efficiency w.r.t Medium++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &MCefficienciesRel17CCe22_mediumMediumPP;
       vectUnc = 0; // no error
     }
     else if (set==14) {//e22_medium trigger SF w.r.t Tight++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &efficienciesRel17CCe22_mediumTightPP;
       vectUnc = &uncertaintiesRel17CCe22_mediumTightPP;
     }
     else if (set==15) {//e22_medium MC trigger efficiency w.r.t Tight++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &MCefficienciesRel17CCe22_mediumTightPP;
       vectUnc = 0; // no error
     }
     else if (set==16) {//e22vh_medium1 trigger SF w.r.t Medium++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &efficienciesRel17CCe22vh_medium1MediumPP;
       vectUnc = &uncertaintiesRel17CCe22vh_medium1MediumPP;
     }
     else if (set==17) {//e22vh_medium1 MC trigger efficiency w.r.t Medium++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &MCefficienciesRel17CCe22vh_medium1MediumPP;
       vectUnc = 0; // no error
     }
     else if (set==18) {//e22vh_medium1 trigger SF w.r.t Tight++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &efficienciesRel17CCe22vh_medium1TightPP;
       vectUnc = &uncertaintiesRel17CCe22vh_medium1TightPP;
     }
     else if (set==19) {//e22vh_medium1 MC trigger efficiency w.r.t Tight++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &MCefficienciesRel17CCe22vh_medium1TightPP;
       vectUnc = 0; // no error
     }
     else if (set==20) {//e20_medium MC trigger efficiency w.r.t Loose++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &MCefficienciesRel17CCe20_mediumLoosePP;
       vectUnc = 0; // no error
     }
     else if (set==21) {//e22_medium MC trigger efficiency w.r.t Loose++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &MCefficienciesRel17CCe22_mediumLoosePP;
       vectUnc = 0; // no error
     }
     else if (set==22) {//e22vh_medium1 MC trigger efficiency w.r.t Loose++ offline
       vectEtaBinning = &m_FineEtabins;
       vectEff = &MCefficienciesRel17CCe22vh_medium1LoosePP;
       vectUnc = 0; // no error
     }
   }
   else if (rel==4) { //release 16.6 numbers estimated from 2011 data, "EPS recommendations" including the low ET region
     vectEtaBinning = &m_FineEtabins;
     if (range==0) { //20-50 GeV region
       if (set==0 || set>4) {
	 std::cout << "egammaSFclass: ERROR : only Medium, Tight and trigger scale factors exist" << std::endl;
	 return make_pair(-1.,-1.);
       }
       else if (set==1) {//Medium
	 if (ET>=15000.) {
	   vectEff = &efficienciesRel166EPSMedium2050;
	   vectUnc = &uncertaintiesRel166EPSMedium2050;
	 }
	 else {
	   vectEtaBinning = &m_Etabins;
	   vectEff = &efficienciesRel166EPSMediumLowET;
	   vectUnc = &uncertaintiesRel166EPSMediumLowET;
	 }
       }
       else if (set==2) {//Tight
	 if (ET>=15000.) {
	   vectEff = &efficienciesRel166EPSTight2050;
	   vectUnc = &uncertaintiesRel166EPSTight2050;
	 }
	 else {
	   vectEtaBinning = &m_Etabins;
	   vectEff = &efficienciesRel166EPSTightLowET;
	   vectUnc = &uncertaintiesRel166EPSTightLowET;
	 }
       }
       else if (set==3) {//Trigger
	 vectEff = &efficienciesRel166EPSTrigger;
	 vectUnc = &uncertaintiesRel166EPSTrigger;
       }
       else if (set==4) {//Reco + track quality requirements
	 vectEff = &efficienciesRel166EPSRecoTrkQual;
	 vectUnc = &uncertaintiesRel166EPSRecoTrkQual;
	 if (ET<15000.) {
	   float eff=1.,unc=0.05;
	   if (fabs(eta)<1.37) {
	     eff=1.;unc=0.02;
	   }
	   return make_pair(eff,unc);
	 }
       }
     }//endif 20-50 GeV
     else {
	 std::cout << "egammaSFclass: ERROR : invalid range" << std::endl;
	 return make_pair(-1.,-1.);
     }
   } 
   else if (rel==3) { //release 16.6 numbers estimated from 2011 data, "EPS recommendations"
     vectEtaBinning = &m_FineEtabins;
     if (range==0) { //20-50 GeV region
       if (set==0 || set>4) {
	 std::cout << "egammaSFclass: ERROR : only Medium, Tight and trigger scale factors exist" << std::endl;
	 return make_pair(-1.,-1.);
       }
       else if (set==1) {//Medium
	 vectEff = &efficienciesRel166EPSMedium2050;
	 vectUnc = &uncertaintiesRel166EPSMedium2050;
       }
       else if (set==2) {//Tight
	 vectEff = &efficienciesRel166EPSTight2050;
	 vectUnc = &uncertaintiesRel166EPSTight2050;
       }
       else if (set==3) {//Trigger
	 vectEff = &efficienciesRel166EPSTrigger;
	 vectUnc = &uncertaintiesRel166EPSTrigger;
       }
       else if (set==4) {//Reco + track quality requirements
	 vectEff = &efficienciesRel166EPSRecoTrkQual;
	 vectUnc = &uncertaintiesRel166EPSRecoTrkQual;
       }
     }//endif 20-50 GeV
     else {
	 std::cout << "egammaSFclass: ERROR : invalid range" << std::endl;
	 return make_pair(-1.,-1.);
     }
   } 
   else if (rel==2) { //release 16.6 numbers estimated from 2010 data
     vectEtaBinning = &m_Etabins;
     if (range==0) { //20-50 GeV region
       if (set==0 || set>2) {//Loose
	 std::cout << "egammaSFclass: ERROR : only Medium and Tight scale factors exist" << std::endl;
	 return make_pair(-1.,-1.);
       }
       else if (set==1) {//Medium
	 vectEff = &efficienciesRel166Data2010Medium2050;
	 vectUnc = &uncertaintiesRel166Data2010Medium2050;
       }
       else if (set==2) {//Tight
	 vectEff = &efficienciesRel166Data2010Tight2050;
	 vectUnc = &uncertaintiesRel166Data2010Tight2050;
       }
     }//endif 20-50 GeV
     else if (range==1) { //30-50 GeV region
       if (set==0 || set>2) {//Loose
	 std::cout << "egammaSFclass: ERROR : only Medium and Tight scale factors exist" << std::endl;
	 return make_pair(-1.,-1.);
       }
       else if (set==1) {//Medium
	 vectEff = &efficienciesRel166Data2010Medium3050;
	 vectUnc = &uncertaintiesRel166Data2010Medium3050;
       }
       else if (set==2) {//Tight
	 vectEff = &efficienciesRel166Data2010Tight3050;
	 vectUnc = &uncertaintiesRel166Data2010Tight3050;
       }
     }//endif 30-50 GeV
     else {
	 std::cout << "egammaSFclass: ERROR : invalid range" << std::endl;
	 return make_pair(-1.,-1.);
     }
   } 
   else if (rel==1) { //release 16 numbers
     vectEtaBinning = &m_Etabins;
     if (range==0) { //20-50 GeV region
       if (set==0 || set>2) {//Loose
	 std::cout << "egammaSFclass: ERROR : only Medium and Tight scale factors exist" << std::endl;
	 return make_pair(-1.,-1.);
       }
       else if (set==1) {//Medium
	 vectEff = &efficienciesRel16Medium2050;
	 vectUnc = &uncertaintiesRel16Medium2050;
       }
       else if (set==2) {//Tight
	 vectEff = &efficienciesRel16Tight2050;
	 vectUnc = &uncertaintiesRel16Tight2050;
       }
     }//endif 20-50 GeV
     else if (range==1) { //30-50 GeV region
       if (set==0 || set>2) {//Loose
	 std::cout << "egammaSFclass: ERROR : only Medium and Tight scale factors exist" << std::endl;
	 return make_pair(-1.,-1.);
       }
       else if (set==1) {//Medium
	 vectEff = &efficienciesRel16Medium3050;
	 vectUnc = &uncertaintiesRel16Medium3050;
       }
       else if (set==2) {//Tight
	 vectEff = &efficienciesRel16Tight3050;
	 vectUnc = &uncertaintiesRel16Tight3050;
       }
     }//endif 30-50 GeV
     else {
	 std::cout << "egammaSFclass: ERROR : invalid range" << std::endl;
	 return make_pair(-1.,-1.);
     }
   }
   else { //release 15 numbers
     vectEtaBinning = &m_Etabins;
     if (range==0) { //20-50 GeV region
       if (set==0) {//Loose
	 vectEff = &efficienciesRel15Loose2050;
	 vectUnc = &uncertaintiesRel15Loose2050;
       }
       else if (set==1) {//Medium
	 vectEff = &efficienciesRel15Medium2050;
	 vectUnc = &uncertaintiesRel15Medium2050;
       }
       else if (set==2) {//Tight
	 vectEff = &efficienciesRel15Tight2050;
	 vectUnc = &uncertaintiesRel15Tight2050;
       }
       else {
	 std::cout << "egammaSFclass: ERROR : invalid set of cuts" << std::endl;
	 return make_pair(-1.,-1.);
       }
     }//endif 20-50 GeV
     else if (range==1) { //30-50 GeV region
       if (set==0) {//Loose
	 vectEff = &efficienciesRel15Loose3050;
	 vectUnc = &uncertaintiesRel15Loose3050;
       }
       else if (set==1) {//Medium
	 vectEff = &efficienciesRel15Medium3050;
	 vectUnc = &uncertaintiesRel15Medium3050;
       }
       else if (set==2) {//Tight
	 vectEff = &efficienciesRel15Tight3050;
	 vectUnc = &uncertaintiesRel15Tight3050;
       }
       else {
	 std::cout << "egammaSFclass: ERROR : invalid set of cuts" << std::endl;
	 return make_pair(-1.,-1.);
       }
     }//endif 30-50 GeV
     else {
	 std::cout << "egammaSFclass: ERROR : invalid range" << std::endl;
	 return make_pair(-1.,-1.);
     }
   }//endif rel15

   //Choice of the eta bin
   int ietabin=-1;
   if (doAbsEta)
     eta = fabs(eta);

   while (ietabin<((int)vectEtaBinning->size()-1) && eta>=vectEtaBinning->at(ietabin+1)) ietabin++;
   if (ietabin<0 || ietabin>=((int)vectEtaBinning->size()-1)) {
     std::cout << "egammaSFclass: ERROR : given eta value outside range of existing scale factors" << std::endl;
     return make_pair(-1.,-1.);
   }


   float effvseta = vectEff->at(ietabin)/100.;
   float uncvseta = 0.;
   if (vectUnc)
     uncvseta = vectUnc->at(ietabin)/100.;

   float eff = effvseta;
   float unc = uncvseta;

   if (etcorrection) {
     std::pair<float,float> corr = etCorrection(ET, set, rel);
     if (corr.first <= 0 || eff <= 0)
       unc = 1.;
     else
       unc = eff*corr.first * sqrt( unc*unc/(eff*eff) + corr.second*corr.second/(corr.first*corr.first) );
     eff *= corr.first;
   }

   return make_pair(eff,unc);
}


//Returns the ET-correction factor (and uncertainty) to the scale factor for the correspond ET bin 
std::pair<float,float> egammaSFclass::etCorrection(float ET, int set, int rel) {
  //for backport of rel16 SF-ET-dependence to rel15
  if (rel==0) rel=1;  
  
  std::vector<float> * vectCorr=0;
  std::vector<float> * vectUncCorr=0;
  std::vector<float> * vectETBinning=0;
  vectETBinning = &m_ETbinsFullRange;
  
  if (rel==1) {
    vectETBinning = &m_ETbins;
    if (set==1) {
      vectCorr = &ETCorrectionsMediumRel16;
      vectUncCorr = &uncertaintiesETCorrectionsMediumRel16;
    }
    else if (rel==2) { // tight
      vectCorr = &ETCorrectionsTightRel16;
      vectUncCorr = &uncertaintiesETCorrectionsTightRel16;
    }
  }
  else if (rel==2) {
    vectETBinning = &m_ETbins;
    if (set==1) { //Medium
      vectCorr = &ETCorrectionsMediumRel166Data2010;
      vectUncCorr = &uncertaintiesETCorrectionsMediumRel166Data2010;
    }
    else if (set==2) {
      vectCorr = &ETCorrectionsTightRel166Data2010;
      vectUncCorr = &uncertaintiesETCorrectionsTightRel166Data2010;
    }
  } 
  else if (rel==3) {
    vectETBinning = &m_ETbins;
    if (set==1) { //Medium
      vectCorr = &ETCorrectionsMediumRel166EPS;
      vectUncCorr = &uncertaintiesETCorrectionsMediumRel166EPS;
    }
    else if (set==2) {
      vectCorr = &ETCorrectionsTightRel166EPS;
      vectUncCorr = &uncertaintiesETCorrectionsTightRel166EPS;
    }
  }
  else if (rel==4) {
    if (set==1) { //Medium
      vectCorr = &ETCorrectionsMediumRel166EPSFullRange;
      vectUncCorr = &uncertaintiesETCorrectionsMediumRel166EPSFullRange;
    }
    else if (set==2) {
      vectCorr = &ETCorrectionsTightRel166EPSFullRange;
      vectUncCorr = &uncertaintiesETCorrectionsTightRel166EPSFullRange;
    }
  }
  else if (rel==5) {// Rel17CC
    if (set==5) { // Loose++
      vectETBinning = &m_ETbinsFullRange;
      vectCorr = &ETCorrectionsRel17CCLoosePP;
      vectUncCorr = &uncertaintiesETCorrectionsRel17CCLoosePP;
    }
    else if (set==6) { // Medium++
      vectETBinning = &m_ETbinsFullRange;
      vectCorr = &ETCorrectionsRel17CCMediumPP;
      vectUncCorr = &uncertaintiesETCorrectionsRel17CCMediumPP;
    }
    else if (set==7) { // Tight++
      vectETBinning = &m_ETbinsFullRange;
      vectCorr = &ETCorrectionsRel17CCTightPP;
      vectUncCorr = &uncertaintiesETCorrectionsRel17CCTightPP;
    }
    else if (set==8) { // e20_medium trigger SF w.r.t Medium++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &efficienciesRel17CCe20_mediumMediumPPET;
      vectUncCorr = &uncertaintiesRel17CCe20_mediumMediumPPET;
    }
    else if (set==9) { // e20_medium trigger MC eff w.r.t Medium++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &MCefficienciesRel17CCe20_mediumMediumPPET;
      vectUncCorr = 0;
    }
    else if (set==10) { // e20_medium trigger SF w.r.t Tight++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &efficienciesRel17CCe20_mediumTightPPET;
      vectUncCorr = &uncertaintiesRel17CCe20_mediumTightPPET;
    }
    else if (set==11) { // e20_medium trigger MC eff w.r.t Tight++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &MCefficienciesRel17CCe20_mediumTightPPET;
      vectUncCorr = 0;
    }
    else if (set==12) { // e22_medium trigger SF w.r.t Medium++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &efficienciesRel17CCe22_mediumMediumPPET;
      vectUncCorr = &uncertaintiesRel17CCe22_mediumMediumPPET;
    }
    else if (set==13) { // e22_medium trigger MC eff w.r.t Medium++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &MCefficienciesRel17CCe22_mediumMediumPPET;
      vectUncCorr = 0;
    }
    else if (set==14) { // e22_medium trigger SF w.r.t Tight++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &efficienciesRel17CCe22_mediumTightPPET;
      vectUncCorr = &uncertaintiesRel17CCe22_mediumTightPPET;
    }
    else if (set==15) { // e22_medium trigger MC eff w.r.t Tight++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &MCefficienciesRel17CCe22_mediumTightPPET;
      vectUncCorr = 0;
    }
    else if (set==16) { // e22vh_medium1 trigger SF w.r.t Medium++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &efficienciesRel17CCe22vh_medium1MediumPPET;
      vectUncCorr = &uncertaintiesRel17CCe22vh_medium1MediumPPET;
    }
    else if (set==17) { // e22vh_medium1 trigger MC eff w.r.t Medium++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &MCefficienciesRel17CCe22vh_medium1MediumPPET;
      vectUncCorr = 0;
    }
    else if (set==18) { // e22vh_medium1 trigger SF w.r.t Tight++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &efficienciesRel17CCe22vh_medium1TightPPET;
      vectUncCorr = &uncertaintiesRel17CCe22vh_medium1TightPPET;
    }
    else if (set==19) { // e22vh_medium1 trigger MC eff w.r.t Tight++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &MCefficienciesRel17CCe22vh_medium1TightPPET;
      vectUncCorr = 0;
    }
    else if (set==20) { // e20_medium trigger MC eff w.r.t Loose++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &MCefficienciesRel17CCe20_mediumLoosePPET;
      vectUncCorr = 0;
    }
    else if (set==21) { // e22_medium trigger MC eff w.r.t Loose++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &MCefficienciesRel17CCe22_mediumLoosePPET;
      vectUncCorr = 0;
    }
    else if (set==22) { // e22vh_medium1 trigger MC eff w.r.t Loose++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &MCefficienciesRel17CCe22vh_medium1LoosePPET;
      vectUncCorr = 0;
    }
  }
  else if (rel==6) {// Rel17 Moriond G4FullSim
    if (set==1) { // Medium
      vectETBinning = &m_ETbinsFullRange;
      vectCorr = &ETCorrectionsRel17MoriondMedium;
      vectUncCorr = &uncertaintiesETCorrectionsRel17MoriondMedium;
    }
    else if (set==5) { // Loose++
      vectETBinning = &m_ETbinsFullRange;
      vectCorr = &ETCorrectionsRel17MoriondLoosePP;
      vectUncCorr = &uncertaintiesETCorrectionsRel17MoriondLoosePP;
    }
    else if (set==6) { // Medium++
      vectETBinning = &m_ETbinsFullRange;
      vectCorr = &ETCorrectionsRel17MoriondMediumPP;
      vectUncCorr = &uncertaintiesETCorrectionsRel17MoriondMediumPP;
    }
    else if (set==7) { // Tight++
      vectETBinning = &m_ETbinsFullRange;
      vectCorr = &ETCorrectionsRel17MoriondTightPP;
      vectUncCorr = &uncertaintiesETCorrectionsRel17MoriondTightPP;
    }
    else if (set==8) { // e20_medium trigger SF w.r.t Medium++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &efficienciesRel17Morionde20_mediumMediumPPET;
      vectUncCorr = &uncertaintiesRel17Morionde20_mediumMediumPPET;
    }
    else if (set==9) { // e20_medium trigger MC eff w.r.t Medium++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &MCefficienciesRel17Morionde20_mediumMediumPPET;
      vectUncCorr = 0;
    }
    else if (set==10) { // e20_medium trigger SF w.r.t Tight++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &efficienciesRel17Morionde20_mediumTightPPET;
      vectUncCorr = &uncertaintiesRel17Morionde20_mediumTightPPET;
    }
    else if (set==11) { // e20_medium trigger MC eff w.r.t Tight++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &MCefficienciesRel17Morionde20_mediumTightPPET;
      vectUncCorr = 0;
    }
    else if (set==12) { // e22_medium trigger SF w.r.t Medium++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &efficienciesRel17Morionde22_mediumMediumPPET;
      vectUncCorr = &uncertaintiesRel17Morionde22_mediumMediumPPET;
    }
    else if (set==13) { // e22_medium trigger MC eff w.r.t Medium++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &MCefficienciesRel17Morionde22_mediumMediumPPET;
      vectUncCorr = 0;
    }
    else if (set==14) { // e22_medium trigger SF w.r.t Tight++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &efficienciesRel17Morionde22_mediumTightPPET;
      vectUncCorr = &uncertaintiesRel17Morionde22_mediumTightPPET;
    }
    else if (set==15) { // e22_medium trigger MC eff w.r.t Tight++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &MCefficienciesRel17Morionde22_mediumTightPPET;
      vectUncCorr = 0;
    }
    else if (set==16) { // e22vh_medium1 trigger SF w.r.t Medium++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &efficienciesRel17Morionde22vh_medium1MediumPPET;
      vectUncCorr = &uncertaintiesRel17Morionde22vh_medium1MediumPPET;
    }
    else if (set==17) { // e22vh_medium1 trigger MC eff w.r.t Medium++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &MCefficienciesRel17Morionde22vh_medium1MediumPPET;
      vectUncCorr = 0;
    }
    else if (set==18) { // e22vh_medium1 trigger SF w.r.t Tight++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &efficienciesRel17Morionde22vh_medium1TightPPET;
      vectUncCorr = &uncertaintiesRel17Morionde22vh_medium1TightPPET;
    }
    else if (set==19) { // e22vh_medium1 trigger MC eff w.r.t Tight++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &MCefficienciesRel17Morionde22vh_medium1TightPPET;
      vectUncCorr = 0;
    }
    else if (set==20) { // e20_medium trigger MC eff w.r.t Loose++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &MCefficienciesRel17Morionde20_mediumLoosePPET;
      vectUncCorr = 0;
    }
    else if (set==21) { // e22_medium trigger MC eff w.r.t Loose++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &MCefficienciesRel17Morionde22_mediumLoosePPET;
      vectUncCorr = 0;
    }
    else if (set==22) { // e22vh_medium1 trigger MC eff w.r.t Loose++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &MCefficienciesRel17Morionde22vh_medium1LoosePPET;
      vectUncCorr = 0;
    }
    else if (set==27) { // e20_medium trigger SF w.r.t Loose++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &efficienciesRel17Morionde20_mediumLoosePPET;
      vectUncCorr = &uncertaintiesRel17Morionde20_mediumLoosePPET;
    }
    else if (set==28) { // e22_medium trigger SF w.r.t Loose++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &efficienciesRel17Morionde22_mediumLoosePPET;
      vectUncCorr = &uncertaintiesRel17Morionde22_mediumLoosePPET;
    }
    else if (set==29) { // e22vh_medium1 trigger SF w.r.t Loose++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &efficienciesRel17Morionde22vh_medium1LoosePPET;
      vectUncCorr = &uncertaintiesRel17Morionde22vh_medium1LoosePPET;
    }

    else if (set==23 || set == 24 || set == 25 || set == 26) { // Forward Loose+Tight: just make sure, it's not used below 20 GeV
      if (ET < 20000.) {
	std::cout << "egammaSFclass: ERROR : Out of Et range for forward electrons" << std::endl;
	return make_pair(-1.,-1.);
      } else {
	return make_pair(1.,0.);
      }
    }
  }
  else if (rel==7) {// Rel17 Moriond AFII
    if (set==1) { // Medium
      vectETBinning = &m_ETbinsFullRange;
      vectCorr = &ETCorrectionsRel17MoriondAFIIMedium;
      vectUncCorr = &uncertaintiesETCorrectionsRel17MoriondAFIIMedium;
    }
    else if (set==5) { // Loose++
      vectETBinning = &m_ETbinsFullRange;
      vectCorr = &ETCorrectionsRel17MoriondAFIILoosePP;
      vectUncCorr = &uncertaintiesETCorrectionsRel17MoriondAFIILoosePP;
    }
    else if (set==6) { // Medium++
      vectETBinning = &m_ETbinsFullRange;
      vectCorr = &ETCorrectionsRel17MoriondAFIIMediumPP;
      vectUncCorr = &uncertaintiesETCorrectionsRel17MoriondAFIIMediumPP;
    }
    else if (set==7) { // Tight++
      vectETBinning = &m_ETbinsFullRange;
      vectCorr = &ETCorrectionsRel17MoriondAFIITightPP;
      vectUncCorr = &uncertaintiesETCorrectionsRel17MoriondAFIITightPP;
    }
    else if (set==8) { // e20_medium trigger SF w.r.t Medium++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &efficienciesRel17MoriondAFIIe20_mediumMediumPPET;
      vectUncCorr = &uncertaintiesRel17MoriondAFIIe20_mediumMediumPPET;
    }
    else if (set==9) { // e20_medium trigger MC eff w.r.t Medium++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &MCefficienciesRel17MoriondAFIIe20_mediumMediumPPET;
      vectUncCorr = 0;
    }
    else if (set==10) { // e20_medium trigger SF w.r.t Tight++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &efficienciesRel17MoriondAFIIe20_mediumTightPPET;
      vectUncCorr = &uncertaintiesRel17MoriondAFIIe20_mediumTightPPET;
    }
    else if (set==11) { // e20_medium trigger MC eff w.r.t Tight++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &MCefficienciesRel17MoriondAFIIe20_mediumTightPPET;
      vectUncCorr = 0;
    }
    else if (set==12) { // e22_medium trigger SF w.r.t Medium++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &efficienciesRel17MoriondAFIIe22_mediumMediumPPET;
      vectUncCorr = &uncertaintiesRel17MoriondAFIIe22_mediumMediumPPET;
    }
    else if (set==13) { // e22_medium trigger MC eff w.r.t Medium++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &MCefficienciesRel17MoriondAFIIe22_mediumMediumPPET;
      vectUncCorr = 0;
    }
    else if (set==14) { // e22_medium trigger SF w.r.t Tight++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &efficienciesRel17MoriondAFIIe22_mediumTightPPET;
      vectUncCorr = &uncertaintiesRel17MoriondAFIIe22_mediumTightPPET;
    }
    else if (set==15) { // e22_medium trigger MC eff w.r.t Tight++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &MCefficienciesRel17MoriondAFIIe22_mediumTightPPET;
      vectUncCorr = 0;
    }
    else if (set==16) { // e22vh_medium1 trigger SF w.r.t Medium++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &efficienciesRel17MoriondAFIIe22vh_medium1MediumPPET;
      vectUncCorr = &uncertaintiesRel17MoriondAFIIe22vh_medium1MediumPPET;
    }
    else if (set==17) { // e22vh_medium1 trigger MC eff w.r.t Medium++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &MCefficienciesRel17MoriondAFIIe22vh_medium1MediumPPET;
      vectUncCorr = 0;
    }
    else if (set==18) { // e22vh_medium1 trigger SF w.r.t Tight++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &efficienciesRel17MoriondAFIIe22vh_medium1TightPPET;
      vectUncCorr = &uncertaintiesRel17MoriondAFIIe22vh_medium1TightPPET;
    }
    else if (set==19) { // e22vh_medium1 trigger MC eff w.r.t Tight++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &MCefficienciesRel17MoriondAFIIe22vh_medium1TightPPET;
      vectUncCorr = 0;
    }
    else if (set==20) { // e20_medium trigger MC eff w.r.t Loose++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &MCefficienciesRel17MoriondAFIIe20_mediumLoosePPET;
      vectUncCorr = 0;
    }
    else if (set==21) { // e22_medium trigger MC eff w.r.t Loose++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &MCefficienciesRel17MoriondAFIIe22_mediumLoosePPET;
      vectUncCorr = 0;
    }
    else if (set==22) { // e22vh_medium1 trigger MC eff w.r.t Loose++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &MCefficienciesRel17MoriondAFIIe22vh_medium1LoosePPET;
      vectUncCorr = 0;
    }
    else if (set==27) { // e20_medium trigger SF w.r.t Loose++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &efficienciesRel17MoriondAFIIe20_mediumLoosePPET;
      vectUncCorr = &uncertaintiesRel17MoriondAFIIe20_mediumLoosePPET;
    }
    else if (set==28) { // e22_medium trigger SF w.r.t Loose++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &efficienciesRel17MoriondAFIIe22_mediumLoosePPET;
      vectUncCorr = &uncertaintiesRel17MoriondAFIIe22_mediumLoosePPET;
    }
    else if (set==29) { // e22vh_medium1 trigger SF w.r.t Loose++ offline
      vectETBinning = &m_ETbinsTrigger;
      vectCorr = &efficienciesRel17MoriondAFIIe22vh_medium1LoosePPET;
      vectUncCorr = &uncertaintiesRel17MoriondAFIIe22vh_medium1LoosePPET;
    }

  }


  if (vectCorr == 0) { // catch all missing cases
    std::cout << "egammaSFclass: ERROR : ET-correction factors not implemented for given selection" << std::endl;
    return make_pair(-1.,-1.);
  }

  int iETbin=-1;
  while (iETbin < int(vectETBinning->size()-1)
	 && ET>=vectETBinning->at(iETbin+1))
    iETbin++;
  if (iETbin<0 || iETbin>= int(vectETBinning->size()-1)) {
    std::cout << "egammaSFclass: ERROR : given ET value (" 
	      << ET << ") outside range of existing ET-correction factors" << std::endl;
    return make_pair(-1.,-1.);
  }

  float eff=vectCorr->at(iETbin)/100.;
  float unc=0;
  if (vectUncCorr)
    unc=vectUncCorr->at(iETbin)/100.;
  return make_pair(eff, unc);
}

std::pair<float,float> egammaSFclass::scaleFactorForward(float eta, int set)
{
  if (set == 0) {
    // Forward loose
    float abseta = std::abs(eta);
    if (2.5 <= abseta && abseta <= 3.16)
      return make_pair(1.010, 0.039);
    else if (3.35 <= abseta && abseta <= 4.9)
      return make_pair(1.020, 0.059);
    else {
      std::cout << "egammaSFclass: ERROR : Out of eta range for forward electrons" << std::endl;
      return make_pair(-1.,-1.);
    }

  } else if (set == 2) {
    // Forward tight
    float abseta = std::abs(eta);
    if (2.5 <= abseta && abseta <= 3.16)
      return make_pair(0.893, 0.047);
    else if (3.35 <= abseta && abseta <= 4.9)
      return make_pair(0.940, 0.061);
    else {
      std::cout << "egammaSFclass: ERROR : Out of eta range for forward electrons" << std::endl;
      return make_pair(-1.,-1.);
    }

  } else {
    std::cout << "egammaSFclass: ERROR : Forward electrons only have Loose or Tight cuts" << std::endl;
    return make_pair(-1.,-1.);
  }

}

void egammaSFclass::copyToVector(const float *myarray, int n, std::vector<float> &dest, double renorm)
{
  for (int i = 0; i < n; i++)
    dest.push_back(myarray[i]*renorm);
}

