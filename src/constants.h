const double pi = 3.14159265;

const uint32_t AmbiguityResolution_Photon_Mask = 1 << 23;
const uint32_t TrackBlayer_Electron_Mask = 1 << 16;


// https://twiki.cern.ch/twiki/bin/view/AtlasProtected/LArCleaningAndObjectQuality#Details_object_quality_flag_and
// Grepability
// 0x8000000 == 134217728 == 0b1000 0000 0000 0000 0000 0000 0000

const uint32_t LARBITS_OUTOFTIME_CLUSTER = 1 << 26;
const uint32_t LARBITS_PHOTON_CLEANING = 1 << 27;

// 0x00085a6 == 34214 == 0b1000 0101 1010 0110
const uint32_t OQ_BAD_BITS = 0x00085a6;

const uint32_t nontight_mask = 0x45fc01;

const uint32_t PDGID_GLUON = 21;
