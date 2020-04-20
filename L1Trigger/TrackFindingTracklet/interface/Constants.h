#ifndef CONSTANTS_H
#define CONSTANTS_H

// ------------------------------------------------------------------------------------
// Different options for root output (imath), tracklet vs hybrid, HLS for KF
// ------------------------------------------------------------------------------------

//Uncomment if you want root output
//#define USEROOT

//Comment out to run the tracklet-only algorithm
#ifdef CMSSW_GIT_HASH
#define USEHYBRID
#endif

//Uncomment to use HLS version of KF. Also read TrackFindingTMTT/README_HLS.txt
//#ifdef USEHYBRID
//#define USE_HLS
//#endif

//Uncomment to run the HLS version of the KF if using the Hybrid (instead of the C++ KF).
//(Please also follow the instructions in L1Trigger/TrackFindingTMTT/README_HLS.txt).
//#define USE_HLS

// ------------------------------------------------------------------------------------
// Turn on/off debug info
// ------------------------------------------------------------------------------------

static const bool printDebugKF = false;  // if true print lots of debugging statements related to the KF fit
static const bool debug1 = false;        // print detailed debug information about tracklet tracking
static const bool writetrace = false;    // print out details about parsing configuration files

//Debug options in TC (should be false for 'normal' operation)
static const bool dumppars = false;
static const bool dumpproj = false;

static const bool warnNoMem = false;  // if true will print out warnings about missing projection memories
static const bool warnNoDer = false;  // if true will print out warnings about missing track fit derivatives

// ------------------------------------------------------------------------------------
// Overall configuration switches
// ------------------------------------------------------------------------------------

static const bool bookHistos =
    false;  //set to true/false to turn on/off histogram booking internal to the tracking (class "HistImp")

static unsigned int nHelixPar = 4;  // 4 or 5 param helix fit.
static bool hourglassExtended =
    false;  // turn on displaced tracking (controlled for CMSSW running via L1TrackNtupleMaker_cfg.py)

//Gemetry extensions -- used only by stand-alone code
static std::string geomext = hourglassExtended ? "hourglassExtended" : "hourglass";

static const bool geomTkTDR =
    false;  // false => newest T14/T15 tracker, true => "TDR" (T5/T6 tracker, D21/D11/D17 CMS geometries)

// Program flow (should be true for normal operation)
// enables the tracklet finding in these layer/disk combinations
static const bool doL1L2 = true;
static const bool doL2L3 = true;
static const bool doL3L4 = true;
static const bool doL5L6 = true;

static const bool doD1D2 = true;
static const bool doD3D4 = true;

static const bool doL1D1 = true;
static const bool doL2D1 = true;

static const bool doL3L4L2 = true;  // only run if hourglassExtended is true
static const bool doL5L6L4 = true;  // only run if hourglassExtended is true
static const bool doL2L3D1 = true;  // only run if hourglassExtended is true
static const bool doD1D2L2 = true;  // only run if hourglassExtended is true

static const int TMUX = 6;  //Only used for link capacity - not affecting tracking

static std::string fitpatternfile = "../data/fitpattern.txt";  //list of the different hit patterns for fits

//If this string is non-empty we will write ascii file with processed events
static const std::string skimfile = "";

// ------------------------------------------------------------------------------------
// Write various lookup tables and autogenerated code (from iMath)
// ------------------------------------------------------------------------------------

static const bool writeVerilog = false;  //Write out auto-generated Verilog mudules used by TCs
static const bool writeHLS = false;      //Write out auto-generated HLS mudules used by TCs

static const bool writeInvTable = false;  //Write out tables of drinv and invt in tracklet calculator for Verilog module
static const bool writeHLSInvTable = false;  //Write out tables of drinv and invt in tracklet calculator for HLS module

// For HLS testing: produce data/MemPrints/*/*.dat files of input/output data of processing modules.
static const bool writememLinks = false;     //Write files for dtc links
static const bool writemem = false;          //Note that for 'full' detector this will open
                                             //a LOT of files, and the program will run excruciatingly slow
static const unsigned int writememsect = 3;  //writemem only for this sector (note that the files will have _4 extension

static const bool writeVMRTables = false;      //write tables used by VMRouter
static const bool writeTripletTables = false;  //Train and write the TED and TRE tables. N.B.: the tables
                                               //cannot be applied while they are being trained, i.e.,
                                               //this flag effectively turns off the cuts in
                                               //TrackletEngineDisplaced and TripletEngine

static const bool writeTETables = false;  //LUTs used in TE
static const bool writeVMTables = false;  //LUTs used for VM consistency in TE
static const bool writeMETables = false;  //LUTS used by ME
static const bool writeMCcuts = false;    //cuts used by MC

static const bool writeFitDerTable = false;  //Write out track derivative tables

static const bool writestubs = false;      // write input stubs in the normal format
static const bool writestubs_in2 = false;  // write input stubs in hardware-ready format
static const bool padding = true;
static const bool writeoutReal = false;

// ------------------------------------------------------------------------------------
// Write out internal processing status. Scripts to make plots in PlotMacros
// ------------------------------------------------------------------------------------

static const bool writeIL = false;
static const bool writeStubsLayer = false;
static const bool writeStubsLayerperSector = false;
static const bool writeAllStubs = false;
static const bool writeVMOccupancyME = false;
static const bool writeVMOccupancyTE = false;
static const bool writeSeeds = false;
static const bool writeTE = false;
static const bool writeTED = false;
static const bool writeTRE = false;
static const bool writeTrackletProcessor = false;
static const bool writeTrackletCalculator = false;
static const bool writeTrackletCalculatorDisplaced = false;
static const bool writeTrackletPars = false;
static const bool writeAllProjections = false;
static const bool writeVMProjections = false;
static const bool writeTrackProjOcc = false;
static const bool writeME = false;
static const bool writeMatchCalculator = false;
static const bool writeResiduals = false;
static const bool writeFitTrack = false;
static const bool writeChiSq = false;

static const bool writeTC = false;  //if true write out which memories track projections will fill

static const bool writeNMatches = false;
static const bool writeHitEff = false;

static const bool writeCabling = false;
static const bool writeHitPattern = false;
static const bool writeTrackletParsOverlap = false;
static const bool writeTrackletParsDisk = false;

static const bool writeAllCT = false;  //write out .dat file containing all output tracks in bitwise format
static const bool writeifit = false;

static const bool writeVariance = false;  //write out residuals for variand matrix determination
static const bool writeResEff =
    false;                            //write files for making resolution & efficiency plots for standable code version
static const bool writePars = false;  //write files for making plots of track parameters

static const bool writeMatchEff = true;  //write files for making plots with truth matched efficiency

// ------------------------------------------------------------------------------------
// Options for chisq fit
// ------------------------------------------------------------------------------------

static const bool useMSFit = false;
static const bool tcorrection = true;
static const bool exactderivatives = false;            //for both the integer and float
static const bool exactderivativesforfloating = true;  //only for the floating point
static const bool useapprox = true;  //use approximate postion based on integer representation for floating point
static const bool usephicritapprox = false;  //use floating point approximate version of phicrit cut if true

// ------------------------------------------------------------------------------------
// Parameters for bit sizes
// ------------------------------------------------------------------------------------

static const int alphashift = 12;
static const int nbitsalpha = 4;      //bits used to store alpha
static const int alphaBitsTable = 2;  //For number of bits in track derivative table
static const int nrinvBitsTable = 3;  //number of bits for tabulating rinv dependence

static const int MEBinsBits = 3;
static const int MEBins = (1 << MEBinsBits);

static const int MEBinsDisks = 8;  //on each side

static const int Nphibits = 2;                  //Number of bits required to label the phi VM
static const int L1Nphi = (1 << Nphibits) - 1;  //Number of odd layer VMs
static const int Nzbits = 3;                    //Number of bits required to label the z VM
static const int L1Nz = (1 << Nzbits);          //Number of z VMs in odd layers
//static const int VMzbits=4;        //Number of bits for the z position in VM
static const int L2Nphi = (1 << Nphibits);  //Number of even layer VMs
static const int L2Nz = (1 << Nzbits);      //Number of z VMs in even layers
//static const int VMrbits=2;        //Number of bits for r position 'in VM'
static const int VMphibits = 3;  //Number of bits for phi position in VM

//Constants for defining stub representations
static const int nbitsrL123 = 7;
static const int nbitsrL456 = 7;

static const int nbitszL123 = 12;
static const int nbitszL456 = 8;

static const int nbitsphistubL123 = 14;
static const int nbitsphistubL456 = 17;

static const int nrbitsdisk = 12;
static const int nzbitsdisk = 7;

static const int nrbitsprojdisk = 12;
static const int nrbitsprojderdisk = 9;

static const int nbitsphiprojL123 = nbitsphistubL123;
static const int nbitsphiprojL456 = nbitsphistubL456;

static const int nbitszprojL123 = 12;
static const int nbitszprojL456 = 8;

static const int nbitsphiprojderL123 = 8 + 2;
static const int nbitsphiprojderL456 = 8 + 2;

static const int nbitszprojderL123 = 8 + 2;
static const int nbitszprojderL456 = 7 + 2;

//vm stubs
static const int nfinephibarrelinner = 2;
static const int nfinephibarrelouter = 3;

static const int nfinephidiskinner = 2;  //too small!
static const int nfinephidiskouter = 3;

static const int nfinephioverlapinner = 2;
static const int nfinephioverlapouter = 3;

//Bits used to store track parameter in tracklet
static const int nbitsrinv = 14;
static const int nbitsphi0 = 18;
static const int nbitsd0 = 13;
static const int nbitst = 14;
static const int nbitsz0 = 10;

//track and tracklet parameters
const int rinv_shift = -8;  // Krinv = 2^shift * Kphi/Kr
const int phi0_shift = 1;   // Kphi0 = 2^shift * Kphi
const int t_shift = -10;    // Kt    = 2^shift * Kz/Kr
const int z0_shift = 0;     // Kz0   = 2^shift * kz

//projections are coarsened from global to stub precision

//projection to R parameters
const int PS_phiL_shift = 0;  // phi projections have global precision in ITC
const int SS_phiL_shift = 0;
const int PS_zL_shift = 0;  // z projections have global precision in ITC
const int SS_zL_shift = 0;

const int PS_phiderL_shift = -5;  // Kderphi = 2^shift * Kphi/Kr
const int SS_phiderL_shift = -5;
const int PS_zderL_shift = -7;  // Kderz = 2^shift * Kz/Kr
const int SS_zderL_shift = -7;

//projection to Z parameters
const int PS_phiD_shift = 3;
const int SS_phiD_shift = 3;
const int PS_rD_shift = 1;  // a bug?! coarser by a factor of two then stubs??
const int SS_rD_shift = 1;

const int PS_phiderD_shift = -4;  //Kderphidisk = 2^shift * Kphi/Kz
const int SS_phiderD_shift = -4;
const int PS_rderD_shift = -6;  //Kderrdisk = 2^shift * Kr/Kz
const int SS_rderD_shift = -6;

//numbers needed for matches & fit, unclear what they are.
static const int idrinvbits = 19;
static const int phi0bitshift = 1;
static const int rinvbitshift = 13;
static const int tbitshift = 9;
static const int z0bitshift = 0;
static const int phiderbitshift = 7;
static const int zderbitshift = 6;
static const int t2bits = 23;
static const int t3shift = 8;
static const int rinvbitshiftdisk = 13;
static const int rprojdiskbitshift = 6;
static const int phiderdiskbitshift = 20;
static const int rderdiskbitshift = 7;

static const int phiresidbits = 12;
static const int zresidbits = 9;
static const int rresidbits = 7;

//Trackfit
static const int fitrinvbitshift = 9;  //6 OK?
static const int fitphi0bitshift = 6;  //4 OK?
static const int fittbitshift = 10;    //4 OK? //lower number gives rounding problems
static const int fitz0bitshift = 8;    //6 OK?

//r correction bits
static const int rcorrbits = 6;

static const int chisqphifactbits = 14;
static const int chisqzfactbits = 14;

// ------------------------------------------------------------------------------------
// Geometry
// ------------------------------------------------------------------------------------

//These define the length scale for both r and z
static const double zlength = 120.0;
static const double rmaxdisk = 120.0;

// these assume either "TDR" tracker geometry (T5 or T6), or otherwise most recent T14 (or equivalently, T15) tracker
// T5: http://cms-tklayout.web.cern.ch/cms-tklayout/layouts/recent-layouts/OT616_200_IT404/layout.html
// T14: http://cms-tklayout.web.cern.ch/cms-tklayout/layouts/recent-layouts/OT616_200_IT404/layout.html

static const double rmeanL1 = geomTkTDR ? (rmaxdisk * 858) / 4096 : (rmaxdisk * 851) / 4096;
static const double rmeanL2 = geomTkTDR ? (rmaxdisk * 1279) / 4096 : (rmaxdisk * 1269) / 4096;
static const double rmeanL3 = geomTkTDR ? (rmaxdisk * 1795) / 4096 : (rmaxdisk * 1784) / 4096;
static const double rmeanL4 = geomTkTDR ? (rmaxdisk * 2347) / 4096 : (rmaxdisk * 2347) / 4096;
static const double rmeanL5 = geomTkTDR ? (rmaxdisk * 2937) / 4096 : (rmaxdisk * 2936) / 4096;
static const double rmeanL6 = geomTkTDR ? (rmaxdisk * 3783) / 4096 : (rmaxdisk * 3697) / 4096;

static const double zmeanD1 = (zlength * 2239) / 2048;
static const double zmeanD2 = (zlength * 2645) / 2048;
static const double zmeanD3 = (zlength * 3163) / 2048;
static const double zmeanD4 = (zlength * 3782) / 2048;
static const double zmeanD5 = (zlength * 4523) / 2048;

static const double rmindiskvm = 22.5;
static const double rmaxdiskvm = 67.0;

static const double rmaxdiskl1overlapvm = 45.0;
static const double rmindiskl2overlapvm = 40.0;
static const double rmindiskl3overlapvm = 50.0;

static const double half2SmoduleWidth = 4.57;

// need separate lookup values for inner two vs outer three disks for 2S modules

// T5:  http://cms-tklayout.web.cern.ch/cms-tklayout/layouts/cmssw-models/I_OT613_200_IT4025/index.html
// T14: http://cms-tklayout.web.cern.ch/cms-tklayout/layouts/cmssw-models/Q_OT616_200_IT613/index.html

static const double rDSSinner_mod1 = geomTkTDR ? 69.2345 : 68.9391;
static const double rDSSinner_mod2 = geomTkTDR ? 80.0056 : 78.7750;
static const double rDSSinner_mod3 = geomTkTDR ? 87.3444 : 85.4550;
static const double rDSSinner_mod4 = geomTkTDR ? 98.2515 : 96.3150;
static const double rDSSinner_mod5 = geomTkTDR ? 104.9750 : 102.3160;

static const double rDSSouter_mod1 = geomTkTDR ? 67.6317 : 66.4903;
static const double rDSSouter_mod2 = geomTkTDR ? 78.1300 : 76.7750;
static const double rDSSouter_mod3 = geomTkTDR ? 86.4293 : 84.4562;
static const double rDSSouter_mod4 = geomTkTDR ? 97.1316 : 94.9920;
static const double rDSSouter_mod5 = geomTkTDR ? 104.9750 : 102.3160;

static const double halfstrip =
    2.5;  //we want the center of the two strip positions in a module, not just the center of a module

static const double rDSSinner[10] = {rDSSinner_mod1 - halfstrip,
                                     rDSSinner_mod1 + halfstrip,
                                     rDSSinner_mod2 - halfstrip,
                                     rDSSinner_mod2 + halfstrip,
                                     rDSSinner_mod3 - halfstrip,
                                     rDSSinner_mod3 + halfstrip,
                                     rDSSinner_mod4 - halfstrip,
                                     rDSSinner_mod4 + halfstrip,
                                     rDSSinner_mod5 - halfstrip,
                                     rDSSinner_mod5 + halfstrip};
static const double rDSSouter[10] = {rDSSouter_mod1 - halfstrip,
                                     rDSSouter_mod1 + halfstrip,
                                     rDSSouter_mod2 - halfstrip,
                                     rDSSouter_mod2 + halfstrip,
                                     rDSSouter_mod3 - halfstrip,
                                     rDSSouter_mod3 + halfstrip,
                                     rDSSouter_mod4 - halfstrip,
                                     rDSSouter_mod4 + halfstrip,
                                     rDSSouter_mod5 - halfstrip,
                                     rDSSouter_mod5 + halfstrip};

static const double drmax = rmaxdisk / 32.0;

static const double dzmax = zlength / 32.0;

static const double drdisk = rmaxdisk;

static const double rmean[6] = {rmeanL1, rmeanL2, rmeanL3, rmeanL4, rmeanL5, rmeanL6};

static const double zmean[5] = {zmeanD1, zmeanD2, zmeanD3, zmeanD4, zmeanD5};

//Number of VM regions for TE by inner/outer layer for each seed
static const unsigned int nvmte_[3][12] = {
    {4, 4, 4, 4, 4, 4, 2, 2, 4, 4, 8, 4}, {8, 4, 8, 8, 4, 4, 4, 4, 8, 8, 4, 4}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 2}};
//number of bits for VM regions for TE by inner/outer layer for each seed
static const unsigned int nbitsvmte[3][12] = {
    {2, 2, 2, 2, 2, 2, 1, 1, 2, 2, 3, 2}, {3, 2, 3, 3, 2, 2, 2, 2, 3, 3, 2, 2}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1}};

//Number of allstub memories (phi regions) for each of the layers+disks
static const unsigned int nallstubs_[11] = {8, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4};
static const unsigned int nbitsallstubs_[11] = {3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2};
//Number of vm for ME memories per phi region for each of the layers+disks
static const unsigned int nvmme[11] = {4, 8, 8, 8, 8, 8, 8, 4, 4, 4, 4};
static const unsigned int nbitsvmme[11] = {2, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2};

static const unsigned int nallstubslayers[6] = {8, 4, 4, 4, 4, 4};   // FIXME should retire
static const unsigned int nvmtelayers[6] = {4, 8, 4, 8, 4, 8};       // FIXME should retire
static const unsigned int nvmteextralayers[6] = {0, 4, 4, 0, 0, 0};  // FIXME should retire

static const unsigned int nallprojlayers[6] = {8, 4, 4, 4, 4, 4};  // FIXME should retire
static const unsigned int nvmmelayers[6] = {4, 8, 8, 8, 8, 8};     // FIXME should retire

static const unsigned int nallstubsdisks[5] = {4, 4, 4, 4, 4};  // FIXME should retire
static const unsigned int nvmtedisks[5] = {4, 4, 4, 4, 4};      // FIXME should retire

static const unsigned int nallprojdisks[5] = {4, 4, 4, 4, 4};  // FIXME should retire
static const unsigned int nvmmedisks[5] = {8, 4, 4, 4, 4};     // FIXME should retire
//for seeding in L1D1 L2D1
static const unsigned int nallstubsoverlaplayers[3] = {8, 4, 4};  // FIXME should retire
static const unsigned int nvmteoverlaplayers[3] = {2, 2, 2};      // FIXME should retire

static const unsigned int nallstubsoverlapdisks[2] = {4, 4};  // FIXME should retire
static const unsigned int nvmteoverlapdisks[2] = {4, 4};      // FIXME should retire

static const double rcrit = 55.0;

static const double rmaxL1 = rmeanL1 + drmax;
static const double rmaxL2 = rmeanL2 + drmax;
static const double rmaxL3 = rmeanL3 + drmax;
static const double rmaxL4 = rmeanL4 + drmax;
static const double rmaxL5 = rmeanL5 + drmax;
static const double rmaxL6 = rmeanL6 + drmax;

static const double zminD1 = zmeanD1 - dzmax;
static const double zmaxD1 = zmeanD1 + dzmax;
static const double zminD2 = zmeanD2 - dzmax;
static const double zmaxD2 = zmeanD2 + dzmax;
static const double zminD3 = zmeanD3 - dzmax;
static const double zmaxD3 = zmeanD3 + dzmax;
static const double zminD4 = zmeanD4 - dzmax;
static const double zmaxD4 = zmeanD4 + dzmax;
static const double zminD5 = zmeanD5 - dzmax;
static const double zmaxD5 = zmeanD5 + dzmax;

static const double two_pi = 2 * M_PI;

static const double ptcut = 1.91;                            //Minimum pt
static const double rinvcut = 0.01 * 0.3 * 3.8 / ptcut;      //0.01 to convert to cm-1
static const double ptcutte = 1.8;                           //Minimum pt in TE
static const double rinvcutte = 0.01 * 0.3 * 3.8 / ptcutte;  //0.01 to convert to cm-1 in TE
static const double bendcut = 1.25;
static const double bendcutdisk = 1.25;
static const double z0cut = 15.0;
static const double mecut = 2.0;
static const double mecutdisk = 1.5;

static const unsigned int NSector = 9;

static const double rinvmax = 0.01 * 0.3 * 3.8 / 2.0;  //0.01 to convert to cm-1

static const double dphisectorHG =
  2 * M_PI / NSector + 2 * fmax(std::abs(asin(0.5 * rinvmax * rmean[0]) - asin(0.5 * rinvmax * rcrit)),
				std::abs(asin(0.5 * rinvmax * rmean[5]) - asin(0.5 * rinvmax * rcrit)));

static const double phicritmin = 0.5 * dphisectorHG - M_PI / NSector;
static const double phicritmax = dphisectorHG - 0.5 * dphisectorHG + M_PI / NSector;

static const double dphicritmc = 0.005;  //lose for MC
static const double phicritminmc = phicritmin - dphicritmc;
static const double phicritmaxmc = phicritmax + dphicritmc;

// Obsolete - only used in TrackletCalculatorDisplaced (Ryd - 2020-01-16)
static const int iphicritminmc = 9253;
static const int iphicritmaxmc = 56269;

static const unsigned int NLONGVMBITS = 3;
static const unsigned int NLONGVMRBITS = 3;  //4 bins on each side (+-z)
static const unsigned int NLONGVMBINS = (1 << NLONGVMBITS);
static const unsigned int NLONGVMRBINS = (1 << NLONGVMRBITS);
static const unsigned int NLONGVMODDLAYERBITS = 6;
static const unsigned int NLONGVMODDDISKBITS = 6;

// ------------------------------------------------------------------------------------
// Truncation cuts
// TMUX = 18 & assuming a clock frequency of 240 MHz (conservative) ==>
// 18*(240 MHz)/(40 MHz) = 108
// ------------------------------------------------------------------------------------

static const unsigned int MAXOFFSET = 10000;  //set to 0 to enable regular truncation or 10000 to disable it.

static const unsigned int MAXSTUBSLINK = 108 + MAXOFFSET;    //Max stubs per link
static const unsigned int MAXLAYERROUTER = 108 + MAXOFFSET;  //Max stubs handled by layer router
static const unsigned int MAXDISKROUTER = 108 + MAXOFFSET;   //Max stubs handled by disk router
static const unsigned int MAXVMROUTER = 108 + MAXOFFSET;     //Max stubs handled by VM router
static const unsigned int MAXTE = 108 + MAXOFFSET;           //Maximum number of stub pairs to try in TE
static const unsigned int MAXTRE = 108 + MAXOFFSET;          //Maximum number of stub pairs to try in TRE
static const unsigned int MAXTC = 108 + MAXOFFSET;           //Maximum number of tracklet parameter calculations
static const unsigned int MAXPROJROUTER = 108 + MAXOFFSET;   //Maximum number of projections to route
static const unsigned int MAXME = 108 + MAXOFFSET;           //Maximum number of stub-projection matches to try
static const unsigned int MAXMC = 108 + MAXOFFSET;           //Maximum number of match calculations
static const unsigned int MAXMP = 108 + MAXOFFSET;           //Maximum number of match calculations
static const unsigned int MAXFIT = 108 + MAXOFFSET;          //Maximum number of track fits

// ------------------------------------------------------------------------------------
// additional track parameter constants
// ------------------------------------------------------------------------------------

// phi sector definition
static const double dphisector = 2 * M_PI / NSector;

//Minimal ranges for track parameters
static const double maxrinv = 0.006;
static const double maxd0 = 10.;

//These are constants defining global coordinate system
static const double kphi = dphisectorHG / (1 << nbitsphistubL123);
static const double kphi1 = dphisectorHG / (1 << nbitsphistubL456);
static const double kz = 2 * zlength / (1 << nbitszL123);
//static const double kr=2*drmax/(1<<nbitsrL456);
static const double kr = rmaxdisk / (1 << nrbitsdisk);
static const double kd0 = 2 * maxd0 / (1 << nbitsd0);

//constants derivative from the above
static double krinvpars, kphi0pars, kd0pars, ktpars, kz0pars;
static double kphiproj123, kphiproj456, kzproj, kphider, kzder;
static double krprojshiftdisk, kphiprojdisk, krprojderdisk;
static double krdisk, krprojderdiskshift, kzpars;

// ------------------------------------------------------------------------------------
// Duplicate Removal
// ------------------------------------------------------------------------------------

static const int minIndStubs = 3;  // not used with merge removal

// possible options
// "ichi" (pairwise, keep track with best ichisq)
// "nstub" (pairwise, keep track with more stubs)
// "grid" (TMTT-like removal)
// "" (no removal)
// "merge" (hybrid dup removal)

#ifdef USEHYBRID
static const std::string RemovalType = "merge";
// "CompareBest" (recommended) Compares only the best stub in each track for each region (best = smallest phi residual) and will merge the two tracks if stubs are shared in three or more regions
// "CompareAll" Compares all stubs in a region, looking for matches, and will merge the two tracks if stubs are shared in three or more regions
static const std::string MergeComparison = "CompareBest";
static const bool doKF = true;  //true => use KF (USEHYBRID is defined)
#else
static const std::string RemovalType = "ichi";
static const bool doKF = false;  //false => use chi2 fit (USEHYBRID is not defined)
#endif

static const bool fakefit =
    false;  //if true, run a dummy fit, producing TTracks directly from output of tracklet pattern reco stage. (Not working for Hybrid)

// ------------------------------------------------------------------------------------
// bit definitions internal to the tracklet algorithm
// ------------------------------------------------------------------------------------

//These are the number of bits used for the fine phi position fo the TE VM stubs by seedindex
static const unsigned int nfinephi_[3][12] = {{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2},   //inner
                                              {3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3},   //outer
                                              {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 3}};  //outermost

//These are the number of bits used for the VM regions in the TE by seedindex
static const unsigned int nphireg_[3][12] = {{5, 4, 4, 4, 4, 4, 4, 3, 4, 4, 5, 4},   //inner
                                             {5, 4, 5, 5, 4, 4, 4, 4, 4, 4, 4, 4},   //outer
                                             {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 4}};  //outermost

static const unsigned int zbitstab[3][12] = {
    {7, 7, 7, 7, 3, 3, 7, 7, 0, 0, 7, 0}, {7, 7, 7, 7, 3, 3, 3, 3, 0, 0, 7, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 7}};

static const unsigned int rbitstab[3][12] = {
    {4, 4, 4, 4, 8, 8, 3, 3, 0, 0, 4, 0}, {4, 4, 4, 4, 7, 7, 7, 7, 0, 0, 4, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 7, 4}};

static const unsigned int lutwidthtab[3][12] = {{10, 11, 11, 11, 11, 11, 11, 11, 0, 0, 11, 0},
                                                {6, 6, 6, 6, 10, 10, 10, 10, 0, 0, 6, 0},
                                                {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 6}};

static const unsigned int lutwidthtabextended[3][12] = {{11, 11, 21, 21, 21, 21, 11, 11, 0, 0, 21, 0},
                                                        {6, 6, 6, 6, 10, 10, 10, 10, 0, 0, 6, 0},
                                                        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 6}};

//projection layers by seed index. For each seeding index (row) the list of layers that
//we consider projections to
static const int projlayers[12][4] = {
    {3, 4, 5, 6},  //0 L1L2
    {1, 4, 5, 6},  //1 L2L3
    {1, 2, 5, 6},  //2 L3L4
    {1, 2, 3, 4},  //3 L5L6
    {1, 2},        //4 D1D2
    {1},           //5 D3D4
    {},            //6 L1D1
    {1},           //7 L2D1
    {1, 5, 6},     //8 L2L3L4
    {1, 2, 3},     //9 L4L5L6
    {1},           //10 L2L3D1
    {1}            //11 D1D2L2
};

//projection disks by seed index. For each seeding index (row) the list of diks that
//we consider projections to
static const int projdisks[12][5] = {
    {1, 2, 3, 4},  //0 L1L2
    {1, 2, 3, 4},  //1 L2L3
    {1, 2},        //2 L3L4
    {},            //3 L5L6
    {3, 4, 5},     //4 D1D2
    {1, 2, 5},     //5 D3D4
    {2, 3, 4, 5},  //6 L1D1
    {2, 3, 4},     //7 L2D1
    {1, 2},        //8 L2L3L4
    {},            //9 L4L5L6
    {2, 3, 4},     //10 L2L3D1
    {3, 4}         //11 D1D2L2
};

//rphi cuts for layers - the column is the seedindex
static const double rphimatchcut[6][12] = {
    {0.0, 0.1, 0.07, 0.08, 0.07, 0.05, 0.0, 0.05, 0.08, 0.15, 0.125, 0.15},  //Layer 1
    {0.0, 0.0, 0.06, 0.08, 0.05, 0.0, 0.0, 0.0, 0.0, 0.1, 0.0, 0.0},         //Layer 2
    {0.1, 0.0, 0.0, 0.08, 0.0, 0.0, 0.0, 0.0, 0.0, 0.08, 0.0, 0.0},          //Layer 3
    {0.19, 0.19, 0.0, 0.05, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},         //Layer 4
    {0.4, 0.4, 0.08, 0.0, 0.0, 0.0, 0.0, 0.0, 0.08, 0.0, 0.0, 0.0},          //Layer 5
    {0.5, 0.0, 0.19, 0.0, 0.0, 0.0, 0.0, 0.0, 0.2, 0.0, 0.0, 0.0}            //Layer 6
};

//z cuts for layers - the column is the seedindex
static const double zmatchcut[6][12] = {
    {0.0, 0.7, 5.5, 15.0, 1.5, 2.0, 0.0, 1.5, 1.0, 8.0, 1.0, 1.5},   //Layer 1
    {0.0, 0.0, 3.5, 15.0, 1.25, 0.0, 0.0, 0.0, 0.0, 7.0, 0.0, 0.0},  //Layer 2
    {0.7, 0.0, 0.0, 9.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.0, 0.0, 0.0},    //Layer 3
    {3.0, 3.0, 0.0, 7.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},    //Layer 4
    {3.0, 3.0, 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.5, 0.0, 0.0, 0.0},    //Layer 5
    {4.0, 0.0, 9.5, 0.0, 0.0, 0.0, 0.0, 0.0, 4.5, 0.0, 0.0, 0.0}     //Layer 6
};

//rphi cuts for PS modules in disks - the column is the seedindex
static const double rphicutPS[5][12] = {
    {0.2, 0.2, 0.0, 0.0, 0.0, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},     //disk 1
    {0.2, 0.2, 0.0, 0.0, 0.0, 0.1, 0.1, 0.1, 0.0, 0.0, 0.15, 0.0},    //disk 2
    {0.25, 0.2, 0.0, 0.0, 0.15, 0.0, 0.2, 0.15, 0.0, 0.0, 0.0, 0.2},  //disk 3
    {0.5, 0.2, 0.0, 0.0, 0.2, 0.0, 0.3, 0.5, 0.0, 0.0, 0.0, 0.0},     //disk 4
    {0.0, 0.0, 0.0, 0.0, 0.25, 0.1, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0}     //disk 5
};

//r cuts for PS modules in disks - the column is the seedindex
static const double rcutPS[5][12] = {
    {0.5, 0.5, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},  //disk 1
    {0.5, 0.5, 0.0, 0.0, 0.0, 0.5, 0.5, 0.5, 0.0, 0.0, 0.5, 0.0},  //disk 2
    {0.5, 0.5, 0.0, 0.0, 0.5, 0.0, 0.6, 0.8, 0.0, 0.0, 0.0, 0.4},  //disk 3
    {0.5, 0.5, 0.0, 0.0, 0.8, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0},  //disk 4
    {0.0, 0.0, 0.0, 0.0, 1.0, 0.5, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0}   //disk 5
};

//rphi cuts for 2S modules in disks = the column is the seedindex
static const double rphicut2S[5][12] = {
    {0.5, 0.5, 0.8, 0.0, 0.0, 0.0, 0.0, 0.0, 0.2, 0.0, 0.0, 0.0},    //disk 1
    {0.5, 0.5, 0.8, 0.0, 0.0, 0.0, 0.5, 0.15, 0.3, 0.0, 0.68, 0.0},  //disk 2
    {0.5, 0.5, 0.0, 0.0, 0.15, 0.0, 0.2, 0.25, 0.0, 0.0, 0.8, 0.1},  //disk 3
    {0.5, 0.5, 0.0, 0.0, 0.2, 0.0, 0.25, 0.5, 0.0, 0.0, 0.6, 0.4},   //disk 4
    {0.0, 0.0, 0.0, 0.0, 0.4, 0.2, 0.4, 0.0, 0.0, 0.0, 0.0, 0.8}     //disk 5
};

//r cuts for 2S modules in disks -the column is the seedindex
static const double rcut2S[5][12] = {
    {3.8, 3.8, 3.8, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0},  //disk 1
    {3.8, 3.8, 3.8, 0.0, 0.0, 0.0, 3.8, 3.4, 3.0, 0.0, 3.0, 0.0},  //disk 2
    {3.6, 3.8, 0.0, 0.0, 3.6, 0.0, 3.6, 3.8, 0.0, 0.0, 3.8, 3.0},  //disk 3
    {3.6, 3.8, 0.0, 0.0, 3.6, 0.0, 3.5, 3.8, 0.0, 0.0, 3.0, 3.0},  //disk 4
    {0.0, 0.0, 0.0, 0.0, 3.6, 3.4, 3.7, 0.0, 0.0, 0.0, 0.0, 3.0}   //disk 5
};

#endif
