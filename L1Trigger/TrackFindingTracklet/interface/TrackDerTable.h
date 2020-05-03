#ifndef L1Trigger_TrackFindingTracklet_interface_TrackDerTable_h
#define L1Trigger_TrackFindingTracklet_interface_TrackDerTable_h

#include <iostream>
#include <fstream>
#include <cassert>
#include <cmath>
#include <vector>
#include <string>

#include "L1Trigger/TrackFindingTracklet/interface/TrackDer.h"

namespace Trklet {

  class Settings;
  class Globals;

  class TrackDerTable {
  public:
    TrackDerTable(const Settings* settings);
    
    ~TrackDerTable() {}

    TrackDer* getDerivatives(int index) { return &derivatives_[index]; }

    TrackDer* getDerivatives(unsigned int layermask,unsigned int diskmask,unsigned int alphaindex,unsigned int rinvindex);
    
    int getIndex(unsigned int layermask, unsigned int diskmask);
    
    void addEntry(unsigned int layermask, unsigned int diskmask, int multiplicity, int nrinv);
    
    void readPatternFile(std::string fileName);
    
    int getEntries() const { return nextLayerDiskValue_; }
    
    void fillTable(const Settings* settings);

    static void invert(double M[4][8], unsigned int n);
    
    static void invert(std::vector<std::vector<double> >& M, unsigned int n);
    
    static void getVarianceMatrix(Globals* globals, bool layer[6], bool disk[5], int ptbin, std::vector<std::vector<double> >& V);
    
    static void calculateDerivatives(const Settings* settings,
				     unsigned int nlayers,
				     double r[6],
				     unsigned int ndisks,
				     double z[5],
				     double alpha[5],
				     double t,
				     double rinv,
				     double D[4][12],
				     int iD[4][12],
				     double MinvDt[4][12],
				     int iMinvDt[4][12],
				     double sigma[12],
				     double kfactor[12]);
    
    static void calculateDerivativesMS(const Settings* settings,
				       Globals* globals,
				       unsigned int nlayers,
				       double r[6],
				       unsigned int ndisks,
				       double z[5],
				       double alpha[5],
				       double t,
				       double rinv,
				       double D[4][12],
				       int iD[4][12],
				       double MinvDt[4][12],
				       int iMinvDt[4][12],
				       double sigma[12],
				       double kfactor[12],
				       int ptbin);

    
    static double gett(const Settings* settings, int diskmask, int layermask);

    
  private:

    const Settings* settings_;
    
    std::vector<int> LayerMem_;
    std::vector<int> DiskMem_;    
    std::vector<int> LayerDiskMem_;
    
    unsigned int LayerMemBits_;
    unsigned int DiskMemBits_;
    unsigned int LayerDiskMemBits_;
    unsigned int alphaBits_;
    
    unsigned int Nlay_;
    unsigned int Ndisk_;
    
    std::vector<TrackDer> derivatives_;
    
    int nextLayerValue_;
    int nextDiskValue_;
    int nextLayerDiskValue_;
    int lastMultiplicity_;
    
  };

};
#endif
