#ifndef L1Trigger_TrackFindingTracklet_interface_MatchProcessor_h
#define L1Trigger_TrackFindingTracklet_interface_MatchProcessor_h

#include "L1Trigger/TrackFindingTracklet/interface/ProcessBase.h"
#include "L1Trigger/TrackFindingTracklet/interface/MatchEngineUnit.h"
#include "L1Trigger/TrackFindingTracklet/interface/ProjectionTemp.h"
#include "L1Trigger/TrackFindingTracklet/interface/CircularBuffer.h"
#include "L1Trigger/TrackFindingTracklet/interface/FullMatchMemory.h"
#include "L1Trigger/TrackFindingTracklet/interface/TrackletProjectionsMemory.h"
#include "L1Trigger/TrackFindingTracklet/interface/VMStubsMEMemory.h"
#include "L1Trigger/TrackFindingTracklet/interface/AllStubsMemory.h"

#include <vector>

class L1TStub;

namespace Trklet {

  class Settings;
  class Globals;
  class MemoryBase;
  class Stub;
  class Tracklet;

  class MatchProcessor:public ProcessBase {
    
  public:
    
    MatchProcessor(string name, const Settings* settings, Globals* global, unsigned int iSector);

    void addOutput(MemoryBase* memory,string output);

    void addInput(MemoryBase* memory,string input);

    void execute();
    
    bool matchCalculator(Tracklet* tracklet,Stub* fpgastub, L1TStub* stub);
    
  private:
    
    int layer_;
    int disk_;
    bool barrel_;

    int nvm_;  //VMs in sector
    int nvmbits_; //# of bits for VMs in sector
    int nvmbins_; //VMs in in phi region
    
    int fact_;
    int icorrshift_;
    int icorzshift_;
    int phi0shift_;
    
    double phimin_;
    double phimax_;
    double phioffset_;
    
    unsigned int phimatchcut_[12];
    unsigned int zmatchcut_[12];
    
    unsigned int rphicutPS_[12];
    unsigned int rphicut2S_[12];
    unsigned int rcutPS_[12];
    unsigned int rcut2S_[12];
    
    double phifact_;
    double rzfact_;
    
    int nrbits_;
    int nphiderbits_;
    
    AllStubsMemory* allstubs_;
    std::vector<VMStubsMEMemory*> vmstubs_;
    std::vector<TrackletProjectionsMemory*> inputprojs_;
    
    int ialphafactinner_[10];
    int ialphafactouter_[10];
    
    //FIXME should index by iSeed
    std::vector<FullMatchMemory*> fullmatches_;
    
    //used in the layers
    std::vector<bool> table_;
    
    //used in the disks
    std::vector<bool> tablePS_;
    std::vector<bool> table2S_;
    
    unsigned int nMatchEngines_;  
    std::vector<MatchEngineUnit> matchengines_;
    
    CircularBuffer<ProjectionTemp> inputProjBuffer_;
  };

};
#endif
