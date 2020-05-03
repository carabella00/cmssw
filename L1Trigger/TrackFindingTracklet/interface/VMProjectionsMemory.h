#ifndef L1Trigger_TrackFindingTracklet_interface_VMProjectionsMemory_h
#define L1Trigger_TrackFindingTracklet_interface_VMProjectionsMemory_h

#include "L1Trigger/TrackFindingTracklet/interface/MemoryBase.h"
#include "L1Trigger/TrackFindingTracklet/interface/Tracklet.h"

#include <vector>

namespace Trklet {

  class Settings;

  class VMProjectionsMemory : public MemoryBase {
  public:
    VMProjectionsMemory(std::string name, const Settings* const settings, unsigned int iSector, double phimin, double phimax);
    
    void addTracklet(Tracklet* tracklet, unsigned int allprojindex);
    
    unsigned int nTracklets() const { return tracklets_.size(); }
    
    Tracklet* getFPGATracklet(unsigned int i) const { return tracklets_[i].first; }
    int getAllProjIndex(unsigned int i) const { return tracklets_[i].second; }
    
    void writeVMPROJ(bool first);
    
    void clean() { tracklets_.clear(); }
    
    int layer() const { return layer_; }
    int disk() const { return disk_; }
    
  private:
    double phimin_;
    double phimax_;
    int layer_;
    int disk_;
    std::vector<std::pair<Tracklet*, unsigned int> > tracklets_;
  };

};
#endif
