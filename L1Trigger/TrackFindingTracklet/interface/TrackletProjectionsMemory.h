#ifndef L1Trigger_TrackFindingTracklet_interface_TrackletProjectionsMemory_h
#define L1Trigger_TrackFindingTracklet_interface_TrackletProjectionsMemory_h

#include "L1Trigger/TrackFindingTracklet/interface/MemoryBase.h"

#include <string>
#include <vector>

namespace Trklet {

  class Settings;
  class Tracklet;

  class TrackletProjectionsMemory : public MemoryBase {
  public:
    TrackletProjectionsMemory(std::string name, const Settings* const settings, unsigned int iSector);

    void addProj(Tracklet* tracklet);

    unsigned int nTracklets() const { return tracklets_.size(); }

    Tracklet* getFPGATracklet(unsigned int i) const { return tracklets_[i]; }

    void clean();

    void writeTPROJ(bool first);

    int layer() const { return layer_; }
    int disk() const { return disk_; }

  private:
    std::vector<Tracklet*> tracklets_;

    int layer_;
    int disk_;
  };

};  // namespace Trklet
#endif
