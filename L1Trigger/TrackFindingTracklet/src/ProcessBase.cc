#include "L1Trigger/TrackFindingTracklet/interface/ProcessBase.h"
#include "L1Trigger/TrackFindingTracklet/interface/MemoryBase.h"
#include "L1Trigger/TrackFindingTracklet/interface/Settings.h"
#include "L1Trigger/TrackFindingTracklet/interface/Globals.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"

using namespace trklet;
using namespace std;

ProcessBase::ProcessBase(string name, const Settings* const settings, Globals* global, unsigned int iSector)
  : name_(name), settings_(settings), globals_(global) {
  iSector_ = iSector;
  double dphi = 2 * M_PI / settings_->NSector();
  double dphiHG = 0.5 * settings_->dphisectorHG() - M_PI / settings_->NSector();
  phimin_ = iSector_ * dphi - dphiHG;
  phimax_ = phimin_ + dphi + 2 * dphiHG;
  phimin_ -= M_PI / settings_->NSector();
  phimax_ -= M_PI / settings_->NSector();
  if (phimin_ > M_PI) {
    phimin_ -= 2 * M_PI;
    phimax_ -= 2 * M_PI;
  }
}

unsigned int ProcessBase::nbits(unsigned int power) {
  if (power == 2)
    return 1;
  if (power == 4)
    return 2;
  if (power == 8)
    return 3;
  if (power == 16)
    return 4;
  if (power == 32)
    return 5;

  throw cms::Exception("LogicError")  << __FILE__ << " " << __LINE__ << "nbits: power = " << power;
  return -1;
}

void ProcessBase::initLayerDisk(unsigned int pos, int& layer, int& disk) {
  string subname = name_.substr(pos, 2);
  layer = 0;
  disk = 0;
  if (subname.substr(0,1)=="L")
    layer = stoi(subname.substr(1,1));
  else if (subname.substr(0,1)=="D")
    disk = stoi(subname.substr(1,1));
  else
    throw cms::Exception("BadConfig") << __FILE__ << " " << __LINE__
				      << " " << name_ << " subname = " << subname << " " << layer << " " << disk;
}

void ProcessBase::initLayerDisk(unsigned int pos, int& layer, int& disk, int& layerdisk) {
  initLayerDisk(pos, layer, disk);

  layerdisk = layer - 1;
  if (disk > 0)
    layerdisk = 5 + disk;
}

unsigned int ProcessBase::initLayerDisk(unsigned int pos) {
  int layer, disk;
  initLayerDisk(pos, layer, disk);

  if (disk > 0)
    return 5 + disk;
  return layer - 1;
}

void ProcessBase::initLayerDisksandISeed(unsigned int& layerdisk1, unsigned int& layerdisk2, unsigned int& iSeed) {
  layerdisk1 = 99;
  layerdisk2 = 99;

  if (name_.substr(0,3)=="TE_") {
    if (name_[3] == 'L') {
      layerdisk1 = name_[4] - '1';
    }
    else if (name_[3] == 'D') {
      layerdisk1 = 6 + name_[4] - '1';
    }
    if (name_[11] == 'L') {
      layerdisk2 = name_[12] - '1';
    }
    else if (name_[11] == 'D') {
      layerdisk2 = 6 + name_[12] - '1';
    }
    else if (name_[12] == 'L') {
      layerdisk2 = name_[13] - '1';
    }
    else if (name_[12] == 'D') {
      layerdisk2 = 6 + name_[13] - '1';
    }  
  }
  
  if (name_.substr(0,3)=="TC_") {
    if (name_[3] == 'L') {
      layerdisk1 = name_[4] - '1';
    }
    else if (name_[3] == 'D') {
      layerdisk1 = 6 + name_[4] - '1';
    }
    if (name_[5] == 'L') {
      layerdisk2 = name_[6] - '1';
    }
    else if (name_[5] == 'D') {
      layerdisk2 = 6 + name_[6] - '1';
    }
  }
  
  if (layerdisk1 == 0 && layerdisk2 == 1)
    iSeed = 0;
  else if (layerdisk1 == 1 && layerdisk2 == 2)
    iSeed = 1;
  else if (layerdisk1 == 2 && layerdisk2 == 3)
    iSeed = 2;
  else if (layerdisk1 == 4 && layerdisk2 == 5)
    iSeed = 3;
  else if (layerdisk1 == 6 && layerdisk2 == 7)
    iSeed = 4;
  else if (layerdisk1 == 8 && layerdisk2 == 9)
    iSeed = 5;
  else if (layerdisk1 == 0 && layerdisk2 == 6)
    iSeed = 6;
  else if (layerdisk1 == 1 && layerdisk2 == 6)
    iSeed = 7;
  else {
    assert(0);
  }
    
}

unsigned int ProcessBase::getISeed(std::string name) {
  //assumes here that namme is on the form XX_L1L2_XXX where L1L2 gives iSeed=0

  std::size_t pos = name.find("_");
  std::string name1 = name.substr(pos + 1);

  pos = name1.find("_");
  std::string name2 = name1.substr(0, pos);

  if (name2 == "L1L2")
    return 0;
  if (name2 == "L2L3")
    return 1;
  if (name2 == "L3L4")
    return 2;
  if (name2 == "L5L6")
    return 3;
  if (name2 == "D1D2")
    return 4;
  if (name2 == "D3D4")
    return 5;
  if (name2 == "L1D1")
    return 6;
  if (name2 == "L2D1")
    return 7;

  if (name2 == "L1L2XX")
    return 0;
  if (name2 == "L2L3XX")
    return 1;
  if (name2 == "L3L4XX")
    return 2;
  if (name2 == "L5L6XX")
    return 3;
  if (name2 == "D1D2XX")
    return 4;
  if (name2 == "D3D4XX")
    return 5;
  if (name2 == "L1D1XX")
    return 6;
  if (name2 == "L2D1XX")
    return 7;
  if (name2 == "L3L4L2")
    return 8;
  if (name2 == "L5L6L4")
    return 9;
  if (name2 == "L2L3D1")
    return 10;
  if (name2 == "D1D2L2")
    return 11;

  edm::LogPrint("Tracklet") << getName() << " name name1 name2 " << name << " - " << name1 << " - " << name2;
  assert(0);
  return 0;
}
