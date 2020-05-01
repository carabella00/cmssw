#ifndef L1Trigger_TrackFindingTracklet_interface_Track_h
#define L1Trigger_TrackFindingTracklet_interface_Track_h

#include <iostream>
#include <fstream>
#include <cassert>
#include <cmath>
#include <vector>
#include <map>

#include "L1Trigger/TrackFindingTracklet/interface/Settings.h"
#include "L1Trigger/TrackFindingTracklet/interface/Util.h"
#include "L1Trigger/TrackFindingTracklet/interface/slhcevent.h"

namespace Trklet {
  
  class Track {
  public:
    Track(int irinv,
	  int iphi0,
	  int id0,
	  int it,
	  int iz0,
	  int ichisqrphi,
	  int ichisqrz,
	  double chisqrphi,
	  double chisqrz,
	  int hitpattern,
	  std::map<int, int> stubID,
	  std::vector<L1TStub*> l1stub,
	  int seed);
    
    ~Track() {}
    
    void setDuplicate(bool flag) { duplicate_ = flag; }
    void setSector(int nsec) { sector_ = nsec; }
    void setStubIDpremerge(std::vector<std::pair<int, int>> stubIDpremerge) { stubIDpremerge_ = stubIDpremerge; }
    void setStubIDprefit(std::vector<std::pair<int, int>> stubIDprefit) { stubIDprefit_ = stubIDprefit; }
    
    int irinv() const { return irinv_; }
    int iphi0() const { return iphi0_; }
    int id0() const { return id0_; }
    int iz0() const { return iz0_; }
    int it() const { return it_; }
    int ichisq() const { return ichisqrphi_ + ichisqrz_; }
    
    std::map<int, int> stubID() const { return stubID_; }
    std::vector<L1TStub*> stubs() const { return l1stub_; }
    std::vector<std::pair<int, int>> stubIDpremerge() const { return stubIDpremerge_; }
    std::vector<std::pair<int, int>> stubIDprefit() const { return stubIDprefit_; }
    
    int hitpattern() const { return hitpattern_; }
    int seed() const { return seed_; }
    int duplicate() const { return duplicate_; }
    int sector() const { return sector_; }
    
    double pt(const Settings* settings, double bfield = 3.811202) const { return (0.3 * bfield / 100.0) / (irinv_ * settings->krinvpars()); }

    double phi0(const Settings* settings) const;

    double eta(const Settings* settings) const { return asinh(it_ * settings->ktpars()); }
    double tanL(const Settings* settings) const { return it_ * settings->ktpars(); }
    double z0(const Settings* settings) const { return iz0_ * settings->kz0pars(); }
    double rinv(const Settings* settings) const { return irinv_ * settings->krinvpars(); }
    double d0(const Settings* settings) const { return id0_ * settings->kd0pars(); }  //Fix when fit for 5 pars
    double chisq() const { return chisqrphi_ + chisqrz_; }
    
    double chisqrphi() const { return chisqrphi_; }
    double chisqrz() const { return chisqrz_; }
    
    int nPSstubs() const {
      int npsstubs = 0;
      for (unsigned int i = 0; i < l1stub_.size(); i++) {
	if (l1stub_[i]->layer() < 3)
	  npsstubs++;
      }
      return npsstubs;
    }
    
  private:
    int irinv_;
    int iphi0_;
    int id0_;
    int iz0_;
    int it_;
    int ichisqrphi_;
    int ichisqrz_;
    
    double chisqrphi_;
    double chisqrz_;
    
    int hitpattern_;
    
    std::vector<std::pair<int, int>> stubIDpremerge_;
    std::vector<std::pair<int, int>> stubIDprefit_;
    std::map<int, int> stubID_;
    std::vector<L1TStub*> l1stub_;
    
    unsigned int nstubs_;
    int seed_;
    bool duplicate_;
    int sector_;
  };

};
#endif
