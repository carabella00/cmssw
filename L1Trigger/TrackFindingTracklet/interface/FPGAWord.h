#ifndef L1Trigger_TrackFindingTracklet_interface_FPGAWord_h
#define L1Trigger_TrackFindingTracklet_interface_FPGAWord_h

#include <iostream>
#include <fstream>
#include <sstream>
#include <bitset>
#include <cassert>
#include <cmath>

#include "FWCore/MessageLogger/interface/MessageLogger.h"

using namespace std;

class FPGAWord {
public:
  FPGAWord() {
    value_ = -1;
    nbits_ = -1;
  }

  FPGAWord(int value, int nbits, bool positive = true, int line = -1, const char* file = 0) {
    set(value, nbits, positive, line, file);
  }

  void set(int value, int nbits, bool positive = true, int line = -1, const char* file = 0) {
    value_ = value;
    nbits_ = nbits;
    positive_ = positive;
    if (positive) {
      if (value < 0) {
        edm::LogProblem("Tracklet") << "FPGAWord got negative value:" << value << " (" << file << ":" << line << ")";
      }
      assert(value >= 0);
    }
    if (nbits >= 22) {
      edm::LogPrint("Tracklet") << "FPGAWord got too many bits:" << nbits << " (" << file << ":" << line << ")";
    }
    assert(nbits < 22);
    if (nbits <= 0) {
      edm::LogPrint("Tracklet") << "FPGAWord got too few bits:" << nbits << " (" << file << ":" << line << ")";
    }
    assert(nbits > 0);
    if (positive) {
      if (value >= (1 << nbits)) {
        if (file != 0) {
          edm::LogProblem("Tracklet") << "value too large:" << value << " " << (1 << nbits) << " (" << file << ":" << line << ")";
        }
      }
      assert(value < (1 << nbits));
    } else {
      if (value > (1 << (nbits - 1))) {
        edm::LogProblem("Tracklet") << "value too large:" << value << " " << (1 << (nbits - 1)) << " (" << file << ":" << line << ")";
      }
      assert(value <= (1 << (nbits - 1)));
      if (value < -(1 << (nbits - 1))) {
        edm::LogProblem("Tracklet") << "value too negative:" << value << " " << -(1 << (nbits - 1)) << " (" << file << ":" << line << ")";
      }
      assert(value >= -(1 << (nbits - 1)));
    }
  }

  ~FPGAWord() {}

  std::string str() const {
    const int nbit = nbits_;

    if (!(nbit > 0 && nbit < 22))
      edm::LogVerbatim("Tracklet") << "nbit: " << nbit;
    if (nbit == -1)
      return "?";
    if (nbit == 0)
      return "~";

    assert(nbit > 0 && nbit < 22);

    std::ostringstream oss;
    if (nbit == 1)
      oss << (bitset<1>)value_;
    if (nbit == 2)
      oss << (bitset<2>)value_;
    if (nbit == 3)
      oss << (bitset<3>)value_;
    if (nbit == 4)
      oss << (bitset<4>)value_;
    if (nbit == 5)
      oss << (bitset<5>)value_;
    if (nbit == 6)
      oss << (bitset<6>)value_;
    if (nbit == 7)
      oss << (bitset<7>)value_;
    if (nbit == 8)
      oss << (bitset<8>)value_;
    if (nbit == 9)
      oss << (bitset<9>)value_;
    if (nbit == 10)
      oss << (bitset<10>)value_;
    if (nbit == 11)
      oss << (bitset<11>)value_;
    if (nbit == 12)
      oss << (bitset<12>)value_;
    if (nbit == 13)
      oss << (bitset<13>)value_;
    if (nbit == 14)
      oss << (bitset<14>)value_;
    if (nbit == 15)
      oss << (bitset<15>)value_;
    if (nbit == 16)
      oss << (bitset<16>)value_;
    if (nbit == 17)
      oss << (bitset<17>)value_;
    if (nbit == 18)
      oss << (bitset<18>)value_;
    if (nbit == 19)
      oss << (bitset<19>)value_;
    if (nbit == 20)
      oss << (bitset<20>)value_;
    if (nbit == 21)
      oss << (bitset<21>)value_;
    if (nbit == 22)
      oss << (bitset<22>)value_;

    return oss.str();
  }

  //return the nbits starting with the lsb. lsb=0 means the least significant bit
  unsigned int bits(unsigned int lsb, unsigned int nbit) {
    assert(lsb + nbit <= (unsigned int)nbits());
    return (value_ >> lsb) & ((1 << nbit) - 1);
  }

  int value() const { return value_; }
  int nbits() const { return nbits_; }  //FEXME should be unsigned

  bool atExtreme() const {
    if (positive_)
      return (value_ == 0) || (value_ == (1 << nbits_) - 1);
    //return ((value_==((1<<nbits_)-1))|(value_==((1<<(nbits_-1))-1)));
    return ((value_ == (-(1 << (nbits_ - 1)))) || (value_ == ((1 << (nbits_ - 1)) - 1)));
  }

  bool operator==(const FPGAWord& other) const {
    return (value_ == other.value_) && (nbits_ == other.nbits_) && (positive_ == other.positive_);
  }

private:
  int value_;
  int nbits_;  //FIXME should be unsigned
  bool positive_;
};

#endif
