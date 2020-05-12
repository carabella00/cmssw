#ifndef L1Trigger_TrackFindingTracklet_interface_IMATH_TrackletCalculatorDisk_h
#define L1Trigger_TrackFindingTracklet_interface_IMATH_TrackletCalculatorDisk_h

#include "Settings.h"
#include "imath.h"

//
// Constants used:
//   dphisector
//   rmaxL6
//   zmaxD5
//   rmaxdisk
//   kr, kphi1, kz
//
//   rmean[], zmean[]

class IMATH_TrackletCalculatorDisk {
public:
  IMATH_TrackletCalculatorDisk(const Trklet::Settings* settings, imathGlobals* globals, int i1, int i2)
      : settings_(settings), globals_(globals) {
#ifndef CMSSW_GIT_HASH
    edm::LogVerbatim("Tracklet") << "=============================================";
    char s[1024];
    sprintf(s, "IMATH Tracklet Calculator for Disk %i %i dphisector = %f", i1, i2, settings->dphisector());
    edm::LogVerbatim("Tracklet") << s;
    sprintf(s, "rmaxL6 = %f, zmaxD5 = %f", settings->rmax(5), settings->zmax(4));
    edm::LogVerbatim("Tracklet") << s;
    sprintf(s, "      stub Ks: kr, kphi1, kz = %g, %g, %g", settings->kr(), settings->kphi1(), settings->kz());
    edm::LogVerbatim("Tracklet") << s;
    sprintf(s,
            "  tracklet Ks: krinvpars, kphi0pars, ktpars, kzpars = %g, %g, %g, %g",
            settings->kphi1() / settings->kr() * pow(2, settings->rinv_shift()),
            settings->kphi1() * pow(2, settings->phi0_shift()),
            settings->kz() / settings->kr() * pow(2, settings->t_shift()),
            settings->kz() * pow(2, settings->z0_shift()));
    edm::LogVerbatim("Tracklet") << s;
    sprintf(s,
            "layer proj Ks: kphiproj456, kphider, kzproj, kzder = %g, %g, %g, %g",
            settings->kphi1() * pow(2, settings->SS_phiL_shift()),
            settings->kphi1() / settings->kr() * pow(2, settings->SS_phiderL_shift()),
            settings->kz() * pow(2, settings->PS_zL_shift()),
            settings->kz() / settings->kr() * pow(2, settings->PS_zderL_shift()));
    edm::LogVerbatim("Tracklet") << s;
    sprintf(s,
            " disk proj Ks: kphiprojdisk, kphiprojderdisk, krprojdisk, krprojderdisk = %g, %g, %g, %g",
            settings->kphi1() * pow(2, settings->SS_phiD_shift()),
            settings->kphi1() / settings->kr() * pow(2, settings->SS_phiderD_shift()),
            settings->kr() * pow(2, settings->PS_rD_shift()),
            settings->kr() / settings->kz() * pow(2, settings->PS_rderD_shift()));
    edm::LogVerbatim("Tracklet") << s;
    edm::LogVerbatim("Tracklet") << "=============================================";
#endif

    z1mean.set_fval(settings->zmean(abs(i1) - 1));
    z2mean.set_fval(settings->zmean(abs(i2) - 1));

    if (i2 < 0) {  //t is negative
      z1mean.set_fval(-settings->zmean(abs(i1) - 1));
      z2mean.set_fval(-settings->zmean(abs(i2) - 1));
      invt.set_mode(var_inv::mode::neg);
      invt.initLUT(0.);
    }

    valid_phiL_0.add_cut(&t_layer_cut);
    valid_phiL_1.add_cut(&t_layer_cut);
    valid_phiL_2.add_cut(&t_layer_cut);

    valid_der_phiL.add_cut(&t_layer_cut);

    valid_zL_0.add_cut(&t_layer_cut);
    valid_zL_1.add_cut(&t_layer_cut);
    valid_zL_2.add_cut(&t_layer_cut);

    valid_der_zL.add_cut(&t_layer_cut);

    valid_phiD_0.add_cut(&t_disk_cut_left);
    valid_phiD_1.add_cut(&t_disk_cut_left);
    valid_phiD_2.add_cut(&t_disk_cut_left);

    valid_der_phiD.add_cut(&t_disk_cut_left);

    valid_rD_0.add_cut(&t_disk_cut_left);
    valid_rD_1.add_cut(&t_disk_cut_left);
    valid_rD_2.add_cut(&t_disk_cut_left);

    valid_der_rD.add_cut(&t_disk_cut_left);

    valid_phiD_0.add_cut(&t_disk_cut_right);
    valid_phiD_1.add_cut(&t_disk_cut_right);
    valid_phiD_2.add_cut(&t_disk_cut_right);

    valid_der_phiD.add_cut(&t_disk_cut_right);

    valid_rD_0.add_cut(&t_disk_cut_right);
    valid_rD_1.add_cut(&t_disk_cut_right);
    valid_rD_2.add_cut(&t_disk_cut_right);

    valid_der_rD.add_cut(&t_disk_cut_right);
  }

  const Trklet::Settings* settings_;
  imathGlobals* globals_;

  //max values
  double dr_max = 20.;
  double delta0_max = 0.005;
  double a2a_max = 0.1;
  double x8_max = 1.;
  double x22_max = 0.3;
  double x13_max = 300.;
  double deltaZ_max = 8.;
  double der_phiD_max = 0.002;

  // constants
  //
  var_param plus2{globals_, "plus2", 2., 10};
  var_param plus1{globals_, "plus1", 1., 10};
  var_param minus1{globals_, "minus1", -1, 10};
  //
  //
  var_param z1mean{globals_, "z1mean", "Kz", settings_->zmax(4), settings_->kz()};
  var_param z2mean{globals_, "z2mean", "Kz", settings_->zmax(4), settings_->kz()};

  //inputs
  var_def r1{globals_, "r1", "Kr", settings_->rmax(5), settings_->kr()};
  var_def r2{globals_, "r2", "Kr", settings_->rmax(5), settings_->kr()};
  var_def z1{globals_, "z1", "Kz", settings_->dzmax(), settings_->kz()};
  var_def z2{globals_, "z2", "Kz", settings_->dzmax(), settings_->kz()};

  var_def phi1{globals_, "phi1", "Kphi", settings_->dphisector() / 0.75, settings_->kphi1()};
  var_def phi2{globals_, "phi2", "Kphi", settings_->dphisector() / 0.75, settings_->kphi1()};

  var_def rproj0{globals_, "rproj0", "Kr", settings_->rmax(5), settings_->kr()};
  var_def rproj1{globals_, "rproj1", "Kr", settings_->rmax(5), settings_->kr()};
  var_def rproj2{globals_, "rproj2", "Kr", settings_->rmax(5), settings_->kr()};

  var_def zproj0{globals_, "zproj0", "Kz", settings_->zmax(4), settings_->kz()};
  var_def zproj1{globals_, "zproj1", "Kz", settings_->zmax(4), settings_->kz()};
  var_def zproj2{globals_, "zproj2", "Kz", settings_->zmax(4), settings_->kz()};

  //calculations

  //tracklet
  var_add z1abs{globals_, "z1abs", &z1, &z1mean, settings_->zmax(4)};
  var_add z2abs{globals_, "z2abs", &z2, &z2mean, settings_->zmax(4)};

  var_subtract dr{globals_, "dr", &r2, &r1, dr_max};

  //R LUT
  var_inv drinv{globals_, "drinv", &dr, 0, 18, 23, 0, var_inv::mode::pos};

  var_subtract dphi{globals_, "dphi", &phi2, &phi1, settings_->dphisector() / 4.};
  var_subtract dz{globals_, "dz", &z2abs, &z1abs, 50.};

  var_mult delta0{globals_, "delta0", &dphi, &drinv, 4 * delta0_max};
  var_mult deltaZ{globals_, "deltaZ", &dz, &drinv, deltaZ_max};
  var_mult delta1{globals_, "delta1", &r1, &delta0};
  var_mult delta2{globals_, "delta2", &r2, &delta0};
  var_mult a2a{globals_, "a2a", &delta1, &delta2, 4 * a2a_max};
  var_nounits a2b{globals_, "a2b", &a2a};
  var_subtract a2{globals_, "a2", &plus2, &a2b, 3.};
  var_neg a2n{globals_, "a2n", &a2};
  var_shift a{globals_, "a", &a2, 1};

  var_add Rabs{globals_, "Rabs", &r1, &r2};
  var_timesC R6{globals_, "R6", &Rabs, 1. / 6., 12};

  var_mult x4{globals_, "x4", &R6, &delta0};
  var_mult x6a{globals_, "x6a", &delta2, &x4, 0.16};
  var_nounits x6b{globals_, "x6b", &x6a};
  var_add x6m{globals_, "x6m", &minus1, &x6b, 2.};
  var_mult phi0a{globals_, "phi0a", &delta1, &x6m, settings_->dphisector()};

  var_mult z0a{globals_, "z0a", &r1, &deltaZ, 240.};
  var_mult z0b{globals_, "z0b", &z0a, &x6m, 240.};

  var_add phi0{globals_, "phi0", &phi1, &phi0a, 2 * settings_->dphisector()};
  var_mult rinv{globals_, "rinv", &a2n, &delta0, 4 * settings_->maxrinv()};
  var_mult t{globals_, "t", &a, &deltaZ, 15.8};
  var_add z0{globals_, "z0", &z1abs, &z0b, 40.};

  var_adjustK rinv_final{
      globals_, "rinv_final", &rinv, settings_->kphi1() / settings_->kr() * pow(2, settings_->rinv_shift())};
  var_adjustK phi0_final{globals_, "phi0_final", &phi0, settings_->kphi1() * pow(2, settings_->phi0_shift())};
  var_adjustK t_final{globals_, "t_final", &t, settings_->kz() / settings_->kr() * pow(2, settings_->t_shift())};
  var_adjustK z0_final{globals_, "z0_final", &z0, settings_->kz() * pow(2, settings_->z0_shift())};

  //projection to r
  //
  var_shift x2{globals_, "x2", &delta0, 1};

  var_mult x1_0{globals_, "x1_0", &x2, &rproj0};
  var_mult x1_1{globals_, "x1_1", &x2, &rproj1};
  var_mult x1_2{globals_, "x1_2", &x2, &rproj2};

  var_mult x8_0{globals_, "x8_0", &x1_0, &a2n, x8_max};
  var_mult x8_1{globals_, "x8_1", &x1_1, &a2n, x8_max};
  var_mult x8_2{globals_, "x8_2", &x1_2, &a2n, x8_max};

  var_mult x12_0{globals_, "x12_0", &x8_0, &x8_0};
  var_mult x12_1{globals_, "x12_1", &x8_1, &x8_1};
  var_mult x12_2{globals_, "x12_2", &x8_2, &x8_2};

  var_nounits x12A_0{globals_, "x12A_0", &x12_0};
  var_nounits x12A_1{globals_, "x12A_1", &x12_1};
  var_nounits x12A_2{globals_, "x12A_2", &x12_2};

  var_timesC x20_0{globals_, "x20_0", &x12A_0, 1. / 6.};
  var_timesC x20_1{globals_, "x20_1", &x12A_1, 1. / 6.};
  var_timesC x20_2{globals_, "x20_2", &x12A_2, 1. / 6.};

  var_add x10_0{globals_, "x10_0", &plus1, &x20_0};
  var_add x10_1{globals_, "x10_1", &plus1, &x20_1};
  var_add x10_2{globals_, "x10_2", &plus1, &x20_2};

  var_mult x22_0{globals_, "x22_0", &x8_0, &x10_0, 2 * x22_max};
  var_mult x22_1{globals_, "x22_1", &x8_1, &x10_1, 2 * x22_max};
  var_mult x22_2{globals_, "x22_2", &x8_2, &x10_2, 2 * x22_max};

  var_subtract phiL_0{globals_, "phiL_0", &phi0_final, &x22_0, -1, phi0_final.get_nbits() + 1};
  var_subtract phiL_1{globals_, "phiL_1", &phi0_final, &x22_1, -1, phi0_final.get_nbits() + 1};
  var_subtract phiL_2{globals_, "phiL_2", &phi0_final, &x22_2, -1, phi0_final.get_nbits() + 1};

  var_shift x3{globals_, "x3", &rinv, 1};
  var_neg der_phiL{globals_, "der_phiL", &x3};

  var_adjustK phiL_0_final{globals_, "phiL_0_final", &phiL_0, settings_->kphi1() * pow(2, settings_->SS_phiL_shift())};
  var_adjustK phiL_1_final{globals_, "phiL_1_final", &phiL_1, settings_->kphi1() * pow(2, settings_->SS_phiL_shift())};
  var_adjustK phiL_2_final{globals_, "phiL_2_final", &phiL_2, settings_->kphi1() * pow(2, settings_->SS_phiL_shift())};

  var_adjustK der_phiL_final{globals_,
                             "der_phiL_final",
                             &der_phiL,
                             settings_->kphi1() / settings_->kr() * pow(2, settings_->SS_phiderL_shift())};

  var_mult x11_0{globals_, "x11_0", &rproj0, &t};
  var_mult x11_1{globals_, "x11_1", &rproj1, &t};
  var_mult x11_2{globals_, "x11_2", &rproj2, &t};

  var_mult x23_0{globals_, "x23_0", &x11_0, &x10_0, 800};
  var_mult x23_1{globals_, "x23_1", &x11_1, &x10_1, 800};
  var_mult x23_2{globals_, "x23_2", &x11_2, &x10_2, 800};

  var_add zL_0{globals_, "zL_0", &z0, &x23_0};
  var_add zL_1{globals_, "zL_1", &z0, &x23_1};
  var_add zL_2{globals_, "zL_2", &z0, &x23_2};

  var_adjustK zL_0_final{globals_, "zL_0_final", &zL_0, settings_->kz() * pow(2, settings_->PS_zL_shift())};
  var_adjustK zL_1_final{globals_, "zL_1_final", &zL_1, settings_->kz() * pow(2, settings_->PS_zL_shift())};
  var_adjustK zL_2_final{globals_, "zL_2_final", &zL_2, settings_->kz() * pow(2, settings_->PS_zL_shift())};

  var_adjustK der_zL_final{
      globals_, "der_zL_final", &t_final, settings_->kz() / settings_->kr() * pow(2, settings_->PS_zderL_shift())};

  //projection to z
  //
  var_inv invt{globals_, "invt", &t_final, 0., 18, 26, 1, var_inv::mode::pos, 13};

  var_mult x7{globals_, "x7", &x2, &a2};

  var_subtract x5_0{globals_, "x5_0", &zproj0, &z0};
  var_subtract x5_1{globals_, "x5_1", &zproj1, &z0};
  var_subtract x5_2{globals_, "x5_2", &zproj2, &z0};

  var_mult x13_0{globals_, "x13_0", &x5_0, &invt, x13_max};
  var_mult x13_1{globals_, "x13_1", &x5_1, &invt, x13_max};
  var_mult x13_2{globals_, "x13_2", &x5_2, &invt, x13_max};

  var_mult x25_0{globals_, "x25_0", &x13_0, &x7, settings_->dphisector()};
  var_mult x25_1{globals_, "x25_1", &x13_1, &x7, settings_->dphisector()};
  var_mult x25_2{globals_, "x25_2", &x13_2, &x7, settings_->dphisector()};

  var_add phiD_0{globals_, "phiD_0", &phi0, &x25_0, 2 * settings_->dphisector()};
  var_add phiD_1{globals_, "phiD_1", &phi0, &x25_1, 2 * settings_->dphisector()};
  var_add phiD_2{globals_, "phiD_2", &phi0, &x25_2, 2 * settings_->dphisector()};

  var_adjustK phiD_0_final{globals_, "phiD_0_final", &phiD_0, settings_->kphi1() * pow(2, settings_->SS_phiD_shift())};
  var_adjustK phiD_1_final{globals_, "phiD_1_final", &phiD_1, settings_->kphi1() * pow(2, settings_->SS_phiD_shift())};
  var_adjustK phiD_2_final{globals_, "phiD_2_final", &phiD_2, settings_->kphi1() * pow(2, settings_->SS_phiD_shift())};

  var_mult der_phiD{globals_, "der_phiD", &x7, &invt, 2 * der_phiD_max};

  var_adjustK der_phiD_final{globals_,
                             "der_phiD_final",
                             &der_phiD,
                             settings_->kphi1() / settings_->kr() * pow(2, settings_->SS_phiderD_shift())};

  var_mult x26_0{globals_, "x26_0", &x25_0, &x25_0};
  var_mult x26_1{globals_, "x26_1", &x25_1, &x25_1};
  var_mult x26_2{globals_, "x26_2", &x25_2, &x25_2};

  var_nounits x26A_0{globals_, "x26A_0", &x26_0};
  var_nounits x26A_1{globals_, "x26A_1", &x26_1};
  var_nounits x26A_2{globals_, "x26A_2", &x26_2};

  var_timesC x9_0{globals_, "x9_0", &x26A_0, 1. / 6.};
  var_timesC x9_1{globals_, "x9_1", &x26A_1, 1. / 6.};
  var_timesC x9_2{globals_, "x9_2", &x26A_2, 1. / 6.};

  var_subtract x27_0{globals_, "x27_0", &plus1, &x9_0};
  var_subtract x27_1{globals_, "x27_1", &plus1, &x9_1};
  var_subtract x27_2{globals_, "x27_2", &plus1, &x9_2};

  var_mult rD_0{globals_, "rD_0", &x13_0, &x27_0, settings_->rmaxdisk()};
  var_mult rD_1{globals_, "rD_1", &x13_1, &x27_1, settings_->rmaxdisk()};
  var_mult rD_2{globals_, "rD_2", &x13_2, &x27_2, settings_->rmaxdisk()};

  var_adjustK rD_0_final{globals_, "rD_0_final", &rD_0, settings_->kr() * pow(2, settings_->PS_rD_shift())};
  var_adjustK rD_1_final{globals_, "rD_1_final", &rD_1, settings_->kr() * pow(2, settings_->PS_rD_shift())};
  var_adjustK rD_2_final{globals_, "rD_2_final", &rD_2, settings_->kr() * pow(2, settings_->PS_rD_shift())};

  var_adjustK der_rD_final{
      globals_, "der_rD_final", &invt, settings_->kr() / settings_->kz() * pow(2, settings_->PS_rderD_shift())};

  var_cut rinv_final_cut{globals_, &rinv_final, -settings_->rinvcut(), settings_->rinvcut()};
  var_cut z0_final_cut{globals_, &z0_final, -settings_->z0cut(), settings_->z0cut()};

  var_cut z1abs_cut{globals_, &z1abs, -settings_->zmax(4), settings_->zmax(4)};
  var_cut z2abs_cut{globals_, &z2abs, -settings_->zmax(4), settings_->zmax(4)};
  var_cut dr_cut{globals_, &dr, -dr_max, dr_max};
  var_cut dphi_cut{globals_, &dphi, -settings_->dphisector() / 4., settings_->dphisector() / 4.};
  var_cut dz_cut{globals_, &dz, -50., 50.};
  var_cut delta0_cut{globals_, &delta0, -delta0_max, delta0_max};
  var_cut deltaZ_cut{globals_, &deltaZ, -deltaZ_max, deltaZ_max};
  var_cut a2a_cut{globals_, &a2a, -a2a_max, a2a_max};
  var_cut a2_cut{globals_, &a2, -3., 3.};
  var_cut x6a_cut{globals_, &x6a, -0.02, 0.02};
  var_cut x6m_cut{globals_, &x6m, -2., 2.};
  var_cut phi0a_cut{globals_, &phi0a, -settings_->dphisector(), settings_->dphisector()};
  var_cut z0a_cut{globals_, &z0a, -205, 205};
  var_cut phi0_cut{globals_, &phi0, -2 * settings_->dphisector(), 2 * settings_->dphisector()};
  var_cut rinv_cut{globals_, &rinv, -settings_->maxrinv(), settings_->maxrinv()};
  var_cut t_cut{globals_, &t, -7.9, 7.9};
  var_cut z0_cut{globals_, &z0, -20., 20.};
  var_cut x8_0_cut{globals_, &x8_0, -x8_max, x8_max};
  var_cut x8_1_cut{globals_, &x8_1, -x8_max, x8_max};
  var_cut x8_2_cut{globals_, &x8_2, -x8_max, x8_max};
  var_cut x22_0_cut{globals_, &x22_0, -x22_max, x22_max};
  var_cut x22_1_cut{globals_, &x22_1, -x22_max, x22_max};
  var_cut x22_2_cut{globals_, &x22_2, -x22_max, x22_max};
  var_cut x23_0_cut{globals_, &x23_0, -200, 200};
  var_cut x23_1_cut{globals_, &x23_1, -200, 200};
  var_cut x23_2_cut{globals_, &x23_2, -200, 200};
  var_cut x13_0_cut{globals_, &x13_0, -x13_max, x13_max};
  var_cut x13_1_cut{globals_, &x13_1, -x13_max, x13_max};
  var_cut x13_2_cut{globals_, &x13_2, -x13_max, x13_max};
  var_cut x25_0_cut{globals_, &x25_0, -settings_->dphisector(), settings_->dphisector()};
  var_cut x25_1_cut{globals_, &x25_1, -settings_->dphisector(), settings_->dphisector()};
  var_cut x25_2_cut{globals_, &x25_2, -settings_->dphisector(), settings_->dphisector()};
  var_cut phiD_0_cut{globals_, &phiD_0, -2 * settings_->dphisector(), 2 * settings_->dphisector()};
  var_cut phiD_1_cut{globals_, &phiD_1, -2 * settings_->dphisector(), 2 * settings_->dphisector()};
  var_cut phiD_2_cut{globals_, &phiD_2, -2 * settings_->dphisector(), 2 * settings_->dphisector()};
  var_cut der_phiD_cut{globals_, &der_phiD, -der_phiD_max, der_phiD_max};
  var_cut rD_0_cut{globals_, &rD_0, -settings_->rmaxdisk(), settings_->rmaxdisk()};
  var_cut rD_1_cut{globals_, &rD_1, -settings_->rmaxdisk(), settings_->rmaxdisk()};
  var_cut rD_2_cut{globals_, &rD_2, -settings_->rmaxdisk(), settings_->rmaxdisk()};

  var_cut t_disk_cut_left{globals_, &t, -7.9, -1};
  var_cut t_disk_cut_right{globals_, &t, 1, 7.9};
  var_cut t_layer_cut{globals_, &t, -2.5, 2.5};

  // the following flags are used to apply the cuts in TrackletCalculator
  // and in the output Verilog
  var_flag valid_trackpar{globals_, "valid_trackpar", &rinv_final, &phi0_final, &t_final, &z0_final};

  var_flag valid_phiL_0{globals_, "valid_phiL_0", &phiL_0_final};
  var_flag valid_phiL_1{globals_, "valid_phiL_1", &phiL_1_final};
  var_flag valid_phiL_2{globals_, "valid_phiL_2", &phiL_2_final};

  var_flag valid_zL_0{globals_, "valid_zL_0", &zL_0_final};
  var_flag valid_zL_1{globals_, "valid_zL_1", &zL_1_final};
  var_flag valid_zL_2{globals_, "valid_zL_2", &zL_2_final};

  var_flag valid_der_phiL{globals_, "valid_der_phiL", &der_phiL_final};
  var_flag valid_der_zL{globals_, "valid_der_zL", &der_zL_final};

  var_flag valid_phiD_0{globals_, "valid_phiD_0", &phiD_0_final};
  var_flag valid_phiD_1{globals_, "valid_phiD_1", &phiD_1_final};
  var_flag valid_phiD_2{globals_, "valid_phiD_2", &phiD_2_final};

  var_flag valid_rD_0{globals_, "valid_rD_0", &rD_0_final};
  var_flag valid_rD_1{globals_, "valid_rD_1", &rD_1_final};
  var_flag valid_rD_2{globals_, "valid_rD_2", &rD_2_final};

  var_flag valid_der_phiD{globals_, "valid_der_phiD", &der_phiD_final};
  var_flag valid_der_rD{globals_, "valid_der_rD", &der_rD_final};
};

#endif
