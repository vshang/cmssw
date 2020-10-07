/**
 * This defines the KF matrices and the operations performance on them.
 *
 *  All variable names & equations come from Fruhwirth KF paper
 * http://dx.doi.org/10.1016/0168-9002%2887%2990887-4
 *
 * Author: Ian Tomalin
 */
 
#ifndef __KalmanMatricesHLS5__
#define __KalmanMatricesHLS5__

// Defines StateHLS & KFstateHLS. Also defines finite bit integers & floats.
#ifdef CMSSW_GIT_HASH
#include "L1Trigger/TrackFindingTMTT/interface/HLS/HLSutilities.h"
#include "L1Trigger/TrackFindingTMTT/interface/HLS/StubHLS.h"
#include "L1Trigger/TrackFindingTMTT/interface/HLS/KFstateHLS.h"
#include "L1Trigger/TrackFindingTMTT/interface/HLS/HLSconstants.h"
#include "L1Trigger/TrackFindingTMTT/interface/HLS/KalmanMatricesHLS.h"
#include "L1Trigger/TrackFindingTMTT/interface/HLS/KalmanMatricesHLS4.h"
#else
#include "HLSutilities.h"
#include "StubHLS.h"
#include "KFstateHLS.h"
#include "HLSconstants.h"
#include "KalmanMatricesHLS.h"
#include "KalmanMatricesHLS5.h"
#endif
 
#ifdef CMSSW_GIT_HASH
namespace TMTT {

namespace KalmanHLS {
#endif

// Utility for -1/r from ROM used by MatrixH<5>. 

class MinusOneOverR {
public:
  enum {BRED = 4}; // To save ROM resources, reduce granularity in r by this number of bits.
  enum {MAXN = 1 << (BSR - BRED)}; // pow(2,BSR) // Max. value of [r / pow(2,BRED)].
  enum {BIR = -8}; // Chosen using CheckCalc output.
  typedef AP_FIXED(BSR+1,BIR) TIR;
public:

  MinusOneOverR() {
    for (unsigned int n = 2; n < MAXN; n++) { // Don't bother initializing first two, as would overflow bit range.
      get[n] = -1./(float((n << BRED)) + 0.5*(1 << BRED)); // Round to nearest half integer
    }
  }
public:
  TIR get[MAXN];
};


// Calculate matrix of derivatives of predicted stub coords w.r.t. helix params.

template <>
class MatrixH<5> {
public:
  enum {BH=BSR+1};
  enum {BHD=MinusOneOverR::BIR};
  typedef AP_FIXED(BH,BH)  TH;  // One extra bit, since "-r" can be -ve.
  typedef AP_FIXED(BH,BHD) THD; // For d0 elements (which have size -1/r)
  typedef AP_UFIXED(1,1)   T1;
  MatrixH(const StubHLS::TR& r) : _00(-r),         _12(r),         _04(this->setH04(r)),
                                           _01(1), _02(0), _03(0),
                                  _10(0),  _11(0),         _13(1), _14(0) {
#ifdef PRINT_SUMMARY
    CHECK_AP::checkCalc("H04", _04, -1./double(r), 0.03);
#endif
  }

  // Set element _04 to -1/r.
  static THD setH04(const StubHLS::TR& r);

public:
  TH       _00, _12;
  THD                          _04;
  const T1      _01, _02, _03, 
           _10, _11,      _13, _14;
};

// S = H * C

template <>
class MatrixS<5> {
public:
  enum {BH=MatrixH<5>::BH, BHD=MatrixH<5>::BHD,
  // Calculate number of integer bits required for all non-zero elements of S.
  // (Assumes that some elements of C & H are zero and that all non-trivial elements of H have BH or BHD bits).
        BS00=MAX3(BH+BHC00, BHC01, BHD+BHC04) - BODGE<5>::S,  // H00*C00 + H01*C10 + (H02*C20 + H03*C30 = zero) + H04*C40.
        BS01=MAX3(BH+BHC01, BHC11, BHD+BHC14) - BODGE<5>::S,  // H00*C01 + H01*C11 + (H02*C21 + H03*C31 = zero) + H04*C41.
        BS12=MAX2(BH+BHC22, BHC23)            - BODGE<5>::S,  // (H10*C02 + H11*C12 = zero) + H12*C22 + H13*C32 + (H14*C42 = zero).
	BS13=MAX2(BH+BHC23, BHC33)            - BODGE<5>::S,  // (H10*C03 + H11*C13 = zero) + H12*C23 + H13*C33 + (H14*C43 = zero);
        BS04=MAX3(BH+BHC04, BHC14, BHD+BHC44) - BODGE<5>::S  // H00*C04 + H01*C14 + (H02*C24 + H03*C34 = zero) + H04*C44.
       };
  typedef AP_FIXED(B27,BS00)  TS00;
  typedef AP_FIXED(B27,BS01)  TS01;
  typedef AP_FIXED(B27,BS12)  TS12;
  typedef AP_FIXED(B27,BS13)  TS13;
  typedef AP_FIXED(B27,BS04)  TS04;
  typedef AP_FIXED(BCORR,0)   T0;     // Neglect correlation between r-phi & r-z planes for now.

public:

  MatrixS(const MatrixH<5>& H, const MatrixC<5>& C);

public:
  
  TS00 _00;
  TS01 _01;
  TS12 _12;
  TS13 _13;
  TS04 _04;
  T0          _02, _03,
    _10, _11,           _14;
};

// Covariance matrix of helix params.

template <>
class MatrixC<5> {
public:
  typedef AP_UFIXED(1,1)      T0; // HLS doesn't like zero bit variables.

  // Determine input helix coviaraiance matrix.
  MatrixC(const KFstateHLS<5>& stateIn) :
             _00(stateIn.cov_00), _11(stateIn.cov_11), _22(stateIn.cov_22), _33(stateIn.cov_33), _44(stateIn.cov_44),
	     _01(stateIn.cov_01), _23(stateIn.cov_23), _04(stateIn.cov_04), _14(stateIn.cov_14),
             _02(0), _03(0), _12(0), _13(0), _42(0), _43(0),
             _10(_01), _32(_23), _40(_04), _41(_14), _20(_02), _30(_03), _21(_12), _31(_13), _24(_42), _34(_43) {}

  // Calculate output helix covariance matrix: C' = C - K*H*C = C - K*S.
  MatrixC(const MatrixC<5>& C, const MatrixK<5>& K, const MatrixS<5>& S);

public:
  // Elements that are finite
  // VHDL interface wierdly uses signed 25 bits for these, but makes more sense to use unsigned 24 instead.
  AP_UFIXED(B24,BHC00-1) _00; // One less integer bit as no sign required.
  AP_UFIXED(B24,BHC11-1) _11;
  AP_UFIXED(B24,BHC22-1) _22;
  AP_UFIXED(B24,BHC33-1) _33;
  AP_UFIXED(B24,BHC44-1) _44;
  KFstateHLS<5>::TC01 _01; // (inv2R, phi0) -- other off-diagonal elements assumed negligeable.
  KFstateHLS<5>::TC23 _23; // (tanL,  z0)   -- other off-diagonal elements assumed negligeable.
  KFstateHLS<5>::TC04 _04; // (inv2R, d0)   -- other off-diagonal elements assumed negligeable.
  KFstateHLS<5>::TC14 _14; // (phi0,  d0)   -- other off-diagonal elements assumed negligeable.
  // Elements that are zero.
  const T0 _02, _03, _12, _13, _42, _43;
  // Elements below the diagonal of this symmetric matrix.
  const KFstateHLS<5>::TC01 &_10;
  const KFstateHLS<5>::TC23 &_32;
  const KFstateHLS<5>::TC04 &_40;
  const KFstateHLS<5>::TC14 &_41;
  const T0 &_20, &_30, &_21, &_31, &_24, &_34;
};

// S(transpose) = C*H(transpose)

template <>
class MatrixS_transpose<5> {
public:
  typedef MatrixS<5>::TS00  TS00;
  typedef MatrixS<5>::TS01  TS01;
  typedef MatrixS<5>::TS12  TS12;
  typedef MatrixS<5>::TS13  TS13;
  typedef MatrixS<5>::TS04  TS04;
  typedef MatrixS<5>::T0    T0;
  MatrixS_transpose(const MatrixS<5>& S) : _00(S._00), _10(S._01), _21(S._12), _31(S._13), _40(S._04),
					   _01(S._10), _11(S._11), _20(S._02), _30(S._03), _41(S._14) {}
public:
  const TS00&  _00;
  const TS01&  _10;
  const TS12&  _21;
  const TS13&  _31;
  const TS04&  _40;
  const T0&       _01,
                  _11,
             _20,
             _30,
                  _41;
};

// Covariance matrix of predicted residuals R = V + H*C*Ht = V + H*St.

template <>
class MatrixR<5> {
public:
  enum {BH=MatrixH<5>::BH, BHD=MatrixH<5>::BHD,
	BS00=MatrixS<5>::BS00,
	BS01=MatrixS<5>::BS01,
	BS12=MatrixS<5>::BS12,
	BS13=MatrixS<5>::BS13,
	BS04=MatrixS<5>::BS04,
        // Calculate number of integer bits required for elements of R.
        BR00 = MAX2(MatrixV::BVPP, MAX3(BH+BS00, BS01, BHD+BS04)) - BODGE<5>::R, // H00*St00 + H01*St10 + (H02*St20 + H03*St30 = zero) + H04*St40
	BR11 = MAX2(MatrixV::BVZZ, MAX2(BH+BS12, BS13))           - BODGE<5>::R, // (H10*St01 + H11*St11 = zero) + H12*St21 + H13*St31 + (H14*St41 = zero)
	BR01 = 0                                         // (H00*St01 + H01*St11 + H02*St21 + H03*St31 + H04*St41 = zero)
       };  
  typedef SW_UFIXED(B34,BR00)   TR00;
  typedef SW_UFIXED(B34,BR11)   TR11;
  typedef SW_UFIXED(BCORR,BR01) TR01;

public:
  MatrixR(const MatrixV& V, const MatrixH<5>& H, const MatrixS_transpose<5>& St);

public:
  TR00  _00;
  TR11  _11;
  TR01  _01;
  TR01& _10; // Matrix symmetric, so can use reference.
};

// Kalman gain matrix K = S(transpose)*R(inverse).

template <>
class MatrixK<5> {
public:
  enum {BS00=MatrixS<5>::BS00,
	BS01=MatrixS<5>::BS01,
	BS12=MatrixS<5>::BS12,
	BS13=MatrixS<5>::BS13,
	BS04=MatrixS<5>::BS04,
        BIR00=MatrixInverseR<5>::BIR00,
        BIR11=MatrixInverseR<5>::BIR11,
        BIR01=MatrixInverseR<5>::BIR01,
        BK00=(BS00+BIR00) - BODGE<5>::K,  // St00*Rinv00 (+ St01*Rinv10 = zero)
        BK10=(BS01+BIR00) - BODGE<5>::K,  // St10*Rinv00 (+ St11*Rinv10 = zero)
        BK21=(BS12+BIR11) - BODGE<5>::K,  // (St20*Rinv01 = zero) + St21*Rinv11
        BK31=(BS13+BIR11) - BODGE<5>::K,  // (St30*Rinv01 = zero) + St31*Rinv11
        BK40=(BS04+BIR00) - BODGE<5>::K}; // St40*Rinv00 (+ St41*Rinv10 = zero)
  typedef SW_FIXED(B35,BK00)  TK00;
  typedef SW_FIXED(B35,BK10)  TK10;
  typedef SW_FIXED(B35,BK21)  TK21;
  typedef SW_FIXED(B35,BK31)  TK31;
  typedef SW_FIXED(B35,BK40)  TK40;
  typedef SW_FIXED(BCORR,0)   T0; // Neglect correlation between r-phi & r-z
  MatrixK(const MatrixS_transpose<5>& St, const MatrixInverseR<5>& RmatInv);
public:
  // Additional types used to cast this matrix to a lower precision one for updated helix param calculation.
  typedef SW_FIXED(B27,BK00)  TK00_short;
  typedef SW_FIXED(B27,BK10)  TK10_short;
  typedef SW_FIXED(B27,BK21)  TK21_short;
  typedef SW_FIXED(B27,BK31)  TK31_short;
  typedef SW_FIXED(B27,BK40)  TK40_short;
public:
  TK00  _00;
  TK10  _10;
  TK21       _21;
  TK31       _31;
  TK40  _40;
  T0         _01,
             _11,
        _20,
        _30,
             _41;
};

// Hit residuals: res = m - H*x. 

template <>
class VectorRes<5> {
public:
  // Use higher granularity for residuals than for stubs.
  // BODGE<5>::RES should be slightly larger than BODGE_V as hits can be several sigma wrong.
  // Add one extra fractional bit relative to stub, to avoid additional precision loss.
  typedef AP_FIXED(B18-BODGE<5>::RES+1,BSP-BODGE<5>::RES) TRP;
  typedef AP_FIXED(B18-BODGE<5>::RES+1,BSZ-BODGE<5>::RES) TRZ;

public:
  VectorRes(const VectorM& m, const MatrixH<5>& H, const VectorX<5>& x);

public:
  TRP _0;
  TRZ _1;
};

// Vector of helix params.

template <>
class VectorX<5> {
public:
  // Determine input helix params.
  VectorX(const KFstateHLS<5>& stateIn) : _0(stateIn.inv2R), _1(stateIn.phi0), _2(stateIn.tanL), _3(stateIn.z0), _4(stateIn.d0) {} 

  // Calculate output helix params: x' = x + K*res
  VectorX(const VectorX<5>& x, const MatrixK<5>& K, const VectorRes<5>& res);

public:
  KFstateHLS<5>::TR _0;
  KFstateHLS<5>::TP _1;  
  KFstateHLS<5>::TT _2;  
  KFstateHLS<5>::TZ _3;  
  KFstateHLS<5>::TD _4;  
};

#ifdef CMSSW_GIT_HASH
}

}
#endif

#endif
