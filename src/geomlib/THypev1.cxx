//*************************************************************************
//* ====================
//*  THype Class
//* ====================
//*
//* (Description)
//*   A class to implement a hyperboloidal surface object.
//* (Requires)
//*     TVSurface
//* (Provides)
//*     class THype
//* (Update Recored)
//*   2003/10/03  K.Fujii       Original version. Currently fXc is
//*                             supposed to be at the origin.
//*************************************************************************
//
#include <iostream>
#include "TCircle.h"
#include "THypev1.h"
#include "TVTrack.h"

using namespace std;

#if __GNUC__ < 4 && !defined(__STRICT_ANSI__)
#else
const Double_t THypev1::kTol = 1.e-5; // tolerance
#endif

//_____________________________________________________________________
//  -----------------------------------
//  THype Class
//  -----------------------------------

ClassImp(THypev1)

//_____________________________________________________________________
//  -----------------------------------
//  Calculate S
//  -----------------------------------
//
Double_t THypev1::CalcS(const TVector3 &xx) const
{
   TVector3 xxc = xx - fXc;
   Double_t s   = xxc.X() * xxc.X() + xxc.Y() * xxc.Y()
                   - fR0 * fR0 - xxc.Z() * xxc.Z() * fTanA * fTanA;
   return s;
}

//_____________________________________________________________________
//  -----------------------------------
//  Calculate (@S/@x)
//  -----------------------------------
//
TMatrixD THypev1::CalcDSDx(const TVector3 &xx) const
{
   TVector3 xxc = xx - fXc;
   TMatrixD dsdx(1,3);
   dsdx(0,0) =  2.* xxc.X();
   dsdx(0,1) =  2.* xxc.Y();
   dsdx(0,2) = -2.* xxc.Z() * fTanA * fTanA;
   return dsdx;
}
