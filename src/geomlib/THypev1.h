#ifndef THYPEV1_H
#define THYPEV1_H 1
//*************************************************************************
//* ====================
//*  THypev1 Class
//* ====================
//*
//* (Description)
//*   A class to implement a hyperboloidal surface object.
//* (Requires)
//*     TVSurface
//* (Provides)
//*     class THype
//* (Update Recored)
//*   2003/10/03  K.Fujii       Original version.
//*   2005/02/23  K.Fujii       Added GetSortingPolicy().
//*
//*************************************************************************
//
#include "TVSurface.h"
#include "TVector3.h"
#include "TMatrixD.h"

class TVTrack;

//_____________________________________________________________________
//  -----------------------------------
//  Hype Class
//  -----------------------------------

class THypev1 : public TVSurface {
public:
    THypev1(Double_t r = 1., Double_t hlen  = 1., Double_t tana = 0.,
            Double_t xc = 0., Double_t yc = 0, Double_t zc = 0.)
        : fR0(r), fHalfLen(hlen), fXc(xc,yc,zc), fTanA(tana) {}

    virtual ~THypev1() {}

    virtual Double_t CalcS   (const TVector3 &xx) const;
    virtual TMatrixD CalcDSDx(const TVector3 &xx) const;

    inline virtual       Bool_t     IsOnSurface(const TVector3 &xx) const;
    inline virtual       Bool_t     IsOutside  (const TVector3 &xx) const;

    inline virtual       Double_t   GetSortingPolicy()              const;

    inline virtual       Double_t   GetR0     () const { return fR0;   } 
    inline virtual const TVector3 & GetXc     () const { return fXc;   } 
    inline virtual       Double_t   GetTanA   () const { return fTanA; }
    inline virtual       Double_t   GetLength () const;
    inline virtual       Double_t   GetZmin   () const;
    inline virtual       Double_t   GetZmax   () const;

private:
    Double_t fR0;          // radius at z = 0.
    Double_t fHalfLen;     // half length
    TVector3 fXc;          // center
    Double_t fTanA;        // tan(stereo angle)
#if __GNUC__ < 4 && !defined(__STRICT_ANSI__)
    const Double_t kTol = 1.e-5; // tolerance
#else
    static const Double_t kTol; // tolerance
#endif

    ClassDef(THypev1,1)      // hype class
};

//=======================================================
// inline functions
//=======================================================

Double_t THypev1::GetLength () const 
{ 
    return 2*fHalfLen;
} 

Double_t THypev1::GetZmin() const
{
    return fXc.Z() - fHalfLen;
}

Double_t THypev1::GetZmax() const
{
    return fXc.Z() + fHalfLen;
}

Bool_t THypev1::IsOnSurface(const TVector3 &xx) const
{
    TVector3 xxc = xx - fXc;
    Double_t r   = xxc.Perp();
    Double_t z   = xxc.Z();
    Double_t s   = r*r - fTanA*fTanA*z*z - fR0*fR0;

    return (TMath::Abs(s) < kTol && xx.Z() >= GetZmin() && xx.Z() <= GetZmax());
} 

Bool_t THypev1::IsOutside(const TVector3 &xx) const
{
    Double_t r  = (xx-fXc).Perp();
    Double_t z  = xx.Z();
    Double_t R2 = fR0*fR0 + fTanA*fTanA*z*z;

    return (r*r > R2 || z < GetZmin() || z > GetZmax());
} 

Double_t THypev1::GetSortingPolicy() const
{
    return GetR0();
}

#endif

