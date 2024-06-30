//*************************************************************************
//* ====================
//*  EXKalState Class
//* ====================
//*
//* (Description)
//*   Sample state vector class used in Kalman Filter.
//* (Requires)
//*     TVKalState
//* (Provides)
//*     class EXKalState
//* (Update Recored)
//*   2003/09/30  Y.Nakashima       Original version.
//*
//*************************************************************************
//
#include <iostream>
#include "EXKalState.h"
#include "EXKalSite.h"

using namespace std;

//_____________________________________________________________________
//  -----------------------------------
//  Base Class for Kalman state vector
//  -----------------------------------
//
// --------------------------------------------------------------------
// Ctors and Dtor
//
EXKalState::EXKalState(Int_t p) 
    : TVKalState(p)
{
}

EXKalState::EXKalState(const TKalMatrix &sv, Int_t type, Int_t p) 
    : TVKalState(sv,type,p)
{
}

EXKalState::EXKalState(const TKalMatrix &sv, const TKalMatrix &c,
                       Int_t type, Int_t p) 
    : TVKalState(sv,c,type,p)
{
}

EXKalState::EXKalState(const TKalMatrix &sv, const TVKalSite &site,
                       Int_t type, Int_t p) 
    : TVKalState(sv,site,type,p)
{
    //std::printf("========Call EXKalState 3  A Type=%d\n",type);
    //kPredicted == 0;
    //kFiltered == 1;
    //kSmoothed ==2;
    //kInvFiltered ==3;
}

EXKalState::EXKalState(const TKalMatrix &sv, const TKalMatrix &c,
                       const TVKalSite &site, Int_t type, Int_t p) 
    : TVKalState(sv,c,site,type,p)
{
}

//
// --------------------------------------------------------------------
// Implementation of base-class pure virtuals
//

EXKalState * EXKalState::MoveTo(TVKalSite  &to,
                                TKalMatrix &F,
                                TKalMatrix *QPtr) const
{
    //std::printf("================Call MoveTo EXKalState2");
    //std::printf("\n");
    Int_t p = GetDimension();
    //std::cout<<"====> "<< this->GetNrows()<<std::endl;
    //std::printf("================Call MoveTo EXKalState P=%d\n",p);
    for (Int_t i=0; i<p; i++) 
        F(i,i) = 1.;

    if (QPtr) {
        QPtr->Zero();
        return (new EXKalState(*this,to,TVKalSite::kPredicted,p));//返回预测状态量，由于状态转移矩阵为单位阵，新的预测量就是现在的状态量
    } else {
        return 0;
    }
}

EXKalState & EXKalState::MoveTo(TVKalSite  &to,
                                TKalMatrix &F,
                                TKalMatrix &Q) const
{
    //std::printf("================Call MoveTo EXKalState1");
    return *MoveTo(to, F, &Q);
}

void EXKalState::DebugPrint() const
{
    cout << " A= " << (*this)(0,0) << endl
        << "  B= " << (*this)(1,0) << endl
        << "  C= " << (*this)(2,0) << endl
        << "  D= " << (*this)(3,0) << endl;
    GetCovMat().DebugPrint(" cov Mat   = ");
}

//
// --------------------------------------------------------------------
// Derived class methods
//

TKalMatrix EXKalState::CalcXAt(Double_t t)
{
    TKalMatrix m(2,1);
    //system equation:
    //y = A*x+B; 
    //z = C*x+D;
    m(0,0) = (*this)(0,0)*t + (*this)(1,0) ;
    m(1,0) = (*this)(2,0)*t + (*this)(3,0) ;
    return m;
}

TKalMatrix EXKalState::CalcXDerivAt(Double_t t)
{
    TKalMatrix H(2,4);
    H(0,0) = t;
    H(0,1) = 1;
    H(0,2) = 0.;
    H(0,3) = 0.;
    H(1,0) = 0.;
    H(1,1) = 0.;
    H(1,2) = t;
    H(1,3) = 1;
    return H;
}
