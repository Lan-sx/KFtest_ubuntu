/*********************************************************************
 * Author           : Lan-sx
 * Email            : shexin@ihep.ac.cn
 * Last modified    : 2024-04-29 14:16
 * Filename         : EXKalTest.cxx
 * Description      : Example import from KalTest, remove ClassDef(),ClassImp()
 * Update           : 
 * ******************************************************************/

#include "EXKalTest.h"
#include "EXHit.h"      // class inherit from TMartixD
#include "EXKalState.h"
#include "EXKalSite.h"
#include "EXKalSystem.h"
#include "TPlane.h"
#include "TCylinder.h"
#include "THypev1.h"

#include "TFile.h"
#include "TNtupleD.h"
#include "TMath.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "Lansxlogon.h"

#define __TSTEP__     1.
#define __SLOPE__     0.2
#define __INTERCEPT__ 2.
#define __SIGMAX__    6
//#define __DEBUG__

using namespace std;

int main (Int_t argc, Char_t **argv)
{
    gROOT->SetBatch();
    TApplication app("EXKalTest", &argc, argv, 0, 0);
    
    myStyle();
    const double A =0.1,B=2,C=-100;
    // -------------------------------------------------------------------
    //  Define Hists and Plots
    // -------------------------------------------------------------------
    TFile hfile("h.root","RECREATE","KalTest");

    TNtupleD *hTrackMonitor = new TNtupleD("track", "", "ndf:chi2:cl");

    auto grhits = new TGraphErrors(50);
    grhits->SetName("grdata");
    auto grpre = new TGraphErrors(50);
    grpre->SetName("grpre");
    auto grfilter= new TGraphErrors(50);
    grfilter->SetName("grfilter");

    LansxFormat::FormatAll(grhits,"%a %d %e %f",kBlue,kBlue,kFullCircle,1);
    LansxFormat::FormatAll(grpre,"%a %d %e %f",kOrange-3,kOrange-3,kFullCircle,1);
    LansxFormat::FormatAll(grfilter,"%a %d %e %f",kGreen,kGreen,kFullCircle,1);
    // -------------------------------------------------------------------
    //  Start event loop
    // -------------------------------------------------------------------

    Int_t ntrks = 5000;
    for (Int_t itrk=0; itrk<ntrks; itrk++) 
    {
        if(itrk % 100 == 0 ) cerr<<itrk<<" tracks done!"<<endl;
        // ---------------------------
        //  Create hits
        // ---------------------------

        Int_t nhits = 50;
        TObjArray hits(nhits);

        for (Int_t ihit=0; ihit<nhits; ihit++) {
            Double_t t  = __TSTEP__*ihit;
            Double_t dx = __SIGMAX__;
            //Double_t x  = __SLOPE__*t + __INTERCEPT__;
            Double_t x  = A*t*t+B*t+C;
            x += gRandom->Gaus(0.,dx);
            grhits->SetPoint(ihit,t,x);grhits->SetPointError(ihit,0,dx);
            hits.Add(new EXHit(t,&x,&dx));
            //cerr << "t = " << t << " x = " << x << " dx = " << dx << endl;
        }

        // ---------------------------
        //  Create Kalman Filter
        // ---------------------------
        EXKalSystem kalsys;
        // ---------------------------
        // Prepare hit iterrator 
        // ---------------------------

        TIter next(&hits);

        // ---------------------------
        //  Create a dummy site: sited 
        // ---------------------------

        if (!hits.GetEntries()) continue;

        EXHit       hitd        = *(EXHit *)next();
        hitd(0,1)   = 1.e4;      // give a huge error to x

        next.Reset();                        // rewind iterator

        EXKalSite  &sited       = *(new EXKalSite(hitd,1,3));
        // sited.Lock();                     // dummy site should not be used
        sited.SetOwner();                    // site owns states
        
        //sited.PrintDe();
        // ---------------------------
        //  Set dummy state to sited
        // ---------------------------

        TKalMatrix  svd(3,1);
        svd(0,0) = -0.2;
        svd(1,0) = 1;
        svd(2,0) = 0.1;

        TKalMatrix C(3,3);   
        C(0,0) = 1.e8;
        C(1,1) = 1.e8;
        C(2,2) = 1.e8;

        sited.Add(new EXKalState(svd,C,TVKalSite::kPredicted,3));//entry 0
        sited.Add(new EXKalState(svd,C,TVKalSite::kFiltered,3)); //entry 1

        //sited.GetState(TVKalSite::kPredicted).DebugPrint();
        //sited.GetState(TVKalSite::kFiltered).DebugPrint();
        //std::printf("----------------------------------------------------------------\n");

        //std::printf("========Debug: sited entries()==%d\n",sited.GetEntries());
        // ---------------------------
        //  Add sited to the system
        // ---------------------------
        kalsys.Add(&sited);

        //std::printf("========Debug: system entries()==%d\n",kalsys.GetEntries());
        // ---------------------------
        //  Start Kalman Filter
        // ---------------------------

        Int_t  loop = 1;
        EXHit *hitPtr;
        Int_t    ndf  = 0;
        Double_t chi2 = 0.;
        Double_t cl   = 0.;

        while ((hitPtr = (EXHit *)next())) {
            EXHit &hit = *hitPtr;
            EXKalSite & site = *(new EXKalSite(hit,1,3));
            if (kalsys.AddAndFilter(site)) 
            {
                ndf  = kalsys.GetNDF();
                chi2 = kalsys.GetChi2();
                cl   = TMath::Prob(chi2,ndf);

                if(itrk==0)
                {
                    auto curSite = dynamic_cast<EXKalSite&>(kalsys.GetCurSite());
                    auto apre = curSite.GetState(TVKalSite::kPredicted)(0,0);
                    auto bpre = curSite.GetState(TVKalSite::kPredicted)(1,0);
                    auto cpre = curSite.GetState(TVKalSite::kPredicted)(2,0);
                    auto afilter= curSite.GetState(TVKalSite::kFiltered)(0,0);
                    auto bfilter= curSite.GetState(TVKalSite::kFiltered)(1,0);
                    auto cfilter= curSite.GetState(TVKalSite::kFiltered)(2,0);

                    auto ti = curSite.GetCurHit()->GetT();

                    grpre->SetPoint(loop-1,ti,apre*ti*ti+bpre*ti+cpre);
                    grfilter->SetPoint(loop-1,ti,afilter*ti*ti+bfilter*ti+cfilter);

                    //cout<<"Loop: "<<loop<<"\t"<<hitPtr->GetX(0)<<"\t"<<kkpre*ti+bbpre<<"\t"<<kkfilter*ti+bbfilter<<endl;
                }
            } else {
                cerr << " site discarded!" << endl;
                delete &site;
            }
            loop++;
        }
        
        //cout<<"====================Filtered : "<<endl;
        //kalsys.GetState(TVKalSite::kFiltered).DebugPrint();
        //

        //cout<<"====================Smoothed : "<<endl;
        //kalsys.SmoothBackTo(2);
        //kalsys.GetState(TVKalSite::kSmoothed).DebugPrint();

        hTrackMonitor->Fill(ndf,chi2,cl);
    }

    grhits->Write();
    auto ftrue = new TF1("fture","0.1*x*x+2*x-100",0,50);
    auto myc = new TCanvas("myc","myc",800,600);
    myc->SetGrid();
    grhits->SetTitle(";t;x");
    grhits->Draw("AP");
    grpre->Draw("P");
    grfilter->Draw("P");
    ftrue->Draw("SAME");
    auto lg = new TLegend(0.2,0.75,0.5,0.92);
    lg->SetFillStyle(000);
    lg->AddEntry(grhits,"Data","p");
    lg->AddEntry(grpre,"Predicted","p");
    lg->AddEntry(grfilter,"Filtered","p");
    lg->AddEntry(ftrue,"x=0.1t^{2}+2*t-100","l");
    lg->Draw();

    myc->Write();
    hTrackMonitor->Write();
    hfile.Close();
    return 0;
}
