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
#include "TGraphErrors.h"
#include "Lansxlogon.h"

#define __TSTEP__     1.
#define __SLOPE__     0.2
#define __INTERCEPT__ 2.
#define __SIGMAX__    0.1
//#define __DEBUG__

using namespace std;

int main (Int_t argc, Char_t **argv)
{
    gROOT->SetBatch();
    TApplication app("EXKalTest", &argc, argv, 0, 0);
    
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

    Int_t ntrks = 1;
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
            Double_t x  = __SLOPE__*t + __INTERCEPT__;
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
                                 
        //auto hit0 = *(EXHit *)(hits.At(0));
        //auto hit1 = *(EXHit *)(hits.At(1));
        //std::printf("================Debug print===========\n");
        //std::printf("=========> hitd %f \n",hitd.GetT());
        //std::printf("=========> hit0 %f \n",hit0.GetT());
        //std::printf("=========> hit1 %f \n",hit1.GetT());

        next.Reset();                        // rewind iterator

        EXKalSite  &sited       = *(new EXKalSite(hitd));

        // sited.Lock();                     // dummy site should not be used
        sited.SetOwner();                    // site owns states

        // ---------------------------
        //  Set dummy state to sited
        // ---------------------------

        Int_t i1st  = 0;
        Int_t ilst  = hits.GetEntries() - 1;
        EXHit &h1st = *(EXHit *)hits.At(i1st);
        EXHit &hlst = *(EXHit *)hits.At(ilst);

        TKalMatrix  svd(2,1);
        svd(0,0) = (hlst.GetX(0) - h1st.GetX(0))/(hlst.GetT() - h1st.GetT());
        svd(1,0) = h1st.GetX(0);

        TKalMatrix C(2,2);   
        C(0,0) = 1.e8;
        C(1,1) = 1.e8;

        sited.Add(new EXKalState(svd,C,TVKalSite::kPredicted));//entry 0
        sited.Add(new EXKalState(svd,C,TVKalSite::kFiltered)); //entry 1

        sited.GetState(TVKalSite::kPredicted).DebugPrint();
        sited.GetState(TVKalSite::kFiltered).DebugPrint();
        std::printf("----------------------------------------------------------------\n");

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
        double predictX =0.;
        double filterX =0.;

        while ((hitPtr = (EXHit *)next())) {
#ifdef __DEBUG__ 
            cerr << "----------------------" << "loop = " << loop        << endl;
#endif
            EXHit &hit = *hitPtr;
            EXKalSite & site = *(new EXKalSite(hit));
            if (kalsys.AddAndFilter(site)) {
                if(loop < 3)
                {
                    kalsys.GetState(TVKalSite::kPredicted).DebugPrint();
                    kalsys.GetState(TVKalSite::kFiltered).DebugPrint();
                }
#ifdef __DEBUG__
                //std::printf("========Debug: loop%d > site entries=%d system entries=%d\n",loop,site.GetEntries(),kalsys.GetEntries());
                site.DebugPrint();
                kalsys.GetState(TVKalSite::kFiltered).DebugPrint();
#endif
                ndf  = kalsys.GetNDF();
                chi2 = kalsys.GetChi2();
                cl   = TMath::Prob(chi2,ndf);
                auto curSite = dynamic_cast<EXKalSite&>(kalsys.GetCurSite());
                auto kkpre = curSite.GetState(TVKalSite::kPredicted)(0,0);
                auto bbpre = curSite.GetState(TVKalSite::kPredicted)(1,0);
                auto kkfilter= curSite.GetState(TVKalSite::kFiltered)(0,0);
                auto bbfilter= curSite.GetState(TVKalSite::kFiltered)(1,0);

                auto ti = curSite.GetCurHit()->GetT();

                grpre->SetPoint(loop-1,ti,kkpre*ti+bbpre);
                grfilter->SetPoint(loop-1,ti,kkfilter*ti+bbfilter);
                
                cout<<"Loop: "<<loop<<"\t"<<hitPtr->GetX(0)<<"\t"<<kkpre*ti+bbpre<<"\t"<<kkfilter*ti+bbfilter<<endl;
                //cout<<"=====> "<<curSite.GetCurHit()->GetT()<<"\t"<<hit.GetT()<<endl;
                //cout<<"====> "<<hit.GetT()<<"\t"<<curSite
                //auto kk = kalsys.GetCurSite().GetState(TVKalSite::kPredicted)(0,0);
                //auto bb = kalsys.GetCurSite().GetState(TVKalSite::kPredicted)(1,0);
#ifdef __DEBUG__ 
                cerr << " ndf = "  << ndf
                     << " chi2 = " << chi2 << endl;
                cerr << " cl = " << cl << endl;
                kalsys.GetState(TVKalSite::kPredicted).Print();
#endif
            } else {
                cerr << " site discarded!" << endl;
                delete &site;
            }
            loop++;
        }

        kalsys.SmoothBackTo(2);
#ifdef __DEBUG__ 
        cerr << " Smoothed-------------------" << endl;
        kalsys.GetState(TVKalSite::kSmoothed).DebugPrint();
#endif   

        hTrackMonitor->Fill(ndf,chi2,cl);

    }
    grhits->Write();
    grpre->Write();
    grfilter->Write();
    auto myc = new TCanvas("myc","myc",800,600);
    grhits->Draw("AP");
    grpre->Draw("P");
    grfilter->Draw("P");
    myc->Write();
    hTrackMonitor->Write();
    hfile.Close();
    return 0;
}
