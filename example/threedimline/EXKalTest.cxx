/*********************************************************************
 * Author           : Lan-sx
 * Email            : shexin@ihep.ac.cn
 * Last modified    : 2024-06-29 11:00
 * Filename         : EXKalTest.cxx
 * Description      : Example import from KalTest, remove ClassDef(),ClassImp()
 * Update           : Kalman fit for 3d straight line
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
#include "TGraph.h"
#include "TGraph2D.h"
#include "TLegend.h"
#include "Lansxlogon.h"

//#define __DEBUG__

using namespace std;

int main (Int_t argc, Char_t **argv)
{
    gROOT->SetBatch();
    TApplication app("EXKalTest", &argc, argv, 0, 0);
    
    LansxFormat::myStyle();
    const double A =0.2,B=-0.1,C=-0.0001,D=15.5;
    // -------------------------------------------------------------------
    //  Define Hists and Plots
    // -------------------------------------------------------------------
    TFile hfile("h.root","RECREATE","KalTest");

    TNtupleD *hTrackMonitor = new TNtupleD("track", "", "ndf:chi2:cl");

    auto grhitsxy = new TGraph(50);
    grhitsxy->SetName("grxydata");
    auto grhitsxz = new TGraph(50);
    grhitsxz->SetName("grxzdata");
    auto grhits3d = new TGraph2D(50);
    grhits3d->SetName("gr3ddata");
    //auto grpre = new TGraphErrors(50);
    //grpre->SetName("grpre");
    //auto grfilter= new TGraphErrors(50);
    //grfilter->SetName("grfilter");

    LansxFormat::FormatAll(grhitsxy,"%a %d %e %f",kBlue,kBlue,kFullCircle,1);
    LansxFormat::FormatAll(grhitsxz,"%a %d %e %f",kBlue,kBlue,kFullCircle,1);
    LansxFormat::FormatAll(grhits3d,"%a %d %e %f",kBlue,kBlue,kFullCircle,1);
    // -------------------------------------------------------------------
    //  Start event loop
    // -------------------------------------------------------------------

    Int_t ntrks = 2000;
    for (Int_t itrk=0; itrk<ntrks; itrk++) 
    {
        if(itrk % 100 == 0 ) cerr<<itrk<<" tracks done!"<<endl;
        // ---------------------------
        //  Create hits
        // ---------------------------

        Int_t nhits = 50;
        TObjArray hits(nhits);

        for (Int_t ihit=0; ihit<nhits; ihit++) {
            Double_t xhit  = 6.*ihit;

            Double_t dy = 0.02;
            Double_t yhit  = A*xhit+B;
            yhit += gRandom->Gaus(0.,dy);
            
            Double_t dz = 0.05;
            Double_t zhit = C*xhit+D;
            zhit +=gRandom->Gaus(0.,dz);

            grhitsxy->SetPoint(ihit,xhit,yhit);
            grhitsxz->SetPoint(ihit,xhit,zhit);
            grhits3d->SetPoint(ihit,xhit,yhit,zhit);
            Double_t yzhit[2] = {yhit,zhit};
            Double_t dyzhit[2] = {dy,dz};
            hits.Add(new EXHit(xhit,yzhit,dyzhit,2));
            //cerr << "t = " << t << " x = " << x << " dx = " << dx << endl;
        }
        //auto exhit_item = dynamic_cast<EXHit*>(hits.At(1));
        //exhit_item->DebugPrint();

//#if 0     
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
        hitd(0,1)   = 1.e4;      // give a huge error to y
        hitd(1,1)   = 1.e4;      // give a huge error to z

        next.Reset();                        // rewind iterator

        EXKalSite  &sited       = *(new EXKalSite(hitd,2,4));
        // sited.Lock();                     // dummy site should not be used
        sited.SetOwner();                    // site owns states
        
        //sited.PrintDe();
        // ---------------------------
        //  Set dummy state to sited
        // ---------------------------

        TKalMatrix  svd(4,1);
        svd(0,0) = 0.18;
        svd(1,0) = 0.;
        svd(2,0) = 0.;
        svd(3,0) = 15.;

        TKalMatrix C(4,4);   
        C(0,0) = 1.e8;
        C(1,1) = 1.e8;
        C(2,2) = 1.e8;
        C(3,3) = 1.e8;

        sited.Add(new EXKalState(svd,C,TVKalSite::kPredicted,4));//entry 0
        sited.Add(new EXKalState(svd,C,TVKalSite::kFiltered,4)); //entry 1

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
            EXKalSite & site = *(new EXKalSite(hit,2,4));
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
                    auto dpre = curSite.GetState(TVKalSite::kPredicted)(3,0);
                    auto afilter= curSite.GetState(TVKalSite::kFiltered)(0,0);
                    auto bfilter= curSite.GetState(TVKalSite::kFiltered)(1,0);
                    auto cfilter= curSite.GetState(TVKalSite::kFiltered)(2,0);
                    auto dfilter= curSite.GetState(TVKalSite::kFiltered)(3,0);

                    auto ti = curSite.GetCurHit()->GetT();
                    cout<<"xi "<<ti<<"\t"<<apre<<"\t"<<bpre<<"\t"<<cpre<<"\t"<<dpre<<endl;
                    cout<<"xi "<<ti<<"\t"<<afilter<<"\t"<<bfilter<<"\t"<<cfilter<<"\t"<<dfilter<<endl;
                    //grpre->SetPoint(loop-1,ti,apre*ti*ti+bpre*ti+cpre);
                    //grfilter->SetPoint(loop-1,ti,afilter*ti*ti+bfilter*ti+cfilter);

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
        kalsys.SmoothBackTo(2);
        //kalsys.GetState(TVKalSite::kSmoothed).DebugPrint();

        hTrackMonitor->Fill(ndf,chi2,cl);
//#if 0
//#endif
    }
    grhitsxy->Write();
    grhitsxz->Write();
    grhits3d->Write();
    hTrackMonitor->Write();
    hfile.Close();
    return 0;
}
