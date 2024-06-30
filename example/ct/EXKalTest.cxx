/*********************************************************************
 * Author           : Lan-sx
 * Email            : shexin@ihep.ac.cn
 * Last modified    : 2024-04-27 17:05
 * Filename         : EXKalTest.cxx
 * Description      : Example import from Kaltest 
 * Update           : remove ClassDef() ClassImp(), run on WSL2 successfully
 * ******************************************************************/
#include "TKalDetCradle.h"    // from KalTrackLib
#include "TKalTrackState.h"   // from KalTrackLib
#include "TKalTrackSite.h"    // from KalTrackLib
#include "TKalTrack.h"        // from KalTrackLib

#include "EXKalTest.h"
#include "EXKalDetector.h"
#include "EXEventGen.h"
#include "EXHit.h"

#include "TNtupleD.h"         // from ROOT
#include "TFile.h"            // from ROOT
#include "TGraph.h"
#include "TGraph2D.h"
#include "Lansxlogon.h"

#include <iostream>

static const Bool_t gkDir = kIterBackward;
//static const Bool_t gkDir = kIterForward;
//#define __STRAIGHT_TRACK__
#define __MS_OFF__

using namespace std;

int main (Int_t argc, Char_t **argv)
{
    gROOT->SetBatch();
    TApplication app("EXKalTest", &argc, argv, 0, 0);

    TFile hfile("hh_debug.root","RECREATE","KalTest");

#ifdef __STRAIGHT_TRACK__
    TNtupleD *hTrackMonitor = new TNtupleD("track", "", "ndf:chi2:cl:phi0:tanl:cpa");
#else
    TNtupleD *hTrackMonitor = new TNtupleD("track", "", "ndf:chi2:cl:cpa");
#endif


    Double_t pt      = 1.;
    Double_t t0in    = 0.;
    Int_t    nevents = 1;

    switch (argc) {
    case 4: 
        nevents = atoi(argv[1]);
        pt      = atof(argv[2]);
        t0in    = atof(argv[3]);
        break;
    case 3: 
        nevents = atoi(argv[1]);
        pt      = atof(argv[2]);
        break;
    case 2:
        nevents = atoi(argv[1]);
        break;
    case 1:
        break;
    default:
        cerr << "Too many command line arguments!" << endl;
        return -1;
        //abort();
    }

    // ===================================================================
    //  Prepare a detector
    // ===================================================================

    TObjArray     kalhits;    // hit buffer
    TKalDetCradle cradle;     // detctor system
    EXKalDetector detector;   // CT detector
    std::printf("==============> test: %d entries~\n",detector.GetEntries());

    cradle.Install(detector); // install detector into its cradle
#ifdef __MS_OFF__
    cradle.SwitchOffMS();     // switch off multiple scattering
#endif

    // ===================================================================
    //  Prepare a Event Generator
    // ===================================================================

    EXEventGen gen(cradle, kalhits);
    gen.SetT0(t0in);

    // ===================================================================
    //  Event loop
    // ===================================================================

    for (Int_t eventno = 0; eventno < nevents; eventno++) { 
        cerr << "------ Event " << eventno << " ------" << endl;

        // ---------------------------
        //  Reset hit data
        // ---------------------------

        kalhits.Delete();

        // ============================================================
        //  Generate a partcle
        // ============================================================

#ifdef __STRAIGHT_TRACK__
        TStraightTrack hel = gen.GenerateStraightTrack(pt);
#else
        THelicalTrack hel = gen.GenerateHelix(pt);
#endif

        // ============================================================
        //  Swim the particle in detector
        // ============================================================
        cout<<"==========> "<<hel.GetRho()<<"\t"<<hel.GetKappa()<<"\t"<<hel.GetPtoR()<<endl;
        gen.Swim(hel);
        std::printf("=====================Swim Done!\n");
        //std::printf("================test=========\n");
        //std::printf("===rho = %.5f [m] \t alpha = %.5f Kappa=%.5f\n",hel.GetRho(),hel.GetPtoR(),hel.GetKappa());
        //std::printf("================test=========\n");
        // ============================================================
        //  Do Kalman Filter
        // ============================================================

        auto grhitsxy = new TGraph(kalhits.GetEntries());
        auto grhits3d = new TGraph2D(kalhits.GetEntries());
        LansxFormat::FormatAll(grhitsxy,"%a %d%e%f",kBlue,kBlue,kFullCircle,1);
        LansxFormat::FormatAll(grhits3d,"%a %d%e%f",kBlue,kBlue,kFullCircle,1);

        Int_t i1, i2, i3;
        if (gkDir == kIterBackward) {
            i3 = 0;
            i1 = kalhits.GetEntries() - 1;
            i2 = i1 / 2;
        } else {
            i1 = 0;
            i3 = kalhits.GetEntries() - 1;
            i2 = i3 / 2;
        }

        for(int k=0;k<kalhits.GetEntries();++k)
        {
            EXHit hitd = *dynamic_cast<EXHit *>(kalhits.At(k));
            TVector3 xx = hitd.GetMeasLayer().HitToXv(hitd);
            grhitsxy->SetPoint(k,xx.X(),xx.Y());
            grhits3d->SetPoint(k,xx.X(),xx.Y(),xx.Z());
        }

        // ---------------------------
        //  Create a dummy site: sited
        // ---------------------------

        EXHit hitd = *dynamic_cast<EXHit *>(kalhits.At(i1));
        hitd(0,1) = 1.e6;   // give a huge error to d
        hitd(1,1) = 1.e6;   // give a huge error to z

        TKalTrackSite &sited = *new TKalTrackSite(hitd);
        sited.SetOwner();   // site owns states

        // ---------------------------
        // Create initial helix
        // ---------------------------

        EXHit   &h1 = *dynamic_cast<EXHit *>(kalhits.At(i1));   // first hit
        EXHit   &h2 = *dynamic_cast<EXHit *>(kalhits.At(i2));   // last hit
        EXHit   &h3 = *dynamic_cast<EXHit *>(kalhits.At(i3));   // middle hit
        TVector3 x1 = h1.GetMeasLayer().HitToXv(h1);
        TVector3 x2 = h2.GetMeasLayer().HitToXv(h2);
        TVector3 x3 = h3.GetMeasLayer().HitToXv(h3);
        THelicalTrack helstart(x1, x2, x3, h1.GetBfield(), gkDir); // initial helix 

        // ---------------------------
        //  Set dummy state to sited
        // ---------------------------

        static TKalMatrix svd(kSdim,1);
        svd(0,0) = 0.;
        svd(1,0) = helstart.GetPhi0();
        svd(2,0) = helstart.GetKappa();
        svd(3,0) = 0.;
        svd(4,0) = helstart.GetTanLambda();
        if (kSdim == 6) svd(5,0) = 0.;
        
#ifdef __STRAIGHT_TRACK__
        svd(2,0) = 0.;
#endif
        static TKalMatrix C(kSdim,kSdim);
        for (Int_t i=0; i<kSdim; i++) {
            C(i,i) = 1.e4;   // dummy error matrix
        }

        sited.Add(new TKalTrackState(svd,C,sited,TVKalSite::kPredicted));
        sited.Add(new TKalTrackState(svd,C,sited,TVKalSite::kFiltered));

        // ---------------------------
        //  Add sited to the kaltrack
        // ---------------------------

        TKalTrack kaltrack;    // a track is a kal system
        kaltrack.SetOwner();   // kaltrack owns sites
        kaltrack.Add(&sited);  // add the dummy site to the track

        // ---------------------------
        //  Prepare hit iterrator
        // ---------------------------
        TIter next(&kalhits, gkDir);   // come in to IP  outer->inner?

        // ---------------------------------------------------------------------------------
        //  Start Kalman Filter
        // ---------------------------------------------------------------------------------

        EXHit *hitp = 0;
        while ((hitp = dynamic_cast<EXHit *>(next()))) {     // loop over hits
            TKalTrackSite  &site = *new TKalTrackSite(*hitp); // create a site for this hit
            if (!kaltrack.AddAndFilter(site)) {               // add and filter this site
                cerr << " site discarded!" << endl;           
                delete &site;                                  // delete this site, if failed
            }
        }
        //kaltrack.SmoothBackTo(1);                          // smooth back.

        // ============================================================
        //  Monitor Fit Result
        // ============================================================

        Int_t    ndf  = kaltrack.GetNDF();
        Double_t chi2 = kaltrack.GetChi2();
        Double_t cl   = TMath::Prob(chi2, ndf);
        Double_t cpa  = kaltrack.GetCurSite().GetCurState()(2, 0);

#ifdef __STRAIGHT_TRACK__
        Double_t tanl1  = kaltrack.GetCurSite().GetCurState()(4, 0);
        Double_t phi0  = kaltrack.GetCurSite().GetCurState()(1, 0);
        hTrackMonitor->Fill(ndf, chi2, cl, phi0, tanl1, cpa);
#else
        hTrackMonitor->Fill(ndf, chi2, cl, cpa);
        //std::printf("###############%d : %.6f\n",eventno,cpa);
#endif
        grhitsxy->Write();
        grhits3d->Write();
    }

    std::printf("===================Code End!\n");


    hTrackMonitor->Write();
    hfile.Close();

    return 0;
}
