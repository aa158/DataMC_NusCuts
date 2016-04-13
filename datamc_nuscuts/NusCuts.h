#pragma once

#include "CAFAna/Core/Cut.h"
#include "CAFAna/Cuts/TimingCuts.h"
#include "StandardRecord/StandardRecord.h"

namespace ana
{

  const double NDL = -191.;
  const double NDR =  192.;
  const double NDB = -187.;
  const double NDT =  194.;
  const double NDF =    0.;
  const double NDE = 1270.;

  const double FDL = -758.;
  const double FDR =  765.;
  const double FDB = -749.;
  const double FDT =  765.;
  const double FDF =    0.;
  const double FDE = 5962.;

  /// Data Quality cuts from docdb 14241
  const Cut kNusEventQuality({
      "slc.ncontplanes", "vtx.nelastic",
      "sel.nuecosrej.hitsperplane",
      "shw.nshwlid", "shw.shwlid.gap"
    },
    [](const caf::StandardRecord* sr)
    {
      if(sr->sel.nuecosrej.hitsperplane >= 8) return false;
      if(sr->shw.nshwlid == 0)                return false;
      if(sr->shw.shwlid[0].gap >= 100.)       return false;
      if(sr->slc.ncontplanes <= 2)            return false;
      if(sr->vtx.nelastic == 0)               return false;
      return true;
    });
  
  /// FD Fiducial volume from docdb 14241
  const Cut kNusFDFiducial({"vtx.nelastic", "vtx.elastic.*"},
    [](const caf::StandardRecord* sr)
    {
      // can't use the assertion for N-1 plots
      //assert(sr->vtx.nelastic > 0 && "Must apply DQ cuts");
      if(sr->vtx.nelastic == 0) return false;

      const TVector3 vtx = sr->vtx.elastic[0].vtx;
      if(vtx.X() < -720.0) return false;
      if(vtx.X() >  720.0) return false;
      if(vtx.Y() < -720.0) return false;
      if(vtx.Y() >  300.0) return false;
      if(vtx.Z() <   50.0) return false;
      if(vtx.Z() > 5450.0) return false;
      return true;
    });

  const Cut kNusFDFidLoose({"vtx.nelastic", "vtx.elastic.*"},
    [](const caf::StandardRecord* sr)
    {
      double buf = 10.;
      if(sr->vtx.nelastic == 0) return false; 

      const TVector3 vtx = sr->vtx.elastic[0].vtx;
      if(vtx.X() < (FDL + buf)) return false;
      if(vtx.X() > (FDR - buf)) return false;
      if(vtx.Y() < (FDB + buf)) return false;
      if(vtx.Y() > (FDT - buf)) return false;
      if(vtx.Z() < (FDF + buf)) return false;
      if(vtx.Z() > (FDE - buf)) return false;
      return true;
    });

  /// Containment variable for NC events from docdb 14241
  const Cut kNusContain({
      "sel.nuecosrej.startfront", "sel.nuecosrej.stopfront",
      "sel.nuecosrej.startback",  "sel.nuecosrej.stopback",
      "sel.nuecosrej.starteast", "sel.nuecosrej.stopeast",
      "sel.nuecosrej.startwest", "sel.nuecosrej.stopwest",
      "sel.nuecosrej.starttop", "sel.nuecosrej.stoptop",
      "sel.nuecosrej.startbottom", "sel.nuecosrej.stopbottom"
    },
    [](const caf::StandardRecord* sr)
    {
      const caf::SRNueCosRej& cr = sr->sel.nuecosrej;

      if(std::min(cr.starteast,   cr.stopeast)   < 10) return false;
      if(std::min(cr.startwest,   cr.stopwest)   < 10) return false;
      if(std::min(cr.starttop,    cr.stoptop)    < 10) return false;
      if(std::min(cr.startbottom, cr.stopbottom) < 10) return false;
      if(std::min(cr.startfront,  cr.stopfront)  < 10) return false;
      if(std::min(cr.startback,   cr.stopback)   < 10) return false;
      return true;
    });

  /// Cut that is more CC rejection than NC selection from docdb 14241
  const Cut kNusNCSel({
      "slc.nhit", "trk.nkalman", "trk.kalman.len",
      "sel.remid.pid", "sel.elecid.ann"
    },
    [](const caf::StandardRecord* sr)
    {

      if(sr->slc.nhit >= 200)            return false;
      if(sr->trk.nkalman == 0)           return false;
      if(sr->trk.kalman[0].len >= 400.0) return false;
      if(sr->sel.remid.pid >= 0.6)       return false;
      if(sr->sel.elecid.ann >= 0.5)      return false;
      return true;
    });
  
  /// Decaf cut motivated by doc-db 15152
  const Cut kNusFDDecafCut({"sel.cosrej.numucontpid"},
			   [](const caf::StandardRecord* sr)
			   {
			     if(sr->sel.cosrej.numucontpid <= 0.42) return false;
			     return true;
			   });
  
  /// Cosmic rejection for the NC sample from docdb 14241
  const Cut kNusCosRej({
      "sel.cosrej.numucontpid", "sel.nuecosrej.partptp",
      "sel.elecid.nshwlid", "sel.elecid.shwlid.ismuon",
      "slc.calE", "slc.nhit"
    },
    [](const caf::StandardRecord* sr)
    {
      double numucontpid = sr->sel.cosrej.numucontpid;
      double partptp = sr->sel.nuecosrej.partptp;
      if(sr->sel.elecid.nshwlid == 0) return false;
      if(numucontpid <= 0.5) return false;
      if(partptp >= 0.8) return false;
      for(int i = 0, n = sr->sel.elecid.nshwlid; i < n; ++i) {
        if(sr->sel.elecid.shwlid[i].ismuon == 1) return false;
      }
      double nhit = (double)sr->slc.nhit;
      if(nhit <= 0.) return false;
      if(sr->slc.calE/nhit <= 0.02) return false;
      return true;
    });
  
  // Cuts that define the full cuts for the Far Detector
  const Cut kNusFDPresel = kNusEventQuality && kNusFDFiducial && kNusContain;
  const Cut kNusFD = kNusFDPresel && kNusNCSel && kNusCosRej;


  //************** Near Detector cuts ****************************
  
  /// ND Fiducial volume from docdb XXXXX
  // TODO: Add doc-db number
  const Cut kNusNDFiducial({"vtx.nelastic", "vtx.elastic.*"},
    [](const caf::StandardRecord* sr)
    {
      assert(sr->vtx.nelastic > 0 && "Must apply DQ cuts");
      
      
      if(sr->vtx.elastic[0].vtx.X() < -140.0) return false; 
      if(sr->vtx.elastic[0].vtx.X() >  140.0) return false;
      if(sr->vtx.elastic[0].vtx.Y() < -140.0) return false;  
      if(sr->vtx.elastic[0].vtx.Y() >  140.0) return false; 
      if(sr->vtx.elastic[0].vtx.Z() <  100.0) return false; 
      if(sr->vtx.elastic[0].vtx.Z() > 1000.0) return false; 
      return true;
    });

  /// Loose containment for the ND for decafs motivated by doc-db XXXXX
  // TODO: Add doc-db number
  const Cut kNusNDFidLoose({"vtx.nelastic", "vtx.elastic.*"},
    [](const caf::StandardRecord* sr)
    {
      double buf = 10.;
      if(sr->vtx.nelastic == 0) return false;

      const TVector3 vtx = sr->vtx.elastic[0].vtx;
      if(vtx.X() < (NDL + buf)) return false;
      if(vtx.X() > (NDR - buf)) return false;
      if(vtx.Y() < (NDB + buf)) return false;
      if(vtx.Y() > (NDT - buf)) return false;
      if(vtx.Z() < (NDF + buf)) return false;
      if(vtx.Z() > (NDE - buf)) return false;
      
      return true;
    });
  
  /// Extra containment cut for the ND
  // TODO: Add doc-db number
  const Cut kNusNDHarshTrk({"trk.nkalman",
      "trk.kalman.start.fX", "trk.kalman.stop.fX",
      "trk.kalman.start.fY", "trk.kalman.stop.fY",
      "trk.kalman.start.fZ", "trk.kalman.stop.fZ"
    },
    [](const caf::StandardRecord* sr)
    {
      double buf = 25.;
      if(sr->trk.nkalman == 0) return false;
      
      for(unsigned int i = 0, n = sr->trk.nkalman; i < n; ++i) {
        TVector3 start = sr->trk.kalman[i].start;
        TVector3 stop  = sr->trk.kalman[i].stop;
        if(std::min(start.X(), stop.X()) < (NDL + buf)) return false;
        if(std::max(start.X(), stop.X()) > (NDR - buf)) return false;
        if(std::min(start.Y(), stop.Y()) < (NDB + buf)) return false;
        if(std::max(start.Y(), stop.Y()) > (NDT - buf)) return false;
        if(std::min(start.Z(), stop.Z()) < (NDF + buf)) return false;
        if(std::max(start.Z(), stop.Z()) > (NDE - buf)) return false;
      }
      
      return true;
    });
  
  // Cuts that define the full cuts for the Near Detector
  const Cut kNusNDPresel = kNusEventQuality && kNusNDFiducial && kNusContain && kNusNDHarshTrk;
  const Cut kNusND = kNusNDPresel && kNusNCSel;
  
}
