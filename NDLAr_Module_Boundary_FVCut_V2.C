#include "CAFAna/Core/Spectrum.h" 
#include "CAFAna/Core/SpectrumLoader.h"
#include "CAFAna/Core/Var.h" 
#include "CAFAna/Core/Selection.h"
#include "CAFAna/AnalysisInfo.h"
#include "CAFAna/Cuts/TruthCuts.h"
#include "CAFAna/StandardRecord/Proxy/SRProxy.h"
#include "CAFAna/StandardRecord/Proxy/SRTrueInteractionProxy.h"
#include "TCanvas.h" 
#include "TGraph.h" 
#include "TGraphErrors.h" 
#include "TLegend.h" 
#include "TH1D.h" 
#include "TFile.h" 
#include "TAxis.h" 
#include "TStyle.h" 
#include <cmath> 
#include <iostream> 
#include <vector> 
#include <string> 
#include <algorithm>
#include "NDLArGeometry.h"
#include "CAFAna/Core/Cut.h"


//define fiducial volume for entire detector outer-boundaries 
//NDLArFV is 25cm 'in' from each face of the following detector bounding box:
double NDLArXLo = -346.9;
double NDLArXHi = 346.9;
double NDLArYLo = -215.5;
double NDLArYHi = 81.7;
double NDLArZLo = 418.2;
double NDLArZHi = 913.3;

//Function to check that intercation coordinates are within outer FV
bool IsInOuterFiducialVolume(double x, double y, double z,
                             const std::vector<double>& outer_margins) {
  double mx = outer_margins.size()>0 ? outer_margins[0] : 0.0;
  double my = outer_margins.size()>1 ? outer_margins[1] : 0.0;
  double mz = outer_margins.size()>2 ? outer_margins[2] : 0.0;
  
  return (x >= NDLArXLo + mx && x <= NDLArXHi - mx &&
          y >= NDLArYLo + my && y <= NDLArYHi - my &&
          z >= NDLArZLo + mz && z <= NDLArZHi - mz);
}


//Function to check whether the interaction coordinates are within the combined FV of all 70 TPCs, 
//With the specified margins
bool IsInCombinedNDLArFV(double x, double y, double z, 
                        const std::vector<double>& box_margins,
                        const std::vector<double>& outer_margins) {
    
  //Check coordinates are in outer FV first, if not return false
  if (!IsInOuterFiducialVolume(x, y, z, outer_margins)) {
      return false;
  }

  // Check individual TPC module box boundaries
  double b_mx = box_margins.size() > 0 ? box_margins[0] : 0.0;
  double b_my = box_margins.size() > 1 ? box_margins[1] : 0.0;
  double b_mz = box_margins.size() > 2 ? box_margins[2] : 0.0;

  //All event coordinates that pass the outer FV cut are then checked against the box's FV
  return NDLArGeo::IsInBoxesFiducialVolume(x, y, z, b_mx, b_my, b_mz);

}



void NDLAr_Module_Boundary_FVCut(){

// Define FV margin vectors: {margin_x, margin_y, margin_z} in cm
const std::vector<double> kBoxMargins   = {5.0, 5.0, 5.0}; // Per-box margins
const std::vector<double> kAVMargins = {25.0, 25.0, 25.0}; // Margins to cut in on in avtive volume

const Cut kIsInFiducialVolume = ([](const caf::SRinteractionProxy* sr) {
    
    // Get the true vertex position from the SRinteractionProxy
    double x = sr->vtx.x;
    double y = sr->vtx.y;
    double z = sr->vtx.z;

    // Check if the vertex is inside the detetcor-encompasing FV margins, and inside the individual TPC boxes with specified FV margins
    return IsInCombinedNDLArFV(x, y, z, kBoxMargins, kAVMargins);
});

}
