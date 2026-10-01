#include "duneanaobj/StandardRecord/Proxy/SRProxy.h"

#include "CAFAna/Core/SpectrumLoader.h"
#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Core/Binning.h"
#include "CAFAna/Core/Var.h"
#include "CAFAna/Core/Cut.h"
#include "CAFAna/Core/HistAxis.h"
#include "CAFAna/Core/TruthMatching.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TF1.h"
#include "TH2.h"
#include "TPad.h"
#include "TLegend.h"
#include "CAFAna/Core/Cut.h"


#include <iostream>
using namespace ana;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//truth overlap function (PARTICLE LEVEL)

const caf::SRTrueParticleProxy* GetTruthWithLargestOverlap(const caf::SRRecoParticleProxy * part){
  if (part->truth.empty()) return nullptr;
  const caf::SRProxy * sr = part->Ancestor<caf::SRProxy>();
  size_t tidx = 0;
  float maxOverlap = 0;
  for (size_t i = 0; i < part->truthOverlap.size(); i++){
    if (part->truthOverlap[i] > maxOverlap){
        tidx = i;
        maxOverlap = part->truthOverlap[i];
    }
  };
  return caf::FindParticle(sr->mc, part->truth[tidx]);
 
};

//////////////////////////////////////////////////////////////////////////////////////////////////
//Creating selection lambda functions

//we want:
// an FV function 
// a NueCC function
// a 1 e produced function

namespace mysel {

//Lets attempt to enumerate the sample:
//rename enum and separate reco + truth
enum RecoEnum {kRejected =0, kContainedNu_eCC1e};
enum TruthEnum {kRejectedTrue =0 , kContainedNu_eCCTruth};

//interaction produces particles
bool has_particles(const caf::SRInteractionProxy* sr){
    if (sr->part.dlp.empty()) {
    return false;
  }
  return true;

}

//defining Fiducial Volume dimensions
double NDLArXLo = -346.9;
double NDLArXHi = 346.9;
double NDLArYLo = -215.5;
double NDLArYHi = 81.7;
double NDLArZLo = 418.2;
double NDLArZHi = 913.3;

double FVCrop = 25;

// interaction is within fiducial volume 
bool InFV (const caf::SRInteractionProxy* sr, double Crop= FVCrop){

  double x_position = sr->vtx.x;
  double y_position = sr->vtx.y;
  double z_position = sr->vtx.z;

  if((x_position <= (NDLArXLo+Crop)) || (x_position >= (NDLArXHi-Crop))){
      return false;
  };

  if((y_position <= (NDLArYLo+Crop)) || (y_position >= (NDLArYHi-Crop))){
    return false;
  };

  if((z_position <= (NDLArZLo+Crop)) || (z_position >= (NDLArZHi-Crop))){
    return false;
  };

  return true;
};


//primary particles are contained 

//not urgent: make particle track containemnt (looking at 1 particle)

bool AllParticlesContained(const caf::SRInteractionProxy* sr)
{
  for (const auto& p : sr->part.dlp){
    if (p.primary == 1 && p.contained !=1){
      return false;
    }
  }
  return true;
   
}
//has a primary electron 
//plot no of electrons hist
bool HasPrimaryElectron(const caf::SRInteractionProxy* sr)
{
    for (const auto& p : sr->part.dlp)
    {
        if (p.primary == 1 && p.pdg == 11)
            return true;
    }

    return false;
}


//has exactly one primary electron 
bool ExactlyOnePrimaryElectron(const caf::SRInteractionProxy* sr)
{
    int nElectron = 0;

    for (const auto& p : sr->part.dlp)
    {
        if (p.primary == 1 && p.pdg == 11)
            ++nElectron;
    }

    return nElectron == 1;
}



RecoEnum ApplySelection(const caf::SRInteractionProxy* sr)
{
    //Are particles being reconstructed?
    if (!has_particles(sr)){
      return kRejected;
    }
    //is the vertex within the FV?
    if (!InFV(sr)){
      return kRejected;
    }

    //Are all primary reconstructed particles contained?
    if (!AllParticlesContained(sr)){
      return kRejected;
    }

    //Does the reconstructed interaction produce a primary electron?
    /*
    if(!HasPrimaryElectron(sr)){
      return kRejected;
    }
      */
    if(!ExactlyOnePrimaryElectron(sr)){
          return kRejected;
    }


    //if else it is a contained Nu_e CC me interaction where m is an integer
      return kContainedNu_eCC1e;


}

//truth equivilant to this selection's FV function

bool InFVTrue (const caf::SRTrueInteractionProxy* sr){

  double x_position = sr->vtx.x;
  double y_position = sr->vtx.y;
  double z_position = sr->vtx.z;

  if((x_position <= (NDLArXLo+25)) || (x_position >= (NDLArXHi-25))){
    return false;
  }

  if((y_position <= (NDLArYLo+25)) || (y_position >= (NDLArYHi-25))){
    return false;
  }

  if((z_position <= (NDLArZLo+25)) || (z_position >= (NDLArZHi-25))){
    return false;
  }

  return true;
}




//selecting the interaction

bool IsNu_eCCTrue (const caf::SRTrueInteractionProxy* sr){

  //is the interaction a CC?
  if (sr->iscc == 0){
    return false;
  }

    //are any particles produced in the ineraction?x
  if (sr->nprim == 0){
    return false;
  }

  //are any of the particle produced an electron?
  for (const auto& p : sr->prim){
    if (p.pdg == 11){
        return true;
    }
  }
  return false;
}


//attempting to apply the selection

TruthEnum ApplySelectionTruth(const caf::SRTrueInteractionProxy* sr)
{

  //is the vertex within the FV?
  if (!InFVTrue(sr)){
    return kRejectedTrue;
  }

  if(!IsNu_eCCTrue(sr)){
    return kRejectedTrue;
  }

  //if else it is a contained Nu_e CC e
  return kContainedNu_eCCTruth;

};


}//end of namespace 

///////////////////////////////////////////////////////////////////////////////////////////////////////
//MAIN BODY OF CODE BELOW

void sel_1(){

//Files which we will be running on:
std::string fname = "/pnfs/dune/persistent/physicsgroups/dunendsim/abooth/nd-production/MiniProdN5/run-cafmaker/MiniProdN5p3_NDComplex_FHC.caf.full.sanddrift.spineonly/CAF/0000000/*.root";
  
SpectrumLoader loader(fname);

/////////////////////////////////////////////////////////////////////////////////////////////////
//Vars

//e.g creating a var, histaxis, spectrum and plotting while implamenting the selection
//creating a number of electrons Var:
const Var eNum([](const caf::SRInteractionProxy* sr){               
  int number = 0;
  for (const auto& p : sr->part.dlp){
      if (p.primary == 1 &&p.pdg == 11){
        ++number;
      }
    }
  return number;
});


//creating a series of vars  (this is still under development)


const RecoPartVar eRecoE([](const caf::SRRecoParticleProxy* p)->double{ 
      return p->E;
});

//var for True electron energy
const TruthPartVar eTrueE([](const caf::SRTrueParticleProxy* p)->double{
    return p->p.E;
});

//var for truth matched energy (for selection effeciency vs true energy)

const RecoPartVar MatchedTrueElectronEnergy(
[](const caf::SRRecoParticleProxy* p) -> double
{
    const auto* truth = GetTruthWithLargestOverlap(p);
    if(!truth) return -1.0;
    return double(truth->p.E);
});


/////////////////////////////////////////////////////////////////////////////////////////////
//histaxis

//creating binning and hist axis:
//number of e
const Binning eRecoNumBinning = Binning::Simple(5, -0.5, 4.5);
const HistAxis ERecoNumAxis("Number of electrons", eRecoNumBinning, eNum);

//Reco electron energy
const Binning eRecoEBinning = Binning::Simple(150, 0, 2);
const RecoPartHistAxis ERecoEAxis("Reco eletron energy", eRecoEBinning, eRecoE);

//true electron energy
const Binning eTrueEBinning = Binning::Simple(150, 0, 2);
const TruthPartHistAxis ETrueEAxis("True Electron energy", eTrueEBinning, eTrueE);

//truth matched electron energy

const RecoPartHistAxis MatchedTrueEAxis("Matched true electron energy",eTrueEBinning, MatchedTrueElectronEnergy);

////////////////////////////////////////////////////////////////////////////////////////////////
//Cuts

//selection application cuts
//reco
const Cut kRecoNu_eCC1e([](const caf::SRInteractionProxy* sr){
   return mysel::ApplySelection(sr) == mysel::kContainedNu_eCC1e;
});
//true
const TruthCut kTrueNu_eCC([](const caf::SRTrueInteractionProxy* sr){
  return mysel::ApplySelectionTruth(sr) == mysel::kContainedNu_eCCTruth;
});



//particle ID cuts

const TruthPartCut kTruthElectron([](const caf::SRTrueParticleProxy* p)
{
    return abs(p->pdg) == 11;
});

const RecoPartCut kRecoElectron([](const caf::SRRecoParticleProxy* p)
{
    return abs(p->pdg) == 11;
});


//electron truth match cut (PARTICLE LEVEL)
const RecoPartCut kTruthMatchElectron([](const caf::SRRecoParticleProxy * part){
  if (part->truth.empty()) return false;
  const caf::SRProxy * sr = part->Ancestor<caf::SRProxy>();
  size_t tidx = 0;
  float maxOverlap = 0;
  for (size_t i = 0; i < part->truthOverlap.size(); i++){
    if (part->truthOverlap[i] > maxOverlap){
        tidx = i;
        maxOverlap = part->truthOverlap[i];
    }
  }
  int pdg = caf::FindParticle(sr->mc, part->truth[tidx])->pdg;
  return pdg == 11;
});     


//Nu_eCC interaction truth matching cut (INTERACTION LEVEL)
const Cut kTruthMatchInteraction([](const caf::SRInteractionProxy* sr)
{
    if(sr->truth.empty()) return false;

    // Find best-matched true interaction
    size_t tidx = 0;
    float maxOverlap = 0.0;

    for(size_t i = 0; i < sr->truthOverlap.size(); ++i){
        if(sr->truthOverlap[i] > maxOverlap){
            maxOverlap = sr->truthOverlap[i];
            tidx = i;
        }
    }

    // Get the ancestor SRProxy
    const caf::SRProxy* proxy = sr->Ancestor<caf::SRProxy>();

    // Get the matched true interaction
    const caf::SRTrueInteractionProxy* truthInt =
        &proxy->mc.nu[sr->truth[tidx]];

    return truthInt->iscc &&
           abs(truthInt->pdg) == 12;

});

///////////////////////////////////////////////////////////////////////////////////////
//Spectrum
//e number
Spectrum sRecoElectronNum(loader.Interactions(RecoType::kDLP)[kRecoNu_eCC1e],ERecoNumAxis);

//reco e energy
Spectrum sRecoElectronEnergy(loader.Interactions(RecoType::kDLP)[kRecoNu_eCC1e].RecoParticles(RecoType::kDLP)[kRecoElectron],ERecoEAxis);

//true e energy
Spectrum sTrueElectronEnergy(loader.NuTruths()[kTrueNu_eCC].TruthParticles(TruePType::kPrim)[kTruthElectron],ETrueEAxis);

//truth matched spectrums:
Spectrum sTruthMatchElectronEnergy(loader.Interactions(RecoType::kDLP)[kRecoNu_eCC1e && kTruthMatchInteraction].RecoParticles(RecoType::kDLP)[kRecoElectron && kTruthMatchElectron],ERecoEAxis);
//this selection is binned in true electron energy, even though it is filled from the selected reco particles.
Spectrum sSelectedTrueE(loader.Interactions(RecoType::kDLP)[kRecoNu_eCC1e && kTruthMatchInteraction].RecoParticles(RecoType::kDLP)[kRecoElectron && kTruthMatchElectron],MatchedTrueEAxis);


//////////////////////////////////////////////////////////////////////////////////////////////////////
//Plotting
//NOTE: I am still developing the truth matching and plotting of reco / truth energy
//Much still needs investigating


loader.Go();

TFile fout("MuonEnergySpectra.root", "RECREATE");

sEnergyNuMuCC0Pi0P0K.SaveTo(&fout, "sEnergyNuMuCC0Pi0P0K");

fout.Close();

std::cout << "Spectra saved!" << std::endl;

const double pot = sRecoElectronNum.POT();

std::cout << "POT = " << pot << std::endl;





//electron number
new TCanvas;

TH1* hElectronNumber = sRecoElectronNum.ToTH1(pot, kBlue);

hElectronNumber->SetTitle("Primary electron multiplicity for Nu_eCC1e;N_{e};Events");

hElectronNumber->Draw("hist");

gPad->SaveAs("Number_of_primary_e_with1e.pdf");




//Reco Electron energy with truth matched component
new TCanvas;

TH1* hElectronRecoEnergy = sRecoElectronEnergy.ToTH1(pot, kGreen);
TH1* hElectronTruthMatchEnergy = sTruthMatchElectronEnergy.ToTH1(pot, kRed+1);

hElectronRecoEnergy->SetTitle("Reco primary electron Energy for nueCCInclusive;Reco E_{e}");

hElectronRecoEnergy->Draw("hist");
hElectronTruthMatchEnergy->Draw("hist same");

// create legend
TLegend* leg = new TLegend(0.55,0.65,0.88,0.85);

leg->AddEntry(hElectronRecoEnergy,
              "All reconstructed electrons",
              "L");

leg->AddEntry(hElectronTruthMatchEnergy,
              "Reco electrons truth matched to #nu_{e}CC true electrons",
              "L");

leg->Draw();


gPad->SaveAs("NueCCInclusive_RecoPrimaryElectronE.pdf");





//True Electron energy
new TCanvas;

TH1* hElectronTrueEnergy = sTrueElectronEnergy.ToTH1(pot, kRed);

hElectronTrueEnergy->SetTitle("Truth primary electron Energy for nueCCInclusive;Truth E_{e}");

hElectronTrueEnergy->Draw("hist");

gPad->SaveAs("NueCCInclusive_TruePrimaryElectronE.pdf");


}