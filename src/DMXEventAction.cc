//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// --------------------------------------------------------------
//   GEANT 4 - Underground Dark Matter Detector Advanced Example
//
//      For information related to this code contact: Alex Howard
//      e-mail: alexander.howard@cern.ch
// --------------------------------------------------------------
// Comments
//
//                  Underground Advanced
//               by A. Howard and H. Araujo
//                    (27th November 2001)
//
// History/Additions:
// 16 Jan 2002  Added analysis
//
//
// EventAction program
// --------------------------------------------------------------

#include "DMXEventAction.hh"

// pass parameters for messengers:
#include "DMXRunAction.hh"
#include "DMXPrimaryGeneratorAction.hh"

// note DMXScintHit.hh is included in DMXEventAction.hh

#include "DMXEventActionMessenger.hh"

#include "G4AnalysisManager.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4UnitsTable.hh"
#include "G4RunManager.hh"
#include "G4Threading.hh"
#include <fstream>
#include <iomanip>

#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
DMXEventAction::DMXEventAction()
  : runAct(0),genAction(0),hitsfile(0)
{

  // create messenger
  eventMessenger = new DMXEventActionMessenger(this);

  // defaults for messenger
  drawColsFlag = "standard";
  drawTrksFlag = "all";
  saveHitsFlag = 1;

  printModulo = 1;

  // hits collections
  scintillatorCollID = -1;

  energy_pri=0;
  seeds=NULL;

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
DMXEventAction::~DMXEventAction() {

  if (hitsfile)
    {
      hitsfile->close();
      delete hitsfile;
    }
  delete eventMessenger;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void DMXEventAction::BeginOfEventAction(const G4Event* evt)
{

  //thread-local run action
  if (!runAct)
    runAct =
      dynamic_cast<const DMXRunAction*>
      (G4RunManager::GetRunManager()->GetUserRunAction());

  if (!genAction)
    genAction = dynamic_cast<const DMXPrimaryGeneratorAction*>
      (G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());


  // grab seeds
  seeds = genAction->GetEventSeeds();

  // grab energy of primary
  energy_pri = genAction->GetEnergyPrimary();

  event_id = evt->GetEventID();

  // print this information event by event (modulo n)
  if (event_id%printModulo == 0)
    {
      G4cout << "\n---> Begin of event: " << event_id << G4endl;
      G4cout << "\n     Primary Energy: " << G4BestUnit(energy_pri,"Energy")
	     << G4endl;
      //      HepRandom::showEngineStatus();
    }


  // get ID for scintillator hits collection
  if (scintillatorCollID==-1) {
    G4SDManager *SDman = G4SDManager::GetSDMpointer();
    scintillatorCollID = SDman->GetCollectionID("scintillatorCollection");
  }

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void DMXEventAction::EndOfEventAction(const G4Event* evt) {

  // check that the hits collection has been defined
  if(scintillatorCollID<0) return;

  // address hits collection
  DMXScintHitsCollection* SHC = NULL;
  G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
  if(HCE) {
    SHC = (DMXScintHitsCollection*)(HCE->GetHC(scintillatorCollID));
  }

  // event summary
  totEnergy         = 0.;
  totEnergyGammas   = 0.;
  totEnergyNeutrons = 0.;
  firstParticleE    = 0.;
  particleEnergy    = 0.;
  firstLXeHitTime   = 0.;

  firstParticleName = "";
  particleName      = "";


  // particle flags
  gamma_ev          = false;
  neutron_ev        = false;
  positron_ev       = false;
  electron_ev       = false;
  proton_ev         = false;
  other_ev          = false;
  start_gamma       = false;
  start_neutron     = false;


  // scintillator hits
  if(SHC) {
    S_hits = SHC->entries();

    for (G4int i=0; i<S_hits; i++) {
      if(i==0) {
	firstParticleName = (*SHC)[0]->GetParticle();
	firstLXeHitTime   = (*SHC)[0]->GetTime();
	firstParticleE = (*SHC)[0]->GetParticleEnergy();
	if (event_id%printModulo == 0 && S_hits > 0) {
	  G4cout << "     First hit in LXe: " << firstParticleName << G4endl;
	  G4cout << "     Number of hits in LXe: " << S_hits << G4endl;
	}
      }
      hitEnergy         = (*SHC)[i]->GetEdep();
      totEnergy        += hitEnergy;

      particleName      = (*SHC)[i]->GetParticle();
      particleEnergy    = (*SHC)[i]->GetParticleEnergy();

      if(particleName == "gamma") {
	gamma_ev = true;
	start_gamma = true;
	start_neutron = false;
      }
      else if(particleName == "neutron")
	neutron_ev = true;
      else if(particleName == "e+")
	positron_ev = true;
      else if(particleName == "e-")
	electron_ev = true;
      else if(particleName == "proton")
	proton_ev = true;
      else {
	other_ev = true;
	start_gamma = false;
	start_neutron = true;
      }

      if(start_gamma && !start_neutron)
	totEnergyGammas += hitEnergy;
      if(start_neutron && !start_gamma)
	totEnergyNeutrons += hitEnergy;
    }

    if (event_id%printModulo == 0)
      G4cout << "     Total energy in LXe: "
	     << G4BestUnit(totEnergy,"Energy") << G4endl;

  }


  // write out event summary
  if(saveHitsFlag)
    writeScintHitsToFile();

  // draw trajectories
  if(drawColsFlag=="standard" && drawTrksFlag!="none")
    drawTracks(evt);

  // print this event by event (modulo n)
  if (event_id%printModulo == 0)
    G4cout << "\n---> End of event: " << event_id << G4endl << G4endl;

}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void DMXEventAction::writeScintHitsToFile()
{

  G4String filename="hits.out";
  if (runAct)
    filename=runAct->GetsavehitsFile();



  //First time it is inkoved
  if (!hitsfile)
    {
      //check that we are in a worker: returns -1 in a master and -2 in sequential
      //one file per thread is produced ending with ".N", with N= thread number
      if (G4Threading::G4GetThreadId() >= 0)
	{
	  std::stringstream sss;
	  sss << filename.c_str() << "." << G4Threading::G4GetThreadId();
	  filename = sss.str();
	  //G4cout << "Filename is: " << filename << G4endl;
	}

      hitsfile = new std::ofstream;
      hitsfile->open(filename);
    }

  if(S_hits) {

    if(hitsfile->is_open()) {


      (*hitsfile) << std::setiosflags(std::ios::fixed)
		  << std::setprecision(4)
		  << std::setiosflags(std::ios::left)
		  << std::setw(6)
		  << event_id << "\t"
		  << energy_pri/MeV << "\t"
		  << totEnergy/MeV << "\t"
		  << S_hits  << "\t"
		  << firstParticleName << "\t"
		  << (gamma_ev    ? "gamma " : "")
		  << (neutron_ev  ? "neutron " : "")
		  << (positron_ev ? "positron " : "")
		  << (electron_ev ? "electron " : "")
		  << (other_ev    ? "other " : "")
		  << G4endl;

      if (event_id%printModulo == 0)
	G4cout << "     Event summary in file " << filename << G4endl;
    }

    G4AnalysisManager* man = G4AnalysisManager::Instance();
    G4int firstparticleIndex = 0;
    if(firstParticleName == "gamma") firstparticleIndex = 1;
    else if (firstParticleName == "neutron") firstparticleIndex = 2;
    else if(firstParticleName == "electron") firstparticleIndex = 3;
    else if(firstParticleName == "positron") firstparticleIndex = 4;
    else{
      firstparticleIndex = 5;
    }

    //Fill ntuple #2
    man->FillNtupleDColumn(2,0,event_id);
    man->FillNtupleDColumn(2,1,energy_pri/keV);
    man->FillNtupleDColumn(2,2,totEnergy);
    man->FillNtupleDColumn(2,3,S_hits);
    man->FillNtupleDColumn(2,4,firstLXeHitTime);
    man->FillNtupleDColumn(2,5,firstparticleIndex);
    man->FillNtupleDColumn(2,6,firstParticleE);
    man->FillNtupleDColumn(2,7,gamma_ev);
    man->FillNtupleDColumn(2,8,neutron_ev);
    man->FillNtupleDColumn(2,9,positron_ev);
    man->FillNtupleDColumn(2,10,electron_ev);
    man->FillNtupleDColumn(2,11,other_ev);
    man->AddNtupleRow(2);

  }

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void DMXEventAction::drawTracks(const G4Event* evt) {

  if(G4VVisManager::GetConcreteInstance()) {
    G4UImanager::GetUIpointer()->ApplyCommand("/vis/scene/notifyHandlers");
    G4TrajectoryContainer* trajContainer = evt->GetTrajectoryContainer();
    G4int n_trajectories = 0;

    if(trajContainer) n_trajectories = trajContainer->entries();
    for (G4int i=0; i<n_trajectories; i++) {
      G4Trajectory* trj = (G4Trajectory*)(*trajContainer)[i];
      if (drawTrksFlag == "all")
	trj->DrawTrajectory();
      else if ((drawTrksFlag == "charged") && (trj->GetCharge() != 0.))
	trj->DrawTrajectory();
      else if ((drawTrksFlag == "noscint")
	       && (trj->GetParticleName() != "opticalphoton"))
	trj->DrawTrajectory();
    }

    // G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/update");
  }

}
