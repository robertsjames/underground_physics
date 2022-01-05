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
// History:
// 17 Jan 2002 Alex Howard Added Analysis
// 23 Oct 2009 Luciano Pandola Removed un-necessary calls from EndOfRun()
//
// RunAction program
// --------------------------------------------------------------

#include "DMXRunActionMessenger.hh"
#include "DMXRunAction.hh"

#include "G4Run.hh"
#include "G4ios.hh"
#include "G4AnalysisManager.hh"

#include <fstream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

DMXRunAction::DMXRunAction()
{
  runMessenger = new DMXRunActionMessenger(this);
  savehitsFile = "hits.out";
  savehistFile = "dmx";
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

DMXRunAction::~DMXRunAction()
{
  delete runMessenger;
  runMessenger = 0;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DMXRunAction::BeginOfRunAction(const G4Run* aRun)
{
  //Master mode or sequential
  if (IsMaster())
    G4cout << "### Run " << aRun->GetRunID() << " starts (master)." << G4endl;
  else
    G4cout << "### Run " << aRun->GetRunID() << " starts (worker)." << G4endl;

  // Book histograms and ntuples
  Book();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DMXRunAction::EndOfRunAction(const G4Run*)
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DMXRunAction::Book()
{
  // Get/create analysis manager
  G4AnalysisManager* man = G4AnalysisManager::Instance();
  man->SetDefaultFileType("xml");

  // Open an output file
  man->OpenFile(savehistFile);
  man->SetFirstHistoId(1);
  man->SetFirstNtupleId(1);


  // ---- primary ntuple ------
  // id==1
  man->CreateNtuple("tree1", "Particle Source Energy");
  man->CreateNtupleDColumn("energy");
  man->FinishNtuple();

  // ---- secondary ntuple ------
  //id==2
  man->CreateNtuple("tree2", "Scintillation Hits Info");
  man->CreateNtupleDColumn("Event");
  man->CreateNtupleDColumn("e_prim");
  man->CreateNtupleDColumn("tot_e");
  man->CreateNtupleDColumn("s_hits");
  man->CreateNtupleDColumn("xe_time");
  man->CreateNtupleDColumn("firstpart");
  man->CreateNtupleDColumn("firstparte");
  man->CreateNtupleDColumn("gamma");
  man->CreateNtupleDColumn("neutron");
  man->CreateNtupleDColumn("posi");
  man->CreateNtupleDColumn("elec");
  man->CreateNtupleDColumn("other");
  man->FinishNtuple();

  return;

}
