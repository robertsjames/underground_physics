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
// DetectorConstruction program
// --------------------------------------------------------------

#include "DMXDetectorConstruction.hh"
#include "DMXDetectorMessenger.hh"

#include "DMXScintSD.hh"


#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4Isotope.hh"
#include "G4UnitsTable.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"

#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpBoundaryProcess.hh"

#include "G4FieldManager.hh"
#include "G4UniformElectricField.hh"
#include "G4TransportationManager.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4EqMagElectricField.hh"
#include "G4ClassicalRK4.hh"
#include "G4ChordFinder.hh"

#include "G4SDManager.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4UserLimits.hh"

#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
DMXDetectorConstruction::DMXDetectorConstruction()
{
  // create commands for interactive definition of time cuts:
  detectorMessenger = new DMXDetectorMessenger(this);

  theUserLimitsForRoom     = 0;
  theUserLimitsForDetector = 0;
  // default time cut = infinite
  //  - note also number of steps cut in stepping action = MaxNoSteps
  theMaxTimeCuts      = DBL_MAX;
  theMaxStepSize      = DBL_MAX;
  theDetectorStepSize = DBL_MAX;
  theRoomTimeCut      = 1000. * nanosecond;
  theMinEkine         = 250.0*eV; // minimum kinetic energy required in volume
  theRoomMinEkine     = 250.0*eV; // minimum kinetic energy required in volume

  //Zero the G4Cache objects to contain logical volumes
  LXeSD.Put(0);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
DMXDetectorConstruction::~DMXDetectorConstruction()
{
  delete theUserLimitsForRoom;
  delete theUserLimitsForDetector;
  delete detectorMessenger;
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void DMXDetectorConstruction::DefineMaterials()
{

#include "DMXDetectorMaterial.icc"

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4VPhysicalVolume* DMXDetectorConstruction::Construct() {

  DefineMaterials();

  // DefineField();

  // make colours
  G4Colour  white   (1.0, 1.0, 1.0) ;
  G4Colour  grey    (0.5, 0.5, 0.5) ;
  G4Colour  lgrey   (.85, .85, .85) ;
  G4Colour  red     (1.0, 0.0, 0.0) ;
  G4Colour  blue    (0.0, 0.0, 1.0) ;
  G4Colour  cyan    (0.0, 1.0, 1.0) ;
  G4Colour  magenta (1.0, 0.0, 1.0) ;
  G4Colour  yellow  (1.0, 1.0, 0.0) ;
  G4Colour  orange  (.75, .55, 0.0) ;
  G4Colour  lblue   (0.0, 0.0, .75) ;
  G4Colour  lgreen  (0.0, .75, 0.0) ;
  G4Colour  green   (0.0, 1.0, 0.0) ;
  G4Colour  brown   (0.7, 0.4, 0.1) ;


  //  un-used colours:
  //  G4Colour  black   (0.0, 0.0, 0.0) ;



  // Universe - room wall - CONCRETE ****************************************

  //NB: measured INSIDE of lab, therefore have to add twice wall thickness
  G4double wallThick   = 24.*cm;
  G4double worldWidth  = 470.0*cm + 2.*wallThick; // "x"
  G4double worldLength = 690.0*cm + 2.*wallThick; // "y"
  G4double worldHeight = 280.0*cm + 2.*wallThick; // "z"

  G4Box* world_box = new G4Box
     ("world_box", 0.5*worldWidth, 0.5*worldLength, 0.5*worldHeight );
  world_log  = new G4LogicalVolume(world_box, world_mat, "world_log");
  world_phys = new G4PVPlacement(0, G4ThreeVector(0.,0.,0.),
     "world_phys", world_log, NULL, false,0);

  //  G4VisAttributes* world_vat= new G4VisAttributes(white);
  world_log->SetVisAttributes(G4VisAttributes::GetInvisible());
  //world_vat->SetVisibility(true);
  //world_vat->SetVisibility(false);
  //world_log->SetVisAttributes(world_vat);


  // Lab Space - AIR ********************************************************

  G4double labWidth  = worldWidth  - 2.*wallThick; //X
  G4double labLength = worldLength - 2.*wallThick; //Y
  G4double labHeight = worldHeight - 2.*wallThick; //Z

  G4Box* lab_box = new G4Box
     ("lab_box", 0.5*labWidth, 0.5*labLength, 0.5*labHeight );
  lab_log  = new G4LogicalVolume(lab_box, lab_mat, "lab_log");
  lab_phys = new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), "lab_phys",
     lab_log, world_phys, false,0);

  G4VisAttributes* lab_vat= new G4VisAttributes(white);
  //  lab_log->SetVisAttributes(G4VisAttributes::GetInvisible());
  //  lab_vat->SetVisibility(true);
  lab_vat->SetVisibility(false);
  lab_log->SetVisAttributes(lab_vat);

  // Now start with detector assembly:

  // first LN2 cooling container: *******************************************

  G4double PosZ = -25.3*cm; // extra z-pos to correspond with height in lab

  G4double jacketHeight       = 680.0*mm;
  G4double jacketflangeHeight = 53.0*mm;
  G4double LN2jacketMetalThick = 2.0*mm;


  // outer vacuum jacket volume: stainless steel ****************************

  G4double jacketRadius     = 127.5*mm;
  //  G4double jacketHeight     = 680.0*mm; // defined above to get full-height
  G4double jacketMetalThick = LN2jacketMetalThick;

  G4Tubs* jacket_tube = new G4Tubs("jacket_tube",
     0.*cm, jacketRadius, 0.5*jacketHeight, 0.*deg, 360.*deg);
  jacket_log  = new G4LogicalVolume(jacket_tube, jacket_mat, "jacket_log");
  jacket_phys = new G4PVPlacement(0, G4ThreeVector(0.,0.,PosZ),
     "jacket_phys", jacket_log, lab_phys, false,0);

  G4VisAttributes* jacket_vat = new G4VisAttributes(grey);
  // jacket_log->SetVisAttributes(G4VisAttributes::GetInvisible());
  jacket_log->SetVisAttributes(jacket_vat);


  // outer vacuum jacket flanges: stainless steel *************************

  G4double jacketflangeRadius  = 127.5*mm;
  // G4double jacketflangeHeight = 53.0*mm; // defined above to get full-height
  G4double topjacketflangePosZ = 0.5*(jacketHeight+jacketflangeHeight);

  G4Tubs* jacketflange_tube = new G4Tubs("jacketflange_tube",
     0.*cm, jacketflangeRadius, 0.5*jacketflangeHeight, 0.*deg, 360.*deg);
  jacketflange_log     = new G4LogicalVolume
    (jacketflange_tube, jacketflange_mat, "jacketflange_log");
  topjacketflange_phys = new G4PVPlacement
    (0, G4ThreeVector(0.,0.,topjacketflangePosZ + PosZ),
     "topjacketflange_phys", jacketflange_log, lab_phys, false,0);
  bottomjacketflange_phys = new G4PVPlacement
    (0, G4ThreeVector(0.,0.,-topjacketflangePosZ + PosZ),
     "bottomjacketflange_phys", jacketflange_log, lab_phys, false,0);

  // jacketflange_log->SetVisAttributes(G4VisAttributes::GetInvisible());
  jacketflange_log->SetVisAttributes(jacket_vat);

  // vacuum **************************************************************

  G4double vacuumRadius = jacketRadius - jacketMetalThick;
  G4double vacuumHeight = jacketHeight - jacketMetalThick;

  G4Tubs* vacuum_tube = new G4Tubs("vacuum_tube",
     0.*cm, vacuumRadius, 0.5*vacuumHeight, 0.*deg, 360.*deg);
  vacuum_log  = new G4LogicalVolume(vacuum_tube, vacuum_mat, "vacuum_log");
  vacuum_phys = new G4PVPlacement(0, G4ThreeVector(0.,0.,0.),
     "vacuum_phys", vacuum_log, jacket_phys, false,0);

  // G4VisAttributes* vacuum_vat= new G4VisAttributes(lgrey);
  vacuum_log->SetVisAttributes(G4VisAttributes::GetInvisible());


  // copper cooling jacket volume: **************************************

  G4double vesselHeight     = 320.0*mm;

  // inner vessel jacket volume: stainless steel ************************

  //  G4double vesselHeight = 320.0*mm; // - moved earlier
  G4double vesselMetalThick      = jacketMetalThick;
  G4double vesselRadius          = 75.0*mm + vesselMetalThick;
  G4double vesselflangeRadius    = 101.5*mm;
  G4double vesselflangeThick     = 40.0*mm;
  G4double PMTvesselHeight       = 152.0*mm;
  G4double vesselVPos            = 7.0*cm;
  G4double TotalvesselHeight     = PMTvesselHeight + vesselHeight;

  G4Tubs* vessel_tube    = new G4Tubs("vessel_tube",
     0.*cm, vesselRadius, 0.5*vesselHeight, 0.*deg, 360.*deg);

  vessel_log  = new G4LogicalVolume(vessel_tube, vessel_mat, "vessel_log");
  vessel_phys = new G4PVPlacement(0, G4ThreeVector(0.,0.,vesselVPos),
     "vessel_phys", vessel_log, vacuum_phys, false,0);


  // flanges: 1=upper half (diff. inner diam.) 2=lower half
  G4Tubs* vesseltop_flange1 = new G4Tubs("vesseltop_flange1",
     0.*cm, vesselflangeRadius, 0.25*vesselflangeThick, 0.*deg, 360.*deg);
  vesseltop_log1  = new G4LogicalVolume
    (vesseltop_flange1, vessel_mat, "vesseltop_log1");
  vesseltop_phys1 = new G4PVPlacement
    (0,
     G4ThreeVector(0.,0.,0.5*(vesselHeight+0.5*vesselflangeThick)+vesselVPos),
     "vesseltop_phys1", vesseltop_log1, vacuum_phys, false,0);

  G4Tubs* vesseltop_flange2 = new G4Tubs("vesseltop_flange2",vesselRadius,
    vesselflangeRadius, 0.25*vesselflangeThick, 0.*deg, 360.*deg);
  vesseltop_log2  = new G4LogicalVolume
    (vesseltop_flange2, vessel_mat, "vesseltop_log2");
  vesseltop_phys2 = new G4PVPlacement
    (0,
     G4ThreeVector(0.,0.,0.5*(vesselHeight-0.5*vesselflangeThick)+vesselVPos),
     "vesseltop_phys2", vesseltop_log2, vacuum_phys, false,0);


  G4Tubs* vesselbottom_flange1 = new G4Tubs
    ("vesselbottom_flange1",vesselRadius, vesselflangeRadius,
     0.25*vesselflangeThick, 0.*deg, 360.*deg);
  vesselbottom_log1  = new G4LogicalVolume
    (vesselbottom_flange1, vessel_mat, "vesselbottom_log1");
  vesselbottom_phys1 = new G4PVPlacement(0,
     G4ThreeVector(0.,0.,-0.5*(vesselHeight-0.5*vesselflangeThick)+vesselVPos),
     "vesselbottom_phys1", vesselbottom_log1, vacuum_phys, false,0);

  G4Tubs* vesselbottom_flange2 = new G4Tubs
    ("vesselbottom_flange2",0.*cm, vesselflangeRadius,
     0.25*vesselflangeThick, 0.*deg, 360.*deg);
  vesselbottom_log2  = new G4LogicalVolume
    (vesselbottom_flange2, vessel_mat, "vesselbottom_log2");
  vesselbottom_phys2 = new G4PVPlacement(0,
     G4ThreeVector(0.,0.,-0.5*(vesselHeight+0.5*vesselflangeThick)+vesselVPos),
     "vesselbottom_phys2", vesselbottom_log2, vacuum_phys, false,0);


  G4VisAttributes* vessel_vat     = new G4VisAttributes(grey);
  //  vessel_log->SetVisAttributes(G4VisAttributes::GetInvisible());
   vessel_vat->SetForceSolid(true);
  vessel_log->SetVisAttributes(vessel_vat);
  vesseltop_log1->SetVisAttributes(vessel_vat);
  vesselbottom_log1->SetVisAttributes(vessel_vat);
  vesseltop_log2->SetVisAttributes(vessel_vat);
  vesselbottom_log2->SetVisAttributes(vessel_vat);



  // *********************************************************************
  // LXe height = 15.75 mm, gXe height = 6.00 mm
  // NB: Increased liquid height by 1mm - to take away problem with
  // over-lapping volumes/ring pronounced from liquid phase..........
  // *********************************************************************

  // detector volume: gas phase ******************************************

  G4double mirrorVPos     = 21.3*cm;
  G4double gasGap         = 6.0*mm;
  G4double DetectorRadius = vesselRadius - vesselMetalThick;
  G4double GXeHeight      = TotalvesselHeight - mirrorVPos + gasGap;
  G4double GXeVPos        = 0.5*vesselHeight - 0.5*GXeHeight;

  G4Tubs* GXe_tube = new G4Tubs("GXe_tube",
     0.*cm, DetectorRadius, 0.5*GXeHeight, 0.*deg, 360.*deg);
  GXe_log  = new G4LogicalVolume(GXe_tube, GXe_mat, "GXe_log");
  GXe_phys = new G4PVPlacement(0, G4ThreeVector(0.,0.,GXeVPos),
     "GXe_phys", GXe_log, vessel_phys, false,0);

  G4VisAttributes* GXe_vat = new G4VisAttributes(cyan);
  GXe_vat->SetForceSolid(true);
  GXe_vat->SetVisibility(true);
  GXe_log->SetVisAttributes(GXe_vat);


  // liquid phase *******************************************************

  G4double LXeHeight         = mirrorVPos - gasGap;
  G4double LXeTubeHeight     = LXeHeight - PMTvesselHeight;
  G4double LXeVPos           = -0.5*TotalvesselHeight + 0.5*LXeHeight;

  G4Tubs* LXe_tube = new G4Tubs("GXe_tube",
     0.*cm, DetectorRadius, 0.5*LXeTubeHeight, 0.*deg, 360.*deg);

  LXe_log  = new G4LogicalVolume(LXe_tube, LXe_mat, "LXe_log");
  LXe_phys = new G4PVPlacement(0, G4ThreeVector(0.*cm, 0.*cm, LXeVPos),
    "LXe_phys", LXe_log, vessel_phys, false, 0);

  // attributes
  G4VisAttributes* LXe_vat = new G4VisAttributes(lblue);
  LXe_vat->SetForceSolid(true);
  LXe_vat->SetVisibility(true);
  LXe_log->SetVisAttributes(LXe_vat);


  // Gas phase vessel lagging - for optical properties:

  G4double laggingThickness = 10.*micrometer;
  G4double laggingRadius    = DetectorRadius - laggingThickness;

  G4Tubs* gaslag_tube = new G4Tubs("gaslag_tube", laggingRadius,
     DetectorRadius, 0.5*GXeHeight, 0.*deg, 360.*deg);
  gaslag_log  = new G4LogicalVolume(gaslag_tube, vessel_mat, "gaslag_log");
  gaslag_phys = new G4PVPlacement(0, G4ThreeVector(0.*cm, 0.*cm, 0.*cm),
    "gaslag_phys", gaslag_log, GXe_phys, false, 0);

  // attributes
  G4VisAttributes* gaslag_vat = new G4VisAttributes(lgreen);
  // gaslag_vat->SetForceSolid(true);
  gaslag_vat->SetVisibility(true);
  gaslag_log->SetVisAttributes(gaslag_vat);


  // liquid phase vessel lagging - for optical properties:

  G4double lagTubeRadius = DetectorRadius - laggingThickness;
  G4double lagTubeHeight = LXeHeight - PMTvesselHeight;

  G4Tubs* liqLag_tube = new G4Tubs("liqlag_tube", lagTubeRadius,
     DetectorRadius, 0.5*lagTubeHeight, 0.*deg, 360.*deg);

  liqLag_log  = new G4LogicalVolume(liqLag_tube, vessel_mat, "liqLag_log");
  liqLag_phys = new G4PVPlacement(0, G4ThreeVector(0.*cm, 0.*cm, 0.*cm),
    "liqLag_phys", liqLag_log, LXe_phys, false, 0);

  // attributes
  G4VisAttributes* liqLag_vat = new G4VisAttributes(magenta);
  // liqLag_vat->SetForceSolid(true);
  liqLag_vat->SetVisibility(true);
  liqLag_log->SetVisAttributes(liqLag_vat);


  // Vessel Wall Optical Surface definition:
  G4OpticalSurface* OpVesselSurface = new G4OpticalSurface
    ("VesselSurface", unified, polished, dielectric_metal);

  // created optical lagging onto vessel - to avoid clash between over-lapping
  // liquid and gas phase - so removed below:
  /*
  G4LogicalBorderSurface* VesselSurface;
  VesselSurface = new G4LogicalBorderSurface
    ("Vessel", liqPhase_phys, vessel_phys, OpVesselSurface);
  */

  std::vector<G4double> vessel_PP   = { 6.5*eV, 7.50*eV };
  std::vector<G4double> vessel_REFL = { 0.2, 0.2 };
  G4MaterialPropertiesTable* vessel_mt = new G4MaterialPropertiesTable();
  vessel_mt->AddProperty("REFLECTIVITY", vessel_PP, vessel_REFL);
  OpVesselSurface->SetMaterialPropertiesTable(vessel_mt);

  // G4LogicalBorderSurface* VesselTopSurface =
  new G4LogicalBorderSurface
    ("VesselTop", GXe_phys, vesseltop_phys1, OpVesselSurface);

  //G4LogicalBorderSurface* VesselBottomSurface =
  new G4LogicalBorderSurface
    ("VesselBottom", LXe_phys, vesselbottom_phys2, OpVesselSurface);

  //G4LogicalBorderSurface* GasVesselSurface =
  new G4LogicalBorderSurface
    ("GasVessel", GXe_phys, gaslag_phys, OpVesselSurface);

  //G4LogicalBorderSurface* LiquidVesselSurface =
  new G4LogicalBorderSurface
    ("LiquidVessel", LXe_phys, liqLag_phys, OpVesselSurface);



  // ......................................................................
  // attach user limits ...................................................


  G4cout << G4endl << "User Limits: " << G4endl
	 << "\t theMaxTimeCuts:     " << G4BestUnit(theMaxTimeCuts,"Time")
	 << G4endl
	 << "\t theRoomTimeCut:     " << G4BestUnit(theRoomTimeCut,"Time")
	 << G4endl
	 << "\t theMaxStepSize:     " << G4BestUnit(theMaxStepSize,"Length")
	 << G4endl
	 << "\t theMinEKine:        " << G4BestUnit(theMinEkine,"Energy")
	 << G4endl
	 << "\t minRoomMinEKine:    " << G4BestUnit(theRoomMinEkine,"Energy")
	 << G4endl << G4endl;

  if (theUserLimitsForRoom != 0) delete theUserLimitsForRoom;
  if (theUserLimitsForDetector != 0) delete theUserLimitsForDetector;

  theUserLimitsForRoom = new G4UserLimits(theMaxStepSize,   // step length max
					  DBL_MAX,          // track length max
					  theRoomTimeCut,   // Time cut
					  theRoomMinEkine); // min energy

  theUserLimitsForDetector = new G4UserLimits(theDetectorStepSize,
					      DBL_MAX, // Track Max
					      theMaxTimeCuts,
					      theMinEkine);

      world_log->SetUserLimits(theUserLimitsForRoom);
        lab_log->SetUserLimits(theUserLimitsForRoom);
     jacket_log->SetUserLimits(theUserLimitsForRoom);
     vacuum_log->SetUserLimits(theUserLimitsForRoom);
     vessel_log->SetUserLimits(theUserLimitsForRoom);
        GXe_log->SetUserLimits(theUserLimitsForDetector);
	//        LXe_log->SetUserLimits(theUserLimitsForXenon);
        LXe_log->SetUserLimits(theUserLimitsForDetector);

 return world_phys;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DMXDetectorConstruction::ConstructSDandField()
{
  // ......................................................................
  // sensitive detectors ..................................................
  // ......................................................................

  if (LXeSD.Get() == 0)
    {
      G4String name="/DMXDet/LXeSD";
      DMXScintSD* aSD = new DMXScintSD(name);
      LXeSD.Put(aSD);
    }
  G4SDManager::GetSDMpointer()->AddNewDetector(LXeSD.Get());
  if (LXe_log)
    SetSensitiveDetector(LXe_log,LXeSD.Get());

  return;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// specific method to G4UserLimits:= SetUserMinEkine
void DMXDetectorConstruction::SetRoomEnergyCut(G4double val)
{
  // set minimum charged particle energy cut - NB: for ROOM
  theRoomMinEkine = val;
  if (theUserLimitsForRoom != 0)
    {
      theUserLimitsForRoom->SetUserMinEkine(val);
      G4cout << " Changing Room energy cut to: " << G4BestUnit(val,"Energy")
	     << G4endl;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// specific method to G4UserLimits:= SetUserMinEkine
void DMXDetectorConstruction::SetEnergyCut(G4double val)
{
  // set minimum charged particle energy cut - NB: for Xenon Detector
  theMinEkine = val;
  if (theUserLimitsForDetector != 0)
    {
      theUserLimitsForDetector->SetUserMinEkine(val);
      G4cout << "Changing Detector energy cut to: " << G4BestUnit(val,"Energy")
	     << G4endl;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// specific method to G4UserLimits:= SetUserMaxTime
void DMXDetectorConstruction::SetRoomTimeCut(G4double val)
{
  // set room time cut:
  theRoomTimeCut = val;
  if (theUserLimitsForRoom != 0)
    {
      theUserLimitsForRoom->SetUserMaxTime(val);
      G4cout << " Changing Room Time cut to: " << G4BestUnit(val,"Time")
	     << G4endl;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// specific method to G4UserLimits:= SetUserMaxTime
void DMXDetectorConstruction::SetTimeCut(G4double val)
{
  // set detector time cut:
  theMaxTimeCuts = val;
  if (theUserLimitsForDetector != 0)
    {
      theUserLimitsForDetector->SetUserMaxTime(val);
      G4cout << " Changing Detector Time cut to: " << G4BestUnit(val,"Time")
	     << G4endl;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
