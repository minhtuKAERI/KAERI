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
/// file:SteppingAction.cc
/// brief:
#include "SteppingAction.hh"
#include "EventAction.hh"
#include "DetectorConstruction.hh"
#include "DNAGeometry.hh"
#include "DNAHit.hh"
#include "UtilityFunctions.hh"
#include "ChromosomeMapper.hh"
#include "OctreeNode.hh"
#include "PlacementVolumeInfo.hh"

#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4AffineTransform.hh"
#include "G4StepPoint.hh"
#include "G4VTouchable.hh"
#include "Randomize.hh"
#include "G4TouchableHistory.hh"

#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(EventAction* evt)
  : fEventAction(evt)
{
  const auto* det = dynamic_cast<const DetectorConstruction*>(
    G4RunManager::GetRunManager()->GetUserDetectorConstruction());

  fDNAGeometry = det->GetDNAGeometry();

  fEventAction = (EventAction*)G4RunManager::GetRunManager()->GetUserEventAction();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
  G4double edep = aStep->GetTotalEnergyDeposit();
  if(edep == 0)
  {
    return;
  }

  // Work out if we are in the Logical volume assigned to DNA
  const G4VTouchable* touch = aStep->GetPreStepPoint()->GetTouchable();

  const G4VProcess* process = aStep->GetPostStepPoint()->GetProcessDefinedStep();
  if(!process) { return; }
  const std::vector<const G4Track*>* secondaries = aStep->GetSecondaryInCurrentStep();

  G4String processName = process->GetProcessName();

  G4bool directIonisation = false;

  // remember the chem volume is the "real" volume.
  G4LogicalVolume* dnaVol = fDNAGeometry->GetDNAChemVolumePointer();
  G4int depth             = touch->GetHistoryDepth();
  G4bool isInDNARegion    = false;

  if(depth >= 0)
  {
    fEventAction->AddCellEdep(edep);
  }

  while((depth >= 0) && !isInDNARegion)
  {
    if(touch->GetVolume(depth) != nullptr)
    {
      isInDNARegion = isInDNARegion ||
                      (touch->GetVolume(depth)->GetLogicalVolume() == dnaVol);
    }
    else
    {
      G4cout << "Null pointer volume in mother with depth of "
             << touch->GetHistoryDepth() << " looking for volume at depth "
             << depth << G4endl;
      assert(depth < 0);
    }
    depth--;
  }

  if(isInDNARegion)
  {
    // IMPORTANT: Pass directIonisation to ProcessDNARegionHit
    // This flag will be applied to ALL DNA hits found in this step
    // if(processName == "GenericIon_G4DNAIonisation" || processName == "e-_G4DNAIonisation")
    // {
    //   //G4cout<<"process name: "<<processName<<G4endl;
    //   // const G4Material* mat =aStep->GetPostStepPoint()->GetMaterial();
    //   // const G4String& name = mat->GetName();
    //   // G4cout<<"base ionized "<<name<<G4endl;
    //   directIonisation = true;
    // }
    // else if (secondaries) 
    // {
    //   for (const G4Track* sec: *secondaries)
    //   {
    //     const G4VProcess* creatorPoc = sec->GetCreatorProcess();
    //     if (creatorPoc) 
    //     {
    //       const G4String& spname = creatorPoc->GetProcessName();
    //       //G4cout<<"secondaries process: "<<spname<<G4endl;
    //       if(spname.find("G4DNAIonisation")!=std::string::npos)
    //       {
    //         directIonisation = true;
    //       }
    //     }
    //   }
    // }

    ProcessDNARegionHit(aStep, false);
    AccumulateEnergyInBases(aStep, edep);
  }

  G4ThreeVector prepos = aStep->GetPreStepPoint()->GetPosition();
  if(fDNAGeometry->GetChromosomeMapper()->GetChromosome(prepos) != nullptr)
  {
    DoChromosomeRegionHit(prepos, edep, false);
  }

  if(aStep->GetTrack()->GetTrackID() == 1)
  {
    fEventAction->SetPrimStopPos(prepos);
    G4double tl = aStep->GetStepLength();
    fEventAction->AddTrackLengthCell(tl);
    if(fDNAGeometry->GetChromosomeMapper()->GetChromosome(prepos) != nullptr)
    {
      fEventAction->AddTrackLengthChro(tl);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::ProcessDNARegionHit(const G4Step* aStep, G4bool directIonisation)
{
  G4double edep             = aStep->GetTotalEnergyDeposit();
  const G4VTouchable* touch = aStep->GetPreStepPoint()->GetTouchable();
  G4ThreeVector pos =
    aStep->GetPreStepPoint()->GetPosition() +
    G4UniformRand() * (aStep->GetPostStepPoint()->GetPosition() -
                       aStep->GetPreStepPoint()->GetPosition());
  G4VPhysicalVolume* pv = aStep->GetPreStepPoint()->GetPhysicalVolume();
  G4LogicalVolume* lv   = pv->GetLogicalVolume();
  //G4cout<<"processDNAREgion1 - directIon check DNAhit = "<<directIonisation<<G4endl;

  std::vector<G4VPhysicalVolume*> molvec;
  OctreeNode* parentNode;
  G4AffineTransform transform = touch->GetHistory()->GetTopTransform();
  parentNode                  = fDNAGeometry->GetTopOctreeNode(lv);
  if(parentNode == nullptr)
  {  // See if parent has an octree.
    parentNode = fDNAGeometry->GetTopOctreeNode(pv->GetMotherLogical());
    if(parentNode != nullptr)
    {  // If mother is a node, update it to be the PV also
      // These methods are unsafe if the pv is the top volume so they
      // need to be in this if statement.
      // lv = pv->GetMotherLogical();
      pv = touch->GetVolume(1);
      transform =
        touch->GetHistory()->GetTransform(touch->GetHistoryDepth() - 1);
    }
  }

  G4ThreeVector localPos = transform.TransformPoint(pos);

  if(parentNode != nullptr)
  {  // Get Molecule Vector
    molvec = parentNode->SearchOctree(
      localPos, fDNAGeometry->GetDirectInteractionRange());
  }

  // debug output for finding logical volumes.
  if(fDNAGeometry->GetVerbosity() > 1)
  {
    if(parentNode != nullptr)
    {
      G4String lvname = pv->GetLogicalVolume()->GetName();
      G4cout << "Found octree for logical volume: " << lvname << G4endl;
    }
    else
    {
      // Don't bother doing error printouts if the particle
      // is at the DNA mother volume level.
      // Remember the Chemistry Volume is the "real" world
      if(pv->GetLogicalVolume() != fDNAGeometry->GetDNAChemVolumePointer())
      {
        G4String lvname = pv->GetLogicalVolume()->GetName();
        G4String motherlvname;
        if(lv != nullptr)
        {
          motherlvname = lv->GetName();
        }
        G4cout << "Could not find octree for logical volume." << G4endl
               << "Particle in LV: " << lvname
               << ", with mother LV: " << motherlvname << G4endl;
      }
    }
  }

  //G4cout<<"processdnaregionhit molvec="<<molvec.empty()<<G4endl;
  if(!molvec.empty())
  {
    // NOTE: This loop can be optimised to avoid checking pairs twice
    G4VPhysicalVolume* closest = molvec[0];
    G4double dist              = (closest->GetTranslation() - localPos).mag();
    if(molvec.size() > 1)
    {
      G4double newdist;
      for(auto it = molvec.begin() + 1; it != molvec.end(); it++)
      {
        newdist = ((*it)->GetTranslation() - localPos).mag();
        if(newdist < dist)
        {
          dist    = newdist;
          closest = (*it);
        }
      }
    }
    //G4cout<<"dist = "<<dist/angstrom<<G4endl;
    G4bool isDNAHit = (dist < fDNAGeometry->GetDirectInteractionRange());
    //G4cout<<"stepping action - dnahit = "<<isDNAHit<<G4endl;
    //G4cout<<"processdnaregionhit DI="<<directIonisation<<G4endl;
    if(isDNAHit)
    {
      //G4cout<<"processDNAREgion2 - directIon check DNAhit = "<<directIonisation<<G4endl;
      DoChromosomeDNAHit(pos, localPos, edep, dist, closest->GetName(),
                         pv->GetName(), directIonisation);
      DoChromosomeRegionHit(pos, edep, true);
    }
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::DoChromosomeRegionHit(const G4ThreeVector& pos,
                                           const G4double& edep,
                                           const G4bool& dnahit)
{
  // TODO: chromosome should probably be based on position of center of
  //       parent volume, not of current pos
  G4String chromosome =
    fDNAGeometry->GetChromosomeMapper()->GetCurrentChromosomeKey(pos);
  if(!(chromosome.empty()))
  {
    fEventAction->AddChromosomeEdep(chromosome, edep, dnahit);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::DoChromosomeDNAHit(
  const G4ThreeVector& pos, const G4ThreeVector& localPos, const G4double& edep,
  const G4double& dist, const G4String& pvName, const G4String& motherName, 
  const G4bool& directIonisation)
{
  G4int placementIdx = std::stoi(utility::Split(motherName, '-').at(1));
  std::vector<G4String> pvNameVec = utility::Split(pvName, '-');
  molecule mol = utility::GetMoleculeEnum(pvNameVec.at(0));
  G4int chainIdx = std::stoi(pvNameVec.at(1));
  G4int strandIdx = std::stoi(pvNameVec.at(2));
  int64_t baseIdx = std::stoll(pvNameVec.at(3));

  // Convert to local chain and base index
  chainIdx = fDNAGeometry->GetGlobalChain(placementIdx, chainIdx);
  baseIdx += fDNAGeometry->GetStartIdx(placementIdx, chainIdx);

  if(fDNAGeometry->GetStrandsFlipped(placementIdx))
  {
      strandIdx = (strandIdx + 1) % 2;
  }
  
  G4String chromosome =
      fDNAGeometry->GetChromosomeMapper()->GetCurrentChromosomeKey(pos);

  // NEW: Accumulate energy for this base
  //BaseKey baseKey = {placementIdx, chainIdx, strandIdx, baseIdx};
  //fEventAction->AccumulateBaseEnergy(baseKey, edep);

  // EXISTING: CREATE DNAHit
  const DNAHit* dnaHit =
      new DNAHit(mol, placementIdx, chainIdx, strandIdx, baseIdx, pos, localPos,
                  edep, dist, chromosome, directIonisation);

  // RECORD DIRECT IONIZATION if flag is set
  // if(directIonisation)
  // {
  //   dnaHit->SetDirectIonisation(true);
  //   fEventAction->AddDirectBaseOxidationEvent(baseIdx, strandIdx);
    
  //   if(fDNAGeometry->GetVerbosity() > 0)
  //   {
  //     G4cout << "*** DIRECT IONIZATION HIT RECORDED ***" << G4endl;
  //     G4cout << "    BP: " << baseIdx << " Strand: " << strandIdx << G4endl;
  //   }
  // }

  fEventAction->AddDNAHit(dnaHit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


// 3. NEW METHOD: Accumulate energy in all nearby bases . Find all bases within 0.8 nm and tell Event action to
// accumulate energy for each 
void SteppingAction::AccumulateEnergyInBases(const G4Step* aStep, G4double edep)
{
    const G4double oxidationRadius = 0.8 * nm;
    
    G4ThreeVector pos = aStep->GetPreStepPoint()->GetPosition() +
        G4UniformRand() * (aStep->GetPostStepPoint()->GetPosition() -
                          aStep->GetPreStepPoint()->GetPosition());
    
    const G4VTouchable* touch = aStep->GetPreStepPoint()->GetTouchable();
    G4VPhysicalVolume* pv = aStep->GetPreStepPoint()->GetPhysicalVolume();
    G4LogicalVolume* lv = pv->GetLogicalVolume();
    
    // Get octree and transform (same logic as ProcessDNARegionHit)
    OctreeNode* parentNode;
    G4AffineTransform transform = touch->GetHistory()->GetTopTransform();
    parentNode = fDNAGeometry->GetTopOctreeNode(lv);
    
    if(parentNode == nullptr)
    {
        // See if parent has an octree
        parentNode = fDNAGeometry->GetTopOctreeNode(pv->GetMotherLogical());
        if(parentNode != nullptr)
        {
            pv = touch->GetVolume(1);
            transform = touch->GetHistory()->GetTransform(touch->GetHistoryDepth() - 1);
        }
    }
    
    if(parentNode == nullptr) return;
    
    G4ThreeVector localPos = transform.TransformPoint(pos);
    
    // Search for ALL molecules within oxidation radius
    std::vector<G4VPhysicalVolume*> molvec = 
        parentNode->SearchOctree(localPos, oxidationRadius);
    
    // Accumulate energy for all bases within range
    for(auto pv_mol : molvec)
    {
        G4double dist = (pv_mol->GetTranslation() - localPos).mag();
        if(dist < oxidationRadius)
        {
            // Extract base identifier (same logic as DoChromosomeDNAHit)
            G4String pvName = pv_mol->GetName();
            G4String motherName = pv->GetName();
            
            G4int placementIdx = std::stoi(utility::Split(motherName, '-').at(1));
            std::vector<G4String> pvNameVec = utility::Split(pvName, '-');
            
            G4int chainIdx = std::stoi(pvNameVec.at(1));
            G4int strandIdx = std::stoi(pvNameVec.at(2));
            int64_t baseIdx = std::stoll(pvNameVec.at(3));
            
            // Convert to global indices
            chainIdx = fDNAGeometry->GetGlobalChain(placementIdx, chainIdx);
            baseIdx += fDNAGeometry->GetStartIdx(placementIdx, chainIdx);
            
            if(fDNAGeometry->GetStrandsFlipped(placementIdx))
            {
                strandIdx = (strandIdx + 1) % 2;
            }
            
            BaseKey baseKey = {placementIdx, chainIdx, strandIdx, baseIdx};
            
            // Accumulate energy for this base
            fEventAction->AccumulateBaseEnergy(baseKey, edep);
        }
    }
}