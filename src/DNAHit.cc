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
/// \file DNAHit.cc
/// \brief Hit class for a hit interacting with a DNA molecule

#include <utility>

#include "DNAHit.hh"

G4ThreadLocal G4Allocator<DNAHit>* MolecularDNAHitAllocator = nullptr;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// PRIMARY CONSTRUCTOR: Physical hit (with energy deposit, no radical)
DNAHit::DNAHit(const molecule& mol, G4int placement_idx, G4int chain,
               G4int strand, int64_t bp, const G4ThreeVector& pos,
               const G4ThreeVector& localpos, G4double energy, G4double d,
               G4String chromo, G4bool directIonisation)
  : fMoleculeEnum(mol)
  , fPlacementIdx(placement_idx)
  , fChainIdx(chain)
  , fStrandIdx(strand)
  , fBasePairIdx(bp)
  , fPosition(pos)
  , fLocalPosition(localpos)
  , fEnergy(energy)
  , fDistance(d)
  , fChromosome(std::move(chromo))
  , fRadical(nullptr)
  , fDirectIonisation(directIonisation)
{
  // Computed quantities based on strand and molecule type
  if(fStrandIdx == 0)
  {
    if((fMoleculeEnum == SUGAR) || (fMoleculeEnum == PHOSPHATE))
    {
      fStrand1Rad = fRadical;
      fStrand1Energy = fEnergy;
    }
    else if((fMoleculeEnum == CYTOSINE) || (fMoleculeEnum == GUANINE) ||
            (fMoleculeEnum == ADENINE) || (fMoleculeEnum == THYMINE))
    {
      fBase1Rad = fRadical;
      fBP1Energy = fEnergy;
    }
    else
    {
      G4Exception("DNAHit", "ERR_UNKNOWN_MOLECULE", JustWarning,
                  "Physical hit with unknown molecule");
    }
  }
  else if(fStrandIdx == 1)
  {
    if((fMoleculeEnum == SUGAR) || (fMoleculeEnum == PHOSPHATE))
    {
      fStrand2Rad = fRadical;
      fStrand2Energy = fEnergy;
    }
    else if((fMoleculeEnum == CYTOSINE) || (fMoleculeEnum == GUANINE) ||
            (fMoleculeEnum == ADENINE) || (fMoleculeEnum == THYMINE))
    {
      fBase2Rad = fRadical;
      fBP2Energy = fEnergy;
    }
    else
    {
      G4Exception("DNAHit", "ERR_UNKNOWN_MOLECULE", JustWarning,
                  "Physical hit with unknown molecule");
    }
  }
  fDirectIonisation = directIonisation;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// SECONDARY CONSTRUCTOR: Chemical hit (with radical, no energy deposit)
DNAHit::DNAHit(const molecule& mol, G4int placement_idx, G4int chain,
               G4int strand, int64_t bp, const G4ThreeVector& pos,
               const G4ThreeVector& localpos, G4String chromo,
               const G4MolecularConfiguration* radical, G4bool directIonisation)
  : fMoleculeEnum(mol)
  , fPlacementIdx(placement_idx)
  , fChainIdx(chain)
  , fStrandIdx(strand)
  , fBasePairIdx(bp)
  , fPosition(pos)
  , fLocalPosition(localpos)
  , fEnergy(0)
  , fDistance(0)
  , fChromosome(std::move(chromo))
  , fRadical(radical)
  , fDirectIonisation(directIonisation)
{
  // Computed quantities based on strand and molecule type
  if(fStrandIdx == 0)
  {
    if((fMoleculeEnum == SUGAR) || (fMoleculeEnum == PHOSPHATE))
    {
      fStrand1Rad = fRadical;
    }
    else if((fMoleculeEnum == CYTOSINE) || (fMoleculeEnum == GUANINE) ||
            (fMoleculeEnum == ADENINE) || (fMoleculeEnum == THYMINE))
    {
      fBase1Rad = fRadical;
    }
    else
    {
      G4Exception("DNAHit", "ERR_UNKNOWN_MOLECULE", JustWarning,
                  "Chemical hit with unknown molecule");
    }
  }
  else if(fStrandIdx == 1)
  {
    if((fMoleculeEnum == SUGAR) || (fMoleculeEnum == PHOSPHATE))
    {
      fStrand2Rad = fRadical;
    }
    else if((fMoleculeEnum == CYTOSINE) || (fMoleculeEnum == GUANINE) ||
            (fMoleculeEnum == ADENINE) || (fMoleculeEnum == THYMINE))
    {
      fBase2Rad = fRadical;
    }
    else
    {
      G4Exception("DNAHit", "ERR_UNKNOWN_MOLECULE", JustWarning,
                  "Chemical hit with unknown molecule");
    }
  }
  fDirectIonisation = directIonisation;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// This method is used to combine the computed quantities of two hits.
// It keeps the settable parameters of the current hit, but the computed
// parameters come from both, letting the hit reflect what happened on each
// base pair.
void DNAHit::AddHit(const DNAHit& right)
{
  this->fStrand1Energy += right.GetStrand1Energy();
  this->fStrand2Energy += right.GetStrand2Energy();
  this->fBP1Energy += right.GetBP1Energy();
  this->fBP2Energy += right.GetBP2Energy();
  if(right.GetStrand1Rad() != nullptr)
  {
    this->fStrand1Rad = right.GetStrand1Rad();
  }
  if(right.GetBase1Rad() != nullptr)
  {
    this->fBase1Rad = right.GetBase1Rad();
  }
  if(right.GetStrand2Rad() != nullptr)
  {
    this->fStrand2Rad = right.GetStrand2Rad();
  }
  if(right.GetBase2Rad() != nullptr)
  {
    this->fBase2Rad = right.GetBase2Rad();
  }
//   // IMPORTANT: Preserve direct ionization flag
//   if(right.GetDirectIonisation())
//   {
//     this->fDirectIonisation = false;
//   }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DNAHit::DNAHit(const DNAHit& right)
{
  this->SetPlacementIdx(right.GetPlacementIdx());
  this->SetMolecule(right.GetMolecule());
  this->SetChainIdx(right.GetChainIdx());
  this->SetStrandIdx(right.GetStrandIdx());
  this->SetBasePairIdx(right.GetBasePairIdx());
  this->SetPosition(right.GetPosition());
  this->SetLocalPosition(right.GetLocalPosition());
  this->SetEnergy(right.GetEnergy());
  this->SetDistance(right.GetDistance());
  this->SetChromosome(right.GetChromosome());
  this->SetRadical(right.GetRadical());
  //this->SetDirectIonisation(right.GetDirectIonisation());

  // Computed Quantities, no setters.
  this->fStrand1Energy = right.GetStrand1Energy();
  this->fStrand2Energy = right.GetStrand2Energy();
  this->fBP1Energy = right.GetBP1Energy();
  this->fBP2Energy = right.GetBP2Energy();
  this->fStrand1Rad = right.GetStrand1Rad();
  this->fBase1Rad = right.GetBase1Rad();
  this->fStrand2Rad = right.GetStrand2Rad();
  this->fBase2Rad = right.GetBase2Rad();

  // IMPORTANT: Copy direct ionization flag
  this->fDirectIonisation = right.GetDirectIonisation();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const DNAHit& DNAHit::operator=(const DNAHit& right)
{
  this->SetMolecule(right.GetMolecule());
  this->SetPlacementIdx(right.GetPlacementIdx());
  this->SetChainIdx(right.GetChainIdx());
  this->SetStrandIdx(right.GetStrandIdx());
  this->SetBasePairIdx(right.GetBasePairIdx());
  this->SetPosition(right.GetPosition());
  this->SetLocalPosition(right.GetLocalPosition());
  this->SetEnergy(right.GetEnergy());
  this->SetDistance(right.GetDistance());
  this->SetChromosome(right.GetChromosome());
  this->SetRadical(right.GetRadical());

  this->fStrand1Energy = right.GetStrand1Energy();
  this->fStrand2Energy = right.GetStrand2Energy();
  this->fBP1Energy = right.GetBP1Energy();
  this->fBP2Energy = right.GetBP2Energy();
  this->fStrand1Rad = right.GetStrand1Rad();
  this->fBase1Rad = right.GetBase1Rad();
  this->fStrand2Rad = right.GetStrand2Rad();
  this->fBase2Rad = right.GetBase2Rad();

  // IMPORTANT: Copy direct ionization flag
  this->fDirectIonisation = right.GetDirectIonisation();

  return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int DNAHit::operator==(const DNAHit& right) const
{
  return (this == &right) ? 1 : 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......