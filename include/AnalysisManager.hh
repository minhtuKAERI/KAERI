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
#ifndef MOLECULAR_ANALYSIS_MANAGER_HH
#define MOLECULAR_ANALYSIS_MANAGER_HH

#include "globals.hh"
#include "G4RootAnalysisManager.hh"
#include "G4AnalysisManager.hh"
#include "G4ThreeVector.hh"
#include "G4OH.hh"
#include "G4Electron_aq.hh"
#include "G4Hydrogen.hh"
//#include "MolecularList.hh"
#include <cstdlib>
#include <ctime>
#include <map>
#include <vector>
#include <random>

class ChromosomeHit;

class DNAHit;

class AnalysisMessenger;

class DNAGeometry;




struct BasePairDamageRecord
{
  G4double fStrand1Energy    = 0;
  G4double fStrand2Energy    = 0;
  G4double fBp1Energy        = 0;
  G4double fBp2Energy        = 0;
  G4bool fBp1IndirectDmg     = false;
  G4bool fBp2IndirectDmg     = false;
  G4bool fStrand1IndirectDmg = false;
  G4bool fStrand2IndirectDmg = false;
  G4bool fBp1IndirectEvt     = false;
  G4bool fBp2IndirectEvt     = false;
  G4bool fStrand1IndirectEvt = false;
  G4bool fStrand2IndirectEvt = false;
  G4bool fbp1DirectDmg       = false;
  G4bool fbp2DirectDmg       = false;
  G4bool fbp1DirectOx        = false;
  G4bool fbp2DirectOx        = false;
  G4bool fStrand1DirectDmg   = false;
  G4bool fStrand2DirectDmg   = false;
  G4bool fbp1InducedBreak    = false;
  G4bool fbp2InducedBreak    = false;
  G4int  fbp1Product         = 0;   
  G4int  fbp2Product         = 0;
  G4int  fbp1ProductOx         = 0;   
  G4int  fbp2ProductOx         = 0;
  G4int  fBaseName           = 0;   
/*
fbp*Product value
0 = no base damage
1 = other base damage
2 = thymine glycol (OH)
8 = 8-oxo-guanine (OH)
21 = tg (direct)
81 = 8og (direct)

fBaseName
  NSPECIFIED  = 0,
  SUGAR       = 1,
  PHOSPHATE   = 2,
  GUANINE     = 3,
  ADENINE     = 4,
  CYTOSINE    = 5,
  THYMINE     = 6,
  Count       = 7


*/
  
};

//exp
enum complexityEnum
{
  SSB,
  SSBplus,
  twoSSB,
  DSB,
  DSBplus,
  DSBplusplus,
  NoneComplexity
};

enum sourceEnum
{
  SSBd,
  SSBi,
  SSBm,
  DSBh,
  DSBm,
  DSBd,
  DSBi,
  undefined
};
//

struct DamageClassification
{
  complexityEnum fComplexity  = NoneComplexity;
  sourceEnum fSource      = undefined;
  G4int fbaseDmg        = 0;
  G4int fStrandDmg      = 0;
  G4int fDirectBreaks   = 0;
  G4int fIndirectBreaks = 0;
  G4int fInducedBreaks  = 0;
};

class DamageRecord
{
 public:
  DamageRecord(const G4String&, int64_t, G4int, G4int);  // dousatsu
  ~DamageRecord();

  void AddBasePairDamage(BasePairDamageRecord* bp, const G4ThreeVector& pos)
  {
    fDamageRecords.push_back(bp);
    fPositions.push_back(pos);
  };

  void AddEmptyBPDamage(int64_t ii);

  void AddStrandHit(const G4MoleculeDefinition* mol);

  void AddBaseHit(const G4MoleculeDefinition* mol);
  
  void AddProductBaseBp1 (G4int ii) 
  { 
    if (ii==0 && (fProductBaseBp1==81 || fProductBaseBp1==21) ) return;
    if (ii==81 || ii==21) {fProductBaseBp1 = ii; return; }
    else {fProductBaseBp1 = ii;}
    //if (ii==0) return;
    //if (fProductBaseBp1 == 0) {fProductBaseBp1 = ii; }
  };

  void AddProductBaseBp2 (G4int ii) 
  { 
    if (ii==0 && (fProductBaseBp2==81 || fProductBaseBp2==21) ) return;
    if (ii==81 || ii==21) {fProductBaseBp2 = ii; return; }
    else {fProductBaseBp2 = ii;}
    // if (ii==0) return;
    // if (fProductBaseBp2 == 0) {fProductBaseBp2 = ii; }
  };

  void AddProductBaseBp1Ox (G4int ii) { fProductBaseBp1Ox = ii; };
  
  void AddProductBaseBp2Ox (G4int ii) { fProductBaseBp2Ox = ii; };

  void AddProductBaseName (G4int ii) { fProductBaseName = ii; };

  void PrintRecord(const G4String&, const G4double& dsbDistance = 10);

  DamageClassification* GetClassification(const G4double& dsbDistance = 10);

  inline const G4String& GetName() const { return fName; };

  int64_t GetSize() const
  {
    return fDamageRecords.size();
  };  // dousatsu

  inline G4int GetOHBaseHits() const { return fOHBase; };

  inline G4int GetEaqBaseHits() const { return fEaqBase; };

  inline G4int GetHBaseHits() const { return fHBase; };

  inline G4int GetOHStrandHits() const { return fOHStrand; };

  inline G4int GetEaqStrandHits() const { return fEaqStrand; };

  inline G4int GetHStrandHits() const { return fHStrand; };

  inline G4int GetPlacementIdx() const { return fStartPlacement; };

  inline G4int GetChainIdx() const { return fChainIdx; };

  inline G4int GetStrandIdx() const { return fStrandIdx; };
  
  inline G4int GetProductBaseBp1() const { return fProductBaseBp1; };
  
  inline G4int GetProductBaseBp2() const { return fProductBaseBp2; };

  inline G4int GetProductBaseBp1Ox() const { return fProductBaseBp1Ox; };
  
  inline G4int GetProductBaseBp2Ox() const { return fProductBaseBp2Ox; };

  inline G4int GetProductBaseName() const { return fProductBaseName; };

  inline int64_t GetStartBPIdx() const
  {
    return fStartIndex;
  };  // dousatsu
  void AddTestDamage(G4int, G4int, G4int, G4int);

  G4ThreeVector GetMeanPosition() const;

  G4double GetMeanDistance() const;

  G4double GetEnergy() const;

 private:
  G4String fName;
  int64_t fStartIndex;
  G4int fStartPlacement, fChainIdx, fStrandIdx = 0;
  G4int fOHBase = 0, fOHStrand = 0, fHBase = 0, fHStrand = 0, fEaqBase = 0,
        fEaqStrand = 0, fProductBaseBp1 = 0, fProductBaseBp2 = 0, fProductBaseBp1Ox = 0, fProductBaseBp2Ox = 0  , fProductBaseName = 0;
  std::vector<BasePairDamageRecord*> fDamageRecords;
  std::vector<G4ThreeVector> fPositions;
  
    const G4MoleculeDefinition* fOH = G4OH::Definition();
const G4MoleculeDefinition* fe_aq = G4Electron_aq::Definition();
const G4MoleculeDefinition* fH = G4Hydrogen::Definition();
  




  static const char* fDirectDamageChar;
  static const char* fIndirectDamageChar;
  static const char* fHitNoDamageChar;
  static const char* fNotHitChar;
  static const char* fBothDamageChar;

  const char* GetChar(const G4bool&, const G4bool&, const G4double&);
};

struct Node
{
  int64_t fkey;
  DNAHit* fdata;
  Node* fleft;
  Node* fright;
  Node* fparent;
};

// A binary tree is used to order the DNA hit objects
// It creates internal copies of all DNA Hits passed to it and then
// deletes them
class BinaryTree
{
 public:
  BinaryTree();

  virtual ~BinaryTree();

  void Insert(const DNAHit*);

  DNAHit* Search(int64_t);

  void Destroy_tree();

  // return left-most node
  DNAHit* First() const;

  // Return next node with higher key
  DNAHit* Next(const DNAHit*) const;

 private:
  void Destroy_tree_(Node*);

  void Insert_(DNAHit*, Node*);

  DNAHit* Search_(int64_t, Node*);

  DNAHit* First_(Node*) const;

  DNAHit* Next_(int64_t, Node*) const;

  Node* fRoot;
};

class AnalysisManager
{
 public:
  AnalysisManager();

  virtual ~AnalysisManager();

  void Initialize();

  void ProcessDNAHitsVector(const std::vector<const DNAHit*>&);//,const std::map<G4int, std::set<G4int>>& oxidations);

  void ProcessChromosomeHitMap(const std::map<uint32_t, ChromosomeHit*>&);

  void ProcessPrimary(const G4ThreeVector&, const G4double&, const G4double&);

  void ProcessCellEdep(const G4double&);  // dousatsu
  void Close();

  inline void SetSaveStrands(const G4bool strand) { fSaveStrands = strand; };

  inline void SetStrandDirectory(const G4String& dir) { fStrandDirectory =
      dir; };

  inline void SetFragmentGap(G4int gap) { fFragmentGap = gap; };

  inline void SetDSBDistance(G4int gap) { fDSBDistance = gap; };

  inline void SetChainToSave(G4int i) { fChainToSave = i; };

  inline void SetFileName(const G4String& name) { fFileName = name; };

  void TestClassification();
  
   // set up likelihood of base lesions
  
  

 private:
 
  G4AnalysisManager* fAnalysisManager;
  G4bool fSaveStrands = false;
  G4String fStrandDirectory = "./";
  G4String fFileName = "molecular-dna";
  G4int fFragmentGap = 100;
  G4int fDSBDistance = 10;
  G4int fChainToSave = -1;
  DNAGeometry* fpDNAGeometry = nullptr;
  AnalysisMessenger* fpAnalysisMessenger;
  G4double prob_val = 1 ;
  
  const G4MoleculeDefinition* fOH = G4OH::Definition();
  const G4MoleculeDefinition* fe_aq = G4Electron_aq::Definition();
  const G4MoleculeDefinition* fH = G4Hydrogen::Definition();

  const G4double PROB_8OG = 0.183;
  const G4double PROB_TG = 0.9;
  const G4double PROB_DHT = 0.187;  // add your probability
  const G4double PROB_DHC = 0.209;  // add your probability
};

#endif  // MOLECULAR_ANALYSIS_MANAGER_HH
