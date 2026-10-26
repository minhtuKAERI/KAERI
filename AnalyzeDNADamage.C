#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <algorithm>
#include <cmath>
#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"

using namespace std;

// Structure to hold damage site information
struct DamageSite {
    int chain;
    int strand;
    int bp1;
    int bp2;
    int dirBreaks;
    int inBreaks;
    double basePair;
};

void AnalyzeDNADamage() {
    
    // ============================================================
    // Configuration Parameters
    // ============================================================
    //system ("hadd -O -f molecular-dna.root molecular-dna_t0.root molecular-dna_t1.root molecular-dna_t2.root molecular-dna_t3.root molecular-dna_t4.root molecular-dna_t5.root molecular-dna_t6.root molecular-dna_t7.root molecular-dna_t8.root molecular-dna_t9.root");//
    const char* inputFile  = "molecular-dna.root";
    const char* outputCSV = "damage_analysis.csv";
    
    int nCells = 1;
    double r3 = 900e-9 * 317e-9 * 900e-9;  // Chromosome volume (m³)
    double Nbp = 13.195594 * nCells;       // DNA length in Mbp
    double mass = 997 * 4 * 3.141592 * r3 / 3;  // Mass in kg
    double eVtoJ = 1.60218e-19;
    double EnergyDeposited_eV = 0;
    double acc_edep = 0;
    double acc_edep2 = 0;

    
    // ============================================================
    // Open ROOT file and get damage tree
    // ============================================================
    TFile* f = TFile::Open(inputFile);
    if (!f || f->IsZombie()) {
        cout << "Error: Cannot open file " << inputFile << endl;
        return;
    }
    
    TTree* tree = (TTree*)f->Get("tuples/damage");
    if (!tree) {
        cout << "Error: Cannot find tuples/damage tree" << endl;
        f->Close();
        return;
    }
    
    // ============================================================
    // Declare variables for tree branches
    // ============================================================
    Char_t Primary;
    Double_t Energy, ED;
    Int_t EB, ES, OHB, OHS, HB, HS;
    Int_t B1PD, B2PD, BD, BaseName;
    Double_t Position_x_um, Position_y_um, Position_z_um;
    Int_t Chain, Strand, DirectBreaks, IndirectBreaks;
    Double_t BasePair;
    char* type = new char[256];
    
    // ============================================================
    // Set branch addresses
    // ============================================================
    tree->SetBranchAddress("Primary", &Primary);
    tree->SetBranchAddress("Energy", &Energy);
    tree->SetBranchAddress("EaqBaseHits", &EB);
    tree->SetBranchAddress("EaqStrandHits", &ES);
    tree->SetBranchAddress("OHBaseHits", &OHB);
    tree->SetBranchAddress("OHStrandHits", &OHS);
    tree->SetBranchAddress("HBaseHits", &HB);
    tree->SetBranchAddress("HStrandHits", &HS);
    tree->SetBranchAddress("TypeClassification", type);
    tree->SetBranchAddress("Bp1Product", &B1PD);
    tree->SetBranchAddress("Bp2Product", &B2PD);
    tree->SetBranchAddress("Position_x_um", &Position_x_um);
    tree->SetBranchAddress("Position_y_um", &Position_y_um);
    tree->SetBranchAddress("Position_z_um", &Position_z_um);
    tree->SetBranchAddress("Chain", &Chain);
    tree->SetBranchAddress("Strand", &Strand);
    tree->SetBranchAddress("BasePair", &BasePair);
    tree->SetBranchAddress("BaseDamage", &BD);
    tree->SetBranchAddress("EnergyDeposited_eV", &ED);
    tree->SetBranchAddress("DirectBreaks", &DirectBreaks);
    tree->SetBranchAddress("IndirectBreaks", &IndirectBreaks);
    tree->SetBranchAddress("BaseName", &BaseName);
    
    // ============================================================
    // Initialize counters for individual base products
    // ============================================================
    int count_8og = 0;       // 8-oxo-G total
    int count_tg = 0;        // Thymine glycol total
    int count_81 = 0;        // Type 81 total
    int count_21 = 0;        // Type 21 total
    
    // Energy deposition for dose calculation
    double totalEnergyDeposit_eV = 0;
    
    // ============================================================
    // Cluster classification counters
    // ============================================================
    int cluster_ssb_8og = 0;           // SSB + 1x 8-oxo-G
    int cluster_ssb_tg = 0;            // SSB + 1x thymine glycol
    int cluster_ssb_both = 0;          // SSB + both lesions
    int cluster_2x8og = 0;             // 2x 8-oxo-G
    int cluster_multi8og = 0;          // >2x 8-oxo-G
    int cluster_2xtg = 0;              // 2x thymine glycol
    int cluster_multitg = 0;           // >2x thymine glycol
    int cluster_8og_tg = 0;            // 1x 8-oxo-G + 1x TG
    int cluster_multi_both = 0;        // Multiple 8-oxo-G + TG
    int cluster_multi_ssb_8og = 0;     // Multiple SSB + 8-oxo-G
    int cluster_multi_ssb_tg = 0;      // Multiple SSB + TG
    int cluster_multi_ssb_both = 0;    // Multiple SSB + both lesions

    // Counter for the DSB criterion: number of tree entries whose
    // TypeClassification is "DSB", "DSB+", or "DSB++" — these are excluded
    // from the SSB+lesion clustering pass since they're already known
    // double-strand breaks, not candidates for that classification. A
    // single DSB event typically produces two tree entries (one per
    // strand) a few bp apart, both labeled DSB*, so positions are
    // collected here and deduped by clustering within 10 bp after the
    // first pass rather than counting raw rows.
    vector<double> dsbPositions;
    
    // ============================================================
    // First pass: Count individual base products and energy
    // ============================================================
    vector<DamageSite> damageList;
    Long64_t nEntries = tree->GetEntries();
    
    cout << "Processing " << nEntries << " entries..." << endl;
    
    for (Long64_t i = 0; i < nEntries; i++) {
        tree->GetEntry(i);

        // Skip entries the tree already classifies as a double-strand
        // break (DSB / DSB+ / DSB++) — these aren't candidates for the
        // SSB+lesion clustering below, they're already known DSBs.
        string typeStr = (string)type;
        if (typeStr == "DSB" || typeStr == "DSB+" || typeStr == "DSB++") {
            dsbPositions.push_back(BasePair);
            continue;
        }
        
        // Count base products
        if (B1PD == 8) count_8og++;
        if (B2PD == 8) count_8og++;
        if (B1PD == 2) count_tg++;
        if (B2PD == 2) count_tg++;
        if (B1PD == 81) count_81++;
        if (B2PD == 81) count_81++;
        if (B1PD == 21) count_21++;
        if (B2PD == 21) count_21++;
        
        // Accumulate energy
        totalEnergyDeposit_eV += ED;
        
        // Store damage sites for clustering analysis
        if (B1PD > 0 || B2PD > 0 || DirectBreaks > 0 || IndirectBreaks > 0) {
            damageList.push_back({Chain, Strand, B1PD, B2PD, DirectBreaks, IndirectBreaks, BasePair});
        }
    }

    // Dedupe DSB positions: cluster entries within 10 bp of each other and
    // count each cluster once, so a DSB event with two nearby DSB-labeled
    // rows is only counted as one removed site.
    sort(dsbPositions.begin(), dsbPositions.end());
    int dsb_criterion_removed = 0;
    for (size_t i = 0; i < dsbPositions.size(); ) {
        dsb_criterion_removed++;
        size_t j = i + 1;
        while (j < dsbPositions.size() && abs(dsbPositions[j] - dsbPositions[i]) <= 10) j++;
        i = j;
    }
    
    // ============================================================
    // Calculate dose
    // ============================================================
    // Measure the Deposited Energy in the whole volume that includes DNA chain (chromosome)

    tree = (TTree *) f->Get("tuples/chromosome_hits");
    Long64_t nentries = tree->GetEntries();
    tree->SetBranchAddress("e_chromosome_kev",&EnergyDeposited_eV);
    nentries = tree->GetEntries();
    for(int i = 0;i<nentries;i++){
      tree->GetEntry(i);
      acc_edep += EnergyDeposited_eV *1e3;
      acc_edep2 += EnergyDeposited_eV *EnergyDeposited_eV *1e6;
    }
    tree->SetBranchAddress("e_dna_kev",&EnergyDeposited_eV);
    nentries = tree->GetEntries();
    for(int i = 0;i<nentries;i++){
      tree->GetEntry(i);
      acc_edep += EnergyDeposited_eV *1e3;
      acc_edep2 += EnergyDeposited_eV *EnergyDeposited_eV *1e6;
    }

    // Close the root file to free space
    //f->Close();

    // Calculate the absorbed dose
    double dose_Gy = acc_edep * eVtoJ / mass;
    
    // Read trees and leaves from root file, and give values to variables
    tree = (TTree*) f->Get("tuples/primary_source");
    Float_t number = (Float_t) tree->GetEntries();
    
    
    
    // ============================================================
    // Sort damage list by base pair position
    // ============================================================
    sort(damageList.begin(), damageList.end(), [](const DamageSite& a, const DamageSite& b) {
        return a.basePair < b.basePair;
    });
    
    // ============================================================
    // Second pass: Classify clusters (within 10 bp)
    // ============================================================
    set<int> visitedIndices;
    
    for (size_t i = 0; i < damageList.size(); i++) {
        if (visitedIndices.count(i)) continue;
        
        visitedIndices.insert(i);
        
        int breakCount = (damageList[i].dirBreaks > 0 ? 1 : 0) + (damageList[i].inBreaks > 0 ? 1 : 0);
        int c8ogCount = (damageList[i].bp1 == 8 || damageList[i].bp1 == 81 ? 1 : 0) + 
                        (damageList[i].bp2 == 8 || damageList[i].bp2 == 81 ? 1 : 0);
        int ctgCount = (damageList[i].bp1 == 2 || damageList[i].bp1 == 21 ? 1 : 0) + 
                       (damageList[i].bp2 == 2 || damageList[i].bp2 == 21 ? 1 : 0);
        
        // Look for nearby damage within 10 bp
        for (size_t j = i + 1; j < damageList.size(); j++) {
            if (abs(damageList[i].basePair - damageList[j].basePair) <= 10) {
                visitedIndices.insert(j);
                
                // Count 8-oxo-G
                if (damageList[j].bp1 == 8 || damageList[j].bp1 == 81) c8ogCount++;
                if (damageList[j].bp2 == 8 || damageList[j].bp2 == 81) c8ogCount++;
                
                // Count thymine glycol
                if (damageList[j].bp1 == 2 || damageList[j].bp1 == 21) ctgCount++;
                if (damageList[j].bp2 == 2 || damageList[j].bp2 == 21) ctgCount++;
                
                // Count breaks (DSB-classified entries were already
                // excluded from damageList in the first pass above via
                // TypeClassification, so anything remaining here is
                // treated as a candidate SSB)
                if (damageList[j].dirBreaks > 0 || damageList[j].inBreaks > 0) {
                    breakCount++;
                }
            } else {
                break; // Beyond 10 bp range
            }
        }
        
        // ============================================================
        // Classify the cluster
        // ============================================================
        if (breakCount > 1) {
            // Multiple SSBs
            if (c8ogCount > 0 && ctgCount > 0) cluster_multi_ssb_both++;
            else if (c8ogCount > 0) cluster_multi_ssb_8og++;
            else if (ctgCount > 0) cluster_multi_ssb_tg++;
        }
        else if (breakCount == 1) {
            // Single SSB
            if (c8ogCount == 1 && ctgCount == 0) cluster_ssb_8og++;
            else if (c8ogCount == 0 && ctgCount == 1) cluster_ssb_tg++;
            else if (c8ogCount > 0 && ctgCount > 0) cluster_ssb_both++;
        }
        else if (breakCount == 0) {
            // No SSB
            if (c8ogCount == 0 && ctgCount == 0) continue;
            else if (c8ogCount == 0 && ctgCount > 0) {
                if (ctgCount == 2) cluster_2xtg++;
                else if (ctgCount > 2) cluster_multitg++;
            }
            else if (c8ogCount > 0 && ctgCount == 0) {
                if (c8ogCount == 2) cluster_2x8og++;
                else if (c8ogCount > 2) cluster_multi8og++;
            }
            else if (c8ogCount == 1 && ctgCount == 1) cluster_8og_tg++;
            else if (c8ogCount > 1 || ctgCount > 1) cluster_multi_both++;
        }
    }
    
    // ============================================================
    // Write results to CSV file - Single header format
    // ============================================================
    
    // Check if file exists
    bool fileExists = false;
    ifstream checkFile(outputCSV);
    if (checkFile.good()) {
        fileExists = true;
    }
    checkFile.close();
    
    // Open file in append mode
    ofstream csvFile;
    if (fileExists) {
        csvFile.open(outputCSV, ios::app);  // Append mode
    } else {
        csvFile.open(outputCSV);  // Create new file
        // Write header only for new file
        csvFile << "num_bp,dose_Gy,num_par";
        csvFile << "8oxoG,TG,type81,type21,";
        csvFile << "cluster_SSB_8oxoG,cluster_SSB_TG,cluster_SSB_both,";
        csvFile << "cluster_2x8oxoG,cluster_multi8oxoG,cluster_2xTG,cluster_multiTG,";
        csvFile << "cluster_8oxoG_TG,cluster_multi_both,";
        csvFile << "cluster_multiSSB_8oxoG,cluster_multiSSB_TG,cluster_multiSSB_both,";
        csvFile << "dsb_criterion_removed\n";
    }
    
    // Write data row
    csvFile << Nbp * 1e6 << "," << dose_Gy << ","<< number << "," ;
    csvFile << count_8og << "," << count_tg << "," << count_81 << "," << count_21 << ",";
    csvFile << cluster_ssb_8og << "," << cluster_ssb_tg << "," << cluster_ssb_both << ",";
    csvFile << cluster_2x8og << "," << cluster_multi8og << "," << cluster_2xtg << "," << cluster_multitg << ",";
    csvFile << cluster_8og_tg << "," << cluster_multi_both << ",";
    csvFile << cluster_multi_ssb_8og << "," << cluster_multi_ssb_tg << "," << cluster_multi_ssb_both << ",";
    csvFile << dsb_criterion_removed << "\n";
    
    csvFile.close();
    
    // ============================================================
    // Print summary to console
    // ============================================================
    cout << "\n========================================\n";
    cout << "DNA Damage Analysis Complete\n";
    cout << "========================================\n";
    cout << "Number of base pairs: " << Nbp * 1e6 << "\n";
    cout << "Dose: " << dose_Gy << " Gy\n\n";
    
    cout << "Individual Base Products:\n";
    cout << "  8-oxo-G: " << count_8og << "\n";
    cout << "  TG: " << count_tg << "\n";
    cout << "  Type 81: " << count_81 << "\n";
    cout << "  Type 21: " << count_21 << "\n\n";
    
    cout << "Total Clusters: " << (cluster_ssb_8og + cluster_ssb_tg + cluster_ssb_both + 
                                    cluster_2x8og + cluster_multi8og + cluster_2xtg + 
                                    cluster_multitg + cluster_8og_tg + cluster_multi_both +
                                    cluster_multi_ssb_8og + cluster_multi_ssb_tg + 
                                    cluster_multi_ssb_both) << "\n\n";

    cout << "Sites removed by DSB criterion\n";
    cout << "(TypeClassification == DSB/DSB+/DSB++): " << dsb_criterion_removed << "\n\n";
    
    cout << "Results saved to: " << outputCSV << "\n";
    cout << "========================================\n";
    
    // Cleanup
    delete[] type;
    f->Close();
}
