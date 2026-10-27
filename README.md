# Quantifying High-LET Heavy-Ion-Induced Non-DSB Clustered DNA Lesions Using Geant4-DNA

This repository contains the simulation code and post-processing scripts associated with the manuscript:

> **Quantifying High-LET Heavy-Ion-Induced Non-DSB Clustered DNA Lesions Using Geant4-DNA**  
> Tu Minh Khuong, Sangyong Lim  
> *Scientific Reports* (under revision)

---

## Overview

This work extends the Geant4-DNA molecularDNA framework to enable lesion-specific quantification of oxidative non-DSB clustered DNA damage. The implementation adds an operational direct/quasi-direct base-oxidation scoring channel alongside the existing indirect (·OH-mediated) pathway, and classifies the resulting 8-oxo-7,8-dihydro-2′-deoxyguanosine (8-oxo-dG) and thymidine glycol (dTg) lesions within non-DSB cluster configurations.

Simulations were performed using a reduced-scale human fibroblast chromatin model for carbon, silicon, argon, and iron ions spanning LET values of 13–148 keV/μm.

---

## Repository Structure

```
├── src/                        # Modified Geant4-DNA source files
│   ├── SteppingAction.cc       # Direct/quasi-direct base oxidation scoring
│   │                           # (AccumulateEnergyInBases: 0.8 nm octree search)
│   ├── EventAction.cc          # Per-event energy accumulation and 10 eV threshold
│   ├── AnalysisManager.cc      # Extended damage ntuple with base product codes
│   │                           # (Bp1Product, Bp2Product, BaseName columns)
│   │                           # and direct oxidation flag integration
│   └── DNAHit.cc               # Identify the direct ionized base
├── human-cell-mini.mac         # Geant4 macro file
├── analysis/                   # Post-processing scripts
│   ├── AnalyzeDNADamage.C      # ROOT-based cluster classification (total yields)
└── README.md
```

## Modifications to Geant4-DNA molecularDNA

The following modifications were made to the Geant4-DNA molecularDNA example:

### SteppingAction.cc
- Added `AccumulateEnergyInBases()` method implementing the 0.8 nm octree search for direct/quasi-direct base oxidation scoring
- Energy deposited within 0.8 nm of each base volume is accumulated per base across the full event

### EventAction.cc
- Added `fBaseEnergyMap` to accumulate energy per base (identified by placement, chain, strand, and base pair index)
- Added `fDirectlyOxidizedBases` set populated at end of event using 10 eV threshold
- Added `AccumulateBaseEnergy()` and `GetDirectlyOxidizedBases()` methods

### AnalysisManager.cc
- Added three new columns to the `damage` ntuple:
  - `Bp1Product` — base product code for strand 1 (0: none, 2: dTg indirect, 8: 8-oxo-dG indirect, 21: dTg direct, 81: 8-oxo-dG direct)
  - `Bp2Product` — base product code for strand 2 (same encoding)
  - `BaseName` — integer identifier of the base molecule type
- Integrated direct oxidation flag from `EventAction::GetDirectlyOxidizedBases()` into damage record processing
- Added stochastic product assignment for ·OH-induced base damage (p(8-oxo-dG) = 0.18 for guanine, p(dTg) = 0.90 for thymine)

---

## Dependencies

### Geant4
- **Version:** geant4-11-02-patch-02 [MT] (21 June 2024), multi-threaded mode
- **Download:** https://geant4.org/download/
- **References:**
  - NIM A 506 (2003), 250–303
  - IEEE-TNS 53 (2006), 270–278
  - NIM A 835 (2016), 186–225

### Geant4-DNA molecularDNA example
The simulation is based on the molecularDNA example distributed with Geant4-DNA. 
Incerti, Sébastien, et al. "The geant4-dna project." International Journal of Modeling, Simulation, and Scientific Computing 1.02 (2010): 157-178.

## Scope of this repository

This repository provides a **partial code release** associated with the manuscript.

It contains the source-code modifications and post-processing scripts developed for lesion-specific quantification of direct/quasi-direct and indirect non-DSB clustered base damage. It does **not** contain the complete Geant4-DNA `molecularDNA` example, the full simulation framework, all project-specific configuration files, or a standalone build system.

It is **not** intended, by itself, to provide a fully compilable or independently
runnable version of the complete simulation workflow.

## License

This software is built upon and distributed with modifications to the Geant4 toolkit, which is subject to the **Geant4 Software License, Version 1.0 (28 June 2006)**.

Copyright (c) Copyright Holders of the Geant4 Collaboration, 1994–2006.  
See http://cern.ch/geant4/license for the full license text.

In accordance with the Geant4 Software License:

- Redistributions of this software, in whole or in part, with or without modification, must reproduce the Geant4 copyright notice and license conditions in the software, user documentation, and any other materials provided with the redistributed software.
- User documentation must include the following notice: *"This product includes software developed by Members of the Geant4 Collaboration (http://cern.ch/geant4)."*
- The names "Geant4" and "The Geant4 toolkit" may not be used to endorse or promote software or products derived therefrom without prior written permission.

