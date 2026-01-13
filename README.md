
# Zebiak and Cane (ZC) Model – Execution Guide

This repository contains a simplified implementation of the Zebiak and Cane (1987) coupled ocean–atmosphere model, together with preprocessing scripts and namelists required to run standalone ocean, atmosphere, and fully coupled simulations.

This document summarizes the directory structure and the standard workflow to compile and run the model.

---

## Directory Structure
```
CODES/
Core source codes and Makefile

Forcing/
Scripts and files for external forcing and background fields
├── OCN/ Ocean mean fields and preprocessing scripts
├── WIND/ Wind forcing files and preprocessing scripts
├── Tzm/ Background vertical temperature gradient fields
└── SST/ Background SST fields

INPUT/
Grid definition files and generation scripts
├── OCN/ Ocean grid
└── ATM/ Atmospheric grid

RUN/
Namelists and helper scripts for model execution
├── OGCM/ Standalone ocean model runs
└── CGCM/ Coupled ocean–atmosphere model runs
```
---
## Workflow Overview
### 0. Clone the Repository
First, clone the Zebiak–Cane model repository to your local machine:

```bash
git clone https://github.com/shokido/ZC_model
cd ZC_model
```

### 1. Compile the Source Codes
Move to the CODES directory and compile the model.
``` bash
cd CODES
```
Edit the Makefile if necessary, paying particular attention to the netCDF include and library paths.
``` bash
make -f Makefile
```
After successful compilation, the following executables should be created:

- exec_solver_ocn_dyn.out
  - Runs the ocean dynamical model only, forced by prescribed winds.
- exec_solver_ocn_full.out
  - Runs the ocean dynamical model including the SST equation, forced by winds. This is typically used for hindcast experiments to evaluate ocean model performance.
- exec_solver_atm.out
  - Computes atmospheric response using a Gill-type model, forced by prescribed SST.
- exec_solver_couple_full.out
  - Runs the fully coupled Zebiak–Cane model after specifying background fields.
