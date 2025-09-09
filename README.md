# fNIRS_ICU: Motor-Imagery Analysis
This repository contains MATLAB code to reproduce the motor-imagery (MI) analyses used in our study. It includes a step-by-step SOP, pinned software versions, and QC checklists to allow the pipeline to be run end-to-end.
Repo: [https://github.com/TheOwenLab/fNIRS_ICU](https://github.com/TheOwenLab/fNIRS_ICU)

## System Requirements
- MATLAB: R2020a–R2025a (tested on R2023b and R2024b)  
- OS: Windows 10/11, macOS 12+, Ubuntu 20.04/22.04  
- fNIRS toolbox: version included in this repo under `/fNIRS_toolbox` (commit pinned)  
- SNIRF I/O: included via the toolbox (Shared Near Infrared File format)  

## Quick Overview
1. Clone the repository:  
   `git clone https://github.com/TheOwenLab/fNIRS_ICU.git`  
   `cd fNIRS_ICU`  
2. Start MATLAB in the repo root:  
   `addpath(genpath(pwd));`  
3. Run the MI pipeline on your dataset from `Example_Motor_Imagery_codes/`:  
   - `main_MI_step_1.m` — preprocessing + channel QC  
   - `main_MI_step_2.m` — GLM at subject level  
   - `main_MI_step_3.m` — subject-level summaries/figures  
   - `main_MI_step_4_HbO_HbR_activatedChannels.m` — activation tables  
   - `main_step_5_leave_one_out_HbO_HbR.m` — optional group-level analysis  

## Expected Outputs
Saved under `/outputs/<participant_id>/`:  
- `*preproc.mat`: preprocessed HbO/HbR time-series  
- `*design.mat`: design matrix used for GLM  
- `QC*.png`: channel coverage, motion flags, short-channel status  
- `GLM*.csv`: beta estimates & stats per channel/condition  
- `activated_channels_*.csv`: HbO/HbR activation table  

## Repository Layout

fNIRS_ICU/
Example_Motor_Imagery_codes/
main_MI_step_1.m
main_MI_step_2.m
main_MI_step_3.m
main_MI_step_4_HbO_HbR_activatedChannels.m
main_step_5_leave_one_out_HbO_HbR.m
InspectQualityOfShortChannels.m
Age_healthy_MI.m
GoodSC_Controls_MI.m
IndividualSplineValuesControls_MI.m


## Standard Operating Procedure (SOP)

### 4.1 Environment Setup
In MATLAB (R2020a–R2025a), from the repo root:  
`addpath(genpath(pwd)); % makes /fNIRS_toolbox and code available`  
`ver % record MATLAB and toolbox versions in logs`

### 4.2 Required Inputs
- SNIRF files for each participant (recommended) or compatible formats  
- Events (`*_events.tsv`) with onset, duration, trial_type for task runs  
- `participants.csv` with at least `participant_id`, `age`, `sex`  
- Short-channel mapping via probe definition in SNIRF; verified in QC  

### 4.3 Parameter Block
Default parameters (used in manuscript; modify as needed):  
- `bandpass = [0.01, 0.2];` Hz  
- `motion.tMotion = 0.5;` s  
- `motion.tMask = 1.0;` s  
- `motion.stdThresh = 50;` z-threshold  
- `motion.ampThresh = 0.5;` OD threshold  
- `spline.order = 3;`  
- `hpf_cutoff = 0.01;` Hz  
- `lpf_cutoff = 0.2;` Hz  
- `hrf = 'canonical';` GLM HRF shape  
- `alpha = 0.05;` significance threshold  
- `multicomp = 'FDR';` correction across channels  

### 4.4 Step-by-Step
**Step 1 — Inspecting Short Channels (`InspectQualityOfShortChannels.m`)**  
- Load raw fNIRS data  
- Generate PSD plots and raw time series for SCs  
- Classify SCs as usable if cardiac peak visible and free of artifacts  
- Save list in `GoodSC_Controls_MI.m`  

**Step 2 — Preprocessing & QC (`main_MI_step_1.m`)**  
- Convert SNIRF → Intensity → Optical Density → HbO/HbR  
- Channel pruning, motion detection/correction, SC regression  
- Save preprocessed data and QC artifacts  
QC outputs: `QC_coverage.png`, `QC_motion.png`, `QC_shortchannels.csv`  

**Step 3 — Subject-level GLM (`main_MI_step_2.m`)**  
- Build design matrix from `*_events.tsv`  
- Fit GLM to HbO/HbR  
- Output: `GLM_betas*.csv`, `GLM_stats*.csv`  

**Step 4 — Subject Summaries (`main_MI_step_3.m`)**  
- Beta maps, summary figures, QC tables  

**Step 5 — Activated Channels (`main_MI_step_4_HbO_HbR_activatedChannels.m`)**  
- Apply FDR correction (α=0.05)  
- Output: `activated_channels*.csv` + summary figure  

**Step 6 — Group Summary / Leave-One-Out (`main_step_5_leave_one_out_HbO_HbR.m`)**  
- Aggregate subject betas  
- Compute group means and stability check  

## What to Edit for Your Study
- `Age_healthy_MI.m`: enter participant ages or link `participants.csv`  
- `GoodSC_Controls_MI.m`: rerun to identify usable SCs  
- Update dataset paths in all `main_*.m` scripts  

## Reproducibility Checklist
- MATLAB version + OS recorded in `outputs/log_session.txt`  
- Toolbox commit pinned in `/fNIRS_toolbox`  
- `params_used.mat` saved with all hyperparameters  
- SC quality documented in `GoodSC_Controls_MI.m`  
- QC artifacts saved (coverage, motion, SC lists)  
- All tables exported as CSV with ISO-8601 timestamps  
- Logs written to `outputs/runlog.txt`  

## Citation
If you use this pipeline or derivatives, please cite:  
**TheOwenLab/fNIRS_ICU** (commit `<hash>`), accessed `<date>`.  
