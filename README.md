fNIRS_ICU: Motor-Imagery Analysis
This repository contains MATLAB code to reproduce the motor-imagery (MI) analyses used in our study. It includes a step-by-step SOP, pinned software versions, and QC checklists to allow the pipeline to be run end-to-end.

Repo: https://github.com/TheOwenLab/fNIRS_ICU

1.	System requirements
•	MATLAB: R2020a–R2025a (tested on R2023b and R2024b). 
•	OS: Windows 10/11, macOS 12+, Ubuntu 20.04/22.04
•	fNIRS toolbox: version included in this repo under /fNIRS_toolbox (commit pinned)
•	SNIRF I/O: included via the toolbox (Shared Near Infrared File format)

2.	Quick Overview
•	Download or clone the repository: git clone https://github.com/TheOwenLab/fNIRS_ICU.git
•	Enter the repo: cd fNIRS_ICU
•	Start MATLAB in the repo root and run: addpath(genpath(pwd));
•	Run the MI pipeline on your dataset from Example_Motor_Imagery_codes:
o	main_MI_step_1 (preprocessing + channel QC)
o	main_MI_step_2 (GLM at subject level)
o	main_MI_step_3 (subject-level summaries/figures)
o	main_MI_step_4_HbO_HbR_activatedChannels
o	main_step_5_leave_one_out_HbO_HbR (optional group-level analysis)

Expected outputs (saved under /outputs/<participant_id>/):
•	*_preproc.mat: preprocessed HbO/HbR time-series
•	*_design.mat: design matrix used for GLM
•	QC_*png: channel coverage, motion flags, short-channel status
•	GLM_*csv: beta estimates & stats per channel/condition
•	activated_channels_*csv: HbO/HbR activation table

3.	Repository layout
fNIRS_ICU/
•	Example_Motor_Imagery_codes/
o	main_MI_step_1.m
o	main_MI_step_2.m
o	main_MI_step_3.m
o	main_MI_step_4_HbO_HbR_activatedChannels.m
o	main_step_5_leave_one_out_HbO_HbR.m
o	InspectQualityOfShortChannels.m
o	Age_healthy_MI.m
o	GoodSC_Controls_MI.m
o	IndividualSplineValuesControls_MI.m
4.	Standard Operating Procedure (SOP)
4.1 Environment setup
In MATLAB (R2020a–R2025a), from the repo root:
addpath(genpath(pwd)); % makes /fNIRS_toolbox and code available
ver % record MATLAB and toolbox versions in logs

4.2 Required inputs
•	SNIRF files for each participant (recommended) or compatible formats supported by /fNIRS_toolbox
•	Events (*_events.tsv) with onset, duration, trial_type for task runs
•	Participants table (participants.csv) with at least participant_id, age, sex
•	Short-channel mapping: provided via probe definition in SNIRF; verified by QC step

4.3 Parameter block (edit in code or pass as struct)
Default parameters (mirroring those used in the manuscript; modify as needed):
bandpass = [0.01, 0.2]; Hz; task band for MI
motion.tMotion = 0.5; s; event size for motion detection
motion.tMask = 1.0; s; padding around motion
motion.stdThresh = 50; z-threshold
motion.ampThresh = 0.5; OD threshold
spline.order = 3; spline interpolation order
hpf_cutoff = 0.01; high-pass (Hz)
lpf_cutoff = 0.2; low-pass (Hz)
hrf = ‘canonical’; GLM HRF shape
alpha = 0.05; significance threshold
multicomp = ‘FDR’; FDR correction across channels

4.4 Step-by-step
Step 1 — Inspecting Short channels (InspectQualityOfShortChannels.m)
•	Load in raw fNIRS data
•	For each participant, predefined short-separation (SC) channels are visualized individually
•	For each SC, the script generates:
• Power spectral density (PSD) plots using Welch’s method to assess frequency-domain characteristics.
• Raw time series plots to evaluate baseline stability and noise artifacts.
•	SCs are classified as “good” if the PSD shows a distinct physiological peak around the cardiac frequency (~1 Hz) and the time series is free of saturation or excessive high-frequency noise.
•	Document the set of usable SCs for each participant in GoodSC_Controls_MI.m, which is then referenced during preprocessing and SC regression.

Step 2 — Preprocessing & QC (main_MI_step_1.m)
•	Load SNIRF → Intensity → Optical Density → HbO/HbR (modified Beer–Lambert).
•	Channel pruning: drop channels with SNR below threshold or missing short-channel pairing.
•	Motion detection & correction: detect spikes (stdThresh/ampThresh), apply spline interpolation, temporal filtering.
•	Short-channel regression: regress nuisance (SC) from long channels if SC is available and passes QC.
•	Save preprocessed HbO/HbR and QC artifacts.

QC outputs to check:
•	QC_coverage.png: expected coverage consistent with probe layout.
•	QC_motion.png: spikes detected and corrected segments.
•	QC_shortchannels.csv: list of “good” short channels used.

Step 3 — Subject-level GLM (main_MI_step_2.m)
•	Builds design matrix from *_events.tsv.
•	Fits GLM to HbO and HbR using selected HRF.
•	Outputs GLM_betas_[HbO|HbR].csv and GLM_stats_[HbO|HbR].csv.

Step 4 — Subject-level summaries (main_MI_step_3.m)
•	Produces per-subject figures/tables: beta maps per condition (HbO/HbR), QC summary table (coverage, motion %, SC regression used).

Step 5 — Activated channels (main_MI_step_4_HbO_HbR_activatedChananels.m)
•	Applies FDR-corrected thresholding (alpha=0.05) to GLM results.
•	Exports activated_channels_[HbO|HbR].csv and summary figure.

Step 6 — Group summary / leave-one-out (main_step_5_leave_one_out_HbO_HbR.m)
•	Aggregates subject betas; computes group means and a leave-one-out stability check.
•	Outputs group-level CSVs and figures.

What to edit in your own study
•	Age_healthy_MI.m: enter participant ages (or point to your participants.csv).
•	GoodSC_Controls_MI.m: run to identify usable short channels; the script saves a list consumed by Step 1.
•	Replace placeholder data paths in the scripts with your dataset (keep folder/file conventions) and update paths at the top of each main_*.m script.

Reproducibility checklist
•	MATLAB version and OS recorded in outputs/log_session.txt (auto-written at runtime).
•	Toolbox commit pinned (this repo’s /fNIRS_toolbox is used; no external updates).
•	params_used.mat saved with all processing hyperparameters.
•	Assess SC quality and note good channels in GoodSC
•	QC artifacts (coverage, motion, SC lists) saved.
•	All tables exported as CSV with ISO-8601 timestamps.
•	Script stdout captured to outputs/runlog.txt.


Citation
If you use this pipeline or derivatives, please cite the paper and this repository:
<Paper citation once available> 
TheOwenLab/fNIRS_ICU (commit <HASH>), accessed <YYYY-MM-DD>.

![Uploading image.png…]()
