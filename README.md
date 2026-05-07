# pABR_click
MATLAB code and data sample for parallel (Poisson) ABR (pABR) click responses, including grand average and ICI-dependent adaptation analysis

## What it does
`draw_pABR_clicks.m` takes preprocessed single-channel ABR data and computes:
- Grand average ABR and cochlear microphonic (CM) across stimulus levels
- ICI-dependent response: how the ABR changes as a function of the 
  inter-click interval (ICI), reflecting auditory adaptation

## Requirements
- MATLAB (tested on R2025)
- Statistics and Machine Learning Toolbox (for `expcdf`)

## How to run
1. Download all files into the same folder
2. Open `draw_pABR_clicks.m` in MATLAB
3. Set `path_read` to the folder containing `B323011003.mat`
4. Run:
```matlab
data_out = draw_pABR_clicks();
```

## Example dataset
`B323011003.mat` contains an example pABR click recording from a single 
session. The file contains:
- `y` — neural recording [n_epochs x n_samples] in uV
- `fs` — sampling rate (Hz)
- `list` — stimulus parameters:
  - `list.level` — target SPL (dB SPL) for each epoch
  - `list.i_onset` — click onset indices within each epoch.
    Sign encodes polarity: positive = condensation, negative = rarefaction.
    Use `abs()` for indexing, `sign()` for polarity.

## Output
`data_out` structure contains:
- `data_out.t` — time axis (ms)
- `data_out.all_levels` — stimulus levels (dB SPL)
- `data_out.y_abr` — grand average ABR [n_levels x n_samples]
- `data_out.y_ici` — ICI-dependent response {n_levels}[n_ici x n_samples]
- `data_out.icis` — ICI values (ms) {n_levels}[1 x n_ici]

## Reference
The analysis methodology and application is described in detail in Chapter 3 of my thesis:
Wang, Y. (2026). *Impacts of Cochlear Synaptopathy on Temporal Dynamics 
of Auditory Neural Responses in the Budgerigar* (Doctoral dissertation, 
University of Rochester). ProQuest Dissertations & Theses.
https://www.proquest.com/docview/3311187178

The stimulus generation follows parallel-ABR developed by Polonenko & Maddox (2019) :
Polonenko, M.J. & Maddox, R.K. (2019). A parallel processing framework 
for auditory brainstem responses. *bioRxiv*.
https://doi.org/10.1101/731760

## Contact
Yingxuan Wang
University of Minnesota
wan04032@umn.edu

## License
MIT License — free to use and modify with attribution.
Copyright (c) 2026 Yingxuan Wang
