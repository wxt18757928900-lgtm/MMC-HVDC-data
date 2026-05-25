# MMC-HVDC Wideband Field Data and Fitting Scripts

This repository contains the processed field-measured frequency-response data and the MATLAB fitting scripts used for the wideband equipment models in the associated MMC-HVDC converter-station paper.

## Contents

- `data/arm_reactor_impedance.csv`: measured arm-reactor terminal impedance.
- `data/transformer_dm_impedance.csv`: measured converter-transformer differential-mode impedance.
- `data/transformer_cm_impedance.csv`: measured converter-transformer common-mode impedance.
- `matlab/fit_arm_reactor.m`: admittance-domain vector fitting and branch extraction for the arm reactor.
- `matlab/fit_transformer_dm.m`: Foster-I differential-mode fitting for the converter transformer.
- `matlab/fit_transformer_cm.m`: Foster-II common-mode fitting for the converter transformer.
- `matlab/run_all.m`: runs the three fitting scripts.

Each CSV file has three columns:

```text
frequency_hz, phase_deg, magnitude_ohm
```

The data are processed frequency-response records extracted from field measurements of installed primary equipment. Station-identifying information and unrelated PSCAD project files are not included.

## MATLAB Requirements

The scripts were prepared for MATLAB. Depending on the script, the following toolboxes may be required:

- RF Toolbox, for `rationalfit` and `freqresp` in `fit_arm_reactor.m`;
- Optimization Toolbox, for `lsqnonlin` in the transformer fitting scripts;
- Signal Processing Toolbox, for `findpeaks`.

## How to Run

From the repository root, run:

```matlab
run('matlab/run_all.m')
```

or run each script separately:

```matlab
run('matlab/fit_arm_reactor.m')
run('matlab/fit_transformer_dm.m')
run('matlab/fit_transformer_cm.m')
```

The fitted branch parameters and verification figures are written to `results/`.

## Citation

If you use these data or scripts, please cite the associated paper:

```text
Propagation Characteristics of High-Frequency Signals in an MMC-HVDC Converter Station Based on Field-Measured Wideband Modeling
```

The final DOI and bibliographic information should be added after publication.

