# Self-Driving Polymer Nanoparticle Data & Instrument Control

This repository contains all datasets and source code used in the **self-driving active learning campaign**, as well as non-self-driving experiments described in: **“Self-driving lab for the data-driven design of single-chain polymer nanoparticles.**
---

## Repository Structure

- **Active Learning Experiment**: All data generated during self-driving experiments.
- **Homopolymer / Copolymer Experiment**: Data from earlier experiments.
- **DLS Presets**: Wyatt DYNAMICS-specific parameter files used in automation.
- **templates**: Spreadsheet templates for autonomous instrument control.
- **data**: Contains all experimental data files.

> All experiment files contain a **typo**:  
> ⚠️ _poly(ethylene glycol) methyl ether acrylate (PEGMEA)_ is represented as **“mPEG”**.
> _Note: For homopolymer / copolymer data, lightbox and DLS template spreadsheets were omitted,
as the instrumentation was not configured in a self-driving fashion at that time._

---

## Source Code Overview

The source code provides a **Python wrapper** for the **Wyatt Technology DynaPro III Plate Reader DLS instrument**, integrating with the **DYNAMICS 7** software and **DYNAMICS SDK 3.0.1** (originally written in C#).

### Wrapper Modules

| File | Description |
|------|-------------|
| `wrap_dls_methods.py` | Global tasks: initialization, shutdown, door control, experiment run |
| `wrap_dls_connect.py` | Connects to the API and control server |
| `wrap_dls_functions.py` | Wraps commands from the API library |
| `wrap_dls_events.py` | Handles asynchronous events and instrument state changes |
| `wrap_dls_helpers.py` | Helper methods for simplifying execution of API functions |
| `wrap_dls_dicts.py` | Handles C# enum-type variables |

### Experiment Scripts

| File | Purpose |
|------|---------|
| `sampling.py` | Samples polymer formulations |
| `preproc.py` | Pre-processes data for instruction file generation |
| `postproc.py` | Post-processes polymer characterization data |

> **Figure Generation:**  
> Post-processing includes figure creation for PCA, regression, and SHAP analysis.  
> _Note: Manuscript figures were generated using different code._

---

## Setting Up an Instrument Control File

To begin running a self-driving experiment:

1. **Copy** the `Experiment Tracker Template` (from the `templates` folder) into your data folder.

2. **Set the first 3 parameters** (immediately below the experiment title) as follows:
   - `start cycle` – The initial active learning cycle (e.g., `1` for seed library).
   - `end cycle` – The final cycle (e.g., `4` for seed + 3 active learning generations).
   - `experiment name` – Descriptive name of the experiment.

3. Ensure the **experiment tracker file is stored** in the **same directory** as the `data folder` specified inside the file.

4. The **filename (without extension)** must be **identical** to the `experiment name`.

5. **Do not end the filename** with a number following an underscore.  
   ✅ Acceptable: `exp`, `exp_alpha_1`, `123_exp`, `exp_1_alpha`  
   ❌ Not acceptable: `exp_1`, `exp_alpha_1`

---

## How It Works

The control file is loaded by a parser that:
- Reads parameters and triggers the `active_learning_with_dls()` method.
- Detects updates between generations (especially the `cycle` parameter).
- Enables intelligent pre-processing (e.g., selecting new polymers vs re-sampling).
- ⚠️ Checkpointing functionality within active learning generations is under development.

---

## 📬 Contact

For issues or questions regarding this dataset or software, please contact:

- **Alexander Suponya** – [asuponya3@gatech.edu](mailto:asuponya3@gatech.edu)  
- **Adam Gormley** – [adam.gormley@rutgers.edu](mailto:adam.gormley@rutgers.edu)

---

## 📄 Citation

If you use this repository, please cite the corresponding publication:  
**“Self-driving lab for the data-driven design of single-chain polymer nanoparticles.”**
