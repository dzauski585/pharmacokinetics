# Anesthetic PK/PD Simulation

A Python-based pharmacokinetic/pharmacodynamic (PK/PD) simulation tool for modeling propofol, remifentanil, and dexmetetomidine concentrations and their effects during anesthesia.

## Pharacokinetics Review

### PK Software

There are commercially available target controlled infusion systems. They traditionally are for soley propofol or a combination of propofol and remifentantil. These ssytems are used in most countries with one notable exception, the United States. The models available for use are the Marsh and Schneider models for only propofol and the Minto model for a combination of propofol and remifentanil. The simulation software described in this readme will have addiitional medications and models available. 

Depending on the model used there are different goal targets for the medication. The Marsh model uses plasma concentration while the others use effect site concentration. The Ce requires an additional model component; the effect site and therefore a keo value desciribn the movement of the drug to this compartment. Cp and Ce are both dependent on hemodynamic factors which are largely assumed to be normal for the purposes of simulation. It should also be noted that Ce is difficult to impossible to measure. It should also be noted that to allow for increases in Ce the Cp must be overshot to ensure that drug delivery continues into the fourth compartment. 

### Three Compartment Model

The classic model is made up of three compartments V1,2,3. V1 is the large central vascular compartment where are bolus and infusions start. V2 is the vessel rich areas like muscles, and V3 is the vessel poor areas such as fat. To model this effect a tri-exponential curve is uses to show the decline in Cp after bolus injection. The Cp at anytime may be derived from the following equation and graph: C(t) = Ae^-alphaT + Be^betaT + Ce^-gammaT. A, B, C represent phase coefficients which sum to the Cp after a bolus. Alpha, Beta, and Gamma are phase rate constants. e is the natural logarithm. T is time and C is concentration. 

![alt text](image-1.png)
![alt text](image.png)

### Marsh Model vs Schneider Model

The Marsh model was first and did not include elderly or obese patients in their cohort. The model was later updated with a larger keo value to increase accuracy. The schnider model incoporates lean body mass and age. The LBM formula used is the James formula. Unfortunatey the model allows for exponential increases in k10 with total body weight and it decreases with lbm. These are opposing actions so the companies using this model impose a lbm limit to utilize the maximum k10. Both models do not accuratlye represent the increase obesity epidemic around the world.

![alt text](image-2.png)

### Eleveld Model

This model is meant to be a unified model that can accurately deliver propofol or remifentanil to a patient of any age or weight. This model is not in use by any commercial TCI pump but is in use in many research applications. It was born out of the open TCI initiative and uses data from many subjects. The model is far more complex.

## Overview

This software implements the multiple models for anesthetic agents to simulate:

- Plasma concentrations (Cp)
- Effect site concentrations (Ce)
- Bispectral Index (BIS) values
- Wake-up times following different dosing regimens

The simulation allows for customized patient demographics (weight, height, age, sex) and flexible dosing regimens including boluses and variable-rate infusions.

### Models

#### Propofol

- Schneider
- Marsh
- Eleveld
- Roberts (10-8-6)

#### Remifentanil

- Minto

#### Dexmetetomdine


## Features

- Patient-specific PK/PD parameter calculation based on weight, height, age, and sex
- Support for complex dosing regimens including boluses and tiered infusions
- Calculation of key clinical metrics:
  - Time to surgical anesthesia (BIS < 60)
  - Time to therapeutic effect site concentration
  - Predicted wake-up time
  - Total drug consumption
- Detailed visualization of simulation results

## Dependencies

- NumPy
- Matplotlib
- SciPy

## Installation

```bash
pip install numpy matplotlib scipy
```

## Usage

### Basic Usage

```python
from propofol_simulation import main

# Run with default settings
main()
```

### Customizing Parameters

Edit the parameters in the `main()` function:

```python
# Patient parameters
weight = 70      # kg
height = 170     # cm
age = 40         # years
sex = 'male'     # 'male' or 'female'

# Bolus and infusion parameters
dosing_regimen = [
    {'type': 'bolus', 'time': 0, 'dose': 2.0, 'duration': 0.1},  # mg/kg bolus at t=0
    {'type': 'infusion', 'time': 0.1, 'rate': 166.67, 'duration': 10},  # mcg/kg/min for 10 min
    {'type': 'infusion', 'time': 10.1, 'rate': 133.33, 'duration': 10},  # mcg/kg/min for 10 min
    {'type': 'infusion', 'time': 20.1, 'rate': 100, 'duration': 99.9},  # mcg/kg/min until 120 min
]

# Simulation parameters
total_time = 180  # minutes
time_step = 0.1   # minutes
```

## Functions

### `get_schneider_parameters(weight, height, age, sex)`

Calculates PK/PD parameters based on the Schneider model using patient demographics.

**Parameters:**

- `weight`: Patient weight in kg
- `height`: Patient height in cm
- `age`: Patient age in years
- `sex`: Patient sex ('male' or 'female')

**Returns:**

- Dictionary containing all PK/PD parameters

### `run_simulation(model_params, weight, dosing_regimen, total_time, dt)`

Runs the PK/PD simulation using the provided parameters.

**Parameters:**

- `model_params`: Dictionary of PK/PD parameters from `get_schneider_parameters()`
- `weight`: Patient weight in kg
- `dosing_regimen`: List of dosing events (boluses and infusions)
- `total_time`: Total simulation time in minutes
- `dt`: Time step for simulation in minutes

**Returns:**

- Tuple containing: time points, infusion rates, bolus times, plasma concentrations, effect site concentrations, BIS values, and wake-up time

### `plot_results(t, infusion_rates, bolus_times, cp, ce, bis, wake_time)`

Creates a visualization of the simulation results.

**Parameters:**

- Results from `run_simulation()`

**Returns:**

- Matplotlib figure object

## Dosing Regimen Format

The dosing regimen is specified as a list of dictionaries:

```python
[
    {
        'type': 'bolus',     # 'bolus' or 'infusion'
        'time': 0,           # Start time in minutes
        'dose': 2.0,         # Dose in mg/kg (for bolus)
        'duration': 0.1      # Duration in minutes
    },
    {
        'type': 'infusion',  # 'bolus' or 'infusion'
        'time': 0.1,         # Start time in minutes
        'rate': 166.67,      # Rate in mcg/kg/min (for infusion)
        'duration': 10       # Duration in minutes
    }
]
```

## Example Output

The simulation produces a three-panel plot showing:

1. Infusion rates and bolus times
2. Plasma and effect site concentrations
3. BIS values over time

A summary of key clinical metrics is also printed to the console:

- Time to surgical anesthesia (BIS < 60)
- Time to therapeutic concentration
- Wake-up time prediction
- Total drug consumption

## Clinical Applications

This simulation can be used for:

- Educational purposes in anesthesiology
- Pre-procedure planning for complex cases
- Research on optimizing propofol dosing regimens
- Better understanding of PK/PD relationships

## Limitations

- This is a simulation tool and should not be used as the sole basis for clinical decisions
- The model is based on population averages and may not perfectly predict individual patient responses
- Always consult with a qualified anesthesiologist and/or certified registered nurse anesthetist for clinical applications

## License

[Include license information here]

## Citations

Al-Rifai, Ziad & Mulvey, David. (2015). Principles of total intravenous anaesthesia: basic pharmacokinetics and model descriptions. BJA Education. 16. 10.1093/bjaceaccp/mkv021.
