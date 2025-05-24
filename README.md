# Flying Base Station Optimization using Genetic Algorithms

This project implements a Genetic Algorithm to optimize the **3D location** and **transmission power** of **Flying Base Stations (FBSs)**. The goal is to maximize user connectivity using a realistic SINR evaluation framework based on **QuaDRiGa** propagation models.

---

## Project Structure
your_project_folder/
│
├── optimize_base_station_ga.m          # Main script to run the optimization
│
├── optimizeBaseStation.m               # Core GA loop
├── evaluatePopulation.m                # Fitness evaluation logic (based on SINR)
├── SINREvaluation.m                    # SINR and connectivity computation
│
├── initializePopulation_ppp.m          # Population initialization using Matern Point Process
├── initializePopulation_uniform.m      # Population initialization using uniformly distributed sampling
│
├── selectParents.m                     # Tournament parent selection
├── mutate.m                            # Gaussian + Bernoulli mutation
├── crossover_blend.m                   # Blend crossover operator
│
├── reflectToBounds.m                   # Ensures individuals stay within bounds
├── clampToBounds.m                     # Clamps values between bounds
├── mergeParams.m                       # Merge user-defined and default parameters
├── maternPP.m                          # Generates repulsive FBS layouts
│
└── README.md                           # You are here

---

## Call Order and Execution Flow

1. **`optimize_base_station_ga.m`**
   - Initializes antennas and GA parameters
   - Calls `optimizeBaseStation(...)`

2. **`optimizeBaseStation.m`**
   - Initializes population via `initializePopulation_ppp(...)`
   - Evaluates fitness using `evaluatePopulation(...)`
   - Applies GA operators: `selectParents`, `crossover_blend`, `mutate`
   - Tracks and visualizes best solutions

3. **`evaluatePopulation.m`**
   - Computes fitness based on number of connected users and total power
   - Calls `SINREvaluation(...)`

4. **`SINREvaluation.m`**
   - Computes SINR using QuaDRiGa `power_map` functions
   - Returns number of connected users and total transmitted power

---

## Dependencies

- MATLAB R2021b or later
- [QuaDRiGa](https://quadriga-channel-model.de) channel simulator (installed and licensed)
- Required Toolboxes:
  - Optimization Toolbox
  - Statistics and Machine Learning Toolbox

---
## Important Notes

Before running the main script, **please ensure** that:

- All project folders (containing the `.m` files)
- The QuaDRiGa source directory

are added to your MATLAB path. You can do this from the MATLAB command window:

```matlab
addpath(genpath(pwd));                % Adds all project subfolders
addpath(genpath('path/to/quadriga')); % Replace with the actual path to your QuaDRiGa source


## How to Run

Run the optimization from the MATLAB command window:

```matlab
>> optimize_base_station_ga
```
To modify the optimization behavior, edit the params struct inside optimize_base_station_ga.m:

```
params = struct(...
    'populationSize', 10, ...
    'numGenerations', 10, ...
    'crossoverProb', 0.3, ...
    'mutationProb', 0.7, ...
    'numBS', numBS, ...
    'bounds', repmat([0 1500; 0 1500; 20 150; 7 10.5; 0 1], numBS, 1), ...
    'mutationScale', 0.17, ...
    'spaceLimit', 1500, ... 
    'verbose', 1 ...
);
```

