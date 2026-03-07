# Flying Base Station Optimization with Genetic Algorithms

This repository optimizes Flying Base Station (FBS) deployment using MATLAB + QuaDRiGa.  
The current workflow supports:
- Single-objective GA (`optimizeBaseStation.m`)
- Optional multi-objective GA / NSGA-II style search (`optimizeBaseStationMoga.m`)
- SINR-based evaluation with optional Macro Base Stations (MBS) and cached MBS power maps

Multiband behavior is intentionally not documented here (per current scope).

## Main Files

- `optimize_base_station_ga.m`: main entry script for experiment setup and GA execution
- `optimizeBaseStation.m`: single-objective GA loop
- `optimizeBaseStationMoga.m`: multi-objective GA loop (users vs power)
- `evaluatePopulation.m`: computes fitness/objectives for each individual
- `SINREvaluation.m`: computes per-user SINR/connectivity/rate metrics
- `precompute_mbs_power_maps.m`: precomputes and caches MBS power maps
- `initializePopulation_uniform.m`, `selectParents.m`, `crossover_blend.m`, `mutate.m`: GA operators

## Optimization Variables (Per FBS)

Each FBS is encoded with 6 decision variables:
1. `x` position (m)
2. `y` position (m)
3. `z` altitude (m)
4. transmit power (W)
5. power on/off status (binary)
6. reserved binary flag (kept for compatibility; multiband details omitted)

So total chromosome size is `6 * numBS`.

## Execution Flow

1. `optimize_base_station_ga.m` defines:
   - Area dimensions and MBS layout
   - Antenna objects
   - Cached MBS power maps (`precompute_mbs_power_maps`)
   - GA parameters (`params`)
2. For single-objective mode (`runMoga = false`):
   - calls `optimizeBaseStation(...)`
3. For multi-objective mode (`runMoga = true`):
   - calls `optimizeBaseStationMoga(...)`

## Fitness and Metrics

`evaluatePopulation.m` gets metrics from `SINREvaluation.m` and supports:
- `targetIdx = 1`: connectivity/power composite score
- `targetIdx = 2`: maximize average connected-user spectral efficiency (`avg_rate_connected_bpsHz`)

Single-objective GA currently uses `targetIdx = 2` (average rate).

Tracked metrics include:
- total connected users
- FBS-connected users
- MBS-connected users
- total transmitted power
- average connected-user rate (bps/Hz)

## Run Instructions

From MATLAB:

```matlab
optimize_base_station_ga
```

## Key Parameters (`optimize_base_station_ga.m`)

- `numBS`: number of FBS nodes
- `populationSize`, `numGenerations`, `initialPopulationSize`
- `crossoverProb`, `mutationProb`, `mutationScale`
- `bounds`: repeated per-FBS variable limits
- `fitnessWeights.beta`, `fitnessWeights.gamma`, `fitnessWeights.fbsWeight`, `fitnessWeights.fbsExponent`
- `maxUsers`, `sinrThreshold`
- `enableLogging`, `enablePerformancePlotting`, `plotTrajectory`
- `spaceLimit`: `[W, H]`
- `mbsCache`: output of `precompute_mbs_power_maps`

## Logging

If `enableLogging = true`, `optimizeBaseStation.m` writes a CSV under `logs/` with run summary fields such as:
- GA weights/probabilities
- number of FBSs
- final FBS coordinates/status/powers
- total power
- connected-user split (FBS vs MBS)
- average rate

## Dependencies

- MATLAB (tested with modern releases)
- QuaDRiGa (required for `power_map` calls)
- Parallel Computing Toolbox recommended (`parfor` in SINR evaluation)

## Notes

- `plotTrajectory` currently supports only `numBS == 1`.
- The repository also contains RL-related scripts/notebooks; they are separate from the GA path described above.
