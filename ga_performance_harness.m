% GA_PERFORMANCE_HARNESS  Structured test harness for GA performance assessment
%
% This script runs the GA across multiple configurations with repeated
% trials per configuration, collects convergence histories and final
% metrics, and saves everything to disk for post-hoc analysis.
%
% ====================================================================
% USAGE:
%   1. Adjust the "User Configuration" section below.
%   2. Run:  ga_performance_harness
%   3. Results are saved to ./ga_perf_results/<timestamp>/
%   4. Use ga_postanalysis.m on the saved results for visualization.
% ====================================================================
%
% OUTPUT STRUCTURE (saved as results.mat):
%   allRuns     - struct array, one entry per (config x trial)
%   configs     - struct array describing each configuration
%   metadata    - timestamp, MATLAB version, total wall time
%
% NOTE: This file does NOT modify any existing project files.

clearvars; close all; clc;

%% ==================== User Configuration ====================

% --- Number of independent trials per configuration (for statistics) ---
nTrials = 3;

% --- Mode: 'quick', 'medium', 'full' ---
%   quick  : small pop, few gens  (debugging / sanity check)
%   medium : moderate settings     (comparing a handful of configs)
%   full   : large pop, many gens  (overnight sweep)
runMode = 'medium';

switch runMode
    case 'quick'
        defaultPop   = 10;  defaultGens = 5;   defaultInitPop = 15;
    case 'medium'
        defaultPop   = 15;  defaultGens = 15;  defaultInitPop = 25;
    case 'full'
        defaultPop   = 20;  defaultGens = 30;  defaultInitPop = 40;
end

% --- Scenario constants (match your main script) ---
coverageFreq = 2e9;
capacityFreq = 2.6e9;
fbsAntenna = [setup_antenna(coverageFreq), setup_antenna(capacityFreq)];
mbsAntenna = setup_antenna();

numMbs = 1;
W = 2000; H = 1500;
margin = 100;
ISD = 500;
[xs, ys] = generate_hex_sites(W, H, ISD, margin, numMbs);
mbs_height = 25;
mbs_power  = 20;
[mbs_params, antennaObjectMbs, containsMbs, ~] = ...
    pack_mbs_params(xs, ys, mbs_height, mbs_power, mbsAntenna);
tempForX = mbs_params(1,:);
mbs_params(1,:) = mbs_params(2,:);
mbs_params(2,:) = tempForX;

subset = struct('xmin', 0, 'xmax', W, 'ymin', 0, 'ymax', H);
cache = precompute_mbs_power_maps( ...
    mbs_params, antennaObjectMbs, subset, ...
    '3GPP_38.901_UMa_LOS', 'quick', 1.5, './cache_mbs_maps');

%% ==================== Configuration Grid ====================
% Each "config" is a named struct describing one hyperparameter setting.
% The sweep below generates the full factorial grid.

% --- Hyperparameters to sweep ---
populationSizes  = [20, 30];
crossoverProbs   = [0.3,0.5];
mutationProbs    = [0.3];
mutationScales   = [0.15];
numBSOptions     = [1, 2, 3];

% --- Fitness weight configs (named presets) ---
fitnessPresets = {
    struct('name', 'balanced',      'beta', 0.8, 'gamma', 0.0, 'fbsWeight', 0.4, 'fbsExponent', 1);
    struct('name', 'users_only',    'beta', 1.0, 'gamma', 0.0, 'fbsWeight', 0.0, 'fbsExponent', 1);
    struct('name', 'power_aware',   'beta', 0.7, 'gamma', 0.0, 'fbsWeight', 0.2, 'fbsExponent', 0.5);
    struct('name', 'fbs_biased',    'beta', 0.5, 'gamma', 0.0, 'fbsWeight', 0.2, 'fbsExponent', 0.5);
};

% --- Build the full configuration list ---
configs = [];
cfgIdx = 0;

for nBS = numBSOptions
    for pi = 1:numel(populationSizes)
        for ci = 1:numel(crossoverProbs)
            for mi = 1:numel(mutationProbs)
                for si = 1:numel(mutationScales)
                    for fi = 1:numel(fitnessPresets)
                        cfgIdx = cfgIdx + 1;
                        fp = fitnessPresets{fi};

                        cfg = struct();
                        cfg.id               = cfgIdx;
                        cfg.numBS            = nBS;
                        cfg.populationSize   = populationSizes(pi);
                        cfg.crossoverProb    = crossoverProbs(ci);
                        cfg.mutationProb     = mutationProbs(mi);
                        cfg.mutationScale    = mutationScales(si);
                        cfg.fitnessPreset    = fp.name;
                        cfg.beta             = fp.beta;
                        cfg.gamma            = fp.gamma;
                        cfg.fbsWeight        = fp.fbsWeight;
                        cfg.fbsExponent      = fp.fbsExponent;
                        cfg.numGenerations   = defaultGens;

                        if isempty(configs)
                            configs = cfg;
                        else
                            configs(end+1) = cfg; %#ok<SAGROW>
                        end
                    end
                end
            end
        end
    end
end

totalConfigs = numel(configs);
totalRuns    = totalConfigs * nTrials;
fprintf('=== GA Performance Harness ===\n');
fprintf('Mode: %s | Configs: %d | Trials/config: %d | Total runs: %d\n', ...
    runMode, totalConfigs, nTrials, totalRuns);

%% ==================== Output Directory ====================
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
outDir = fullfile(pwd, 'ga_perf_results', timestamp);
if ~exist(outDir, 'dir'); mkdir(outDir); end
fprintf('Results will be saved to: %s\n\n', outDir);

%% ==================== Main Execution Loop ====================
allRuns = [];
runCounter = 0;
wallStart = tic;

for c = 1:totalConfigs
    cfg = configs(c);
    fprintf('\n--- Config %d/%d: nBS=%d pop=%d xo=%.2f mut=%.2f mutS=%.2f fitness=%s ---\n', ...
        c, totalConfigs, cfg.numBS, cfg.populationSize, ...
        cfg.crossoverProb, cfg.mutationProb, cfg.mutationScale, cfg.fitnessPreset);

    for t = 1:nTrials
        runCounter = runCounter + 1;

        % Set reproducible seed
        seed = c * 1000 + t;
        rng(seed);

        % Build params struct for this run
        numBS = cfg.numBS;
        params = struct( ...
            'enablePerformancePlotting', false, ...
            'enableLogging',            false, ...
            'plotTrajectory',           false, ...
            'initialPopulationSize',    max(defaultInitPop, cfg.populationSize), ...
            'populationSize',           cfg.populationSize, ...
            'numGenerations',           cfg.numGenerations, ...
            'crossoverProb',            cfg.crossoverProb, ...
            'mutationProb',             cfg.mutationProb, ...
            'mutationScale',            cfg.mutationScale, ...
            'fitnessWeights',           struct('beta', cfg.beta, 'gamma', cfg.gamma, ...
                                               'fbsWeight', cfg.fbsWeight, 'fbsExponent', cfg.fbsExponent), ...
            'maxUsers',                 1000, ...
            'sinrThreshold',            5, ...
            'logFile',                  '', ...
            'numBS',                    numBS, ...
            'mbsBandId',                0, ...
            'bounds',                   repmat([0 W; 0 H; 20 150; 7 10.5; 0 1; 0 1], numBS, 1), ...
            'spaceLimit',               [W, H], ...
            'mbsCache',                 cache, ...
            'verbose',                  0, ...
            'targetIdx',                1 ...  % 1 -> connectivity, 2 -> avg sum rate
        );

        % Run the GA
        runTimer = tic;
        try
            [bestInd, bestFit, history] = optimizeBaseStation( ...
                fbsAntenna, containsMbs, antennaObjectMbs, mbs_params, params);
            status = 'ok';
            errmsg = '';
        catch ME
            bestInd = NaN(1, 6*numBS);
            bestFit = NaN;
            history = struct();
            status  = 'error';
            errmsg  = ME.message;
            fprintf(2, '  ERROR in trial %d: %s\n', t, errmsg);
        end
        elapsed = toc(runTimer);

        % Store result
        run = struct();
        run.configId       = c;
        run.trial          = t;
        run.seed           = seed;
        run.status         = status;
        run.errmsg         = errmsg;
        run.elapsed_sec    = elapsed;
        run.bestFitness    = bestFit;
        run.bestIndividual = bestInd;

        % Extract convergence trace
        if isfield(history, 'bestFitness')
            run.bestFitnessTrace = history.bestFitness(:)';
            run.avgFitnessTrace  = history.avgFitness(:)';
            run.stdFitnessTrace  = history.stdFitness(:)';
        else
            run.bestFitnessTrace = NaN;
            run.avgFitnessTrace  = NaN;
            run.stdFitnessTrace  = NaN;
        end

        % Extract connected users trace
        if isfield(history, 'bestConnectedUsers')
            run.bestUsersTrace = history.bestConnectedUsers(:)';
            run.avgUsersTrace  = history.avgConnectedUsers(:)';
        else
            run.bestUsersTrace = NaN;
            run.avgUsersTrace  = NaN;
        end

        % Extract operator stats
        if isfield(history, 'crossovers')
            run.crossoversTrace = history.crossovers(:)';
            run.mutationsTrace  = history.mutations(:)';
        else
            run.crossoversTrace = NaN;
            run.mutationsTrace  = NaN;
        end

        % Variable trajectories
        if isfield(history, 'bestIndividuals')
            run.bestIndividualsTrace = history.bestIndividuals;
        end
        if isfield(history, 'avgVariables')
            run.avgVariablesTrace = history.avgVariables;
        end

        % Generation timing
        if isfield(history, 'time') && isfield(history.time, 'generation')
            run.genTimesTrace = history.time.generation(:)';
        end

        % Raw physical metrics (scale-independent, safe for cross-config comparison)
        if isfield(history, 'rawMetrics')
            run.numUsers         = history.rawMetrics.numUsers;
            run.fbsUsers         = history.rawMetrics.fbsUsers;
            run.mbsUsers         = history.rawMetrics.mbsUsers;
            run.transmittedPower = history.rawMetrics.transmittedPower;
            run.avgRate          = history.rawMetrics.avgRate;
            run.fbsFreqFlags     = history.rawMetrics.fbsFreqFlags;
        else
            run.numUsers         = NaN;
            run.fbsUsers         = NaN;
            run.mbsUsers         = NaN;
            run.transmittedPower = NaN;
            run.avgRate          = NaN;
            run.fbsFreqFlags     = NaN;
        end

        if isempty(allRuns)
            allRuns = run;
        else
            allRuns(end+1) = run; %#ok<SAGROW>
        end

        fprintf('  Trial %d/%d: fitness=%.4f | users=%d (fbs=%d mbs=%d) | power=%.2fW | rate=%.2f | time=%.1fs | %s\n', ...
            t, nTrials, bestFit, ...
            run.numUsers, run.fbsUsers, run.mbsUsers, ...
            run.transmittedPower, run.avgRate, elapsed, status);
    end

    % --- Periodic checkpoint save (every 10 configs) ---
    if mod(c, 10) == 0
        save(fullfile(outDir, 'checkpoint.mat'), 'allRuns', 'configs', 'nTrials', 'runMode', '-v7.3');
        fprintf('  [Checkpoint saved at config %d/%d]\n', c, totalConfigs);
    end
end

totalWallTime = toc(wallStart);

%% ==================== Save Final Results ====================
metadata = struct();
metadata.timestamp     = timestamp;
metadata.matlabVersion = version;
metadata.runMode       = runMode;
metadata.nTrials       = nTrials;
metadata.totalConfigs  = totalConfigs;
metadata.totalRuns     = totalRuns;
metadata.totalWallTime = totalWallTime;
metadata.scenarioW     = W;
metadata.scenarioH     = H;

save(fullfile(outDir, 'results.mat'), 'allRuns', 'configs', 'metadata', 'nTrials', '-v7.3');

% --- Also export a summary CSV for quick inspection ---
summaryTable = table();
for r = 1:numel(allRuns)
    row = struct();
    row.configId      = allRuns(r).configId;
    row.trial         = allRuns(r).trial;
    row.seed          = allRuns(r).seed;
    row.numBS         = configs(allRuns(r).configId).numBS;
    row.populationSize= configs(allRuns(r).configId).populationSize;
    row.crossoverProb = configs(allRuns(r).configId).crossoverProb;
    row.mutationProb  = configs(allRuns(r).configId).mutationProb;
    row.mutationScale = configs(allRuns(r).configId).mutationScale;
    row.fitnessPreset = string(configs(allRuns(r).configId).fitnessPreset);
    row.beta          = configs(allRuns(r).configId).beta;
    row.gamma         = configs(allRuns(r).configId).gamma;
    row.fbsWeight     = configs(allRuns(r).configId).fbsWeight;
    row.bestFitness   = allRuns(r).bestFitness;
    row.numUsers      = allRuns(r).numUsers;
    row.fbsUsers      = allRuns(r).fbsUsers;
    row.mbsUsers      = allRuns(r).mbsUsers;
    row.txPower       = allRuns(r).transmittedPower;
    row.avgRate       = allRuns(r).avgRate;
    row.elapsed_sec   = allRuns(r).elapsed_sec;
    row.status        = string(allRuns(r).status);

    summaryTable = [summaryTable; struct2table(row)]; %#ok<AGROW>
end
writetable(summaryTable, fullfile(outDir, 'summary.csv'));

fprintf('\n=== Harness Complete ===\n');
fprintf('Total wall time: %.1f minutes (%.1f hours)\n', totalWallTime/60, totalWallTime/3600);
fprintf('Results saved to: %s\n', outDir);
fprintf('  - results.mat  (full data)\n');
fprintf('  - summary.csv  (flat table)\n');
fprintf('\nNext step: run ga_postanalysis(''%s'') for visualization.\n', outDir);
