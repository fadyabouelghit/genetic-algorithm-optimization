%GA_ABLATION_CAPACITY  Paired ablation: does GA control of the MBS capacity
% slot help connectivity? For every (code, seed) the GA runs TWICE — once
% with gaControlsMbsCapacity=false, once with true — under randomizeGA=true
% with gaSeed=seed, so both arms share the exact same GA randomness (init
% population + operators) and differ ONLY in the capacity-slot control.
% original: 2 codes x 20 seeds x 2 flags = 80 runs;
% connonly_4scen: 4 codes x 20 seeds x 2 flags = 160 runs (see variant below).
%
% Results land in ga_studies/<studyName>_<DATE>_<TIME>/runs.csv;
% each row carries code, gaSeed and gaControlsMbsCapacity (plus the per-tier
% connected-user split), so paired (seed, flag) comparisons can be queried
% directly from the CSV.
%
% Dry run (print the plan without training):
%   dryRun = true; ga_ablation_capacity

if ~exist('dryRun', 'var'); dryRun = false; end

% ---- Study variant ----------------------------------------------------
%   'original'       : as-published ablation, default fitness shaping
%                      (fbsWeight 0.4, gamma from Z code, beta 0.8); 2 codes.
%   'connonly_4scen' : pure connectivity (fbsWeight 0, gamma 0, beta 1.0)
%                      over the macro x FBS scenario grid
%                      (1-1-1, 2-1-1, 1-2-1, 2-2-1 = {1,2} macro x {1,2} FBS).
variant = 'connonly_4scen';     % 'original' | 'connonly_4scen'
useZ2Variants = false;          % (original only) true -> Z=2 codes (gamma 0.1)

seeds = 1:20;
switch variant
    case 'original'
        if useZ2Variants
            codes = {'1-1-2', '2-1-2'};
        else
            codes = {'1-1-1', '2-1-1'};
        end
        studyName = 'capacity_ablation_v1';
        pureConn  = false;
    case 'connonly_4scen'
        codes     = {'1-1-1', '2-1-1', '1-2-1', '2-2-1'};
        studyName = 'capacity_ablation_connonly_4scen_v1';
        pureConn  = true;
    otherwise
        error('ga_ablation_capacity:badVariant', 'Unknown variant "%s".', variant);
end
capacityFlags = [false, true];

% Build aligned (code, overrides) lists for ga_study 'list' mode: row i of
% the overrides struct array pairs with runCodes{i}.
N = numel(codes) * numel(seeds) * numel(capacityFlags);
runCodes  = cell(1, N);
overrides = repmat(struct( ...
    'randomizeGA',           true, ...
    'gaSeed',                NaN, ...
    'gaControlsMbsCapacity', false, ...
    'gaControlsFbsBand',     false, ...
    'populationSize',        25, ...
    'numGenerations',        25, ...
    'targetIdx',             1), 1, N);
idx = 1;
for c = 1:numel(codes)
    for s = seeds
        for f = capacityFlags
            runCodes{idx} = codes{c};
            overrides(idx).gaSeed = s;
            overrides(idx).gaControlsMbsCapacity = f;
            idx = idx + 1;
        end
    end
end

% Pure-connectivity objective: drop the power-penalty (gamma) and FBS-share
% (fbsWeight) terms so fitness == normalized connected-user count (beta=1).
if pureConn
    [overrides.gamma]     = deal(0);
    [overrides.fbsWeight] = deal(0);
    [overrides.beta]      = deal(1.0);
end

study = ga_study.makeStudy(studyName, ...
    'codes',     runCodes, ...
    'overrides', overrides, ...
    'mode',      'list', ...
    'force',     true);

% --- Plan table: one row per run with its (code, seed, flag) triple --------
planCode = strings(study.numRuns, 1);
planSeed = zeros(study.numRuns, 1);
planFlag = zeros(study.numRuns, 1);
for i = 1:study.numRuns
    planCode(i) = string(study.runs{i}.code);
    planSeed(i) = study.runs{i}.overrides.gaSeed;
    planFlag(i) = double(study.runs{i}.overrides.gaControlsMbsCapacity);
end
plan = table((1:study.numRuns)', planCode, planSeed, planFlag, ...
    'VariableNames', {'run_index', 'code', 'gaSeed', 'gaControlsMbsCapacity'});
disp(plan);
fprintf('[ga_ablation_capacity] variant "%s" -> study "%s": %d planned runs\n', ...
    variant, study.name, study.numRuns);

if dryRun
    fprintf('[ga_ablation_capacity] dry run -- not executing the study.\n');
else
    ga_study.run(study);
end
