%GA_POP_SWEEP  Population-size sensitivity, seed-paired across levels.
%
% Two framings, selected by budgetMatched below:
%   false (default): popLevels x fixed 25 generations — answers "what does a
%                    larger population buy per generation?"
%   true           : (pop, gens) pairs at a fixed ~500-evaluation budget —
%                    answers "what is the best allocation of fixed compute?"
%
% The SAME gaSeed list is reused at every level, so per-seed comparisons
% across levels are paired (same GA randomness, only the level differs).
% initialPopulationSize is tied to populationSize (no gen-0 oversampling)
% so the population effect is not confounded by a fixed oversampled pool.
% NOTE: this differs from earlier studies (which used initPop=30); compare
% levels within this study, not against historical runs.
%
% Dry run (print the plan, train nothing):
%   dryRun = true; ga_pop_sweep

if ~exist('dryRun', 'var'); dryRun = false; end
budgetMatched = false;

% ---- Study variant ----------------------------------------------------
%   'original'       : as-published sweep, default fitness shaping
%                      (fbsWeight 0.4, gamma from Z code, beta 0.8); 2 codes.
%   'connonly_4scen' : pure connectivity (fbsWeight 0, gamma 0, beta 1.0)
%                      over the macro x FBS scenario grid
%                      (1-1-1, 2-1-1, 1-2-1, 2-2-1 = {1,2} macro x {1,2} FBS).
variant = 'connonly_4scen';     % 'original' | 'connonly_4scen'

seeds = 1:15;
switch variant
    case 'original'
        codes    = {'1-1-1', '2-1-1'};
        suffix   = '';
        pureConn = false;
    case 'connonly_4scen'
        codes    = {'1-1-1', '2-1-1', '1-2-1', '2-2-1'};
        suffix   = '_connonly_4scen';
        pureConn = true;
    otherwise
        error('ga_pop_sweep:badVariant', 'Unknown variant "%s".', variant);
end

if budgetMatched
    levels = struct('pop', {10, 20, 25, 50}, 'gens', {50, 25, 20, 10});
    name   = ['pop_sweep_budget'   suffix '_v1'];
else
    levels = struct('pop', {10, 20, 30, 40}, 'gens', {25, 25, 25, 25});
    name   = ['pop_sweep_fixedgen' suffix '_v1'];
end

N = numel(codes) * numel(seeds) * numel(levels);
runCodes  = cell(1, N);
overrides = repmat(struct( ...
    'randomizeGA',           true, ...
    'gaSeed',                NaN, ...
    'populationSize',        NaN, ...
    'initialPopulationSize', NaN, ...
    'numGenerations',        NaN, ...
    'gaControlsMbsCapacity', false, ...
    'gaControlsFbsBand',     false, ...
    'targetIdx',             1), 1, N);
idx = 1;
for c = 1:numel(codes)
    for s = seeds
        for L = 1:numel(levels)
            runCodes{idx} = codes{c};
            overrides(idx).gaSeed                = s;
            overrides(idx).populationSize        = levels(L).pop;
            overrides(idx).initialPopulationSize = levels(L).pop;
            overrides(idx).numGenerations        = levels(L).gens;
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

study = ga_study.makeStudy(name, ...
    'codes', runCodes, 'overrides', overrides, 'mode', 'list', 'force', true);
fprintf('[ga_pop_sweep] variant "%s" -> "%s": %d planned runs\n', ...
    variant, study.name, study.numRuns);

if dryRun
    for i = 1:study.numRuns
        o = study.runs{i}.overrides;
        fprintf('%3d  %s  seed=%2d  pop=%2d  gens=%2d\n', i, ...
            study.runs{i}.code, o.gaSeed, o.populationSize, o.numGenerations);
    end
    fprintf('[ga_pop_sweep] dry run -- not executing the study.\n');
else
    ga_study.run(study);
end
