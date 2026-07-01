%GA_MACRO_CAPACITY_SWEEP  Gamma sweep x macro capacity-carrier ablation (LIST mode).
%
% Conference-paper-2 main campaign. One explicit run list (not a grid):
%
%   Capacity ON  (gaControlsMbsCapacity = true)  -> the GA may switch the macro
%       capacity carrier on; swept over the full gamma list. This is the
%       threshold curve: carriers drop off as the power price gamma rises.
%   Capacity OFF (gaControlsMbsCapacity = false) -> capacity carrier locked off,
%       FBS-only baseline; only a few gamma anchors (it is flat in gamma).
%
% Cost: the new power-aware cost with powerIncludesMacroCapacity = true, so an
% activated capacity carrier costs its 20 W. Pure connectivity-vs-power shaping
% (beta = 1, fbsWeight = 0). FBS band locked to coverage. No .mat dumps.
%
% Scenarios: 1-1-1 (1 FBS, 1 macro) and 2-1-1 (2 FBS, 1 macro).
% Seeds: 15 per cell. pop 15, gens 20.
%
% Dry run (print the full run list, train nothing):
%   dryRun = true; ga_macro_capacity_sweep

if ~exist('dryRun', 'var'); dryRun = false; end

% ---- Knobs (edit here to trim) ----------------------------------------
codes_all  = {'1-1-1', '2-1-1'};
gammas_on  = [0 0.02 0.06 0.08 0.10 0.15 0.20];   % capacity ON: full sweep
gammas_off = [0 0.08 0.20];                            % capacity OFF: flat baseline anchors
seeds      = 1:15;
studyName  = 'macro_capacity_sweep_v1';

% ---- Build the explicit run list --------------------------------------
nOnCells  = numel(codes_all) * numel(gammas_on);
nOffCells = numel(codes_all) * numel(gammas_off);
N = (nOnCells + nOffCells) * numel(seeds);

runCodes  = cell(1, N);
overrides = repmat(struct( ...
    'randomizeGA',                true, ...
    'gaSeed',                     NaN, ...
    'populationSize',             15, ...
    'numGenerations',             20, ...
    'gaControlsMbsCapacity',      true, ...
    'gaControlsFbsBand',          false, ...
    'powerIncludesMacroCapacity', true, ...
    'targetIdx',                  1, ...
    'gamma',                      0, ...
    'fbsWeight',                  0, ...
    'beta',                       1), 1, N);

idx = 1;
for c = 1:numel(codes_all)
    % Capacity ON -- full gamma sweep
    for g = gammas_on
        for s = seeds
            runCodes{idx} = codes_all{c};
            overrides(idx).gaControlsMbsCapacity = true;
            overrides(idx).gamma  = g;
            overrides(idx).gaSeed = s;
            idx = idx + 1;
        end
    end
    % Capacity OFF -- baseline anchors only
    for g = gammas_off
        for s = seeds
            runCodes{idx} = codes_all{c};
            overrides(idx).gaControlsMbsCapacity = false;
            overrides(idx).gamma  = g;
            overrides(idx).gaSeed = s;
            idx = idx + 1;
        end
    end
end

fprintf('[macro_capacity_sweep] %d runs: %d ON cells + %d OFF cells x %d seeds\n', ...
    N, nOnCells, nOffCells, numel(seeds));

% ---- Dry run or launch ------------------------------------------------
if dryRun
    for i = 1:N
        o = overrides(i);
        fprintf('%3d  %s  cap=%d  gamma=%.2f  seed=%2d\n', ...
            i, runCodes{i}, o.gaControlsMbsCapacity, o.gamma, o.gaSeed);
    end
    fprintf('[macro_capacity_sweep] dry run -- nothing executed.\n');
else
    study = ga_study.makeStudy(studyName, ...
        'codes', runCodes, 'overrides', overrides, 'mode', 'list', 'force', true);
    fprintf('[macro_capacity_sweep] launching study "%s" (%d runs)\n', ...
        study.name, study.numRuns);
    ga_study.run(study);
end
