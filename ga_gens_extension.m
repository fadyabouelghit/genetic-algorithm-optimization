%GA_GENS_EXTENSION  Confirm the convergence plateau holds beyond 25 gens.
%
% The 31-seed robustness study (pop 25, 25 gens) showed median plateau at
% generation 13-16, but 19-35%% of runs were still improving in the final
% generations. This small study runs 40 generations at otherwise identical
% settings to check whether fitness resumes climbing after generation 25.
% Analysis: overlay median/IQR best-fitness curves; if gens 26-40 add less
% than ~0.5%% median fitness, the 25-generation budget is confirmed.
%
% Dry run (print the plan, train nothing):
%   dryRun = true; ga_gens_extension

if ~exist('dryRun', 'var'); dryRun = false; end

% ---- Study variant ----------------------------------------------------
%   'original'       : as-published gens extension, default fitness shaping
%                      (fbsWeight 0.4, gamma from Z code, beta 0.8); 2 codes.
%   'connonly_4scen' : pure connectivity (fbsWeight 0, gamma 0, beta 1.0)
%                      over the macro x FBS scenario grid
%                      (1-1-1, 2-1-1, 1-2-1, 2-2-1 = {1,2} macro x {1,2} FBS).
variant = 'connonly_4scen';     % 'original' | 'connonly_4scen'

seeds = 1:10;
switch variant
    case 'original'
        codes     = {'1-1-1', '2-1-1'};
        studyName = 'gens_extension_v1';
        pureConn  = false;
    case 'connonly_4scen'
        codes     = {'1-1-1', '2-1-1', '1-2-1', '2-2-1'};
        studyName = 'gens_extension_connonly_4scen_v1';
        pureConn  = true;
    otherwise
        error('ga_gens_extension:badVariant', 'Unknown variant "%s".', variant);
end

N = numel(codes) * numel(seeds);
runCodes  = cell(1, N);
overrides = repmat(struct( ...
    'randomizeGA',           true, ...
    'gaSeed',                NaN, ...
    'populationSize',        25, ...
    'numGenerations',        40, ...
    'gaControlsMbsCapacity', false, ...
    'gaControlsFbsBand',     false, ...
    'targetIdx',             1), 1, N);
idx = 1;
for c = 1:numel(codes)
    for s = seeds
        runCodes{idx} = codes{c};
        overrides(idx).gaSeed = s;
        idx = idx + 1;
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
    'codes', runCodes, 'overrides', overrides, 'mode', 'list', 'force', true);
fprintf('[ga_gens_extension] variant "%s" -> "%s": %d planned runs\n', ...
    variant, study.name, study.numRuns);

if dryRun
    for i = 1:study.numRuns
        o = study.runs{i}.overrides;
        fprintf('%3d  %s  seed=%2d  gens=%2d\n', i, ...
            study.runs{i}.code, o.gaSeed, o.numGenerations);
    end
    fprintf('[ga_gens_extension] dry run -- not executing the study.\n');
else
    ga_study.run(study);
end
