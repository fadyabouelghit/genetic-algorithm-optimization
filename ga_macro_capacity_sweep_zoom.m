%GA_MACRO_CAPACITY_SWEEP_ZOOM  Supplementary runs to zoom the gamma knee.
%
% The first campaign (macro_capacity_sweep_v1) saturated: the 1-1-1 carrier is
% ON through gamma 0.02 and fully OFF by 0.06, with no samples in between. This
% adds the missing knee points and runs the full 2-1-1 sweep. It REUSES the
% existing 1-1-1 ON runs at gamma {0, 0.02, 0.06, 0.08}, so only the new cells
% run here.
%
% Writes to a fresh macro_capacity_sweep_v1_<date> folder (same prefix), and the
% figure pollers merge ALL macro_capacity_sweep_v1_* folders into one curve.
%
% Dry run (print the list, train nothing):
%   dryRun = true; ga_macro_capacity_sweep_zoom

if ~exist('dryRun', 'var'); dryRun = false; end

seeds     = 1:15;
studyName = 'macro_capacity_sweep_v1';   % same prefix so the pollers merge it

% Per-code plan: {code, capacityON, gamma-list}
%   1-1-1 ON  -> only the 3 NEW knee points (0,0.02,0.06,0.08 already exist)
%   2-1-1 ON  -> full zoomed sweep (nothing ran yet)
%   *   OFF   -> baseline anchors aligned to the ON grid
plan = {
    '1-1-1', true,  [0.03 0.04 0.05]
    '1-1-1', false, [0 0.04 0.08]
    '2-1-1', true,  [0 0.02 0.03 0.04 0.05 0.06 0.08]
    '2-1-1', false, [0 0.04 0.08]
};

N = 0;
for r = 1:size(plan,1); N = N + numel(plan{r,3}) * numel(seeds); end

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
for r = 1:size(plan,1)
    code  = plan{r,1};
    capOn = plan{r,2};
    for g = plan{r,3}
        for s = seeds
            runCodes{idx} = code;
            overrides(idx).gaControlsMbsCapacity = capOn;
            overrides(idx).gamma  = g;
            overrides(idx).gaSeed = s;
            idx = idx + 1;
        end
    end
end

fprintf('[zoom] %d supplementary runs planned\n', N);

if dryRun
    for i = 1:N
        o = overrides(i);
        fprintf('%3d  %s  cap=%d  gamma=%.3f  seed=%2d\n', ...
            i, runCodes{i}, o.gaControlsMbsCapacity, o.gamma, o.gaSeed);
    end
    fprintf('[zoom] dry run -- nothing executed.\n');
else
    study = ga_study.makeStudy(studyName, ...
        'codes', runCodes, 'overrides', overrides, 'mode', 'list', 'force', true);
    fprintf('[zoom] launching "%s" (%d runs)\n', study.name, study.numRuns);
    ga_study.run(study);
end
