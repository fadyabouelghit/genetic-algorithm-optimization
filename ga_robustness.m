%GA_ROBUSTNESS  GA statistical-robustness campaign (many independent seeds).
%
% Re-runs the GA from scratch nReps times per scenario code under
% randomizeGA=true (each repeat draws an independent shuffled seed), so the
% spread of best-fitness / connected-user outcomes measures how robust the
% search is to GA randomness alone. The user map is fixed in SINREvaluation,
% so only the GA's own randomness (init pool + operators) varies across reps.
%
% This is the dedicated driver for the campaign that previously lived only as
% an inline block in ga_study_example.m (study 5). Pick the variant below.
%
%   'original'       : as-published sweep over the Y=1 codes with both Z
%                      configs ({'1-1-1','1-1-2','2-1-1','2-1-2'}) and default
%                      fitness shaping (fbsWeight 0.4, gamma from Z, beta 0.8).
%   'connonly_4scen' : pure connectivity (fbsWeight 0, gamma 0, beta 1.0)
%                      over the macro x FBS scenario grid
%                      (1-1-1, 2-1-1, 1-2-1, 2-2-1 = {1,2} macro x {1,2} FBS).
%
% Dry run (print the plan, train nothing):
%   dryRun = true; ga_robustness

if ~exist('dryRun', 'var'); dryRun = false; end

variant = 'connonly_4scen';     % 'original' | 'connonly_4scen'
nReps   = 31;                   % independent random GA runs per code

switch variant
    case 'original'
        baseCodes = {'1-1-1', '1-1-2', '2-1-1', '2-1-2'};
        studyName = 'robustness_v1';
        pureConn  = false;
    case 'connonly_4scen'
        baseCodes = {'1-1-1', '2-1-1', '1-2-1', '2-2-1'};
        studyName = 'robustness_connonly_4scen_v1';
        pureConn  = true;
    otherwise
        error('ga_robustness:badVariant', 'Unknown variant "%s".', variant);
end

overrides = struct( ...
    'randomizeGA',    true, ...
    'populationSize', 25, ...
    'numGenerations', 25);

% Pure-connectivity objective: drop the power-penalty (gamma) and FBS-share
% (fbsWeight) terms so fitness == normalized connected-user count (beta=1).
if pureConn
    overrides.gamma     = 0;
    overrides.fbsWeight = 0;
    overrides.beta      = 1.0;
end

study = ga_study.makeStudy(studyName, ...
    'codes',     repmat(baseCodes, 1, nReps), ...
    'overrides', overrides, ...
    'force',     true);
fprintf('[ga_robustness] variant "%s" -> "%s": %d planned runs (%d codes x %d reps)\n', ...
    variant, study.name, study.numRuns, numel(baseCodes), nReps);

if dryRun
    fprintf('[ga_robustness] dry run -- not executing the study.\n');
else
    ga_study.run(study);
end
