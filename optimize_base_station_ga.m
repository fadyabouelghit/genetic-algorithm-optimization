% Main GA script for Flying Base Station Optimization

clearvars;
clear;
clc;

% ---------- Experiment I/O setup ----------
% timestamp = datestr(now,'yyyymmdd_HHMMSS');
% outdir = fullfile(pwd, ['exp_alpha_', timestamp]);
% if ~exist(outdir, 'dir'); mkdir(outdir); end
% figdir = fullfile(outdir, 'figs'); if ~exist(figdir,'dir'); mkdir(figdir); end

% Setup antennas and MBS parameters
fbsAntenna = setup_antenna();
mbsAntenna = setup_antenna();

numMbs = 1;  
W = 2000; H = 1500;      
margin = 100;            % edge distance threshold
ISD = 500;               % inter-site distance
[xs, ys] = generate_hex_sites(W, H, ISD, margin, numMbs);
mbs_height   = 25;
mbs_power    = 20;
[mbs_params, antennaObjectMbs, containsMbs, numMbs] = ...
    pack_mbs_params(xs, ys, mbs_height, mbs_power, mbsAntenna);
tempForX = mbs_params(1,:);
mbs_params(1,:) = mbs_params(2,:);
mbs_params(2,:) = tempForX;

subset = struct('xmin', 0, 'xmax', W, ...
                'ymin', 0, 'ymax', H);

cache = precompute_mbs_power_maps( ...
    mbs_params, ...
    antennaObjectMbs, ...
    subset, ...
    '3GPP_38.901_UMa_LOS', ...   % scenario
    'quick', ...                  % mode
    1.5, ...                      % UE height
    './cache_mbs_maps' ...        % folder for on-disk cache
);

disp(['Generated ', num2str(numMbs), ' MBSs']);
% mbs_x = [300];       
% mbs_y = [350];     
% mbs_height = [25];        
% mbs_power = [20]; 
% mbs_params = [mbs_x; mbs_y; mbs_height; mbs_power];


numBS = 2; % number of FBSs
containsMbs = 1;
% Optimization parameters
params = struct(...
    'enableLogging', true, ...
    'plotTrajectory', false, ...
    'initialPopulationSize', 30, ...
    'populationSize', 15, ...
    'numGenerations', 10, ...
    'crossoverProb', 0.3, ...
    'mutationProb', 0.8, ...
    'mutationScale', 0.2, ...
    'fitnessWeights', struct('beta', 1, 'gamma', 0, 'fbsWeight', 0.4, 'fbsExponent', 1), ...
    'maxUsers', 1000, ...
    'sinrThreshold', 5, ...
    'logFile', '', ...
    'numBS', numBS, ...
    'bounds', repmat([0 W; 0 H; 20 150; 7 10.5; 0 1], numBS, 1), ...
    'spaceLimit', [W,H], ...
    'mbsCache', cache, ...
    'verbose', 1 ...
);

runMoga = false;
[bestInd, bestFit, history] = optimizeBaseStation( ...
    fbsAntenna, containsMbs, antennaObjectMbs, mbs_params, params);

% ---- Multi-objective GA (connectivity vs. power) ----
if runMoga
    fprintf('\n=== Launching Multi-Objective GA (MOGA) ===\n');
    mogaParams = params;
    mogaParams.enableLogging = false;
    mogaParams.numGenerations = max(20, params.numGenerations);
    mogaParams.tournamentSize = 3;
    mogaParams.populationSize = 150; 
    mogaParams.mutationScale = 0.5;
    mogaParams.crossoverProb = 0.3;
    mogaParams.mutationProb = 0.7;
    mogaParams.FbsWeight = 0;
    mogaParams.beta = 1;
    mogaParams.gamma = 0;

    if mogaParams.initialPopulationSize < mogaParams.populationSize
        mogaParams.initialPopulationSize = mogaParams.populationSize;
    end
    [paretoFront, mogaHistory] = optimizeBaseStationMoga( ...
        fbsAntenna, containsMbs, antennaObjectMbs, mbs_params, mogaParams);

    if ~isempty(paretoFront.users)
        paretoTable = table( ...
            paretoFront.users(:), ...
            paretoFront.power(:), ...
            paretoFront.avgRate(:), ...
            'VariableNames', {'ConnectedUsers', 'TransmittedPowerW', 'AvgRatebpsHz'});
        disp('Final Pareto-optimal trade-offs (users vs. power):');
        disp(paretoTable);
    else
        warning('Pareto front is empty. Check GA settings or SINR configuration.');
    end
end

% % ---- RL knobs / action sets ----
% params.defaultSigmaScales = [0.2, 0.2, 0.2, 0.2];   % X, Y, Z, Power
% params.defaultMeanBiases  = [0, 0, 0, 0];   % X, Y, Z, Power
% 
% % One enable flag per continuous variable (X, Y, Z, Power)
% params.rl.enableSigma = [false, false, false, false];  % agent controls sigma scale per variable
% params.rl.enableMean  = [false, false, false, false];  % agent controls mean bias per variable
% 
% params.rl.enableCrossoverProb  = true;   % agent can change crossoverProb
% params.rl.enableMutationProb   = true;   % agent can change mutationProb
% 
% % Action sets 
% params.rl.sigma_scales = { [0.01 0.2], [0.01 0.2], [], [] };
% params.rl.mean_shifts  = { [], [], [], [] };
% params.rl.crossover_set = [0.3 0.5 0.7];
% params.rl.mutation_set  = [0.3 0.5 0.7];
% 
% params.rl.q.alpha        = 0.9;   % learning rate
% params.rl.q.gamma        = 0.3;   % discount factor
% params.rl.q.epsilon      = 0.8;   % initial exploration
% params.rl.q.epsilonDecay = 0.96;  % epsilon decay per generation
% params.rl.q.minEpsilon   = 0.05;   % exploration floor
% 
% 
% % ---------- Sweep settings ----------
% alphas = 0.1:0.1:1.0;
% nRuns  = 5;
% 
% % ---------- Storage ----------
% run_alpha   = [];        % double
% run_id      = [];        % double (1..nRuns)
% run_seed    = [];        % double
% run_bestFit = [];        % double
% run_runtime = [];        % double
% run_bestInd_str = {};    % cellstr (stringified vector)
% run_logfile = {};        % cellstr
% run_status  = {};        % 'ok' | 'error'
% run_errmsg  = {};        % error text (if any)
% 
% % ---------- Live plotting setup ----------
% % 1) Convergence figure (each run’s bestFitness over generations)
% figConv = figure('Name','Convergence per Run','NumberTitle','off');
% clf(figConv); hold on; grid on;
% xlabel('Generation'); ylabel('Best Fitness'); title('Convergence per Run');
% 
% % 2) BestFit vs Alpha (filled as stats are computed later)
% figAlpha = figure('Name','BestFit vs Alpha','NumberTitle','off');
% clf(figAlpha); hold on; grid on;
% xlabel('\alpha'); ylabel('Best Fitness'); title('Per-Alpha Summary');
% 
% % ---------- Main sweep ----------
% fprintf('=== Starting alpha sweep (%d alphas × %d runs) ===\n', numel(alphas), nRuns);
% 
% for ai = 1:numel(alphas)
%     a = alphas(ai);
%     fprintf('\n=== Alpha %.2f (%d/%d) ===\n', a, ai, numel(alphas));
% 
%     % Per-alpha accumulators for live summary
%     alpha_bestFits = zeros(nRuns,1);
%     alpha_runtimes = zeros(nRuns,1);
%     alpha_ok       = false(nRuns,1);
% 
%     for r = 1:nRuns
%         seed = ai*1000 + r;    % deterministic per (alpha, run)
%         rng(seed);
% 
%         % Choose where alpha lives: params.alpha vs params.rl.q.alpha
%         params_cur = params;
%         if isfield(params, 'alpha')
%             params_cur.alpha = a;
%         elseif isfield(params,'rl') && isfield(params.rl,'q') && isfield(params.rl.q,'alpha')
%             params_cur.rl.q.alpha = a;
%         else
%             params_cur.alpha = a; % add top-level alpha if neither exists
%         end
% 
%         % Build logfile & start diary
%         logFile = fullfile(outdir, sprintf('alpha_%0.2f_run_%02d.log', a, r));
%         diary(logFile);
%         fprintf('>>> Alpha=%.2f | Run=%02d/%02d | Seed=%d\n', a, r, nRuns, seed);
% 
%         t0 = tic;
%         status = 'ok'; errmsg = '';
%         bestInd = []; bestFit = NaN; history = struct();
%         try
%             [bestInd, bestFit, history] = OptimizeBaseStationRl_prob( ...
%                 fbsAntenna, containsMbs, antennaObjectMbs, mbs_params, params_cur);
%         catch ME
%             status = 'error';
%             errmsg = ME.message;
%             fprintf(2, 'ERROR: %s\n', errmsg);
%         end
%         elapsed = toc(t0);
% 
%         % Diagnostics
%         if isfield(history,'bestFitness')
%             bf_hist = history.bestFitness(:).';
%             fprintf('    Generations=%d | final bestFit=%.6g | mean(bestFitness)=%.6g | runtime=%.2fs\n', ...
%                 numel(bf_hist), bestFit, mean(bf_hist,'omitnan'), elapsed);
% 
%             % Live convergence plot
%             figure(figConv);
%             plot(bf_hist, '-o', 'DisplayName', sprintf('\\alpha=%.2f run=%d', a, r));
%             legend('-dynamiclegend'); drawnow;
%         else
%             fprintf('    No history.bestFitness captured. final bestFit=%.6g | runtime=%.2fs\n', bestFit, elapsed);
%         end
%         diary off;
% 
%         % Save per-run MAT snapshot (even on error, for forensics)
%         save(fullfile(outdir, sprintf('alpha_%0.2f_run_%02d.mat', a, r)), ...
%             'a','seed','bestInd','bestFit','history','elapsed','params_cur','status','errmsg');
% 
%         % Append to aggregate arrays
%         run_alpha(end+1,1) = a;                 %#ok<SAGROW>
%         run_id(end+1,1)    = r;                 %#ok<SAGROW>
%         run_seed(end+1,1)  = seed;              %#ok<SAGROW>
%         run_bestFit(end+1,1) = bestFit;         %#ok<SAGROW>
%         run_runtime(end+1,1) = elapsed;         %#ok<SAGROW>
%         run_bestInd_str{end+1,1} = mat2str(bestInd(:).'); %#ok<SAGROW>
%         run_logfile{end+1,1}     = logFile;     %#ok<SAGROW>
%         run_status{end+1,1}      = status;      %#ok<SAGROW>
%         run_errmsg{end+1,1}      = errmsg;      %#ok<SAGROW>
% 
%         % Per-alpha accumulators
%         alpha_bestFits(r) = bestFit;
%         alpha_runtimes(r) = elapsed;
%         alpha_ok(r)       = strcmp(status,'ok');
%     end
% 
%     % Per-alpha console summary
%     okFits = alpha_bestFits(alpha_ok);
%     if ~isempty(okFits)
%         fprintf('--- Alpha=%.2f summary: mean=%.6g | std=%.6g | min=%.6g | max=%.6g | median=%.6g | mean_rt=%.2fs | ok=%d/%d ---\n', ...
%             a, mean(okFits), std(okFits,0), min(okFits), max(okFits), median(okFits), ...
%             mean(alpha_runtimes(alpha_ok)), sum(alpha_ok), nRuns);
%     else
%         fprintf('--- Alpha=%.2f summary: no successful runs (ok=0/%d) ---\n', a, nRuns);
%     end
% end
% 
% % ---------- Build per-run table and write CSV ----------
% runsTbl = table(run_alpha, run_id, run_seed, run_bestFit, run_runtime, run_bestInd_str, run_logfile, run_status, run_errmsg, ...
%     'VariableNames', {'alpha','run','seed','bestFit','runtime_sec','bestInd','logfile','status','errmsg'});
% writetable(runsTbl, fullfile(outdir, 'runs.csv'));
% 
% % ---------- Compute per-alpha statistics ----------
% alphas_u = unique(runsTbl.alpha);
% stats = table('Size',[numel(alphas_u) 10], ...
%     'VariableTypes', {'double','double','double','double','double','double','double','double','double','double'}, ...
%     'VariableNames', {'alpha','mean_bestFit','std_bestFit','min_bestFit','max_bestFit','median_bestFit','mean_runtime_sec','ok_runs','total_runs','error_rate'});
% 
% for i = 1:numel(alphas_u)
%     a = alphas_u(i);
%     rows = runsTbl.alpha == a & strcmp(runsTbl.status,'ok');
%     bf = runsTbl.bestFit(rows);
%     rt = runsTbl.runtime_sec(rows);
%     stats.alpha(i)            = a;
%     stats.mean_bestFit(i)     = mean(bf,'omitnan');
%     stats.std_bestFit(i)      = std(bf,0,'omitnan');
%     stats.min_bestFit(i)      = min(bf,[],'omitnan');
%     stats.max_bestFit(i)      = max(bf,[],'omitnan');
%     stats.median_bestFit(i)   = median(bf,'omitnan');
%     stats.mean_runtime_sec(i) = mean(rt,'omitnan');
%     stats.ok_runs(i)          = sum(rows);
%     stats.total_runs(i)       = sum(runsTbl.alpha == a);
%     stats.error_rate(i)       = 1 - stats.ok_runs(i) / stats.total_runs(i);
% end
% 
% writetable(stats, fullfile(outdir, 'alpha_stats.csv'));
% save(fullfile(outdir, 'summary.mat'), 'runsTbl', 'stats', 'alphas', 'nRuns', 'timestamp', 'outdir');
% 
% % ---------- Plot per-alpha summary (bestFit vs alpha) ----------
% figure(figAlpha); clf(figAlpha); hold on; grid on;
% xlabel('\alpha'); ylabel('Best Fitness'); title('Per-Alpha Summary');
% 
% % Scatter all runs (faded), then overlay stats (bold)
% scatter(runsTbl.alpha, runsTbl.bestFit, 36, 'filled', 'MarkerFaceAlpha', 0.25);
% plot(stats.alpha, stats.mean_bestFit, '-o', 'LineWidth', 2, 'DisplayName', 'Mean bestFit');
% plot(stats.alpha, stats.median_bestFit, '-s', 'LineWidth', 2, 'DisplayName', 'Median bestFit');
% legend('All runs','Mean','Median','Location','best');
% drawnow;
% 
% % Save figures
% saveas(figConv, fullfile(figdir, 'convergence_per_run.png'));
% saveas(figAlpha, fullfile(figdir, 'bestFit_vs_alpha.png'));
% 
% % ---------- Console recap ----------
% disp('--- Sweep complete ---');
% disp(stats);
% fprintf('Artifacts saved under: %s\n', outdir);



% [bestInd, bestFit, history] = OptimizeBaseStationRl_prob( ...
%     fbsAntenna, containsMbs, antennaObjectMbs, mbs_params, params);

% [bestInd, bestFit, history] = OptimizeBaseStationRl( ...
%     fbsAntenna, containsMbs, antennaObjectMbs, mbs_params, params);
