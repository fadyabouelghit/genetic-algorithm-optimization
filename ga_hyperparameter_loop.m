% ============================
% Hyperparameter grids
% ============================
betaGrid       = [0.6 0.7 0.8 0.9 1.0];
gammaGrid      = [0.0 0.1 0.2 0.3 0.4];
fbsWeightGrid  = [0.4 0.6];
fbsExpGrid     = [1 2];
numFbsGrid     = [1 2];   

% Base bounds for a single FBS (x, y, height, power, on/off)
baseBounds = [0   W;      % x
              0   H;      % y
              20  150;    % height
              7   10.5;   % power
              0   1];     % power status or binary var

% Directory for log files (per run)
logDir = './logs_grid';
if ~exist(logDir, 'dir')
    mkdir(logDir);
end


% ============================
% Full grid search
% ============================
results = [];   % optional: in-memory summary
runIdx  = 0;

for nFbs = numFbsGrid
    for betaVal = betaGrid
        for gammaVal = gammaGrid
            for wVal = fbsWeightGrid
                for expVal = fbsExpGrid

                    runIdx = runIdx + 1;

                    % --- Update GA / fitness parameters for this run ---
                    params.numBS = nFbs;
                    params.bounds = repmat(baseBounds, nFbs, 1);

                    params.fitnessWeights.beta       = betaVal;
                    params.fitnessWeights.gamma      = gammaVal;
                    params.fitnessWeights.fbsWeight  = wVal;
                    params.fitnessWeights.fbsExponent= expVal;

                    % adjust looping configs if required with different FBS
                    % numbers
                    % if nFbs == 2
                        % params.populationSize = 20;
                    % else
                        % params.populationSize = 15;
                    % end

                    % --- Unique log file name per run ---
                    logName = sprintf( ...
                        'log_nFbs%d_beta%.2f_gamma%.2f_w%.2f_exp%.1f_run%03d.mat', ...
                        nFbs, betaVal, gammaVal, wVal, expVal, runIdx);

                    params.logFile = fullfile(logDir, logName);

                    fprintf('Run %d: nFbs=%d, beta=%.2f, gamma=%.2f, w=%.2f, exp=%.1f\n', ...
                            runIdx, nFbs, betaVal, gammaVal, wVal, expVal);

                    % --- Call optimizer (logging is enabled in params) ---
                    [bestInd, bestFit, history] = optimizeBaseStation( ...
                        fbsAntenna, containsMbs, antennaObjectMbs, mbs_params, params);

                    % --- Optional: store a compact summary in memory ---
                    results(runIdx).nFbs        = nFbs;
                    results(runIdx).beta        = betaVal;
                    results(runIdx).gamma       = gammaVal;
                    results(runIdx).fbsWeight   = wVal;
                    results(runIdx).fbsExponent = expVal;
                    results(runIdx).bestInd     = bestInd;
                    results(runIdx).bestFit     = bestFit;
                    % You can also store summary metrics extracted from history here

                end
            end
        end
    end
end

% Optionally save the summary struct
save(fullfile(logDir, 'grid_results_summary.mat'), 'results');
