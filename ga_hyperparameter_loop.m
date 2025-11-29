% ============================
% Hyperparameter grids
% ============================
betaGrid       = [0.6 0.7 0.8 0.9 1.0];
gammaGrid      = [0.0 0.1 0.2 0.3 0.4];
fbsWeightGrid  = [0.4 0.6];
fbsExpGrid     = [1 2];
numFbsGrid     = [1 2];   % number of FBSs to test

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
