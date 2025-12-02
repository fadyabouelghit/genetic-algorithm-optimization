function [bestIndividual, bestFitness, history] = optimizeBaseStation(l, containsMbs, antennaObjectMbs, mbs_params, params)
% OPTIMIZEBASESTATION Genetic algorithm optimizer with enhanced diagnostics
%   Inputs:
%   l       - QuaDRiGa layout object
%   params  - Structure containing algorithm parameters:
%       populationSize    : Individuals per generation
%       numGenerations    : Number of generations
%       crossoverProb     : Crossover probability [0-1]
%       mutationProb      : Mutation probability [0-1]
%       bounds            : Parameter bounds matrix [min, max] per parameter
%       verbose           : Display level (0: none, 1: normal, 2: detailed)
%
%   Outputs:
%   bestIndividual      - Optimal parameters found [x, y, z, power]
%   bestFitness         - Best fitness value achieved
%   history             - Structure with optimization history

% Setup diagnostics
    defaultParams.verbose = 1;
    defaultParams.fitnessWeights = struct('beta', 1, 'gamma', 1, 'epsilon', 1e-3, 'fbsWeight', 0, 'fbsExponent', 1);
    defaultParams.initialPopulationSize = [];
    defaultParams.plotTrajectory = false;
    defaultParams.maxUsers = 1000;
    defaultParams.sinrThreshold = 5;
    defaultParams.enableLogging = true;
    defaultParams.enablePerformancePlotting = false;
    defaultParams.logFile = '';
    params = mergeParams(defaultParams, params);
if isempty(params.initialPopulationSize)
    params.initialPopulationSize = params.populationSize;
end
if params.initialPopulationSize < params.populationSize
    error('initialPopulationSize (%d) must be >= populationSize (%d).', ...
        params.initialPopulationSize, params.populationSize);
end

if params.enableLogging
    functionDir = fileparts(mfilename('fullpath'));
    logDir = fullfile(functionDir, 'logs');
    if ~exist(logDir, 'dir')
        mkdir(logDir);
    end
    if isempty(params.logFile)
        timestamp = datestr(now, 'yyyymmdd_HHMMSS');
        params.logFile = sprintf('ga_results_%s.csv', timestamp);
    elseif ~endsWith(params.logFile, '.csv', 'IgnoreCase', true)
        params.logFile = [params.logFile '.csv'];
    end
    [~, baseName, ext] = fileparts(params.logFile);
    if isempty(ext)
        ext = '.csv';
    end
    params.logFile = fullfile(logDir, [baseName ext]);
end
% 1 -> connectivity / 2 -> avg sum rate
targetIdx = 1; 

% Initialize history tracking
history = struct(...
    'bestFitness', zeros(params.numGenerations, 1), ...
    'avgFitness', zeros(params.numGenerations, 1), ...
    'stdFitness', zeros(params.numGenerations, 1), ...
    'minFitness', zeros(params.numGenerations,1), ...
    'maxFitness', zeros(params.numGenerations,1), ...
    'bestIndividuals', zeros(params.numGenerations, size(params.bounds, 1)), ...
    'crossovers', zeros(params.numGenerations, 1), ...
    'mutations', zeros(params.numGenerations, 1), ...
    'time', struct('total', 0, 'generation', zeros(params.numGenerations, 1)), ...
    'avgVariables', zeros(params.numGenerations, size(params.bounds,1)), ...  % Average value for each variable
    'currentAvgReward', zeros(params.numGenerations, 1) ...         % Current average reward (same as avgFitness but explicit)
);

totalTimer = tic;

if params.verbose > 0
    fprintf('\n=== Genetic Algorithm Optimization ===\n');
    fprintf('Population: %d, Generations: %d\n', params.populationSize, params.numGenerations);
    fprintf('Crossover: %.1f%%, Mutation: %.1f%%\n', params.crossoverProb*100, params.mutationProb*100);
    disp('Initializing population...');
end

filename = 'population_evolution.xlsx';
function headers = create_headers(n_fbs)
    headers = cell(1, n_fbs*5 + 1);
    headers{1} = 'Generation';
    for bs = 1:n_fbs
        headers{(bs-1)*5+2} = sprintf('BS%d_X', bs);
        headers{(bs-1)*5+3} = sprintf('BS%d_Y', bs);
        headers{(bs-1)*5+4} = sprintf('BS%d_Z', bs);
        headers{(bs-1)*5+5} = sprintf('BS%d_Power', bs);
        headers{(bs-1)*5+6} = sprintf('BS%d_Power_status', bs);
    end
end


population = initializePopulation_uniform(params.initialPopulationSize, params.bounds, params.numBS);

evalParams = params.fitnessWeights;
evalParams.maxUsers = params.maxUsers;
evalParams.sinrThreshold = params.sinrThreshold;

if params.initialPopulationSize > params.populationSize
    [initialFitness, ~] = evaluatePopulation(l, population, params.verbose, params.numBS, ...
        params.spaceLimit ,containsMbs, mbs_params, antennaObjectMbs, params.bounds, ...
        params.mbsCache, targetIdx, evalParams);
    [~, sortIdx] = sort(initialFitness, 'descend');
    keepCount = params.populationSize;
    population = population(sortIdx(1:keepCount), :);
else
    population = population(1:params.populationSize, :);
end

trajectoryPlot.enabled = params.plotTrajectory;
if trajectoryPlot.enabled && params.numBS ~= 1
    warning('plotTrajectory currently supports only numBS == 1. Disabling plot.');
    trajectoryPlot.enabled = false;
end
if trajectoryPlot.enabled
    trajectoryPlot.maxGenerations = params.numGenerations;
    trajectoryPlot.fig = figure('Name', 'FBS Population Trajectory', 'NumberTitle', 'off');
    trajectoryPlot.ax = axes('Parent', trajectoryPlot.fig);
    hold(trajectoryPlot.ax, 'on');
    grid(trajectoryPlot.ax, 'on');
    xlabel(trajectoryPlot.ax, 'X (m)');
    ylabel(trajectoryPlot.ax, 'Y (m)');
    zlabel(trajectoryPlot.ax, 'Z (m)');
    title(trajectoryPlot.ax, 'Population Trajectory (Gen 0)');
    trajectoryPlot.colormap = parula(max(2, params.numGenerations));
    colormap(trajectoryPlot.ax, trajectoryPlot.colormap);
    trajectoryPlot.cb = colorbar(trajectoryPlot.ax);
    trajectoryPlot.cb.Label.String = 'Generation';
    caxis(trajectoryPlot.ax, [1 params.numGenerations]);
end

headers = create_headers(params.numBS);
initialData = [zeros(params.populationSize,1) population]; % Gen 0
% writetable(array2table(initialData, 'VariableNames', headers), filename, 'WriteMode', 'overwrite');

mbs_x = mbs_params(1,:);
mbs_y = mbs_params(2,:);
mbs_height = mbs_params(3,:);
mbs_power = mbs_params(4,:);

globalBestFitness = -inf;
globalBestIndividual = population(1,:);


for gen = 1:params.numGenerations
    genTimer = tic;
    
    % Evaluate fitness
    if params.verbose > 0
        fprintf('\nGeneration %d/%d:\n', gen, params.numGenerations);
        if params.verbose > 1
            disp('Evaluating population fitness...');
        end
    end
    
    % taregtIdx = 1; % optimize for connected users
    % taregtIdx = 2; % optimize for sum rate
    evalParams = params.fitnessWeights;
    evalParams.maxUsers = params.maxUsers;
    evalParams.sinrThreshold = params.sinrThreshold;
    [fitness, evalDetails] = evaluatePopulation(l, population, params.verbose, params.numBS, params.spaceLimit ,containsMbs, mbs_params, antennaObjectMbs, params.bounds, params.mbsCache, targetIdx, evalParams);
    
    if trajectoryPlot.enabled
        updateTrajectoryPlot(trajectoryPlot, population, gen, params.bounds);
    end
    
    % Track population statistics
    history.bestFitness(gen) = max(fitness);
    history.maxFitness(gen) = max(fitness);
    history.minFitness(gen) = min(fitness);
    history.avgFitness(gen) = mean(fitness);
    history.stdFitness(gen) = std(fitness);

    [bestFitness, bestIdx] = max(fitness);
    bestIndividual = population(bestIdx, :);
    history.bestIndividuals(gen, :) = bestIndividual;
    
    if bestFitness > globalBestFitness
        globalBestFitness = bestFitness;
        globalBestIndividual = bestIndividual;
    end
    
    history.avgVariables(gen, :) = mean(population, 1);
    
    
    % Display generation summary
    if params.verbose > 0
        fprintf('Fitness: Best=%.3f, Avg=%.2f, Std=%.2f\n', ...
            history.bestFitness(gen), history.avgFitness(gen), history.stdFitness(gen));
            
        for fb = 1:params.numBS
            startIdx = (fb-1)*5 + 1;
            
            coords = bestIndividual(startIdx : startIdx+2);
            coordCells = cellstr(num2str(coords', '%g'));
            coordStr = strjoin(coordCells, ', ');

            powerVal = bestIndividual(startIdx+3);            
            binaryVal = bestIndividual(startIdx+4);

            fprintf('Best Individual (FBS %d): [%s] Power: %.1f, Power Status: %d\n', ...
                    fb, coordStr, powerVal, binaryVal);
        end
        
        if params.verbose > 1
            fprintf('Fitness values:\n');
            disp(fitness');
        end
    end
    
    % Create new population
    eliteIndividual = bestIndividual;
    newPopulation = [eliteIndividual; zeros(params.populationSize-1, size(params.bounds, 1))];
    crossoverCount = 0;
    mutationCount = 0;
    
    % Generate offspring
    for i = 2:params.populationSize
        % Parent selection
        parents = selectParents(population, fitness, 2);
        
        % Crossover
%         [child1, child2, crossoverFlag] = crossover(parents(1,:), parents(2,:), params.crossoverProb, params.bounds);
%         [child1, child2, crossoverFlag] = crossover_sbx(parents(1,:), parents(2,:), params.crossoverProb, params.bounds);
        [child1, child2, crossoverFlag] = crossover_blend(parents(1,:), parents(2,:), params.crossoverProb, params.bounds);
        crossoverCount = crossoverCount + crossoverFlag;
        
        % Mutation
        [child, mutationFlags] = mutate(child1, params.mutationProb, params.bounds, params.mutationScale, gen, params, []);
        mutationCount = mutationCount + sum(mutationFlags);
        
        % Add to new population
        newPopulation(i,:) = child;
        
     end
    
    % Update operator statistics
    history.crossovers(gen) = crossoverCount;
    history.mutations(gen) = mutationCount;
    history.time.generation(gen) = toc(genTimer);
    
    if params.verbose > 0
        fprintf('Operators: %d crossovers, %d mutations\n', crossoverCount, mutationCount);
        fprintf('Generation time: %.1f seconds\n', history.time.generation(gen));
    end
    
    population = newPopulation;
    genData = [gen*ones(params.populationSize,1) newPopulation];
    % writetable(array2table(genData, 'VariableNames', headers), filename, 'WriteMode', 'append');
    % population
    
    % Early stopping check
    % if gen > 8 && std(history.bestFitness(gen-7:gen)) < 0.01
    %     if params.verbose > 0
    %         fprintf('\nEarly stopping: No improvement in last 8 generations\n');
    %     end
    %     break;
    % end
end

history.time.total = toc(totalTimer);
bestFitness = globalBestFitness;
bestIndividual = globalBestIndividual;

    [~, ~, numUsers, transmittedPower, avg_rate_connected_bpsHz, fbsUsers, mbsUsers] = SINREvaluation( ...
        l, bestIndividual(5:5:end), bestIndividual(1:5:end), bestIndividual(2:5:end), bestIndividual(3:5:end), params.numBS, ...
        bestIndividual(4:5:end), mbs_y, mbs_x, mbs_height, mbs_power, ...
        0, params.spaceLimit(1), 0, params.spaceLimit(2), params.maxUsers, params.sinrThreshold, containsMbs, antennaObjectMbs, params.mbsCache);

if params.verbose > 0
    fprintf('\n=== Optimization Complete ===\n');
    fprintf('Total time: %.1f minutes\n', history.time.total/60);
    fprintf('Best Fitness: %.2f users\n', bestFitness);
    
    fprintf('The best individual performance: \n') 
    fprintf('Total Connected Users: %d\n', numUsers);
    fprintf(' - FBS-connected Users: %d\n', fbsUsers);
    fprintf(' - MBS-connected Users: %d\n', mbsUsers);
    fprintf('Total Transmitted Power: %.2f W\n', transmittedPower);
    fprintf('Avg. Sum Rate: %.2f bps/Hz\n', avg_rate_connected_bpsHz);

    
    fprintf('Final Parameters:\n');
    
    % Create parameter table for multiple BS
    numBS = params.numBS;
    paramNames = cell(5*numBS, 1);
    for bs = 1:numBS
        paramNames{(bs-1)*5 + 1} = sprintf('BS%d X (m)', bs);
        paramNames{(bs-1)*5 + 2} = sprintf('BS%d Y (m)', bs);
        paramNames{(bs-1)*5 + 3} = sprintf('BS%d Z (m)', bs);
        paramNames{(bs-1)*5 + 4} = sprintf('BS%d Power (W)', bs);
        paramNames{(bs-1)*5 + 5} = sprintf('BS%d Power Status', bs);
    end
    
    disp(array2table(bestIndividual', ...
        'VariableNames', {'Value'}, ...
        'RowNames', paramNames));
    
    if params.enablePerformancePlotting
    % Enhanced visualization
    figure('Name', 'Optimization Results', 'Position', [100 100 1200 800]);
    
    % Fitness Progression
    subplot(2,2,1);
    hold on;
    
    % Ensure everything is a row vector
    gens = 1:gen;
    avg = history.avgFitness(1:gen)';
    stddev = history.stdFitness(1:gen)';
    best = history.bestFitness(1:gen)';
    
    % Shaded area for average ± std
    fill([gens, fliplr(gens)], ...
         [avg + stddev, fliplr(avg - stddev)], ...
         [1 0.7 0.7], ...
         'EdgeColor', 'none', ...
         'FaceAlpha', 0.3, ...
         'DisplayName', '±1 Std Dev');
    
    % Plot lines
    plot(gens, best, 'b-o', 'DisplayName', 'Best Fitness');
    plot(gens, avg,  'r--', 'LineWidth', 1.5, 'DisplayName', 'Average Fitness');

    legend('Location', 'best');
    xlabel('Generation');
    ylabel('Connected Users');
    title('Fitness Progression');
    grid on;

    % Operator Activity
    subplot(2,2,2);
    bar([history.crossovers(1:gen), history.mutations(1:gen)], 'stacked');
    legend('Crossovers', 'Mutations', 'Location', 'best');
    xlabel('Generation');
    ylabel('Operator Count');
    title('Genetic Operator Activity');
    grid on;
    
    % Spatial Distribution History
    subplot(2,2,3);
    hold on;
    colors = lines(numBS);
    for bs = 1:numBS
        colBase = (bs-1)*5;
        x = history.bestIndividuals(1:gen, colBase + 1);
        y = history.bestIndividuals(1:gen, colBase + 2);
        plot(x, y, 'o-', 'Color', colors(bs,:), 'DisplayName', sprintf('BS%d', bs));
    end
    xlabel('X Position (m)');
    ylabel('Y Position (m)');
    title('Base Station Position Evolution');
    legend('show');
    grid on;


    % Altitude Distribution
    subplot(2,2,4);
    altitudeData = zeros(gen, numBS);
    for bs = 1:numBS
        colBase = (bs-1)*5;
        altitudeData(:, bs) = history.bestIndividuals(1:gen, colBase + 3);
    end
    boxplot(altitudeData, ...
        'Labels', arrayfun(@(n)sprintf('BS%d',n), 1:numBS, 'UniformOutput', false));
    ylabel('Altitude (m)');
    title('Base Station Altitude Distribution');
    grid on;

    
    % trajectory plot 
    step = 5; 
    num_pairs = params.numBS;
    
    figure; 
    hold on; axis equal; grid on;
    title('Trajectories of Multiple XY Pairs');
    xlabel('X'); ylabel('Y');
    
    colors = lines(num_pairs);
    trajectory_handles = gobjects(num_pairs, 1);
    scatter(mbs_y,mbs_x, 'marker','x','MarkerEdgeColor','black', 'MarkerFaceColor','black')
    
    for k = 1:num_pairs
        col_x = (k-1)*step + 1;
        col_y = (k-1)*step + 2;
    
        x = history.bestIndividuals(:, col_x);
        y = history.bestIndividuals(:, col_y);
        xy = [x y];
    
        [unique_xy, ia, ~] = unique(xy, 'rows', 'stable');
        x_unique = unique_xy(:,1);
        y_unique = unique_xy(:,2);
    
        % Draw arrows with smaller arrowheads
        for i = 1:(size(unique_xy,1)-1)
            quiver(x_unique(i), y_unique(i), ...
                   x_unique(i+1) - x_unique(i), ...
                   y_unique(i+1) - y_unique(i), ...
                   0, 'Color', colors(k,:), 'LineWidth', 1.5, 'MaxHeadSize', 0.2);
        end
    
        % Scatter for legend handle
        trajectory_handles(k) = scatter(x_unique, y_unique, 50, colors(k,:), 'filled');
    
        % Offset for start text
        offset_x = 30; % Adjust as needed for your data scale
        offset_y = 30; % Adjust as needed for your data scale
    
        text(x_unique(1) + offset_x, y_unique(1) + offset_y, sprintf('Start %d', k), ...
            'FontSize', 15, 'FontWeight', 'bold', 'Color', colors(k,:), ...
            'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
    
        % Offset and size for point number labels
        for i = 1:length(x_unique)
            text(x_unique(i) + 10, y_unique(i) + 10, sprintf('%d', i), ...
                'FontSize', 13, 'Color', colors(k,:), ...
                'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
        end
    
        % Mark start and end points
        plot(x_unique(1), y_unique(1), 's', 'MarkerSize', 12, 'LineWidth', 2, 'Color', colors(k,:));
        plot(x_unique(end), y_unique(end), '*', 'MarkerSize', 12, 'LineWidth', 2, 'Color', colors(k,:));
    end
    
    % lgd = legend(trajectory_handles, {'FBS 1', 'FBS 2'}, 'Location', 'northwest');
    % set(lgd, 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial');



    % power vs height plot
    figure;
    heights = history.bestIndividuals(1:gen, 3:5:end);
    powers  = history.bestIndividuals(1:gen, 4:5:end);
    generations = 1:gen;
    
    colors = lines(numBS);
    yyaxis left
    hold on;
    for bs = 1:numBS
        plot(generations, heights(:, bs), '-o', 'Color', colors(bs,:), 'DisplayName', sprintf('Height BS%d', bs));
    end
    ylabel('Height (m)');
    minHeight = min(heights(:));
    maxHeight = max(heights(:));
    ylim([minHeight-5, maxHeight+5]);
    
    yyaxis right
    for bs = 1:numBS
        plot(generations, powers(:, bs), '--', 'Color', colors(bs,:), 'DisplayName', sprintf('Power BS%d', bs));
    end
    ylabel('Transmitted Power (W)');
    minPower = min(powers(:));
    maxPower = max(powers(:));
    ylim([minPower-0.5, maxPower+0.5]);
    
    xlabel('Generation');
    title('Progression of Height and Power over Generations');
    legend('Location', 'bestoutside');
    grid on;

    end

end

if params.enableLogging
    logEntry = struct( ...
        'beta', params.fitnessWeights.beta, ...
        'gamma', params.fitnessWeights.gamma, ...
        'fbsWeight', params.fitnessWeights.fbsWeight, ...
        'fbsExponent', params.fitnessWeights.fbsExponent, ...
        'mutationProb', params.mutationProb, ...
        'crossoverProb', params.crossoverProb, ...
        'numFbs', params.numBS, ...
        'finalXYZ', formatFinalLocations(bestIndividual, params.numBS), ...
        'powerStatus', formatPowerStatuses(bestIndividual, params.numBS), ...
        'totalPower', transmittedPower, ...
        'fbsPowers', formatFbsPowers(bestIndividual, params.numBS), ...
        'fbsConnected', fbsUsers, ...
        'mbsConnected', mbsUsers, ...
        'totalConnected', numUsers, ...
        'avgRate', avg_rate_connected_bpsHz);
    appendRunLog(params.logFile, logEntry);
end

end
function updateTrajectoryPlot(plotStruct, population, generation, bounds)
    ax = plotStruct.ax;
    if isempty(population)
        cla(ax);
        return;
    end

    cla(ax);
    hold(ax, 'on');
    x = population(:,1);
    y = population(:,2);
    z = population(:,3);
    status = population(:,5) >= 0.5;
    colorValues = generation * ones(size(x));

    if any(status)
        scatter3(ax, x(status), y(status), z(status), 60, colorValues(status), 'filled', ...
            'Marker', 'o', 'DisplayName', 'Active FBS');
    end
    if any(~status)
        scatter3(ax, x(~status), y(~status), z(~status), 60, colorValues(~status), ...
            'Marker', 'x', 'LineWidth', 1.5, 'DisplayName', 'Inactive FBS');
    end

    title(ax, sprintf('Population Trajectory (Gen %d)', generation));
    legend(ax, 'Location', 'best');
    axis(ax, 'tight');
    grid(ax, 'on');
    hold(ax, 'off');

    if size(bounds,1) >= 3
        xlim(ax, bounds(1,:));
        ylim(ax, bounds(2,:));
        zlim(ax, bounds(3,:));
    end
end

function finalXYZStr = formatFinalLocations(individual, numBS)
    triplets = strings(numBS, 1);
    for bs = 1:numBS
        idx = (bs-1)*5 + 1;
        coords = individual(idx:idx+2);
        triplets(bs) = sprintf('[%g,%g,%g]', coords);
    end
    finalXYZStr = strjoin(triplets, '; ');
end

function statusStr = formatPowerStatuses(individual, numBS)
    statuses = strings(numBS, 1);
    for bs = 1:numBS
        idx = (bs-1)*5 + 5;
        statuses(bs) = string(individual(idx));
    end
    statusStr = strjoin(statuses, '; ');
end

function powersStr = formatFbsPowers(individual, numBS)
    powers = strings(numBS, 1);
    for bs = 1:numBS
        idx = (bs-1)*5 + 4;
        powers(bs) = sprintf('%g', individual(idx));
    end
    powersStr = strjoin(powers, '; ');
end

function appendRunLog(filename, logEntry)
    logTable = struct2table(logEntry);
    if isfile(filename)
        writetable(logTable, filename, 'WriteMode', 'append', 'WriteVariableNames', false);
    else
        writetable(logTable, filename);
    end
end
