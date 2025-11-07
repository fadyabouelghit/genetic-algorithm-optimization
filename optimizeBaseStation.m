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
params = mergeParams(defaultParams, params);
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


% population = initializePopulation_sobol(params.populationSize, params.bounds);
% population = initializePopulation_ppp(params.populationSize, params.bounds, params.numBS);
population = initializePopulation_uniform(params.populationSize, params.bounds, params.numBS);


headers = create_headers(params.numBS);
initialData = [zeros(params.populationSize,1) population]; % Gen 0
% writetable(array2table(initialData, 'VariableNames', headers), filename, 'WriteMode', 'overwrite');

mbs_x = mbs_params(1,:);
mbs_y = mbs_params(2,:);
mbs_height = mbs_params(3,:);
mbs_power = mbs_params(4,:);


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
    [fitness, evalDetails] = evaluatePopulation(l, population, params.verbose, params.numBS, params.spaceLimit ,containsMbs, mbs_params, antennaObjectMbs, params.bounds, params.mbsCache, targetIdx);
    
    % Track population statistics
    history.bestFitness(gen) = max(fitness);
    history.maxFitness(gen) = max(fitness);
    history.minFitness(gen) = min(fitness);
    history.avgFitness(gen) = mean(fitness);
    history.stdFitness(gen) = std(fitness);

    [bestFitness, bestIdx] = max(fitness);
    bestIndividual = population(bestIdx, :);
    history.bestIndividuals(gen, :) = bestIndividual;
    
    history.avgVariables(gen, :) = mean(population, 1);
    
    
    % Display generation summary
    if params.verbose > 0
        fprintf('Fitness: Best=%.2f, Avg=%.2f, Std=%.2f\n', ...
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
    newPopulation = [bestIndividual; zeros(params.populationSize-1, size(params.bounds, 1))];
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
    population
    
    % Early stopping check
    % if gen > 8 && std(history.bestFitness(gen-7:gen)) < 0.01
    %     if params.verbose > 0
    %         fprintf('\nEarly stopping: No improvement in last 8 generations\n');
    %     end
    %     break;
    % end
end

history.time.total = toc(totalTimer);
[bestFitness, bestIdx] = max(fitness);
bestIndividual = population(bestIdx, :);

if params.verbose > 0
    fprintf('\n=== Optimization Complete ===\n');
    fprintf('Total time: %.1f minutes\n', history.time.total/60);
    fprintf('Best Fitness: %.2f users\n', bestFitness);
    
    fprintf('The best individual performance: \n') 
    [~, ~, numUsers, transmittedPower, avg_rate_connected_bpsHz] = SINREvaluation(...,
                                                        l, bestIndividual(5:5:end), bestIndividual(1:5:end), bestIndividual(2:5:end), bestIndividual(3:5:end), params.numBS, ...,
                                                        bestIndividual(4:5:end), mbs_y, mbs_x, mbs_height, mbs_power, ...,
                                                        0, params.spaceLimit(1), 0, params.spaceLimit(2), 1000, 5, containsMbs, antennaObjectMbs, params.mbsCache);
    fprintf('Total Connected Users: %d\n', numUsers);
    fprintf('Total Transmitted Power: %.2f W\n', transmittedPower);
    fprintf('Avg. Sum Rate: %.2f W\n', avg_rate_connected_bpsHz);

    
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
        x = history.bestIndividuals(1:gen, (bs-1)*4 + 1);
        y = history.bestIndividuals(1:gen, (bs-1)*4 + 2);
        plot(x, y, 'o-', 'Color', colors(bs,:), 'DisplayName', sprintf('BS%d', bs));
    end
    xlabel('X Position (m)');
    ylabel('Y Position (m)');
    title('Base Station Position Evolution');
    legend('show');
    grid on;


    % Altitude Distribution
    subplot(2,2,4);
    boxplot(reshape(history.bestIndividuals(1:gen, 3:4:end), [], numBS), ...
        'Labels', arrayfun(@(n)sprintf('BS%d',n), 1:numBS, 'UniformOutput', false));
    ylabel('Altitude (m)');
    title('Base Station Altitude Distribution');
    grid on;

    
    % trajectory plot 
    step = 5; 
    num_pairs = floor(size(history.bestIndividuals, 2) / step);
    
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
    % Extract height and power from history
    heights = history.bestIndividuals(:, 3);  % Assuming 3rd column = height
    powers  = history.bestIndividuals(:, 4);  % Assuming 4th column = power
    generations = 1:length(heights);          % X-axis: generation index
    
    yyaxis left
    plot(generations, heights, '-o', 'LineWidth', 2);
    ylabel('Height (m)');
    ylim([min(heights)-5, max(heights)+5]);  % Optional padding
    
    yyaxis right
    plot(generations, powers, '-s', 'LineWidth', 2);
    ylabel('Transmitted Power (W)');
    ylim([min(powers)-0.5, max(powers)+0.5]);  % Optional padding
    
    xlabel('Generation');
    title('Progression of Height and Power over Generations');
    legend('Height', 'Power');
    grid on;

end

end
