function [paretoFront, history] = optimizeBaseStationMoga(l, containsMbs, antennaObjectMbs, mbs_params, params)
% OPTIMIZEBASESTATIONMOGA Multi-objective GA for connectivity vs power
%   Runs an NSGA-II style evolutionary loop to maximize connected users
%   while minimizing transmitted power. The resulting Pareto fronts per
%   generation are stored in history and the final front is returned.

    defaultParams = struct(...
        'verbose', 1, ...
        'tournamentSize', 3, ...
        'fitnessWeights', struct('beta', 1, 'gamma', 1, 'epsilon', 1e-3), ...
        'initialPopulationSize', [], ...
        'maxUsers', 1000, ...
        'sinrThreshold', 5 ...
    );
    params = mergeParams(defaultParams, params);
    if isempty(params.initialPopulationSize)
        params.initialPopulationSize = params.populationSize;
    end
    if params.initialPopulationSize < params.populationSize
        error('initialPopulationSize (%d) must be >= populationSize (%d).', ...
            params.initialPopulationSize, params.populationSize);
    end

    populationSize = params.populationSize;
    numParams = size(params.bounds, 1);
    targetIdx = 1; % for population evaluation: 1 - connectivity & power, 2 - avg rate

    population = initializePopulation_uniform(params.initialPopulationSize, params.bounds, params.numBS);

    evalParams = params.fitnessWeights;
    evalParams.maxUsers = params.maxUsers;
    evalParams.sinrThreshold = params.sinrThreshold;

    [~, details] = evaluatePopulation(l, population, params.verbose, params.numBS, ...
        params.spaceLimit, containsMbs, mbs_params, antennaObjectMbs, params.bounds, params.mbsCache, targetIdx, evalParams);
    objectives = buildObjectives(details);

    [fronts, ranks] = fastNonDominatedSort(objectives);
    crowding = computeCrowdingDistances(objectives, fronts);

    if params.initialPopulationSize > populationSize
        survivorIdx = environmentalSelection(fronts, crowding, populationSize);
        population = population(survivorIdx, :);
        details = sliceDetails(details, survivorIdx);
        objectives = objectives(survivorIdx, :);
        [fronts, ranks] = fastNonDominatedSort(objectives);
        crowding = computeCrowdingDistances(objectives, fronts);
    end

    history = initHistory(params.numGenerations);

    for gen = 1:params.numGenerations
        if params.verbose > 0
            fprintf('\n[MOGA] Generation %d/%d\n', gen, params.numGenerations);
        end

        [offspring, crossovers, mutations] = createOffspring(population, params, ranks, crowding, gen);

        evalParams = params.fitnessWeights;
        evalParams.maxUsers = params.maxUsers;
        evalParams.sinrThreshold = params.sinrThreshold;
        [~, offspringDetails] = evaluatePopulation(l, offspring, params.verbose, params.numBS, ...
            params.spaceLimit, containsMbs, mbs_params, antennaObjectMbs, params.bounds, params.mbsCache, targetIdx, evalParams);
        offspringObjectives = buildObjectives(offspringDetails);

        combinedPopulation = [population; offspring];
        combinedDetails = stackDetails(details, offspringDetails);
        combinedObjectives = [objectives; offspringObjectives];

        [combinedFronts, combinedRanks] = fastNonDominatedSort(combinedObjectives);
        combinedCrowding = computeCrowdingDistances(combinedObjectives, combinedFronts);

        survivorIdx = environmentalSelection(combinedFronts, combinedCrowding, populationSize);
        population = combinedPopulation(survivorIdx, :);
        details = sliceDetails(combinedDetails, survivorIdx);
        objectives = combinedObjectives(survivorIdx, :);

        [fronts, ranks] = fastNonDominatedSort(objectives);
        crowding = computeCrowdingDistances(objectives, fronts);

        history.crossovers(gen) = crossovers;
        history.mutations(gen) = mutations;
        history.fronts{gen} = captureFront(population, details, fronts);
        history.populationObjectives{gen} = [details.numUsers, details.transmittedPower, details.activeFbs];
        history.ranks{gen} = ranks;

        if params.verbose > 0
            users = details.numUsers;
            power = details.transmittedPower;
            fprintf('  Population stats | users(min/mean/max)=%.0f/%.1f/%.0f | power(min/mean/max)=%.2f/%.2f/%.2f W\n', ...
                min(users), mean(users), max(users), min(power), mean(power), max(power));
            fprintf('  Operators: %d crossovers, %d mutations\n', crossovers, mutations);

            currentFront = history.fronts{gen};
            if isempty(currentFront.users)
                fprintf('  No feasible Pareto points captured this generation.\n');
            else
                [~, bestIdx] = max(currentFront.users - currentFront.power);
                fprintf('  Pareto front size: %d | best trade-off users=%.0f, power=%.2f W, active FBS=%d\n', ...
                    numel(currentFront.users), currentFront.users(bestIdx), currentFront.power(bestIdx), currentFront.activeFbs(bestIdx));
            end
        end
    end

    paretoFront = history.fronts{params.numGenerations};

    if params.verbose > 0
        plotParetoHistory(history);
    end
end

function history = initHistory(numGenerations)
    history = struct(...
        'fronts', {cell(numGenerations, 1)}, ...
        'populationObjectives', {cell(numGenerations, 1)}, ...
        'ranks', {cell(numGenerations, 1)}, ...
        'crossovers', zeros(numGenerations, 1), ...
        'mutations', zeros(numGenerations, 1) ...
    );
end

function objectives = buildObjectives(details)
    objectives = [-details.numUsers, details.transmittedPower, details.activeFbs];
end

function [offspring, crossoverCount, mutationCount] = createOffspring(population, params, ranks, crowding, gen)
    popSize = params.populationSize;
    numParams = size(population, 2);
    offspring = zeros(popSize, numParams);
    offspringIdx = 1;
    crossoverCount = 0;
    mutationCount = 0;
    tournamentSize = max(2, params.tournamentSize);

    while offspringIdx <= popSize
        p1Idx = binaryTournamentSelection(ranks, crowding, tournamentSize);
        p2Idx = binaryTournamentSelection(ranks, crowding, tournamentSize);

        parent1 = population(p1Idx, :);
        parent2 = population(p2Idx, :);

        [child1, child2, crossoverFlag] = crossover_blend(parent1, parent2, params.crossoverProb, params.bounds);
        crossoverCount = crossoverCount + crossoverFlag;

        [child1, mutationFlags1] = mutate(child1, params.mutationProb, params.bounds, params.mutationScale, gen, params, []);
        [child2, mutationFlags2] = mutate(child2, params.mutationProb, params.bounds, params.mutationScale, gen, params, []);
        mutationCount = mutationCount + sum(mutationFlags1) + sum(mutationFlags2);

        offspring(offspringIdx, :) = child1;
        if offspringIdx + 1 <= popSize
            offspring(offspringIdx + 1, :) = child2;
        end
        offspringIdx = offspringIdx + 2;
    end
end

function idx = binaryTournamentSelection(ranks, crowding, tournamentSize)
    popSize = numel(ranks);
    candidates = randi(popSize, tournamentSize, 1);
    best = candidates(1);
    for i = 2:numel(candidates)
        contender = candidates(i);
        if ranks(contender) < ranks(best)
            best = contender;
        elseif ranks(contender) == ranks(best)
            if crowding(contender) > crowding(best)
                best = contender;
            elseif crowding(contender) == crowding(best) && rand() < 0.5
                best = contender;
            end
        end
    end
    idx = best;
end

function [fronts, ranks] = fastNonDominatedSort(objectives)
    popSize = size(objectives, 1);
    dominationCounts = zeros(popSize, 1);
    dominatedSets = cell(popSize, 1);
    ranks = inf(popSize, 1);
    fronts = {};
    firstFront = [];

    for p = 1:popSize
        dominatedSets{p} = [];
        for q = 1:popSize
            if p == q
                continue;
            end
            if dominates(objectives(p,:), objectives(q,:))
                dominatedSets{p}(end+1) = q; %#ok<AGROW>
            elseif dominates(objectives(q,:), objectives(p,:))
                dominationCounts(p) = dominationCounts(p) + 1;
            end
        end
        if dominationCounts(p) == 0
            ranks(p) = 1;
            firstFront(end+1) = p; %#ok<AGROW>
        end
    end

    fronts{1} = firstFront;
    i = 1;
    while ~isempty(fronts{i})
        nextFront = [];
        for p = fronts{i}
            for q = dominatedSets{p}
                dominationCounts(q) = dominationCounts(q) - 1;
                if dominationCounts(q) == 0
                    ranks(q) = i + 1;
                    nextFront(end+1) = q; %#ok<AGROW>
                end
            end
        end
        i = i + 1;
        fronts{i} = nextFront; %#ok<AGROW>
    end

    if isempty(fronts{end})
        fronts(end) = [];
    end
end

function flag = dominates(a, b)
    flag = all(a <= b) && any(a < b);
end

function distances = computeCrowdingDistances(objectives, fronts)
    popSize = size(objectives, 1);
    numObjectives = size(objectives, 2);
    distances = zeros(popSize, 1);

    for f = 1:numel(fronts)
        front = fronts{f};
        if numel(front) == 0
            continue;
        elseif numel(front) <= 2
            distances(front) = inf;
            continue;
        end

        frontObjectives = objectives(front, :);
        frontDistances = zeros(numel(front), 1);

        for m = 1:numObjectives
            [sortedVals, sortedIdx] = sort(frontObjectives(:, m));
            sortedFront = front(sortedIdx);
            frontDistances(sortedIdx(1)) = inf;
            frontDistances(sortedIdx(end)) = inf;

            objectiveRange = sortedVals(end) - sortedVals(1);
            if objectiveRange == 0
                continue;
            end

            for k = 2:numel(sortedFront)-1
                prevVal = sortedVals(k-1);
                nextVal = sortedVals(k+1);
                frontDistances(sortedIdx(k)) = frontDistances(sortedIdx(k)) + ...
                    (nextVal - prevVal) / objectiveRange;
            end
        end

        distances(front) = frontDistances;
    end
end

function survivorIdx = environmentalSelection(fronts, crowding, populationSize)
    survivorIdx = [];
    for f = 1:numel(fronts)
        front = fronts{f};
        if isempty(front)
            continue;
        end
        remainingSlots = populationSize - numel(survivorIdx);
        if numel(front) <= remainingSlots
            survivorIdx = [survivorIdx; front(:)]; %#ok<AGROW>
        else
            [~, sortIdx] = sort(crowding(front), 'descend');
            chosen = front(sortIdx(1:remainingSlots));
            survivorIdx = [survivorIdx; chosen(:)]; %#ok<AGROW>
            break;
        end
    end
end

function stacked = stackDetails(a, b)
    fields = fieldnames(a);
    stacked = struct();
    for i = 1:numel(fields)
        stacked.(fields{i}) = [a.(fields{i}); b.(fields{i})];
    end
end

function sliced = sliceDetails(details, idx)
    fields = fieldnames(details);
    sliced = struct();
    for i = 1:numel(fields)
        values = details.(fields{i});
        sliced.(fields{i}) = values(idx, :);
    end
end

function frontStruct = captureFront(population, details, fronts)
    if isempty(fronts) || isempty(fronts{1})
        frontStruct = struct('individuals', [], 'users', [], 'power', [], 'avgRate', [], 'activeFbs', []);
        return;
    end
    idx = fronts{1};
    frontStruct = struct(...
        'individuals', population(idx, :), ...
        'users', details.numUsers(idx), ...
        'power', details.transmittedPower(idx), ...
        'avgRate', details.avgRate(idx), ...
        'activeFbs', details.activeFbs(idx));
end

function plotParetoHistory(history)
    numGenerations = numel(history.fronts);
    colors = parula(max(numGenerations, 2));
    figure('Name', 'MOGA Pareto Front', 'Position', [100 100 900 600]);
    hold on; grid on;

    for gen = 1:numGenerations
        front = history.fronts{gen};
        if isempty(front.users)
            continue;
        end
        scatter3(front.power, front.users, front.activeFbs, 50, colors(gen, :), 'filled', ...
            'DisplayName', sprintf('Gen %d', gen));
    end

    finalFront = history.fronts{end};
    if ~isempty(finalFront.users)
        [sortedPower, order] = sortrows([finalFront.power(:), finalFront.users(:), finalFront.activeFbs(:)], [1 3]);
        plot3(sortedPower(:,1), sortedPower(:,2), sortedPower(:,3), '-k', 'LineWidth', 1.5, 'DisplayName', 'Final Front');
    end

    xlabel('Transmitted Power (W)');
    ylabel('Connected Users');
    zlabel('Active FBS Count');
    title('Pareto Front Evolution (Power vs Connectivity vs Active FBS)');
    legend('Location', 'best');
    view(35, 25);
end
