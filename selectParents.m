function parents = selectParents(population, fitness, nParents)
% SELECTPARENTS Performs tournament selection of parents for crossover

    tournamentSize = 3;
    parents = zeros(nParents, size(population,2));
    for i = 1:nParents
        candidates = randperm(size(population,1), tournamentSize);
        [~, idx] = max(fitness(candidates));
        parents(i,:) = population(candidates(idx),:);
    end
end
