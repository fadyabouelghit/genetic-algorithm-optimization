function population = initializePopulation_uniform(popSize, bounds, n_fbs)

rng(43);
numParams = size(bounds,1);
population = zeros(popSize, numParams);

for i = 1:numParams
    lb = bounds(i,1);
    ub = bounds(i,2);

    if mod(i,5) == 0  % 5th parameter in each FBS block is power_status
        % Binary sampling: 0 or 1
        population(:,i) = randi([0 1], popSize, 1);
    else
        % Uniform sampling in [lb, ub]
        population(:,i) = lb + (ub - lb) * rand(popSize, 1);
    end
end

end
