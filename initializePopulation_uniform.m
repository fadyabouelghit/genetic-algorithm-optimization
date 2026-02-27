function population = initializePopulation_uniform(popSize, bounds, n_fbs)

rng(43);
numParams = size(bounds,1);
population = zeros(popSize, numParams);
numFbsParams = 5 * n_fbs;
hasFreqFlag = (numParams == numFbsParams + 1);

for i = 1:numParams
    
    % if i == 1 || i == 2
    %     lb = 0;
    %     ub = 300;
    % elseif i == 3
    %     lb = 20;
    %     ub = 50;
    % else
        % lb = bounds(i,1);
        % ub = bounds(i,2);
    % end

    lb = bounds(i,1);
    ub = bounds(i,2);

    if hasFreqFlag && i == numFbsParams + 1
        % Global binary frequency selector gene: fbsFreqFlag
        population(:,i) = randi([0 1], popSize, 1);
    elseif i <= numFbsParams && mod(i,5) == 0  % 5th parameter in each FBS block is power_status
        % Binary sampling: 0 or 1
        population(:,i) = randi([0 1], popSize, 1);
    else
        % Uniform sampling in [lb, ub]
        population(:,i) = lb + (ub - lb) * rand(popSize, 1);
    end
end

end
