function population = initializePopulation_uniform(popSize, bounds, n_fbs)

rng(43);
numParams = size(bounds,1);
population = zeros(popSize, numParams);
blockSize = 6;
expectedParams = blockSize * n_fbs;
if numParams ~= expectedParams
    error('initializePopulation_uniform expects %d params (6 per FBS), got %d.', expectedParams, numParams);
end

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

    posInBlock = mod(i-1, blockSize) + 1;
    if posInBlock == 5 || posInBlock == 6
        % Binary sampling: 0 or 1
        population(:,i) = randi([0 1], popSize, 1);
    else
        % Uniform sampling in [lb, ub]
        population(:,i) = lb + (ub - lb) * rand(popSize, 1);
    end
end

end
