function population = initializePopulation_ppp(popSize, bounds, n_fbs)
    
    rng(42);
    % Parameters
    minDistance = 300; % Minimum distance between BSs in meters
    maxAttempts = 150; % Max attempts per BS placement
    areaSize = [bounds(1,2)-bounds(1,1), bounds(2,2)-bounds(2,1)]; % [width, height]

    population = zeros(popSize, n_fbs*5);

    for i = 1:popSize
        [x, y] = maternPP(areaSize, n_fbs, minDistance, maxAttempts);
        z = bounds(3,1) + (bounds(3,2)-bounds(3,1))*rand(n_fbs,1);
        power = bounds(4,1) + (bounds(4,2)-bounds(4,1))*rand(n_fbs,1);

        % assign power status bernoulli p=0.5 
        power_status = randi([0, 1], 1, n_fbs);


        % Store in individual vector
        individual = zeros(1, n_fbs*5);
        for bs = 1:n_fbs
            idx = (bs-1)*5 + 1;
            individual(idx:idx+4) = [x(bs), y(bs), z(bs), power(bs), power_status(bs)];
        end

        population(i,:) = individual;
    end

    % Shift to actual bounds
    population(:,1:4:end) = population(:,1:4:end) + bounds(1,1);
    population(:,2:4:end) = population(:,2:4:end) + bounds(2,1);
end
