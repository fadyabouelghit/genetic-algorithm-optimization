function [child1, child2, crossoverFlag] = crossover_blend(parent1, parent2, prob, bounds)
    
    crossoverFlag = 0;
    alpha = 0.5;
    blockSize = 6;
    if mod(length(parent1), blockSize) ~= 0
        error('crossover_blend expects chromosome size to be a multiple of 6.');
    end
    n_fbs = floor(length(parent1) / blockSize);  % Each BS has 6 parameters

    if rand() < prob
        crossoverFlag = 1;
        child1 = parent1;
        child2 = parent2;

        % Loop through each base station block
        for bs = 1:n_fbs
            % Continuous parameters (X, Y, Z, Power)
            idx = (bs-1)*blockSize + 1 : (bs-1)*blockSize + 4;
            % Binary parameters
            status_idx = (bs-1)*blockSize + 5;
            freq_idx = (bs-1)*blockSize + 6;

            % Blend crossover for continuous parameters
            gamma = (1 + 2*alpha) * rand(1,4) - alpha;
            child1(idx) = (1 - gamma) .* parent1(idx) + gamma .* parent2(idx);
            child2(idx) = gamma .* parent1(idx) + (1 - gamma) .* parent2(idx);
            
            if parent1(status_idx) == parent2(status_idx)
                child1(status_idx) = parent1(status_idx);
                child2(status_idx) = parent1(status_idx);
            else
                child1(status_idx) = randi([0, 1]);
                child2(status_idx) = 1 - child1(status_idx);
            end

            if parent1(freq_idx) == parent2(freq_idx)
                child1(freq_idx) = parent1(freq_idx);
                child2(freq_idx) = parent1(freq_idx);
            else
                child1(freq_idx) = randi([0, 1]);
                child2(freq_idx) = 1 - child1(freq_idx);
            end
        end

        for bs = 1:n_fbs
            idx = (bs-1)*blockSize + 1 : (bs-1)*blockSize + 4;
            child1(idx) = clampToBounds(child1(idx), bounds(idx,:));
            child2(idx) = clampToBounds(child2(idx), bounds(idx,:));
        end
    else
        child1 = parent1;
        child2 = parent2;
    end
end
