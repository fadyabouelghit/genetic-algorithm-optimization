function [child1, child2, crossoverFlag] = crossover_blend(parent1, parent2, prob, bounds)
    
    crossoverFlag = 0;
    alpha = 0.5;
    n_fbs = length(parent1) / 5;  % Each BS has 5 parameters

    if rand() < prob
        crossoverFlag = 1;
        child1 = parent1;
        child2 = parent2;

        % Loop through each base station block
        for bs = 1:n_fbs
            % Continuous parameters (X, Y, Z, Power)
            idx = (bs-1)*5 + 1 : (bs-1)*5 + 4;
            % Binary parameter (Power Status)
            bin_idx = (bs-1)*5 + 5;

            % Blend crossover for continuous parameters
            gamma = (1 + 2*alpha) * rand(1,4) - alpha;
            child1(idx) = (1 - gamma) .* parent1(idx) + gamma .* parent2(idx);
            child2(idx) = gamma .* parent1(idx) + (1 - gamma) .* parent2(idx);
            
            if parent1(bin_idx) == parent2(bin_idx) % XOR crossover with random replacement
                child1(bin_idx) = parent1(bin_idx);
                child2(bin_idx) = parent1(bin_idx);
            else
                child1(bin_idx) = randi([0, 1]);
                child2(bin_idx) = 1 - child1(bin_idx);  % Opposite value
            end

        end

        % Clamp continuous parameters to bounds
        child1(1:end-1) = clampToBounds(child1(1:end-1), bounds(1:end-1,:));
        child2(1:end-1) = clampToBounds(child2(1:end-1), bounds(1:end-1,:));
    else
        child1 = parent1;
        child2 = parent2;
    end
end
