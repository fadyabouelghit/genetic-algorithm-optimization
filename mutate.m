function [mutated, mutationFlags] = mutate(individual, prob, bounds, mutationScale, gen, params)
    
    n_fbs = length(individual)/5;
    mutationFlags = false(size(individual));
    
    for bs = 1:n_fbs
        if rand() < prob
            idx = (1:5) + (bs-1)*5;
            bin_idx = bs * 5;
            
            paramRanges = bounds(idx,2) - bounds(idx,1);
            
            % optional diminishing sigma based on generation index
            % sigma = (params.numGenerations - gen)/(params.numGenerations) * mutationScale * paramRanges;
            sigma = mutationScale * paramRanges;
            
            perturbation = sigma(1:4) .* randn(4,1);
            
            mutatedBS = individual(idx(1:end-1)) + perturbation';
            mutatedBS = reflectToBounds(mutatedBS, bounds(idx(1:end-1),:));
            
            individual(idx(1:end-1)) = mutatedBS; % continuous variable mutations 
            individual(bin_idx) = 1 - individual(bin_idx); % binary variable mutations/flipping 

            mutationFlags(idx) = true;
        end

    end

    mutated = individual;
end
