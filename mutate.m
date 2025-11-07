function [mutated, mutationFlags] = mutate(individual, prob, bounds, mutationScale, gen, params, adaptiveParams)
    
    n_fbs = length(individual)/5;
    mutationFlags = false(size(individual));
    
    for bs = 1:n_fbs
        if rand() < prob
            idx = (1:5) + (bs-1)*5;
            bin_idx = bs * 5;
            
            paramRanges = bounds(idx,2) - bounds(idx,1);
            
            if exist('adaptiveParams','var') && ~isempty(adaptiveParams)
                sigmaLearn = zeros(4,1);
                muLearn    = zeros(4,1);
            
                % Loop over each continuous variable (X, Y, Z, Power)
                for varIdx = 1:4
                    if params.rl.enableSigma(varIdx)
                        sigmaLearn(varIdx) = adaptiveParams(varIdx) * paramRanges(varIdx);
                    else
                        sigmaLearn(varIdx) = mutationScale * paramRanges(varIdx);
                    end
            
                    if params.rl.enableMean(varIdx)
                        muLearn(varIdx) = adaptiveParams(4 + varIdx) * paramRanges(varIdx);
                    else
                        muLearn(varIdx) = 0;
                    end
                end
            else
                % Fallback to static mutation for all variables
                sigmaLearn = mutationScale * paramRanges(1:4);
                muLearn    = zeros(4,1);
            end
            
            perturbation = muLearn(:) + sigmaLearn(:) .* randn(4,1);



            mutatedBS = individual(idx(1:end-1)) + perturbation';
            mutatedBS = reflectToBounds(mutatedBS, bounds(idx(1:end-1),:));
            
            individual(idx(1:end-1)) = mutatedBS; % continuous variable mutations 
            individual(bin_idx) = 1 - individual(bin_idx); % binary variable mutations/flipping 

            mutationFlags(idx) = true;
        end

    end

    mutated = individual;
end

%% prior to change with crossover and mutation probability inclusion in q learning

% function [mutated, mutationFlags] = mutate(individual, prob, bounds, mutationScale, gen, params, adaptiveParams)
% 
%     n_fbs = length(individual)/5;
%     mutationFlags = false(size(individual));
% 
%     for bs = 1:n_fbs
%         if rand() < prob
%             idx = (1:5) + (bs-1)*5;
%             bin_idx = bs * 5;
% 
%             paramRanges = bounds(idx,2) - bounds(idx,1);
% 
%             % optional diminishing sigma based on generation index
%             % sigma = (params.numGenerations - gen)/(params.numGenerations) * mutationScale * paramRanges;
% 
%             if isempty(adaptiveParams)
%                 % normal sigma
%                 sigma = mutationScale * paramRanges;
%             else
%                 % adaptive sigma from RL output
%                 sigma = adaptiveParams(1,1) * paramRanges;
%                 mean_perturbations =  adaptiveParams(1,2)*paramRanges;
%                 mean_perturbations = mean_perturbations(1:4);
%             end
% 
%             if isempty(adaptiveParams)
%                 perturbation = sigma(1:4) .* randn(4,1);
% 
%             else
%                 perturbation = mean_perturbations + sigma(1:4) .* randn(4,1);
%             end
% 
%             mutatedBS = individual(idx(1:end-1)) + perturbation';
%             mutatedBS = reflectToBounds(mutatedBS, bounds(idx(1:end-1),:));
% 
%             individual(idx(1:end-1)) = mutatedBS; % continuous variable mutations 
%             individual(bin_idx) = 1 - individual(bin_idx); % binary variable mutations/flipping 
% 
%             mutationFlags(idx) = true;
%         end
% 
%     end
% 
%     mutated = individual;
% end
