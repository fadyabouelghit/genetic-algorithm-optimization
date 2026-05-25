function [mutated, mutationFlags] = mutate(individual, prob, bounds, mutationScale, gen, params, adaptiveParams)

    blockSize = 6;
    if isfield(params, 'numBS')
        n_fbs = params.numBS;
    else
        n_fbs = floor(length(individual)/blockSize);
    end
    fbsCount = blockSize * n_fbs;
    n_mbs = length(individual) - fbsCount;
    if n_mbs < 0
        error('mutate: chromosome shorter than 6*n_fbs (got %d, n_fbs=%d).', length(individual), n_fbs);
    end
    if isfield(params, 'gaControlsFbsBand'),     gaControlsFbsBand     = params.gaControlsFbsBand;     else, gaControlsFbsBand     = true; end
    if isfield(params, 'gaControlsMbsCapacity'), gaControlsMbsCapacity = params.gaControlsMbsCapacity; else, gaControlsMbsCapacity = true; end
    mutationFlags = false(size(individual));

    for bs = 1:n_fbs
        if rand() < prob
            idx = (1:blockSize) + (bs-1)*blockSize;
            status_idx = (bs-1)*blockSize + 5;
            freq_idx = (bs-1)*blockSize + 6;

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



            cont_idx = idx(1:4);
            mutatedBS = individual(cont_idx) + perturbation';
            mutatedBS = reflectToBounds(mutatedBS, bounds(cont_idx,:));

            individual(cont_idx) = mutatedBS; % continuous variable mutations
            individual(status_idx) = 1 - individual(status_idx); % power_status flip
            if gaControlsFbsBand
                individual(freq_idx) = 1 - individual(freq_idx); % fbsFreqFlag flip
            end

            mutationFlags(idx) = true;
        end

    end

    % MBS capacity-flag suffix: each gene is its own atomic block.
    % For base MBSs this drives the capacity slot on/off; for fixed sites
    % the gene is bit-flipped but ignored at evaluation time. Skipped when
    % the GA does not control MBS capacity so the RNG stream matches the
    % toggle-off behaviour.
    if gaControlsMbsCapacity
        for j = 1:n_mbs
            if rand() < prob
                mbsIdx = fbsCount + j;
                individual(mbsIdx) = 1 - round(individual(mbsIdx));
                mutationFlags(mbsIdx) = true;
            end
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
