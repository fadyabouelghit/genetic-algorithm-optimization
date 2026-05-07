function [child1, child2, crossoverFlag] = crossover_blend(parent1, parent2, prob, bounds, n_fbs, controls)

    crossoverFlag = 0;
    alpha = 0.5;
    blockSize = 6;

    if nargin < 5 || isempty(n_fbs)
        if mod(length(parent1), blockSize) ~= 0
            error(['crossover_blend: chromosome length is not a multiple of 6 and ' ...
                'n_fbs was not supplied; pass n_fbs explicitly.']);
        end
        n_fbs = floor(length(parent1) / blockSize);
    end
    fbsCount = blockSize * n_fbs;
    n_mbs = length(parent1) - fbsCount;
    if n_mbs < 0
        error('crossover_blend: chromosome shorter than 6*n_fbs (got %d, n_fbs=%d).', length(parent1), n_fbs);
    end
    if nargin < 6 || isempty(controls)
        controls = struct();
    end
    if ~isfield(controls, 'fbsBand'), controls.fbsBand = true; end
    if ~isfield(controls, 'mbsBand'), controls.mbsBand = true; end

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

            if controls.fbsBand
                if parent1(freq_idx) == parent2(freq_idx)
                    child1(freq_idx) = parent1(freq_idx);
                    child2(freq_idx) = parent1(freq_idx);
                else
                    child1(freq_idx) = randi([0, 1]);
                    child2(freq_idx) = 1 - child1(freq_idx);
                end
            else
                % Toggle off: hold at 0, no RNG consumed.
                child1(freq_idx) = 0;
                child2(freq_idx) = 0;
            end
        end

        for bs = 1:n_fbs
            idx = (bs-1)*blockSize + 1 : (bs-1)*blockSize + 4;
            child1(idx) = clampToBounds(child1(idx), bounds(idx,:));
            child2(idx) = clampToBounds(child2(idx), bounds(idx,:));
        end

        % MBS frequency-flag suffix: same equal-or-randomize logic per gene.
        % Skip the loop entirely when GA does not control MBS bands so the RNG
        % stream matches the pre-feature behavior.
        if controls.mbsBand
            for j = 1:n_mbs
                mIdx = fbsCount + j;
                p1g = round(parent1(mIdx));
                p2g = round(parent2(mIdx));
                if p1g == p2g
                    child1(mIdx) = p1g;
                    child2(mIdx) = p1g;
                else
                    child1(mIdx) = randi([0, 1]);
                    child2(mIdx) = 1 - child1(mIdx);
                end
            end
        end
    else
        child1 = parent1;
        child2 = parent2;
    end
end
