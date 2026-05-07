function population = initializePopulation_uniform(popSize, bounds, n_fbs, n_mbs, controls)

rng(43);
numParams = size(bounds,1);
population = zeros(popSize, numParams);
blockSize = 6;

if nargin < 4 || isempty(n_mbs)
    n_mbs = numParams - blockSize * n_fbs;
end
if nargin < 5 || isempty(controls)
    controls = struct();
end
if ~isfield(controls, 'fbsBand'), controls.fbsBand = true; end
if ~isfield(controls, 'mbsBand'), controls.mbsBand = true; end

expectedParams = blockSize * n_fbs + n_mbs;
if numParams ~= expectedParams
    error('initializePopulation_uniform expects %d params (6 per FBS + 1 per MBS), got %d.', ...
        expectedParams, numParams);
end

fbsCount = blockSize * n_fbs;

for i = 1:fbsCount

    lb = bounds(i,1);
    ub = bounds(i,2);

    posInBlock = mod(i-1, blockSize) + 1;
    if posInBlock == 5 || posInBlock == 6
        % Binary sampling: 0 or 1 (kept regardless of toggle so RNG stream is
        % identical to the legacy initializer; the freq col is zeroed below
        % when the toggle is off).
        population(:,i) = randi([0 1], popSize, 1);
    else
        % Uniform sampling in [lb, ub]
        population(:,i) = lb + (ub - lb) * rand(popSize, 1);
    end
end

% FBS frequency-flag override: hold at 0 (coverage band) when GA does not control it.
if ~controls.fbsBand
    for bs = 1:n_fbs
        population(:, (bs-1)*blockSize + 6) = 0;
    end
end

% MBS frequency-flag suffix: one binary gene per MBS. Sample only if GA controls it;
% otherwise leave at the pre-allocated 0 (no extra RNG consumed).
if controls.mbsBand
    for j = 1:n_mbs
        population(:, fbsCount + j) = randi([0 1], popSize, 1);
    end
end

end
