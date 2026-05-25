function [fitness, details] = evaluatePopulation(l, population, verbose, n_fbs, spaceLimit, containsMbs, mbs_params, antennaObjectMbs, bounds, mbsCache, targetIdx, weightParams)
% EVALUATEPOPULATION Computes the fitness of each individual in the population
% based on normalized number of connected users and transmission power.
% fitness(i) = β * log(normUsers + ε) - γ * log(normPower + ε)          --- option 1
% fitness(i) = (norm_numUsers ^ beta) * ((1 - norm_power) ^ gamma);     --- option 2


    if nargin < 12 || isempty(weightParams)
        weightParams = struct();
    end
    if ~isfield(weightParams, 'beta'), weightParams.beta = 1; end
    if ~isfield(weightParams, 'gamma'), weightParams.gamma = 1; end
    if ~isfield(weightParams, 'epsilon'), weightParams.epsilon = 1e-3; end
    if ~isfield(weightParams, 'fbsWeight'), weightParams.fbsWeight = 0; end
    if ~isfield(weightParams, 'fbsExponent'), weightParams.fbsExponent = 1; end
    if ~isfield(weightParams, 'maxUsers'), weightParams.maxUsers = 1000; end
    if ~isfield(weightParams, 'sinrThreshold'), weightParams.sinrThreshold = 5; end
    if ~isfield(weightParams, 'gaControlsFbsBand'), weightParams.gaControlsFbsBand = true; end
    if ~isfield(weightParams, 'gaControlsMbsCapacity'), weightParams.gaControlsMbsCapacity = true; end
    if ~isfield(weightParams, 'mbsBandPolicy'), weightParams.mbsBandPolicy = []; end
    beta = weightParams.beta;
    gamma = weightParams.gamma;
    epsilon = weightParams.epsilon;
    fbsWeight = weightParams.fbsWeight;
    fbsExponent = weightParams.fbsExponent;
    maxUsers = weightParams.maxUsers;
    sinrThreshold = weightParams.sinrThreshold;
    gaControlsFbsBand = weightParams.gaControlsFbsBand;
    gaControlsMbsCapacity = weightParams.gaControlsMbsCapacity;
    mbsBandPolicy = weightParams.mbsBandPolicy;

    fitness = zeros(size(population,1), 1);
    numIndividuals = size(population,1);
    details = struct(...
        'numUsers', zeros(numIndividuals, 1), ...
        'transmittedPower', zeros(numIndividuals, 1), ...
        'avgRate', zeros(numIndividuals, 1), ...
        'fbsUsers', zeros(numIndividuals, 1), ...
        'mbsUsers', zeros(numIndividuals, 1), ...
        'activeFbs', zeros(numIndividuals, 1));

    fbsBoundRows = bounds(1:6*n_fbs, :);
    powerBounds = fbsBoundRows(4:6:end, :);
    maxPower = sum(powerBounds(:,2));
    minPower = 0;

    mbs_x = mbs_params(1,:); mbs_y = mbs_params(2,:);
    mbs_height = mbs_params(3,:); mbs_power = mbs_params(4,:);
    blockSize = 6;
    numMbs = containsMbs * size(mbs_params, 2);
    fbsCount = blockSize * n_fbs;
    expectedLen = fbsCount + numMbs;

    % Resolve MBS band policy. When no policy is supplied, treat every MBS
    % site as a "base" MBS (dual-band, GA-gated capacity). Femtos only show
    % up when the caller supplies the policy explicitly.
    if isempty(mbsBandPolicy)
        siteIsBaseMbs = true(1, numMbs);
        siteFixedBand = zeros(1, numMbs);
    else
        siteIsBaseMbs = logical(reshape(mbsBandPolicy.siteIsBaseMbs, 1, []));
        siteFixedBand = double(reshape(mbsBandPolicy.siteFixedBand, 1, []));
        assert(numel(siteIsBaseMbs) == numMbs, ...
            'mbsBandPolicy.siteIsBaseMbs must have %d entries (got %d).', ...
            numMbs, numel(siteIsBaseMbs));
        assert(numel(siteFixedBand) == numMbs, ...
            'mbsBandPolicy.siteFixedBand must have %d entries (got %d).', ...
            numMbs, numel(siteFixedBand));
    end

    for i = 1:size(population,1)
        ind = population(i,:);
        if numel(ind) ~= expectedLen
            error('evaluatePopulation expects %d decision vars (6 per FBS + 1 per MBS), got %d.', ...
                expectedLen, numel(ind));
        end
        fbsBlock = ind(1:fbsCount);
        x = fbsBlock(1:blockSize:end);
        y = fbsBlock(2:blockSize:end);
        z = fbsBlock(3:blockSize:end);
        power = fbsBlock(4:blockSize:end);
        power_status = fbsBlock(5:blockSize:end);
        if gaControlsFbsBand
            fbsFreqFlags = double(fbsBlock(6:blockSize:end) >= 0.5);
        else
            fbsFreqFlags = zeros(1, n_fbs);
        end

        fbsAntennaEval = repmat(l(1), 1, n_fbs);
        if numel(l) >= 2
            capMask = (fbsFreqFlags >= 0.5);
            fbsAntennaEval(capMask) = l(2);
        end

        if numMbs > 0
            mbsGenes = double(ind(fbsCount + (1:numMbs)) >= 0.5);
            if ~gaControlsMbsCapacity
                mbsGenes(:) = 0;
            end
            [mbsSlotMap, mbsSlotBands] = build_mbs_slots(siteIsBaseMbs, siteFixedBand, mbsGenes);
        else
            mbsSlotMap   = zeros(0, 3);
            mbsSlotBands = zeros(1, 0);
        end
        bsBandIds = [fbsFreqFlags, mbsSlotBands];

        [~, ~, numUsers, transmittedPower, avg_rate_connected_bpsHz, fbsUsers, mbsUsers, sum_rate_connected_bpsHz] = SINREvaluation(fbsAntennaEval, power_status, ...
            x, y, z, n_fbs, power, ...
            mbs_y, mbs_x, mbs_height, mbs_power, ...
            0, spaceLimit(1), 0, spaceLimit(2), maxUsers, sinrThreshold, containsMbs, antennaObjectMbs, mbsCache, bsBandIds, mbsSlotMap);

        details.numUsers(i) = numUsers;
        details.transmittedPower(i) = transmittedPower;
        details.avgRate(i) = avg_rate_connected_bpsHz;
        details.fbsUsers(i) = fbsUsers;
        details.mbsUsers(i) = mbsUsers;
        details.activeFbs(i) = sum(power_status >= 0.5);

        if targetIdx == 1
            norm_numUsers = numUsers / maxUsers;
            norm_numUsers = min(max(norm_numUsers, 0), 1);
            powerRange = maxPower - minPower;
            if powerRange <= 0
                powerRange = 1;
            end
            norm_power = (transmittedPower - minPower) / powerRange;
            norm_power = min(max(norm_power, 0), 1);

            base = (norm_numUsers ^ beta) * ((1 - norm_power) ^ gamma);
            totalConnections = max(numUsers, 1);
            fbsShare = fbsUsers / totalConnections;
            fbsShare = min(max(fbsShare, 0), 1);
            fbsTerm = fbsShare ^ fbsExponent;
            fitness(i) = (1 - fbsWeight) * base + fbsWeight * fbsTerm;
        
        elseif targetIdx == 2
            fitness(i) = sum_rate_connected_bpsHz;

        end
    end
end

function [slotMap, slotBands] = build_mbs_slots(siteIsBaseMbs, siteFixedBand, mbsGenes)
% Expand MBS sites into per-band slots.
%   - Base MBS site j -> two slots: [j, 1 (cov), 1 always-on] and
%     [j, 2 (cap), mbsGenes(j) (capacity gated by GA gene)].
%   - Fixed (femto) site j -> one slot: [j, siteFixedBand(j)+1, 1].
% slotBands is the band-id vector (0/1) that lines up with the slot rows
% and is appended to bsBandIds for same-band interference accounting.

    numMbs = numel(siteIsBaseMbs);
    rows = cell(1, numMbs);
    bands = cell(1, numMbs);
    for j = 1:numMbs
        if siteIsBaseMbs(j)
            rows{j} = [ j, 1, 1; ...
                        j, 2, mbsGenes(j) ];
            bands{j} = [0, 1];
        else
            bandIdx = siteFixedBand(j) + 1;
            rows{j}  = [ j, bandIdx, 1 ];
            bands{j} = siteFixedBand(j);
        end
    end
    if isempty(rows)
        slotMap   = zeros(0, 3);
        slotBands = zeros(1, 0);
    else
        slotMap   = vertcat(rows{:});
        slotBands = horzcat(bands{:});
    end
end
