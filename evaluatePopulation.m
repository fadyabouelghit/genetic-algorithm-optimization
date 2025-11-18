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
    beta = weightParams.beta;
    gamma = weightParams.gamma;
    epsilon = weightParams.epsilon;
    fbsWeight = weightParams.fbsWeight;
    fbsExponent = weightParams.fbsExponent;
    maxUsers = weightParams.maxUsers;
    sinrThreshold = weightParams.sinrThreshold;

    fitness = zeros(size(population,1), 1);
    numIndividuals = size(population,1);
    details = struct(...
        'numUsers', zeros(numIndividuals, 1), ...
        'transmittedPower', zeros(numIndividuals, 1), ...
        'avgRate', zeros(numIndividuals, 1), ...
        'fbsUsers', zeros(numIndividuals, 1), ...
        'mbsUsers', zeros(numIndividuals, 1));

    powerBounds = bounds(4:5:end, :);
    maxPower = sum(powerBounds(:,2));
    minPower = 0;

    mbs_x = mbs_params(1,:); mbs_y = mbs_params(2,:);
    mbs_height = mbs_params(3,:); mbs_power = mbs_params(4,:);

    for i = 1:size(population,1)
        ind = population(i,:);
        x = ind(1:5:end); y = ind(2:5:end); z = ind(3:5:end);
        power = ind(4:5:end); power_status = ind(5:5:end);

        [~, ~, numUsers, transmittedPower, avg_rate_connected_bpsHz, fbsUsers, mbsUsers] = SINREvaluation(l, power_status, ...
            x, y, z, n_fbs, power, ...
            mbs_y, mbs_x, mbs_height, mbs_power, ...
            0, spaceLimit(1), 0, spaceLimit(2), maxUsers, sinrThreshold, containsMbs, antennaObjectMbs, mbsCache);

        details.numUsers(i) = numUsers;
        details.transmittedPower(i) = transmittedPower;
        details.avgRate(i) = avg_rate_connected_bpsHz;
        details.fbsUsers(i) = fbsUsers;
        details.mbsUsers(i) = mbsUsers;

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
            fitness(i) = avg_rate_connected_bpsHz;
        
        end
    end
end
