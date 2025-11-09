function [fitness, details] = evaluatePopulation(l, population, verbose, n_fbs, spaceLimit, containsMbs, mbs_params, antennaObjectMbs, bounds, mbsCache, targetIdx)
% EVALUATEPOPULATION Computes the fitness of each individual in the population
% based on normalized number of connected users and transmission power.

    fitness = zeros(size(population,1), 1);
    numIndividuals = size(population,1);
    details = struct(...
        'numUsers', zeros(numIndividuals, 1), ...
        'transmittedPower', zeros(numIndividuals, 1), ...
        'avgRate', zeros(numIndividuals, 1));

    alpha = 0; 
    maxPower = bounds(4,2); 
    minPower = bounds(4,1); 
    maxUsers = 1000;

    mbs_x = mbs_params(1,:); mbs_y = mbs_params(2,:);
    mbs_height = mbs_params(3,:); mbs_power = mbs_params(4,:);

    for i = 1:size(population,1)
        ind = population(i,:);
        x = ind(1:5:end); y = ind(2:5:end); z = ind(3:5:end);
        power = ind(4:5:end); power_status = ind(5:5:end);

        [~, ~, numUsers, transmittedPower, avg_rate_connected_bpsHz] = SINREvaluation(l, power_status, ...
            x, y, z, n_fbs, power, ...
            mbs_y, mbs_x, mbs_height, mbs_power, ...
            0, spaceLimit(1), 0, spaceLimit(2), maxUsers, 5, containsMbs, antennaObjectMbs, mbsCache);

        details.numUsers(i) = numUsers;
        details.transmittedPower(i) = transmittedPower;
        details.avgRate(i) = avg_rate_connected_bpsHz;

        if targetIdx == 1
            norm_numUsers = numUsers / maxUsers;
            % norm_numUsers = numUsers;
            norm_power = (transmittedPower - minPower) / (maxPower - minPower);
            normalized_cost = norm_numUsers - alpha * norm_power;
            fitness(i) = (normalized_cost + alpha) / (1 + alpha);
        
        elseif targetIdx == 2
            fitness(i) = avg_rate_connected_bpsHz;
        
        end
    end
end
