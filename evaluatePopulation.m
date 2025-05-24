function [fitness, details] = evaluatePopulation(l, population, verbose, n_fbs, spaceLimit, containsMbs, mbs_params, antennaObjectMbs, bounds)
% EVALUATEPOPULATION Computes the fitness of each individual in the population
% based on normalized number of connected users and transmission power.

    fitness = zeros(size(population,1), 1);
    details = struct('numUsers', zeros(size(population,1), 1), 'transmittedPower', zeros(size(population,1), 1));

    alpha = 0; 
    maxPower = bounds(4,2); 
    minPower = bounds(4,1); 
    maxUsers = 100;

    mbs_x = mbs_params(1,:); mbs_y = mbs_params(2,:);
    mbs_height = mbs_params(3,:); mbs_power = mbs_params(4,:);

    for i = 1:size(population,1)
        ind = population(i,:);
        x = ind(1:5:end); y = ind(2:5:end); z = ind(3:5:end);
        power = ind(4:5:end); power_status = ind(5:5:end);

        [~, ~, numUsers, transmittedPower] = SINREvaluation(l, power_status, ...
            x, y, z, n_fbs, power, ...
            mbs_x, mbs_y, mbs_height, mbs_power, ...
            0, spaceLimit, 0, spaceLimit, maxUsers, 5, containsMbs, antennaObjectMbs);

        norm_numUsers = numUsers / maxUsers;
        norm_power = (transmittedPower - minPower) / (maxPower - minPower);
        normalized_cost = norm_numUsers - alpha * norm_power;
        fitness(i) = (normalized_cost + alpha) / (1 + alpha);
    end
end
