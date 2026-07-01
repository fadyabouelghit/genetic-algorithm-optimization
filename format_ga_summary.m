function s = format_ga_summary(bestIndividual, bestFitness, history, params)
%FORMAT_GA_SUMMARY  Render the post-GA summary block as a plain text string.
%
% Mirrors the "=== Optimization Complete ===" block printed by
% optimizeBaseStation.m, but as a string. Useful for shipping the result
% via email / push without scraping stdout.

    if ~isfield(history, 'rawMetrics') || isempty(history.rawMetrics)
        rm = struct('numUsers', NaN, 'fbsUsers', NaN, 'mbsUsers', NaN, ...
                    'transmittedPower', NaN, 'avgRate', NaN);
    else
        rm = history.rawMetrics;
    end

    numBS  = params.numBS;
    numMbs = max(0, numel(bestIndividual) - 6 * numBS);

    if isstruct(history) && isfield(history, 'time') && isfield(history.time, 'total')
        totalMin = history.time.total / 60;
    else
        totalMin = NaN;
    end

    lines = strings(0, 1);
    lines(end+1, 1) = "=== Optimization Complete ===";
    lines(end+1, 1) = sprintf('Total time: %.1f minutes', totalMin);
    lines(end+1, 1) = sprintf('Best Fitness: %.2f users', bestFitness);
    lines(end+1, 1) = "The best individual performance:";
    lines(end+1, 1) = sprintf('Total Connected Users: %d', round(rm.numUsers));
    lines(end+1, 1) = sprintf(' - FBS-connected Users: %d', round(rm.fbsUsers));
    lines(end+1, 1) = sprintf(' - MBS-connected Users: %d', round(rm.mbsUsers));
    lines(end+1, 1) = sprintf('Total Transmitted Power: %.2f W', rm.transmittedPower);
    lines(end+1, 1) = sprintf('Avg. Sum Rate: %.2f bps/Hz', rm.avgRate);
    lines(end+1, 1) = "Final Parameters:";

    for bs = 1:numBS
        base = (bs - 1) * 6;
        lines(end+1, 1) = sprintf('    BS%d X (m)           %g',  bs, bestIndividual(base + 1)); %#ok<AGROW>
        lines(end+1, 1) = sprintf('    BS%d Y (m)           %g',  bs, bestIndividual(base + 2)); %#ok<AGROW>
        lines(end+1, 1) = sprintf('    BS%d Z (m)           %g',  bs, bestIndividual(base + 3)); %#ok<AGROW>
        lines(end+1, 1) = sprintf('    BS%d Power (W)       %g',  bs, bestIndividual(base + 4)); %#ok<AGROW>
        lines(end+1, 1) = sprintf('    BS%d Power Status    %d',  bs, round(bestIndividual(base + 5))); %#ok<AGROW>
        lines(end+1, 1) = sprintf('    BS%d fbsFreqFlag     %d',  bs, double(bestIndividual(base + 6) >= 0.5)); %#ok<AGROW>
    end
    for mIdx = 1:numMbs
        lines(end+1, 1) = sprintf('    MBS%d CapacityOn     %d', ...
            mIdx, double(bestIndividual(6 * numBS + mIdx) >= 0.5)); %#ok<AGROW>
    end

    s = strjoin(lines, newline);
end
