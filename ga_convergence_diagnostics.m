function diagnostics = ga_convergence_diagnostics(history, configLabel)
% GA_CONVERGENCE_DIAGNOSTICS  Compute convergence metrics from a single GA run.
%
%   diagnostics = ga_convergence_diagnostics(history)
%   diagnostics = ga_convergence_diagnostics(history, configLabel)
%
%   Inputs:
%       history     - The history struct returned by optimizeBaseStation
%       configLabel - (optional) string label for plots/display
%
%   Outputs:
%       diagnostics - struct with the following fields:
%
%   --- Convergence metrics ---
%       converged           : true if fitness plateau detected
%       convergenceGen      : generation at which convergence was detected
%       improvementRate     : per-generation improvement in best fitness
%       relativeImprovement : (final - initial) / initial best fitness
%       plateauLength       : number of gens at end with no meaningful change
%       earlyGain           : fraction of total improvement in first 30% of gens
%
%   --- Diversity metrics ---
%       finalStdFitness     : population fitness std at final generation
%       diversityRatio      : std / mean of fitness at final generation
%       diversityTrace      : std(fitness) over all generations
%
%   --- Stability metrics ---
%       bestFitnessCV       : coeff of variation of best fitness (last 30%)
%       oscillationCount    : number of sign changes in bestFitness diff
%
%   This function does NOT modify any existing files.

    if nargin < 2 || isempty(configLabel)
        configLabel = 'GA Run';
    end

    nGens = numel(history.bestFitness);
    bestF = history.bestFitness(:);
    avgF  = history.avgFitness(:);
    stdF  = history.stdFitness(:);

    diagnostics = struct();
    diagnostics.label = configLabel;
    diagnostics.nGenerations = nGens;

    %% --- Convergence Detection ---
    % Use a sliding window to detect plateau
    windowSize = max(3, round(nGens * 0.2));  % 20% of generations
    plateauThreshold = 1e-4;  % relative change threshold

    % Find first plateau
    diagnostics.converged = false;
    diagnostics.convergenceGen = nGens;
    diagnostics.plateauLength = 0;

    for g = windowSize:nGens
        windowSlice = bestF(g-windowSize+1:g);
        relChange = abs(windowSlice(end) - windowSlice(1)) / (abs(windowSlice(1)) + 1e-12);
        if relChange < plateauThreshold
            diagnostics.converged = true;
            diagnostics.convergenceGen = g - windowSize + 1;
            break;
        end
    end

    % Plateau at the tail
    tailDiffs = abs(diff(bestF(max(1, nGens-windowSize):nGens)));
    diagnostics.plateauLength = sum(tailDiffs < plateauThreshold * abs(bestF(end)));

    %% --- Improvement Rate ---
    diffs = diff(bestF);
    diagnostics.improvementRate = diffs;  % per-generation
    diagnostics.meanImprovementRate = mean(diffs);

    if abs(bestF(1)) > 1e-12
        diagnostics.relativeImprovement = (bestF(end) - bestF(1)) / abs(bestF(1));
    else
        diagnostics.relativeImprovement = bestF(end) - bestF(1);
    end

    %% --- Early Gain (how front-loaded is progress?) ---
    earlyBound = max(1, round(nGens * 0.3));
    totalGain  = bestF(end) - bestF(1);
    earlyGain  = bestF(earlyBound) - bestF(1);
    if abs(totalGain) > 1e-12
        diagnostics.earlyGain = earlyGain / totalGain;
    else
        diagnostics.earlyGain = NaN;  % no improvement at all
    end

    %% --- Diversity Metrics ---
    diagnostics.finalStdFitness = stdF(end);
    if abs(avgF(end)) > 1e-12
        diagnostics.diversityRatio = stdF(end) / abs(avgF(end));
    else
        diagnostics.diversityRatio = NaN;
    end
    diagnostics.diversityTrace = stdF;

    %% --- Stability Metrics ---
    tailStart = max(1, round(nGens * 0.7));
    tailBest  = bestF(tailStart:end);
    diagnostics.bestFitnessCV = std(tailBest) / (abs(mean(tailBest)) + 1e-12);

    % Count oscillations (sign changes in the diff of best fitness)
    bestDiffs = diff(bestF);
    signChanges = diff(sign(bestDiffs));
    diagnostics.oscillationCount = sum(abs(signChanges) > 0);

    %% --- Variable Drift (how much does the best solution move?) ---
    if isfield(history, 'bestIndividuals')
        bi = history.bestIndividuals;
        % Euclidean distance between consecutive best individuals
        stepDists = sqrt(sum(diff(bi, 1, 1).^2, 2));
        diagnostics.meanStepDistance   = mean(stepDists);
        diagnostics.totalPathLength    = sum(stepDists);
        diagnostics.directDisplacement = norm(bi(end,:) - bi(1,:));
        if diagnostics.totalPathLength > 1e-12
            diagnostics.pathEfficiency = diagnostics.directDisplacement / diagnostics.totalPathLength;
        else
            diagnostics.pathEfficiency = NaN;
        end
    end

    %% --- Operator Effectiveness ---
    if isfield(history, 'crossovers') && isfield(history, 'mutations')
        diagnostics.totalCrossovers = sum(history.crossovers);
        diagnostics.totalMutations  = sum(history.mutations);
        diagnostics.crossoverMutationRatio = diagnostics.totalCrossovers / ...
            (diagnostics.totalMutations + 1e-12);
    end

    %% --- Timing ---
    if isfield(history, 'time')
        diagnostics.totalTime       = history.time.total;
        diagnostics.meanGenTime     = mean(history.time.generation);
        diagnostics.genTimeVariance = var(history.time.generation);
    end

    %% --- Summary Print ---
    fprintf('\n=== Convergence Diagnostics: %s ===\n', configLabel);
    fprintf('Generations: %d\n', nGens);
    fprintf('Converged: %s (gen %d)\n', string(diagnostics.converged), diagnostics.convergenceGen);
    fprintf('Relative improvement: %.4f\n', diagnostics.relativeImprovement);
    fprintf('Early gain (first 30%%): %.1f%%\n', diagnostics.earlyGain * 100);
    fprintf('Plateau length (tail): %d gens\n', diagnostics.plateauLength);
    fprintf('Final diversity ratio (std/mean): %.4f\n', diagnostics.diversityRatio);
    fprintf('Stability CV (last 30%%): %.6f\n', diagnostics.bestFitnessCV);
    fprintf('Oscillation count: %d\n', diagnostics.oscillationCount);
    if isfield(diagnostics, 'pathEfficiency')
        fprintf('Search path efficiency: %.4f\n', diagnostics.pathEfficiency);
    end
    if isfield(diagnostics, 'totalTime')
        fprintf('Total time: %.1f s (%.2f s/gen avg)\n', ...
            diagnostics.totalTime, diagnostics.meanGenTime);
    end
    fprintf('=======================================\n');

end
