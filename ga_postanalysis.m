function ga_postanalysis(resultsDir)
% GA_POSTANALYSIS  Post-hoc visualization and statistical analysis of GA runs
%
%   ga_postanalysis(resultsDir)
%
%   Loads results.mat from the specified directory (output of
%   ga_performance_harness.m) and produces:
%
%   1. Convergence comparison  - overlaid best-fitness traces grouped by config
%   2. Boxplot dashboard       - best fitness distribution per hyperparameter
%   3. Heatmaps                - pairwise hyperparameter interaction effects
%   4. Statistical ranking     - configs ranked by mean fitness (with CI)
%   5. Timing analysis         - cost per generation and per configuration
%   6. Convergence diagnostics - calls ga_convergence_diagnostics for top configs
%
%   All figures are saved to resultsDir/figures/.
%
%   NOTE: This file does NOT modify any existing project files.

    if nargin < 1
        resultsDir = uigetdir(pwd, 'Select results directory');
        if resultsDir == 0; return; end
    end

    fprintf('Loading results from: %s\n', resultsDir);
    data = load(fullfile(resultsDir, 'results.mat'));
    allRuns  = data.allRuns;
    configs  = data.configs;
    nTrials  = data.nTrials;
    metadata = data.metadata;

    figDir = fullfile(resultsDir, 'figures');
    if ~exist(figDir, 'dir'); mkdir(figDir); end

    nConfigs = numel(configs);
    nRuns    = numel(allRuns);

    fprintf('Configs: %d | Trials/config: %d | Total runs: %d\n', nConfigs, nTrials, nRuns);
    fprintf('Total wall time was: %.1f hours\n\n', metadata.totalWallTime / 3600);

    %% ======== 1. Aggregate Statistics Per Config ========
    % Aggregate both fitness (for within-preset comparisons) and raw metrics
    % (for cross-preset, scale-independent comparisons).
    configStats = struct();
    for c = 1:nConfigs
        mask = [allRuns.configId] == c & strcmp({allRuns.status}, 'ok');
        fits = [allRuns(mask).bestFitness];

        configStats(c).configId    = c;
        configStats(c).meanFitness = mean(fits, 'omitnan');
        configStats(c).stdFitness  = std(fits, 0, 'omitnan');
        configStats(c).minFitness  = min(fits);
        configStats(c).maxFitness  = max(fits);
        configStats(c).medianFit   = median(fits, 'omitnan');
        configStats(c).nOk         = sum(mask);

        elapsedArr = [allRuns(mask).elapsed_sec];
        configStats(c).meanTime = mean(elapsedArr, 'omitnan');

        % Raw physical metrics (scale-independent)
        if isfield(allRuns, 'numUsers')
            users = [allRuns(mask).numUsers];
            configStats(c).meanUsers   = mean(users, 'omitnan');
            configStats(c).stdUsers    = std(users, 0, 'omitnan');
            configStats(c).meanFbsU    = mean([allRuns(mask).fbsUsers], 'omitnan');
            configStats(c).meanMbsU    = mean([allRuns(mask).mbsUsers], 'omitnan');
            configStats(c).meanPower   = mean([allRuns(mask).transmittedPower], 'omitnan');
            configStats(c).stdPower    = std([allRuns(mask).transmittedPower], 0, 'omitnan');
            configStats(c).meanRate    = mean([allRuns(mask).avgRate], 'omitnan');
            configStats(c).stdRate     = std([allRuns(mask).avgRate], 0, 'omitnan');
        else
            configStats(c).meanUsers = NaN; configStats(c).stdUsers = NaN;
            configStats(c).meanFbsU = NaN;  configStats(c).meanMbsU = NaN;
            configStats(c).meanPower = NaN; configStats(c).stdPower = NaN;
            configStats(c).meanRate = NaN;  configStats(c).stdRate = NaN;
        end
    end

    %% ======== 2. Rank Configs by Raw Metrics (Connected Users) ========
    % Primary ranking uses connected users — comparable across all configs
    [~, rankOrder] = sort([configStats.meanUsers], 'descend');

    fprintf('=== Top 10 Configurations (by mean connected users) ===\n');
    fprintf('%-5s %-5s %-4s %-5s %-5s %-5s %-12s %-7s %-7s %-8s %-7s %-8s\n', ...
        'Rank', 'CfgID', 'nBS', 'Pop', 'XO', 'Mut', 'Preset', 'Users', 'StdU', 'Power(W)', 'Rate', 'Time(s)');
    for r = 1:min(10, nConfigs)
        c = rankOrder(r);
        cfg = configs(c);
        fprintf('%-5d %-5d %-4d %-5d %-5.2f %-5.2f %-12s %-7.0f %-7.1f %-8.2f %-7.2f %-8.1f\n', ...
            r, c, cfg.numBS, cfg.populationSize, cfg.crossoverProb, cfg.mutationProb, ...
            cfg.fitnessPreset, configStats(c).meanUsers, configStats(c).stdUsers, ...
            configStats(c).meanPower, configStats(c).meanRate, configStats(c).meanTime);
    end

    %% ======== 3. Boxplots: Raw Metrics by Hyperparameter ========
    % Uses connected users (scale-independent) as the primary metric,
    % plus separate panels for power and rate.
    hyperparams = {'numBS', 'populationSize', 'crossoverProb', 'mutationProb', ...
                   'mutationScale', 'fitnessPreset'};
    hpLabels    = {'Num FBS', 'Population Size', 'Crossover Prob', 'Mutation Prob', ...
                   'Mutation Scale', 'Fitness Preset'};

    % Metrics to boxplot: {field on allRuns, ylabel, filename suffix}
    metricsToPlot = {
        'numUsers',         'Connected Users',        'users';
        'transmittedPower', 'Transmitted Power (W)',   'power';
        'avgRate',          'Avg Sum Rate (bps/Hz)',   'rate';
        'bestFitness',      'Fitness (warning: scale varies by preset)', 'fitness';
    };

    for hi = 1:numel(hyperparams)
        hp = hyperparams{hi};

        for mi = 1:size(metricsToPlot, 1)
            metricField  = metricsToPlot{mi, 1};
            metricYLabel = metricsToPlot{mi, 2};
            metricSuffix = metricsToPlot{mi, 3};

            fig = figure('Name', sprintf('Boxplot %s: %s', metricSuffix, hpLabels{hi}), ...
                'Visible', 'off', 'Position', [100 100 800 500]);

            if strcmp(hp, 'fitnessPreset')
                groupLabels = {};
                metricVals  = [];
                for r = 1:nRuns
                    if strcmp(allRuns(r).status, 'ok')
                        c = allRuns(r).configId;
                        groupLabels{end+1} = configs(c).fitnessPreset; %#ok<AGROW>
                        metricVals(end+1)  = allRuns(r).(metricField); %#ok<AGROW>
                    end
                end
                boxplot(metricVals, groupLabels);
            else
                groupVals  = [];
                metricVals = [];
                for r = 1:nRuns
                    if strcmp(allRuns(r).status, 'ok')
                        c = allRuns(r).configId;
                        groupVals(end+1)  = configs(c).(hp); %#ok<AGROW>
                        metricVals(end+1) = allRuns(r).(metricField); %#ok<AGROW>
                    end
                end
                [uVals, ~, groupIdx] = unique(groupVals);
                boxplot(metricVals, groupIdx, 'Labels', arrayfun(@num2str, uVals, 'UniformOutput', false));
            end

            xlabel(hpLabels{hi});
            ylabel(metricYLabel);
            title(sprintf('%s by %s', metricYLabel, hpLabels{hi}));
            grid on;

            saveas(fig, fullfile(figDir, sprintf('boxplot_%s_%s.png', metricSuffix, hp)));
            saveas(fig, fullfile(figDir, sprintf('boxplot_%s_%s.fig', metricSuffix, hp)));
            close(fig);
        end
    end

    %% ======== 4. Heatmap: Pairwise Interactions ========
    % Focus on the most impactful pairs
    pairsToPlot = {
        'crossoverProb', 'mutationProb',  'Crossover vs Mutation';
        'mutationProb',  'mutationScale', 'Mutation Prob vs Scale';
        'beta',          'fbsWeight',     'Beta vs FBS Weight';
        'populationSize','numBS',         'Pop Size vs Num FBS';
    };

    for pi = 1:size(pairsToPlot, 1)
        p1 = pairsToPlot{pi, 1};
        p2 = pairsToPlot{pi, 2};
        pLabel = pairsToPlot{pi, 3};

        % Get unique values
        v1 = unique([configs.(p1)]);
        v2 = unique([configs.(p2)]);

        % Use connected users (scale-independent) for heatmaps
        heatData = NaN(numel(v2), numel(v1));
        heatCounts = zeros(numel(v2), numel(v1));

        for c = 1:nConfigs
            i1 = find(v1 == configs(c).(p1));
            i2 = find(v2 == configs(c).(p2));
            if ~isempty(i1) && ~isempty(i2)
                prev = heatData(i2, i1);
                cur  = configStats(c).meanUsers;
                if isnan(prev)
                    heatData(i2, i1) = cur;
                    heatCounts(i2, i1) = 1;
                else
                    % Running average across configs sharing this cell
                    n = heatCounts(i2, i1);
                    heatData(i2, i1) = (prev * n + cur) / (n + 1);
                    heatCounts(i2, i1) = n + 1;
                end
            end
        end

        fig = figure('Name', ['Heatmap: ' pLabel], 'Visible', 'off', ...
            'Position', [100 100 700 500]);
        imagesc(v1, v2, heatData);
        cb = colorbar;
        ylabel(cb, 'Mean Connected Users');
        set(gca, 'YDir', 'normal');
        xlabel(p1, 'Interpreter', 'none');
        ylabel(p2, 'Interpreter', 'none');
        title(sprintf('Mean Connected Users: %s', pLabel));
        colormap(parula);

        % Add text labels on cells
        for ii = 1:numel(v1)
            for jj = 1:numel(v2)
                if ~isnan(heatData(jj, ii))
                    text(v1(ii), v2(jj), sprintf('%.3f', heatData(jj, ii)), ...
                        'HorizontalAlignment', 'center', 'FontSize', 8);
                end
            end
        end

        saveas(fig, fullfile(figDir, sprintf('heatmap_%s_vs_%s.png', p1, p2)));
        saveas(fig, fullfile(figDir, sprintf('heatmap_%s_vs_%s.fig', p1, p2)));
        close(fig);
    end

    %% ======== 5. Convergence Overlay: Top vs Bottom Configs ========
    nCompare = min(5, nConfigs);
    topIds    = rankOrder(1:nCompare);
    bottomIds = rankOrder(end-nCompare+1:end);

    fig = figure('Name', 'Convergence: Top vs Bottom', 'Visible', 'off', ...
        'Position', [100 100 1000 600]);

    subplot(1,2,1); hold on; title('Top Configs'); grid on;
    xlabel('Generation'); ylabel('Best Fitness');
    cmap = lines(nCompare);
    for ki = 1:nCompare
        c = topIds(ki);
        mask = [allRuns.configId] == c & strcmp({allRuns.status}, 'ok');
        runs = allRuns(mask);
        for r = 1:numel(runs)
            trace = runs(r).bestFitnessTrace;
            if numel(trace) > 1
                plot(1:numel(trace), trace, '-', 'Color', [cmap(ki,:), 0.4], ...
                    'HandleVisibility', 'off');
            end
        end
        plot(NaN, NaN, '-', 'Color', cmap(ki,:), 'LineWidth', 2, ...
            'DisplayName', sprintf('Cfg %d (%.0f users)', c, configStats(c).meanUsers));
    end
    legend('Location', 'best', 'FontSize', 7);

    subplot(1,2,2); hold on; title('Bottom Configs'); grid on;
    xlabel('Generation'); ylabel('Best Fitness');
    for ki = 1:nCompare
        c = bottomIds(ki);
        mask = [allRuns.configId] == c & strcmp({allRuns.status}, 'ok');
        runs = allRuns(mask);
        for r = 1:numel(runs)
            trace = runs(r).bestFitnessTrace;
            if numel(trace) > 1
                plot(1:numel(trace), trace, '-', 'Color', [cmap(ki,:), 0.4], ...
                    'HandleVisibility', 'off');
            end
        end
        plot(NaN, NaN, '-', 'Color', cmap(ki,:), 'LineWidth', 2, ...
            'DisplayName', sprintf('Cfg %d (%.0f users)', c, configStats(c).meanUsers));
    end
    legend('Location', 'best', 'FontSize', 7);

    saveas(fig, fullfile(figDir, 'convergence_top_vs_bottom.png'));
    saveas(fig, fullfile(figDir, 'convergence_top_vs_bottom.fig'));
    close(fig);

    %% ======== 6. Timing Analysis ========
    fig = figure('Name', 'Timing Analysis', 'Visible', 'off', ...
        'Position', [100 100 900 400]);

    subplot(1,2,1);
    times = [configStats.meanTime];
    bar(times(rankOrder));
    xlabel('Config (ranked by fitness)');
    ylabel('Mean Run Time (s)');
    title('Run Time by Config Rank');
    grid on;

    subplot(1,2,2);
    scatter([configs.populationSize], times, 60, [configs.numBS], 'filled');
    colorbar; ylabel(colorbar, 'Num FBS');
    xlabel('Population Size');
    ylabel('Mean Run Time (s)');
    title('Time vs Population Size (color = numBS)');
    grid on;

    saveas(fig, fullfile(figDir, 'timing_analysis.png'));
    saveas(fig, fullfile(figDir, 'timing_analysis.fig'));
    close(fig);

    %% ======== 7. Export Ranked Config Table ========
    rankTable = table();
    for r = 1:nConfigs
        c = rankOrder(r);
        row = struct();
        row.rank          = r;
        row.configId      = c;
        row.numBS         = configs(c).numBS;
        row.populationSize= configs(c).populationSize;
        row.crossoverProb = configs(c).crossoverProb;
        row.mutationProb  = configs(c).mutationProb;
        row.mutationScale = configs(c).mutationScale;
        row.fitnessPreset = string(configs(c).fitnessPreset);
        row.beta          = configs(c).beta;
        row.gamma         = configs(c).gamma;
        row.fbsWeight     = configs(c).fbsWeight;
        row.fbsExponent   = configs(c).fbsExponent;
        row.meanFitness   = configStats(c).meanFitness;
        row.stdFitness    = configStats(c).stdFitness;
        row.medianFitness = configStats(c).medianFit;
        row.minFitness    = configStats(c).minFitness;
        row.maxFitness    = configStats(c).maxFitness;
        % Raw metrics (scale-independent)
        row.meanUsers     = configStats(c).meanUsers;
        row.stdUsers      = configStats(c).stdUsers;
        row.meanFbsUsers  = configStats(c).meanFbsU;
        row.meanMbsUsers  = configStats(c).meanMbsU;
        row.meanPower     = configStats(c).meanPower;
        row.stdPower      = configStats(c).stdPower;
        row.meanRate      = configStats(c).meanRate;
        row.stdRate        = configStats(c).stdRate;
        row.meanTime_sec  = configStats(c).meanTime;
        row.nOkTrials     = configStats(c).nOk;

        rankTable = [rankTable; struct2table(row)]; %#ok<AGROW>
    end
    writetable(rankTable, fullfile(resultsDir, 'ranked_configs.csv'));

    %% ======== 8. Summary Figure (Dashboard) ========
    fig = figure('Name', 'GA Performance Dashboard', 'Visible', 'off', ...
        'Position', [50 50 1400 900]);

    % 2x3 grid — all panels use raw metrics (connected users) for
    % scale-independent comparison across fitness presets.
    okRuns = allRuns(strcmp({allRuns.status}, 'ok'));

    subplot(2,3,1);
    histogram([okRuns.numUsers], 30);
    xlabel('Connected Users'); ylabel('Count');
    title('Connected Users Distribution (All Runs)'); grid on;

    subplot(2,3,2);
    meanU = [configStats.meanUsers];
    stdU  = [configStats.stdUsers];
    errorbar(1:nConfigs, meanU(rankOrder), stdU(rankOrder), '.', 'MarkerSize', 4);
    xlabel('Config Rank'); ylabel('Mean Users \pm Std');
    title('Config Ranking (Connected Users)'); grid on;

    subplot(2,3,3);
    scatter(stdU, meanU, 30, 'filled');
    xlabel('Std Users'); ylabel('Mean Users');
    title('Mean vs Variability (Users)'); grid on;

    subplot(2,3,4);
    nbs_vals = unique([configs.numBS]);
    hold on;
    for ni = 1:numel(nbs_vals)
        nbv = nbs_vals(ni);
        mask = [configs.numBS] == nbv;
        cfgIds = find(mask);
        users = arrayfun(@(c) configStats(c).meanUsers, cfgIds);
        scatter(ones(size(users)) * nbv + randn(size(users))*0.05, users, 20, 'filled');
    end
    xlabel('Num FBS'); ylabel('Mean Connected Users');
    title('Users by Num FBS'); grid on; hold off;

    subplot(2,3,5);
    scatter([okRuns.elapsed_sec], [okRuns.numUsers], 15, 'filled', 'MarkerFaceAlpha', 0.3);
    xlabel('Run Time (s)'); ylabel('Connected Users');
    title('Users vs Compute Cost'); grid on;

    subplot(2,3,6);
    errorRuns = sum(strcmp({allRuns.status}, 'error'));
    text(0.5, 0.7, sprintf('Total runs: %d', nRuns), ...
        'HorizontalAlignment', 'center', 'FontSize', 12, 'Units', 'normalized');
    text(0.5, 0.55, sprintf('Successful: %d (%.1f%%)', nRuns-errorRuns, (1-errorRuns/nRuns)*100), ...
        'HorizontalAlignment', 'center', 'FontSize', 12, 'Units', 'normalized');
    text(0.5, 0.4, sprintf('Wall time: %.1f hours', metadata.totalWallTime/3600), ...
        'HorizontalAlignment', 'center', 'FontSize', 12, 'Units', 'normalized');
    text(0.5, 0.25, sprintf('Best users: %.0f', max([okRuns.numUsers])), ...
        'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold', 'Units', 'normalized');
    text(0.5, 0.1, sprintf('Best rate: %.2f bps/Hz', max([okRuns.avgRate])), ...
        'HorizontalAlignment', 'center', 'FontSize', 12, 'Units', 'normalized');
    title('Summary'); axis off;

    saveas(fig, fullfile(figDir, 'dashboard.png'));
    saveas(fig, fullfile(figDir, 'dashboard.fig'));
    close(fig);

    fprintf('\n=== Post-Analysis Complete ===\n');
    fprintf('Figures saved to: %s\n', figDir);
    fprintf('Ranked configs saved to: %s\n', fullfile(resultsDir, 'ranked_configs.csv'));

end
