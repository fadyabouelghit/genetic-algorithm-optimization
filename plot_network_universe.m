function fig = plot_network_universe(mbs_params, numBaseMbs, W, H, numUsers, varargin)
% PLOT_NETWORK_UNIVERSE Visualize users, macro BSs, and extra fixed BSs.
%
%   Uses the post-swap mbs_params convention from optimize_base_station_ga.m:
%     row 1 = world y, row 2 = world x, row 3 = height, row 4 = power.
%   The first numBaseMbs columns are treated as macro BSs; any remaining
%   columns are treated as extra/additional fixed BSs.
%
%   Name-value options:
%     'Seed'           (default 0)   rng seed used to match SINREvaluation
%     'Title'          (default 'Network universe')
%     'ShowFbsBounds'  (default true) draw the world rectangle
%     'UserPositions'  (default [])  N×2 [x y]; if empty, regenerate from Seed
%     'IsConnected'    (default [])  N×1 logical; if non-empty, split users
%                                    into connected/unconnected categories

    p = inputParser;
    addParameter(p, 'Seed', 0);
    addParameter(p, 'Title', 'Network universe');
    addParameter(p, 'ShowFbsBounds', true);
    addParameter(p, 'UserPositions', []);
    addParameter(p, 'IsConnected', []);
    parse(p, varargin{:});

    mbs_world_x = mbs_params(2, :);
    mbs_world_y = mbs_params(1, :);

    macroX = mbs_world_x(1:numBaseMbs);
    macroY = mbs_world_y(1:numBaseMbs);
    extraX = mbs_world_x(numBaseMbs+1:end);
    extraY = mbs_world_y(numBaseMbs+1:end);

    if isempty(p.Results.UserPositions)
        rngState = rng;
        cleanup = onCleanup(@() rng(rngState));
        rng(p.Results.Seed);
        ux = randi([max(1,0), W], numUsers, 1);
        uy = randi([max(1,0), H], numUsers, 1);
    else
        ux = p.Results.UserPositions(:, 1);
        uy = p.Results.UserPositions(:, 2);
    end

    fig = figure('Name', 'Network universe', 'NumberTitle', 'off');
    hold on; axis equal; grid on; box on;

    if p.Results.ShowFbsBounds
        rectangle('Position', [0 0 W H], 'EdgeColor', [0.4 0.4 0.4], ...
            'LineStyle', '--');
    end

    isConn = p.Results.IsConnected;
    if isempty(isConn)
        hUsers = scatter(ux, uy, 14, [0.20 0.55 0.85], 'filled', ...
            'MarkerFaceAlpha', 0.45, 'DisplayName', sprintf('Users (n=%d)', numel(ux)));
        legendHandles = hUsers;
    else
        isConn = logical(isConn(:));
        hUnconn = scatter(ux(~isConn), uy(~isConn), 14, [0.55 0.55 0.55], 'filled', ...
            'MarkerFaceAlpha', 0.35, ...
            'DisplayName', sprintf('Unconnected (n=%d)', sum(~isConn)));
        hConn = scatter(ux(isConn), uy(isConn), 16, [0.15 0.65 0.30], 'filled', ...
            'MarkerFaceAlpha', 0.65, ...
            'DisplayName', sprintf('Connected to MBS (n=%d)', sum(isConn)));
        legendHandles = [hUnconn, hConn];
    end

    if ~isempty(macroX)
        hMacro = scatter(macroX, macroY, 220, 'p', ...
            'MarkerFaceColor', [0.95 0.75 0.10], ...
            'MarkerEdgeColor', 'k', 'LineWidth', 1.2, ...
            'DisplayName', sprintf('Macro BS (n=%d)', numel(macroX)));
        legendHandles(end+1) = hMacro;
    end

    if ~isempty(extraX)
        hExtra = scatter(extraX, extraY, 140, '^', ...
            'MarkerFaceColor', [0.85 0.30 0.30], ...
            'MarkerEdgeColor', 'k', 'LineWidth', 1.0, ...
            'DisplayName', sprintf('Additional BS (n=%d)', numel(extraX)));
        legendHandles(end+1) = hExtra;
    end

    xlim([0 W]); ylim([0 H]);
    xlabel('x [m]'); ylabel('y [m]');
    title(p.Results.Title);
    legend(legendHandles, 'Location', 'bestoutside');
    hold off;
end
