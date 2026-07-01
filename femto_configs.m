function [extras, descr, fig] = femto_configs(modeName, varargin)
%FEMTO_CONFIGS  Build a struct array of femtocell extras for the GA.
%
% Returns a struct array compatible with the `extras` argument of
% append_fixed_bs.m, ready to be assigned to:
%   - `extraBs` in optimize_base_station_ga.m, or
%   - `exp.extraBs` in ga_experiment.makeConfig / ga_experiment.train.
%
% USAGE
%   % Plain, no figure:
%   extras = femto_configs('coverage_fill_v1');
%
%   % With a preview figure (not saved to disk):
%   extras = femto_configs('coverage_fill_v1', 'showFigure', true);
%
%   % For a different scenario / map:
%   extras = femto_configs('coverage_fill_v1', 'expCode', '2-2-1');
%
%   % Capture description and figure handle as well:
%   [extras, descr, fig] = femto_configs('cell_edge_ring_v1', 'showFigure', true);
%
% INTEGRATION (the two intended use sites)
%
%   --- optimize_base_station_ga.m ---
%       extraBs = femto_configs('coverage_fill_v1');
%       % ... existing call to append_fixed_bs / precompute_mbs_power_maps
%
%   --- ga_experiment style ---
%       exp = ga_experiment.makeConfig('1-1-1');
%       exp.extraBs = femto_configs('coverage_fill_v1', 'expCode', '1-1-1');
%       runDir = ga_experiment.train(exp);
%
% INPUTS
%   modeName : char/string. One of (case-insensitive):
%                'coverage_fill_v1', 'coverage_fill_4kx3k_v1',
%                'residential_random_v1', 'cell_edge_ring_v1',
%                'per_mbs_edge_v1', 'hex_lattice_v1',
%                'hotspot_clusters_v1', 'none'.
%              Modes ending in '_4kx3k_v1' are tuned for the X-2-Z
%              scenario (4000x3000 map, 2 MBSs); use them with
%              expCode = '<x>-2-<z>'.
%
% NAME-VALUE OPTIONS
%   'showFigure'  (logical, default false)
%       Open a preview figure showing MBS + femto positions over the map.
%       The figure is shown only -- nothing is saved to disk.
%   'expCode'     (char, default '1-1-1')
%       X-Y-Z code passed to ga_experiment.makeConfig to resolve geometry
%       (W, H, margin, numMbs, ISD, xs, ys). Use this so the femto layout
%       respects whichever scenario you're targeting.
%   'fcCoverage'  (double, default 2e9)
%       Femto antenna coverage-band centre frequency [Hz].
%   'fcCapacity'  (double, default 2e9)
%       Femto antenna capacity-band centre frequency [Hz]. Set to 2.6e9
%       for a real dual-band split.
%   'numHotspots'      (integer, default 3)
%       Used by 'hotspot_clusters_v1': number of cluster centres.
%   'femtosPerHotspot' (integer, default 5)
%       Used by 'hotspot_clusters_v1': femtos per cluster. Total femto
%       count for that mode is numHotspots * femtosPerHotspot.
%   'femtosPerRing'    (integer, default 6)
%       Used by 'per_mbs_edge_v1': femtos placed around each MBS.
%   'ringRadius'       (double, default 700)
%       Used by 'per_mbs_edge_v1': cell-edge radius [m] for each ring.
%
% OUTPUTS
%   extras : 1×N struct array with fields x, y, height, power, antenna.
%            Empty struct array if N==0 (mode 'none').
%   descr  : short text label for logs/plots.
%   fig    : figure handle if 'showFigure' is true, otherwise [].
%
% HOW TO ADD A NEW MODE
%   Append a new `case 'my_mode_name'` block inside build_femto_extras()
%   and return the same (extras, descr) signature. Nothing else needs to
%   change.

    % --- Parse name-value options ---------------------------------------
    p = inputParser;
    p.addRequired('modeName', @(x) ischar(x) || isstring(x));
    p.addParameter('showFigure', false, @(x) islogical(x) || isnumeric(x));
    p.addParameter('expCode',    '1-1-1', @(x) ischar(x) || isstring(x));
    % Defaults sourced from band_frequencies so the femto antennas track
    % whatever the FBS/MBS antennas are using. Override per-call only when
    % you intentionally want a femto-only carrier (and remember that the
    % band-id contract assumes same physical band per band id).
    freqsDefault = band_frequencies();
    p.addParameter('fcCoverage', freqsDefault.coverage, @(x) isnumeric(x) && isscalar(x) && x > 0);
    p.addParameter('fcCapacity', freqsDefault.capacity, @(x) isnumeric(x) && isscalar(x) && x > 0);
    p.addParameter('numHotspots',      3, @(x) isnumeric(x) && isscalar(x) && x >= 1);
    p.addParameter('femtosPerHotspot', 5, @(x) isnumeric(x) && isscalar(x) && x >= 1);
    p.addParameter('femtosPerRing',    6, @(x) isnumeric(x) && isscalar(x) && x >= 1);
    p.addParameter('ringRadius',     700, @(x) isnumeric(x) && isscalar(x) && x > 0);
    p.parse(modeName, varargin{:});
    opt = p.Results;

    % --- Resolve experiment geometry ------------------------------------
    % Resolve only the geometry fields we need from the X-Y-Z code. We
    % deliberately avoid ga_experiment.makeConfig here: makeConfig itself
    % calls femto_configs to populate exp.extraBs, which would recurse
    % back into us. ga_experiment.scenarios is a leaf method.
    ga_experiment.ensurePaths();
    codeParts = strsplit(char(opt.expCode), '-');
    if numel(codeParts) ~= 3
        error('femto_configs:badExpCode', ...
              'expCode must be X-Y-Z (got %s).', char(opt.expCode));
    end
    scen = ga_experiment.scenarios(str2double(codeParts{2}));
    % margin/ISD mirror the constants in ga_experiment.makeConfig.
    exp = struct( ...
        'numMbs', scen.numMbs, 'xs', scen.xs, 'ys', scen.ys, ...
        'W',      scen.W,      'H',  scen.H, ...
        'margin', 100,         'ISD', 500);

    % --- Default femto antenna stack ------------------------------------
    defaultFemtoAnt = [local_setup_femto_antenna(opt.fcCoverage), ...
                       local_setup_femto_antenna(opt.fcCapacity)];

    % --- Build the chosen layout ----------------------------------------
    [extras, descr] = build_femto_extras(char(modeName), exp, defaultFemtoAnt, opt);

    % --- Optional preview figure ----------------------------------------
    fig = [];
    if logical(opt.showFigure)
        fig = preview_layout(exp, extras, descr, char(modeName));
    end
end


% =============================================================================
%                       FEMTO LAYOUT DEFINITIONS
% =============================================================================

function [extras, descr] = build_femto_extras(modeName, exp, defaultFemtoAnt, opt)
%BUILD_FEMTO_EXTRAS  Return struct array of fixed femtos for a named mode.
%   exp             : geometry struct (numMbs, xs, ys, W, H, margin, ISD)
%   defaultFemtoAnt : 1×nBands antenna stack used as fall-through
%   opt             : parsed name-value options from femto_configs (mode-
%                     specific knobs like numHotspots, femtosPerHotspot)

    switch lower(modeName)

        % -------------------------------------------------------------------
        case 'coverage_fill_v1'
        % -------------------------------------------------------------------
        % Hand-placed femtos targeting the visible MBS coverage holes on the
        % default 2000×1500 map with a single MBS near (1100, 950).
        %
        % Strategy:
        %   - Bottom strip  (y < 400) is the largest dead zone -> 6 femtos.
        %   - Far-left/right edges  -> 2-3 femtos.
        %   - Top corners (above MBS coverage radius) -> 3 femtos.
        %   - 1 in-fill femto at the bottom-right transition.
        %
        % Power:  0.25 W (24 dBm)  -- outdoor lamppost-style HeNB.
        % Height: 5 m              -- wall-mount / lamp post.
        % Antenna: omni at default band stack.
        %
        % Total: 13 femtos.

        descr = '13 outdoor femtos (5 m, 0.25 W) hand-placed in MBS dead zones';

        positions = [ ...
            % --- Bottom strip ---
             300,  200; ...   %  1: bottom-left corner
             700,  250; ...   %  2: bottom inner-left
            1100,  220; ...   %  3: bottom directly below MBS
            1500,  250; ...   %  4: bottom inner-right
            1800,  200; ...   %  5: bottom-right corner
             500,  450; ...   %  6: lower-left transition

            % --- Side edges (mid-height) ---
             150,  800; ...   %  7: far-left mid
            1900,  800; ...   %  8: far-right mid
             150, 1200; ...   %  9: far-left upper

            % --- Top edges & corners ---
             400, 1350; ...   % 10: top-left
            1850, 1350; ...   % 11: top-right
            1100, 1400; ...   % 12: top centre (above MBS, edge of coverage)

            % --- Far bottom-right in-fill ---
            1750,  500];      % 13: bottom-right in-fill

        N = size(positions, 1);
        extras = repmat(empty_extra(), 1, N);
        for k = 1:N
            extras(k).x       = positions(k, 1);
            extras(k).y       = positions(k, 2);
            extras(k).height  = 5;          % outdoor lamppost height
            extras(k).power   = 0.25;       % 24 dBm
            extras(k).antenna = defaultFemtoAnt;
        end

        % -------------------------------------------------------------------
        case 'coverage_fill_4kx3k_v1'
        % -------------------------------------------------------------------
        % Hand-placed dead-zone fill for the X-2-Z scenario:
        %   4000x3000 map with 2 MBSs at (1000, 1000) and (3500, 2200).
        % Targets the four map corners, the long bottom/top strips, and the
        % weak-coverage band running diagonally between the two MBSs.
        %
        % Power: 0.25 W (24 dBm), Height: 5 m, Antenna: defaultFemtoAnt.
        % Total: 18 femtos. Positions are clipped to the margin'd map so
        % calling this with a smaller scenario won't push BSs off-grid.

        descr = '18 outdoor femtos (5 m, 0.25 W) hand-placed for 4000x3000 / 2 MBS';

        positions = [ ...
            % --- Bottom strip (y < 400) ---
             400,  250; ...   %  1: bottom-left corner
            1000,  250; ...   %  2: bottom directly below MBS A
            1700,  300; ...   %  3: bottom mid-left
            2400,  300; ...   %  4: bottom mid-right
            3100,  300; ...   %  5: bottom directly below MBS B
            3800,  250; ...   %  6: bottom-right corner

            % --- Side edges (mid-height) ---
             200, 1500; ...   %  7: far-left mid
            3800, 1100; ...   %  8: far-right lower
             200, 2500; ...   %  9: far-left upper
            3800, 2800; ...   % 10: top-right corner

            % --- Top strip (y > 2700) ---
            3000, 2900; ...   % 11: top right of MBS B
            2200, 2800; ...   % 12: top mid-right
            1400, 2800; ...   % 13: top mid-left
             500, 2800; ...   % 14: top-left

            % --- Diagonal in-fill between MBSs ---
            2200, 1500; ...   % 15: centre fill
            2400, 1100; ...   % 16: between-MBS lower
            2700, 1700; ...   % 17: between-MBS upper
             600, 1900];      % 18: upper-left in-fill

        N = size(positions, 1);
        extras = repmat(empty_extra(), 1, N);
        for k = 1:N
            extras(k).x       = max(exp.margin, min(exp.W - exp.margin, positions(k, 1)));
            extras(k).y       = max(exp.margin, min(exp.H - exp.margin, positions(k, 2)));
            extras(k).height  = 5;
            extras(k).power   = 0.25;
            extras(k).antenna = defaultFemtoAnt;
        end

        % -------------------------------------------------------------------
        case 'residential_random_v1'
        % -------------------------------------------------------------------
        % Stub: uniform random scatter of indoor residential HeNBs.

        descr = '20 indoor residential femtos (3 m, 0.1 W) — uniform random';

        N    = 20;
        seed = 42;
        rng(seed);
        fx = exp.margin + (exp.W - 2*exp.margin) * rand(N, 1);
        fy = exp.margin + (exp.H - 2*exp.margin) * rand(N, 1);

        extras = repmat(empty_extra(), 1, N);
        for k = 1:N
            extras(k).x       = fx(k);
            extras(k).y       = fy(k);
            extras(k).height  = 3;          % residential
            extras(k).power   = 0.1;        % 20 dBm
            extras(k).antenna = defaultFemtoAnt;
        end

        % -------------------------------------------------------------------
        case 'cell_edge_ring_v1'
        % -------------------------------------------------------------------
        % Stub: ring of femtos along an iso-distance contour from the MBS
        % nominal position. Useful for studying cell-edge throughput uplift.

        descr = 'ring of 12 femtos (5 m, 0.25 W) at radius 600 m from MBS centre';

        % MBS nominal centre — pulled from the same hex generator as ga_experiment
        if exp.numMbs == 1
            [xs, ys] = generate_hex_sites(exp.W, exp.H, exp.ISD, exp.margin, exp.numMbs);
        else
            xs = exp.xs; ys = exp.ys;
        end
        mbs_x_ctr = mean(xs);
        mbs_y_ctr = mean(ys);

        N      = 12;
        radius = 600;
        theta  = linspace(0, 2*pi, N+1); theta(end) = [];
        fx     = mbs_x_ctr + radius * cos(theta);
        fy     = mbs_y_ctr + radius * sin(theta);

        % Clamp inside the map (keeps things sane near edges)
        fx = max(exp.margin, min(exp.W - exp.margin, fx));
        fy = max(exp.margin, min(exp.H - exp.margin, fy));

        extras = repmat(empty_extra(), 1, N);
        for k = 1:N
            extras(k).x       = fx(k);
            extras(k).y       = fy(k);
            extras(k).height  = 5;
            extras(k).power   = 0.25;
            extras(k).antenna = defaultFemtoAnt;
        end

        % -------------------------------------------------------------------
        case 'per_mbs_edge_v1'
        % -------------------------------------------------------------------
        % Cell-edge ring around EACH MBS (vs. cell_edge_ring_v1 which uses
        % a single ring around the centroid). Generic across scenarios but
        % especially useful for X-2-Z layouts where two MBSs sit far apart
        % and a centroid-anchored ring lands in no-man's-land.
        %
        % Knobs:
        %   'femtosPerRing' (default 6)
        %   'ringRadius'    (default 700)  -- cell-edge radius [m]
        % Total femtos = numMbs * femtosPerRing.
        %
        % Power: 0.25 W, Height: 5 m. Each femto is clamped inside the
        % margin'd map.

            N      = round(opt.femtosPerRing);
            R      = opt.ringRadius;

            if exp.numMbs == 1
                [mxs, mys] = generate_hex_sites(exp.W, exp.H, exp.ISD, exp.margin, exp.numMbs);
            else
                mxs = exp.xs; mys = exp.ys;
            end

            descr = sprintf('%d MBS rings × %d femtos @ R=%.0f m (5 m, 0.25 W)', ...
                            numel(mxs), N, R);

            theta = linspace(0, 2*pi, N+1); theta(end) = [];
            extras = repmat(empty_extra(), 1, numel(mxs) * N);
            idx = 1;
            for m = 1:numel(mxs)
                for j = 1:N
                    fx = mxs(m) + R * cos(theta(j));
                    fy = mys(m) + R * sin(theta(j));
                    extras(idx).x       = max(exp.margin, min(exp.W - exp.margin, fx));
                    extras(idx).y       = max(exp.margin, min(exp.H - exp.margin, fy));
                    extras(idx).height  = 5;
                    extras(idx).power   = 0.25;
                    extras(idx).antenna = defaultFemtoAnt;
                    idx = idx + 1;
                end
            end

        % -------------------------------------------------------------------
        case 'hex_lattice_v1'
        % -------------------------------------------------------------------
        % Stub: regular hex lattice of femtos with chosen ISD.

        descr = 'hex lattice of femtos at ISD=400 m (3 m, 0.1 W)';

        femtoISD    = 400;
        femtoMargin = 200;
        % Reuse the existing hex generator. It errors if N exceeds capacity,
        % so bisect down until it fits.
        fx = []; fy = [];
        for cap = 50:-1:1
            try
                [fx, fy] = generate_hex_sites(exp.W, exp.H, femtoISD, femtoMargin, cap);
                break;
            catch
            end
        end

        N = numel(fx);
        extras = repmat(empty_extra(), 1, N);
        for k = 1:N
            extras(k).x       = fx(k);
            extras(k).y       = fy(k);
            extras(k).height  = 3;
            extras(k).power   = 0.1;
            extras(k).antenna = defaultFemtoAnt;
        end

        % -------------------------------------------------------------------
        case 'hotspot_clusters_v1'
        % -------------------------------------------------------------------
        % Outdoor pico/femto hotspots — multiple dense clusters of femtos
        % scattered across the map (think shopping plazas, transit hubs,
        % event venues). Cluster centres are placed deterministically with
        % a fixed rng seed and rejected if they fall too close to any MBS.
        %
        % Knobs (passed via name-value to femto_configs):
        %   'numHotspots'      (default 3)
        %   'femtosPerHotspot' (default 5)
        % Total femtos = numHotspots * femtosPerHotspot.
        %
        % Power: 0.25 W (24 dBm), Height: 5 m, Antenna: defaultFemtoAnt.

            K = round(opt.numHotspots);
            N = round(opt.femtosPerHotspot);
            descr = sprintf('%d hotspot cluster(s) × %d femtos (5 m, 0.25 W)', K, N);

            if exp.numMbs == 1
                [mxs, mys] = generate_hex_sites(exp.W, exp.H, exp.ISD, exp.margin, exp.numMbs);
            else
                mxs = exp.xs; mys = exp.ys;
            end

            rngState = rng;
            cleanupRng = onCleanup(@() rng(rngState));
            rng(7);   % deterministic but distinct from residential_random_v1

            minMbsDist = max(300, 0.25 * min(exp.W, exp.H));
            centres    = zeros(K, 2);
            placed     = 0;
            tries      = 0;
            while placed < K && tries < 5000
                tries = tries + 1;
                cx = exp.margin + (exp.W - 2*exp.margin) * rand;
                cy = exp.margin + (exp.H - 2*exp.margin) * rand;
                if all(hypot(mxs(:) - cx, mys(:) - cy) > minMbsDist)
                    placed = placed + 1;
                    centres(placed, :) = [cx, cy];
                end
            end
            % Fallback (small map / many MBS): drop the MBS-distance constraint
            for k = (placed + 1):K
                centres(k, :) = [exp.margin + (exp.W - 2*exp.margin) * rand, ...
                                 exp.margin + (exp.H - 2*exp.margin) * rand];
            end

            jitterRadius = 60;
            extras = repmat(empty_extra(), 1, K * N);
            idx = 1;
            for k = 1:K
                cx = centres(k, 1); cy = centres(k, 2);
                for j = 1:N
                    theta = 2*pi*rand;
                    r     = jitterRadius * sqrt(rand);
                    fx = max(exp.margin, min(exp.W - exp.margin, cx + r*cos(theta)));
                    fy = max(exp.margin, min(exp.H - exp.margin, cy + r*sin(theta)));
                    extras(idx).x       = fx;
                    extras(idx).y       = fy;
                    extras(idx).height  = 5;
                    extras(idx).power   = 0.25;
                    extras(idx).antenna = defaultFemtoAnt;
                    idx = idx + 1;
                end
            end

        % -------------------------------------------------------------------
        case 'none'
        % -------------------------------------------------------------------
        % Sanity baseline: no femtos at all. Useful for A/B comparisons.
            descr  = 'no femtos (baseline)';
            extras = empty_extra(); extras(1) = []; % return 1×0 struct array

        % -------------------------------------------------------------------
        otherwise
        % -------------------------------------------------------------------
            error('femto_configs:unknownMode', ...
                  'Unknown modeName = "%s". Valid: coverage_fill_v1, coverage_fill_4kx3k_v1, residential_random_v1, cell_edge_ring_v1, per_mbs_edge_v1, hex_lattice_v1, hotspot_clusters_v1, none.', ...
                  modeName);
    end
end


% =============================================================================
%                              HELPERS
% =============================================================================

function s = empty_extra()
% Default schema dictated by append_fixed_bs.m
    s = struct('x', 0, 'y', 0, 'height', 3, 'power', 0.1, 'antenna', []);
end

function l = local_setup_femto_antenna(centerFreq)
% Femto-flavoured equivalent of setup_antenna(): omnidirectional 0 dBi.
    if nargin < 1 || isempty(centerFreq), centerFreq = 2e9; end
    s = qd_simulation_parameters;
    s.center_frequency = centerFreq;
    l = qd_layout(s);
    l.no_tx = 1;
    l.tx_array(1,1) = qd_arrayant('omni');
    l.rx_array      = qd_arrayant('omni');
end

function fig = preview_layout(exp, extras, descr, modeName)
% Lightweight preview: MBS positions, femto positions, map bounds.

    fig = figure('Name', sprintf('Femto preview — %s', modeName), ...
                 'NumberTitle', 'off', 'Position', [120 120 900 600]);
    ax = axes('Parent', fig); hold(ax, 'on'); grid(ax, 'on');
    set(ax, 'Box', 'on', 'Layer', 'top');
    xlim(ax, [0 exp.W]); ylim(ax, [0 exp.H]);
    xlabel(ax, 'x [m]'); ylabel(ax, 'y [m]');
    title(ax, sprintf('Femto deployment preview: %s', descr), ...
          'Interpreter', 'none', 'FontSize', 11);

    % Map margin rectangle
    mrec = exp.margin;
    plot(ax, [mrec exp.W-mrec exp.W-mrec mrec mrec], ...
              [mrec mrec exp.H-mrec exp.H-mrec mrec], ...
              ':', 'Color', [0.5 0.5 0.5], 'DisplayName', 'margin');

    % MBS positions
    if exp.numMbs == 1
        [xs, ys] = generate_hex_sites(exp.W, exp.H, exp.ISD, exp.margin, exp.numMbs);
    else
        xs = exp.xs; ys = exp.ys;
    end
    plot(ax, xs, ys, 'p', 'MarkerSize', 18, 'MarkerFaceColor', [1 0.85 0.15], ...
         'MarkerEdgeColor', 'k', 'LineWidth', 1.2, ...
         'DisplayName', sprintf('Macro BS (n=%d)', exp.numMbs));

    % Femto positions
    if ~isempty(extras)
        fx = [extras.x];
        fy = [extras.y];
        scatter(ax, fx, fy, 90, [0.20 0.55 0.85], 'filled', ...
                'MarkerEdgeColor', 'k', 'LineWidth', 0.6, ...
                'DisplayName', sprintf('Femtos (n=%d)', numel(extras)));

        % Coverage radius hint per femto (rough, cosmetic only)
        for k = 1:numel(extras)
            r = femto_radius_hint_m(extras(k).power);
            theta = linspace(0, 2*pi, 60);
            plot(ax, extras(k).x + r*cos(theta), extras(k).y + r*sin(theta), ...
                 '-', 'Color', [0.20 0.55 0.85 0.25], 'HandleVisibility', 'off');
        end
    end

    legend(ax, 'Location', 'eastoutside');
    axis(ax, 'equal');
    xlim(ax, [0 exp.W]); ylim(ax, [0 exp.H]);
end

function r = femto_radius_hint_m(power_W)
% Very rough free-space breakeven radius for visualization only.
% NOT used by SINREvaluation; just shades a cosmetic disk on the preview.
    if power_W <= 0.05      % <= 17 dBm
        r = 25;
    elseif power_W <= 0.12  % residential ~20 dBm
        r = 40;
    elseif power_W <= 0.30  % outdoor ~24 dBm
        r = 70;
    else
        r = 100;
    end
end