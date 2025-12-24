function cache = precompute_mbs_power_maps(mbs_params, antennaObjectMbs, subset, scenario, mode, ue_height, cacheDir, returnCell)
% Precompute (or load) per-MBS power maps and store them in a cache struct.
% Inputs:
%   mbs_params        : 4×N [x; y; z; power]
%   antennaObjectMbs  : 1×N array of antenna objects (template cloned per MBS)
%   subset            : struct with fields xmin, xmax, ymin, ymax
%   scenario          : e.g., '3GPP_38.901_UMa_LOS'
%   mode              : e.g., 'quick'
%   ue_height         : e.g., 1.5
%   cacheDir          : folder to store .mat cache files (optional but recommended)
%
% Output:
%   cache(i) : struct with fields
%       .key, .x, .y, .z, .power, .scenario, .subset (struct),
%       .ue_height, .mode, .x_coords, .y_coords, .map (single)
%
% Notes:
%   - Ensures all sites share identical x_coords/y_coords; warns if not.
%   - Uses on-disk caching keyed by MD5 of the config to avoid recomputing.

    arguments
        mbs_params (4,:) double
        antennaObjectMbs (1,:)    % pass your replicated antenna objects
        subset struct
        scenario char
        mode char
        ue_height double = 1.5
        cacheDir char = ''
        returnCell logical = false
    end

    N = size(mbs_params, 2);
    if ~isempty(cacheDir) && ~exist(cacheDir,'dir')
        mkdir(cacheDir);
    end
    cache = repmat(struct('key',[], 'x',[], 'y',[], 'z',[], 'power',[], ...
                          'scenario',[], 'subset',[], 'ue_height',[], 'mode',[], ...
                          'x_coords',[], 'y_coords',[], 'map',[]), 1, N);

    commonGrid = [];  % to verify identical coordinate grids across sites

    for i = 1:N
        x  = mbs_params(1,i);
        y  = mbs_params(2,i);
        z  = mbs_params(3,i);
        pw = mbs_params(4,i);

        key = make_mbs_key(x,y,z,pw,scenario,mode,ue_height,subset);

        map_file = '';
        if ~isempty(cacheDir)
            map_file = fullfile(cacheDir, [key '.mat']);
        end

        if ~isempty(map_file) && exist(map_file,'file')
            S = load(map_file, 'map', 'x_coords', 'y_coords', 'meta');
            map       = S.map;
            x_coords  = S.x_coords;
            y_coords  = S.y_coords;
        else
            % Compute once
            ant = antennaObjectMbs(i);
            ant.tx_position(:,1) = [x; y; z];

            % Your original call
            [cellMaps, x_coords, y_coords] = ant.power_map( ...
                scenario, mode, 1, ...
                subset.xmin, subset.xmax, subset.ymin, subset.ymax, ...
                ue_height, pw);

            map = sum(cat(3, cellMaps{:}), 3);
            map = single(map);  % memory saver, usually sufficient

            % Save to disk cache
            if ~isempty(map_file)
                meta = struct('x',x,'y',y,'z',z,'power',pw, ...
                              'scenario',scenario,'mode',mode, ...
                              'ue_height',ue_height,'subset',subset);
                try
                    save(map_file, 'map', 'x_coords', 'y_coords', 'meta', '-v7.3');
                catch
                    % If -v7.3 fails for any reason, fallback
                    save(map_file, 'map', 'x_coords', 'y_coords', 'meta');
                end
            end
        end

        % Fill cache entry
        cache(i).key       = key;
        cache(i).x         = x;
        cache(i).y         = y;
        cache(i).z         = z;
        cache(i).power     = pw;
        cache(i).scenario  = scenario;
        cache(i).subset    = subset;
        cache(i).ue_height = ue_height;
        cache(i).mode      = mode;
        cache(i).x_coords  = x_coords;
        cache(i).y_coords  = y_coords;
        cache(i).map       = map;

        % Consistency check for grid alignment
        thisGrid = [numel(x_coords), numel(y_coords), x_coords(1), x_coords(end), y_coords(1), y_coords(end)];
        if isempty(commonGrid)
            commonGrid = thisGrid;
        else
            if any(abs(thisGrid - commonGrid) > 1e-9)
                warning('MBS #%d grid differs from others. Power maps may not align exactly.', i);
            end
        end
    end
    if returnCell
        cache = arrayfun(@(s) s, cache, 'UniformOutput', false);
    end
end

function [map, x_coords, y_coords] = lookup_mbs_power_map(cache, i)
% Retrieve cached map for MBS index i
    map       = cache(i).map;
    x_coords  = cache(i).x_coords;
    y_coords  = cache(i).y_coords;
end

function key = make_mbs_key(x,y,z,pw,scenario,mode,ue_height,subset)
% Deterministic key from full configuration (robust to float formatting)
    payload = struct('x',x,'y',y,'z',z,'pw',pw, ...
                     'scenario',string(scenario), 'mode',string(mode), ...
                     'ue',ue_height, ...
                     'xmin',subset.xmin,'xmax',subset.xmax,'ymin',subset.ymin,'ymax',subset.ymax);
    s = jsonencode(payload);  % stable-ish textual representation
    key = md5_of_string(s);
end

function h = md5_of_string(s)
% MD5 as hex using Java (available in MATLAB)
    md = java.security.MessageDigest.getInstance('MD5');
    md.update(uint8(s));
    d = typecast(md.digest, 'uint8');
    h = lower(reshape(dec2hex(d)',1,[]));
end
