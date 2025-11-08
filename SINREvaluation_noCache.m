function [user_positions, df_users_table, total_connected_users, total_transmitted_pwr] = SINREvaluation_noCache(antenna_object, power_status, tx_x, tx_y, tx_height, no_fbs, tx_power, mbs_x, mbs_y, mbs_height, mbs_power, subset_x_min, subset_x_max, subset_y_min, subset_y_max, num_users, threshold, containsMbs, antennaObjectMbs)
% SINREVALUATION Computes the SINR values for users based on FBS and optional MBS parameters.
%
%   INPUTS:
%     antenna_object     : FBS antenna object (Quadriga-compatible)
%     power_status       : Activation vector for each FBS (1 = active, 0 = inactive)
%     tx_x, tx_y         : X and Y coordinates of each FBS
%     tx_height          : Heights of each FBS
%     no_fbs             : Number of Flying Base Stations (FBSs)
%     tx_power           : Transmission power for each FBS
%     mbs_x, mbs_y       : Coordinates of MBSs (optional)
%     mbs_height         : Heights of MBSs (optional)
%     mbs_power          : Transmission power for each MBS
%     subset_x_min/max,
%     subset_y_min/max   : Bounds for the region of interest
%     num_users          : Number of users
%     threshold          : SINR threshold in dB
%     containsMbs        : Binary flag indicating MBS presence
%     antennaObjectMbs   : Array of MBS antenna objects
%
%   OUTPUTS:
%     user_positions         : Matrix of user coordinates
%     df_users_table         : Table with power values and SINRs per BS
%     total_connected_users  : Total number of connected users (SINR ≥ threshold)
%     total_transmitted_pwr  : Sum of transmission power of all active FBSs

    user_positions = generate_user_positions(subset_x_min, subset_x_max, subset_y_min, subset_y_max, num_users);

    [~, df_users_table] = calculate_power_iterator(antenna_object, no_fbs, power_status, tx_y, tx_x, tx_height, tx_power, ...
        mbs_x, mbs_y, mbs_height, mbs_power, user_positions, ...
        subset_x_min, subset_x_max, subset_y_min, subset_y_max, threshold, containsMbs, antennaObjectMbs);

    total_connected_users = sum(df_users_table.is_connected);
    total_transmitted_pwr = sum(tx_power);
end

function user_positions = generate_user_positions(x_min, x_max, y_min, y_max, num_users)
% GENERATE_USER_POSITIONS Generates uniformly distributed user positions.
%
%   INPUTS:
%     x_min, x_max     : Range of x-coordinates
%     y_min, y_max     : Range of y-coordinates
%     num_users        : Number of users
%
%   OUTPUT:
%     user_positions   : num_users × 2 matrix of (x, y) user positions

    rng(0);  % For reproducibility

    % Adjust zero bounds for indexing
    x_min = max(1, x_min);
    y_min = max(1, y_min);

    user_positions = [randi([x_min, x_max], num_users, 1), ...
                      randi([y_min, y_max], num_users, 1)];
end

function user_positions = generate_user_positions_clustered(x_min, x_max, y_min, y_max, num_users, separation, std1, std2)
% GENERATE_USER_POSITIONS_CLUSTERED Creates two user clusters within defined bounds.
%
%   INPUTS:
%     x_min, x_max, y_min, y_max : Bounds of area
%     num_users                  : Number of users
%     separation, std1, std2     : Cluster separation and standard deviations
%
%   OUTPUT:
%     user_positions             : Clustered user positions matrix

    rng(0);

    adjusted_x_min = x_min + 1;
    adjusted_x_max = x_max - 1;
    adjusted_y_min = y_min + 1;
    adjusted_y_max = y_max - 1;

    mean1 = [adjusted_x_min + separation/2, adjusted_y_min + separation/2];
    mean2 = [adjusted_x_max - separation/2, adjusted_y_max - separation/2];

    cluster1 = mean1 + std1 * randn(num_users/2, 2);
    cluster2 = mean2 + std2 * randn(num_users/2, 2);

    user_positions = round([cluster1; cluster2]);
    user_positions(:, 1) = max(min(user_positions(:, 1), adjusted_x_max), adjusted_x_min);
    user_positions(:, 2) = max(min(user_positions(:, 2), adjusted_y_max), adjusted_y_min);
end

function [power_map, x_coords, y_coords] = calculate_power(antenna_object, tx_x, tx_y, tx_height, tx_power, subset_x_min, subset_x_max, subset_y_min, subset_y_max)
% CALCULATE_POWER Computes received power map for a transmitter using Quadriga.
%
%   OUTPUTS:
%     power_map : Grid of received power values (mW)
%     x_coords, y_coords : Coordinates associated with power_map

    antenna_object.tx_position(:,1) = [tx_x; tx_y; tx_height];

    [map, x_coords, y_coords] = antenna_object.power_map('3GPP_38.901_UMa_LOS', 'quick', 1, ...
        subset_x_min, subset_x_max, subset_y_min, subset_y_max, 1.5, tx_power);

    power_map = sum(cat(3, map{:}), 3); % Total received power in mW
end

function [df_users, df_users_table] = calculate_power_iterator(antenna_object, no_fbs, power_status, tx_y, tx_x, tx_height, tx_power, ...
    mbs_x, mbs_y, mbs_height, mbs_power, user_positions, ...
    subset_x_min, subset_x_max, subset_y_min, subset_y_max, ... 
    threshold, containsMbs, antennaObjectMbs, mbsCache)
% CALCULATE_POWER_ITERATOR Evaluates SINR at each user for all BSs.
%
%   OUTPUTS:
%     df_users       : Raw matrix of power received at each user from each BS
%     df_users_table : Structured table of power and SINR values

    num_users = size(user_positions, 1);
    df_users = zeros(num_users, no_fbs);
    noise_power = 1e-11;

    % FBS loop
    for fbs_id = 1:no_fbs
        if power_status(fbs_id)
            [P,x_coords,y_coords] = calculate_power(antenna_object, tx_x(fbs_id), tx_y(fbs_id), tx_height(fbs_id), ...
                                tx_power(fbs_id), subset_x_min, subset_x_max, subset_y_min, subset_y_max);
        else
            P = zeros(subset_x_max + 1, subset_y_max + 1);
        end
        
        % idx = sub2ind(size(P), user_positions(:, 2), user_positions(:, 1));
        % df_users(:, fbs_id) = P(idx);
        df_users(:,fbs_id) = sample_nearest(P, user_positions);

    end

    % MBS loop
    if containsMbs
        num_mbs = numel(antennaObjectMbs);
        df_users(:, no_fbs+1:no_fbs+num_mbs) = 0;

        for mbs_idx = 1:num_mbs
            P_mbs = calculate_power(antennaObjectMbs(mbs_idx), ...
                mbs_x(mbs_idx), mbs_y(mbs_idx), mbs_height, mbs_power, ...
                subset_x_min, subset_x_max, subset_y_min, subset_y_max);

            % idx = sub2ind(size(P_mbs), user_positions(:, 1), user_positions(:, 2));
            % df_users(:, no_fbs + mbs_idx) = P_mbs(idx);
            df_users(:,no_fbs + mbs_idx) = sample_nearest(P_mbs, user_positions);
        end
    end

    df_users_table = convert_to_dataframe_style(df_users);

    % SINR evaluation
    max_sinr = -inf(num_users, 1);
    max_fbs_index = zeros(num_users, 1);

    bsIterator = size(df_users, 2);

    for bs_id = 1:bsIterator
        sinr_column = zeros(num_users, 1);
        parfor user_id = 1:num_users
            signal = df_users(user_id, bs_id);
            interference = sum(df_users(user_id, :)) - signal;
            sinr = signal / (interference + noise_power);
            sinr_db = 10 * log10(sinr);
            sinr_column(user_id) = sinr_db;

            eta = log2(1 + sinr);

            if sinr_db >= threshold && sinr_db > max_sinr(user_id)
                max_sinr(user_id) = sinr_db;
                max_fbs_index(user_id) = bs_id;
            end
        end
        df_users_table.(['SINR_FBS_' num2str(bs_id)]) = sinr_column;
    end

    df_users_table.is_connected = max_sinr >= threshold;
    df_users_table.FBS_connection_index = max_fbs_index;
end

function df_users_table = convert_to_dataframe_style(df_users)
% CONVERT_TO_DATAFRAME_STYLE Converts power matrix into a table with labeled columns.
%
%   INPUT:
%     df_users         : num_users × num_base_stations matrix
%   OUTPUT:
%     df_users_table   : MATLAB table with named columns for each FBS

    [~, no_fbs] = size(df_users);
    col_names = arrayfun(@(i) sprintf('power_FBS_%d', i), 1:no_fbs, 'UniformOutput', false);
    df_users_table = array2table(df_users, 'VariableNames', col_names);
end

function user_powers = sample_nearest(P, user_positions)
    % P: (ny x nx), ny = 3001, nx = 3001, covering x,y in [0,3000]
    % user_positions: (N x 2) [x y], within [0,3000]
    % Returns: (N x 1) power values from P at nearest grid point

    ix = round(user_positions(:,1)) + 1;  % columns (x)
    iy = round(user_positions(:,2)) + 1;  % rows (y)

    % Clamp to valid indices (handles rare boundary rounding)
    [nx, ny] = size(P);
    ix = min(max(ix, 1), nx);
    iy = min(max(iy, 1), ny);

    % Linear indexing
    user_powers = P(sub2ind([nx, ny], ix, iy));
end
