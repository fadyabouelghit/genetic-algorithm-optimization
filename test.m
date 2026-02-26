

% no_fbs = 2;
% tx_x = [1317.2, 1236.9];
% tx_y = [115.85, 1475.5];
% tx_height = [24.79, 37.236];
% tx_power = [9.41, 8.546];
% power_status = ones(1,no_fbs);

    % 
    % BS1 X (m)           72.579
    % BS1 Y (m)           356.31
    % BS1 Z (m)           56.741
    % BS1 Power (W)         9.62
    % BS1 Power Status         1
    % BS2 X (m)           77.406
    % BS2 Y (m)           283.79
    % BS2 Z (m)           24.139
    % BS2 Power (W)       9.9492
    % BS2 Power Status         1

no_fbs = 1;
tx_x = 198.73;
tx_y = 192.09;
tx_height = 21.899;
tx_power = 9.1958;   
power_status = 1;

% no_fbs = 2;
% tx_x = [72.579 77.406];
% tx_y = [356.31 283.79];
% tx_height = [56.741 24.139];
% tx_power = [9.62 9.949];   
% power_status = [1 1];


fbsAntenna = setup_antenna();
mbsAntenna = setup_antenna();

numMbs = 1;
containsMbs = 1;
W = 2000; H = 1500;      
margin = 100;            
ISD = 500;               
[xs, ys] = generate_hex_sites(W, H, ISD, margin, numMbs);
mbs_height   = 25;
mbs_power    = 20;
[mbs_params, antennaObjectMbs, containsMbs, numMbs] = pack_mbs_params...
    (xs, ys, mbs_height, mbs_power, mbsAntenna);

tempForX = mbs_params(1,:);
mbs_params(1,:) = mbs_params(2,:);
mbs_params(2,:) = tempForX;
% xs 
% ys
% mbs_params(1,:) = [750, 750];
% mbs_params(2,:) = [750, 1500];
cache = precompute_mbs_power_maps( ...
    mbs_params, ...
    antennaObjectMbs, ...
    subset, ...
    '3GPP_38.901_UMa_LOS', ...   % scenario
    'quick', ...                  % mode
    1.5, ...                      % UE height
    './cache_mbs_maps' ...        % folder for on-disk cache
);
% mbs_x = [300];       
% mbs_y = [350];     
% mbs_height = [25];        
% mbs_power = [20];
% mbs_params = [mbs_x; mbs_y; mbs_height; mbs_power]; 


num_users = 1000;
threshold = 5; % in dB

subset_x_min = 0;
subset_x_max = W;
subset_y_min = 0;
subset_y_max = H;

% FBS antenna object 
l = setup_antenna();


containsMbs = 1;
[user_positions, df_users_table, total_connected_users, total_transmitted_pwr, avg_rate_connected_bpsHz, ~, ~] = SINREvaluation(...,
                                                                                                l,power_status, tx_x, tx_y, tx_height, no_fbs, tx_power, ...,
                                                                                                mbs_params(2,:), mbs_params(1,:), mbs_height, mbs_power,...,
                                                                                                subset_x_min, subset_x_max, subset_y_min, subset_y_max, num_users, ..., 
                                                                                                threshold,containsMbs,antennaObjectMbs,cache);
 
% [user_positions, df_users_table, total_connected_users, total_transmitted_pwr] = SINREvaluation_noCache(...,
%                                                                                                 l,power_status, tx_x, tx_y, tx_height, no_fbs, tx_power, ...,
%                                                                                                 xs, ys, mbs_height, mbs_power,...,
%                                                                                                 subset_x_min, subset_x_max, subset_y_min, subset_y_max, num_users, ..., 
%                                                                                                 threshold,containsMbs,antennaObjectMbs);

x = user_positions(:, 1); % Extract x-coordinates
y = user_positions(:, 2); % Extract y-coordinates
% data = [x(:), y(:)]; % Ensures column vectors and proper concatenation
% k = 1; 
% [idx, centroids] = kmeans(data, k);


% [user_positions, df_users_table, total_connected_users, total_transmitted_pwr] = SINREvaluation(...,
%                                                                                                 l,power_status, centroids(2), centroids(1), tx_height, no_fbs, tx_power, ...,
%                                                                                                 mbs_x, mbs_y, mbs_height, mbs_power,...,
%                                                                                                 subset_x_min, subset_x_max, subset_y_min, subset_y_max, num_users, ..., 
%                                                                                                 threshold,containsMbs,[l2,l2]);


fprintf("total connected users: %.2f     total power: %.2f \n", total_connected_users, total_transmitted_pwr*power_status)






connected_indices = df_users_table.is_connected == 1;
unconnected_indices = df_users_table.is_connected == 0;

% Create the scatter plot
figure;
hold on;
scatter(x(connected_indices), y(connected_indices), 'b', 'DisplayName', 'Connected');
scatter(x(unconnected_indices), y(unconnected_indices), 'r', 'DisplayName', 'Unconnected');
scatter(mbs_params(2,:),mbs_params(1,:),'marker','x','MarkerEdgeColor','black', 'MarkerFaceColor','black','LineWidth',3);
xlabel('X Coordinate');
ylabel('Y Coordinate');
title('Connected and Unconnected Users');
legend;
xlim([-500 W+500])
ylim([-500 H+500])
hold off;

fbs1_indices = (df_users_table.FBS_connection_index == 1) & connected_indices;
fbs2_indices = (df_users_table.FBS_connection_index == 2) & connected_indices;
fbs3_indices = (df_users_table.FBS_connection_index == 3) & connected_indices;

% Create the scatter plot
% figure;
% hold on;
% scatter(x(fbs1_indices), y(fbs1_indices), 'g', 'DisplayName', 'Connected to FBS 1');
% scatter(x(fbs2_indices), y(fbs2_indices), 'b', 'DisplayName', 'Connected to FBS 2');
% scatter(x(fbs3_indices), y(fbs3_indices), 'm', 'DisplayName', 'Connected to FBS 3');
% scatter(x(unconnected_indices), y(unconnected_indices), 'r', 'DisplayName', 'Unconnected');
% scatter(xs,ys, 'marker','x','MarkerEdgeColor','black', 'MarkerFaceColor','black');
% xlabel('X Coordinate');
% ylabel('Y Coordinate');
% title('User Connections by FBS');
% xlim([-500 W+500])
% ylim([-500 H+500])
% legend;
% hold off;

