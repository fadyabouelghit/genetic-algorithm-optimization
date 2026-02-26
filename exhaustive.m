
no_fbs = 1;
num_users = 100;
threshold = 5; 

tx_x = 1140;
tx_y = 550;
tx_height = 24;
tx_power = 10.5;   
power_status = 1;

mbs_x = [300];       
mbs_y = [350];     
mbs_height = [25];        
mbs_power = [20];
mbs_params = [mbs_x; mbs_y; mbs_height; mbs_power]; 


subset_x_min = 0;
subset_x_max = 1500;
subset_y_min = 0;
subset_y_max = 1500;

% FBS antenna object 
l = setup_antenna();
l2 = l;

containsMbs = 1;
[user_positions, df_users_table, total_connected_users, total_transmitted_pwr, ~, ~, ~] =  SINREvaluation(...,
                                                                                                l,power_status, tx_x, tx_y, tx_height, no_fbs, tx_power, ...,
                                                                                                mbs_x, mbs_y, mbs_height, mbs_power,...,
                                                                                                subset_x_min, subset_x_max, subset_y_min, subset_y_max, num_users, ..., 
                                                                                                threshold,containsMbs,[l2]);

% fprintf("total connected users: %.2f     total power: %.2f \n", total_connected_users, total_transmitted_pwr)

x = user_positions(:, 1); % Extract x-coordinates
y = user_positions(:, 2); % Extract y-coordinates




%% Exhaustive loop 1 FBS

delta = 10;
x_min = 300; x_max = 1500;
y_min = 300; y_max = 1500;
% height_values = [30, 40, 60, 70];
height_values = [24];

% Generate grid vectors
x_vals = x_min:delta:x_max;
y_vals = y_min:delta:y_max;

position_list = [];

for xi = 1:1:length(x_vals)
    for yi = 1:1:length(y_vals)
        x_grid = x_vals(xi);
        y_grid = y_vals(yi);
        
        for h = height_values
            position_list = [position_list; x_grid, y_grid, h];
        end
    end
end


fprintf("size of generated list: %d \n", size(position_list,1))
fprintf("min of x = %.2f // max of x = %.2f \n", min(position_list(:,1)), max(position_list(:,1)))
fprintf("min of y = %.2f // max of y = %.2f \n", min(position_list(:,2)), max(position_list(:,2)))

n_points = size(position_list, 1);
progress = 0;

if isempty(gcp('nocreate'))
    % parpool;
    parpool('local', 10);  

end

D = parallel.pool.DataQueue;
n_points = size(position_list, 1);
progress = 0;
update_interval = 100;
afterEach(D, @(x) updateProgress(x, n_points, update_interval));

results_matrix = zeros(n_points, 4);  

parfor i = 1:n_points
    x_grid = position_list(i, 1);
    y_grid = position_list(i, 2);
    h = position_list(i, 3);

    [~, ~, connected_users, ~, ~, ~, ~] = SINREvaluation(...,
        l, power_status, x_grid, y_grid, h, no_fbs, tx_power, ...
        mbs_x, mbs_y, mbs_height, mbs_power, ...
        subset_x_min, subset_x_max, subset_y_min, subset_y_max, num_users, ...
        threshold, containsMbs, [l2]);

    results_matrix(i, :) = [x_grid, y_grid, h, connected_users];
    send(D, i); % Send update

end

results_table = array2table(results_matrix, ...
    'VariableNames', {'X', 'Y', 'Height', 'ConnectedUsers'});

delete(gcp('nocreate'));

save('exhaustive_granular_new.mat', 'results_table')
[maxVal, idx] = max(results_table.ConnectedUsers);

% 2D slices on a 3D plane at different heights 
scatter3(results_table.X, results_table.Y, results_table.Height, ...
         50, results_table.ConnectedUsers, 'filled');
xlabel('X');
ylabel('Y');
zlabel('Height');
title('Connected Users by 3D Position');
colorbar;
view(45,30);

% box plot per height entry for all locations aggregated 
grouped = groupsummary(results_table, 'Height', 'mean', 'ConnectedUsers');
plot(grouped.Height, grouped.mean_ConnectedUsers, '-o');
xlabel('Height');
ylabel('Mean Connected Users');
title('Average Users vs. Height');
grid on;
boxplot(results_table.ConnectedUsers, results_table.Height);
xlabel('Height');
ylabel('Connected Users');
title('Distribution of Connected Users by Height');
grid on;

heights = height_values;
for i = 1:length(heights)
    h = heights(i);
    subset = results_table(results_table.Height == h, :);
    [Xq, Yq] = meshgrid(unique(subset.X), unique(subset.Y));
    Z = griddata(subset.X, subset.Y, subset.ConnectedUsers, Xq, Yq);

    figure;
    contourf(Xq, Yq, Z, 20);
    title(['Performance Map at Height = ', num2str(h)]);
    xlabel('X'); ylabel('Y'); colorbar;
end


fprintf('Best results are: ')
fprintf("X = %.3f    Y = %.3f   Z = %.2f \n with connectivity %d", results_table.X(idx), results_table.Y(idx), results_table.Height(idx), results_table.ConnectedUsers(idx))

results_table

% Define the update function
function updateProgress(~, n_points, update_interval)
    persistent progress counter
    if isempty(progress) || isempty(counter)
        progress = 0;
        counter = 0;
    end
    progress = progress + 1;
    counter = counter + 1;
    
    if counter >= update_interval || progress == n_points
        fprintf('Progress: %.2f%% (%d/%d)\n', (progress/n_points)*100, progress, n_points);
        counter = 0; % Reset counter
    end
end
