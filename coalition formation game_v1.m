%% 定义变量
clc; clear;
% 带宽总量
a = 10;
 
num_uav = 20;
num_area = 5;

% 基地到灾难区域的距离
d_m = [0.1, 0.6, 0.5, 0.8, 0.3];
E_const_1 = 10;
E_const_2 = 0.1;
E_const_3 = 0.02;

theta_j = normrnd(0.5, 0.2, [1, num_area]);
gamma_i = rand(1, num_uav);
% 行表示区域,列表示无人机
P_i_j = theta_j' * gamma_i * 5;

% 带宽单位成本
c = 2;

% 20个1~5的数字
uav_dist = round(rand(1,num_uav) * (num_area-1)) + 1;

% end
uavs_per_area = histc(uav_dist, 1:num_area);

filename = "model_args.mat";
save(filename);

%%计算无人机效用

balance_counter = 0;
iter = 0;
%%
while true
    if balance_counter >= num_uav*num_area
        break
    end
    
    idx_uav = mod(iter, num_uav) + 1;
    idx_area = uav_dist(idx_uav);
    
    init_band = ones(1, num_area)' + (10-num_area)/num_area;
    bandwidth = test_projgrad(init_band, P_i_j, uav_dist);
    
    if iter == 0
        non_0_input = uavs_per_area(uavs_per_area ~= 0);
        idx_non_0_input = find(uavs_per_area ~= 0);
        A = -objective(bandwidth(idx_non_0_input),non_0_input, P_i_j, uav_dist);
        E_area_all = 0;
        for i = 1:5
            E_area_all = E_area_all + sum(P_i_j(i, find(uav_dist == i)) + E_const_2 + E_const_3 * (d_m(i))^2);
        end
        Ubo = E_const_1 * A - c * a - E_area_all;
        fprintf("Initial Ubo: %g\n", Ubo);
        
        switch idx_condi
            case 1
                Y_Ubo1 = [Y_Ubo1, Ubo];
            case 2
                Y_Ubo2 = [Y_Ubo2, Ubo];
            case 3
                Y_Ubo3 = [Y_Ubo3, Ubo];
        end
        
    end
    

    for new_idx_area = 1:num_area
        
        idx_area = uav_dist(idx_uav);
        if new_idx_area == idx_area
            continue;
        end
        
        % 计算U_uav_from_old
        uavs_from = find(uav_dist == idx_area);
        U_uav_from_old = zeros(1, num_uav);
        for idx_uav_from_old = uavs_from
            U_area_from_old = E_const_1 * (bandwidth(idx_area) / uavs_per_area(idx_area)) * ...
                sum(log2(1 + (P_i_j(idx_area, uavs_from)) / ((bandwidth(idx_area) / uavs_per_area(idx_area)))));
            
            E_area_from_old = sum(P_i_j(idx_area, uavs_from) + E_const_2 + E_const_3 * (d_m(idx_area))^2);
            
            % 无人机效用
            U_uav_from_old(idx_uav_from_old) = (U_area_from_old - c * bandwidth(idx_area) - E_area_from_old) / uavs_per_area(idx_area);
            
        end
        
        
        % 计算U_uav_to_old
        uavs_to = find(uav_dist == new_idx_area);
        U_uav_to_old = zeros(1, num_uav);
        for idx_uav_to_old = uavs_to
            U_area_to_old = E_const_1 * (bandwidth(new_idx_area) / uavs_per_area(new_idx_area)) * ...
                sum(log2(1 + (P_i_j(new_idx_area, uavs_to)) / ((bandwidth(new_idx_area) / uavs_per_area(new_idx_area)))));

            E_area_to_old = sum(P_i_j(new_idx_area, uavs_to) + E_const_2 + E_const_3 * (d_m(new_idx_area))^2);
            
            % 无人机效用
            U_uav_to_old(idx_uav_to_old) = (U_area_to_old - c * bandwidth(new_idx_area) - E_area_to_old) / uavs_per_area(new_idx_area);
            
        end
        
        
        
        new_uav_dist = uav_dist;
        new_uav_dist(idx_uav) = new_idx_area;
        new_uavs_per_area = histc(new_uav_dist, 1:num_area);
        new_bandwidth = test_projgrad(bandwidth, P_i_j, new_uav_dist);
        
        
        % 计算U_uav_from_new
        uavs_from_new = find(new_uav_dist == idx_area);
        U_uav_from_new = zeros(1, num_uav);
        for idx_uav_from_new = uavs_from_new
            U_area_from_new = E_const_1 * (new_bandwidth(idx_area) / new_uavs_per_area(idx_area)) * ...
                sum(log2(1 + (P_i_j(idx_area, uavs_from_new)) / ((new_bandwidth(idx_area) / new_uavs_per_area(idx_area)))));
          
            E_area_from_new = sum(P_i_j(idx_area, uavs_from_new) + E_const_2 + E_const_3 * (d_m(idx_area))^2);
            
            % 无人机效用
            U_uav_from_new(idx_uav_from_new) = (U_area_from_new - c * new_bandwidth(idx_area) - E_area_from_new) / new_uavs_per_area(idx_area);
            
        end
        
        % 计算U_uav_to_new
        uavs_to_new = find(new_uav_dist == new_idx_area);
        U_uav_to_new = zeros(1, num_uav);
        for idx_uav_to_new = uavs_to_new
            U_area_to_new = E_const_1 * (new_bandwidth(new_idx_area) / new_uavs_per_area(new_idx_area)) * ...
                sum(log2(1 + (P_i_j(new_idx_area, uavs_to_new)) / ((new_bandwidth(new_idx_area) / new_uavs_per_area(new_idx_area)))));
          
            E_area_to_new = sum(P_i_j(new_idx_area, uavs_to_new) + E_const_2 + E_const_3 * (d_m(new_idx_area))^2);
            
            % 无人机效用
            U_uav_to_new(idx_uav_to_new) = (U_area_to_new - c * new_bandwidth(new_idx_area) - E_area_to_new) / new_uavs_per_area(new_idx_area);
            
        end










