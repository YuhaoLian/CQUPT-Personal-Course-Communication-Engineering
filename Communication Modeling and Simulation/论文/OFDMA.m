% %% Petri网络在OFDMA系统中的应用
% 
% % 定义 Petri 网络
% % 资源分配变迁
% allocate_res = [0 0 1 0; 0 1 0 0; 1 0 0 0; 0 0 0 0];
% % 数据传输变迁
% transfer_data = [0 0 0 1; 0 0 0 1; 0 0 0 1; 0 0 0 0];
% % 初始标记
% marking = [0 0 0 1];
% 
% % 定义仿真参数
% num_users = 8;      % 用户数量
% num_slots = 8;      % 仿真时隙数量
% tx_power = 23;      % 发射功率
% channel_gain = [0.3 0.4 0.2 0.3 0.4 0.2 0.1 0.5];    % 信道增益
% noise_power = 1;    % 噪声功率
% bandwidth = 10e6;   % 带宽
% 
% % 定义仿真变量
% user_bits = zeros(1, num_users);    % 用户传输数据量
% user_rate = zeros(1, num_users);    % 用户传输率
% 
% % 执行仿真
% for i = 1:num_slots
%     % 执行资源分配变迁
%     marking = marking * allocate_res;
%     % 计算每个用户的传输速率
%     user_rate = tx_power .* channel_gain ./ (noise_power + user_bits .* bandwidth);
%     % 执行数据传输变迁
%     transfer_bits = min(user_rate) * bandwidth;    % 计算每个用户的传输数据量
%     user_bits = user_bits + transfer_bits;         % 更新用户传输数据量
%     marking = marking * transfer_data;
%     % 输出每个时隙的用户传输率
%     disp(['Slot ' num2str(i) ' user rate: ' num2str(user_rate)]);
% end


%% Petri网络在OFDMA系统中的应用（可视化版）

% 定义 Petri 网络
% 资源分配变迁
allocate_res = [0 0 1 0; 0 1 0 0; 1 0 0 0; 0 0 0 0];
% 数据传输变迁
transfer_data = [0 0 0 1; 0 0 0 1; 0 0 0 1; 0 0 0 0];
% 初始标记
marking = [0 0 0 1];

% 定义仿真参数
num_users = 3;      % 用户数量
num_slots = 5;      % 仿真时隙数量
tx_power = 23;      % 发射功率
channel_gain = [0.3 0.4 0.2];    % 信道增益
noise_power = 1;    % 噪声功率
bandwidth = 10e6;   % 带宽

% 定义仿真变量
user_bits = zeros(1, num_users);    % 用户传输数据量
user_rate = zeros(1, num_users);    % 用户传输率
user_bits_history = zeros(num_users, num_slots);   % 用户传输数据量历史记录

% 执行仿真
for i = 1:num_slots
    % 执行资源分配变迁
    marking = marking * allocate_res;
    % 计算每个用户的传输速率
    user_rate = tx_power .* channel_gain ./ (noise_power + user_bits .* bandwidth);
    % 执行数据传输变迁
    transfer_bits = min(user_rate) * bandwidth;    % 计算每个用户的传输数据量
    user_bits = user_bits + transfer_bits;         % 更新用户传输数据量
    marking = marking * transfer_data;
    % 记录每个用户的传输数据量
    user_bits_history(:, i) = user_rate';
    % 输出每个时隙的用户传输率
    disp(['Slot ' num2str(i) ' user rate: ' num2str(user_rate)]);
end

% 画出每个用户的传输数据量随时间的变化趋势图
figure;
for i = 1:num_users
    subplot(num_users, 1, i);
    plot(1:num_slots, user_bits_history(i, :), '-o');
    xlabel('Time slot');
    ylabel('Data rate');
    title(['User ' num2str(i) ' data transmission']);
end

