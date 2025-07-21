clear,clc;

data = load('q_true.mat');

beta_array_xingyan = 5/96 * [1, 1, 1, 1, 1];
beta_array_zjy_025 = 0.025 * [1, 1, 1, 1, 1];
beta_array_zjy_05 = 0.05 * [1, 1, 1, 1, 1];
beta_array_zjy_035 = 0.035 * [1, 1, 1, 1, 1];
beta_array_zjy_015 = 0.015 * [1, 1, 1, 1, 1];

beta_array_tjq = 0.6 * [1, 1, 1, 1, 1];

Q0 = [-0.0221, 0.6157, 0.5597, 0.5541];
Q1 = Q0;
Q2 = [-0.0222, 0.3065, -0.2452, 0.9195];
Q3 = [-0.0222, -0.0092, 0.4026, 0.9151];
Q4 = [-0.0995, -0.2643, -0.2142, 0.9351];
Q5 = [-0.0221, -0.7604, 0.2304, 0.6068];
Q6 = Q5;
Q_mat = [dwh(Q0); dwh(Q1); dwh(Q2); dwh(Q3); dwh(Q4); dwh(Q5); dwh(Q6)];
QQQ_xing = zjyC36(Q_mat, 750, 4, beta_array_xingyan);
QQQ_tan = tjqC25(Q_mat, 750, 4, beta_array_tjq );
QQQ_025 = zjyC36(Q_mat, 750, 4, beta_array_zjy_025);
QQQ_015  = zjyC36(Q_mat, 750, 4, beta_array_zjy_015 );
QQQ_035  =zjyC36(Q_mat, 750, 4, beta_array_zjy_035);
QQQ_05 = zjyC36(Q_mat, 750, 4, beta_array_zjy_05);




q_fit_tan = QQQ_tan;
q_fit_zhu_025 = QQQ_025;
q_fit_zhu_015= QQQ_015;
q_fit_zhu_035 = QQQ_035;
q_fit_zhu_05 = QQQ_05;
q_fit_xing = QQQ_xing;


% 初始化球面距离数组
N = 3000;  % 获取数据的行数，即时间步长
spherical_distances_tan = zeros(N, 1);  % 存储 q_fit_tan 的球面距离
spherical_distances_zhu_025 = zeros(N, 1);   % 存储 q_fit_zhu 的球面距离
spherical_distances_zhu_015 = zeros(N, 1); 
spherical_distances_zhu_035 = zeros(N, 1); 
spherical_distances_zhu_05 = zeros(N, 1); 
spherical_distances_xing = zeros(N, 1); 
% 计算球面距离
for i = 1:N
    % 计算 q_fit_tan 的球面距离
    dot_product_tan = dot(data.QQQ(i,:), q_fit_tan(i,:));
    dot_product_tan = min(max(dot_product_tan, -1), 1);
    spherical_distances_tan(i) = 2 * acos(abs(dot_product_tan));

    % 计算 q_fit_zhu 的球面距离
    dot_product_zhu_025 = dot(data.QQQ(i,:), q_fit_zhu_025(i,:));
    dot_product_zhu_025 = min(max(dot_product_zhu_025, -1), 1);
    spherical_distances_zhu_025(i) = 2 * acos(abs(dot_product_zhu_025)); 

    dot_product_zhu_015 = dot(data.QQQ(i,:), q_fit_zhu_015(i,:));
    dot_product_zhu_015 = min(max(dot_product_zhu_015, -1), 1);
    spherical_distances_zhu_015(i) = 2 * acos(abs(dot_product_zhu_015)); 

    dot_product_zhu_035 = dot(data.QQQ(i,:), q_fit_zhu_035(i,:));
    dot_product_zhu_035 = min(max(dot_product_zhu_035, -1), 1);
    spherical_distances_zhu_035(i) = 2 * acos(abs(dot_product_zhu_035)); 

    dot_product_zhu_05 = dot(data.QQQ(i,:), q_fit_zhu_05(i,:));
    dot_product_zhu_05 = min(max(dot_product_zhu_05, -1), 1);
    spherical_distances_zhu_05(i) = 2 * acos(abs(dot_product_zhu_05)); 

        % 计算 q_fit_xing 的球面距离
    dot_product_xing = dot(data.QQQ(i,:), q_fit_xing(i,:));
    dot_product_xing = min(max(dot_product_xing, -1), 1);
    spherical_distances_xing(i) = 2 * acos(abs(dot_product_xing));
end

% 绘制球面距离随时间变化的图像
figure;
hold on;  % 保持当前图形'0.25' ,'0.2 ' ,'0.3' ,'0.5' , 'Tan','Xing'
plot_colors = { ...
    [0, 0.4470, 0.7410], ... % 蓝色
    [0.8500, 0.3250, 0.0980], ... % 橙色
    [0.9290, 0.6940, 0.1250], ... % 金黄色
    [0.4940, 0.1840, 0.5560], ... % 紫色
    [0.4660, 0.6740, 0.1880], ... % 绿色
    [0.3010, 0.7450, 0.9330] ... % 青色
};
t = linspace(0, 4, 3000);
plot(t, spherical_distances_zhu_025, 'Color', plot_colors{1},  'LineWidth', 2, 'DisplayName', '0.25');
plot(t, spherical_distances_zhu_015, 'Color', plot_colors{2},  'LineWidth', 2, 'DisplayName', '0.2');
plot(t, spherical_distances_zhu_035, 'Color', plot_colors{3},  'LineWidth', 2, 'DisplayName', '0.3');
plot(t, spherical_distances_zhu_05, 'Color', plot_colors{4}, 'LineWidth', 2, 'DisplayName', '0.5');
plot(t, spherical_distances_tan,  'Color', plot_colors{5}, 'LineWidth', 2, 'DisplayName', 'Tan');
plot(t, spherical_distances_xing, 'Color', plot_colors{6},  'LineWidth', 2, 'DisplayName', 'Xing');
xlabel('Rotational motion time');
ylabel('Distance error');
% title('Spherical Distance Between True Quaternion and Fitted Quaternions');
  plot_labels = {'0.25' ,'0.15 ' ,'0.35' ,'0.5' , 'Tan','Xing'};
legend(plot_labels, 'FontName', 'SimSun', 'Interpreter', 'latex', 'Location', 'northeast', ...
       'FontSize', 16, 'FontWeight', 'bold');  % 显示图例
set(gca, 'FontSize', 16, 'LineWidth', 1.5);
xlim([0, 4]);
xticks(0:0.5:4);
hold off;  % 释放当前图形