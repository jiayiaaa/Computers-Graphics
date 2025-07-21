%% 此程序用来生成6个球面图，每个图包含真实路径（红色）和对比路径
clear; clc; close all;

data = load('q_true.mat');

% 样条函数和参数
spline_funcs = {@zjyC36, @zjyC36, @zjyC36, @zjyC36, @zjyC36, @tjqC25, @zjyC36};


spline_params = {
   [0.025740000000000  , 0.026400000000000 ,  0.020600000000000  , 0.024460000000000, 0.026000000000000], ...
    0.025 * [1, 1, 1, 1, 1], ...
    0.015 * [1, 1, 1, 1, 1], ...
    0.035 * [1, 1, 1, 1, 1], ...
    0.05 * [1, 1, 1, 1, 1], ...
    0.6 * [1, 1, 1, 1, 1], ...
    5/96 * [1, 1, 1, 1, 1] ...
};
plot_colors = {
    [0, 0.4470, 0.7410], [0.8500, 0.3250, 0.0980], [0.9290, 0.6940, 0.1250], ...
    [0.4940, 0.1840, 0.5560], [0.4660, 0.6740, 0.1880], [0.3010, 0.7450, 0.9330], ...
    [0.6350, 0.0780, 0.1840]
};
plot_labels = {'true', '025', '02', '03', '05', 'Tan', 'Xingyan'};

% 初始四元数控制点
Q0 = [-0.0221, 0.6157, 0.5597, 0.5541];
Q2 = [-0.0222, 0.3065, -0.2452, 0.9195];
Q3 = [-0.0222, -0.0092, 0.4026, 0.9151];
Q4 = [-0.0995, -0.2643, -0.2142, 0.9351];
Q5 = [-0.0221, -0.7604, 0.2304, 0.6068];
Q_mat = [Q0; Q0; Q2; Q3; Q4; Q5; Q5];

%% 生成真实路径数据
v = [0, 1, 0]; % 初始球面点

% 生成真实路径（第一个样条函数）
true_func = spline_funcs{1};
true_params = spline_params{1};
QQQ_true = true_func(Q_mat, 750, 4, true_params);
QQQ_true = QQQ_true ./ vecnorm(QQQ_true, 2, 2);

s_true = [];
for n = 1:size(QQQ_true, 1)
    rotated_point = xz(v, QQQ_true(n, :));
    s_true = cat(1, s_true, rotated_point);
end

%% 循环生成6个对比图
for i = 2:numel(spline_funcs)  % 只处理2-7号样条函数，生成6幅图
    % 生成当前路径数据
    current_func = spline_funcs{i};
    current_params = spline_params{i};
    QQQ_current = current_func(Q_mat, 750, 4, current_params);
    QQQ_current = QQQ_current ./ vecnorm(QQQ_current, 2, 2);
    
    s_current = [];
    for n = 1:size(QQQ_current, 1)
        rotated_point = xz(v, QQQ_current(n, :));
        s_current = cat(1, s_current, rotated_point);
    end
    
    % 创建新图形
    figure;
    axis equal;
    [X, Y, Z] = sphere(100);
    surf(X, Y, Z, 'EdgeColor', 'none');
    shading interp;
    alpha(0.2);
    hold on;
    mesh(X, Y, Z, 'FaceColor', 'none', 'EdgeColor', [0.5, 0.5, 0.5], 'LineStyle', ':');
    set(gcf, 'Color', [1 1 1]);
    set(gca, 'XColor', 'none', 'YColor', 'none', 'ZColor', 'none');
    view(10, 5);
    
    % 绘制真实路径（红色）
    plot3(s_true(:,1), s_true(:,2), s_true(:,3), 'r-', 'LineWidth', 2);
    
    % 绘制当前路径（原颜色）
    plot3(s_current(:,1), s_current(:,2), s_current(:,3), 'b--', 'LineWidth', 2);
    
    % 标注关键点（使用当前路径）
    key_points = s_current([1, 751, 1500, 2250, 3000], :);
    plot3(key_points(:,1), key_points(:,2), key_points(:,3), '.k', 'MarkerSize', 20);
    % text(key_points(3,1), key_points(3,2), key_points(3,3), '$q_3$',...
         % 'FontSize', 24, 'Interpreter', 'latex', 'FontName', 'Times New Roman');
    
    axis off;
    grid off;
    set(gcf, 'Position', [200, 190, 600, 600]);
    set(gca, 'Position', [0 0 1 1]);
    hold off;
end