clear; clc; close all;

% 定义 beta_array 和四元数控制点
beta_array_xingyan = 5/96 * [1, 1, 1, 1, 1];
beta_array_zjy_025 = 0.025 * [1, 1, 1, 1, 1];
beta_array_zjy_05 = 0.05 * [1, 1, 1, 1, 1];
% beta_array_zjy_02 = 0.02 * [1, 1, 1, 1, 1];
% beta_array_zjy_03 = 0.03 * [1, 1, 1, 1, 1];
beta_array_zjy_015 = 0.015 * [1, 1, 1, 1, 1];
beta_array_zjy_035 = 0.035 * [1, 1, 1, 1, 1]; 

beta_array_tjq = 0.6 * [1, 1, 1, 1, 1];

Q0 = [-0.0221, 0.6157, 0.5597, 0.5541];
Q1 = Q0;
Q2 = [-0.0222, 0.3065, -0.2452, 0.9195];
Q3 = [-0.0222, -0.0092, 0.4026, 0.9151];
Q4 = [-0.0995, -0.2643, -0.2142, 0.9351];
Q5 = [-0.0221, -0.7604, 0.2304, 0.6068];
Q6 = Q5;
Q_mat = [dwh(Q0); dwh(Q1); dwh(Q2); dwh(Q3); dwh(Q4); dwh(Q5); dwh(Q6)];
ture = [0.025740000000000  , 0.026400000000000 ,  0.020600000000000  , 0.024460000000000, 0.026000000000000];

% 样条函数和绘图参数
spline_funcs = {@zjyC36,@zjyC36,@zjyC36,@zjyC36,@zjyC36,@tjqC25, @zjyC36 };
% spline_params = {ture,beta_array_zjy_025,beta_array_zjy_02,beta_array_zjy_03,beta_array_zjy_05,beta_array_tjq,beta_array_xingyan };
spline_params = {ture,beta_array_zjy_025,beta_array_zjy_015,beta_array_zjy_035,beta_array_zjy_05,beta_array_tjq,beta_array_xingyan };
% plot_colors = { [1, 0, 0], [0.9290, 0.6940, 0.1250],[0.4940, 0.1840, 0.5560] };
plot_colors = { ...
    [1,0,0],...%red
    [0, 0.4470, 0.7410], ... % 蓝色
    [0.8500, 0.3250, 0.0980], ... % 橙色
    [0.9290, 0.6940, 0.1250], ... % 金黄色
    [0.4940, 0.1840, 0.5560], ... % 紫色
    [0.4660, 0.6740, 0.1880], ... % 绿色
    [0.3010, 0.7450, 0.9330] ... % 青色
};

 % plot_labels = {'true','025' ,'02 ' ,'03' ,'05' , 'Tan','Xing'};
  plot_labels = {'true','0.25' ,'0.15 ' ,'0.35' ,'0.5' , 'Tan','Xing'};
% plot_labels = {'0.025*[1, 1, 1, 1, 1]','0.05*[1, 1, 1, 1, 1]','0.1*[1, 1, 1, 1, 1]', '0.2*[1, 1, 1, 1, 1]'};
% 初始化用于作用的单位向量
init_vectors = {[1, 0, 0], [0, 0, 1]};
vector_labels = {'nose', 'wing'};

% 预计算每种样条曲线的四元数和旋转角度
results = struct();

for i = 1:numel(spline_funcs)
    % 生成四元数序列并归一化
    QQQ = spline_funcs{i}(Q_mat, 1000, 4, spline_params{i});
    QQQ([1000, 2000, 3000], :) = [];  % 去除跳跃点
    QQQ = QQQ ./ vecnorm(QQQ, 2, 2);  % 归一化
    
    % 计算旋转角度和角速度
    results(i).angles = 2 * acos(QQQ(:, 1))';  % 矢量化计算旋转角度
    results(i).angular_velocity = diff(results(i).angles) / 0.001;
    
    % 计算对不同初始向量旋转后的球面速度
    for j = 1:numel(init_vectors)
        transformed_vectors = quatrotate_vector(init_vectors{j}, QQQ);
        results(i).speed_on_sphere{j} = compute_speed_on_sphere(transformed_vectors);
    end
end

% 绘制旋转角度曲线
figure;
hold on;
for i = 1:numel(spline_funcs)
    t = linspace(0, 4, length(results(i).angles));
    if isequal(spline_params{i}, ture)  % 检查是否为 0.05 的参数
        plot(t, results(i).angles, '-', 'Color', plot_colors{i}, 'LineWidth', 2); % 使用虚线
    else
        plot(t, results(i).angles, '--', 'Color', plot_colors{i}, 'LineWidth', 2); % 使用默认实线
    end
end
legend(plot_labels, 'FontName', 'SimSun', 'Interpreter', 'latex', 'Location', 'northeast', ...
       'FontSize', 16, 'FontWeight', 'bold');
xlabel('Rotational motion time', 'FontSize', 10);
ylabel('Rotation angle ', 'FontSize', 10);
set(gca, 'FontSize', 16, 'LineWidth', 1.5);
xlim([0, 4]);
xticks(0:0.5:4);
title('Rotation Angle Curve');
hold off;


% 绘制角速度曲线
figure;
hold on;
for i = 1:numel(spline_funcs)
    t = linspace(0, 4, length(results(i).angular_velocity));
    if isequal(spline_params{i}, ture)  % 检查是否为 0.05 的参数
        plot(t, results(i).angular_velocity, '-', 'Color', plot_colors{i}, 'LineWidth', 2); % 使用虚线
    else
        plot(t, results(i).angular_velocity, '--', 'Color', plot_colors{i}, 'LineWidth', 2); % 使用默认实线
    end
end
legend(plot_labels, 'FontName', 'SimSun', 'Interpreter', 'latex', 'Location', 'northeast', ...
       'FontSize', 16, 'FontWeight', 'bold');
xlabel('Rotational motion time', 'FontSize', 10);
ylabel('Angular velocity ', 'FontSize', 10);
set(gca, 'FontSize', 16, 'LineWidth', 1.5);
xlim([0, 4]);
ylim([-1, 1]);
yticks(-1:0.2:1);
xticks(0:0.5:4);
title('Angular Velocity Curve');
hold off;


% 绘制球面速度变化曲线 (作用在 [1,0,0] 和 [0,0,1] )
for j = 1:numel(init_vectors)
    figure;
    hold on;
    for i = 1:numel(spline_funcs)
        t = linspace(0, 4, length(results(i).speed_on_sphere{j}));
        if isequal(spline_params{i}, ture)  % 检查是否为 0.05 的参数
            plot(t, results(i).speed_on_sphere{j}, '-', 'Color', plot_colors{i}, 'LineWidth', 2, ...
                 'DisplayName', plot_labels{i}); % 使用虚线
        else
            plot(t, results(i).speed_on_sphere{j}, '--', 'Color', plot_colors{i}, 'LineWidth', 2, ...
                 'DisplayName', plot_labels{i}); % 使用默认实线
        end
    end
    legend(plot_labels, 'FontName', 'SimSun', 'Interpreter', 'latex', 'Location', 'northeast', ...
           'FontSize', 16, 'FontWeight', 'bold');
    xlabel('Rotational motion time', 'FontSize', 10);
    ylabel('Speed on sphere' , 'FontSize', 10);
    set(gca, 'FontSize', 16, 'LineWidth', 1.5);
    xlim([0, 4]);
    xticks(0:0.5:4);
    title(['Speed on Sphere (Vector ' vector_labels{j} ')']);
    hold off;
end


%% 辅助函数

% 计算在球面上的速度变化
function speed = compute_speed_on_sphere(vectors)
    speed = vecnorm(diff(vectors), 2, 2)' / 0.001;  % 使用矢量化计算速度
end

% 将四元数序列应用到初始向量上
function transformed_vectors = quatrotate_vector(v, QQQ)
    % 使用矩阵运算实现四元数旋转，避免循环
    v_q = [zeros(size(QQQ, 1), 1), repmat(v, size(QQQ, 1), 1)];
    q_conj = [QQQ(:, 1), -QQQ(:, 2:4)];
    transformed_vectors = quatmultiply(quatmultiply(QQQ, v_q), q_conj);
    transformed_vectors = transformed_vectors(:, 2:4);
end

% 示例单位化四元数函数
function q_unit = dwh(q)
    q_unit = q / norm(q);
end

% 四元数乘法
function q = quatmultiply(q1, q2)
    % q1和q2均为 [w, x, y, z] 四元数
    w = q1(:,1).*q2(:,1) - q1(:,2).*q2(:,2) - q1(:,3).*q2(:,3) - q1(:,4).*q2(:,4);
    x = q1(:,1).*q2(:,2) + q1(:,2).*q2(:,1) + q1(:,3).*q2(:,4) - q1(:,4).*q2(:,3);
    y = q1(:,1).*q2(:,3) - q1(:,2).*q2(:,4) + q1(:,3).*q2(:,1) + q1(:,4).*q2(:,2);
    z = q1(:,1).*q2(:,4) + q1(:,2).*q2(:,3) - q1(:,3).*q2(:,2) + q1(:,4).*q2(:,1);
    q = [w, x, y, z];
end
