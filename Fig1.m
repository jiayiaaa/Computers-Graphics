%% 此程序用来一次画出7个球面图,7个飞机3D图
clear; clc; close all;

data = load('q_true.mat');
% 样条函数和参数
spline_funcs = {@zjyC36, @zjyC36, @zjyC36, @zjyC36, @zjyC36, @tjqC25, @zjyC36};
spline_params = {
     [0.025740000000000  , 0.026400000000000 ,  0.020600000000000  , 0.024460000000000, 0.026000000000000], ...  
    0.025 * [1, 1, 1, 1, 1], ...
    0.02 * [1, 1, 1, 1, 1], ...
    0.03 * [1, 1, 1, 1, 1], ...
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
% % 画7个球面图
% 
% % 初始化点
% v = [0, 1, 0]; % 初始球面点
% 
% 
% % 循环绘制每种样条函数
% for i = 1:numel(spline_funcs)
%     % 样条插值并归一化
%     QQQ = spline_funcs{i}(Q_mat, 750, 4, spline_params{i});
%     QQQ = QQQ ./ vecnorm(QQQ, 2, 2); % 归一化四元数
% 
%     % 计算旋转轨迹
%     s = [];
%     for n = 1:size(QQQ, 1)
%         rotated_point = xz(v, QQQ(n, :));
%         s = cat(1, s, rotated_point);
%     end
% figure;
% axis equal;
% [X, Y, Z] = sphere(100); % 生成球体数据
% 
% % 设置单一颜色背景的球体
% colormap(gray); % 使用灰度映射作为单色背景
% surf(X, Y, Z, 'EdgeColor', 'none'); % 取消边缘线框
% shading interp;
% alpha(0.2); % 设置透明度
% 
% % 绘制经度和纬度的虚线
% hold on;
% mesh(X, Y, Z, 'FaceColor', 'none', 'EdgeColor', [0.5, 0.5, 0.5], 'LineStyle', ':');
% set(gcf, 'Color', [1 1 1]); % 设置窗口背景为白色
% set(gca, 'XColor', 'none', 'YColor', 'none', 'ZColor', 'none'); % 隐藏坐标轴
% 
% view(10, 5);
%     % 绘制轨迹
%     plot3(s(:, 1), s(:, 2), s(:, 3), '-', 'Color',[0.1 0.5,0.9],'LineWidth', 2);
% 
% 
%     % 特殊点标注
%     key_points = s([1, 751, 1500 , 2250,3000], :); % 提取关键点
%     plot3(key_points(:, 1), key_points(:, 2), key_points(:, 3), '.k', 'MarkerSize', 20);
%     text(key_points(3, 1), key_points(3, 2), key_points(3, 3), '$q_3$', ...
%         'FontSize', 24, 'Interpreter', 'latex', 'FontName', 'Times New Roman');
%      axis off;
%     %legend show;
%     grid off;
%       % 设置图形窗口为方形
%     set(gcf, 'Position', [200, 200, 600, 600]); % 设置图形窗口大小为 600x600 像素
%     set(gca, 'Position', [0 0 1 1]); % 扩展图像至整个窗口
% 
%     hold off;
% end
%% 画飞机3D
stlFile = 'airplane.stl';  % 替换为你的STL文件路径
fv = stlread(stlFile);

% 获取模型顶点和面
verts = fv.Points;
faces = fv.ConnectivityList;
for i = 1:numel(spline_funcs)
    % 样条插值并归一化
    QQQ = spline_funcs{i}(Q_mat, 5, 4, spline_params{i});
      QQQ(5,:)=[];
  QQQ(9,:)=[];
  QQQ(13,:)=[];
    QQQ = QQQ ./ vecnorm(QQQ, 2, 2); % 归一化四元数
figure
p = patch('Faces', faces, 'Vertices', verts, 'FaceColor', 'interp'); % 设置面颜色为插值
p.EdgeColor = 'none';           % 移除边缘线条
p.FaceAlpha = 1;                % 完全不透明
colormap("jet");                % 使用 jet 颜色映射，可以尝试其他映射

% 设置光源，增加多个光源以减少阴影区域
light('Position', [1, 1, 1], 'Style', 'infinite');
light('Position', [-1, -1, -1], 'Style', 'infinite');
light('Position', [0, 0, 1], 'Style', 'local');  % 本地光源

% 使用 Phong 照明
lighting phong;

% 增加环境光的强度
material([0.3, 0.7, 0.4]);  % [环境光, 漫反射, 高光] 的强度设置

axis equal;
axis off;                    % 去掉坐标轴
view(180, 50);               % 调整视角

scale = 5;
q_rotate = [cos(pi/4), 0, 0, sin(pi/4)]; 
for n = 1:length(QQQ)
    % 变换顶点
    Ax = arrayfun(@(i) xz(verts(i, :), q_rotate), 1:size(verts, 1), 'UniformOutput', false);
    Ax = cell2mat(Ax(:));
    Ax = arrayfun(@(i) xz(Ax(i, :), QQQ(n, :)), 1:size(Ax, 1), 'UniformOutput', false);
    Ax = cell2mat(Ax(:));
    %g = repmat([n * scale, 0, 0], size(verts, 1), 1);
    g = repmat([n * scale, 0, scale*3.35 ], size(verts, 1), 1);
    v_transformed = Ax + 6*g;
    if n == 1 || n == 5 || n == 9 || n == 13|| n == 17% First transformed teapot is positive (colored red)
        color = repmat([1 0 0], size(verts, 1), 1);
        patch('Faces', faces, 'Vertices', v_transformed, 'FaceVertexCData', color, 'FaceColor','interp', 'FaceAlpha', 0.95, 'EdgeColor', 'none');
     else
         color = repmat([1 1 0], size(verts, 1), 1);
         patch('Faces', faces, 'Vertices', v_transformed, 'FaceVertexCData', color, 'FaceColor','interp', 'FaceAlpha', 0.95, 'EdgeColor', 'none');
    end
       hold on;
end
hold off;
axis off;
set(gca, 'Position', [0 0 1 1]); % 扩展图像至整个窗口
set(gca, 'LooseInset', get(gca, 'TightInset')); % 去除多余边界 
end
 figure
p = patch('Faces', faces, 'Vertices', verts, 'FaceColor', 'interp'); % 设置面颜色为插值
p.EdgeColor = 'none';           % 移除边缘线条
p.FaceAlpha = 1;                % 完全不透明
colormap("jet");                % 使用 jet 颜色映射，可以尝试其他映射

% 设置光源，增加多个光源以减少阴影区域
light('Position', [1, 1, 1], 'Style', 'infinite');
light('Position', [-1, -1, -1], 'Style', 'infinite');
light('Position', [0, 0, 1], 'Style', 'local');  % 本地光源

% 使用 Phong 照明
lighting phong;

% 增加环境光的强度
material([0.3, 0.7, 0.4]);  % [环境光, 漫反射, 高光] 的强度设置

axis equal;
axis off;                    % 去掉坐标轴
view(180, 50);               % 调整视角

scale = 5;
q_rotate = [cos(pi/4), 0, 0, sin(pi/4)]; 
for n = 1:length(QQQ)
    % 变换顶点
    Ax = arrayfun(@(i) xz(verts(i, :), q_rotate), 1:size(verts, 1), 'UniformOutput', false);
    Ax = cell2mat(Ax(:));
    Ax = arrayfun(@(i) xz(Ax(i, :), QQQ(n, :)), 1:size(Ax, 1), 'UniformOutput', false);
    Ax = cell2mat(Ax(:));
    %g = repmat([n * scale, 0, 0], size(verts, 1), 1);
    g = repmat([n * scale, 0, scale*3.35 ], size(verts, 1), 1);
    v_transformed = Ax + 6*g;
    if n == 1 || n == 5 || n == 9 || n == 13|| n == 17% First transformed teapot is positive (colored red)
        color = repmat([1 0 0], size(verts, 1), 1);
        patch('Faces', faces, 'Vertices', v_transformed, 'FaceVertexCData', color, 'FaceColor','interp', 'FaceAlpha', 0.95, 'EdgeColor', 'none');
     % else
     %     color = repmat([1 1 0], size(verts, 1), 1);
     %     patch('Faces', faces, 'Vertices', v_transformed, 'FaceVertexCData', color, 'FaceColor','interp', 'FaceAlpha', 0.95, 'EdgeColor', 'none');
    end
end
