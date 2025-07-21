clear,clc;close all;
%% 定义变量
syms x t P_0 P_1 P_2 P_3 Q_0 Q_1 Q_2 Q_3 beta_ip1 beta_i     real;%real确保在实数域
%% 定义矩阵:bernstein基函数矩阵,关联矩阵,控制点矩阵
B_matrix = arrayfun(@(i) B(i, 6), 0:6);

M = [  0 , 1    ,  0  , 0;
        - beta_i, 1   ,    beta_i, 0;
    - 2 * beta_i + 1/12, 5/6    , 2 * beta_i + 1/12, 0;
    0, 1/2, 1/2, 0;
    0,  2 * beta_ip1 + 1/12, 5/6   , - 2 * beta_ip1 + 1/12;
    0,    beta_ip1, 1    ,     - beta_ip1;
    0,   0 , 1     ,  0  ];
P_matrix = [P_0 P_1 P_2 P_3]';
r_x = (B_matrix*M*P_matrix);%曲线方程
BB = (B_matrix*M)';%BB为4行1列调配函数矩阵
%这里我考虑以后可能会用到V和P的控制点之间的关系,最后是V用P表示
V_matrix = M*P_matrix;
for i = 0:6
    % 使用 eval 函数创建变量名，并将对应的矩阵行赋值给变量
    eval(sprintf('v_%d = V_matrix(%d, :);', i, i+1));
end
B0 = BB(1);
B1 = BB(2);
B2 = BB(3);
B3 = BB(4);
BB_flipped = flipud(BB);

% 计算累积和
BB_flipped_cumsum = cumsum(BB_flipped);

% 再次上下翻转以恢复原始顺序
BB_Allocation = flipud(BB_flipped_cumsum);
BB_Allocation_fun_matrix = simplify(BB_Allocation);
C_tilde_0 = BB_Allocation_fun_matrix(1);
C_tilde_1 = BB_Allocation_fun_matrix(2);
C_tilde_2 = BB_Allocation_fun_matrix(3);
C_tilde_3 = BB_Allocation_fun_matrix(4);

%% 飞机姿态拟合
Q0 =[  -0.022100000000000   0.615700000000000   0.559700000000000   0.554100000000000];

Q1 = Q0;
Q2 =[-0.022200000000000   0.306500000000000  -0.245200000000000   0.919500000000000];


Q3 =[-0.022200000000000  -0.009200000000000   0.402600000000000   0.915100000000000];


Q4 =[-0.099516867564659  -0.264250973283681  -0.214237782302162   0.935077530915483];


Q5 =[-0.022096619278684  -0.760408224603820   0.230426734728430   0.606790401451533];


Q6 =[-0.022096619278684  -0.760408224603820   0.230426734728430   0.606790401451533];
% Q5 =[-0.022200000000000   0.306500000000000  -0.245200000000000   0.919500000000000];
% Q6 =[-0.022200000000000   0.306500000000000  -0.245200000000000   0.919500000000000];
% 
% Q4 =[-0.022200000000000  -0.009200000000000   0.402600000000000   0.915100000000000];
% 
% 
% Q3 =[-0.099516867564659  -0.264250973283681  -0.214237782302162   0.935077530915483];
% 
% 
% 
% 
% 
% Q2 =[-0.022096619278684  -0.760408224603820   0.230426734728430   0.606790401451533];
[verts, faces, cindex] = teapotGeometry;
% stlFile = 'airplane.stl';  % 替换为你的STL文件路径
% fv = stlread(stlFile);
% 
% % 获取模型顶点和面
% verts = fv.Points;
% faces = fv.ConnectivityList;

figure
p = patch('Faces',faces,'Vertices',verts,'FaceVertexCData',cindex,'FaceColor','interp');
p.EdgeColor = 'none';      % 移除边缘线条
p.FaceColor = 'none';
p.FaceAlpha = 1;
p.LineStyle = 'none';
colormap("summer");
l = light('Position',[0.3 -3 0],'Style','infinite');
lighting gouraud

axis equal off;
axis off;%去掉坐标轴
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
 %view(180, 50);               % 调整视角
 view(0,19.462825236670692)

beta_array =  (0.025)*[1 1 1 1 1];
 node = [0,0.5,1.5,2,4];
 
syms x;
x=linspace(0,1,3);
for i=1:length(x)
        beta_i_val = beta_array(1); % 当前段的第一个beta值
    beta_ip1_val = beta_array(2); % 当前段的第二个beta值
        % 替换C_tilde中的beta参数
    C_tilde_0_subs = subs(C_tilde_0, {beta_i, beta_ip1}, {beta_i_val, beta_ip1_val});
    C_tilde_1_subs = subs(C_tilde_1, {beta_i, beta_ip1}, {beta_i_val, beta_ip1_val});
    C_tilde_2_subs = subs(C_tilde_2, {beta_i, beta_ip1}, {beta_i_val, beta_ip1_val});
    C_tilde_3_subs = subs(C_tilde_3, {beta_i, beta_ip1}, {beta_i_val, beta_ip1_val});

 % 将替换后的函数转换为可计算的函数句柄, 直接应用于 x(i)
    C_tilde_0_val = double(subs(C_tilde_0_subs, x(i)));
    C_tilde_1_val = double(subs(C_tilde_1_subs, x(i)));
    C_tilde_2_val = double(subs(C_tilde_2_subs, x(i)));
    C_tilde_3_val = double(subs(C_tilde_3_subs, x(i)));

     m11= quaternion_mys(quatmultiply(quatconj(Q0),Q1),C_tilde_1_val);%quatconj()取共轭,((Q0)^(-1)*Q1)^(B1(t))
     m21= quaternion_mys(quatmultiply(quatconj(Q1),Q2),C_tilde_2_val);%quatmultiply( , )四元数乘法
     m31= quaternion_mys(quatmultiply(quatconj(Q2),Q3),C_tilde_3_val);
     QQ1(i,:)=quatmultiply(quatmultiply(quatmultiply(Q0,m11),m21),m31); % 对应的四元数值,第i段单位四元数样条曲线
end
 AQQ1=QQ1;
 % 应该求出AQ1有5行的四元数
 x=linspace(0,1,5);
 for i=1:length(x)
         beta_i_val = beta_array(2); % 当前段的第一个beta值
    beta_ip1_val = beta_array(3); % 当前段的第二个beta值
        % 替换C_tilde中的beta参数
    C_tilde_0_subs = subs(C_tilde_0, {beta_i, beta_ip1}, {beta_i_val, beta_ip1_val});
    C_tilde_1_subs = subs(C_tilde_1, {beta_i, beta_ip1}, {beta_i_val, beta_ip1_val});
    C_tilde_2_subs = subs(C_tilde_2, {beta_i, beta_ip1}, {beta_i_val, beta_ip1_val});
    C_tilde_3_subs = subs(C_tilde_3, {beta_i, beta_ip1}, {beta_i_val, beta_ip1_val});

 % 将替换后的函数转换为可计算的函数句柄, 直接应用于 x(i)
    C_tilde_0_val = double(subs(C_tilde_0_subs, x(i)));
    C_tilde_1_val = double(subs(C_tilde_1_subs, x(i)));
    C_tilde_2_val = double(subs(C_tilde_2_subs, x(i)));
    C_tilde_3_val = double(subs(C_tilde_3_subs, x(i)));
     m12= quaternion_mys(quatmultiply(quatconj(Q1),Q2),C_tilde_1_val);
     m22= quaternion_mys(quatmultiply(quatconj(Q2),Q3),C_tilde_2_val);
     m32= quaternion_mys(quatmultiply(quatconj(Q3),Q4),C_tilde_3_val);
     QQ2(i,:)=quatmultiply(quatmultiply(quatmultiply(Q1,m12),m22),m32); % 对应的四元数值 i+1段
 end
  % 应该求出AQ2有5行的四元数
  AQQ2=QQ2;
  x=linspace(0,1,4);
  for i=1:length(x)
          beta_i_val = beta_array(3); % 当前段的第一个beta值
    beta_ip1_val = beta_array(4); % 当前段的第二个beta值
        % 替换C_tilde中的beta参数
    C_tilde_0_subs = subs(C_tilde_0, {beta_i, beta_ip1}, {beta_i_val, beta_ip1_val});
    C_tilde_1_subs = subs(C_tilde_1, {beta_i, beta_ip1}, {beta_i_val, beta_ip1_val});
    C_tilde_2_subs = subs(C_tilde_2, {beta_i, beta_ip1}, {beta_i_val, beta_ip1_val});
    C_tilde_3_subs = subs(C_tilde_3, {beta_i, beta_ip1}, {beta_i_val, beta_ip1_val});

 % 将替换后的函数转换为可计算的函数句柄, 直接应用于 x(i)
    C_tilde_0_val = double(subs(C_tilde_0_subs, x(i)));
    C_tilde_1_val = double(subs(C_tilde_1_subs, x(i)));
    C_tilde_2_val = double(subs(C_tilde_2_subs, x(i)));
    C_tilde_3_val = double(subs(C_tilde_3_subs, x(i)));
      m13= quaternion_mys(quatmultiply(quatconj(Q2),Q3),C_tilde_1_val);
      m23= quaternion_mys(quatmultiply(quatconj(Q3),Q4),C_tilde_2_val);
      m33= quaternion_mys(quatmultiply(quatconj(Q4),Q5),C_tilde_3_val);
      QQ3(i,:)=quatmultiply(quatmultiply(quatmultiply(Q2,m13),m23),m33); % 对应的四元数值 i+2段
   end
  % 应该求出AQ3有5行的四元数 
   AQQ3=QQ3;
   x=linspace(0,1,8);
  for i=1:length(x)
          beta_i_val = beta_array(4); % 当前段的第一个beta值
    beta_ip1_val = beta_array(5); % 当前段的第二个beta值
 % 将替换后的函数转换为可计算的函数句柄, 直接应用于 x(i)
     C_tilde_0_subs = subs(C_tilde_0, {beta_i, beta_ip1}, {beta_i_val, beta_ip1_val});
    C_tilde_1_subs = subs(C_tilde_1, {beta_i, beta_ip1}, {beta_i_val, beta_ip1_val});
    C_tilde_2_subs = subs(C_tilde_2, {beta_i, beta_ip1}, {beta_i_val, beta_ip1_val});
    C_tilde_3_subs = subs(C_tilde_3, {beta_i, beta_ip1}, {beta_i_val, beta_ip1_val});

    C_tilde_0_val = double(subs(C_tilde_0_subs, x(i)));
    C_tilde_1_val = double(subs(C_tilde_1_subs, x(i)));
    C_tilde_2_val = double(subs(C_tilde_2_subs, x(i)));
    C_tilde_3_val = double(subs(C_tilde_3_subs, x(i)));

     m14= quaternion_mys(quatmultiply(quatconj(Q3),Q4),C_tilde_1_val);%quatconj()取共轭,((Q0)^(-1)*Q1)^(B1(t))
     m24= quaternion_mys(quatmultiply(quatconj(Q4),Q5),C_tilde_2_val);%quatmultiply( , )四元数乘法
     m34= quaternion_mys(quatmultiply(quatconj(Q5),Q6),C_tilde_3_val);
     QQ4(i,:)=quatmultiply(quatmultiply(quatmultiply(Q3,m14),m24),m34); % 对应的四元数值,第i段单位四元数样条曲线
end 
 AQQ4=QQ4;
  QQQ=cat(1,AQQ1,AQQ2,AQQ3,AQQ4);
  QQQ(5,:)=[];
  QQQ(9,:)=[];
  QQQ(13,:)=[];
scale = 5;
q_rotate = [cos(pi/4), 0, 0, sin(pi/4)]; 
% for n = 1:length(QQQ)
%     % 变换顶点
%     Ax = arrayfun(@(i) xz(verts(i, :), q_rotate), 1:size(verts, 1), 'UniformOutput', false);%改变飞机机头朝向，%%鱼头草谁，谁大草谁
%     Ax = cell2mat(Ax(:));
%     Ax = arrayfun(@(i) xz(Ax(i, :), QQQ(n, :)), 1:size(Ax, 1), 'UniformOutput', false);
%     Ax = cell2mat(Ax(:));
%     %g = repmat([n * scale, 0, 0], size(verts, 1), 1);
%     g = repmat([n * scale, 0, scale*3.35 ], size(verts, 1), 1);
%     v_transformed = Ax + 6*g;
%     if n == 1 || n == 3 || n == 7 || n == 10|| n == 17% First transformed teapot is positive (colored red)
%         color = repmat([1 0 0], size(verts, 1), 1);
%         patch('Faces', faces, 'Vertices', v_transformed, 'FaceVertexCData', color, 'FaceColor','interp', 'FaceAlpha', 0.95, 'EdgeColor', 'none');
%      else
%          color = repmat([1 1 0], size(verts, 1), 1);
%          patch('Faces', faces, 'Vertices', v_transformed, 'FaceVertexCData', color, 'FaceColor','interp', 'FaceAlpha', 0.95, 'EdgeColor', 'none');
%     end
%        hold on;
% end
for n = 1:length(QQQ)
    % 变换顶点
    Ax = arrayfun(@(i) xz(verts(i, :), QQQ(n, :)), 1:size(verts, 1), 'UniformOutput', false);
    Ax = cell2mat(Ax(:));
    %g = repmat([n * scale, 0, 0], size(verts, 1), 1);
    g = repmat([n * scale, 0, scale*1.35], size(verts, 1), 1);
    v_transformed = Ax + g;
       % if n == 1 || n == 5 || n == 9 || n == 13|| n == 17% First transformed teapot is positive (colored red)
      if n == 1 || n == 3 || n == 7 || n == 10|| n == 17
        color = repmat([1 0 0], size(verts, 1), 1);
        patch('Faces', faces, 'Vertices', v_transformed, 'FaceVertexCData', color, 'FaceColor','interp', 'FaceAlpha', 0.95, 'EdgeColor', 'none');
    else
        color = (1:size(verts, 1))';
        patch('Faces', faces, 'Vertices', v_transformed, 'FaceVertexCData', color, 'FaceColor','interp', 'FaceAlpha', 0.95, 'EdgeColor', 'none');
    end
       hold on;

end
axis off;
set(gca, 'Position', [0 0 1 1]); % 扩展图像至整个窗口
set(gca, 'LooseInset', get(gca, 'TightInset')); % 去除多余边界

%% 辅助函数
function plotRotatedTeapots(verts, faces, QQ, scale)
    for n = 1:size(QQ, 1)
        Ax = arrayfun(@(i) rotateVectorByQuaternionMatrix(verts(i, :), QQ(n, :)), 1:size(verts, 1), 'UniformOutput', false);
        Ax = cell2mat(Ax(:));
        g = repmat([n * scale, 0, scale * 1.35 * cos(n * 0.668)], size(verts, 1), 1);
        v_transformed = Ax + g;
        if mod(n, 4) == 1
            color = repmat([1 0 0], size(verts, 1), 1);
            patch('Faces', faces, 'Vertices', v_transformed, 'FaceVertexCData', color, 'FaceColor', 'interp', 'FaceAlpha', 0.95, 'EdgeColor', 'none');
        else
            color = (1:size(verts, 1))';
            patch('Faces', faces, 'Vertices', v_transformed, 'FaceVertexCData', color, 'FaceColor', 'interp', 'FaceAlpha', 0.95, 'EdgeColor', 'none');
        end
        hold on;
    end
end

function expression = B(i, k)
    % 计算 Bernstein 基函数的表达式
    syms x real;
        expression = nchoosek(k, i) * x^i * (1 - x)^(k - i);
end

function q_unit = normalizeQuaternion(q)
    %转化单位四元数
    norm_q = sqrt(sum(q.^2));
    q_unit = q / norm_q;
end


function resultQuaternion = quatmultiply(q1, q2)
    %四元数乘积
    s1 = q1(:, 1);
    v1 = q1(:, 2:4);
    s2 = q2(:, 1);
    v2 = q2(:, 2:4);
    resultQuaternion = [s1 .* s2 - dot(v1, v2, 2), s1 .* v2 + s2 .* v1 + cross(v1, v2, 2)];
end

function q_conj = quatconj(q)
    %四元数共轭,单位四元数的逆=共轭
    q_conj = [q(1), -q(2:4)];
end

function q_log = quatlog(q)
    theta = acos(q(1));
    if theta == 0
        q_log = [0, 0, 0, 0];
    else
        q_log = [0, q(2:4) * theta / sin(theta)];
    end
end

function q_exp = quatexp(v)
    % v is assumed to be a 1x3 vector [v1, v2, v3]
    theta = norm(v);
    
    if theta == 0
        q_exp = [1, 0, 0, 0];
    else
        q_exp = [cos(theta), v * sin(theta) / theta];
    end
end


function v_rot = rotateVectorByQuaternionMatrix(v, q)
    
    q = normalizeQuaternion(q);
   
    w = q(1);
    x = q(2);
    y = q(3);
    z = q(4);
  
    R = [1 - 2*y^2 - 2*z^2,    2*x*y + 2*w*z,    2*x*z - 2*w*y;
         2*x*y - 2*w*z,    1 - 2*x^2 - 2*z^2,    2*y*z + 2*w*x;
         2*x*z + 2*w*y,    2*y*z - 2*w*x,    1 - 2*x^2 - 2*y^2
    ];
    
    % Perform the rotation
    v_rot = v * R;
end
    
function q_pow = quaternionPower(q, n)
 
    q = normalizeQuaternion(q);
 
    w = q(1);
    x = q(2);
    y = q(3);
    z = q(4);
 
    theta = acos(w);
 
    sin_theta = sqrt(1 - w^2);
 
    n_theta = n * theta;
 
    w_new = cos(n_theta);
    if sin_theta > 0
        x_new = x / sin_theta * sin(n_theta);
        y_new = y / sin_theta * sin(n_theta);
        z_new = z / sin_theta * sin(n_theta);
    else
        x_new = 0;
        y_new = 0;
        z_new = 0;
    end
 
    q_pow = [w_new, x_new, y_new, z_new];
end

function q = axisAngleToQuaternion(axis, theta)
    % axisAngleToQuaternion - 将旋转轴和旋转角度转换为四元数
    %
    % 输入:
    %   axis - 旋转轴 [x, y, z] (应该是一个单位向量)
    %   theta - 旋转角度 (单位为弧度)
    %
    % 输出:
    %   q - 四元数 [w, x, y, z]，表示旋转
    
    % 确保旋转轴是单位向量
    axis = axis / norm(axis);
    
    % 计算四元数的分量
    w = cos(theta );
    sin_half_theta = sin(theta );
    x = axis(1) * sin_half_theta;
    y = axis(2) * sin_half_theta;
    z = axis(3) * sin_half_theta;
    
    % 返回四元数 [w, x, y, z]
    q = [w, x, y, z];
end
function q = eulerToQuaternion(pitch, roll, yaw)
    % 将角度转换为弧度
    roll = deg2rad(roll);
    pitch = deg2rad(pitch);
    yaw = deg2rad(yaw);
    
    % 计算四元数分量
    qw = cos(yaw/2) * cos(pitch/2) * cos(roll/2) + sin(yaw/2) * sin(pitch/2) * sin(roll/2);
    qx = sin(roll/2) * cos(pitch/2) * cos(yaw/2) - cos(roll/2) * sin(pitch/2) * sin(yaw/2);
    qy = cos(yaw/2) * sin(pitch/2) * cos(roll/2) + sin(yaw/2) * cos(pitch/2) * sin(roll/2);
    qz = cos(roll/2) * cos(pitch/2) * sin(yaw/2) - sin(roll/2) * sin(pitch/2) * cos(yaw/2);
    
    % 输出单位四元数
    q = [qw; qx; qy; qz];
end


