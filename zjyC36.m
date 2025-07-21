function QQQ = zjyC36(Qua_mat, t_step, total_seg_num, beta_array)
    % 检查输入的尺寸和一致性
    [Q_m, Q_n] = size(Qua_mat);
    if Q_m ~= (total_seg_num + 3)
        error('四元数矩阵与整体段数不符');
    end
    if Q_n ~= 4
        error('这不是四元数矩阵');
    end
    if length(beta_array) ~= (total_seg_num + 1)
        error('给的参数值数目不匹配');
    end

    % 定义 B_matrix 和 M 矩阵
    B_matrix = arrayfun(@(i) B(i, 6), 0:6);
    syms beta_i beta_ip1 real;  % 定义符号变量
    M = [0, 1, 0, 0;
        -beta_i, 1, beta_i, 0;
        -2 * beta_i + 1/12, 5/6, 2 * beta_i + 1/12, 0;
        0, 1/2, 1/2, 0;
        0, 2 * beta_ip1 + 1/12, 5/6, -2 * beta_ip1 + 1/12;
        0, beta_ip1, 1, -beta_ip1;
        0, 0, 1, 0];
    
    % 计算调配函数矩阵
    P_matrix = sym('P', [4, 1]);  % 定义控制点的符号矩阵
    BB = (B_matrix * M)';
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
    % 存储整体四元数结果
    QQQ = [];
    x = linspace(0, 1, t_step);

    % 计算每一段的四元数曲线
    for seg = 1:total_seg_num
        % 获取当前段的 beta_i 和 beta_ip1
        beta_i_val = beta_array(seg);
        beta_ip1_val = beta_array(seg + 1);
        
        C_tilde_0_subs = subs(C_tilde_0, {beta_i, beta_ip1}, {beta_i_val, beta_ip1_val});
    C_tilde_1_subs = subs(C_tilde_1, {beta_i, beta_ip1}, {beta_i_val, beta_ip1_val});
    C_tilde_2_subs = subs(C_tilde_2, {beta_i, beta_ip1}, {beta_i_val, beta_ip1_val});
    C_tilde_3_subs = subs(C_tilde_3, {beta_i, beta_ip1}, {beta_i_val, beta_ip1_val});
        
        % 存储该段的四元数结果
        QQ_seg = zeros(length(x), 4);
        
        % 获取控制点（从 Qua_mat 中取相邻的四个四元数）
        Q0 = Qua_mat(seg, :);
        Q1 = Qua_mat(seg + 1, :);
        Q2 = Qua_mat(seg + 2, :);
        Q3 = Qua_mat(seg + 3, :);
        
        for i = 1:length(x)
              C_tilde_0_val = double(subs(C_tilde_0_subs, x(i)));
    C_tilde_1_val = double(subs(C_tilde_1_subs, x(i)));
    C_tilde_2_val = double(subs(C_tilde_2_subs, x(i)));
    C_tilde_3_val = double(subs(C_tilde_3_subs, x(i)));
            % 计算四元数曲线
            m1 = quaternion_mys(quatmultiply(quatconj(Q0), Q1), C_tilde_1_val);
            m2 = quaternion_mys(quatmultiply(quatconj(Q1), Q2), C_tilde_2_val);
            m3 = quaternion_mys(quatmultiply(quatconj(Q2), Q3), C_tilde_3_val);
            QQ_seg(i, :) = quatmultiply(quatmultiply(quatmultiply(Q0, m1), m2), m3);
        end
        % 将结果添加到整体四元数结果中
        QQQ = [QQQ; QQ_seg];
    end
end

% 辅助函数 - 示例B函数
function expression = B(i, k)
    syms x real;
    expression = nchoosek(k, i) * x^i * (1 - x)^(k - i);
end
