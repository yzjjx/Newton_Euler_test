%% 用于生成符号化的^{i-1}_{i}T矩阵
% 输入DH参数和i，输出矩阵，默认输入为弧度
syms theta_1 d_2 L_2 theta_3 real

DH_matrix = [
    0     0   0     theta_1;
    pi/2  0   d_2   0;
    0     0   L_2   theta_3];

get_T = @(alp,a,d,th)[
    cos(th)              -sin(th)                 0              a;
    sin(th)*cos(alp)      cos(th)*cos(alp)       -sin(alp)      -sin(alp)*d;
    sin(th)*sin(alp)      cos(th)*sin(alp)        cos(alp)       cos(alp)*d;
    0                     0                       0              1
];

[rows_DH,cols_DH] = size(DH_matrix);
T_total = eye(4);

% 循环代入参数,将1:row_DH替换为1，得到T_i就是^0_1T；替换为2，得到T_i就是^1_2T
for i = 1:rows_DH
    alp_i = DH_matrix(i, 1);
    a_i   = DH_matrix(i, 2);
    d_i   = DH_matrix(i, 3);
    th_i  = DH_matrix(i, 4);

    % 将参数代入到get_T里面
    T_i = get_T(alp_i, a_i, d_i, th_i);
    % 累乘
    T_total = T_total * T_i;
end

