%% 牛顿欧拉外推+内推
% 给定初始的速度和角速度，给定关节数量，通过递推获得每一个关节的速度和角速度
% 首先确定关节的数量
N = 3;
Z = [0;0;1];

% 初始化速度和角速度数组，因为w从0开始，matlab的下标从1开始，整体向右移动
% sym保证了d_theta_1能进去，因为zeros()默认是double类型
w   = sym(zeros(3,N+1));
dw  = sym(zeros(3,N+1));
dv  = sym(zeros(3,N+1));
dvc = sym(zeros(3,N+1));
Pc  = sym(zeros(3,N+1));
F   = sym(zeros(3,N+1));
N_n = sym(zeros(3,N+1));
f   = sym(zeros(3,N+2));
n   = sym(zeros(3,N+2));
tau = sym(zeros(1,N+2));
% 定义导数,因为不涉及到对cos和sin的求导，因此暂时不需要用到diff函数
syms d_theta_1 d_theta_2 d_theta_3 real
% 关节3的转角恒为0度
d_theta = [d_theta_1, d_theta_2, 0];

syms dd_theta_1 dd_theta_2 dd_theta3 real
dd_theta = [dd_theta_1,dd_theta_2,0];

% 设定边界条件,注意此时递推下标统一向右平移一位
w(:,1) = [0;0;0];
dw(:,1) = [0;0;0];
syms g real;
dv(:,1) = [0;g;0];

f(:,N+2) = [0;0;0];
n(:,N+2) = [0;0;0];
% 定义连杆质量
syms m1 m2 m3 real
m = [0,m1,m2,0];

% 定义惯性张量矩阵
I = cell(1, N+1);
I{2} = [0, 0, 0;
        0, 0, 0;
        0, 0, 0];
I{3} = [0, 0, 0;
        0, 0, 0;
        0, 0, 0];
I{4} = [0, 0, 0;
        0, 0, 0;
        0, 0, 0];
%% 在矩阵T中提取旋转矩阵
syms theta_1 theta_2 L_1 L_2 real

DH_matrix = [
    0     0     0     theta_1;
    0     L_1   0     theta_2;
    0     L_2   0     0];

get_T = @(alp,a,d,th)[
    cos(th)              -sin(th)                 0              a;
    sin(th)*cos(alp)      cos(th)*cos(alp)       -sin(alp)      -sin(alp)*d;
    sin(th)*sin(alp)      cos(th)*sin(alp)        cos(alp)       cos(alp)*d;
    0                     0                       0              1
];

[rows_DH,cols_DH] = size(DH_matrix);
T_total = eye(4);

% 存放旋转矩阵、平移矩阵的元胞数组
R_store = cell(1,rows_DH);
R_store_2 = cell(1,rows_DH);
P_store = cell(1,rows_DH);

% 循环代入参数,将1:row_DH替换为1，得到T_i就是^0_1T；替换为2，得到T_i就是^1_2T
for i = 1:rows_DH
    alp_i = DH_matrix(i, 1);
    a_i   = DH_matrix(i, 2);
    d_i   = DH_matrix(i, 3);
    th_i  = DH_matrix(i, 4);

    % 将参数代入到get_T里面
    T_i = get_T(alp_i, a_i, d_i, th_i);
    % 提取旋转矩阵的部分
    R_i = T_i(1:3,1:3);
    % 矩阵转置并且保存
    R_store{i} = R_i.';
    R_store_2{i} = R_i;
    % 提取位置矩阵的部分
    P_i = T_i(1:3,4);
    % 矩阵转置并且保存
    P_store{i} = P_i;
end

%% 定义Pc为连杆质心的位置矢量
Pc(:,2) = [L_1; 0; 0];
Pc(:,3) = [L_2; 0; 0];

%% 外推计算
for i = 2:N+1
    % Newton-Euler外推计算
    w(:,i)   = R_store{i-1}*w(:,i-1)+d_theta(i-1)*Z;
    dw(:,i)  = R_store{i-1}*dw(:,i-1)+cross((R_store{i-1}*w(:,i-1)),d_theta(i-1)*Z)+dd_theta(i-1)*Z;
    dv(:,i)  = R_store{i-1}*(cross(dw(:,i-1),P_store{i-1})+cross(w(:,i-1),cross(w(:,i-1),P_store{i-1}))+dv(:,i-1));
    dvc(:,i) = cross(dw(:,i),Pc(:,i))+ cross(w(:,i),cross(w(:,i),Pc(:,i)))+dv(:,i);
    F(:,i)   = m(i)*dvc(:,i);
    N_n(:,i) = I{i}*dw(:,i)+cross(w(:,i),(I{i}*w(:,i)));
end

%% 内推计算
for i = N+1:-1:2
    f(:,i-1)   = R_store_2{i-1}*f(:,i)+F(:,i-1);
    n(:,i-1)   = N_n(:,i-1)+R_store_2{i-1}*n(:,i)+cross(Pc(:,i-1),F(:,i-1))+cross(P_store{i-1},R_store_2{i-1}*f(:,i));
    tau(:,i-1) = (n(:,i-1).')*Z;
end

%% 对关节相互作用力f内推的验证
f_1 = [cos(theta_2) -sin(theta_2)  0;
       sin(theta_2)  cos(theta_2)  0;
       0             0             1]...
       *[m2*L_1*sin(theta_2)*dd_theta_1-m2*L_1*cos(theta_2)*d_theta_1^2+m2*g*(sin(theta_1)*cos(theta_2)+sin(theta_2)*cos(theta_1))-m2*L_2*(d_theta_1+d_theta_2)^2;
                m2*L_1*cos(theta_2)*dd_theta_1+m2*L_1*sin(theta_2)*d_theta_1^2+m2*g*(cos(theta_1)*cos(theta_2)-sin(theta_2)*sin(theta_1))+m2*L_2*(dd_theta_1+dd_theta_2);0]+[-m1*L_1*d_theta_1^2+m1*g*sin(theta_1);m1*L_1*dd_theta_1+m1*g*cos(theta_1);0];
f_test=simplify(f_1-f(:,2));

%% 对关节相互作用力矩n内推的验证
n_2=[0;0;m2*L_1*L_2*cos(theta_2)*dd_theta_1+m2*L_1*L_2*sin(theta_2)*d_theta_1^2+m2*L_2*g*cos(theta_1+theta_2)+m2*L_2^2*(dd_theta_1+dd_theta_2)];
n_1=[0;0;m2*L_1*L_2*cos(theta_2)*dd_theta_1+m2*L_1*L_2*sin(theta_2)*d_theta_1^2+m2*L_2*g*cos(theta_1+theta_2)+m2*L_2^2*(dd_theta_1+dd_theta_2)]+...
    [0;0;m1*L_1^2*dd_theta_1+m1*L_1*g*cos(theta_1)]+...
    [0;0;m2*L_1^2*dd_theta_1-m2*L_1*L_2*sin(theta_2)*(d_theta_1+d_theta_2)^2+m2*L_1*g*sin(theta_2)*sin(theta_1+theta_2)+m2*L_1*L_2*cos(theta_2)*(dd_theta_1+dd_theta_2)+m2*L_1*g*cos(theta_2)*cos(theta_1+theta_2)];
n_test=simplify(n_1-n(:,2));
n_test_2=simplify(n_2-n(:,3));

%% 关节力矩验证
% 关节力矩就是n的Z方向分量
tau_1 = n_1(3,1);
tau_2 = n_2(3,1);
tau_test = simplify(tau(:,2)-tau_1);
tau_test_2 = simplify(tau(:,3)-tau_2);