%% 符号表达式的尝试
syms t;
x = sin(t^2);
v = diff(x,t);%% 外推计算连杆速度与加速度

%% 符号表达式求导尝试
syms theta_1(t)
x = sin(theta_1(t));
res = diff(x,t);
res2 = diff(theta_1(t),t);
% 对导数替换，将diff(theta_1(t),t)这种未知的导数替换为d_theta_1
syms d_theta_1
res_clean = subs(res,diff(theta_1(t),t),d_theta_1);
res2_clean = subs(res2,diff(theta_1),d_theta_1);