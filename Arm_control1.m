function Tasks=Arm_control1(T,dim)
    
        samples = 50 * T;% 样本数量
        x = rand(samples, 2);% 生成随机样本
        [Idx, C] = kmeans(x, T); % K均值聚类
        task_para = C;% 任务参数
       
        for t = 1:T
            Tasks(t).dims = dim;% 设置每个任务的维度
            Tasks(t).Lb = zeros(1, Tasks(t).dims);%下界
            Tasks(t).Ub = ones(1, Tasks(t).dims);%上界
            Amax = task_para(t, 1);%最大角度
            Lmax = task_para(t, 2);%最大长度
            Tasks(t).fnc = @(x)fitness_arm(x, Amax, Lmax);% 每个任务的适应度函数
        end

end

function [Objs, Cons] = fitness_arm(angles_var, Amax, Lmax)% 计算运动臂控制问题的适应度和约束
Objs = [];
for i = 1:size(angles_var, 1)
    angles = angles_var(i, :); % 获取角度
    angular_range = Amax / length(angles);% 计算角度范围
    lengths = ones(1, length(angles)) * Lmax / length(angles); % 计算长度
    target = 0.5 * ones(1, 2); % 设置目标
    command = (angles - 0.5) * angular_range * pi * 2; % 计算指令
    ef = fw_kinematics(command, lengths);% 计算末端执行器位置
    fitness = sum((ef - target) .* (ef - target))^0.5;% 计算适应度
    Objs(i, :) = fitness;% 存储适应度
end
Cons = zeros(size(angles_var, 1), 1); % 无约束
end

function [joint_xy] = fw_kinematics(p, lengths)% 平面运动臂的正运动学
mat = eye(4);% 初始化单位矩阵
p = [p, 0];% 添加末端位置
n_dofs = length(p);% 自由度数量
joint_xy = zeros(1, 2);% 初始化关节坐标
lengths = [0, lengths];% 添加末端长度
for i = 1:n_dofs
    m = [cos(p(i)), -sin(p(i)), 0, lengths(i);
        sin(p(i)), cos(p(i)), 0, 0;
        0, 0, 1, 0;
        0, 0, 0, 1; ];% 计算变换矩阵
    mat = mat * m;% 累积变换矩阵
    v = mat * ([0, 0, 0, 1]');% 计算末端坐标
    joint_xy = v';% 更新关节坐标
end
joint_xy = joint_xy(1:2);% 去除末端坐标的Z轴分量
end