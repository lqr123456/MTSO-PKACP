function Tasks=Arm_control1(T,dim)
    
        samples = 50 * T;
        x = rand(samples, 2);
        [Idx, C] = kmeans(x, T); 
        task_para = C;
       
        for t = 1:T
            Tasks(t).dims = dim;
            Tasks(t).Lb = zeros(1, Tasks(t).dims);
            Tasks(t).Ub = ones(1, Tasks(t).dims);
            Amax = task_para(t, 1);
            Lmax = task_para(t, 2);
            Tasks(t).fnc = @(x)fitness_arm(x, Amax, Lmax);
        end

end

function [Objs, Cons] = fitness_arm(angles_var, Amax, Lmax)
Objs = [];
for i = 1:size(angles_var, 1)
    angles = angles_var(i, :); 
    angular_range = Amax / length(angles);
    lengths = ones(1, length(angles)) * Lmax / length(angles); 
    target = 0.5 * ones(1, 2); 
    command = (angles - 0.5) * angular_range * pi * 2; 
    ef = fw_kinematics(command, lengths);
    fitness = sum((ef - target) .* (ef - target))^0.5;
    Objs(i, :) = fitness;
end
Cons = zeros(size(angles_var, 1), 1); 
end

function [joint_xy] = fw_kinematics(p, lengths)
mat = eye(4);
p = [p, 0];
n_dofs = length(p);
joint_xy = zeros(1, 2);
lengths = [0, lengths];
for i = 1:n_dofs
    m = [cos(p(i)), -sin(p(i)), 0, lengths(i);
        sin(p(i)), cos(p(i)), 0, 0;
        0, 0, 1, 0;
        0, 0, 0, 1; ];
    mat = mat * m;
    v = mat * ([0, 0, 0, 1]');
    joint_xy = v';
end
joint_xy = joint_xy(1:2);
end
