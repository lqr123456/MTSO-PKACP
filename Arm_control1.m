function Tasks=Arm_control1(T,dim)
    
        samples = 50 * T;% ��������
        x = rand(samples, 2);% �����������
        [Idx, C] = kmeans(x, T); % K��ֵ����
        task_para = C;% �������
       
        for t = 1:T
            Tasks(t).dims = dim;% ����ÿ�������ά��
            Tasks(t).Lb = zeros(1, Tasks(t).dims);%�½�
            Tasks(t).Ub = ones(1, Tasks(t).dims);%�Ͻ�
            Amax = task_para(t, 1);%���Ƕ�
            Lmax = task_para(t, 2);%��󳤶�
            Tasks(t).fnc = @(x)fitness_arm(x, Amax, Lmax);% ÿ���������Ӧ�Ⱥ���
        end

end

function [Objs, Cons] = fitness_arm(angles_var, Amax, Lmax)% �����˶��ۿ����������Ӧ�Ⱥ�Լ��
Objs = [];
for i = 1:size(angles_var, 1)
    angles = angles_var(i, :); % ��ȡ�Ƕ�
    angular_range = Amax / length(angles);% ����Ƕȷ�Χ
    lengths = ones(1, length(angles)) * Lmax / length(angles); % ���㳤��
    target = 0.5 * ones(1, 2); % ����Ŀ��
    command = (angles - 0.5) * angular_range * pi * 2; % ����ָ��
    ef = fw_kinematics(command, lengths);% ����ĩ��ִ����λ��
    fitness = sum((ef - target) .* (ef - target))^0.5;% ������Ӧ��
    Objs(i, :) = fitness;% �洢��Ӧ��
end
Cons = zeros(size(angles_var, 1), 1); % ��Լ��
end

function [joint_xy] = fw_kinematics(p, lengths)% ƽ���˶��۵����˶�ѧ
mat = eye(4);% ��ʼ����λ����
p = [p, 0];% ���ĩ��λ��
n_dofs = length(p);% ���ɶ�����
joint_xy = zeros(1, 2);% ��ʼ���ؽ�����
lengths = [0, lengths];% ���ĩ�˳���
for i = 1:n_dofs
    m = [cos(p(i)), -sin(p(i)), 0, lengths(i);
        sin(p(i)), cos(p(i)), 0, 0;
        0, 0, 1, 0;
        0, 0, 0, 1; ];% ����任����
    mat = mat * m;% �ۻ��任����
    v = mat * ([0, 0, 0, 1]');% ����ĩ������
    joint_xy = v';% ���¹ؽ�����
end
joint_xy = joint_xy(1:2);% ȥ��ĩ�������Z�����
end