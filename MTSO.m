function task = MTSO(Tasks, N, T, rmp,nRepeat,ntasks)

    if mod(N, 2) ~= 0 
        N = N + 1;
    end

    no_of_tasks = ntasks; 
    if no_of_tasks <= 1 
        error('At least 2 tasks required for MFEA');
    end

    
    D = zeros(1, no_of_tasks);
    for i = 1:no_of_tasks
        D(i) = Tasks.D(i);
    end
    D_multitask = max(D);

    %options = optimoptions(@fminunc, 'Display', 'off', 'Algorithm', 'quasi-newton', 'MaxIter', 2); % settings for individual learning
    
    %bestobj = Inf(1, no_of_tasks);
 for rep = 1:nRepeat
    %initial 
    tic
    vec_flag=[1,-1];
    Threshold=0.23;
    Thresold2= 0.45;
    C1=0.5;
    C2=.05;
    C3=2;
    Z=3000;
    CLS=ch(1,Z);
    HMCR=0.95;
    PAR=0.4;
    u=0.005;
    VarSize=[1 D_multitask]; 
    pCR=0.2; 
    beta_min=0.2;   % Lower Bound of Scaling Factor
    beta_max=0.8;
    
for k=1:no_of_tasks
    
    Task(k).fobj=Tasks.Fnc{k};
    UB(k,:)=Tasks.Ub{k}(1:D_multitask);
    LB(k,:)=Tasks.Lb{k}(1:D_multitask);
    range(k,:)=UB(k,:)-LB(k,:); 

    Task(k).X = repmat(LB(k,:),N,1)+rand(N,D_multitask).* repmat((UB(k,:)-LB(k,:)),N,1);  
    for i=1:N
        %Task(k).fitness(i)=feval(Task(k).fobj,Task(k).X(i,:),Tasks.Amax(k), Tasks.Lmax(k));  
        Task(k).fitness(i)=Task(k).fobj(Task(k).X(i,:));
    end
    [GYbest(k), gbest(k)] = min(Task(k).fitness);
    Task(k).Xfood = Task(k).X(gbest(k),:);
    %Diving the swarm into two equal groups males and females
    Nm=round(N/2);%eq.(2&3)
    Nf=N-Nm;
    Task(k).Xm=Task(k).X(1:Nm,:);
    Task(k).Xf=Task(k).X(Nm+1:N,:);
    Task(k).fitness_m=Task(k).fitness(1:Nm);
    Task(k).fitness_f=Task(k).fitness(Nm+1:N);
    [Task(k).fitnessBest_m, Task(k).gbest1] = min(Task(k).fitness_m);
    Task(k).Xbest_m = Task(k).Xm(Task(k).gbest1,:);
    [Task(k).fitnessBest_f, Task(k).gbest2] = min(Task(k).fitness_f);
    Task(k).Xbest_f = Task(k).Xf(Task(k).gbest2,:);
    
%     [xxx, y] = sort(Task(k).fitness); 
%     Task(k).X = Task(k).X(y); 
%     for j = 1:N
%         Task(k).X(j).factorial_ranks(k) = j; 
%     end
%     bestobj(k) = xxx(1); 
end


for t = 1:T
    for k=1:no_of_tasks
        Temp=exp(-((t)/T));  %eq.(4)
        Q=C1*exp(((t-T)/(T)));%eq.(5)
        if Q>1        Q=1;    end
        % Exploration Phase (no Food)
        if Q<Threshold
            for i=1:Nm
                
                    rand_leader_index = floor(Nm*rand()+1);
                    Task(k).X_randm = Task(k).Xm(rand_leader_index, :);
                    flag_index = floor(2*rand()+1);
                    Flag=vec_flag(flag_index);
                    Am=exp(-Task(k).fitness_m(rand_leader_index)/(Task(k).fitness_m(i)+eps));%eq.(7)
                    Task(k).Xnewm(i,:)=Task(k).X_randm+Flag*C2*Am*((UB(k,:)-LB(k,:))*rand+LB(k,:));%eq.(6)
                
            end
            for i=1:Nf
                %for j=1:1:D_multitask
                    
                    rand_leader_index = floor(Nf*rand()+1);
                    Task(k).X_randf = Task(k).Xf(rand_leader_index, :);
                    flag_index = floor(2*rand()+1);
                    Flag=vec_flag(flag_index);
                    Af=exp(-Task(k).fitness_f(rand_leader_index)/(Task(k).fitness_f(i)+eps));%eq.(9)
                    Task(k).Xnewf(i,:)=Task(k).X_randf+Flag*C2*Af*((UB(k,:)-LB(k,:))*rand+LB(k,:));%eq.(8)
                    
                %end
            end
            
        else %Exploitation Phase (Food Exists)
            if Temp>Thresold2  %hot
                for i=1:Nm
                    flag_index = floor(2*rand()+1);
                    Flag=vec_flag(flag_index);
                    for j=1:1:D_multitask
                        Task(k).Xnewm(i,j)=Task(k).Xfood(j)+C3*Flag*Temp*rand*(Task(k).Xfood(j)-Task(k).Xm(i,j));%eq.(10)
                    end
                end
                for i=1:Nf
                    flag_index = floor(2*rand()+1);
                    Flag=vec_flag(flag_index);
                    for j=1:1:D_multitask
                        Task(k).Xnewf(i,j)=Task(k).Xfood(j)+Flag*C3*Temp*rand*(Task(k).Xfood(j)-Task(k).Xf(i,j));%eq.(10)
                    end
                end
                
            else %cold
                if rand>0.6 %fight
                    for i=1:Nm
                        for j=1:1:D_multitask
                            FM=exp(-(Task(k).fitnessBest_f)/(Task(k).fitness_m(i)+eps));%eq.(13)
                            Task(k).Xnewm(i,j)=Task(k).Xm(i,j) +C3*FM*rand*(Q*Task(k).Xbest_f(j)-Task(k).Xm(i,j));%eq.(11)
                        end
                    end
                    for i=1:Nf
                        for j=1:1:D_multitask
                            FF=exp(-(Task(k).fitnessBest_m)/(Task(k).fitness_f(i)+eps));%eq.(14)
                            Task(k).Xnewf(i,j)=Task(k).Xf(i,j)+C3*FF*rand*(Q*Task(k).Xbest_m(j)-Task(k).Xf(i,j));%eq.(12)
                        end
                    end
                    
                else%mating
                    for i=1:Nm
                        for j=1:1:D_multitask
                            Mm=exp(-Task(k).fitness_f(i)/(Task(k).fitness_m(i)+eps));%eq.(17)
                            Task(k).Xnewm(i,j)=Task(k).Xm(i,j) +C3*rand*Mm*(Q*Task(k).Xf(i,j)-Task(k).Xm(i,j));%eq.(15
                        end
                    end
                    for i=1:Nf
                        for j=1:1:D_multitask
                            Mf=exp(-Task(k).fitness_m(i)/(Task(k).fitness_f(i)+eps));%eq.(18)
                            Task(k).Xnewf(i,j)=Task(k).Xf(i,j) +C3*rand*Mf*(Q*Task(k).Xm(i,j)-Task(k).Xf(i,j));%eq.(16)
                        end
                    end
                    
                    flag_index = floor(2*rand()+1);
                    egg=vec_flag(flag_index);
                    if egg==1;
                        [Task(k).GYworst, Task(k).gworst] = max(Task(k).fitness_m);
                        
                        rand_leader_index1 = floor(Nm*rand()+1);
                        rand_leader_index2 = floor(Nm*rand()+1);
                        X_randm1 = Task(k).Xm(rand_leader_index1, :);
                        X_randm2 = Task(k).Xm(rand_leader_index2, :);
                        beta=unifrnd(beta_min,beta_max,VarSize);
                        m=Task(k).Xbest_m+beta.*(X_randm1-X_randm2);
%                 rand_leader_index = floor(Nm*rand()+1);
%                 X_randm = Xm(rand_leader_index, :);
                        z=zeros(size(X_randm1));
                        j0=randi([1 D_multitask]);
                        for j=1:D_multitask
                            if j==j0 || rand<=pCR
                                z(j)=Task(k).Xbest_m(j);
                            else
                                z(j)=m(j);
                            end
                        end
                        Task(k).Xnewm(Task(k).gworst,:)=z;%eq.(19)
                        %Task(k).Xnewm(Task(k).gworst,:)=repmat(LB(k,:),1,1)+rand(1,D_multitask).* repmat((UB(k,:)-LB(k,:)),1,1); %eq.(19)
                        
                        [Task(k).GYworst, Task(k).gworst] = max(Task(k).fitness_f);
                        
                        rand_leader_index1 = floor(Nf*rand()+1);
                        rand_leader_index2 = floor(Nf*rand()+1);
                        X_randf1 = Task(k).Xf(rand_leader_index1, :);
                        X_randf2 = Task(k).Xf(rand_leader_index2, :);
                        beta=unifrnd(beta_min,beta_max,VarSize);
                        m=Task(k).Xbest_f+beta.*(X_randf1-X_randf2);
%                 rand_leader_index = floor(Nf*rand()+1);
%                 X_randf = Xf(rand_leader_index, :);
                        a=zeros(size(X_randf1));
                        j0=randi([1 D_multitask]);
                        for j=1:D_multitask
                            if j==j0 || rand<=pCR
                                a(j)=Task(k).Xbest_f(j);
                            else
                                a(j)=m(j);
                            end
                        end
                        Task(k).Xnewf(Task(k).gworst,:)=a;%eq.(20)
                        %Task(k).Xnewf(Task(k).gworst,:)=repmat(LB(k,:),1,1)+rand(1,D_multitask).* repmat((UB(k,:)-LB(k,:)),1,1);%eq.(20)
                    end
                end
            end
        end
        
        for j=1:Nm
            Flag4ub=Task(k).Xnewm(j,:)>UB(k);
            Flag4lb=Task(k).Xnewm(j,:)<LB(k);
            Task(k).Xnewm(j,:)=(Task(k).Xnewm(j,:).*(~(Flag4ub+Flag4lb)))+UB(k,:).*Flag4ub+LB(k,:).*Flag4lb;
            y = feval(Task(k).fobj,Task(k).Xnewm(j,:));
            if y<Task(k).fitness_m(j)
                Task(k).fitness_m(j)=y;
                Task(k).Xm(j,:)= Task(k).Xnewm(j,:);
            end
        end
        [Task(k).Ybest1,Task(k).gbest1] = min(Task(k).fitness_m);
    
        for j=1:Nf
            Flag4ub=Task(k).Xnewf(j,:)>UB(k);
            Flag4lb=Task(k).Xnewf(j,:)<LB(k);
            Task(k).Xnewf(j,:)=(Task(k).Xnewf(j,:).*(~(Flag4ub+Flag4lb)))+UB(k,:).*Flag4ub+LB(k,:).*Flag4lb;
            y = feval(Task(k).fobj,Task(k).Xnewf(j,:));
            if y<Task(k).fitness_f(j)
                Task(k).fitness_f(j)=y;
                Task(k).Xf(j,:)= Task(k).Xnewf(j,:);
            end
        end
        [Task(k).Ybest2,Task(k).gbest2] = min(Task(k).fitness_f);
    
        if Task(k).Ybest1<Task(k).fitnessBest_m
            Task(k).Xbest_m = Task(k).Xm(Task(k).gbest1,:);
            Task(k).fitnessBest_m=Task(k).Ybest1;
        end
        if Task(k).Ybest2<Task(k).fitnessBest_f
            Task(k).Xbest_f = Task(k).Xf(Task(k).gbest2,:);
            Task(k).fitnessBest_f=Task(k).Ybest2;
        end
    
        if Task(k).Ybest1<Task(k).Ybest2
            Task(k).gbest_t(t)=min(Task(k).Ybest1);
        else
        Task(k).gbest_t(t)=min(Task(k).Ybest2);
        end
        
        if Task(k).fitnessBest_m<Task(k).fitnessBest_f
            Task(k).GYbest=Task(k).fitnessBest_m;
            Task(k).Xfood=Task(k).Xbest_m;
        else
            Task(k).GYbest=Task(k).fitnessBest_f;
            Task(k).Xfood=Task(k).Xbest_f;
        end
        
        r=5*(1-t/T)^2;
        Task(k).XX=[Task(k).Xm;Task(k).Xf];
        rand_index1 = floor(N*rand()+1);
        rand_index2= floor(N*rand()+1);
        XX_rand1 = Task(k).XX(rand_index1, :);
        XX_rand2 = Task(k).XX(rand_index2, :);
        XFood=Task(k).Xfood+r*CLS(t)*(XX_rand1-XX_rand2);
        
        for z=1:D_multitask
            Flag4ub=XFood(z)>UB(k,z);
            Flag4lb=XFood(z)<LB(k,z);
            XFood(z)=(XFood(z).*(~(Flag4ub+Flag4lb)))+UB(k,z).*Flag4ub+LB(k,z).*Flag4lb;
        end

        % Flag4ub=XFood>UB(k,:);
        % Flag4lb=XFood<LB(k,:);
        % XFood=repmat(LB(k,:),1,1)+rand(1,D_multitask).* repmat((UB(k,:)-LB(k,:)),1,1);
        

        %XFood=Xfood+r*CLS(floor(Z*rand+1))*(XX_rand1-XX_rand2);
        F=feval(Task(k).fobj,XFood);
        if F<Task(k).GYbest
            if Task(k).fitnessBest_m<Task(k).fitnessBest_f
                Task(k).Xbest_m=XFood;
                Task(k).Xm(Task(k).gbest1,:)=XFood;
                Task(k).fitnessBest_m=F;
                Task(k).Ybest1=F;
                Task(k).fitness_m(Task(k).gbest1)=F;
            else
                Task(k).Xbest_f=XFood;
                Task(k).Xf(Task(k).gbest2,:)=XFood;
                Task(k).fitnessBest_f=F;
                Task(k).Ybest2=F;
                Task(k).fitness_f(Task(k).gbest2)=F;
            end
            Task(k).Xfood=XFood;
            Task(k).GYbest=F;
        end
        Task(k).gbest_t(t)=Task(k).GYbest;
        
        
    end
    
    
    n=N/5;
    d=D_multitask;
    RMP=rand();
    
    for k=1:no_of_tasks
        Task(k).XX=[Task(k).Xm;Task(k).Xf]; 
        Task(k).fitness=[Task(k).fitness_m,Task(k).fitness_f];
        [x,y]=sort(Task(k).fitness);
        Task(k).fitness1=x;
        Task(k).XX1=Task(k).XX(y,:);
        Task(k).XX_elite=Task(k).XX1(1:n,:);
        Task(k).fitness_elite=Task(k).fitness1(1:n);
        
       
    
        Task(k).XXX=(Task(k).XX_elite - LB(k,:))./range(k,:);
    end
    
    for k=1:no_of_tasks
        if RMP<rmp
            if rand<HMCR
                
                k0 = ceil(rand*no_of_tasks);
                while k0 == k
                    k0 = ceil(rand*no_of_tasks);
                end
                a = ceil(rand*n);
                b = ceil(rand*d);
                Task(k).XXX(a,b) = Task(k0).XXX(a,b);
                Task(k).XX_elite1=(Task(k).XXX.*range(k,:))+LB(k,:);
                for i=1:n

                    for z=1:D_multitask
                        Flag4ub=Task(k).XX_elite1(i,z)>UB(k,z);
                        Flag4lb=Task(k).XX_elite1(i,z)<LB(k,z);
                        Task(k).XX_elite1(i,z)=(Task(k).XX_elite1(i,z).*(~(Flag4ub+Flag4lb)))+UB(k,z).*Flag4ub+LB(k,z).*Flag4lb;
                    end

                    %if Task(k).XX_elite1(i,:)>UB(k,:)||Task(k).XX_elite1(i,:)<LB(k,:)
                    % Flag4ub=Task(k).XX_elite1(i,:)>UB(k,:);
                    % Flag4lb=Task(k).XX_elite1(i,:)<LB(k,:);
                    % Task(k).XX_elite1(i,:)=repmat(LB(k,:),1,1)+rand(1,D_multitask).* repmat((UB(k,:)-LB(k,:)),1,1);
                    %end

                    f(i)=feval(Task(k).fobj,Task(k).XX_elite1(i,:));   
                end
                
                Task(k).fit=[Task(k).fitness_elite,f];
                Task(k).XX_elite2=[Task(k).XX_elite;Task(k).XX_elite1];
                
                [fitk,index]=sort(Task(k).fit);
                Task(k).fit1=fitk;
                Task(k).XX_elite3=Task(k).XX_elite2(index,:);
                Task(k).XX_elite=Task(k).XX_elite3(1:n,:);
                
                Task(k).fitness1(1:n)=Task(k).fit1(1:n);
                Task(k).XX1(1:n,:)=Task(k).XX_elite;

                
                Task(k).XX(y,:)=Task(k).XX1;
                Task(k).fitness(y)=Task(k).fitness1;
            else
                if rand<PAR
                    [Task(k).GYworst2, Task(k).gworst2] = max(Task(k).fitness);
                    Task(k).XXnew= Task(k).XX(Task(k).gworst2,:)+2*u*rand-u;
                    
                    for z=1:D_multitask
                        Flag4ub=Task(k).XXnew(z)>UB(k,z);
                        Flag4lb=Task(k).XXnew(z)<LB(k,z);
                        Task(k).XXnew(z)=(Task(k).XXnew(z).*(~(Flag4ub+Flag4lb)))+UB(k,z).*Flag4ub+LB(k,z).*Flag4lb;
                    end

                    %if Task(k).XXnew>UB(k,:)||Task(k).XXnew<LB(k,:)
                        % Flag4ub=Task(k).XXnew>UB(k,:);
                        % Flag4lb=Task(k).XXnew<LB(k,:);
                        % Task(k).XXnew=repmat(LB(k,:),1,1)+rand(1,D_multitask).* repmat((UB(k,:)-LB(k,:)),1,1);
                    %end

                    f=feval(Task(k).fobj,Task(k).XXnew);
                    if f<Task(k).GYworst2
                        Task(k).XX(Task(k).gworst2,:)=Task(k).XXnew;
                        Task(k).fitness(Task(k).gworst2) = f;
                    end
                end
            end
        else
            %%LOBL strategy
            h=(1+(t/T)^0.5)^10;
            [Task(k).GYworst2, Task(k).gworst2] = max(Task(k).fitness);
            Task(k).Xnew = (UB(k,:)+LB(k,:))/2+(UB(k,:)+LB(k,:))/(2*h)-Task(k).XX(Task(k).gworst2,:)/h;
            %update positions

            for z=1:D_multitask
                flag4ub = Task(k).Xnew(z)>UB(k,z);
                flag4lb = Task(k).Xnew(z)<LB(k,z);
                Task(k).Xnew(z) = (Task(k).Xnew(z).*(~(flag4ub+flag4ub)))+LB(k,z).*flag4lb+UB(k,z).*flag4ub;
            end



            y_new = feval(Task(k).fobj,Task(k).Xnew);
            if y_new<Task(k).GYworst2
                Task(k).XX(Task(k).gworst2,:) = Task(k).Xnew;
                Task(k).fitness(Task(k).gworst2) = y_new;
            end
            
        end
     
    [GYbest(k), gbest(k)] = min(Task(k).fitness);
    Task(k).Xfood = Task(k).XX(gbest(k),:);
    %Diving the swarm into two equal groups males and females
    Nm=round(N/2);%eq.(2&3)
    Nf=N-Nm;
    Task(k).Xm=Task(k).XX(1:Nm,:);
    Task(k).Xf=Task(k).XX(Nm+1:N,:);
    Task(k).fitness_m=Task(k).fitness(1:Nm);
    Task(k).fitness_f=Task(k).fitness(Nm+1:N);
    [Task(k).fitnessBest_m, Task(k).gbest1] = min(Task(k).fitness_m);
    Task(k).Xbest_m = Task(k).Xm(Task(k).gbest1,:);
    [Task(k).fitnessBest_f, Task(k).gbest2] = min(Task(k).fitness_f);
    Task(k).Xbest_f = Task(k).Xf(Task(k).gbest2,:);
    
    Task(k).GYbest=GYbest(k);
    Task(k).gbest_t(t)=Task(k).GYbest;
    
    %evBestFitness = zeros(no_of_tasks*nRepeat,T); 
    evBestFitness(k+no_of_tasks*(rep-1),t)=Task(k).GYbest;
    bestIndData{rep,k}=Task(k).Xfood;
    end
    
    
    
end

% for k=1:no_of_tasks
%     task=struct('Task',Task);
%     task(nRepeat).Task(k).GYbest = Task(k).GYbest;  
%     task(nRepeat).Task(k).Xfood =Task(k).Xfood;
%     task(nRepeat).Task(k).Gbest_t=Task(k).gbest_t;
% end

 end
% task.wallClockTime=toc;
task.wallClockTime=toc;
task.bestFitness=evBestFitness;
task.bestIndData=bestIndData;
%task.totalEvals=TotalEvaluations;

end
