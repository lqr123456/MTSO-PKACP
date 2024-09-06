function dataEMTEA = EMTEA(tasks,popSize,nGen,selPocess,pIL,nRepeat,dq,initPop,ntasks)
% EMTEA, Dongrui WU (drwu@hust.edu.cn), 8/20/2018

tic;

dataDisp=cell(1,3);
dataDisp{1}='PKACP'; dataDisp{2}='EMTEA';

mu = 2;     % Index of Simulated Binary Crossover (tunable)
mum = 5;    % Index of polynomial mutation
%nTasks=length(tasks);
nTasks=ntasks;
while mod(popSize,nTasks)
    popSize = popSize + 1;
end

dimTasks=zeros(1,nTasks); population=cell(1,nTasks);
for i=1:nTasks
    dimTasks(i) = tasks(i).dims;
end

callsPerIndividual=zeros(1,popSize);
bestFitness = zeros(nRepeat,nGen,nTasks);    % best fitness found
totalEvals=zeros(nRepeat,nGen);   % total number of task evaluations so far
options = optimoptions(@fminunc,'Display','off','Algorithm','quasi-newton','MaxIter',2);  % settings for individual learning
nTransfer=10;
fCosts=cell(1,2);
M=cell(nTasks,nTasks);
rnvecs=cell(1,2);

for r = 1:nRepeat
    dataDisp{3}=r;
    dq.send(dataDisp);
    
    % Initialize the first generation randomly
    gen=1;
    for idxTask=1:nTasks
        for i = 1 : popSize
            population{idxTask}(i) = Chromosome();
            population{idxTask}(i) = initialize(population{idxTask}(i),dimTasks(idxTask));
            if nargin>=9
                population{idxTask}(i).rnvec(1:dimTasks(idxTask))=initPop{idxTask,r}(i,1:dimTasks(idxTask));
            else
                initPop{idxTask,r}(i,1:dimTasks(idxTask))=population{idxTask}(i).rnvec;
            end
            [population{idxTask}(i),callsPerIndividual(i)] = evaluate_SOO(population{idxTask}(i),tasks(idxTask),pIL,options);
        end
        totalEvals(r,gen)=totalEvals(r,gen)+sum(callsPerIndividual);
        fCosts{idxTask}=[population{idxTask}.factorial_costs];
        [fCosts{idxTask},idsCost]=sort(fCosts{idxTask});
        population{idxTask}=population{idxTask}(idsCost); % sort the chromosomes according to their costs
        bestFitness(r,gen,idxTask)=fCosts{idxTask}(1);
        rnvecs{idxTask}=reshape([population{idxTask}.rnvec],dimTasks(idxTask),popSize)';
        rnvecs{idxTask}=cat(2,rnvecs{idxTask},ones(popSize,1));
    end
    % Compute M
    for i=1:nTasks-1
        for j=i+1:nTasks
            if dimTasks(i)>dimTasks(j)
                rnvecs{j}=cat(2,rnvecs{j},zeros(popSize,dimTasks(i)-dimTasks(j)));
            elseif dimTasks(i)<dimTasks(j)
                rnvecs{i}=cat(2,rnvecs{i},zeros(popSize,dimTasks(j)-dimTasks(i)));
            end
            M{i,j}=rnvecs{i}\rnvecs{j};
            M{j,i}=rnvecs{j}\rnvecs{i};
        end
    end
    
    for gen=2:nGen
        totalEvals(r,gen)=totalEvals(r,gen-1);
        for idxTask=1:nTasks
            prevTask=idxTask-1;
            if idxTask==1; prevTask=nTasks; end
            tempPopulation=population{idxTask}(end:-1:1);
            
            % Transfer some chromosomes from the previous task for reproduction
            if ~rem(gen,10)
                % replace bad chromosomes in the current population by good
                % chromosomes from the previous population
                if dimTasks(idxTask)>dimTasks(prevTask) % the previous task has a smaller dimensionality
                    for i=1:nTransfer
                        temp=[population{prevTask}(i).rnvec 1 zeros(1,dimTasks(idxTask)-dimTasks(prevTask))]...
                            *M{prevTask,idxTask};
                        tempPopulation(i).rnvec=temp(1:dimTasks(idxTask));
                    end
                else
                    for i=1:nTransfer
                        temp=[population{prevTask}(i).rnvec 1]*M{prevTask,idxTask};
                        tempPopulation(i).rnvec=temp(1:dimTasks(idxTask));
                    end
                end
            end
            
            idsOrder = randperm(popSize);
            for i = 1 : popSize/2
                p1 = idsOrder(i);
                p2 = idsOrder(i+popSize/2);
                u = rand(1,dimTasks(idxTask));
                cf = zeros(1,dimTasks(idxTask));
                cf(u<=0.5)=(2*u(u<=0.5)).^(1/(mu+1));
                cf(u>0.5)=(2*(1-u(u>0.5))).^(-1/(mu+1));
                child(2*i-1)=Chromosome();
                child(2*i)=Chromosome();
                child(2*i-1) = crossover(child(2*i-1),tempPopulation(p1),tempPopulation(p2),cf);
                child(2*i) = crossover(child(2*i),tempPopulation(p2),tempPopulation(p1),cf);
                child(2*i-1)=mutate(child(2*i-1),child(2*i-1),dimTasks(idxTask),mum);
                child(2*i)=mutate(child(2*i),child(2*i),dimTasks(idxTask),mum);
                %                 % variable swap
                %                 swap_indicator = (rand(1,dimTasks(idxTask)) >= 0.5);
                %                 temp = child(2*i-1).rnvec(swap_indicator);
                %                 child(2*i-1).rnvec(swap_indicator) = child(2*i).rnvec(swap_indicator);
                %                 child(2*i).rnvec(swap_indicator) = temp;
            end
            
            for i = 1 : popSize
                [child(i),callsPerIndividual(i)] = evaluate_SOO(child(i),tasks(idxTask),pIL,options);
            end
            totalEvals(r,gen)=totalEvals(r,gen)+sum(callsPerIndividual);
            intpopulation(1:popSize)=child;
            intpopulation(popSize+1:2*popSize)=population{idxTask};
            [~,idsCost]=sort([intpopulation.factorial_costs]);
            intpopulation=intpopulation(idsCost);
            bestFitness(r,gen,idxTask)=bestFitness(r,gen-1,idxTask);
            if intpopulation(1).factorial_costs<=bestFitness(r,gen,idxTask)
                bestFitness(r,gen,idxTask)=intpopulation(1).factorial_costs;
                bestChromosome(r,idxTask)=intpopulation(1);
            end
            
            if strcmp(selPocess,'elitist')
                if length(intpopulation)>=popSize
                    population{idxTask}=intpopulation(1:popSize);
                else
                    population{idxTask}(1:length(intpopulation))=intpopulation;
                    for i=length(intpopulation)+1:popSize
                        population{idxTask}(i) = Chromosome();
                        population{idxTask}(i) = initialize(population{idxTask}(i),dimTasks(idxTask));
                    end
                end
            elseif strcmp(selPocess,'roulette wheel')
                for i = 1:length(intpopulation)
                    intpopulation(i).scalar_fitness=1/i;
                end
                for i=1:popSize
                    population{idxTask}(i)=intpopulation(RouletteWheelSelection([intpopulation.scalar_fitness]));
                end
            end
        end
    end
end
dataEMTEA.wallClockTime=toc;
dataEMTEA.bestFitness=bestFitness;
dataEMTEA.bestIndData=bestChromosome;
dataEMTEA.totalEvals=totalEvals;
dataEMTEA.initPop=initPop;