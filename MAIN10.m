clc; clearvars; close all; rng('default'); %warning off all;
popSize=200; % Population size
PopSize=100;
nGen=500; % Number of generations
selProcess = 'elitist'; % Choose either 'elitist' or 'roulette wheel
pIL = 0; % robability of individual learning (BFGA quasi-Newton Algorithm) --> Indiviudal Learning is an IMPORTANT component of the MFEA.
rmp=0.3; % Random mating probability
nRepeat = 10; % Number of repeats; should be 30 in final submission; must >1 to avoid no-display problem
pTransfer=0.4; % Portion of chromosomes to transfer from one task to another
eMin=0.01; % Threshold of accumulated survival rate of divergents
Rmp=0.5;
dim = 20;
ntasks=5;

dqWorker = parallel.pool.DataQueue; %dqWorker 
afterEach(dqWorker, @(data) fprintf('%d-%s-%d ', data{1},data{2},data{3})); % print progress of parfor
dqClient=parallel.pool.DataQueue; afterEach(dqClient,@showResults);

tasks = Arm_control(ntasks,dim);
Tasks = Arm_control1(ntasks,dim);

initPop=cell(10,nRepeat);
SOEA1=SOEA(Tasks(1),popSize,nGen,selProcess,pIL,nRepeat,dqWorker);
SOEA2=SOEA(Tasks(2),popSize,nGen,selProcess,pIL,nRepeat,dqWorker);
SOEA3=SOEA(Tasks(3),popSize,nGen,selProcess,pIL,nRepeat,dqWorker);
SOEA4=SOEA(Tasks(4),popSize,nGen,selProcess,pIL,nRepeat,dqWorker);
SOEA5=SOEA(Tasks(5),popSize,nGen,selProcess,pIL,nRepeat,dqWorker);
% SOEA6=SOEA(Tasks(6),popSize,nGen,selProcess,pIL,nRepeat,dqWorker);
% SOEA7=SOEA(Tasks(7),popSize,nGen,selProcess,pIL,nRepeat,dqWorker);
% SOEA8=SOEA(Tasks(8),popSize,nGen,selProcess,pIL,nRepeat,dqWorker);
% SOEA9=SOEA(Tasks(9),popSize,nGen,selProcess,pIL,nRepeat,dqWorker);
% SOEA10=SOEA(Tasks(10),popSize,nGen,selProcess,pIL,nRepeat,dqWorker);

% SOEA11=SOEA(Tasks(11),popSize,nGen,selProcess,pIL,nRepeat,dqWorker);
% SOEA12=SOEA(Tasks(12),popSize,nGen,selProcess,pIL,nRepeat,dqWorker);
% SOEA13=SOEA(Tasks(13),popSize,nGen,selProcess,pIL,nRepeat,dqWorker);
% SOEA14=SOEA(Tasks(14),popSize,nGen,selProcess,pIL,nRepeat,dqWorker);
% SOEA15=SOEA(Tasks(15),popSize,nGen,selProcess,pIL,nRepeat,dqWorker);
% SOEA16=SOEA(Tasks(16),popSize,nGen,selProcess,pIL,nRepeat,dqWorker);
% SOEA17=SOEA(Tasks(17),popSize,nGen,selProcess,pIL,nRepeat,dqWorker);
% SOEA18=SOEA(Tasks(18),popSize,nGen,selProcess,pIL,nRepeat,dqWorker);
% SOEA19=SOEA(Tasks(19),popSize,nGen,selProcess,pIL,nRepeat,dqWorker);
% SOEA20=SOEA(Tasks(20),popSize,nGen,selProcess,pIL,nRepeat,dqWorker);

    for r=1:nRepeat
        initPop{1,r}=SOEA1.initPop{r}; initPop{2,r}=SOEA2.initPop{r}; % initial population
        initPop{3,r}=SOEA3.initPop{r}; initPop{4,r}=SOEA4.initPop{r};
        initPop{5,r}=SOEA5.initPop{r}; %initPop{6,r}=SOEA6.initPop{r};
        % initPop{7,r}=SOEA7.initPop{r}; initPop{8,r}=SOEA8.initPop{r};
        % initPop{9,r}=SOEA9.initPop{r}; initPop{10,r}=SOEA10.initPop{r};
        
%         initPop{11,r}=SOEA11.initPop{r}; initPop{12,r}=SOEA12.initPop{r}; % initial population
%         initPop{13,r}=SOEA13.initPop{r}; initPop{14,r}=SOEA14.initPop{r};
%         initPop{15,r}=SOEA15.initPop{r}; initPop{16,r}=SOEA16.initPop{r};
%         initPop{17,r}=SOEA17.initPop{r}; initPop{18,r}=SOEA18.initPop{r};
%         initPop{19,r}=SOEA19.initPop{r}; initPop{20,r}=SOEA20.initPop{r};
    end

%%
data2=MFEA(Tasks,popSize,nGen,selProcess,rmp,pIL,nRepeat,dqWorker,initPop,ntasks); % The provided MFEA benchmark algorithm
Mean_MFEA=mean(data2.bestFitness(1:end,end))

S1=mean(data2.bestFitness(1:5,end)); S2=mean(data2.bestFitness(6:10,end));
S3=mean(data2.bestFitness(11:15,end)); S4=mean(data2.bestFitness(16:20,end));
S5=mean(data2.bestFitness(21:25,end));S=[S1,S2,S3,S4,S5];
% S6=mean(data2.bestFitness(51:60,end));
% S7=mean(data2.bestFitness(61:70,end)); S8=mean(data2.bestFitness(71:80,end));
% S9=mean(data2.bestFitness(81:90,end)); S10=mean(data2.bestFitness(91:100,end));
%S=[S1,S2,S3,S4,S5,S6,S7,S8,S9,S10];

STD_MFEA=std(S)

%%
data3=MFEARR(Tasks,popSize,nGen,selProcess,rmp,pIL,nRepeat,dqWorker,eMin,initPop,ntasks);
Mean_MFEARR=mean(data3.bestFitness(1:end,end))

S1=mean(data3.bestFitness(1:5,end)); S2=mean(data3.bestFitness(6:10,end));
S3=mean(data3.bestFitness(11:15,end)); S4=mean(data3.bestFitness(16:20,end));
S5=mean(data3.bestFitness(21:25,end));S=[S1,S2,S3,S4,S5];

% S1=mean(data3.bestFitness(1:10,end)); S2=mean(data3.bestFitness(11:20,end));
% S3=mean(data3.bestFitness(21:30,end)); S4=mean(data3.bestFitness(31:40,end));
% S5=mean(data3.bestFitness(41:50,end)); S6=mean(data3.bestFitness(51:60,end));
% S7=mean(data3.bestFitness(61:70,end)); S8=mean(data3.bestFitness(71:80,end));
% S9=mean(data3.bestFitness(81:90,end)); S10=mean(data3.bestFitness(91:100,end));
% S=[S1,S2,S3,S4,S5,S6,S7,S8,S9,S10];

STD_MFEARR=std(S)
%%
data6=EBSGA(Tasks,PopSize,nGen,selProcess,rmp,pIL,nRepeat,dqWorker,initPop,ntasks);
A6=squeeze(mean(data6.bestFitness,1));
Mean_EBSGA=mean(A6(end,:))

S1=mean(data6.bestFitness(1,end,1:5)); S2=mean(data6.bestFitness(2,end,1:5));
S3=mean(data6.bestFitness(3,end,1:5)); S4=mean(data6.bestFitness(4,end,1:5));
S5=mean(data6.bestFitness(5,end,1:5));S=[S1,S2,S3,S4,S5];

% S1=mean(data6.bestFitness(1,end,1:10)); S2=mean(data6.bestFitness(2,end,1:10));
% S3=mean(data6.bestFitness(3,end,1:10)); S4=mean(data6.bestFitness(4,end,1:10));
% S5=mean(data6.bestFitness(5,end,1:10)); S6=mean(data6.bestFitness(6,end,1:10));
% S7=mean(data6.bestFitness(7,end,1:10)); S8=mean(data6.bestFitness(8,end,1:10));
% S9=mean(data6.bestFitness(9,end,1:10)); S10=mean(data6.bestFitness(10,end,1:10));
% S=[S1,S2,S3,S4,S5,S6,S7,S8,S9,S10];
%S=[S1,S2];
STD_EBSGA=std(S)

%%
data7=GMFEA(Tasks,popSize,nGen,selProcess,rmp,pIL,nRepeat,dqWorker,initPop,ntasks);
Mean_GMFEA=mean(data7.bestFitness(1:end,end))

S1=mean(data7.bestFitness(1:5,end)); S2=mean(data7.bestFitness(6:10,end));
S3=mean(data7.bestFitness(11:15,end)); S4=mean(data7.bestFitness(16:20,end));
S5=mean(data7.bestFitness(21:25,end));S=[S1,S2,S3,S4,S5];

% S1=mean(data7.bestFitness(1:10,end)); S2=mean(data7.bestFitness(11:20,end));
% S3=mean(data7.bestFitness(21:30,end)); S4=mean(data7.bestFitness(31:40,end));
% S5=mean(data7.bestFitness(41:50,end)); S6=mean(data7.bestFitness(51:60,end));
% S7=mean(data7.bestFitness(61:70,end)); S8=mean(data7.bestFitness(71:80,end));
% S9=mean(data7.bestFitness(81:90,end)); S10=mean(data7.bestFitness(91:100,end));
% S=[S1,S2,S3,S4,S5,S6,S7,S8,S9,S10];
%S=[S1,S2];
STD_GMFEA=std(S)

%%
data8=EMTEA(Tasks,PopSize,nGen,selProcess,pIL,nRepeat,dqWorker,initPop,ntasks);
A8=squeeze(mean(data8.bestFitness,1));
Mean_EMFEA=mean(A8(end,:))

S1=mean(data8.bestFitness(1,end,1:5)); S2=mean(data8.bestFitness(2,end,1:5));
S3=mean(data8.bestFitness(3,end,1:5)); S4=mean(data8.bestFitness(4,end,1:5));
S5=mean(data8.bestFitness(5,end,1:5));S=[S1,S2,S3,S4,S5];

% S1=mean(data8.bestFitness(1,end,1:10)); S2=mean(data8.bestFitness(2,end,1:10));
% S3=mean(data8.bestFitness(3,end,1:10)); S4=mean(data8.bestFitness(4,end,1:10));
% S5=mean(data8.bestFitness(5,end,1:10)); S6=mean(data8.bestFitness(6,end,1:10));
% S7=mean(data8.bestFitness(7,end,1:10)); S8=mean(data8.bestFitness(8,end,1:10));
% S9=mean(data8.bestFitness(9,end,1:10)); S10=mean(data8.bestFitness(10,end,1:10));
% S=[S1,S2,S3,S4,S5,S6,S7,S8,S9,S10];

STD_EMFEA=std(S)

%%
data9=MTEAbest(Tasks,PopSize,nGen,selProcess,pIL,nRepeat,pTransfer,dqWorker,initPop,ntasks); % Our algorithm, v2
A9=squeeze(mean(data9.bestFitness,1));
Mean_MFEAbest=mean(A9(end,:))

S1=mean(data9.bestFitness(1,end,1:5)); S2=mean(data9.bestFitness(2,end,1:5));
S3=mean(data9.bestFitness(3,end,1:5)); S4=mean(data9.bestFitness(4,end,1:5));
S5=mean(data9.bestFitness(5,end,1:5));S=[S1,S2,S3,S4,S5];

% S1=mean(data9.bestFitness(1,end,1:10)); S2=mean(data9.bestFitness(2,end,1:10));
% S3=mean(data9.bestFitness(3,end,1:10)); S4=mean(data9.bestFitness(4,end,1:10));
% S5=mean(data9.bestFitness(5,end,1:10)); S6=mean(data9.bestFitness(6,end,1:10));
% S7=mean(data9.bestFitness(7,end,1:10)); S8=mean(data9.bestFitness(8,end,1:10));
% S9=mean(data9.bestFitness(9,end,1:10)); S10=mean(data9.bestFitness(10,end,1:10));
% S=[S1,S2,S3,S4,S5,S6,S7,S8,S9,S10];

STD_MTEAbest=std(S)

%%
data1=MTSO(tasks, PopSize,nGen, Rmp,nRepeat,ntasks);
Mean_MTSO=mean(data1.bestFitness(1:end,end))

S1=mean(data1.bestFitness(1:5,end)); S2=mean(data1.bestFitness(6:10,end));
S3=mean(data1.bestFitness(11:15,end)); S4=mean(data1.bestFitness(16:20,end));
S5=mean(data1.bestFitness(21:25,end));S=[S1,S2,S3,S4,S5];

% S1=mean(data1.bestFitness(1:10,end)); S2=mean(data1.bestFitness(11:20,end));
% S3=mean(data1.bestFitness(21:30,end)); S4=mean(data1.bestFitness(31:40,end));
% S5=mean(data1.bestFitness(41:50,end)); S6=mean(data1.bestFitness(51:60,end));
% S7=mean(data1.bestFitness(61:70,end)); S8=mean(data1.bestFitness(71:80,end));
% S9=mean(data1.bestFitness(81:90,end)); S10=mean(data1.bestFitness(91:100,end));
% S=[S1,S2,S3,S4,S5,S6,S7,S8,S9,S10];

STD_MTSO=std(S)

%%
% A1_normalized = normalize(data1.bestFitness, 1, 'range');
% A11=mean(A1_normalized(1:end,:)); 
A1=mean(data1.bestFitness(1:end,:)); 
A2= mean(data2.bestFitness(1:end,:)); 
A3=mean(data3.bestFitness(1:end,:));
A6_20=squeeze(mean(data6.bestFitness,1)); 
%A6=mean(A6_20(end,:));
A6=mean(A6_20,2);
A7=mean(data7.bestFitness(1:end,:));
A8_20=squeeze(mean(data8.bestFitness,1));
%A8=mean(A8_20(end,:));
A8=mean(A8_20,2);
A9_20=squeeze(mean(data9.bestFitness,1));
%A9=mean(A9_20(end,:));
A9=mean(A9_20,2);
 

 figure; 
    semilogy(A2(1,:),'b-','linewidth',1);hold on;
    semilogy(A3(1,:),'k-','linewidth',1);
    semilogy(A6(:,1),'m-','linewidth',1);
    semilogy(A7(1,:),'y-','linewidth',1);
    semilogy(A8(:,1),'g-','linewidth',1);
    semilogy(A9(:,1),'c-','linewidth',1);
    semilogy(A1(1,:),'r-','linewidth',1);
    axis tight;
    legend('MFEA','MFEARR','EBSGA','GMFEA','EMTEA','MTEA','MTSO','location','northeast');
    title('Plaanar Kinematic Arm Control Problem');
    xlabel('Iteration'); ylabel('The average value of 5 tasks in 20 dimensions');
   
    set(gcf,'Position',[400 200 500 400]);

