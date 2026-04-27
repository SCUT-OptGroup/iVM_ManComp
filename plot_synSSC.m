%% ***************************************************************
% filename: plot_realSSC
%% ***************************************************************
% to test the performance of iRVM with difference mu_max for solving the problem
%%
%% min <XX',L> + lambda||XX'||_1, s.t. X'X=I_q
%%
%% **************************************************************
clear; clear all; clc

close all;

restoredefaultpath;
addpath(genpath('Data'));
addpath(genpath('SparSC'));
addpath(genpath('functions'));

%% Set problem size and lambda

n = 300; q = 6; p = 2000;

lambda = 5e-5;

Imat = eye(n);

OPTIONS_iRVM.maxiter = 5000;

OPTIONS_iRVM.printyes = 1;

OPTIONS_iRVM.tol = 5e-8;

gamma = 1e-5;  beta = 1e-2;

mu_list = [1  10   50   100   200   300   400   500   600   700   800   900   1000];

nmu = length(mu_list);

ntest = 10;

%% Zeros array to collect each iteration

NMI = zeros(nmu,ntest); obj = zeros(nmu,ntest);  time = zeros(nmu,ntest);

display = 0;

for i = 1:nmu
    i
    mu_max = mu_list(i);
       
    for k = 1:ntest
        k
        randstate = 33333+k*80000;
        randn('state',double(randstate));
        rand('state',double(randstate));
        
        %% ****************** upload the data ****************************
        
        [true_labs, X] = syn_SSC(p,n,q,display,randstate); % synthetic data
    
        simX = func_simlarity(X);   % calculate the similarity matrix
        
        rsum = sum(simX,2);
        
        dvec = 1./rsum.^(1/2);
        
        Lap = Imat - dvec.*(simX.*dvec');
        
        X0 = orth(randn(n,q)); % generate an initial point X0\in M
        
        Lf = 5e-1*svds(Lap,1);

        OPTIONS_iRVM.normLap = norm(Lap,'fro');
        
        tstart = clock;
        
        iRVM = SSC_iRVM(X0,gamma,Lf,mu_max,beta,OPTIONS_iRVM,lambda,Imat,Lap,n,q);
                
        time(i,k)= etime(clock,tstart);
        
       %% Results
        [labs,~] = kmeans(iRVM.X, q);
        NMI(i,k)   = Cal_NMI(true_labs,labs)
        obj(i,k)   = iRVM.Theta;       
    end
      
end
aNMI = mean(NMI,2); aobj = mean(obj,2); atime = mean(time,2);

save('nu_max','mu_list','aNMI','aobj','atime');

%% Plot the figures
% figure(1)          % define figure
% set(gcf,'position',[250 300 600 250])
% 
% subplot(1,2,1);
% plot(mu_list,atime,'r-','linewidth',2);
% xlabel('\mu_{max}');
% ylabel('time (s)');
% hold on
% 
% 
% subplot(1,2,2);
% plot(mu_list,aNMI,'b-','linewidth',2);
% xlabel('\mu_{max}');
% ylabel('NMI');







% subplot(2,2,3);
% plot(RXaxis,Rtime,'-*');figure(1)          % define figure
set(gcf,'position',[50 50 1000 270])
subplot(1,3,1);     % subplot(x,y,n)
plot(mu_list,aobj,'r-','linewidth',2);
xlabel('\mu_{max}');
ylabel('obj');
hold on

subplot(1,3,2);
plot(mu_list,atime,'b-','linewidth',2);
xlabel('\mu_{max}');
ylabel('time (s)');
hold on


subplot(1,3,3);
plot(mu_list,aNMI,'c-','linewidth',2);
xlabel('\mu_{max}');
ylabel('NMI');
% xlabel('\mu_{max}');
% ylabel('time');
% 
% subplot(2,2,4);
% plot(RXaxis,RIter,'-*');
% xlabel('\mu_{max}');
% ylabel('iter');

% %%
% figure(2)          % define figure
% set(gcf,'position',[250 300 600 180])
% subplot(1,2,1);     % subplot(x,y,n)
% plot(RXaxis,Robj,'-*');
% xlabel('\mu_{max}');
% ylabel('obj');
% 
% 
% subplot(1,2,2);
% plot(RXaxis,Rtime,'-*');
% xlabel('\mu_{max}');
% ylabel('time');