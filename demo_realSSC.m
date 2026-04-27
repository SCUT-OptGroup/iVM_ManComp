%% ***************************************************************
% Filename: demo_reSSC
%% ***************************************************************
% to test the performance of iRVM, RADMM, RiAL_RDG for solving the problem
%%
%% min <XX',L> + lambda||XX'||_1, s.t. X'X=I_q
%%
%% **************************************************************
clear; clear all; clc;

restoredefaultpath;
addpath(genpath('Data'));
addpath(genpath('functions'));
addpath(genpath('SparSC'));

ntest = 10;

iRVM_obj = zeros(1,ntest); iRVM_time = zeros(1,ntest);  

iRVM_score = zeros(1,ntest);  iRVM_spa = zeros(1,ntest); 

iRVM_nsub = zeros(1,ntest); iRVM_nSN = zeros(1,ntest);

iALM_obj = zeros(1,ntest); iALM_time = zeros(1,ntest);  

iALM_score = zeros(1,ntest); iALM_spa = zeros(1,ntest); 

admm_obj = zeros(1,ntest); admm_time = zeros(1,ntest);  

admm_score = zeros(1,ntest); admm_spa = zeros(1,ntest); 

iRVM_iter = zeros(1,ntest); admm_iter = zeros(1,ntest);

iALM_iter = zeros(1,ntest); iALM_iter_in = zeros(1,ntest);

%% ***************** upload test problems **************************

num = 6 %[1 2 3 5 7 8] 1:Data_Buettner; 2:Data_Deng;  3:Data_Schlitzer; 4:Data_Macosko;
% 5:Data_Pollen;  6:Data_Tasic; 7:Data_Ting; 8:Data_Treutlin

[si_X,true_labs] = info_data(num);

[n, q, Lap] = Lapla_fun(si_X,true_labs);
n
q
norm_Lap = norm(Lap,'fro');

%% ****************** set the stop condition *********************

for j=1:ntest
    j
    randstate = 9999+j*2000;
    randn('state',double(randstate));
    rand('state',double(randstate));
    
    state = rng;
    
    lambda = 1e-5;
    
    Imat = eye(n);
    
    Lf = 5e-1*svds(Lap,1);
    
    X0 = orth(randn(n,q)); % generate an initial point X0\in M
    
    %% ********************* iRVM *************************************
    
    OPTIONS_iRVM.maxiter = 5000;  OPTIONS_iRVM.tol = 1e-8;

    OPTIONS_iRVM.printyes = 1;    OPTIONS_iRVM.normLap = norm_Lap;

    gamma = 1e-5;  mu_max = 5e+2; beta = 1e-2;  

    rng(state);

    tstart = clock;

    iRVM = SSC_iRVM(X0,gamma,Lf,mu_max,beta,OPTIONS_iRVM,lambda,Imat,Lap,n,q);

    iRVM_time(j) = etime(clock,tstart);

    iRVM_obj(j) = iRVM.Theta;  

    iRVM_iter(j) = iRVM.iter;

    temp_XXtv = iRVM.XXt(:);

    rs_zidx = abs(temp_XXtv)<=1e-4*max(1,max(abs(temp_XXtv)));

    iRVM_spa(j) = sum(rs_zidx)/n^2;

    [labs,~] = kmeans(iRVM.X,q);

    iRVM_score(j) = Cal_NMI(true_labs, labs);

    iRVM_nsub(j)=sum(iRVM.tnls)/nnz(iRVM.tnls);

    iRVM_nSN(j)=sum(iRVM.tsubit)/nnz(iRVM.tnls);

%% ********************* RiALM *************************************
    
    OPTIONS_iALM.maxiter = 100;  OPTIONS_iALM.tol = 1e-8; OPTIONS_iALM.printyes = 1;
    
    sigma_iALM = 1.5;  b = 1.5; err = 1e-3;
    
    rng(state);
    
    tstart = clock;

    RiALM = SSC_RiALM(X0,Lap,sigma_iALM,b,err,OPTIONS_iALM,lambda,n,q);
    
    iALM_time(j) = etime(clock,tstart);
    
    iALM_obj(j) = RiALM.Theta;

    iALM_iter(j) = RiALM.iter;

    iALM_iter_in(j) = sum(RiALM.iter_in)/nnz(RiALM.iter_in);
    
    temp_XXtv = RiALM.XXt(:);
    
    rs_zidx = abs(temp_XXtv)<=1e-4*max(1,max(abs(temp_XXtv)));
    
    iALM_spa(j) = sum(rs_zidx)/n^2;
    
    [labs,~] = kmeans(RiALM.X,q);
    
    iALM_score(j) = Cal_NMI(true_labs, labs);
    %
    %% ********************* RADMM *************************************

    OPTIONS_admm.maxiter = 1e5; OPTIONS_admm.tol = 1e-10; OPTIONS_admm.printyes = 1;

    eta_RADMM = 1e-1;  sigma_RADMM = 50; gamma = 1e-8;

    rng(state);

    tstart = clock;

    RADMM = SSC_RADMM(X0,Lap,eta_RADMM,sigma_RADMM,gamma,OPTIONS_admm,lambda,n,q);

    temp_XXtv = RADMM.XXt(:);

    rs_zidx = abs(temp_XXtv)<=1e-4*max(1,max(abs(temp_XXtv)));

    admm_spa(j) = sum(rs_zidx)/n^2;

    admm_time(j) = etime(clock,tstart);

    admm_obj(j) = RADMM.Theta;

    admm_iter(j) = RADMM.iter;

    [labs,~] = kmeans(RADMM.X,q);

    admm_score(j) = Cal_NMI(true_labs, labs);
    
end

iRVM_aNMI = mean(iRVM_score); iRVM_aobj = mean(iRVM_obj);  iRVM_aspa = mean(iRVM_spa);

iRVM_atime = mean(iRVM_time); iRVM_ansub = mean(iRVM_nsub); iRVM_anSN = mean(iRVM_nSN);

iRVM_aiter = mean(iRVM_iter); 

fprintf('\n Method: %s, obj#: %3.4e, spar#: %5.4e, nsubj#: %3.2e, nSN#: %3.2e, CPU time: %f, NMI: %f, Iter: %3.1f', 'iRVM', iRVM_aobj, iRVM_aspa,iRVM_ansub,iRVM_anSN,iRVM_atime,iRVM_aNMI,iRVM_aiter);

iALM_aNMI = mean(iALM_score); iALM_aobj = mean(iALM_obj); 

iALM_aspa = mean(iALM_spa); iALM_atime = mean(iALM_time);

iALM_aiter = mean(iALM_iter); iALM_aiter_in = mean(iALM_iter_in);

fprintf('\n Method: %s, obj#: %3.4e, spar#: %5.4e, CPU time: %f, NMI: %f, Iter: %3.1f, Iter_inner: %3.1f', 'iALM', iALM_aobj,iALM_aspa,iALM_atime,iALM_aNMI,iALM_aiter,iALM_aiter_in);


admm_aNMI = mean(admm_score); admm_aobj = mean(admm_obj); 

admm_aspa = mean(admm_spa); admm_atime = mean(admm_time);

admm_aiter = mean(admm_iter);

fprintf('\n Method: %s, obj#: %3.4e, spar#: %5.4e,CPU time: %f, NMI: %f, Iter: %3.1f \n', 'ADMM', admm_aobj, admm_aspa,admm_atime,admm_aNMI,admm_aiter);

