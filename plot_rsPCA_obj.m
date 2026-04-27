%% ***************************************************************
% filename: plot_rsPCA_obj
%% ***************************************************************
% to test the performance of the iRVM for solving the problem
%%
%% min -tr(X^THX) + lambda||X||_1 + rho||offdiag(X^THX)||_1, s.t. U^TU=I_q
%%
%% **************************************************************

restoredefaultpath;
addpath(genpath('functions'));
%addpath(genpath('r21PCA_draw'));
addpath(genpath('r21PCA'));


%% Set problem size and lambda

m = 50; n = 1000; q = 5;  

lambda = 2.05;   rho = 0.5;

ntest = 10;

iRVM_obj = cell(ntest,1); iRVM_Iter = cell(ntest,1); iRVM_Time = cell(ntest,1); miRVM_iter = 5000;

ALM_obj = cell(ntest,1); ALM_Iter = cell(ntest,1); ALM_Time = cell(ntest,1); mALM_iter = 100;

ADMM_obj = cell(ntest,1); ADMM_Iter = cell(ntest,1); ADMM_Time = cell(ntest,1); mADMM_iter = 2e5;

for jj = 1 : ntest
    
    jj

    randstate = 88888+jj*10000;
    randn('state',double(randstate));
    rand('state',double(randstate));
    
    B = randn(m,n);
    
    B = B - repmat(mean(B,1),m,1);
    
    B = normc(B);
    
    A = B'*B;  E = ones(q)-eye(q);
    
    normA = norm(A,'fro');
      
    Lf = 5e-1*svds(B,1)^2;  %
    
    X0 = orth(randn(n,q)); % generate an initial point X0\in M

    % ********************* iRVM *************************************

    OPTIONS_iRVM.maxiter = 5000;  OPTIONS_iRVM.tol = 1e-8;
    
    OPTIONS_iRVM.printyes = 1;   OPTIONS_iRVM.normA = normA;
    
    iRVM_gamma = 1e-5;  mu_max = 5e2;  beta = 1e-2;
    
    iRVM = plot_PCA_iRVM(X0,Lf,iRVM_gamma,mu_max,beta,OPTIONS_iRVM,lambda,rho,m,n,q,E,A,B);
    
    iRVM_obj{jj} = iRVM.Tobj; iRVM_Iter{jj} = iRVM.Iter;  iRVM_Time{jj} = iRVM.Ttime;

    miRVM_iter = min(miRVM_iter,nnz(iRVM.Iter));

    %% ********************* RiAl_RGD *************************************

    OPTIONS_iALM.maxiter = 100;  OPTIONS_iALM.tol =1e-7;
    
    OPTIONS_iALM.printyes =1;   OPTIONS_iALM.normA = norm(A,'fro');
    
    sigma_RiAl_RGD = 1.5;  b = 1.5; err = 1.5;
    
    iALM = plot_PCA_RiALM(X0,sigma_RiAl_RGD,b,err,OPTIONS_iALM,lambda,rho,m,n,q,E,A,B);
    
    ALM_obj{jj} = iALM.Tobj; ALM_Iter{jj} = iALM.Iter; ALM_Time{jj} = iALM.Ttime;

    mALM_iter = min(mALM_iter,nnz(iALM.Iter));

    %% ********************* RADMM *************************************

    OPTIONS_RADMM.maxiter = 2e5;  OPTIONS_RADMM.tol = 1e-10;
    
    OPTIONS_RADMM.printyes = 1;  OPTIONS_RADMM.normA = norm(A,'fro');
    
    eta_RADMM = 5e-5;  sigma_RADMM = 5e1; gamma_RADMM = 1e-8;
    
    RADMM = plot_PCA_RADMM(X0,eta_RADMM,sigma_RADMM,gamma_RADMM,OPTIONS_RADMM,lambda,rho,m,n,q,E,A,B);
    
    ADMM_obj{jj} = RADMM.Tobj; ADMM_Iter{jj} = RADMM.Iter; ADMM_Time{jj} = RADMM.Ttime;

    mADMM_iter = min(mADMM_iter,nnz(RADMM.Iter));

end

miRVM_iter = miRVM_iter+1; mALM_iter = mALM_iter + 1; mADMM_iter = mADMM_iter + 1;

iRVMobj = zeros(ntest,miRVM_iter);  iRVMTime = zeros(ntest,miRVM_iter); aiRVM_Iter = iRVM_Iter{1}(1:miRVM_iter);
for i = 1 : ntest
    iRVMobj(i,:)=iRVM_obj{i}(1:miRVM_iter); iRVMTime(i,:)=iRVM_Time{i}(1:miRVM_iter); 
end

ALMobj = zeros(ntest,mALM_iter);  ALMTime = zeros(ntest,mALM_iter); aALM_Iter = ALM_Iter{1}(1:mALM_iter);
for i = 1 : ntest
    ALMobj(i,:)=ALM_obj{i}(1:mALM_iter); ALMTime(i,:)=ALM_Time{i}(1:mALM_iter); 
end

ADMMobj = zeros(ntest,mADMM_iter);  ADMMTime = zeros(ntest,mADMM_iter); aADMM_Iter = ADMM_Iter{1}(1:mADMM_iter);
for i = 1 : ntest
    ADMMobj(i,:)=ADMM_obj{i}(1:mADMM_iter); ADMMTime(i,:)=ADMM_Time{i}(1:mADMM_iter); 
end

aiRVMobj = mean(iRVMobj,1); aiRVMTime = mean(iRVMTime,1);

aALMobj = mean(ALMobj,1); aALMTime = mean(ALMTime,1);

aADMMobj = mean(ADMMobj,1); aADMMTime = mean(ADMMTime,1);


%% Plot the figures
figure(1)          % define figure
set(gcf,'position',[50 50 1000 450])
subplot(2,3,1);     % subplot(x,y,n)
plot(aiRVM_Iter,aiRVMobj,'r-','linewidth',2);
xlabel('iter');
ylabel('obj');
xlim([0 miRVM_iter]);
legend('RiVMPL');
hold on


subplot(2,3,2);
plot(aALM_Iter,aALMobj,'b-','linewidth',2);
xlabel('iter');
ylabel('obj');
xlim([0 mALM_iter]);
legend('RiALM');
hold on

subplot(2,3,3);
plot(aADMM_Iter,aADMMobj,'c-','linewidth',2);
xlabel('iter');
ylabel('obj');
xlim([0 mADMM_iter]);
legend('RADMM');


subplot(2,3,4);     % subplot(x,y,n)
plot(aiRVM_Iter,aiRVMTime,'r-','linewidth',2);
xlabel('iter');
ylabel('time (s)');
legend('RiVMPL');
xlim([0 miRVM_iter]);
hold on

subplot(2,3,5);
plot(aALM_Iter,aALMTime,'b-','linewidth',2);
xlabel('iter');
ylabel('time (s)');
xlim([0 mALM_iter]);
legend('RiALM');
hold on

subplot(2,3,6);
plot(aADMM_Iter,aADMMTime,'c-','linewidth',2);
xlabel('iter');
ylabel('time (s)');
xlim([0 mADMM_iter]);
legend('RADMM');



