%% ***************************************************************
% filename: plot_rsPCA_SNcg_APG
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

ntest = 1;

iRVM_obj = cell(ntest,1); iRVM_SNcg = cell(ntest,1); iRVM_Iter = cell(ntest,1); iRVM_Time = cell(ntest,1); miRVM_iter = 5000;

iRVMapg_obj = cell(ntest,1); iRVM_apg = cell(ntest,1); iRVMapg_Iter = cell(ntest,1); iRVMapg_Time = cell(ntest,1); miRVMapg_iter = 5000;

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

    OPTIONS_iRVM.maxiter = 1000;  OPTIONS_iRVM.tol = -1e-8;
    
    OPTIONS_iRVM.printyes = 1;   OPTIONS_iRVM.normA = normA;
    
    iRVM_gamma = 1e-5;  mu_max = 5e2;  beta = 1e-2;
    
    iRVM = plot_PCA_iRVM1(X0,Lf,iRVM_gamma,mu_max,beta,OPTIONS_iRVM,lambda,rho,m,n,q,E,A,B);
    
    iRVM_obj{jj} = iRVM.Tobj; iRVM_SNcg{jj} = iRVM.SNcg; 
    
    iRVM_Iter{jj} = iRVM.Iter;  iRVM_Time{jj} = iRVM.Ttime;

    miRVM_iter = min(miRVM_iter,nnz(iRVM.Iter));

    % ********************* iRVMapg *************************************

    OPTIONS_iRVM.maxiter = 1000;  OPTIONS_iRVM.tol = -1e-8;
    
    OPTIONS_iRVM.printyes = 1;   OPTIONS_iRVM.normA = normA;
    
    iRVM_gamma = 1e-5;  mu_max = 5e2;  beta = 1e-2;
    
    iRVMapg = plot_r21PCA_iRVMapg(X0,Lf,iRVM_gamma,mu_max,beta,OPTIONS_iRVM,lambda,rho,m,n,q,E,A,B);
    
    iRVMapg_obj{jj} = iRVMapg.Tobj; iRVM_apg{jj} = iRVMapg.APG; 
    
    iRVMapg_Iter{jj} = iRVMapg.Iter;  iRVMapg_Time{jj} = iRVMapg.Ttime;

    miRVMapg_iter = min(miRVMapg_iter,nnz(iRVMapg.Iter));

   
end

iRVMobj = zeros(ntest,miRVM_iter); iRVMSNcg = zeros(ntest,miRVM_iter);  iRVMTime = zeros(ntest,miRVM_iter); aiRVM_Iter = iRVM_Iter{1}(1:miRVM_iter);
for i = 1 : ntest
    iRVMobj(i,:)=iRVM_obj{i}(1:miRVM_iter); iRVMSNcg(i,:)=iRVM_SNcg{i}(1:miRVM_iter); iRVMTime(i,:)=iRVM_Time{i}(1:miRVM_iter); 
end

iRVMApg_obj = zeros(ntest,miRVMapg_iter); iRVMApg = zeros(ntest,miRVMapg_iter);  iRVMapgTime = zeros(ntest,miRVMapg_iter); aiRVMapg_Iter = iRVMapg_Iter{1}(1:miRVMapg_iter);
for i = 1 : ntest
    iRVMApg_obj(i,:)=iRVMapg_obj{i}(1:miRVMapg_iter); iRVMApg(i,:)=iRVM_apg{i}(1:miRVMapg_iter); iRVMapgTime(i,:)=iRVMapg_Time{i}(1:miRVMapg_iter); 
end

aiRVMobj = mean(iRVMobj,1); aiRVMSNcg = mean(iRVMSNcg,1); aiRVMTime = mean(iRVMTime,1);

aiRVMApg_obj = mean(iRVMApg_obj,1); aiRVMapg = mean(iRVMApg,1); aiRVMapgTime = mean(iRVMapgTime,1);


%% Plot the figures
% 
% figure(1)          % define figure
% set(gcf,'position',[50 50 1000 270])
% 
% subplot(1,3,1);
% plot(aiRVM_Iter,aiRVMSNcg,'r-','linewidth',2);
% hold on
% plot(aiRVMapg_Iter,aiRVMapg,'b-','linewidth',2);
% xlabel('Iter');
% ylabel('subiter');
% legend('SNcg','APG');
% 
% subplot(1,3,2);
% plot(aiRVM_Iter,aiRVMTime,'r-','linewidth',2);
% hold on
% plot(aiRVMapg_Iter,aiRVMapgTime,'b-','linewidth',2);
% xlabel('Iter');
% ylabel('time (s)');
% legend('SNcg','APG');
% 
% 
% subplot(1,3,3);
% plot(aiRVM_Iter,aiRVMobj,'r-','linewidth',2);
% hold on
% plot(aiRVMapg_Iter,aiRVMApg_obj,'b-','linewidth',2);
% xlabel('Iter');
% ylabel('obj');
% legend('SNcg','APG');







figure(1)          % define figure
set(gcf,'position',[50 50 1000 370])
subplot(1,2,1);     % subplot(x,y,n)
plot(aiRVM_Iter,aiRVMSNcg,'r-','linewidth',2);
hold on
plot(aiRVMapg_Iter,aiRVMapg,'b-','linewidth',2);
xlabel('iter');
ylabel('subiter');
legend('RiVMPL-SNcg','RiVMPL-DFG');

subplot(1,2,2);
plot(aiRVM_Iter,aiRVMTime,'r-','linewidth',2);
hold on
plot(aiRVMapg_Iter,aiRVMapgTime,'b-','linewidth',2);
xlabel('iter');
ylabel('time (s)');
legend('RiVMPL-SNcg','RiVMPL-DFG');




