%% ***************************************************************
% filename: plot_rsPCA_rho
%% **************************************************************
% to test the performance of the iRVM for solving the problem
%%
%% min -tr(X^THX) + lambda||X||_1 + rho||offdiag(X^THX)||_1, s.t. U^TU=I_q
%%
%% **************************************************************

clear; clear all; clc;

restoredefaultpath;
addpath(genpath('functions'));
addpath(genpath('r21PCA'));

%% *********************** Initialization *************************

rhogroup = [0.05 0.1  0.2  0.3  0.4  0.5  0.6  0.7  0.8  0.9  1.0];

ntest = 10; 

nrho = length(rhogroup);

iRVM_obj = zeros(nrho,ntest); iRVM_time = zeros(nrho,ntest); iRVM_rspa = zeros(nrho,ntest); iRVM_voffd = zeros(nrho,ntest); 

iALM_obj = zeros(nrho,ntest); iALM_time = zeros(nrho,ntest); iALM_rspa = zeros(nrho,ntest); iALM_voffd = zeros(nrho,ntest); 

admm_obj = zeros(nrho,ntest); admm_time = zeros(nrho,ntest); admm_rspa = zeros(nrho,ntest); admm_voffd = zeros(nrho,ntest); 
% 
% iRVM_obj = [zeros(1,ntest);iRVM_obj]; iRVM_time = [zeros(1,ntest);iRVM_time]; iRVM_rspa = [zeros(1,ntest);iRVM_rspa]; iRVM_voffd = [zeros(1,ntest);iRVM_voffd];
% 
% iALM_obj = [zeros(1,ntest);iALM_obj]; iALM_time = [zeros(1,ntest);iALM_time]; iALM_rspa = [zeros(1,ntest);iALM_rspa]; iALM_voffd = [zeros(1,ntest);iALM_voffd];
% 
% admm_obj = [zeros(1,ntest);admm_obj]; admm_time = [zeros(1,ntest);admm_time]; admm_rspa = [zeros(1,ntest);admm_rspa]; admm_voffd = [zeros(1,ntest);admm_voffd];

%% *********** Set the size and parameters of problems ************

m = 50;  n = 1000;  q = 5;

lambda = 2.05;

%% ****************** upload the data ****************************
for i = 1 : nrho

    rho = rhogroup(i);

    for jj = 1:ntest
        i
        jj
        randstate = 88888+(jj-1)*10000;
        randn('state',double(randstate));
        rand('state',double(randstate));
        
        B = randn(m,n);
        
        B = B - repmat(mean(B,1),m,1);
        
        B = normc(B);
        
        A = B'*B;  E = ones(q)-eye(q);
        
        normA = norm(A,'fro');
          
        Lf = 5e-1*svds(B,1)^2;  %
        
        X0 = orth(randn(n,q)); % generate an initial point X0\in M
        
        %% ********************* iRVM *************************************

        OPTIONS_iRVM.maxiter = 5000;  OPTIONS_iRVM.tol = 5e-8;  

        OPTIONS_iRVM.printyes = 1;   OPTIONS_iRVM.normA = normA;

        gamma = 1e-5;  mu_max = 5e2;  beta = 1e-2;

        tstart = clock;

        iRVM = r21PCA_iRVM(X0,Lf,gamma,mu_max,beta,OPTIONS_iRVM,lambda,rho,m,n,q,E,A,B);

        iRVM_time(i,jj) = etime(clock,tstart);

        rs_zidx = iRVM.r21norm<=1e-4*max(1,max(iRVM.r21norm));

        rspar = sum(rs_zidx)/n;

        iRVM_obj(i,jj) = iRVM.Theta; iRVM_rspa(i,jj) = rspar; iRVM_voffd(i,jj) = iRVM.offd_vio;

    %% ********************* RiALM *************************************

        OPTIONS_iALM.maxiter = 100;  OPTIONS_iALM.tol = 1e-6; 

        OPTIONS_iALM.printyes = 1;   OPTIONS_iALM.normA = normA;

        sigma_RiALM = 1.5; b = 1.5; err0 = 1.5;

        tstart = clock;

        RiALM = r21PCA_RiALM(X0,sigma_RiALM,b,err0,OPTIONS_iALM,lambda,rho,m,n,q,E,A,B);

        iALM_time(i,jj) = etime(clock,tstart);

        rs_zidx = RiALM.r21norm<=1e-4*max(1,max(RiALM.r21norm));

        rspar = sum(rs_zidx)/n;

        iALM_rspa(i,jj) = rspar; iALM_obj(i,jj) = RiALM.Theta; iALM_voffd(i,jj) = RiALM.offd_vio;

      %% ********************* RADMM *************************************

        OPTIONS_RADMM.maxiter = 2e5;  OPTIONS_RADMM.tol = 1e-9;

        OPTIONS_RADMM.printyes = 1;  OPTIONS_RADMM.normA = normA;

        eta_RADMM = 5e-5; sigma_RADMM = 5e1; gamma = 1e-8;

        tstart = clock;

        RADMM = r21PCA_RADMM(X0,eta_RADMM,sigma_RADMM,gamma,OPTIONS_RADMM,lambda,rho,m,n,q,E,A,B);

        admm_time(i,jj) = etime(clock,tstart);

        rs_zidx = RADMM.r21norm<=1e-4*max(1,max(RADMM.r21norm));

        rspar = sum(rs_zidx)/n;

        admm_rspa(i,jj) = rspar; admm_obj(i,jj) = RADMM.Theta; admm_voffd(i,jj) = RADMM.offd_vio;
end

end

aiRVM_obj = mean(iRVM_obj,2); aiRVM_time = mean(iRVM_time,2); aiRVM_rspa = mean(iRVM_rspa,2); aiRVM_voffd = mean(iRVM_voffd,2);

aadmm_obj = mean(admm_obj,2); aadmm_time = mean(admm_time,2); aadmm_rspa = mean(admm_rspa,2); aadmm_voffd = mean(admm_voffd,2);

aiALM_obj = mean(iALM_obj,2); aiALM_time = mean(iALM_time,2); aiALM_rspa = mean(iALM_rspa,2); aiALM_voffd = mean(iALM_voffd,2);

%% Plot the figures
figure(1)          % define figure
set(gcf,'position',[50 50 1000 370])
subplot(1,2,1);     % subplot(x,y,n)
plot(rhogroup,aiRVM_voffd,'r-','linewidth',2);
hold on
plot(rhogroup,aiALM_voffd,'b-','linewidth',2);
hold on
plot(rhogroup,aadmm_voffd,'c-','linewidth',2);
xlabel('\rho');
ylabel('infeasibility');
legend('RiVMPL','RiALM','RADMM');

subplot(1,2,2);
plot(rhogroup,aiRVM_time,'r-','linewidth',2);
hold on
plot(rhogroup,aiALM_time,'b-','linewidth',2);
hold on
plot(rhogroup,aadmm_time,'c-','linewidth',2);
xlabel('\rho');
ylabel('time (s)');
legend('RiVMPL','RiALM','RADMM');