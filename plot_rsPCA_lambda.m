%% ***************************************************************
% filename: plot_rsPCA_lambda
%% **************************************************************
% to test the performance of the iRVM for solving the problem
%%
%% min -tr(X^THX) + lambda||X||_1 + rho||offdiag(X^THX)||_1, s.t. U^TU=I_q
%%
%% **************************************************************

%clear all;

restoredefaultpath;
addpath(genpath('functions'));
addpath(genpath('r21PCA'));

%% *********************** Initialization *************************

lamgroup = [0.6  0.8  1.0  1.2  1.4  1.6  1.8  2.0  2.2  2.4  2.6 2.8 3.0 3.2];

ntest = 10; 

nlam = length(lamgroup);

iRVM_obj = zeros(nlam,ntest); iRVM_time = zeros(nlam,ntest); iRVM_rspa = zeros(nlam,ntest); 

iALM_obj = zeros(nlam,ntest); iALM_time = zeros(nlam,ntest); iALM_rspa = zeros(nlam,ntest); 

admm_obj = zeros(nlam,ntest); admm_time = zeros(nlam,ntest); admm_rspa = zeros(nlam,ntest);

% iRVM_obj = [iRVM_obj; zeros(3,ntest)]; iRVM_time = [iRVM_time; zeros(3,ntest)]; iRVM_rspa = [iRVM_rspa; zeros(3,ntest)];
% 
% iALM_obj = [iALM_obj; zeros(3,ntest)]; iALM_time = [iALM_time; zeros(3,ntest)]; iALM_rspa = [iALM_rspa; zeros(3,ntest)];
% 
% admm_obj = [admm_obj; zeros(3,ntest)]; admm_time = [admm_time; zeros(3,ntest)]; admm_rspa = [admm_rspa; zeros(3,ntest)];

%% *********** Set the size and parameters of problems ************

m = 50;  n = 1000;  q = 5;

rho = 0.5;

%% ****************** upload the data ****************************
for i = 1 : nlam

    lambda = lamgroup(i);

    for jj = 1:ntest
        i
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
        
        %% ********************* iRVM *************************************

        OPTIONS_iRVM.maxiter = 5000;  OPTIONS_iRVM.tol = 5e-8;  

        OPTIONS_iRVM.printyes = 1;   OPTIONS_iRVM.normA = normA;

        gamma = 1e-5;  mu_max = 5e2;  beta = 1e-2;

        tstart = clock;

        iRVM = r21PCA_iRVM(X0,Lf,gamma,mu_max,beta,OPTIONS_iRVM,lambda,rho,m,n,q,E,A,B);

        iRVM_time(i,jj) = etime(clock,tstart);

        rs_zidx = iRVM.r21norm<=1e-4*max(1,max(iRVM.r21norm));

        rspar = sum(rs_zidx)/n;

        iRVM_obj(i,jj) = iRVM.Theta; iRVM_rspa(i,jj) = rspar; 

    %% ********************* RiALM *************************************

        OPTIONS_iALM.maxiter = 100;  OPTIONS_iALM.tol = 1e-6; 

        OPTIONS_iALM.printyes =1;   OPTIONS_iALM.normA = normA;

        sigma_RiALM = 1.5; b = 1.5; err0 = 1.5;

        tstart = clock;

        RiALM = r21PCA_RiALM(X0,sigma_RiALM,b,err0,OPTIONS_iALM,lambda,rho,m,n,q,E,A,B);

        iALM_time(i,jj) = etime(clock,tstart);

        rs_zidx = RiALM.r21norm<=1e-4*max(1,max(RiALM.r21norm));

        rspar = sum(rs_zidx)/n;

        iALM_rspa(i,jj) = rspar; iALM_obj(i,jj) = RiALM.Theta; 

      %% ********************* RADMM *************************************

        OPTIONS_RADMM.maxiter = 2e5;  OPTIONS_RADMM.tol = 1e-9;

        OPTIONS_RADMM.printyes = 1;  OPTIONS_RADMM.normA = normA;

        eta_RADMM = 5e-5; sigma_RADMM = 5e1; gamma = 1e-8;

        tstart = clock;

        RADMM = r21PCA_RADMM(X0,eta_RADMM,sigma_RADMM,gamma,OPTIONS_RADMM,lambda,rho,m,n,q,E,A,B);

        admm_time(i,jj) = etime(clock,tstart);

        rs_zidx = RADMM.r21norm<=1e-4*max(1,max(RADMM.r21norm));

        rspar = sum(rs_zidx)/n;

        admm_rspa(i,jj) = rspar; admm_obj(i,jj) = RADMM.Theta; 
end

end

aiRVM_obj = mean(iRVM_obj,2); aiRVM_time = mean(iRVM_time,2); aiRVM_rspa = mean(iRVM_rspa,2); 

aadmm_obj = mean(admm_obj,2); aadmm_time = mean(admm_time,2); aadmm_rspa = mean(admm_rspa,2); 

aiALM_obj = mean(iALM_obj,2); aiALM_time = mean(iALM_time,2); aiALM_rspa = mean(iALM_rspa,2); 

%% Plot the figures
figure(1)          % define figure
set(gcf,'position',[50 50 1000 370])
subplot(1,2,1);     % subplot(x,y,n)
plot(lamgroup,aiRVM_rspa,'r-','linewidth',2);
hold on
plot(lamgroup,aiALM_rspa,'b-','linewidth',2);
hold on
plot(lamgroup,aadmm_rspa,'c-','linewidth',2);
xlabel('\lambda');
ylabel('sparsity');
legend('RiVMPL','RiALM','RADMM');

subplot(1,2,2);
plot(lamgroup,aiRVM_time,'r-','linewidth',2);
hold on
plot(lamgroup,aiALM_time,'b-','linewidth',2);
hold on
plot(lamgroup,aadmm_time,'c-','linewidth',2);
xlabel('\lambda');
ylabel('time (s)');
legend('RiVMPL','RiALM','RADMM');

