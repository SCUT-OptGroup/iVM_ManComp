%% ***************************************************************
% filename: demo_r21PCA
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

ntest = 10;

iRVM_obj = zeros(1,ntest); iRVM_time = zeros(1,ntest);  iRVM_rspa = zeros(1,ntest);

iRVM_voffd = zeros(1,ntest); iRVM_nsub = zeros(1,ntest); iRVM_nSN = zeros(1,ntest);

iALM_obj = zeros(1,ntest); iALM_time = zeros(1,ntest);

iALM_rspa = zeros(1,ntest); iALM_voffd = zeros(1,ntest);

admm_obj = zeros(1,ntest); admm_time = zeros(1,ntest);

admm_rspa = zeros(1,ntest); admm_voffd = zeros(1,ntest);

iRVM_iter = zeros(1,ntest); admm_iter = zeros(1,ntest);

iALM_iter = zeros(1,ntest); iALM_iter_in = zeros(1,ntest);

%% *********** Set the size and parameters of problems ************

m = 50;  n = 1000;  q = 9;

lambda = 2.0;   rho = 0.5;

%% ****************** upload the data ****************************
for jj = 1:ntest
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
    %iRVM = r21PCA_iRVMapg(X0,Lf,gamma,mu_max,beta,OPTIONS_iRVM,lambda,rho,m,n,q,E,A,B);

    iRVM_time(jj) = etime(clock,tstart);

    rs_zidx = iRVM.r21norm<=1e-4*max(1,max(iRVM.r21norm));

    rspar = sum(rs_zidx)/n;

    iRVM_obj(jj) = iRVM.Theta

    iRVM_iter(jj) = iRVM.iter;

    iRVM_rspa(jj) = rspar; iRVM_voffd(jj) = iRVM.offd_vio;

    iRVM_nsub(jj)=sum(iRVM.tnls)/nnz(iRVM.tnls); 

    iRVM_nSN(jj)=sum(iRVM.tsubit)/nnz(iRVM.tnls);
    
  %% ********************* RiALM *************************************

    OPTIONS_iALM.maxiter = 100;  OPTIONS_iALM.tol = 1e-6; 

    OPTIONS_iALM.printyes = 1;   OPTIONS_iALM.normA = normA;

    sigma_RiALM = 1.5; b = 1.5; err0 = 1.5;

    tstart = clock;

    RiALM = r21PCA_RiALM(X0,sigma_RiALM,b,err0,OPTIONS_iALM,lambda,rho,m,n,q,E,A,B);

    iALM_time(jj) = etime(clock,tstart);

    rs_zidx = RiALM.r21norm<=1e-4*max(1,max(RiALM.r21norm));

    rspar = sum(rs_zidx)/n;

    iALM_rspa(jj) = rspar;

    iALM_obj(jj) = RiALM.Theta

    iALM_iter(jj) = RiALM.iter;

    iALM_iter_in(jj) = sum(RiALM.iter_in)/nnz(RiALM.iter_in);

    iALM_voffd(jj) = RiALM.offd_vio;

  %% ********************* RADMM *************************************

    OPTIONS_RADMM.maxiter = 2e5;  OPTIONS_RADMM.tol = 1e-9;

    OPTIONS_RADMM.printyes = 1;  OPTIONS_RADMM.normA = normA;

    eta_RADMM = 5e-5; sigma_RADMM = 5e1; gamma = 1e-8;

    tstart = clock;

    RADMM = r21PCA_RADMM(X0,eta_RADMM,sigma_RADMM,gamma,OPTIONS_RADMM,lambda,rho,m,n,q,E,A,B);

    admm_time(jj) = etime(clock,tstart);

    rs_zidx = RADMM.r21norm<=1e-4*max(1,max(RADMM.r21norm));

    rspar = sum(rs_zidx)/n;

    admm_rspa(jj) = rspar;

    admm_obj(jj) = RADMM.Theta; 
    
    admm_iter(jj) = RADMM.iter;
    
    admm_voffd(jj) = RADMM.offd_vio;
end

aiRVM_obj = mean(iRVM_obj); aiRVM_time = mean(iRVM_time); aiRVM_nsub = mean(iRVM_nsub); aiRVM_nSN = mean(iRVM_nSN);

aiRVM_rspa = mean(iRVM_rspa); aiRVM_offd = mean(iRVM_voffd);

iRVM_aiter = mean(iRVM_iter); 

fprintf('\n Method: %s, obj#: %3.4e, row spar#: %f, off-diag val#: %3.2e,CPU time: %.2f, nls#: %.2f, nSSN#: %.2f, Iter: %.1f','iRVM',aiRVM_obj,aiRVM_rspa,aiRVM_offd,aiRVM_time,aiRVM_nsub,aiRVM_nSN,iRVM_aiter);

aadmm_obj = mean(admm_obj); aadmm_time = mean(admm_time);

aadmm_rspa = mean(admm_rspa); aadmm_offd = mean(admm_voffd);

admm_aiter = mean(admm_iter);

fprintf('\n Method: %s, obj#: %3.4e, row spar#: %f, off-diag val#: %3.2e,CPU time: %.2f, Iter: %.1f','RADMM',aadmm_obj,aadmm_rspa,aadmm_offd,aadmm_time,admm_aiter);

aiALM_obj = mean(iALM_obj); aiALM_time = mean(iALM_time);

aiALM_rspa = mean(iALM_rspa); aiALM_offd = mean(iALM_voffd);

iALM_aiter = mean(iALM_iter); iALM_aiter_in = mean(iALM_iter_in);

fprintf('\n Method: %s, obj#: %3.4e, row spar#: %f, off-diag val#: %3.2e,CPU time: %.2f, Iter: %.1f, Iter_inner: %.1f','RiALM',aiALM_obj,aiALM_rspa,aiALM_offd,aiALM_time,iALM_aiter,iALM_aiter_in);
