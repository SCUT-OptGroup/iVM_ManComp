%% ***************************************************************
% filename: demo_SDec_iter
%% **************************************************************
% to test the performance of the iRVM for solving the problem
%%
%% min ||A-XX^+A||_1, s.t. X^TJ_2nX = J_2q, where X^+=J_2q^TX^TJ_2n
%%
%% **************************************************************

restoredefaultpath;

addpath(genpath('functions'));

%addpath(genpath('SDec_Fro'));

addpath(genpath('Sdecom'));


for data_type = [1, 2]
    m = 50; n = 500;  r = 20;
    for q = [5, 10, 15]

        %% *********************** Initialization *************************

        ntest = 10;
        
        iRVM_obj = zeros(1,ntest); iRVM_time = zeros(1,ntest); iRVM_nls = zeros(1,ntest); iRVM_nssn = zeros(1,ntest);
        
        iALM_obj = zeros(1,ntest); iALM_time = zeros(1,ntest); 
        
        admm_obj = zeros(1,ntest); admm_time = zeros(1,ntest);

        iRVM_iter = zeros(1,ntest); admm_iter = zeros(1,ntest);

        iALM_iter = zeros(1,ntest); iALM_iter_in = zeros(1,ntest);
        
        %% ********** Set the size of problems **************************
        
        Imat2q = eye(2*q);  
        
        J2q= [zeros(q,q)  eye(q); -eye(q)  zeros(q)];
        
        for jj = 1:ntest
            %jj
            %% ****************** upload the data ****************************
            randstate = 88888+jj*10000;
            randn('state',double(randstate));
            rand('state',double(randstate));
        
            state = rng;
            
            %% ************ to generate data and an initial X0\in M ************
            
            switch data_type
                
                case 1
                    
                    A = randn(2*n,2*m); normA = norm(A,'fro');  A = A/normA;
                    
                    A1 = A(1:n,:);   A2 = A(n+1:end,:);                     
                    
                case 2
        
                    M0 = symplecticStiefelfactory(n,r,0);
        
                    T = M0.rand();
        
                    C = randn(2*r,2*m); C = C/norm(C,'fro');
        
                    A = T*C; 
        
                    A1 = A(1:n,:); A2 = A(n+1:end,:); normA = norm(A,'fro');
                    
                    % W = randn(2*r,2*r); W = W'*W+0.8*eye(2*r); 
                    % 
                    % E = expm([W(r+1:end,:); -W(1:r,:)]);
                    % 
                    % C = randn(2*r,2*m); C = C/norm(C,'fro');
                    % 
                    % A = [E(1:r,:);zeros(n-r,2*r);E(r+1:end,:);zeros(n-r,2*r)]*C;
                    % 
                    % A1 = A(1:n,:); A2 = A(n+1:end,:); normA = norm(A,'fro');  
            end
        
            M = symplecticStiefelfactory(n,q,0);
        
            X0 = M.rand();
            
            %W = randn(2*q,2*q); W = W'*W+0.1*eye(2*q);  E = expm([W(q+1:end,:); -W(1:q,:)]);
             
            %X0 = [E(1:q,:);zeros(n-q,2*q); E(q+1:end,:);zeros(n-q,2*q)];
            
            % J2nX0 = [X0(n+1:end,:); -X0(1:n,:)];
            % 
            % norm(X0'*J2nX0-J2q,'fro')
            
            %% ********************* iRVM *************************************
        
            OPTIONS_iRVM.maxiter = 5000;  OPTIONS_iRVM.tol = 1e-7;
        
            OPTIONS_iRVM.printyes = 0;    OPTIONS_iRVM.normA = normA;
        
            Lf = 1e-5;  mu_max = 5e+2;  beta = 1e-2;
        
            gamma = 1e-5;
        
            rng(state);
        
            tstart = clock;
        
            iRVM = SDec_iRVM(X0,Lf,gamma,mu_max,beta,OPTIONS_iRVM,m,n,q,A,A1,A2,Imat2q);
        
            iRVM_time(jj) = etime(clock,tstart);
        
            iRVM_obj(jj) = iRVM.Theta;

            iRVM_iter(jj) = iRVM.iter;
        
            niter = nnz(iRVM.tnls);
        
            iRVM_nls(jj) = sum(iRVM.tnls)/niter; iRVM_nssn(jj) = sum(iRVM.tsubit)/niter;
        
            %% ********************* RADMM *************************************
        
            OPTIONS_RADMM.maxiter = 1e5;  OPTIONS_RADMM.tol = 1e-7;
        
            OPTIONS_RADMM.printyes = 0;   OPTIONS_RADMM.normA = normA;
        
            eta_RADMM = 1e-1; sigma_RADMM = 5e1;  gamma_RADMM = 1e-8;
        
            rng(state);
        
            tstart = clock;
        
            RADMM = SDec_RADMM(X0,eta_RADMM,sigma_RADMM,gamma_RADMM,OPTIONS_RADMM,m,n,q,A,A1,A2,Imat2q);
        
            admm_time(jj) = etime(clock,tstart);
        
            admm_obj(jj) = RADMM.Theta;

            admm_iter(jj) = RADMM.iter;
        
          %% ********************* RiALM_RGD *************************************
        
            OPTIONS_RiALM.maxiter = 100;  OPTIONS_RiALM.tol = 5e-5;
        
            OPTIONS_RiALM.printyes = 0;   OPTIONS_RiALM.normA = normA;
        
            sigma_RALM = 1.5;  b = 1.5;  err = 1.5;
        
            rng(state);
        
            tstart = clock;
        
            RiALM = SDec_RiALM(X0,sigma_RALM,b,err,OPTIONS_RiALM,m,n,q,A,A1,A2,Imat2q);
        
            iALM_time(jj) = etime(clock,tstart);
        
            iALM_obj(jj) = RiALM.Theta;

            iALM_iter(jj) = RiALM.iter;

            iALM_iter_in(jj) = sum(RiALM.iter_in)/nnz(RiALM.iter_in);
        end

        fprintf('\n**************************************************************************************');
        fprintf('\n**************************************************************************************');
        fprintf('\n Type: %d, n: %d, r: %d', data_type, n, q);
        
        aiRVM_obj = mean(iRVM_obj); aiRVM_time = mean(iRVM_time); aiRVM_nls = mean(iRVM_nls); aiRVM_nssn = mean(iRVM_nssn);

        iRVM_aiter = mean(iRVM_iter);
        
        fprintf('\n Method: %s, obj#: %3.4e, CPU time: %.2f, Nls#: %.2f, Nssn#: %.2f, Iter: %.1f','iRVM',aiRVM_obj,aiRVM_time,aiRVM_nls,aiRVM_nssn,iRVM_aiter);
        
        aiALM_obj = mean(iALM_obj); aiALM_time = mean(iALM_time);

        iALM_aiter = mean(iALM_iter); iALM_aiter_in = mean(iALM_iter_in);
        
        fprintf('\n Method: %s, obj#: %3.4e, CPU time: %.2f, Iter: %.1f, Iter_inner: %.1f','RiALM',aiALM_obj,aiALM_time,iALM_aiter,iALM_aiter_in);
        
        aadmm_obj = mean(admm_obj); aadmm_time = mean(admm_time);

        admm_aiter = mean(admm_iter);
        
        fprintf('\n Method: %s, obj#: %3.4e, CPU time: %.2f, Iter: %.1f','RADMM',aadmm_obj,aadmm_time,admm_aiter);

    end

end

for data_type = [1, 2]
    m = 50; q = 10;  r = 20;
    for n = [400, 600, 700]

        %% *********************** Initialization *************************

        ntest = 10;
        
        iRVM_obj = zeros(1,ntest); iRVM_time = zeros(1,ntest); iRVM_nls = zeros(1,ntest); iRVM_nssn = zeros(1,ntest);
        
        iALM_obj = zeros(1,ntest); iALM_time = zeros(1,ntest); 
        
        admm_obj = zeros(1,ntest); admm_time = zeros(1,ntest);

        iRVM_iter = zeros(1,ntest); admm_iter = zeros(1,ntest);

        iALM_iter = zeros(1,ntest); iALM_iter_in = zeros(1,ntest);
        
        %% ********** Set the size of problems **************************
        
        Imat2q = eye(2*q);  
        
        J2q= [zeros(q,q)  eye(q); -eye(q)  zeros(q)];
        
        for jj = 1:ntest
            %jj
            %% ****************** upload the data ****************************
            randstate = 88888+jj*10000;
            randn('state',double(randstate));
            rand('state',double(randstate));
        
            state = rng;
            
            %% ************ to generate data and an initial X0\in M ************
            
            switch data_type
                
                case 1
                    
                    A = randn(2*n,2*m); normA = norm(A,'fro');  A = A/normA;
                    
                    A1 = A(1:n,:);   A2 = A(n+1:end,:);                     
                    
                case 2
        
                    M0 = symplecticStiefelfactory(n,r,0);
        
                    T = M0.rand();
        
                    C = randn(2*r,2*m); C = C/norm(C,'fro');
        
                    A = T*C; 
        
                    A1 = A(1:n,:); A2 = A(n+1:end,:); normA = norm(A,'fro');
                    
                    % W = randn(2*r,2*r); W = W'*W+0.8*eye(2*r); 
                    % 
                    % E = expm([W(r+1:end,:); -W(1:r,:)]);
                    % 
                    % C = randn(2*r,2*m); C = C/norm(C,'fro');
                    % 
                    % A = [E(1:r,:);zeros(n-r,2*r);E(r+1:end,:);zeros(n-r,2*r)]*C;
                    % 
                    % A1 = A(1:n,:); A2 = A(n+1:end,:); normA = norm(A,'fro');  
            end
        
            M = symplecticStiefelfactory(n,q,0);
        
            X0 = M.rand();
            
            %W = randn(2*q,2*q); W = W'*W+0.1*eye(2*q);  E = expm([W(q+1:end,:); -W(1:q,:)]);
             
            %X0 = [E(1:q,:);zeros(n-q,2*q); E(q+1:end,:);zeros(n-q,2*q)];
            
            % J2nX0 = [X0(n+1:end,:); -X0(1:n,:)];
            % 
            % norm(X0'*J2nX0-J2q,'fro')
            
            %% ********************* iRVM *************************************
        
            OPTIONS_iRVM.maxiter = 5000;  OPTIONS_iRVM.tol = 1e-7;
        
            OPTIONS_iRVM.printyes = 0;    OPTIONS_iRVM.normA = normA;
        
            Lf = 1e-5;  mu_max = 5e+2;  beta = 1e-2;
        
            gamma = 1e-5;
        
            rng(state);
        
            tstart = clock;
        
            iRVM = SDec_iRVM(X0,Lf,gamma,mu_max,beta,OPTIONS_iRVM,m,n,q,A,A1,A2,Imat2q);
        
            iRVM_time(jj) = etime(clock,tstart);
        
            iRVM_obj(jj) = iRVM.Theta;

            iRVM_iter(jj) = iRVM.iter;
        
            niter = nnz(iRVM.tnls);
        
            iRVM_nls(jj) = sum(iRVM.tnls)/niter; iRVM_nssn(jj) = sum(iRVM.tsubit)/niter;
        
            %% ********************* RADMM *************************************
        
            OPTIONS_RADMM.maxiter = 1e5;  OPTIONS_RADMM.tol = 1e-7;
        
            OPTIONS_RADMM.printyes = 0;   OPTIONS_RADMM.normA = normA;
        
            eta_RADMM = 1e-1; sigma_RADMM = 5e1;  gamma_RADMM = 1e-8;
        
            rng(state);
        
            tstart = clock;
        
            RADMM = SDec_RADMM(X0,eta_RADMM,sigma_RADMM,gamma_RADMM,OPTIONS_RADMM,m,n,q,A,A1,A2,Imat2q);
        
            admm_time(jj) = etime(clock,tstart);
        
            admm_obj(jj) = RADMM.Theta;

            admm_iter(jj) = RADMM.iter;
        
          %% ********************* RiALM_RGD *************************************
        
            OPTIONS_RiALM.maxiter = 100;  OPTIONS_RiALM.tol = 5e-5;
        
            OPTIONS_RiALM.printyes = 0;   OPTIONS_RiALM.normA = normA;
        
            sigma_RALM = 1.5;  b = 1.5;  err = 1.5;
        
            rng(state);
        
            tstart = clock;
        
            RiALM = SDec_RiALM(X0,sigma_RALM,b,err,OPTIONS_RiALM,m,n,q,A,A1,A2,Imat2q);
        
            iALM_time(jj) = etime(clock,tstart);
        
            iALM_obj(jj) = RiALM.Theta;

            iALM_iter(jj) = RiALM.iter;

            iALM_iter_in(jj) = sum(RiALM.iter_in)/nnz(RiALM.iter_in);
        end

        fprintf('\n**************************************************************************************');
        fprintf('\n**************************************************************************************');
        fprintf('\n Type: %d, n: %d, r: %d', data_type, n, q);
        
        aiRVM_obj = mean(iRVM_obj); aiRVM_time = mean(iRVM_time); aiRVM_nls = mean(iRVM_nls); aiRVM_nssn = mean(iRVM_nssn);

        iRVM_aiter = mean(iRVM_iter);
        
        fprintf('\n Method: %s, obj#: %3.4e, CPU time: %.2f, Nls#: %.2f, Nssn#: %.2f, Iter: %.1f','iRVM',aiRVM_obj,aiRVM_time,aiRVM_nls,aiRVM_nssn,iRVM_aiter);
        
        aiALM_obj = mean(iALM_obj); aiALM_time = mean(iALM_time);

        iALM_aiter = mean(iALM_iter); iALM_aiter_in = mean(iALM_iter_in);
        
        fprintf('\n Method: %s, obj#: %3.4e, CPU time: %.2f, Iter: %.1f, Iter_inner: %.1f','RiALM',aiALM_obj,aiALM_time,iALM_aiter,iALM_aiter_in);
        
        aadmm_obj = mean(admm_obj); aadmm_time = mean(admm_time);

        admm_aiter = mean(admm_iter);
        
        fprintf('\n Method: %s, obj#: %3.4e, CPU time: %.2f, Iter: %.1f','RADMM',aadmm_obj,aadmm_time,admm_aiter);

    end

end

