%% ***************************************************************
% Filename: plot_gam_SDec
%% ***************************************************************
% to test the performance of iRVM with difference gamma for solving the problem
%%
%% min <XX',L> + lambda||XX'||_1, s.t. X'X=I_q
%%
%% **************************************************************
clear;

restoredefaultpath;
addpath(genpath('Sdecom'));
addpath(genpath('functions'));

% Set problem size and lambda

m = 50; n = 200; q = 5;  r = 20;

Imat2q = eye(2*q);

OPTIONS_iRVM.maxiter = 2000;

OPTIONS_iRVM.printyes = 1;

OPTIONS_iRVM.tol = 1e-7;

OPTIONS_iRVM.normA = 1;

Lf = 1e-5;

mu_max = 5e2; beta = 1e-2;

gam_list = [1e-6 5e-6 1e-5  5e-5  1e-4  5e-4  1e-3  5e-3  1e-2  5e-2  1e-1];

ngam = length(gam_list);

ntest = 10;

%% Zeros array to collect each iteration

obj = zeros(ngam,ntest);  time = zeros(ngam,ntest);

for k = 1 : ngam
    
    gamma = gam_list (k);
    
    for i = 1 : ntest
        
        k
        i
     %% ****************** upload the data ************************
        
        randstate = 8888+i*2000;
        randn('state',double(randstate));
        rand('state',double(randstate));        
        
        %% ************ to generate data and an initial X0\in M ************

        A = randn(2*n,2*m); normA = norm(A,'fro');  A = A/normA;
            
        A1 = A(1:n,:);   A2 = A(n+1:end,:); 

        
        M = symplecticStiefelfactory(n,q,0);

        X0 = M.rand();
        
        tstart = clock;
        
        iRVM = SDec_iRVM(X0,Lf,gamma,mu_max,beta,OPTIONS_iRVM,m,n,q,A,A1,A2,Imat2q);
        
        time(k,i)  = etime(clock,tstart);
        
        obj(k,i) = iRVM.Theta; 
        
    end
end
aobj = mean(obj,2); atime = mean(time,2);

%% Plot the figures
figure(1)           % define figure
set(gcf,'position',[50 50 1000 370])
subplot(1,2,1);     % subplot(x,y,n)x表示显示的行数，y表示列数，n表示第几幅图片
plot(gam_list,aobj,'r-','linewidth',2);
xlabel('$\overline{\gamma}$', 'Interpreter', 'latex');
ylabel('obj');

subplot(1,2,2);
plot(gam_list,atime,'b-','linewidth',2);
xlabel('$\overline{\gamma}$', 'Interpreter', 'latex');
ylabel('time (s)');
