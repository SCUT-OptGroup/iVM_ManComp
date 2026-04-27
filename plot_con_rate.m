%% ***************************************************************
% Filename: demo_plot_con_rate
%% ***************************************************************
%% to plot the convergence tate curve for problem
%%
%% **************************************************************
 clear; clear all; clc

close all;

%restoredefaultpath;
addpath(genpath('data'));
addpath(genpath('SparSC'));
addpath(genpath('r21PCA'));
addpath(genpath('Sdecom'));
addpath(genpath('functions'));

%% Set size and lambda for problem 1

n = 300; q = 6; p = 2000;

lambda = 5e-5;

Imat = eye(n);

display = 0;

randstate = 9999;
randn('state',double(randstate));
rand('state',double(randstate));

%% ****************** upload the data ****************************
        
[true_labs, X] = syn_SSC(p,n,q,display,randstate); % synthetic data

simX = func_simlarity(X);   % calculate the similarity matrix

rsum = sum(simX,2);

dvec = 1./rsum.^(1/2);

Lap = Imat - dvec.*(simX.*dvec');

X0 = orth(randn(n,q)); % generate an initial point X0\in M

OPTIONS_iRVM.maxiter = 5000; OPTIONS_iRVM.tol = 1e-12;

OPTIONS_iRVM.printyes = 1; OPTIONS_iRVM.normLap = norm(Lap,'fro');

gamma = 1e-5;  beta = 1e-2; mu_max = 5e+2;

Lf = 5e-1*svds(Lap,1);

%% ********************* iRVM *************************************

iRVM = SSC_iRVM_con_rate(X0,gamma,Lf,mu_max,beta,OPTIONS_iRVM,lambda,Imat,Lap,n,q);

Xcell_p1 = iRVM.Xcell; Objvalue_p1 = iRVM.objvale;

num1 = nnz(Objvalue_p1);

Dist_p1 = zeros(num1,1); Value_p1 = zeros(num1,1);

for i = 1:num1
    Dist_p1(i) = log(norm(Xcell_p1{i}-Xcell_p1{num1},'fro'));
    Value_p1(i) = log(Objvalue_p1(i)-Objvalue_p1(num1));
end

%% ************************************************************************************
%% Set size and lambda for problem 2

m = 50; n = 1000; q = 5;  

lambda = 2.05;   rho = 0.5;

B = randn(m,n);
    
B = B - repmat(mean(B,1),m,1);

B = normc(B);

A = B'*B;  E = ones(q)-eye(q);

normA = norm(A,'fro');
  
Lf = 5e-1*svds(B,1)^2;  %

X0 = orth(randn(n,q)); % generate an initial point X0\in M

%% ********************* iRVM *************************************

OPTIONS_iRVM.maxiter = 10000;  OPTIONS_iRVM.tol = 1e-12;

OPTIONS_iRVM.printyes = 1;   OPTIONS_iRVM.normA = normA;

iRVM_gamma = 1e-5;  mu_max = 5e2;  beta = 1e-2;

iRVM = r21PCA_iRVM_con_rate(X0,Lf,iRVM_gamma,mu_max,beta,OPTIONS_iRVM,lambda,rho,m,n,q,E,A,B);

Xcell_p2 = iRVM.Xcell; Objvalue_p2 = iRVM.objvale;

num2 = nnz(Objvalue_p2);

Dist_p2 = zeros(num2,1); Value_p2 = zeros(num2,1);

for i = 1:num2
    Dist_p2(i) = log(norm(Xcell_p2{i}-Xcell_p2{num2},'fro'));
    Value_p2(i) = log(Objvalue_p2(i)-Objvalue_p2(num2));
end

%% ************************************************************************************
%% Set size and lambda for problem 3

m = 50; n = 500;  q = 10;  r = 20;

lambda = 1e-4;

Imat2q = eye(2*q);  

J2q= [zeros(q,q)  eye(q); -eye(q)  zeros(q)];

A = randn(2*n,2*m); normA = norm(A,'fro');  A = A/normA;
            
A1 = A(1:n,:);   A2 = A(n+1:end,:);

M = symplecticStiefelfactory(n,q,0);

X0 = M.rand();

%% ********************* iRVM *************************************
    
OPTIONS_iRVM.maxiter = 1500;  OPTIONS_iRVM.tol = 1e-9;

OPTIONS_iRVM.printyes = 1;    OPTIONS_iRVM.normA = normA;

Lf = 1e-5;  mu_max = 5e+2;  beta = 1e-2;

gamma = 1e-5;

iRVM = SDecFL1_iRVM_con_rate(X0,Lf,gamma,mu_max,beta,lambda,OPTIONS_iRVM,m,n,q,A,A1,A2,Imat2q);

Xcell_p3 = iRVM.Xcell; Objvalue_p3 = iRVM.objvale;

num3 = nnz(Objvalue_p3);

Dist_p3 = zeros(num3,1); Value_p3 = zeros(num3,1);

for i = 1:num3
    Dist_p3(i) = log(norm(Xcell_p3{i}-Xcell_p3{num3},'fro'));
    Value_p3(i) = log(Objvalue_p3(i)-Objvalue_p3(num3));
end


%% Plot the figures
num1 = 550; num2 = 6500; num3 = 1400;

figure(1)          % define figure
set(gcf,'position',[50 50 1333 360])
subplot(1,3,1);     % subplot(x,y,n)
plot(1:num1,Dist_p1(1:num1),'r-','linewidth',2);
xlabel('iter');
ylabel('$\log\|x^k - x^*\|_F(\log|\Theta(x^k) - \Theta(x^*)|)$', 'Interpreter', 'latex');
xlim([0 num1+50]);
hold on

subplot(1,3,1);     % subplot(x,y,n)
plot(1:num1,Value_p1(1:num1),'b-','linewidth',2);
hold on

subplot(1,3,2);     % subplot(x,y,n)
plot(1:num2,Dist_p2(1:num2),'r-','linewidth',2);
xlabel('iter');
ylabel('$\log\|x^k - x^*\|_F(\log|\Theta(x^k) - \Theta(x^*)|)$', 'Interpreter', 'latex');
xlim([0 num2+500]);
hold on

subplot(1,3,2);     % subplot(x,y,n)
plot(1:num2,Value_p2(1:num2),'b-','linewidth',2);
hold on

subplot(1,3,3);     % subplot(x,y,n)
plot(1:num3,Dist_p3(1:num3),'r-','linewidth',2);
xlabel('iter');
ylabel('$\log\|x^k - x^*\|_F(\log|\Theta(x^k) - \Theta(x^*)|)$', 'Interpreter', 'latex');
xlim([0 num3+100]);
hold on

subplot(1,3,3);     % subplot(x,y,n)
plot(1:num3,Value_p3(1:num3),'b-','linewidth',2);
xlabel('iter');

% %% Plot the figures
% figure(1)          % define figure
% set(gcf,'position',[50 50 1000 450])
% subplot(2,3,1);     % subplot(x,y,n)
% plot(1:num1,Dist_p1,'r-','linewidth',2);
% xlabel('Iter');
% ylabel('$\log\|x^k - x^*\|_F$', 'Interpreter', 'latex');
% xlim([0 num1]);
% hold on
% 
% figure(1)          % define figure
% set(gcf,'position',[50 50 1000 450])
% subplot(2,3,4);     % subplot(x,y,n)
% plot(1:num1,Value_p1,'r-','linewidth',2);
% xlabel('Iter');
% ylabel('$\log|\Theta(x^k) - \Theta(x^*)|$', 'Interpreter', 'latex');
% xlim([0 num1]);
% hold on
% 
% figure(1)          % define figure
% set(gcf,'position',[50 50 1000 450])
% subplot(2,3,2);     % subplot(x,y,n)
% plot(1:num2,Dist_p2,'b-','linewidth',2);
% xlabel('Iter');
% ylabel('$\log\|x^k - x^*\|_F$', 'Interpreter', 'latex');
% xlim([0 num2]);
% hold on
% 
% figure(1)          % define figure
% set(gcf,'position',[50 50 1000 450])
% subplot(2,3,5);     % subplot(x,y,n)
% plot(1:num2,Value_p2,'b-','linewidth',2);
% xlabel('Iter');
% ylabel('$\log|\Theta(x^k) - \Theta(x^*)|$', 'Interpreter', 'latex');
% xlim([0 num2]);
% hold on
% 
% figure(1)          % define figure
% set(gcf,'position',[50 50 1000 450])
% subplot(2,3,3);     % subplot(x,y,n)
% plot(1:num3,Dist_p3,'c-','linewidth',2);
% xlabel('Iter');
% ylabel('$\log\|x^k - x^*\|_F$', 'Interpreter', 'latex');
% xlim([0 num3]);
% hold on
% 
% figure(1)          % define figure
% set(gcf,'position',[50 50 1000 450])
% subplot(2,3,6);     % subplot(x,y,n)
% plot(1:num3,Value_p3,'c-','linewidth',2);
% xlabel('Iter');
% ylabel('$\log|\Theta(x^k) - \Theta(x^*)|$', 'Interpreter', 'latex');
% xlim([0 num3]);
