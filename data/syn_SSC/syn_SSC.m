%% ***************************************************************
%% Filename: syn_SSC
%% **************************************************************
function [true_lable,Y] = syn_SSC(p,n,C,tag,randstate)

if ~isempty(randstate)
    
    randn('state',double(randstate));
    
    rand('state',double(randstate));
end

%% Generate the center points

Zcent = zeros(2,C);

for j=1:C

    a = cos(2*j*pi/C); b = sin(2*j*pi/C);
    
    Zcent(:,j)=10+10*[a;b];
end


%% Randomly generate the number of points in each class

while 1
    
    C_num = randi([1,n],1,C);
    
    nq = n - C;
    
    C_num = round(nq*C_num/sum(C_num));
    
    [a, idxa] = max(C_num); [b, idxb] = min(C_num);
    
    Num = sum(C_num);
    
    if  Num > nq
        
        C_num(idxa) = a - (Num - nq);
        
    elseif Num < nq
        
        C_num(idxb) = b + (Num - nq);
        
    end
    
    if sum(C_num) == nq
        
        break;
    end
    
end
C_num = C_num + ones(1,C);

Z = zeros(2, n);

for i = 1 : C
    
    if i==1
        
        tempv = Zcent(:,i);
        
        n1 = C_num(i);

        Z(:,1:n1)= repmat(tempv,1,n1);
    
    else
        
        tempv = Zcent(:,i);
        
        Z(:,1+sum(C_num(1:i-1)):sum(C_num(1:i)))= repmat(tempv,1,C_num(i));    
    end
end

Z = Z + randn(2,n);

%% true lable

true_lable =[];

for i = 1:C
    
    true_lable =[true_lable; i*ones(C_num(i),1)];
end

%% draw the scatter plot of the synthetic data
if tag > 0
    j = 1; k = C_num(1);
    for i=1:C-1
        scatter(Z(1,j:k),Z(2,j:k));
        hold on
        j = j + C_num(i); k = k + C_num(i+1);
    end
    scatter(Z(1,j:k),Z(2,j:k));
end


%% Project the data Z to a p-dimensional space

Z = Z';

q = 0.1*p;

P = [rand(2,q)  zeros(2, p-q)];

X = Z*P;

X1 = X + 1.5*randn(n,p);     % adding independent Gaussian noise

Y = X1.*(binornd(1,1-exp(-0.01*X1.^2)));

end
