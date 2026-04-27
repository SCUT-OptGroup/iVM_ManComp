%% *******************************************************************
%  filename: retr_sp
%
%% *******************************************************************
%% Compute the Cayley retraction on the symplectic stiefel manifold
%%
%% ***************************************************************
%%  
function val = retr_sp(X,V,XJ2q,XJJXt,Imat2q,n,q)

H = V - XJJXt*V;

HJ2q = [-H(:,q+1:end)  H(:,1:q)];

J2nH = [H(n+1:end,:); -H(1:n,:)];

J2nV = [V(n+1:end,:); -V(1:n,:)];

W = 0.25*HJ2q'*J2nH - 0.5*XJ2q'*J2nV + Imat2q;

val = (H+2*X)/W - X;
end
%Z = inv(W')*(H+2*X)';

%Z = linsolve(W',(H+2*X)');