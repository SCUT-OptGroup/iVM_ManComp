%% ***************************************************************
% filename: retr_st
%% ***************************************************************
% The retraction function of Stiefel manifold supports for 4 types of retraction method:
% Exponential mapping (1),
% Polar decomposition (2),
% QR decomposition (3), and
% Cayley transform (4)
% Inputs:
%   Xk: the original point on Stiefel manifold
%   VK: the tangent vector in T_{Xk}M
%   method: the retraction method selected, default: Cayley
%   Xnew: the new point on Stiefel manifold

function Xnew = retr_st(Xk,Vk,r,method)

if nargin < 3 || isempty(method)
    method = 4;
end

if method == 1
    Imat = [eye(r);zeros(r)];
    XktVk = Xk'*Vk;
    Q = qr(-Vk+Xk*XktVk,0);
    Xnew = [Xk, Q]*expm([-XkVk -R'; R  zeros(r)])*Imat;
    
elseif method == 2
    
    [P,D] = eig(eye(r)+ Vk'*Vk);
    
    d = diag(D).^(-1/2);
    
    Xnew = (Xk+Vk)*(P*diag(d)*P');
    
elseif method == 3
    
    [Q,RR] = qr(Xk+Vk,0);
    
    diagRR = sign(diag(RR));
    
    ndr = diagRR < 0;
    
    if nnz(ndr)>0
        Q = Q*spdiags(diagRR,0,r,r);
    end
    Xnew = Q;
    
elseif method == 4
    
    n = size(Xk,1);
    
    Imat = eye(n);
    
    IXXtVXt = (Imat-0.5*(Xk*Xk'))*(Vk*Xk');
    
    W = IXXtVXt - IXXtVXt';
    
    Xnew =linsolve(Imat-0.5*W, (Xk+0.5*W*Vk));
    
elseif method == 5
    
    PY = Xk+Vk;
    
    PYq = PY'*PY;
    
    PYq = 0.5*(PYq+PYq');
    
    [U,SIGMA,S] = eig(PYq);
    
    SIGMA =diag(SIGMA);
    
    Xnew = PY*(U*diag(sqrt(1./SIGMA))*S');
end

%       [U,SIGMA,S] = svd(PY);   
%     
%       SIGMA = diag(SIGMA);    
%      
%       Xnew = PY*(U*diag(sqrt(1./SIGMA))*S');