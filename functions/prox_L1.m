%% ***************************************************************
%% Filename: prox_L1
%% ***************************************************************
%  Compute the proximal mapping 
%
%  prox_{theta}(Z)
%
% where theta(Z) = lambda*||Z||_1
%
%% **************************************************************

function [proxZ,eML1,pidx] = prox_L1(Z,lambda)

absZ = abs(Z);

proxZ = sign(Z).*max(absZ-lambda,0);

if nargout>=2

 eML1 = 0.5*norm(proxZ-Z,'fro')^2 + lambda*norm(proxZ(:),1);
 
 pidx = absZ>=lambda;
   
end

end