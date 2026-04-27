%% ***************************************************************
%% Filename: prox_max
%% ***************************************************************
%  Compute the proximal mapping 
%
%  prox_{theta}(Z)
%
% where theta(Z) = lambda*sum_{i,j}max{Z_{i,j},0}
%
%% **************************************************************

function [proxZ,eML1,pidx] = prox_max(Z,lambda)

pidx1 = Z<=0; pidx2 = Z>=lambda;

proxZ = pidx1.*Z + pidx2.*(Z-lambda);

if nargout > 1

    pZidx = proxZ > 0;

 eML1 = 0.5*norm(proxZ-Z,'fro')^2 + lambda*sum(sum(pZidx.*proxZ));
 
 pidx = pidx1 + pidx2;
   
end

end