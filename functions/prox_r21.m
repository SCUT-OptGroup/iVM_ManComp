%% ****************************************************************
%% filename: prox_r21
%% ****************************************************************
%  to compute the proximal mapping and Moreau-envelope of g(X)=lambda*||X||_{2,1}
%  where ||X||_{2,1} denotes the row L21-norm of X
%% nr: the number of rows of Z

function [proxZ,ZldZ,tempv,eMr21] = prox_r21(Z,lambda,nr)

tempv = zeros(nr,1); ldZ = tempv;

Zrnorm = vecnorm(Z,2,2);

pidx = Zrnorm>lambda;

Zpnorm = Zrnorm(pidx);

tempv(pidx) = 1-lambda./Zpnorm; 

proxZ = Z.*tempv;   % multiply each row Zi of Z with u(i)

eMr21 = 0.5*norm(Z-proxZ,'fro')^2+lambda*sum(vecnorm(proxZ,2,2));

if nargout > 1

    ldZ(pidx) = lambda./(Zpnorm.^3);

    ZldZ = Z.*ldZ;
end

end