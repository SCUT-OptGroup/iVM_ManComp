%% ***************************************************************
%% Filename: prox_LF
%% ***************************************************************
%  Compute the proximal mapping of the function $f(Z)=lambda*||Z||_F$ 
%
%% **************************************************************

function [proxZ,eMLF,pidx,JacLF] = prox_LF(Z,lambda)

Zfro = norm(Z,'fro');

pidx = Zfro>=lambda;

proxZ = max(1-lambda/Zfro,0)*Z;

if pidx==1
  
   JacLF.const = 1/Zfro;
   
   JacLF.mat = Z;   
else
   JacLF.const = 0;    
end

if nargout > 1

eMLF = 0.5*norm(proxZ-Z,'fro')^2 + lambda*norm(proxZ,'fro');

end

end

