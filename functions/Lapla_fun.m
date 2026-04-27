%% ***************************************************************
%% filename: Lapla_fun
%% To generate the normalized Laplacian matrix L for a given data matrix in_X
%% **************************************************************

function [n,C,Lap] = Lapla_fun(si_X,true_labs)

n = length(true_labs); C = max(true_labs);

dvec = 1./sum(si_X,2).^(1/2);

Lap = eye(n) - (dvec.*si_X).*dvec';

end