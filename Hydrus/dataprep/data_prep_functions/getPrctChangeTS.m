function dataout=getPrctChangeTS(datain)
%% Evaluates % change between timesteps
% datain can be a column vector or matrix w columns indicating different TS
% dataout has one less row than datain because change in calculated from
% 2nd timestep until the last.

    dataout=100*diff(datain)./datain(1:end-1,:);
end