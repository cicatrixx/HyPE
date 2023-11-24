function [Pbias, NSE]= evalStats(sim, obs)
% PBIAS and NSE as defined in Moriasi et al. 2007 Eqn 1 and 2
% 
% "The optimal value of PBIAS is 0.0, with low-magnitude values indicating
% accurate model simulation. Positive values indicate model underestimation bias, and negative values indicate model overestimation bias (Gupta et al., 1999)."
% "Model simulation can be judged as satisfactory if NSE > 0.50 and RSR <
% 0.70, and if PBIAS  25% for streamflow, PBIAS +/- 55% for sediment."
  
    Pbias=sum(obs-sim)./sum(obs)*100; 
    NSE=1-sum((obs-sim).^2)./sum((obs-mean(obs)).^2);
end