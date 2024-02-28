function [R2eff, p, F, Reffall] = r2cpmg_fourstate_ppm_constrain(nu,P)

% Calculates the relaxation decay in a CPMG experiment.
% This program calculates the full solution to the 
% Bloch-McConnell equations, including exchange.
% See Protein NMR Spectroscopy by Cavanagh et al, pp 291ff.
% Copyright 1999. Kristofer Modig (kristofer.modig@bpc.lu.se), 
% Mikael Akke (mikael.akke@bpc.lu.se)
% 120331
% System description: 4 states, two pathways (select-fit, induced-fit)
% Pc (1), Po (2), PcL (3), PoL (4) (closed, open, closed-ligand, open-ligand)
%
% This version uses shifts in ppm and translates to Hz from transmitter
% frequency.
%
% This version also constrains the value of X to be max 1e5, 
% via the parameter B, defined as 
%
%   Y/1000 = 0 + 100*(exp(B)/(1 + exp(B))
%
%
% Parameters:
%     
% Supposed to be fitted
%   P(1) = kon/1000      % on-rate in mM^-1.s^-1 (assumed equal to both 
%                        % affinity states). 
%   P(2) = X/100         % factor by which binding to lo-affinity conf
%                        % is weaker than binding to hi-affinity conf, X > 1.
%   P(3) = B             % Parameter used to calculate Y/1000 = 0 + 100*(exp(B)/(1 + exp(B)).
%                        % Y is the factor by which open-to-close rate is faster 
%                        % without ligand bound, ie Y = k21/k43;
%   P(4) = ppm2          % Offset of state 2 (open, conf selection path) (ppm)
%   P(5) = ppm3          % Offset of state 3 (closed.lig, induced fit path) (ppm)
%   P(6) = R20            % R2 in absense of exchange. Assumed equal for all
%                        % states.
% Known
%   P(7) = k12 (= 20)    % closed (lo-affinity) to open (hi-affinity) 
%                        % (free protein), known from exchange in free state
%   P(8) = k21 (= 600)   % open (hi-affinity) to closed (lo-affinity) (free protein), known from exchange in free state
%   P(9) = Ptot (= 0.5)  % Total protein conc in mM (note that this really 
%                        % should change with ligand additions!)
%   P(10) = L             % Free ligand conc in mM. (Carl's data cover 0.01-0.09 mM.)
%   P(11) = Kd (= 0.228)    % Kd in mM
%   P(12) = ppm1
%   P(13) = ppm4
%   P(14) = Tc (=0.1)    % Constant relaxation time (s). Same for all
%                        % frequencies (R2effs)!
%   P(15) = sfrq         % Transmitter frequency (MHz)
%
% Assumption that kon1 = kon2 (verified for Gal3C from exchange
% involving lac, L2 and L3, where kon is independent on Kd!
%
% The program may return the equilibrium populations p and the fluxes in 
% the vector F, with:
%
% Fluxes through CS and IF (ligand binding):
%    F(1) = Fcs
%    F(2) = Fif
% Reverse fluxes through CS and IF (ligand release):
%    F(3) = FRcs
%    F(4) = FRif


% Go from ppm to Hz for shifts
Pin([4 5 13 14]) = P([4 5 12 13]).*P(15);

% Rescale some parameters
Pin(1) = P(1)*1000;
Pin(2) = P(2)*100;

% Limits for Y (0 < Y/1000 < 100)
Pin(3) = P(3); %100*exp(P(3))/(1 + exp(P(3)))*1000;

% Set both R20 equal
Pin(6) = P(6);
Pin(7) = P(6);

% The rest of the parameters
Pin([8 9 10 11 12 15]) = P([7 8 9 10 11 14]);

% Everything should be > 0
Pin = abs(Pin);

% Call CPMG-R2 routine
if nargout == 1,
    R2eff = r2cpmg_fourstate_diag(nu,Pin);
else
    [R2eff, p, F, Reffall] = r2cpmg_fourstate_diag(nu,Pin);
end


