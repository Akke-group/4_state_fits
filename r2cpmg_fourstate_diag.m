function [R2eff, p, F, Reffall] = r2cpmg_fourstate_diag(nu,P)


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
% This version calculates the exponential matrix using diagonalization
% berofore the CPMG-loop, which turns out to be 10 times faster. Might not
% be as accurate, but I've found no differences.
%
% Parameters:
%     
% Supposed to be fitted
%   P(1) = kon           % on-rate in M^-1.s^-1 (assumed equal to both 
%                        % affinity states)
%   P(2) = X             % factor by which binding to lo-affinity conf 
%                        % is weaker than binding to hi-affinity conf, X > 1.
%   P(3) = Y             % factor by which open-to-close rate is faster 
%                        % without ligand bound, ie Y = k21/k43;
%   P(4) = off2          % Offset of state 2 (open, conf selection path) (Hz)
%   P(5) = off3          % Offset of state 3 (closed.lig, induced fit path) (Hz)
%   P(6) = R20P          % R2 in absense of exchange for Pc and Po states
%   P(7) = R20PL         % R2 in absense of exchange for PcL and PoL states
%
% Known
%   P(8) = k12          % closed (lo-affinity) to open (hi-affinity) 
%                       % (free protein), known from exchange in free state
%   P(9) = k21          % open (hi-affinity) to closed (lo-affinity) (free protein),
%                       % known from exchange in free state
%   P(10) = Ptot        % Total protein conc 
%   P(11) = L           % Free ligand conc 
%   P(12) = Kd          % Kd 
%   P(13) = off1
%   P(14) = off4
%   P(15) = Tc          % Constant relaxation time (s). Same for all
%                        % frequencies (R2effs)!
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


% Take care of the parameters
kon = P(1);
X = P(2); 
Y = P(3);
off2 = P(4);
off3 = P(5);
R20P = P(6);
R20PL = P(7);
k12 = P(8);
k21 = P(9);
Ptot = P(10); 
L = P(11);
Kd = P(12);
off1 = P(13);
off4 = P(14);
Tc = P(15); 

% The two on-rates are assumed equal
kon1 = kon;  % on-rate in M^-1.s^-1 (kon1.L = k13)
kon2 = kon;  % on-rate in M^-1.s^-1 (kon2.L = k24)


% ================ Define system parameters ==================================
% First, define system size
systemsize = 4; % number of exchanging sites

% Define number of CPMG delays, minimum and maximum delays, and relaxation delay.
% Each CPMG delay is sampled by the constant-relaxation time approach:
% Reff = -(1/T)*ln(Int(T)/Int(0))
% Note the definition used here: tau =  spacing between 180' pulses
% One CPMG block = 4*tau and contains two 180' pulses.
% The CPMG frequency (number of 2*pi rotations per second) is thus
% 1/(4*tau)

% Calculate tau from cpmg frequencies
tau = 1./(4*nu);
nr_tau = length(tau); % number of experimental points

% Use round below: Warning: Tc must be checked before hand!
nr_blocks = round(Tc./(4.*tau));  % Each block needs 4*tau


% Concentration of free protein
Pfree = Ptot/(1+L/Kd);

% Concentration of bound protein
PL = Ptot/(1+Kd/L);

% Calculate rate constants from input parameters
koff4 = kon2*(k12+k21/X)/(k12+k21)*Kd;   % off-rate from open state 
               % (hi-affinity) (koff4 = k42), (koff4 approx. 300 from NMR)
koff3 = X*koff4;  % off-rate from state 3 (lo-affinity) (koff3 = k31), X ? 1
k34 = X*k12/Y;      % k43 is given by other rates if kon1=kon2, X = ratio of koff3 and koff4
k43 = k21/Y; % k43 is presumably slower than k21, ie Y > 1. 


% Define equilibrium concentrations.
p = zeros(systemsize,1); % column vector
p(1) = k21/(k12+k21)*Pfree/Ptot; % ground-state apo, closed (lo-affinity)
p(2) = k12/(k12+k21)*Pfree/Ptot; % high-energy state apo, open (hi-affinity)
p(3) = k43/(k34+k43)*PL/Ptot;    % non-specific encounter complex, closed.L
p(4) = k34/(k34+k43)*PL/Ptot;    % specific, hi-affinity complex, open.L

%if nargout >= 3
    % Fluxes through CS and IF (ligand binding):
    F(1) = 1/(1/(k12*p(1)*Ptot)+1/(kon2*p(2)*Ptot*L));
    F(2) = 1/(1/(kon1*p(1)*Ptot*L)+1/(k34*p(3)*Ptot));
    % Reverse fluxes through CS and IF (ligand release):
    F(3) = 1/(1/(k21*p(2)*Ptot)+1/(koff4*p(4)*Ptot));
    F(4) = 1/(1/(koff3*p(3)*Ptot)+1/(k43*p(4)*Ptot));
%end


% Give R2 relaxation rates for each site. Assume that R2 is the same for
% all states
R2 = zeros(systemsize,1);
R2(1) = R20P;
R2(2) = R20P;
R2(3) = R20PL;
R2(4) = R20PL;


% Define offsets of the states, in Hz relative to the rotating frame reference.
offset = zeros(systemsize,1);
offset(1) = (off1-off2*p(2)/(p(2)+p(1)))/(p(1)/(p(2)+p(1))); 
offset(2) = off2;
offset(3) = off3;
offset(4) = (off4-off3*p(3)/(p(3)+p(4)))/(p(4)/(p(3)+p(4)));

if L == 0
    offset(4) = off4;
end


% Define reaction rates for each transition.
% Note kij is the rate going from i to j.
% The elements of the rate matrix K are defined: Kij is the rate going from
% j to i. We take care of this below when we define K.
kmat = zeros(systemsize) ;
kmat(1,2) = k12 ;
kmat(1,3) = kon1*L ;
kmat(1,4) = 0 ;
kmat(2,1) = k21 ;
kmat(2,3) = 0 ;
kmat(2,4) = kon2*L ;
kmat(3,1) = koff3 ;
kmat(3,2) = 0 ;
kmat(3,4) = k34 ;
kmat(4,1) = 0 ;
kmat(4,2) = koff4 ;
kmat(4,3) = k43 ;

K = zeros(systemsize);
for j = 1:systemsize,
	for l = 1:systemsize
		K(j,l) = kmat(l,j) ;
	end;
end;
for j = 1:systemsize,
	K(j,j) = -(sum(K(:,j))) ;
end;

% Define R relaxation matrix
R = diag(R2);

% Define W chemical shift (referenced to rotating frame) matrix
W = zeros(systemsize);
for j = 1:systemsize,
	W(j,j) = 2*pi*offset(j);
end;

Liouvillian = 1i*W - R + K; % this is the Liouvillian


[U,D] = eig(Liouvillian);
invU = inv(U);

Reffall = zeros(nr_tau,systemsize);
for k = 1:nr_tau,
    % Define initial magnetizations
    % initial magnetizations = equilibrium concentrations
	M = p; % column vector
    time = tau(k);
    for j = 1:nr_blocks(k),
       % One CPMG block = 4*tau and contains two 180' pulses.
	   M = (U*diag(exp(diag(D.*time)))*invU)*M;
	   M = conj(M) ; % 180' pulse (infinitely short).
	   M = (U*diag(exp(diag(D.*(2*time))))*invU)*M; 
	   M = conj(M) ; % 180' pulse (infinitely short).
	   M = (U*diag(exp(diag(D.*time)))*invU)*M; 
    end;
    Reffall(k,:) = -(1/Tc)*log(real(M)./real(p));
end;

R2eff = Reffall(:,1)';



