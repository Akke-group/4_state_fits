% fits the relaxation decay in a CPMG experiment.
% This program calculates the full solution to the
% Bloch-McConnell equations, including exchange.
% See Protein NMR Spectroscopy by Cavanagh et al, pp 291ff.
% Copyright 1999. Kristofer Modig (kristofer.modig@bpc.lu.se),
% Mikael Akke (mikael.akke@bpc.lu.se)
% 120331
% System description: 4 states, two pathways (select-fit, induced-fit)
% Pc (1), Po (2), PcL (3), PoL (4) (closed, open, closed-ligand, open-ligand)
%
%
% Parameters:
%
% Supposed to be fitted
%   kon           % on-rate in M^-1.s^-1 (assumed equal to both
%                        % affinity states)
%   X             % factor by which binding to lo-affinity conf
%                        % is weaker than binding to hi-affinity conf, X > 1.
%   Y             % factor by which open-to-close rate is faster
%                        % without ligand bound, ie Y = k21/k43;
%   off2          % Offset of state 2 (open, conf selection path) (ppm)
%   off3          % Offset of state 3 (closed.lig, induced fit path) (ppm)
%   R20P          % R2 in absense of exchange for Pc and Po states
%   R20PL         % R2 in absense of exchange for PcL and PoL states
%
% Known
%   k12          % closed (lo-affinity) to open (hi-affinity)
%                       % (free protein), known from exchange in free state
%   k21          % open (hi-affinity) to closed (lo-affinity) (free protein),
%                       % known from exchange in free state
%   Ptot        % Total protein conc
%   L           % Free ligand conc
%   Kd          % Kd
%   off1        % Offset of state 1 (ppm)
%   off4        % Offset of state 4 (ppm)
%   Tc          % Constant relaxation time (s). Same for all
%                        % frequencies (R2effs)!
%
%
% It takes several CPMG data sets at different ligand concentrations concentrations.
% In addition, chemical shifts of the apo state and ligand saturated protein state is needed.


Pin = [40000, 2, 0.3] % initial guesses of parameters Kon/1000 rho_off/100 and rho_close
%


makeplot.make = 'y';
makeplot.print = 'y';
makeplot.dat2file = 'y';
makeplot.print_directory = 'FinalMC';
makeplot.plot_dir = '';
freeze = [];


% If Gridsearch have been performed load the saved Matlab data structure and use best grid point as input to fitting
in = load('');
for i = 1:length(in.fitstats);
    chi2(i) = in.fitstats{i}.chi2;
end
[s, Isort] = sort(chi2);
Ifit = Isort(1);
gridsearch = 0;
%Ifit = 1;

%% Load data
cpmgfiles = {'Gal3C_apo_800.cfi'}; % list of cmpg data files


T = [28]; % temperature in C

for i = 1:length(cpmgfiles),
    [R2{i}, R2sig{i}, nucp{i}, resname{i}, carrier_freq{i}] = ...
        cpf_readdata(['cpmg_all_data/', cpmgfiles{i}]);
end

% Total ligand konc
totligC = [];

totligC(exclude) = [];
% Free ligand konc
freeligC = [0, 0.00936, 0.01904, 0.02905, 0.060865, 0.083463, 0.09773, 0.710];

freeligC(exclude) = [];

% parameters
gal3C = 0.5; % protein concentration in mM
Kd = 0.228; % Dissociation constant in mM

%% Load the files containing shifts for apo and holo forms
fid = fopen('gal3_apo_dhnn.txt');
shiftsApo = textscan(fid, '%d%f%f', 'CommentStyle', '#');
fclose(fid);

fid = fopen('gal3_lac_dhnn.txt');
shiftsLac = textscan(fid, '%d%f%f', 'CommentStyle', '#');
fclose(fid);

%% Create huge data vectors for joint fits.
% Read data one file at a time and extract data for the residues showing
% dispersions (as marked in showsDisps).
%
% At the same time, we create the parameter vector allP. It has the
% following structure:
%   Parameters 1 to 3 are fitted parameters that are common to all
%   data sets: kon, X and Y. Matched to P(1:3) in r2cpmg_fourstate_ppm.
%
%   Parameters 4 to 7 are known parameters (not fitted for) that are common
%   to all data sets: k12, k21, Ptot, Kd. Matched to P(7, 8, 9, 11).
%
%   Then comes parameters that are specific to each dispersion, but are not
%   residue specific: L, Tc and sfrq. Matched to P(10, 14, 15)
%
%   Now, we have parameters that are residue specific, but do not vary
%   over the different data sets. They are added residue by residue in the
%   following order:
%      ppm2 (fitted)
%      ppm3 (fitted)
%      ppm1 (known, apo)
%      ppm4 (known, holo)
%   Matched to P(4, 5, 12, 13)
%
%   Finally, we have R20 that is different for each residue and for each
%   sfrq, but assumed equal at all ligand concentrations. Matched to P(6).

% Constrains Y to max 1e5 and min 0
YconstrainMax = 100;
y = 0.1;
allP = Pin;

allP(4) = 33.8; % k12
allP(5) = 735; % k21
allP(6) = gal3C; % Ptot
allP(7) = Kd;
allP(8) = 0.06; % Tc
I_lig = length(freeligC);
Is = length(allP) + 1;
allP(Is:Is+I_lig-1) = freeligC; % Free ligand concentrations


freeze = [freeze, 4:length(allP)];


ResToFit = [144, 147, 155, 174, 176, 183, 184, 187];

%% Res specific starting parameters

gridBestFitRed = zeros(250, 2);


% Chemical shifts of states 2 and 3
gridBestFitRed(144, :) = [117.73, 117.80];
gridBestFitRed(147, :) = [121.80, 122.0];
gridBestFitRed(155, :) = [121.30, 121.40];
gridBestFitRed(173, :) = [123.70, 123.70];
gridBestFitRed(174, :) = [119.3060, 119.6];
gridBestFitRed(176, :) = [130.660, 130.55];
gridBestFitRed(183, :) = [124.50, 124.30];
gridBestFitRed(184, :) = [125.4010, 126.0];

gridBestFitRed(187, :) = [121.71, 120.82];
gridBestFitRed(236, :) = [127.4, 127.3];


error_shift2 = [0.5 * ones(length(ResToFit), 1)];


nlopt = {};

nlopt.pbounds = [2, 0.3, inf; ...
    3, 0.001, inf];


X = {};
Y = {};
Ysig = {};
GeneralNfit = [];
R20_freeze = [];
%
% YconstrainMax = 100;
% y = 10000;
% allP = [];
% allP(1) = 2000/1000;  % kon
% allP(2) = 400/100;   % X
% allP(3) = y/1000;   % Y
% allP(3) = fzero(@(B) YconstrainMax*exp(B)/(1 + exp(B))-allP(3),0); % Constrain Y (see function)
% allP(4) = 33.8;  % k12
% allP(5) = 735.0;  % k21
% allP(6) = gal3C;  % Ptot
% allP(7) = Kd;
% allP(8) = 0.06; % Tc
% I_lig = length(freeligC);
% Is = length(allP)+1;
% allP(Is:Is+I_lig-1) = freeligC;  % Free ligand concentrations

gridrange = [1, 0, 0, 200; ...
    2, 0, 0, 10]; %...
%3 0 10 10000];
%3 1 0 10000;...
%4 0 Pin(4)-1 Pin(4)+1;...
%5 0 Pin(5)-3 Pin(5)+3];

Ndata = 0;
for i = 1:length(ResToFit) % loop over each res to fit
    for j = 1:length(freeligC) % loop over each concentration
        for k = 1:length(carrier_freq{j}) % Loop over each carrier freq for each res
            I = find(str2double(resname{j}(:)) == ResToFit(i));
            Ndata = Ndata + 1; % Data counter
            redFitted(Ndata) = ResToFit(i);
            X{Ndata} = nucp{j}{k}; % Assign X data (nu cpmg)
            Y{Ndata} = R2{j}{k}(I, :); % Assign Y data (R2)
            Ysig{Ndata} = R2sig{j}{k}(I, :); % Assign error in y data

            res_lig_field(Ndata, :) = [ResToFit(i), freeligC(j), carrier_freq{j}(k)];


            GeneralNfit(Ndata, [1:3, 7:9, 11, 14]) = [1:8];
            GeneralNfit(Ndata, 10) = Is + j - 1;

            % Chemical shift
            if j == 1 && k == 1
                apoIx = find(shiftsApo{1} == ResToFit(i));
                lacIx = find(shiftsLac{1} == ResToFit(i));

                allP = [allP, gridBestFitRed(ResToFit(i), 1:2)];
                ppm23Pix = [length(allP) - 1:length(allP)];

                % restrict chemical shift of state 2 and 3
                nlopt.pbounds(size(nlopt.pbounds, 1)+1, :) = [ppm23Pix(1), allP(ppm23Pix(1)) - error_shift2(i), allP(ppm23Pix(1)) + error_shift2(i)];
                nlopt.pbounds(size(nlopt.pbounds, 1)+1, :) = [ppm23Pix(2), double(shiftsLac{3}(lacIx)) - 2, double(shiftsLac{3}(lacIx)) + 2];

                allP = [allP, shiftsApo{3}(apoIx), double(shiftsLac{3}(lacIx))]; %shiftsLac{1}(lacIx)
                ppm14Pix = [length(allP) - 1:length(allP)];
                gridrange(size(gridrange, 1)+1, :) = [ppm23Pix(1), 0, allP(ppm23Pix(1)) - 3, allP(ppm23Pix(1)) + 3];
                gridrange(size(gridrange, 1)+1, :) = [ppm23Pix(2), 0, allP(ppm23Pix(2)) - 3, allP(ppm23Pix(2)) + 3];


                freeze = [freeze, length(allP) - 1:length(allP)];
            end
            GeneralNfit(Ndata, [4, 5]) = ppm23Pix;
            GeneralNfit(Ndata, [12, 13]) = ppm14Pix;

            %Carrier frq
            allP = [allP, carrier_freq{j}(k)];
            freeze = [freeze, length(allP)];
            GeneralNfit(Ndata, 15) = length(allP);


            [tmp, maxix] = max(X{Ndata});


            % R20
            if j == 1
                allP = [allP, Y{Ndata}(maxix)];
                r20Pix = length(allP);
                GeneralNfit(Ndata, 6) = r20Pix; % R20
                R20_freeze = [R20_freeze, length(allP)];
                nlopt.pbounds(size(nlopt.pbounds, 1)+1, :) = [r20Pix, 8, inf];
            elseif j == 2
                allP = [allP, Y{Ndata}(maxix)];
                r20Pix = length(allP);
                GeneralNfit(Ndata, 6) = r20Pix; % R20
                R20_freeze = [R20_freeze, length(allP)];
                nlopt.pbounds(size(nlopt.pbounds, 1)+1, :) = [r20Pix, 8, inf];
            else
                GeneralNfit(Ndata, 6) = r20Pix;
            end


        end
    end
end

for i = 1:size(GeneralNfit, 1)
    paraMap{i} = GeneralNfit(i, :);
end

freeze_grid = sort([freeze, R20_freeze]);
nlopt.paramap = paraMap;

%% Grid search
if gridsearch

    nlopt.paramap = paraMap;

    [Pgrid, chi2min, gridTables] = ...
        nlgridsearch(X, Y, Ysig, @r2cpmg_fourstate_ppm_constrain, allP, freeze_grid, gridrange, 3000, nlopt)

    nlgridinspect(Pgrid, gridrange, gridTables)


else
    Pgrid = allP;
end

nlopt.paramap = paraMap;
nlopt.figno = 0;

% Fitting
for k = 1:length(Ifit)

    % initial parameters
    Pgrid = in.Pfi{Ifit(k)};

    % fit data
    [Pfi{k}, Psig{k}, fitstats{k}, COV{k}, SIM{k}] = nlfit(X, Y, Ysig, @r2cpmg_fourstate_ppm_constrain, Pgrid, freeze, nlopt);

    %% plotting


    if makeplot.make == 'y'
        for i = 1:length(X);
            figure(length(X)+i);
            hold on;
            errorbar(X{i}, Y{i}, Ysig{i}, '.', 'MarkerSize', 13);
            plot(X{i}, r2cpmg_fourstate_ppm_constrain(X{i}, Pfi{k}(GeneralNfit(i, :))));
            tit = sprintf('Res: %3.0f,  L_{free}: %1.4f,  field: %2.2f', res_lig_field(i, 1), res_lig_field(i, 2), res_lig_field(i, 3));
            title(tit);
            xlabel('\nu_{cpmg} (Hz)');
            ylabel('R_2 (s^{-1})');

            if makeplot.print == 'y';
                if ~exist(sprintf('2.0f', k), 'dir');
                    mkdir(sprintf('%2.0f', k));
                end
                file_name = sprintf('%2.0f/Res_%3.0f_L_%1.4f_field_%2.2f.tiff', k, res_lig_field(i, 1), res_lig_field(i, 2), res_lig_field(i, 3));
                print(gcf, '-dtiff', '-r300', file_name);
            end
        end
    end

    %% print data to file
    if makeplot.dat2file == 'y'
        if isfield(makeplot, 'print_directory');
            txt_file = sprintf('%2.0f/X_%3.0f.txt', k, Pgrid(2)*100);
            fileID = fopen(txt_file, 'w+');
            for i = 1:length(X);
                out = [res_lig_field(i, :), Pfi{k}(GeneralNfit(i, :))];
                for j = 1:length(out);
                    fprintf(fileID, '%0.5f   ', out(j));
                end
                fprintf(fileID, '\n');
            end
        end
    end
    close all
end

save('Final')
