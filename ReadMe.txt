


This script fits a four state model accounting for conformational selection and induced fit. The main m-file is the Fit_4_state_global_R20. It takes as input several CPMG data sets at different ligand concentrations. In addition, chemical shifts of the apo state and ligand-saturated protein state are needed.


Pin - [kon  rho_close  rho_off] kon divided by 1000, rho_close divided by 100 and rho_off as described in the article
cpmgfiles - contains the path and file names to the cpmg data
totligC - list of all total ligand concentrations, one for each entry in cpmgfiles
freeligC - free ligand concentration
gal3C - concentration of protein
Kd - dissociation constant
fid - loads the data for chemical shifts (apo/ligand bound) 
allP(4) - exchange rate going from state 1 to state 2
allP(5) - exchange rate going from state 2 to state 1
allP(8) - Constant relaxation time in seconds
ResToFit - residue number of residues to include in fit
gridBestFitRed(residue number) - guesses of chemical shift of states 2 and 3



Matlab 2019b used 