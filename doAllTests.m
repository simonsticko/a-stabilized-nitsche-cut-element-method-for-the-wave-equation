clear;
addpath('performanceTests/');
addpath('assembling');
addpath('analyticFunctions/plane/')
addpath('analyticFunctions/graffWave/')
addpath('analyticFunctions/besselWave/')
addpath('errorCalculation/')
outerProblConvergence();
membraneConvergence();
foldername=sliverTimeStep();
processtimestep(foldername);
domainDepedenceGammaM();
outerProblConvergence();
dispersionError();