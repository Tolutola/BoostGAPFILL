% demo run of BoostGAPFILL
% requires COBRA toolbox

clear
clc
addpath('code')
global def_weight
% load data/iJO1366
% baseModel=Model; clear Model

% iAF1260 E. coli model
% universal matrix is stored in the struct
load data/iAF1260b

%BoostGAPFILL
%BoostGapFill paramters
%options structure
% mFlag:          boolean indicating if a new model with the predicted reactions is to be returned
% blacklist:      list of blacklisted reactions
% maxIter:        maximum number of iterations
% solver:         solver to be used
% mode:           different modes of running BoostGAPFILL ('pattern','pattern-constraints','integrated')
% newMet:         flag to determine whether or not to include reactions with new
% metabolites
% integerSoln:   flag to determine whether to solve the orignal integer
% least squares problem or the relaxed version

options.maxIter=10; %maximum number of interations to run
options.numAlternatives=1; % number of alternative solutions to generate
% solvers
% the following solvers can be used: 'ibm_cplex','gurobi',and MATLAB's 'lsqlin'.
options.solver='ibm_cplex'; %

options.mFlag=1; % return extended Model with new reactions
options.newMet=false;
options.integerSoln=false;
options.numSolns=1; % number of solution sets wanted
Blacklist=[]; % no blacklisted reaction
options.biomass_threshold=0.05; % threshold for viable biomass growth
options.def_weight=10; % used when integrating with FASTGAPFILL
options.maxTime=2*60*60; % used to limit the time for each iteration of the optimization problem in BoostGapFill.
% important to set when running in integrated mode and when integer solutions are required
options.rxnThreshold=1e-12; % value between 0 and 1 ; used in BoostGAPFILL to alter the number of reactions selected; lower values in larger number of reactions selected.
options.newMetPenalty=3; % paramter that penalizes the addition of new metabolites. higher numbers mean fewer selected reactions with new metabolites
def_weight=options.def_weight;


Rxns2Remove=[]; % no reaction is removed from the model. This variable was useful during testing
noOfRxns2Add=10; % the number of reaction predictions can be directly set

%     % BoostGAPFILL (mode 1)
%
options.mode='pattern'; %run BoostGapFill in the pattern mode
tic
[newReactionsB,RxnsRecoveredB,newMetabolites,newModelsB]= ...
    BoostGapFill(Model,Rxns2Remove,noOfRxns2Add,Blacklist,options);
TimesB=toc;

%     % BoostGAPFILL (mode 2)

options.mode='pattern-constraints';% run BoostGapFill as preprocessing step for FASTGAPFILL
tic
[newReactionsB2,RxnsRecoveredB2,newMetabolites2,newModelsB2]= ...
    BoostGapFill(Model,Rxns2Remove,noOfRxns2Add,Blacklist,options);
TimesB2=toc;


% BoostGAPFILL (mode 3)
options.mode='integrated';%
tic
[newReactionsB3,RxnsRecoveredB3,newMetabolites3,newModelsB3]= ...
    BoostGapFill(Model,Rxns2Remove,noOfRxns2Add,Blacklist,options);
TimesB3=toc;


