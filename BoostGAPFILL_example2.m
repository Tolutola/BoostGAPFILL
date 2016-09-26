% demo run of BoostGAPFILL
% requires COBRA toolbox
% Comparing BoostGAPFILL and FastGAPFILL

clear
clc
addpath('code')
global def_weight


% iAF1260 E. coli model
% universal matrix is stored in the struct
load data/iAF1260b

% experimental design
nTest=10; %
% number of reactions to delete
TestVector=floor(linspace(10,1500,nTest));
Rxns2Remove=cell(nTest,1);
OldRootGaps=cell(nTest,1);% gaps in old model (after deleting reactions)
numOldGaps=zeros(nTest,1); % number of gaps in old model (after deleting reactions)


%BoostGAPFILL
%BoostGapFill paramters
%options structure
% mFlag:                boolean indicating if a new model with the predicted reactions is to be returned
% maxIter:              maximum number of iterations
% solver:               solver to be used% the following solvers can be used: 'ibm_cplex','gurobi' and MATLAB's 'lsqlin'.
% mode:                 different modes of running BoostGAPFILL ('pattern','pattern-constraints','integrated')
% newMet:               flag to determine whether or not to include reactions with new metabolites
% integerSoln:          flag to determine whether to solve the orignal integer
%                       least squares problem or the relaxed version
% numSolns:             number of solution sets required
% biomass_threshold:    threshold for viable biomass growth
% maxTime:              used to limit the time for each iteration of the optimization problem in BoostGapFill.
%                       important to set when running in integrated mode and when integer solutions are required
%  def_weight:          used when integrating with FASTGAPFILL
% rxnThreshold:         value between 0 and 1 ; used in BoostGAPFILL to alter the number of reactions selected;
%                       lower values in larger number of reactions selected.        
% newMetPenalty:        newMetPenalty paramter that penalizes the addition of new metabolites. 
%                       higher numbers mean fewer selected reactions with new metabolites       

options.maxIter=10; %maximum number of interations to run
options.numAlternatives=1; % 
options.solver='ibm_cplex'; %

options.mFlag=1; % return extended Model with new reactions
options.newMet=false;
options.integerSoln=false;
options.numSolns=1; 
Blacklist=[]; % no blacklisted reaction
options.biomass_threshold=0.05; 
options.def_weight=10; %
options.maxTime=2*60*60; 
options.rxnThreshold=1e-12; 
options.newMetPenalty=3; % 
def_weight=options.def_weight;

newReactionsB=cell(nTest,1); %new reaction predictions
RxnsRecoveredB=cell(nTest,1); % new reactions that were part of the deleted set
NewRootGapsB=cell(nTest,1);
newMetabolites=cell(nTest,1);
numNewGapsB=zeros(nTest,1); % number of gaps in new model (after adding predicted reactions)
newModelsB=cell(nTest,1); % new model with predicted reactions

newReactionsB2=cell(nTest,1); %new reaction predictions
RxnsRecoveredB2=cell(nTest,1); % new reactions that were part of the deleted set
NewRootGapsB2=cell(nTest,1);
newMetabolites2=cell(nTest,1);
numNewGapsB2=zeros(nTest,1); % number of gaps in new model (after adding predicted reactions)
newModelsB2=cell(nTest,1); % new model with predicted reactions

newReactionsB3=cell(nTest,1); %new reaction predictions
RxnsRecoveredB3=cell(nTest,1); % new reactions that were part of the deleted set
NewRootGapsB3=cell(nTest,1);
newMetabolites3=cell(nTest,1);
numNewGapsB3=zeros(nTest,1); % number of gaps in new model (after adding predicted reactions)
newModelsB3=cell(nTest,1); % new model with predicted reactions

%FastGAPFILL
newReactionsF=cell(nTest,1); %new reaction predictions
RxnsRecoveredF=cell(nTest,1); % new reactions that were part of the deleted set
NewRootGapsF=cell(nTest,1);
numNewGapsF=zeros(nTest,1); % number of gaps in new model (after adding predicted reactions)
newModelsF=cell(nTest,1); % new model with predicted reactions
Stats=cell(nTest,1); % stats struct (only for FastGAPFILL)

%times
TimesB=cell(nTest,1);
TimesB2=cell(nTest,1);
TimesB3=cell(nTest,1);
TimesF=cell(nTest,1);


%gaps in original iAF1260
[~,origGaps] = gapFind2(Model);


% experiment begins
for j=1:nTest
    
    
    noOfRxns2Add=TestVector(j); %number of reactions to remove
    %
    %     %randomly select reactions
    perm = randperm(size(Model.S,2));    %randomly split S and R into test and train
    Rxns2Remove{j} = Model.rxns(perm(1:noOfRxns2Add));
    %
    %     %ensure that the biomass reaction is not removed
    keepID=find(strcmp(Rxns2Remove{j},Model.rxns(find(Model.c~=0))));
    if ~isempty(keepID)
        Rxns2Remove{j}(keepID)=Model.rxns(perm(noOfRxns2Add+1));
    end
    %
    %     %calculate number of gaps in reduced model
    tempModel=removeRxns(Model,Rxns2Remove{j},[],false);
    [~,OldRootGaps{j}] = gapFind2(tempModel);
    numOldGaps(j)=length(OldRootGaps{j});
    
    % BoostGAPFILL (mode 1)
    
    options.mode='pattern'; %run BoostGapFill in the pattern mode
    tic
    [newReactionsB{j,1},RxnsRecoveredB{j,1},newMetabolites{j},newModelsB{j}]= ...
        BoostGapFill(Model,Rxns2Remove{j},noOfRxns2Add,Blacklist,options);
    TimesB{j}=toc;
    % calculate new number of gaps in extended model
    [~,NewRootGapsB{j}] = gapFind2(newModelsB{j}{:});
    numNewGapsB(j)=length(NewRootGapsB{j});
    %
    %     %
    % BoostGAPFILL (mode 2)
    
    options.mode='pattern-constraints';% run BoostGapFill as preprocessing step for FASTGAPFILL
    %         options.rxnThreshold=testThreshold(j);
    tic
    [newReactionsB2{j,1},RxnsRecoveredB2{j,1},newMetabolites2{j},newModelsB2{j}]= ...
        BoostGapFill(Model,Rxns2Remove{j},noOfRxns2Add,Blacklist,options);
    TimesB2{j}=toc;
    % calculate new number of gaps in extended model
    [~,NewRootGapsB2{j}] = gapFind2(newModelsB2{j}{:});
    numNewGapsB2(j)=length(NewRootGapsB2{j});
    
    % BoostGAPFILL (mode 3)
    options.mode='integrated';%
    tic
    [newReactionsB3{j,1},RxnsRecoveredB3{j,1},newMetabolites3{j},newModelsB3{j}]= ...
        BoostGapFill(Model,Rxns2Remove{j},noOfRxns2Add,Blacklist,options);
    TimesB3{j}=toc;
    % calculate new number of gaps in extended model
    [~,NewRootGapsB3{j}] = gapFind2(newModelsB3{j}{:});
    numNewGapsB3(j)=length(NewRootGapsB3{j});
    
    % FASTGAPFILL
    tic
    [newReactionsF{j,1},RxnsRecoveredF{j,1},Stats{j,1},newModelsF{j}]=testFastGapFill(Model,Rxns2Remove{j},options.newMet);
    TimesF{j}=toc;
    % calculate new number of gaps in extended model
    [~,NewRootGapsF{j}] = gapFind2(newModelsF{j}{:});
    numNewGapsF(j)=length(NewRootGapsF{j});
    
    
    disp (j)
end






