% demo run of BoostGAPFILL
% requires COBRA toolbox
% Comparing BoostGAPFILL and FastGAPFILL

clear
clc
addpath('code')
% changeCobraSolver('ibm_cplex','LP');
% changeCobraSolver('ibm_cplex','MILP');



% iAF1260 E. coli model
% universal matrix is stored in the struct
load data/iAF1260b


%BoostGapFill paramters
maxIter=10; %maximum number of interations to run
solver='ibm_cplex'; %  the MATLAB solver 'lsqlin' can also be used.
mode='pattern';  % run BoostGapFill in the pattern mode
mFlag=1; % return extended Model with new reactions
Blacklist=[]; % no blacklisted reaction

% experimental design
nTest=5; % five time points
% number of reactions to delete
TestVector=floor(linspace(10,1500,nTest)); 
Rxns2Remove=cell(nTest,1);
OldRootGaps=cell(nTest,1);% gaps in old model (after deleting reactions)
numOldGaps=zeros(nTest,1); % number of gaps in old model (after deleting reactions)


%BoostGAPFILL
newReactionsB=cell(nTest,1); %new reaction predictions
RxnsRecoveredB=cell(nTest,1); % new reactions that were part of the deleted set
NewRootGapsB=cell(nTest,1);

numNewGapsB=zeros(nTest,1); % number of gaps in new model (after adding predicted reactions)
newModelsB=cell(nTest,1); % new model with predicted reactions

%FastGAPFILL
newReactionsF=cell(nTest,1); %new reaction predictions
RxnsRecoveredF=cell(nTest,1); % new reactions that were part of the deleted set
NewRootGapsF=cell(nTest,1);
OldRootGapsF=cell(nTest,1);
numNewGapsF=zeros(nTest,1); % number of gaps in new model (after adding predicted reactions)
newModelsF=cell(nTest,1); % new model with predicted reactions
Stats=cell(nTest,1); % stats struct (only for FastGAPFILL)

% experiment begins
for j=1:nTest
    
    rng('shuffle', 'twister')
    noOfRxns2Add=TestVector(j); %number of reactions to remove
    
    %randomly select reactions
    
    perm = randperm(size(Model.S,2));    %randomly split S and R into test and train
    Rxns2Remove{j} = Model.rxns(perm(1:noOfRxns2Add));
    
    %calculate number of gaps in reduced model
    tempModel=removeRxns(Model,Rxns2Remove{j},[],false);
    [~,OldRootGaps{j}] = gapFind2(tempModel);
    numOldGaps(j)=length(OldRootGaps{j});
    
    % BoostGAPFILL    
    [newReactionsB{j,1},RxnsRecoveredB{j,1},newModelsB{j}]= ...
    BoostGapFill(Model,Rxns2Remove{j},noOfRxns2Add,mFlag,Blacklist,maxIter,solver,mode);
    % calculate new number of gaps in extended model
    [~,NewRootGapsB{j}] = gapFind2(newModelsB{j});
    numNewGapsB(j)=length(NewRootGapsB{j});
    
    % FastGAPFILL
    [newReactionsF{j,1},RxnsRecoveredF{j,1},Stats{j,1},newModelsF{j}]=testFastGapFill(Model,Rxns2Remove{j});
    % calculate new number of gaps in extended model
    [~,NewRootGapsF{j}] = gapFind2(newModelsF{j});
    numNewGapsF(j)=length(NewRootGapsF{j});
 
   disp (j)
end





