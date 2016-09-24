% Comparing GapFilling algorithms
clear
clc

changeCobraSolver('ibm_cplex','LP');
changeCobraSolver('ibm_cplex','MILP');

load data/iJO1366
baseModel=Model; clear Model
load Kdictionary

%model to use
% core model
% load data/e_coli_core


% iAF1260b model
load data/iAF1260b
% Model=Ecoli

% % iJO1366 model
% load data/iJO1366


% reactions to remove
nTest=5;
TestVector=linspace(10,1500,nTest);
RVector=linspace(300,1000,nTest);
RFrac=zeros(nTest,1);
newReactions=cell(nTest,1);
RxnsRecovered=cell(nTest,1);
goodRxns=cell(nTest,1);
Rxns2Remove=cell(nTest,1);
NewRootGaps=cell(nTest,1);
OldRootGaps=cell(nTest,1);
Stats=cell(nTest,1);
% noOfRxns2Add=500;

% [allGaps1,rootGaps1,downstreamGaps1] = gapFind(Model);
% noOldGaps=length(rootGaps1);
noNewGaps=zeros(nTest,1);
noOldGaps=zeros(nTest,1);
rootGapsMatrix=cell(nTest,1);
% Rxns2Remove={'12DGR120tipp'}; %remove a trivial reaction
newModels=cell(nTest,1);

%biomass reaction id
bioID=1005;

for j=1:nTest
    %     i=1;
    % for i=1:nTest
    rng('shuffle', 'twister')
    noOfRxns2Add=floor(TestVector(j)); %number of reactions to remove
    
    %randomly select reactions
    
    perm = randperm(size(Model.S,2));    %randomly split S and R into test and train
    Rxns2Remove{j} = Model.rxns(perm(1:noOfRxns2Add));
    
    %ensure that the biomass reaction is not removed
    keepID=find(strcmp(Rxns2Remove{j},Model.rxns(bioID)));
    if ~isempty(keepID)
        Rxns2Remove{j}(keepID)=Model.rxns(perm(noOfRxns2Add+1));
    end
        
    
    %calculate number of gaps in reduced model
    tempModel=removeRxns(Model,Rxns2Remove{j},[],false);
    [~,OldRootGaps{j}] = gapFind2(tempModel);
    noOldGaps(j)=length(OldRootGaps{j});
    
    
    %     noOfRxns2Add=length(Rxns2Remove);
    
    % %     [newReactions{j,1},RxnsRecovered{j,1},newModels{j}]=GapFill(Model,Rxns2Remove{j},noOfRxns2Add,noOfExpts,method,1);
    %     %         [newReactions{j,1},RxnsRecovered{j,1},Stats{j,1},newModels{j}]=testFastGapFill(Model,Rxns2Remove{j});
    %     [newReactions{j,1},RxnsRecovered{j,1},newModels{j}]=BGF(Model,Rxns2Remove{j},noOfRxns2Add,noOfExpts,method,1);
    iterations=25;
    threshold=0.05;
    [solution]=growthExpMatch2(tempModel, 'KEGGFilename.lst', '[c]', iterations, Kdictionary,'GEMLog',threshold,[],Rxns2Remove{j});
    
    [nr,nc]=size(solution.importedRxns);
    tempR=cell(nr,nc);
    
    for p=1:nr
        
        for k=1:nc
            if ~isempty(solution.importedRxns{p,k})
                temp=cellstr(solution.importedRxns{p,k});
                temp=temp{:};
                tempR{p,k}=temp(1:end-2);
                
            end
        end
        nl=sum(cell2mat(cellfun(@(x) ~isempty(x),tempR(p,:),'UniformOutput',false)));
        newReactions{j,1}=union(newReactions{j,1},tempR(p,1:nl));
        
    end
    if length(newReactions{j,1})> noOfRxns2Add
        tempx=newReactions{j,1};
        newReactions{j,1}=tempx(1:noOfRxns2Add);
    end
    goodRxns{j}=intersect(newReactions{j,1},Model.rxns);
    
     
    % calculate new number of gaps in extended model
    %     [~,NewRootGaps{j}] = gapFind2(newModels{j});
    %     noNewGaps(j)=length(NewRootGaps{j});
    
    %         RFrac(j)=length(RxnsRecovered{j})/length(newReactions{j});
    disp (j)
end
% end
save SMILEYFill3 noNewGaps goodRxns newReactions RxnsRecovered noOldGaps newModels baseModel NewRootGaps OldRootGaps Rxns2Remove
% save FastGapFill2 noNewGaps goodRxns newReactions RxnsRecovered noOldGaps newModels NewRootGaps OldRootGaps Rxns2Remove Stats
% save FastGapFill1 newReactions goodRxns RxnsRecovered Stats noOfRxns2Add
% save BoostGapFill2b noNewGaps goodRxns newReactions RxnsRecovered noOldGaps newModels baseModel NewRootGaps OldRootGaps Rxns2Remove
% save BoostGapFill1 noNewGaps goodRxns newReactions RxnsRecovered noOldGaps newModels baseModel
% save FastGapFill newReactions RxnsRecovered Rxns2Remove RFrac

%% SMILEY algorithm
% set these reactions to zero (ARAI, RBK_L1, RMPA, LYXI, RMI, RMK, and
% LACZ)
% they are not present in K12 MG1655

% essential genes (run pFBA)
% core biomass objective function used
% reactions constrained to zero to avoid unrealistic behavior
% CAT,DHPTDNR, DHPTDNRN, FHL (formate hydrogen lyase), SPODM,
% SPODMpp, SUCASPtpp, SUCFUMtpp, SUCMALtpp, and SUCTARTtpp

% [solution]=growthExpMatch(Model, KEGGFilename, compartment, iterations, dictionary, logFile, threshold,KEGGBlackList);


% comparison of OMNI to different levels of gene deletion

