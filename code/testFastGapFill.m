function [newReactions,RxnsRecovered,Stats,newModel]=testFastGapFill(model,Rxns2Remove)
% adaptation of the fastGAPFILL algorithm (to facilitate comparison)
%change the nomenclature of mets
for i =1:length(model.mets)
    temp=model.mets{i};
    model.mets{i}=strcat(temp(1:end-2),'[',temp(end),']');
end

oldRxns=model.rxns;
model.OldMets=model.mets;
model = removeRxns(model,Rxns2Remove,[],false);
% tempModel=model;

% define weights for reactions to be added - the lower the weight the
% higher the priority
weights.MetabolicRxns = 0.1; % Kegg metabolic reactions
weights.ExchangeRxns = 500; % Exchange reactions
weights.TransportRxns = 500; % Transport reactions

% performance of algorithm is best if the weighting parameter is not 1
%
% all reactions are equally weighted
% weights.MetabolicRxns = 1000; % Kegg metabolic reactions
% weights.ExchangeRxns = 1000; % Exchange reactions
% weights.TransportRxns = 1000; % Transport reactions

%% Do not change below here
% Prepare the output table with statistics
cnt = 1;
Stats{cnt,1} = 'Model name';cnt = cnt+1;
Stats{cnt,1} = 'Size S (original model)';cnt = cnt+1;
Stats{cnt,1} = 'Number of compartments';cnt = cnt+1;
Stats{cnt,1} = 'List of compartments';cnt = cnt+1;
Stats{cnt,1} = 'Number of blocked reactions';cnt = cnt+1;
Stats{cnt,1} = 'Number of solvable blocked reactions';cnt = cnt+1;
Stats{cnt,1} = 'Size S (flux consistent)';cnt = cnt+1;
Stats{cnt,1} = 'Size SUX (including solvable blocked reactions)';cnt = cnt+1;
Stats{cnt,1} = 'Number of added reactions (all)';cnt = cnt+1;
Stats{cnt,1} = 'Number of added metabolic reactions ';cnt = cnt+1;
Stats{cnt,1} = 'Number of added transport reactions ';cnt = cnt+1;
Stats{cnt,1} = 'Number of added exchange reactions ';cnt = cnt+1;
Stats{cnt,1} = 'Time preprocessing';cnt = cnt+1;
Stats{cnt,1} = 'Time fastGapFill';


i=1;
cnt = 1;
% remove constraints from exchange reactions
EX = strmatch('EX_',model.rxns);
model.lb(EX)=-100;
model.ub(EX)=100;
clear EX
% get stats
Stats{cnt,i+1} = 'model';cnt = cnt+1;
[a,b] = size(model.S);
Stats{cnt,i+1} = strcat(num2str(a),'x',num2str(b));cnt = cnt+1;
% Number of compartments
[tok,rem] = strtok(model.mets,'\[');
rem = unique(rem);
Stats{cnt,i+1} = num2str(length(rem));cnt = cnt+1;
Rem = rem{1};
for j = 2:length(rem)
    Rem = strcat(Rem,',',rem{j});
end
Stats{cnt,i+1} = Rem;cnt = cnt+1;
clear Rem tok rem;



%Prepare fastGapFill
[consistModel,consistMatricesSUX,BlockedRxns] = prepareFastGapFill2(model);

Stats{cnt,i+1} = num2str(length(BlockedRxns.allRxns));cnt = cnt+1;
Stats{cnt,i+1} = num2str(length(BlockedRxns.solvableRxns));cnt = cnt+1;
[a,b] = size(consistModel.S);
Stats{cnt,i+1} = strcat(num2str(a),'x',num2str(b));cnt = cnt+1;
[a,b] = size(consistMatricesSUX.S);
Stats{cnt,i+1} = strcat(num2str(a),'x',num2str(b));cnt = cnt+1;

% fastGapsFill
epsilon = 1e-4;
[AddedRxns,newModel] = fastGapFill2(consistMatricesSUX,epsilon,weights);

Stats{cnt,i+1} = num2str(length(AddedRxns.rxns));


% Postprocessing
[AddedRxnsExtended] = postProcessGapFillSolutions(AddedRxns,model,BlockedRxns,0);

Stats{cnt,i+1} = num2str(AddedRxnsExtended.Stats.metabolicSol);cnt = cnt+1;
Stats{cnt,i+1} = num2str(AddedRxnsExtended.Stats.transportSol);cnt = cnt+1;
Stats{cnt,i+1} = num2str(AddedRxnsExtended.Stats.exchangeSol);cnt = cnt+1;

newReactions=setdiff(AddedRxns.rxns,model.rxns);
RxnsRecovered=newReactions(ismember(newReactions,oldRxns));



