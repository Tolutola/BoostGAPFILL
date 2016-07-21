function [newReactions,RxnsRecovered,tempModel]=BoostGapFill(Model,Rxns2Remove,noOfRxns2Add,mFlag,blacklist,maxIter,solver,mode)

% main BoostGapFill function
%input
% model:          COBRA model
% Rxns2Remove:    reactions to be deleted
% noOfRxns2Add:   number of reactions to  be added
% mFlag:          boolean indicating if a new model with the predicted reactions is to be returned
% blacklist:      list of blacklisted reactions
% maxIter:        maximum number of iterations
% solver:         solver to be used
% mode:           different modes of running BoostGAPFILL ('pattern','pattern-constraints','integrated')

%output
% newReactions : new reaction predictions
% RxnsRecovered: reactions that match deleted reactions
% tempModel:     new model created after adding predicted reactions

% Authors: Tolutola Oyetunde and Muhan Zhang,  Washington University in St. Louis


% preprocess data
tempModel=removeRxns(Model,Rxns2Remove,[],false);
S = Model.S;    %stoichiometric matrix S
S = full(spones(S));   %make to binary
R = Model.rxns;   %reactions R
US = full(spones(Model.US)); %US with logical values
sUS = Model.US;
unrnames = Model.unrnames;

Testnumber = length(Rxns2Remove);
num_test = length(Testnumber);
RxnsRecovered = zeros(num_test,1);

threshold=noOfRxns2Add;

trainnames=tempModel.rxns;
trainr = sparse(S(:,ismember(R,trainnames)));    %train reactions

tmp = strcmp(repmat(trainnames,1,size(unrnames,2)),repmat(unrnames,size(trainnames,1),1));  %compare train with unr, save to a matrix
check_repeat = find(sum(tmp)==0);     %keep those unr which never appears in the train reactions
unrepeat_US = US(:,check_repeat);     %the unrepeated US, which excludes the train reactions
unrepeat_unrnames = unrnames(1,check_repeat);

%black list
if isempty(blacklist)
    [~,cc] = size(unrepeat_US);
    blacklist = zeros(cc,1);
end


[Lambda,scores] = HLpredict(trainr,unrepeat_US,threshold,maxIter,blacklist,solver);

switch mode
    case 'pattern'

    newReactions = unrepeat_unrnames(Lambda');
    match = strcmp(repmat(Rxns2Remove,1,size(newReactions,2)),repmat(newReactions,size(Rxns2Remove,1),1));
    RxnsRecovered= nnz(match);
  
    if mFlag % if new extendedModel is required

        sUS = sUS(:,check_repeat);
        dS=sUS(:,logical(Lambda'));

        % create new Model with added reactions
        for nrxns=1:noOfRxns2Add
            cVector=dS(:,nrxns);
            cLogic=cVector~=0;
            tempModel = addReaction(tempModel,newReactions{nrxns},Model.mets(cLogic),cVector(cLogic));
        end

    end
    
    % case 'pattern-constraints'
    
    % case 'integrated'
end


end