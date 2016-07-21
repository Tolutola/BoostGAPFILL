function [newReactions,RxnsRecovered,Stats,newModel]=BGF(Model,Rxns2Remove,noOfRxns2Add,noOfExpts,method,mFlag,blacklist)

% preprocess data
tempModel=removeRxns(Model,Rxns2Remove,[],false);
S = Model.S;    %stoichiometric matrix S
S = full(spones(S));   %make to binary
R = Model.rxns;   %reactions R
M = Model.mets;   %metabolites M
US = full(spones(Model.US)); %US with logical values
sUS = Model.US;   %US with stoichiometry values    unrnames = Model.unrnames;
unrnames = Model.unrnames;


%% experiment begins
Testnumber = length(Rxns2Remove);
num_test = length(Testnumber);
RxnsRecovered = zeros(num_test,noOfExpts);


for change_testnumber = 1:num_test
    testnumber = Testnumber(change_testnumber);
    threshold=noOfRxns2Add;
    %     threshold=length(unrnames);
    %     poolobj = parpool(feature('numcores'));
    for ith_experiment = 1:noOfExpts
        
        trainnames=tempModel.rxns;
        trainr = sparse(S(:,ismember(R,trainnames)));    %train reactions
        
        tmp = strcmp(repmat(trainnames,1,size(unrnames,2)),repmat(unrnames,size(trainnames,1),1));  %compare train with unr, save to a matrix
        check_repeat = find(sum(tmp)==0);     %keep those unr which never appears in the train reactions
        unrepeat_US = US(:,check_repeat);     %the unrepeated US, which excludes the train reactions
        unrepeat_unrnames = unrnames(1,check_repeat);
        
        %black list
        if nargin ~= 7
            [~,cc] = size(unrepeat_US);
            blacklist = zeros(cc,1);
        end
        % METHOD
        match_all = strcmp(repmat(Rxns2Remove,1,size(unrepeat_unrnames,2)),repmat(unrepeat_unrnames,size(Rxns2Remove,1),1));
        switch method
            case 1
                [Lambda,scores] = HLpredict(trainr,unrepeat_US,threshold,'Greedy',ith_experiment,blacklist);
            case 2
                [Lambda,scores] = HLpredict(trainr,unrepeat_US,threshold,'CM',ith_experiment,blacklist);
            case 3
                [Lambda,scores] = HLpredict(trainr,unrepeat_US,threshold,'FM',ith_experiment,blacklist);
            case 4
                [Lambda,scores] = HLpredict(trainr,unrepeat_US,threshold,'MATBoost',ith_experiment,10,blacklist);
            case 5
                [Lambda,scores] = HLpredict(trainr,unrepeat_US,threshold,'HCN',ith_experiment,blacklist);
            case 6
                [Lambda,scores] = HLpredict(trainr,unrepeat_US,threshold,'HKatz',ith_experiment,blacklist);
            case 7
                [Lambda,scores] = HLpredict(trainr,unrepeat_US,threshold,'Submodular',ith_experiment,blacklist);
            case 8
                [Lambda,scores] = HLpredict(trainr,unrepeat_US,threshold,'HPLSF',ith_experiment,blacklist);
            case 9
                [Lambda,scores] = HLpredict(trainr,unrepeat_US,threshold,'SPHC',ith_experiment,blacklist);
        end
        
        [~,I] = sort(scores,1,'descend');
        
        %calculate weights for reactions
        %method 1 (very bad)
%         tempWeights=logspace(log10(0.001),log10(1000),length(scores));
%         WeightsPerRxn.var=tempWeights(I);
        %method 2 (fair)
        %         weights=scores;
        %         weights(weights<1e-3)=1e-3;
        %         weights=(0.1./weights).^3;
        %         WeightsPerRxn.var=weights;
        
        %         %method 3 (bad)
        %         WeightsPerRxn.var=1000*ones(length(scores),1);
        %         WeightsPerRxn.var(I(1:noOfRxns2Add))=500;
        
%         method 4 
%         WeightsPerRxn.var=1000*ones(length(scores),1);
%         WeightsPerRxn.var(I(1:noOfRxns2Add))=1e-7;
        
        %method 5
        %use a threshold like 0.9
         WeightsPerRxn.var=1000*ones(length(scores),1);
        WeightsPerRxn.var(scores>=0.9)=10;
                
        WeightsPerRxn.rxns=unrepeat_unrnames';
        
        [newReactions,RxnsRecovered,Stats,newModel]=FastGapFillSub(Model,Rxns2Remove,WeightsPerRxn);
        
        
    end
    
    
end