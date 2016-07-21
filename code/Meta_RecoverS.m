clear all;

%% load a metabolic model
mo = 1;
for mo = 3:3
    mo
switch mo
    case 1
        model = 'iJO1366';
        Testnumber = 50:50:400;
    case 2
        model = 'iAF1260b';
        Testnumber = 50:50:400;
    case 3
        model = 'e_coli_core';
        Testnumber = 5:5:40;
    case 4
        model = 'RECON1';
        Testnumber = 100:100:800;
    case 5
        model = 'iAT_PLT_636';
        Testnumber = 50:50:400;
    case 6
        model = 'iAB_RBC_283';
        Testnumber = 25:25:200;
    case 7
        model = 'iAF692';
        Testnumber = 25:25:200;
    case 8
        model = 'iHN637';
        Testnumber = 25:25:200;
    case 9
        model = 'iIT341';
        Testnumber = 25:25:200;
end
load(sprintf('data/%s.mat',model),'Model');
S = Model.S;    %stoichiometric matrix S
S = full(spones(S));   %make to binary
R = Model.rxns;   %reactions R
M = Model.mets;   %metabolites M
US = full(spones(Model.US)); %US with logical values
sUS = Model.US;   %US with stoichiometry values
unrnames = Model.unrnames;

%% experiment begins
numOfExperiment = 1;        %independent experiment numbers
%method = 1;   %1: Greedy 2: CM 3: FM 4: MATBoost 5:common neighbor 6:Katz 7:Submodular 8:HPLSF  9:SPHC
%Method = [3,4,5,6,8,9]
Method = [2]
for md = 1:length(Method)
method = Method(md)

%Testnumber = [800]  %number of test missing reactions

num_test = length(Testnumber);

Num_of_matched_reactions = zeros(num_test,numOfExperiment);
Num_of_selected_reactions = zeros(num_test,numOfExperiment);
Average_guess_match_num = zeros(num_test,numOfExperiment);
AUC_of_Reaction_Prediction = zeros(num_test,numOfExperiment);
for change_testnumber = 1:num_test
    testnumber = Testnumber(change_testnumber)
    threshold = testnumber;             % select the same number missing hyperlinks from test
    rseed = 2732*(testnumber^2-1);
    %poolobj = parpool(feature('numcores'));
    for ith_experiment = 1:numOfExperiment
        ith_experiment
        
        rand('seed',rseed*(ith_experiment^2-1));   %to reproduce multiple exp results
        
        perm = randperm(size(S,2));    %randomly split S and R into test and train
        testr = sparse(S(:,perm(1:testnumber)));     %test(missing) reactions (as columns in matrix)
        trainr = sparse(S(:,perm(testnumber+1:end)));    %train reactions
        testnames = R(perm(1:testnumber));
        trainnames = R(perm(testnumber+1:end));
        
        tmp = strcmp(repmat(trainnames,1,size(unrnames,2)),repmat(unrnames,size(trainnames,1),1));  %compare train with unr, save to a matrix
        check_repeat = find(sum(tmp)==0);     %keep those unr which never appears in the train reactions
        unrepeat_US = US(:,check_repeat);     %the unrepeated US, which excludes the train reactions
        
        unrepeat_unrnames = unrnames(1,check_repeat);
        
        % METHOD
        match_all = strcmp(repmat(testnames,1,size(unrepeat_unrnames,2)),repmat(unrepeat_unrnames,size(testnames,1),1));
        labels = sum(match_all)>0;
        switch method
            case 1
                [Lambda,scores] = HLpredict(trainr,unrepeat_US,threshold,'Greedy',ith_experiment);
            case 2
                [Lambda,scores] = HLpredict(trainr,unrepeat_US,threshold,'CM',ith_experiment);
            case 3
                [Lambda,scores] = HLpredict(trainr,unrepeat_US,threshold,'FM',ith_experiment);
            case 4
                [Lambda,scores] = HLpredict(trainr,unrepeat_US,threshold,'MATBoost',ith_experiment,10);
            case 5
                [Lambda,scores] = HLpredict(trainr,unrepeat_US,threshold,'HCN',ith_experiment);
            case 6
                [Lambda,scores] = HLpredict(trainr,unrepeat_US,threshold,'HKatz',ith_experiment);
            case 7
                [Lambda,scores] = HLpredict(trainr,unrepeat_US,threshold,'Submodular',ith_experiment);
            case 8
                [Lambda,scores] = HLpredict(trainr,unrepeat_US,threshold,'HPLSF',ith_experiment);
            case 9
                [Lambda,scores] = HLpredict(trainr,unrepeat_US,threshold,'SPHC',ith_experiment);
        end
        MUnames = unrepeat_unrnames(Lambda');
        match = strcmp(repmat(testnames,1,size(MUnames,2)),repmat(MUnames,size(testnames,1),1));
        num_of_matched_reactions = nnz(match);
        Num_of_matched_reactions(change_testnumber,ith_experiment) = num_of_matched_reactions;
        num_of_selected_reactions = nnz(Lambda);
        Num_of_selected_reactions(change_testnumber,ith_experiment) = num_of_selected_reactions;
        %average_guess_match_num = nnz(Lambda)*testnumber/size(Lambda,1);
        average_guess_match_num = testnumber*testnumber/size(Lambda,1);
        Average_guess_match_num(change_testnumber,ith_experiment) = average_guess_match_num;
        
        %matched_reactions = MUnames(logical(sum(match)));    %output the matched reactions
        sUS = sUS(:,check_repeat);
        sS = Model.S;
        predS = full([sS(:,perm(testnumber+1:end)),sUS(:,logical(Lambda'))]);    %train reactions
        %predS = [trainr,unrepeat_US(:,logical(Lambda'))];    %the recovered S matrix (orginal S + dS)
        
        [~,~,~,auc] = perfcurve(labels,scores,true);
        if isnan(auc)
            auc=0.5;
        end
        AUC_of_Reaction_Prediction(change_testnumber,ith_experiment) = auc;
        
    end
    %delete(poolobj)
    
end


Recall = bsxfun(@times,Num_of_matched_reactions,1./Testnumber');
Precision = Num_of_matched_reactions./Num_of_selected_reactions;
average_recall = mean(Recall,2);
std_recall = std(Recall,0,2);
average_precision = mean(Precision,2);
std_precision = std(Precision,0,2);
average_guess_match_num = mean(Average_guess_match_num,2)
std_guess_match_num = std(Average_guess_match_num,0,2);
average_AUC = mean(AUC_of_Reaction_Prediction,2)
std_AUC = std(AUC_of_Reaction_Prediction,0,2);
average_match_num = mean(Num_of_matched_reactions,2)
std_match_num = std(Num_of_matched_reactions,0,2);
sound(sin(2*pi*25*(1:4000)/100));


save(sprintf('data/result/%s_%d.mat',model,method),'average_match_num', 'std_match_num', 'average_guess_match_num', 'std_guess_match_num', 'average_AUC', 'std_AUC','average_recall','std_recall','average_precision','std_precision','Testnumber','threshold');
end
end