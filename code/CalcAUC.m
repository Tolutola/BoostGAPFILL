function [ precision,dA,wdA ] = CalcAUC( train, test, sim, ~ )    %% you can change precision to auc here if you want auc output
    %% output wdA, dA
    selected = round(0.1*size(sim,1)*size(sim,2));    %the threshold number, you can decide whether to remain all sim scores in dA/wdA or only keep the top "selected" number
	sim(isnan(sim)) = 0;
    sim(isinf(sim)) = 0;
    wtrain = train;    %the weighted train adjacency matrix
    wdA = sim;
    train = spones(train);   %the 0/1 train adjacency matrix
    %sim = triu(sim - sim.*train) - diag(diag(sim));      
    sim1 = triu(sim) - diag(diag(sim));         %we keep all sim scores here (include train links)
   
    %notrain = sim.*(train==0);           %to delete train links' scores, however, it will also delete many test links (due to train,test overlap)
    simv = sim1(:);          %keep all predicted scores
    [y,i] = sort(simv,1,'descend');
    thres = y(selected);
    dA = sim1>=thres;
    dA = dA + dA' + diag(diag(sim));
    %wdA = sim1;    %without threshold, i.e, keep weights for all links
	%wdA = wdA - wtrain;  %whether to delete train links scores
    %wdA(wdA<thres)=0; display('keep only selected 1000');   %only keep links with weights >= thres
    %wdA = wdA + wdA';
    
    %% calculate the link prediction precision/auc
    
    non = 1 - train - test - eye(max(size(train,1),size(train,2)));   %only keep links outside train and test as non-exist links
    non(non<0)=0;         %the overlapped links do not belong to non
    test = triu(test);
    non = triu(non);
    
    test_pre = sim .* test;
    non_pre = sim .* non;
    
    test_data =  test_pre( test ~= 0 )';      % the vector of all test scores
    non_data =  non_pre( non ~= 0 )';      %the vector of all non scores
    
    
    % use perfcurve to calculate prediction auc
    for runthis = 1:0
        labels = [ones(1,size(test_data,2)), zeros(1,size(non_data,2))];
        scores = [test_data, non_data];
        [X,Y,T,auc] = perfcurve(labels, scores, 1);
    end
    
    % calculate confusion matrix and link prediction precision
    for runthis = 1:1
        labels = [ones(1,size(test_data,2)), 2*ones(1,size(non_data,2))];       %be careful when test scores = non scores, selected will all come from top entries in labels,
        scores = [test_data, non_data];                                         %then these selected will all be classified as positive, but actually not.
        [y,i] = sort(scores,2,'descend');
        y = y(:,1:selected);
        i = i(:,1:selected);
        thres = y(selected);
        
        g1 = labels(i);
        g2 = ones(1,selected);
        C = confusionmat(g1,g2);
        precision = C(1,1)/selected;           %the precision score at "selected"
        recall = C(1,1)/length(test_data);     %the recall score at "selected"
    end
    
    % to draw the ROC curve
    for runthis = 1:0
        figure(1);
        plot(X,Y,'r','LineWidth',4)
        xlabel('False positive rate')
        ylabel('True positive rate')
        title('ROC');
    end
    
end
