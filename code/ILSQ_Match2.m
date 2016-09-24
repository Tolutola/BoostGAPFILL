function [Lambda,scores,resnorm] = ILSQ_Match2(dA,U1,num_prediction,rr,cc,mask,integer,blacklist,solver)         
%  for solving the (integer) least square optimization
%%
global maxTime penID
l0 = rand(cc,1);
opts = optimset('MaxFunEval',inf,'MaxIter',Inf,'display','off');

dA1 = reshape(dA,rr^2,1);

blacklist=1-blacklist;   
if integer == 0    %integer == 0, solve relaxed continuous lsq
    if nargin<6     %no mask
        scores = lsqlin(U1,dA1,[],[],[],[],zeros(cc,1),ones(cc,1),l0,opts);    %solve the relaxed continuous linear least squared problem
    else         %using mask
        mask = mask(:);
        if strcmp(solver,'ibm_cplex')
            C=bsxfun(@times,U1,mask);
            d=dA1.*mask;
            C=[C;zeros(1,size(C,2))];
            C(end,penID)=1; % penalizing the addition of new metabolites
            d=[d;0];
                       
            [scores,resnorm] = cplexlsqlin(C,d,[],[],[],[],zeros(cc,1),blacklist,l0);    %solve using cplex, more stable
        else
           [scores,resnorm] = lsqlin(bsxfun(@times,U1,mask),dA1.*mask,[],[],[],[],zeros(cc,1),blacklist,l0,opts);  
        end
        
    end
    
    [~,I] = sort(scores,1,'descend');
    if num_prediction > 1
        if num_prediction> length(I), num_prediction=length(I); end
        Lambda(I(1:num_prediction)) = 1;    %only keep hl with top scores
        Lambda = logical(Lambda);
    else
        Lambda(scores>num_prediction) = 1;
        Lambda = logical(Lambda);
    end

else 
    mask = mask(:);
    options=cplexoptimset('MaxTime',maxTime);
    scores = cplexlsqbilp(bsxfun(@times,U1,mask),dA1.*mask,ones(1,size(U1,2)),num_prediction,[],[],[],[],[],options);  
    Lambda = logical(scores);
    
end
