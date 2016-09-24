function [Lambda,scores,resnorm] = ILSQ_Match_Flux(dA,U1,num_prediction,rr,cc,mask,integer,U,vlb,vub,blacklist)         
%  with flux constraints, need additional parameters U, vlb, vub
%%
global tempModel biomass_threshold
dA1 = reshape(dA,rr^2,1);

if nargin == 11
    blacklist = 1 - blacklist;
else
    blacklist = ones(cc,1);
end

mask = mask(:);
U1 = bsxfun(@times,U1,mask);
%U2 = [U1,(zeros(size(U1)))];

solutionFBA=optimizeCbModel(tempModel,[],'one');

dA1 = dA1.*mask;

lb = [zeros(cc,1);vlb;tempModel.lb];
ub = [blacklist;vub;tempModel.ub];
[nm,nr]=size(tempModel.S);

x0 = [rand(cc,1);randn(cc,1);solutionFBA.x];
Aeq=[zeros(nm,cc*2),tempModel.S];
beq=zeros(nm,1);

A=[zeros(1,cc*2),-tempModel.c'];
b=-biomass_threshold;


options = optimoptions('fmincon','Display','iter','MaxFunctionEvaluations',500000,'MaxIterations',20000,'SpecifyObjectiveGradient',true);

    
if integer == 0    %integer == 0, solve relaxed continuous lsq
    
    nonlcon = @(x) nlc(U,nr,x);
    result = fmincon(@obj,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
    scores = result(1:cc);
    
    [~,I] = sort(scores,1,'descend');
    if num_prediction > 1
        Lambda(I(1:num_prediction)) = 1;    %only keep hl with top scores
        Lambda = logical(Lambda);
    else
        Lambda(scores>num_prediction) = 1;
        Lambda = logical(Lambda);
    end
end


function [f,g] = obj(x)   % x = [lambda; v;vknown]
l = x(1:cc);
f = (U1*l-dA1)'*(U1*l-dA1);
g = 2*U1'*U1*l - 2*U1'*dA1;
g = [g;zeros(cc+nr,1)];
end

end



