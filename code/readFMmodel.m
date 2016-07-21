function [Feature,trainlabels,Feature1] = readFMmodel(hltrain,hltest,wtrain,k)

% w=dlmread('data/FMmodelw.txt');
% v=dlmread('data/FMmodelv.txt');
% w0=dlmread('data/FMmodelw0.txt');
% 
% %% for testing if model is read correctly (by comparing the out_1.txt)
% for runthis = 1:1
% k=1;
% g=2;
% a = v(:,k);
% c = v(:,g);
% b = sum((a+c).^2 - a.^2-c.^2)/2 + w0 + w(k) + w(g);
% b1 = a'*c + w0 + w(k) + w(g);
% assert(b-b1<0.01);
% end
% V = [w';v];    %to combine v and w into a single matrix

A = spones(wtrain);
D = eye(size(A));
dis = 10*ones(size(A));
for i = 1:5
    D = spones(D*A);
    tmp = D*i;
    dis = min(dis,tmp);
end
dis(dis==10) = 0;
[V,D] = eigs(dis);
k = min(k,size(D,1));
V(:,k+1:end) = [];
D = D(1:k,1:k);
Z = V*(D^2);
V = Z';

%% generate training examples
%fid=fopen(strcat(datapath,'LatentFactorTrainingExamples.txt'),'w+');
n = 2;   %n times negative samples
trainlabels = [ones(1,size(hltrain,2)),zeros(1,n*size(hltrain,2))]';

%to generate positive training examples from trainr
Feature = [];
for i = 1:size(hltrain,2)
    r = hltrain(:,i);
    example = V(:,logical(r));
    feature = [];
    for j=1:size(V,1)
        feature = [feature,entropy(example(j,:))];
    end
    Feature = [Feature;feature];
    %example = [reshape(example,1,[])];
    %fprintf(fid,strcat(num2str(example),'\r\n'));
end
%to sample randomly from empties in hltest to form negative examples
for nn = 1:n
for i = 1:size(hltrain,2)
    r = hltrain(:,i);
    num_meta = sum(r);
    
    while 1
        c = zeros(size(hltrain,1),1);
        perm = randperm(size(hltrain,1));
        c(perm(1:num_meta)) = 1;
        lia = 0;   %faster, not too much influence because sparsity
        %lia = ismember(c',hltrain','rows');    %check if c in US
        if lia == 0
            i;
            break
        end
    end
    example = V(:,logical(c));
    %example = [0,reshape(example,1,[])];
    feature = [];
    for j=1:size(V,1)
        feature = [feature,entropy(example(j,:))];
    end
    Feature = [Feature;feature];
    
    %fprintf(fid,strcat(num2str(example),'\r\n'));
end
end
%fclose(fid);

%% generate testing examples
Feature1 = [];
%fid2=fopen(strcat(datapath,'LatentFactorTestingExamples.txt'),'w+');
for i = 1:size(hltest,2)
    r = hltest(:,i);
    example = V(:,logical(r));
    feature = [];
    for j=1:size(V,1)
        feature = [feature,entropy(example(j,:))];
    end
    Feature1 = [Feature1;feature];
    %example = [labels(i),reshape(example,1,[])];
    %fprintf(fid2,strcat(num2str(example),'\r\n'));
end
