function [I,scores,hlpredicted] = HLpredict(hltrain,hltest,num_prediction,max_iter,blacklist,solver)
%  Usage: for hyperlink predictions in hypernetworks.  [~,I,scores] = HLpredict(hltrain,hltest,num_prediction)
% with 'MATBoost'--the Alternating Completion-Match algorithm,
%  --Input--
%  -hltrain: existing training hyperlinks, each row is a vertex, each column is a hyperlink.
%           hltrain(i,j)==1 iff the jth hyperlink involves vertex i.
%  -hltest: the test hyperlinks to predict which of them are positive, same format as hltrain
%  -num_prediction: how many number of hyperlinks do you want to keep as positive. If you input
%           value < 1, this will be used as a threshold. hls with scores > threshold will be
%           predicted as positive.
%  -ith_experiment: a positive integer number indicating you are running the ith experiment
%           (for parallelly running multiple experiments purpose, avoid file reading conflict)
%  -max_iter: used specifically for 'MATBoost' alg, the maximum allowed number of iterations
%  --Output--
%  -I: logical column indicator vector, I(i)==1 iff hltest(:,i) is predicted as positive (0 otherwise)
%  -scores: the column vector of scores that are assigned to each test hyperlink, higher scores
%           indicates higher probabilities to be positive hyperlinks
%  -hlpredicted: the predicted hyperlinks, same format as hltrain
%    -modeFlag: indicates whether BoostGapFill is run in the integrated mode or not
%
%  *author: Muhan Zhang, Washington University in St. Louis
%%

A = hltrain * hltrain';   % project the hyperlinks into the vertex space to form the adjacency matrix (with weights)
A = A - diag(diag(A));   % remove self adjacency
B = hltest * hltest';
B = B - diag(diag(B));
B = spones(B);    %the adjacency matrix of test hyperlinks, no use in hyperlink prediction, only for testing link prediction performance
k = 8; % the number of latent factors used in matrix factorization (used in libFM)


[I,scores] = ILSQ_wrapper(A,B,k,hltest,num_prediction,max_iter,blacklist,solver);
hlpredicted = hltest(:,I);

