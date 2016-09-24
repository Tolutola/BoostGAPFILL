function [c,ceq] = nlc(U,nr, x)   % nonlinear constraint flux
[rr,~] = size(x);
rr=rr-nr;
lambda = x(1:rr/2);
v = x(rr/2+1:rr);
lv = lambda.*v;
ceq = U*lv;
c = [];
end
