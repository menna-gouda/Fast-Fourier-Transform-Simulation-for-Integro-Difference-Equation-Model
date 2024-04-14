function err=nll_typeII(r,x)

alpha=x(1);
d=x(2);

[p,c1]=typeII(x);

err=length(r)*(2*log(alpha)-log(c1)) +sum((d+1).*(r/alpha)-0.5.*log(r/alpha)-d.*log(1+0.5.*(r/alpha)));

end