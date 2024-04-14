function err=nll_bessel(r,x)
alpha=x(1);

[p,c1]=bessel(x);

K=besselk(0,r/alpha);
err=length(r)*(2*log(alpha)-log(c1))+sum(-log(K));

end