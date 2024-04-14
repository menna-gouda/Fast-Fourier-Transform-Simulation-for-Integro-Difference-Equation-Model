function err=nll_typeIII(r,x)

alpha=x(1);
d=x(2);

[p,c1]=typeIII(x);

t=(r./(sqrt(8).*alpha)).*sqrt(1+(1+16*(alpha./r).^2)); 
err=length(r)*(2*log(alpha)-log(c1))+...
    sum(-0.5*log(t)+0.5*log((t.^2)+2)+(d+1).*(t-atan(t)+((1./(4.*t)).*((r./alpha).^2))));

end