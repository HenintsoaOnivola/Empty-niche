function K=carryin_capacity(k0,x0,a,x)
        K=k0*(1-((x-x0)/a).^2);
        K(abs(x-x0)>a)=0;
 end