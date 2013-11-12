function power=two_sided_t_pow(ES,sig,bias,n,al)

df= 2*n-2;
t=tinv(1-al/2,df);
omega=((ES+bias)*sqrt(n/2))./sig;
power=(1-tcdf(t-omega,df))+tcdf(-t-omega,df);