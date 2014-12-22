function[Rn]=getRn(n,r,gamma,R)
[an]=get_ansoft(n,gamma,R);
Rn=besselj(n,gamma*r)+an*besselh(n,2,gamma*r);
end