function[an]=get_ansoft(n,gamma,R)
an=-dJdr(n,gamma*R)./dHdr(n,gamma*R);
end
