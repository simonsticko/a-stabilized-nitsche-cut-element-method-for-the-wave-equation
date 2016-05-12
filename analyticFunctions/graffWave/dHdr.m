function[dhdr]=dHdr(n,r)
dhdr=.5*(besselh(n-1,2,r)-besselh(n+1,2,r));
end
