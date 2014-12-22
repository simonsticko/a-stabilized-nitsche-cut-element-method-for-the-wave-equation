function[djdr]=dJdr(n,r)
djdr=.5*(besselj(n-1,r)-besselj(n+1,r));
end