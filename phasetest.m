clear;
close all;
period=2;
x=linspace(-period/2,period/2,200);
k=6*pi;
shift=.1;
f=@(x) cos(k*x)+sin(2*pi*x);
f1=f(x);
f2=f(x-shift);
plot(x,f1,'b',x,f2,'r');
[computedShift]=computeSignalShift(f1,f2,period)