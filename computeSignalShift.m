%Given two functions f1,f2 period in x on with period period, such that
%f2(x)=f1(x-shift) computes the shift based on a fourier transform.
function[computedShift]=computeSignalShift(f1,f2,period)
F1=fft(f1);
F2=fft(f2);
[~,maxBin]=max(abs(F1));
aangle=angle(F2(maxBin)/F1(maxBin));
%Indices goes from 0 to N-1 in the sum, but matlab vectors goes from 1 to
%N.
mathematicalBin=maxBin-1;
computedShift=-period/(2*pi*mathematicalBin)*aangle;