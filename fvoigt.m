function v=voigt(x,gd)
%
% VOIGHT(X,[GAMMA,DELTA]) returns the Voigt lineshape which
% is the convolution of a Lorentzian and a Gaussian, where
% the gaussian fwhm is delta and the lorentzian fwhm is
% gamma.  Uses the complex error function CERF.
% Gamma may be set = 0, delta may not.
gamma=gd(1)/2;
delta=gd(2)/sqrt(4*log(2));
pi=atan(1)*4;
i=sqrt(-1);
v=real(cerf((x+gamma*i)/delta));
v=v/max(v);