% Function to evaluate complex error function for any
% z = x + i*y.  This uses Hui's rational approximation.
% Real part is Voigt function, imaginary part is dispersion
% (ie., complex index of refraction).  Works for both scalars
% and vectors.  Have not yet made function applicable to whole
% z-plane.
% by D. Holmgren, July 23 95.
function f = cerf(z)
% coefficients of rational approximation...      
      a0 = 37.24429446739879 + 0*i;
      a1 = 57.90331938807185 + 0*i;
      a2 = 43.16280063072749 + 0*i;
      a3 = 18.64649990312317 + 0*i;
      a4 = 4.67506018267650 + 0*i;
      a5 = 0.56418958297228 + 0*i;
      b0 = 37.2442945086 + 0*i;
      b1 = 99.9290005933 + 0*i;
      b2 = 118.6763981260 + 0*i;
      b3 = 80.6459493922 + 0*i;
      b4 = 33.5501020941 + 0*i;
      b5 = 8.2863279156 + 0*i;
% evaluate complex error function...      
      x = real(z); y = imag(z); nl = length(x);
% z is a scalar...      
      if nl == 1,
      zh = y - i*x;
      f=(((((a5*zh+a4)*zh+a3)*zh+a2)*zh+a1)*zh+a0)/((((((zh+b5)*zh+b4)*zh+b3)*zh+b2)*zh+b1)*zh+b0);
      end
% z is a vector...      
      if nl ~= 1,
      zh = y - i.*x;
      f=(((((a5 .*zh+a4).*zh+a3).*zh+a2).*zh+a1).*zh+a0)./((((((zh+b5).*zh+b4).*zh+b3).*zh+b2).*zh+b1).*zh+b0);
      end


