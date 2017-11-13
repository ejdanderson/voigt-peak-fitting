function [f,G,fit,out] = fitvoigt(p,y,x,sigma),
%fitvoight Fits an arbitrary number of voigt functions, a gaussian, and a baseline,
%returns chi2 and G (vector of constraints for simps.m)
%
%  usage [f,G,output] = fitvoight([a1,pos1,tau1l,tau1g,a2,pos2,tau2l, ...
%		    tau2g,a3,pos3,tau3g,offset],data,abcissa,std);
%
%  f    chi2
%  G    used by simps
%  fit  fit to data
%  out  output for plotting

    p_length = length(p);
    for i=0:((p_length - 3)/4) - 1
        p_index = (i * 4) + 1;
        amp(i+1) = p(p_index);       % peak amplitude
        pos(i+1) = p(p_index + 1);   % peak position
        tau_l(i+1) = p(p_index + 2); % lorentzian width of peak (fwhm)
        tau_g(i+1) = p(p_index + 3); % gaussian width of peak (fwhm)
    end
    
    
    offset = p(p_length - 2);
    curve = p(p_length - 1);
    most = p(p_length);
  
    % insist a1, a2, a3, tau1l, tau2l > 0
    % @TODO is this really needed? Should be for all amps/taus
	if (amp(1) < 0) | (amp(2) < 0) | (amp(3) < 0), 
        f=1e18;
        fit=zeros(size(y));
        G=[];
        disp('fixing amp'); 
        return;
    end 
  
	if (tau_l(1) < 0) | (tau_l(2) < 0), 
        f=1e25;
        fit=zeros(size(y));
        G=[];
        disp('fixing 1l'); 
        return
    end
  
    length(amp);
    voigts = zeros(length(x), length(amp));
    for i=1:length(amp)
        voigts(:,i) = amp(i) * fvoigt(x - pos(i), [tau_l(i), tau_g(i)]);
    end
    
  %{
  voight1 = a1*voight(x-pos1,[tau1l,tau1g]);
  voight2 = a2*voight(x-pos2,[tau2l,tau2g]);
  voight3 = a3*voight(x-pos3,[tau3l,tau3g]);
  voight4 = a4*voight(x-pos4,[tau4l,tau4g]);
  voight5 = a5*voight(x-pos5,[tau5l,tau5g]);
  %}
   
	baseline= offset + curve * (x - most).^2;
    
    fit = zeros(length(x),1);
	for i=1:length(amp)
        fit = fit + voigts(:, i);
    end
    
    fit = fit + baseline;
    f = sum(((y - fit)./sigma).^2);
    G=[];
    out = {x, y, fit, voigts, baseline};