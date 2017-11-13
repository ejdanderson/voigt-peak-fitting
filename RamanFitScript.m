files = dir('data/ping/line_data/normal_tooth/side_area/e1.txt');

% Region of interest in cm^-1
rois = 400;
roie = 1700;

% Files are loaded into an array of data, if you have additional columns
% this will change which indices are wavenumber and intensity. e.g. winspec
% exports 3 columns, the third being intensity
x_index = 1;
y_index = 2;

% Background Removal Polynomial Order
poly_order = 2;

% Set guesses within for loop, we extract peak information which needs to
% be set

for i=1:length(files)
   
    % This may be slightly different depending on what OS you are using
    file_data = load(strcat(files(i).folder, '/', files(i).name));
    
    % Get the peak scale for guesses off of the largest
    peak_scale = max(file_data(:, y_index)) / 4;

    % Peak Height, Peak position, gaussian fwhm, lorentzian fwhm
    % Peak scaling will clearly vary by sample composition
    guess = [
        1.0*peak_scale 433 20 20 ... % 1
        1.2*peak_scale 578 20 20 ... % 2
        4.0*peak_scale 962 10 5 ... % 3
        1.2*peak_scale 1045 15 15 ...% 4
        1.0*peak_scale 1070 15 15 ...% 5
        1.2*peak_scale 1220 20 20 ...% 6
        1.2*peak_scale 1412 20 20 ...% 7
        1.0*peak_scale 1620 20 20 ...% 8
        0 -2.5e-3 1140]; % Additional baseline fitting
    
    % Free parameters based on guess above add/remove to set free or fix
    free_parameters = [
        1 ... % 2 3 4
        5 ... % 6 7 8
        9 10 11 ... % 12
        13 ... % 14 15 16
        17 ... % 18 19 20 
        21 ... % 22 23 24
        25 ... % 26 27 28 
        29 ... % 30 31 32  
        33 34 35];

    % Guess range to fit. 
    guess_delta = [
        0.5*peak_scale 10 20 20 ... % 1
        1.0*peak_scale 10 20 20 ... % 2
        0.5*peak_scale 10 10 10 ... % 3
        1.0*peak_scale 10 15 15 ... % 4
        1.0*peak_scale 10 15 15 ... % 5
        1.0*peak_scale 10 20 20 ... % 6
        1.0*peak_scale 10 20 20 ... % 7
        1.0*peak_scale 10 20 20 ... % 8
        1.0*peak_scale 1e2 1800];

    
    
    % Background removal
    % Most of the time, the SIMPS background removal is good enough
    %{
    poly_x_vals = (1:size(a, x_index))';
    [poly_background, poly_struct, poly_mu] = polyfit(poly_x_vals, file_data(:, y_index), poly_order);
    poly_y_vals = polyval(p, t, [], mu);
    file_data(:, y_index) = file_data(:, y_index) - poly_y_vals;
    %}
    
    % Use this if your x values are in nm
    % cm = 10^7 * ( -1 ./ a(:, x_index) + 1 / 634.578);
    
    % Use this if your x values are in cm^-1
    x_in_cm = file_data(:, x_index);
    
    % Region of interest
    roi_start = find(x_in_cm >= rois, 1);
    roi_end = find(x_in_cm >= roie, 1);
    roi = [roi_start:roi_end];
    roguess = find(x_in_cm >= 1200, 1);
    
    figure(i)
    clf;
    
    guess(length(guess) - 2) = file_data(roguess, y_index); 

    high_guess = guess + guess_delta;
    low_guess = guess - guess_delta;
    
    
    [answer, g] = simps('fitvoigt', guess,(free_parameters),[],low_guess, high_guess, file_data(roi, y_index), x_in_cm(roi), 1);
    [f, G, fit, out] = fitvoigt(answer, file_data(roi, y_index), x_in_cm(roi), 1);
    
    plot(out{1}, out{2}, out{1}, out{3});

    results(i, :)=answer;
    
    title(files(i).name)
    ylabel('Intensity (arb. u.)')
    xlabel('Raman Shift (cm^-^1)')
end

