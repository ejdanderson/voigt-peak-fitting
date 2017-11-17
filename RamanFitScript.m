files = dir('data/009*dentin*.txt');

% Region of interest in cm^-1
rois = 800;
roie = 1800;

laser_line = 632.46;

% Files are loaded into an array of data, if you have additional columns
% this will change which indices are wavenumber and intensity. e.g. winspec
% exports 3 columns, the third being intensity
x_index = 1;
y_index = 3;

% Background Removal Polynomial Order
poly_order = 2;

% Set guesses within for loop, we extract peak information which needs to
% be set

results = zeros(length(files), 35);
for i=1:length(files)
   
    % This may be slightly different depending on what OS you are using
    file_data = load(strcat(files(i).folder, '/', files(i).name));
    
    % Get the peak scale for guesses off of the largest
    peak_scale = max(file_data(:, y_index)) / 4;

    % Peak Height, Peak position, gaussian fwhm, lorentzian fwhm
    % Peak scaling will clearly vary by sample composition
    % Note that the guess amplitude should be a guess AFTER baseline
    % subtraction
    guess = [
        4.0*peak_scale 962 10 5 ... % 3
        2000 1005 10 5 ... % 3
        3000 1040 15 15 ...% 4
        10000 1075 15 15 ...% 5
        2000 1249 30 20 ...% 6
        2000 1270 25 20 ...% 6
        3000 1460 20 20 ...% 7
        3000 1680 20 20 ...% 8
        0 -2.5e-3 1140]; % Additional baseline fitting
    
    % Guess range to fit. 
    guess_delta = [
        1.0*peak_scale 10 10 10 ... % 3
        1000 10 10 10 ... 
        3000 10 15 15 ... % 4
        5000 10 15 15 ... % 5
        0.5*peak_scale 10 20 20 ... % 6
        0.5*peak_scale 10 20 20 ... % 7
        0.5*peak_scale 10 20 20 ... % 7
        0.5*peak_scale 10 20 20 ... % 8
        0.5*peak_scale 1e2 1800];
    
    free_parameters = [
        1 2 3 4 ... % 2 3 4
        5 6 7 8 ... % 12
        9 10 11 12 ... 
        13 14 15 16 ... % 14 15 16
        17 18 19 20 ... % 18 19 20 
        21 22 23 24 ... % 22 23 24
        25 26 27 28 ...
        29 30 31 32 ...
        33 34 35];
    
    % Background removal
    % Most of the time, the SIMPS background removal is good enough
    %{
    poly_x_vals = (1:size(a, x_index))';
    [poly_background, poly_struct, poly_mu] = polyfit(poly_x_vals, file_data(:, y_index), poly_order);
    poly_y_vals = polyval(p, t, [], mu);
    file_data(:, y_index) = file_data(:, y_index) - poly_y_vals;
    %}
    
    % Use this if your x values are in nm
    x_in_cm = 10^7 * ( -1 ./ file_data(:, x_index) + 1 / laser_line);
    
    % Use this if your x values are in cm^-1
    %x_in_cm = file_data(:, x_index);
    
    % Region of interest
    roi_start = find(x_in_cm >= rois, 1);
    roi_end = find(x_in_cm >= roie, 1);
    roi = [roi_start:roi_end];
    roguess = find(x_in_cm >= 1200, 1);
    
    hf=figure(i)
    
    clf;
    
    guess(length(guess) - 2) = file_data(roguess, y_index); 

    high_guess = guess + guess_delta;
    low_guess = guess - guess_delta;
    
    
    [answer, g] = simps('fitvoigt', guess,(free_parameters),[],low_guess, high_guess, file_data(roi, y_index), x_in_cm(roi), 1);
    [f, G, fit, out] = fitvoigt(answer, file_data(roi, y_index), x_in_cm(roi), 1);
    
    % -------- Plot Fit -------- %
    subplot(2,1,1)
    plot(out{1}, out{2}-out{5}, out{1}, out{3}-out{5});

    title(files(i).name)
    ylabel('Intensity (arb. u.)')
    xlabel('Raman Shift (cm^-^1)')
    
    % -------- Process data -------- %
    results(i, :)=answer; % Export this to a file if you want all data in a csv
    
    % No need to show baseline fit, hence -3 parameters
    table_rows = (length(answer) - 3) / 4;
    table_data = zeros(table_rows, 4);
    for j=1:table_rows
        index = (j - 1) * 4 + 1;
        table_data(j, 1) = answer(index);
        table_data(j, 2) = answer(index + 1);
        table_data(j, 3) = answer(index + 2);
        table_data(j, 4) = answer(index + 3);
    end
  
    % -------- Plot data -------- %
    % MATLAB trickery, produce a subplot, get its position and delete it.
    % Then put the uitable into the subplot position
    sp = subplot(2, 1, 2);
    pos = get(sp, 'Position');
    un = get(sp, 'Units');
    delete(sp);
    cnames={'Amp', 'Position', 'Lorentzian FWHM', 'Gaussian FWHM'};
    t = uitable(hf, 'Data', table_data, 'ColumnName', cnames, 'Units', un, 'Position', pos);
end


