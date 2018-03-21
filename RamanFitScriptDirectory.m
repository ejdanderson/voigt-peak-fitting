%addpath(genpath(strcat(fileparts(which(mfilename)), '/natsortfiles')));
%files = dir('/Users/evan/Dropbox/Teeth/20171114-ascii/*enamel*.txt');
files = dir('F:\Dropbox\Dropbox\Raman\20180305-ascii\data\side-transect\*.txt');
%filenames = natsortfiles({files.name});
    
% Region of interest in cm^-1
rois = 900;
roie = 1050;

show_individual = true;

laser_line = 633.18;

% Files are loaded into an array of data, if you have additional columns
% this will change which indices are wavenumber and intensity. e.g. winspec
% exports 3 columns, the third being intensity
x_index = 1;
y_index = 3;

% Background Removal Polynomial Order
poly_order = 2;

% Set guesses within for loop, we extract peak information which needs to
% be set

    free_parameters = [
        1 2 3 4 ... % 2 3 4
        5 6 7];

results = zeros(length(files), length(free_parameters));
peak_fwhm = zeros(length(files), 2);
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
        250000 962 20 20 ... % 3
        0 -2.5e-3 1140]; % Additional baseline fitting
    %{
               170000	1045 20.6008	29.9878 ... % 3
        170000	1070 23.1763 6.1766e-04 ...% 4
        %}
    
    % Guess range to fit. 
    guess_delta = [
        200000 5 10 10 ... % 3
        0.5*peak_scale 1e2 1800];
    %{
                70000 10 10 10 ... 
        70000 10 15 15 ... % 4
        %}
    

%{
                5 6 7 8 ...
        9 10 11 12 ... 
        %}
    
    % Use this if your x values are in cm^-1
    x_in_cm = 10^7*(1/laser_line - 1./file_data(:, x_index));
    
    % Background removal
    % Most of the time, the SIMPS background removal is good enough
    
    poly_x_vals = (1:size(file_data, x_index))';
    [poly_background, poly_struct, poly_mu] = polyfit(poly_x_vals, file_data(:, y_index), poly_order);
    poly_y_vals = polyval(poly_background, poly_x_vals, [], poly_mu);
    %file_data(:, y_index) = file_data(:, y_index) - poly_y_vals;
    
    
    % Region of interest
    roi_start = find(x_in_cm >= rois, 1);
    roi_end = find(x_in_cm >= roie, 1);
    roi = [roi_start:roi_end];
    roguess = find(x_in_cm >= 1200, 1);
    
    if show_individual
        hf=figure;
        clf;
    end
    
    guess(length(guess) - 2) = file_data(roguess, y_index); 
    
    guess_delta(length(guess) - 2) = file_data(roguess, y_index); 

    high_guess = guess + guess_delta;
    low_guess = guess - guess_delta;
    
    [answer, g] = simps('fitvoigt', guess,(free_parameters),[],low_guess, high_guess, file_data(roi, y_index), x_in_cm(roi), 1);
    [f, G, fit, out] = fitvoigt(answer, file_data(roi, y_index), x_in_cm(roi), 1);
    
    if show_individual
        % -------- Plot Fit -------- %
        subplot(2,1,1)
        plot(out{1}, out{2},out{1}, out{3} );

        title(files(i).name)
        ylabel('Intensity (arb. u.)')
        xlabel('Raman Shift (cm^-^1)')
    end
    
    % -------- Process data -------- %
    results(i, :)=answer; % Export this to a file if you want all data in a csv
    
    % No need to show baseline fit, hence -3 parameters
    table_rows = (length(answer) - 3) / 4;
    table_data = zeros(table_rows, 5);
    
    for j=1:table_rows
        index = (j - 1) * 4 + 1;
        gauss_fwhm = answer(index + 3);
        lorentz_fwhm = answer(index + 2);
        table_data(j, 1) = answer(index);
        table_data(j, 2) = answer(index + 1);
        table_data(j, 3) = lorentz_fwhm;
        table_data(j, 4) = gauss_fwhm;
        %voigt peak width
        table_data(j, 5) = gauss_fwhm*(1-2.0056*1.0593+sqrt((lorentz_fwhm/gauss_fwhm)^2+2*1.0593*lorentz_fwhm/gauss_fwhm+2.0056^2*1.0593^2));
        
        if j==1
            peak_fwhm(i, 1) = answer(index+1);
            peak_fwhm(i, 2) = table_data(j, 5);
        end
    end
  
    if show_individual
        % -------- Plot data -------- %
        % MATLAB trickery, produce a subplot, get its position and delete it.
        % Then put the uitable into the subplot position
    
        sp = subplot(2, 1, 2);
        pos = get(sp, 'Position');
        un = get(sp, 'Units');
        delete(sp);
        cnames={'Amp', 'Position', 'Lorentzian FWHM', 'Gaussian FWHM', 'Voigt FWHM'};
        t = uitable(hf, 'Data', table_data, 'ColumnName', cnames, 'Units', un, 'Position', pos);
    end
    
end


%figure
%uitable('Data', results);
%figure
%uitable('Data', peak_fwhm, 'ColumnName', {'Peak Pos', 'Voigt FWHM'});
figure
scatter(1:length(files), peak_fwhm(:,2));




