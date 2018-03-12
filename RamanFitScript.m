addpath(genpath(strcat(fileparts(which(mfilename)), '\natsortfiles')));
files = dir('F:\Dropbox\Dropbox\Raman\pingdata\side\*.txt');
filenames = natsortfiles({files.name});

% Region of interest in cm^-1
rois = 800;
roie = 1200;

laser_line = 633.18;

show_individual = false;

% Files are loaded into an array of data, if you have additional columns
% this will change which indices are wavenumber and intensity. e.g. winspec
% exports 3 columns, the third being intensity
x_index = 1;
y_index = 2;

% Change this to match the peaks, Amp, Pos, Lorentzian FWHM, Gauss FWHM
% Last 3 are background, should not need to adjust these
guess = [
    4200 962 15 15 ... % 3
    600	1045	20.6008	29.9878 ... % 3
    400	1070 23.1763 6.1766e-04 ...% 4
    ... %10000 1075 15 15 ...% 5
    ... %2000 1249 30 20 ...% 6
    ... %2000 1270 25 20 ...% 6
    ... %3000 1460 20 20 ...% 7
    ... %3000 1680 20 20 ...% 8
    0 -2.5e-3 1140]; % Additional baseline fitting

% Guess range to fit. 
guess_delta = [
    2000 5 10 10 ... % 3
    400 10 10 10 ... 
    300 10 15 15 ... % 4
   ... % 5000 10 15 15 ... % 5
   ... % 0.5*peak_scale 10 20 20 ... % 6
   ... % 0.5*peak_scale 10 20 20 ... % 7
   ... % 0.5*peak_scale 10 20 20 ... % 7
   ... % 0.5*peak_scale 10 20 20 ... % 8
   500 1e2 1800];

free_parameters = [
    1 2 3 4 ... % 2 3 4
    5 6 7 8 ...
    9 10 11 12 ... 
    ... %13 14 15 16 ... % 14 15 16
    ... %17 18 19 20 ... % 18 19 20 
    ... %21 22 23 24 ... % 22 23 24
    ... %25 26 27 28 ...
    ... %29 30 31 32 ...
    13 14 15];


results = zeros(length(files), length(free_parameters));
peak_fwhm = zeros(length(files), 2);
for i=1:length(files)
   
    % This may be slightly different depending on what OS you are using
    file_data = load(strcat(files(i).folder, '/', filenames{i}));
   
    % Use this if your x values are in nm
    x_in_cm = file_data(:, x_index); %10^7 * ( -1 ./ file_data(:, x_index) + 1 / laser_line);
    
    % Use this if your x values are in cm^-1
    %x_in_cm = file_data(:, x_index);
    
    % Region of interest
    roi_start = find(x_in_cm >= rois, 1);
    roi_end = find(x_in_cm >= roie, 1);
    roi = [roi_start:roi_end];
    roguess = find(x_in_cm >= 1200, 1);
    
    if show_individual
        hf=figure(i);  
        clf;
    end
    
    guess(length(guess) - 2) = file_data(roguess, y_index); 
    high_guess = guess + guess_delta;
    low_guess = guess - guess_delta;
  
    [answer, g] = simps('fitvoigt', guess,(free_parameters),[],low_guess, high_guess, file_data(roi, y_index), x_in_cm(roi), 1);
    [f, G, fit, out] = fitvoigt(answer, file_data(roi, y_index), x_in_cm(roi), 1);
    
    if show_individual
        % -------- Plot Fit -------- %
        subplot(2,1,1)
        plot(out{1}, out{2}-out{5}, out{1}, out{3}-out{5});

        title(filenames{i})
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
        table_data(j, 1) = answer(index);
        table_data(j, 2) = answer(index + 1);
        lorentz_fwhm = answer(index + 2);
        table_data(j, 3) = answer(index + 2);
        gauss_fwhm = answer(index + 3);
        table_data(j, 4) = answer(index + 3);
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


figure
uitable('Data', results);
figure
uitable('Data', peak_fwhm, 'ColumnName', {'Peak Pos', 'Voigt FWHM'});
figure
scatter(1:length(files), peak_fwhm(:,2));
title("Side FWHM (Ping's data, Evan's Fit")

