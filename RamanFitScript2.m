files = dir('F:\Dropbox\Dropbox\Raman\20180220-beryl-ascii\epoxy-690_1.txt');
files = dir('F:\Dropbox\Dropbox\Teeth\20171114-ascii\003-dentin-far-from-edge-side-toothB_1.txt');

df = 'F:\Dropbox\Dropbox\Raman\20180305-ascii\halfmicronscan-mat';
A = importdata(df);
x_vals = importdata('F:\Dropbox\Dropbox\Raman\20180305-ascii\xvals.txt');

show_individual = true;

% Region of interest in cm^-1
rois = 800;
roie = 1200;
tic;
laser_line = 633.46;

% Files are loaded into an array of data, if you have additional columns
% this will change which indices are wavenumber and intensity. e.g. winspec
% exports 3 columns, the third being intensity
x_index = 1;
y_index = 3;

% Background Removal Polynomial Order
poly_order = 2;

% Set guesses within for loop, we extract peak information which needs to
% be set

% Peak Height, Peak position, gaussian fwhm, lorentzian fwhm
% Peak scaling will clearly vary by sample composition
% Note that the guess amplitude should be a guess AFTER baseline
% subtraction
guess = [
    780000 962 10 5 ...
    25000 1005 10 10 ...
    100000 1045 15 15 ...
    100000 1070 15 15 ...
    0 -2.5e-3 1140];    

% Guess range to fit. 
guess_delta = [
    780000 10 10 10 ...
    25000 10 10 10 ...
    100000 10 10 10 ... 
    100000 5 10 10 ... 
    500 1e2 1800];

free_parameters = [
    1 2 3 4 ...
    5 6 7 8 ...
    9 10 11 12 ...
    13 14 15 16 ...
    17 18 19];


loop_size = 10;%size(A,1);
results = zeros(length(loop_size), length(free_parameters));
peak_fwhm = zeros(loop_size, 3);

for i=(size(A,1)-10):size(A,1)%loop_size%(size(A,1)-10):size(A,1)
    
    ydata = transpose(A(i,[1:1340]));

    % Background removal
    % Most of the time, the SIMPS background removal is good enough
    %{
    poly_x_vals = (1:size(a, x_index))';
    [poly_background, poly_struct, poly_mu] = polyfit(poly_x_vals, file_data(:, y_index), poly_order);
    poly_y_vals = polyval(p, t, [], mu);
    file_data(:, y_index) = file_data(:, y_index) - poly_y_vals;
    %}
    
    % Use this if your x values are in nm
    x_in_cm = x_vals;
    
    
    % Region of interest
    roi_start = find(x_in_cm >= rois, 1);
    roi_end = find(x_in_cm >= roie, 1);
    roi = [roi_start:roi_end];
    roguess = 900;%[1:1340];%find(x_in_cm >= 1200, 1);
        
    if show_individual
        hf=figure(i);
        clf;
    end
    
    guess(length(guess) - 2) = ydata(roguess); 

    high_guess = guess + guess_delta;
    low_guess = guess - guess_delta;
    
    
    [answer, g] = simps('fitvoigt', guess,(free_parameters),[],low_guess, high_guess, ydata(roi), x_in_cm(roi), 1);
    [f, G, fit, out] = fitvoigt(answer, ydata(roi), x_in_cm(roi), 1);
   
    if show_individual
        % -------- Plot Fit -------- %
        subplot(2,1,1)
        plot(out{1}, out{2}-out{5}, out{1}, out{3}-out{5});
        title('Evans Cool')
        ylabel('Intensity (arb. u.)')
        xlabel('Raman Shift (cm^-^1)')
    end
    
    % -------- Process data -------- %
    results(i, :) = answer; % Export this to a file if you want all data in a csv
    
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
            peak_fwhm(i, 3) = table_data(j, 1);
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
%{
figure
uitable('Data', peak_fwhm, 'ColumnName', {'Peak Pos', 'Voigt FWHM', 'Amp'});
figure
scatter(1:size(A,1), peak_fwhm(:,2));
figure
scatter(1:size(A,1), peak_fwhm(:,1));
figure
scatter(1:size(A,1), peak_fwhm(:,3));
%}

toc

