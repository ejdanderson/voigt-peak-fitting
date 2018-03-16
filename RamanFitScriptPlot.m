% @TODO Protect the user from over-writing results at the bottom of the file?
% It would be easy to run a 10 hour fit which auto-exports the results,
% then run a subset of points which over-writes the previous output.

df = 'F:\Dropbox\Dropbox\Raman\20180305-ascii\tip-halfmicron';
x_in_cm = importdata('F:\Dropbox\Dropbox\Raman\20180305-ascii\xvals.txt');

%df = '~/Dropbox/Raman/20180305-ascii/tip-halfmicron';
%x_in_cm = importdata('~/Dropbox/Raman/20180305-ascii/xvals.txt');

delimiterIn = '\t';
headerlinesIn = 3;
A = importdata(df, delimiterIn, headerlinesIn);

num_steps_x = str2num(A.textdata{2});
num_steps_y = str2num(A.textdata{1});

% Which peaks to create a 2-D map for. These correspond to the peaks
% defined in the guess variable below
peaks_to_plot = [1,3];

fwhm_grids = zeros(num_steps_x, num_steps_y, length(peaks_to_plot));
pp_grids = zeros(num_steps_x, num_steps_y, length(peaks_to_plot));
intensity_grids = zeros(num_steps_x, num_steps_y, length(peaks_to_plot));

% Wether or not to plot the fit for every data point (usually you dont want
% this set if you have more than 20 data points)
% Manually set the x and y ranges in the for loops to test your data
show_individual = true;

% x and y Range to loop over. TO get a sample of data, you might do
% something like: x_range=0:40:399, y_range=1:1
x_range = 0:40:399;%0:num_steps_x-1;
y_range = 1:1;%num_steps_y

% Region of interest in cm^-1 (WHAT PEAKS ARE YOU LOOKING AT?, put them in the ROI)
rois = 900;
roie = 1050;

% Files are loaded into an array of data, if you have additional columns
% this will change which indices are wavenumber and intensity. e.g. winspec
% exports 3 columns, the third being intensity
x_index = 1;
y_index = 3;

% Set guesses within for loop, we extract peak information which needs to
% be set

% Peak Height, Peak position, gaussian fwhm, lorentzian fwhm
% Note that the guess amplitude should be a guess AFTER baseline
% subtraction, you'll likely have to run this once, look at a couple fits
% and make adjustments

guess = [
    18000 962 20 20 ...
    100000 1045 15 15 ...
    100000 1070 15 15 ...
    0 -2.5e-3 1140];    

% Guess range to fit. 
guess_delta = [
    17000 10 10 10 ...
    100000 10 10 10 ...
    100000 10 10 10 ...
    500 10 1800];

free_parameters = [
    1 2 3 4 ...
    5 6 7 8 ... 
    9 10 11 12 ...
    13 14 15];

loop_size = size(A,1);
results = zeros(length(loop_size), length(free_parameters));
peak_fwhm = zeros(loop_size, 3);

for i = 0:num_steps_y-1
    for j=1:num_steps_x
        k = j+num_steps_x*i;
        l = k;
        Spectrum(k,:)=A.data(l,:);
    end
end

% Times how long the rest of the code takes
runTime = tic;

for x = x_range
    for y = y_range
        tic;
        index = (num_steps_y * x) + y;
        index
        if mod(x, 2) == 0
            y_pos = num_steps_y - y + 1;
        else
            y_pos = y;
        end
        
        y_intensity = transpose(Spectrum(index, :));

        % Region of interest
        roi_start = find(x_in_cm >= rois, 1);
        roi_end = find(x_in_cm >= roie, 1);
        roi = [roi_start:roi_end];
        roguess = find(x_in_cm >= 994, 1);

        if show_individual
            hf=figure(index);
            clf;
        end

        % Set baseline param
        guess(length(guess) - 2) = y_intensity(roguess); 
        guess_delta(length(guess) - 2) = y_intensity(roguess)/2;
        
        guess(length(guess)) = max(y_intensity);
        guess_delta(length(guess)) = max(y_intensity)/2;
        
        %Setup guesses
        high_guess = guess + guess_delta;
        low_guess = guess - guess_delta;

        % Run the fitting procedure    
        [answer, g] = simps('fitvoigt', guess,(free_parameters),[],low_guess, high_guess, y_intensity(roi), x_in_cm(roi), 1);
        [f, G, fit, out] = fitvoigt(answer, y_intensity(roi), x_in_cm(roi), 1);
        
        
        if show_individual
            % -------- Plot Fit -------- %
            subplot(2,1,1)
            plot(out{1}, out{2}-out{5}, out{1}, out{3}-out{5});
            title(strcat('x=', int2str(x), ', y=', int2str(y_pos)))
            ylabel('Intensity (arb. u.)')
            xlabel('Raman Shift (cm^-^1)')
        end

        % -------- Process data -------- %
        results(index, :) = answer; % Export this to a file if you want all data in a csv

        % No need to show baseline fit, hence -3 parameters
        table_rows = (length(answer) - 3) / 4;
        table_data = zeros(table_rows, 5);
        
        for j=1:table_rows
            index = (j - 1) * 4 + 1;
            gauss_fwhm = answer(index + 3);
            lorentz_fwhm = answer(index + 2);
            table_data(j, 1) = answer(index); % Intensity
            table_data(j, 2) = answer(index + 1); % Peak Pos
            table_data(j, 3) = lorentz_fwhm;
            table_data(j, 4) = gauss_fwhm;
            %voigt peak width
            table_data(j, 5) = gauss_fwhm*(1-2.0056*1.0593+sqrt((lorentz_fwhm/gauss_fwhm)^2+2*1.0593*lorentz_fwhm/gauss_fwhm+2.0056^2*1.0593^2));

            if ismember(j, peaks_to_plot)
                pp_grids(x+1, y_pos, j) = table_data(j, 2);
                fwhm_grids(x+1, y_pos, j) = table_data(j, 5);
                intensity_grids(x+1, y_pos, j) = table_data(j, 1);
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
        toc
    end % for y = y_range
end % for x = x_range


[fpath, fname] = fileparts(df);


for j=1:length(peaks_to_plot)
    filename_base = fullfile(fpath, fname);
    
    i = peaks_to_plot(j);
    i_str = int2str(i);
    
    figure
    imagesc(transpose(pp_grids(:, :, i)))
    pbaspect([1 num_steps_y/num_steps_x 1])
    colorbar
    title(strcat('Peak Position, peak', i_str));
    savefig(strcat(filename_base, '_peak_pos_peak_', i_str, '.fig'));
    dlmwrite(strcat(filename_base, '_peak_pos_peak_', i_str, '.txt'), pp_grids(:, :, i));

    figure
    imagesc(transpose(fwhm_grids(:, :, i)))
    pbaspect([1 num_steps_y/num_steps_x 1])
    colorbar
    title(strcat('FWHM, peak', i_str));
    savefig(strcat(filename_base, '_fwhm_peak_', i_str, '.fig'))
    dlmwrite(strcat(filename_base, '_fwhm_peak_', i_str, '.txt'), fwhm_grids(:, :, i));

    figure
    imagesc(transpose(intensity_grids(:, :, i)))
    pbaspect([1 num_steps_y/num_steps_x 1])
    colorbar
    title(strcat('Intensity, peak', i_str));
    savefig(strcat(filename_base, '_intensity_peak_', i_str, '.fig'))
    dlmwrite(strcat(filename_base, '_intensity_peak_', i_str, '.txt'), intensity_grids(:, :, i));
end

dlmwrite(fullfile(fpath, strcat(fname, '_all_data.txt')), results)

toc(runTime)

