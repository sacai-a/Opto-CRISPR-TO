%this code is to load a tiff video and the mask, it will calculate the intensity in the mask region. Sa Cai
% Clear workspace and close all figures
clear;
close all;

% Prompt user to select a TIFF video
[tifFileName, tifPathName] = uigetfile('*.tif', 'Select TIFF Video');
if isequal(tifFileName, 0)
    disp('No file selected for TIFF video. Exiting...');
    return;
end

% Read the TIFF stack information
info = imfinfo(fullfile(tifPathName, tifFileName));
nFrames = numel(info);

% Preallocate a matrix to hold frame data
tifData = zeros(info(1).Height, info(1).Width, nFrames, 'uint16'); % Change uint16 as needed based on your image type

% Read the TIFF file
for k = 1:nFrames
    tifData(:, :, k) = imread(fullfile(tifPathName, tifFileName), k);
end

% Prompt user to select a MATLAB mask file
[maskFileName, maskPathName] = uigetfile('*.mat', 'Select Mask File');
if isequal(maskFileName, 0)
    disp('No file selected for mask. Exiting...');
    return;
end

% Load the mask
maskData = load(fullfile(maskPathName, maskFileName)); % Ensure it loads a variable
fields = fieldnames(maskData);
mask = maskData.(fields{1}); % Load the first variable in the .mat file

% Ensure that the mask is binary and the correct size
if size(mask, 1) ~= info(1).Height || size(mask, 2) ~= info(1).Width || size(mask, 3) ~= nFrames
    disp('Mask dimensions do not match the TIFF video dimensions. Exiting...');
    return;
end

% Initialize array to store mean intensities
meanIntensities = zeros(nFrames, 1);

% Process each frame using its corresponding mask
for k = 1:nFrames
    frameData = tifData(:, :, k);
    
    % Extract the mask for the current frame
    currentMask = mask(:, :, k);
    
    % Convert the mask to binary if it isn't already
    currentMask = currentMask > 0; 

    % Calculate mean intensity within the masked area
    meanIntensity = mean(frameData(currentMask), 'omitnan'); % Omit NaN values in case of no overlap
    meanIntensities(k) = meanIntensity; % Store the mean intensity

    % Optionally display the masked frame
    %figure;
    %imshow(frameData, []);
    %hold on;
    %visboundaries(currentMask, 'Color', 'r'); % Visualize the mask
    %title(['Masked Frame ', num2str(k)]);
    %hold off;

   % pause(0.01); % Pause to view the masked frame result before moving to the next frame
   % close; % Close the figure
end

% Plot mean intensities over time
figure;
plot(meanIntensities, 'b-o');
xlabel('Frame Number');
ylabel('Mean Intensity');
title('Mean Intensity of Masked Area Over Time');
grid on;

disp('Processing complete. Mean intensities calculated and plotted.');
% Save results if desired
writetable(array2table(meanIntensities), fullfile(tifPathName, 'mean_intensities_1um.csv'));
disp(['Mean intensities saved as: ', fullfile(tifPathName, 'mean_intensities.csv')]);