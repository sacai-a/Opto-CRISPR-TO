 %this code is to find the max intensity point for each frame, have a given radiuc and save the area as a mask    
%user will pick a retangular region to confine the center region Sa 
% Clear workspace and close all figures 
clear;
close all;

% Prompt user to select the TIFF file
[fileName1, pathName] = uigetfile('*.tif', 'Select the first channel TIFF video');
if isequal(fileName1, 0)
    disp('No file selected for the first channel. Exiting...');
    return;
end

% Read the TIFF stack information for the first channel
info1 = imfinfo(fullfile(pathName, fileName1));
nFrames1 = numel(info1);


% Determine the width and height of the image
imgWidth = info1(1).Width;
imgHeight = info1(1).Height;

% Preallocate matrices for channels
channel1 = zeros(imgHeight, imgWidth, nFrames1, 'uint16'); % First channel


% Read the TIFF files frame by frame for both channels
for k = 1:nFrames1
    channel1(:, :, k) = imread(fullfile(pathName, fileName1), k);
    
end

% Ask user to draw an ROI on the first frame of the first channel
figure;
imshow(channel1(:, :, 1), []);
title('Draw a ROI on the first frame of the first channel');
h = drawrectangle('Color', 'r'); % Allow the user to draw a rectangle
position = h.Position; % Get the position of the rectangle (x, y, width, height)

% Release the rectangle object after use
wait(h);
roiMask = createMask(h); % Create a binary mask from the ROI

% Close the ROI figure
close(gcf);

% Initialize arrays for picked points and intensities
pickedPoints = zeros(nFrames1, 2); % To store (x, y) points for max values
intensities = zeros(nFrames1, 1); 
masks = false(imgHeight, imgWidth, nFrames1); % To store masks for each frame

% Process each frame to find the maximum pixel in the ROI
for k = 1:nFrames1
    % Get the current frame from channel 1
    currentFrame = channel1(:, :, k);
    
    % Apply ROI mask to isolate the ROI
    roiValues = currentFrame(roiMask); % Extract values in the ROI
    
    % Find the maximum value and its index in the ROI
    [maxValue, maxIndex] = max(roiValues);
    
    roiHeight = round(position(4)); % Get the height of the ROI
    roiWidth = round(position(3)); % Get the width of the ROI
    
    % Use the dimensions of the ROI instead of the mask
    [rowROI, colROI] = ind2sub([roiHeight, roiWidth], maxIndex); % Use ROI size

    % Calculate the actual (x, y) coordinates in the original image
    x = colROI + round(position(1)); % Adding the x offset of the ROI
    y = rowROI + round(position(2)); % Adding the y offset of the ROI

    % Store the picked points (x, y) for the maximum value
    pickedPoints(k, :) = [x, y];

    % Ensure x, y are within the image bounds before measuring intensity
    if x < 1 || x > imgWidth || y < 1 || y > imgHeight
        disp(['Calculated coordinates (', num2str(x), ', ', num2str(y), ') are out of bounds; skipping frame ', num2str(k)]);
        continue; % Skip to the next frame if the coordinates are out of bounds
    end
    
    % Define a default radius for measurement
    radius_um = 2; % Radius in micrometers
    pixelSize = 0.1; % Updating this value according to your pixel size (in um/pixel)
    radius_pix = round(radius_um / pixelSize); % Convert radius to pixels

    % Create a mask for the circle
    [X, Y] = meshgrid(1:imgWidth, 1:imgHeight);
    mask = (X - x).^2 + (Y - y).^2 <= radius_pix^2;

    % Store the mask for the current frame
    masks(:, :, k) = mask; % Save mask for later use

end

% Save picked points to a MAT file for later use
save(fullfile(pathName, 'picked_points.mat'), 'pickedPoints');
save(fullfile(pathName, 'masks_2um'), 'masks');

