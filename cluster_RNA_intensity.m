% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code selected ROI of tiff files and then convert it to binary
% images, calculate the areas of different clusters 
% Author: Yanyu Zhu, Date: 6/7/2025
% Output: All the cluster areas distribution as well as the distribution of
% clusters within a certain area range
% clusterSizes and cluster_filter are the two variables needed for next
% steps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
[filename_RNA1,filepath]=uigetfile({'*.tif', 'Tiff File (*.tif)'}, 'Select a file');  % select tiff image file to get cell outline
cd(filepath); 

img1 = imread(filename_RNA1,2);  % Read the TIFF images, RNA channel
% Display the image and allow the user to draw multiple free-hand ROIs
figure;
imshow(img1, []);
caxis([100,500]);
title('Draw multiple ROIs on the Image. Right-click to finish each ROI.');

% Initialize an empty mask for all ROIs
allROIMask = false(size(img1));

% Allow user to draw multiple free-hand ROIs
while true
    roi = drawfreehand();
    if isempty(roi.Position) % If ROI is empty, break the loop
        break;
    end
    % Add the currently drawn ROI mask to the cumulative mask
    allROIMask = allROIMask | createMask(roi);
    
    % Ask the user if they want to draw another ROI
    choice = questdlg('Do you want to draw another ROI?', ...
        'Multiple ROIs', ...
        'Yes', 'No', 'No');
    if strcmp(choice, 'No')
        break;
    end
end

% Process only the pixels within the ROIs
% Calculate the area of the drawn ROI(s)
roiStats = regionprops(allROIMask, 'Area'); % Get the area of all ROIs
totalArea = sum([roiStats.Area]); % Sum the areas to get the total area

%% 
% this section using Xiang Wu's code to determine the threshold for
% binarization as another option. However, we use adpative method here. 
% Img_RNA_Vec_ROI = double(img1(allROIMask));  % Only use pixels within the ROI
% 
% [freq_RNA, Pix_Int_RNA] = ksdensity (Img_RNA_Vec_ROI,1:500);
% 
% RNA_Cas_thres_factor = 1;   
%    % figure(222)
%    % hold on;
%    % histogram(Current_Img_RNA_Vec,'Normalization','probability','BinLimits',[1,500]);
% 
%     % hold on;
%     % plot(Pix_Int_RNA,freq_RNA,'LineWidth',2);
%     % hold on;
% 
%             % Use the 2 most abundant pixel value as the background (weighted average)
%             [max_freq_RNA_a, ind_RNA_a] = max(freq_RNA);
% 
%             freq_RNA(ind_RNA_a) = 0;
% 
%             [max_freq_RNA_b, ind_RNA_b] = max(freq_RNA);
% 
%             background_RNA_int = (max_freq_RNA_a^2 * Pix_Int_RNA(ind_RNA_a) + max_freq_RNA_b^2 * Pix_Int_RNA(ind_RNA_b)) / (max_freq_RNA_a^2 + max_freq_RNA_b^2);
%             background_RNA = background_RNA_int+ std(Img_RNA_Vec_ROI) * RNA_Cas_thres_factor;        
% 
% background_RNA_normal=background_RNA/2^16;
   %%
        img1_binary = imbinarize(img1,"adaptive",'Sensitivity',0.4); % option 1 using adaptive 
       %  img1_binary = imbinarize(img1,background_RNA_normal); % option 2 using Xiang's code 
       %   img1_binary = imbinarize(img1); % option 3 using Otsu's threshold
       img1_binary(~allROIMask) = 0;  % Set pixels outside the ROI to 0

%%
% Create a figure and visualize the matrix
figure;
imagesc(img1_binary);

% Set the colormap to grayscale to clearly distinguish between 0 and 1
colormap(gray);
% Add a color bar
colorbar;
% Optionally, set axis properties for better visualization
axis equal tight;
%set(gca, 'XTick', 1:size(img1_binary, 2), 'YTick', 1:size(img1_binary, 1));
title('Binary Images');
xlabel('X(pixel)');
ylabel('Y(pixel)');

%%

pixel=0.1625 ; % pixel size in um 
% Label each connected component in the binary image
labeledImage = bwlabel(img1_binary);

% Measure the area of each connected component
clusterProperties = regionprops(labeledImage, 'Area','Centroid');

% Extract the area of each cluster into a vector
clusterSizes = [clusterProperties.Area]*pixel*pixel;

% Sort clusters by area and get the largest 30 clusters, can be revised by
% the users
[~, sorted_indices] = sort(clusterSizes, 'descend');
largest_clusters_indices = sorted_indices(1:30);

% Visualize the labeled image for reference
figure;
imshow(label2rgb(labeledImage, @jet, 'k', 'shuffle'));
title('Labeled Clusters');

hold on;
for i = 1:numel(largest_clusters_indices)
    cluster_index = largest_clusters_indices(i);
    centroid = clusterProperties(cluster_index).Centroid;
    area = clusterProperties(cluster_index).Area*pixel*pixel;
    text(centroid(1), centroid(2), sprintf('%0.4g', area), 'Color', 'yellow', 'FontSize', 12, 'FontWeight', 'bold');
end
hold off;




figure(1112);
hold on;
histogram (clusterSizes,'Normalization', 'probability','BinLimits',[1,50],'BinWidth',1);
title ('Normalized distribution of areas (um2) of all the clusters')
xlabel('area(um2)');


cluster_filter=clusterSizes(clusterSizes>50 & clusterSizes < 1200); % min cluster area 50 um2 can be revised
cluster_number = size(cluster_filter,2);
localization_index = 10^6 *cluster_number/totalArea;
figure(1113);
hold on;
histogram (cluster_filter,'Normalization', 'probability','BinWidth',5);
title ('Normalized distribution of areas (um2) of filtered clusters within area range ')
xlabel('area(um2)');