% Read the CSV file into a table
filename = 'C1-Pbody1_time stamp_spots_alltraj.csv';
data = readtable(filename);
%%
data = data(2:end,:); %remove first line that contains only NaNs

%%
% Initialize containers for trajectories
%uniqueTracks = unique(data.TrackID);
uniqueTracks = unique(data.TRACK_ID);
numTracks = length(uniqueTracks);
%%
% Create cell arrays to store X, Y positions, and Distance_to_ROI for each track
xPositions = cell(numTracks, 1);
yPositions = cell(numTracks, 1);
distanceToROI = cell(numTracks, 1);
for i = 1:numTracks
  thisTrack = uniqueTracks(i);
  % Extract rows for the current track
 % subData = data(data.TrackID == thisTrack, :);
  subData = data(data.TRACK_ID == thisTrack, :);
% Sort by frame T so that positions are in chronological order
  % subData = sortrows(subData, 'T');
  subData = sortrows(subData, 'FRAME'); 
% Store X, Y, and Distance_to_ROI in cell arrays
% Make sure these column names match the ones in your CSV file
   xPositions{i}    = subData.POSITION_X;
   yPositions{i}    = subData.POSITION_Y;
  % distanceToROI{i} = subData.Distance_to_ROI;
end 

%%
% Create an array to store distance change for each track
%distanceChange = NaN(numTracks, 1);

for i = 1:numTracks
    if ~isempty(distanceToROI{i})
        distanceChange(i) = distanceToROI{i}(end) - distanceToROI{i}(1);
    else
        distanceChange(i) = NaN; % Handle empty cases
    end
end
%%
% Find the maximum trajectory length to pad shorter trajectories
maxLen = max(cellfun(@length, xPositions));

% Convert cell arrays to matrices, padding with NaNs for uneven lengths
xMatrix = NaN(numTracks, maxLen);
yMatrix = NaN(numTracks, maxLen);
%distanceMatrix = NaN(numTracks, maxLen);

for i = 1:numTracks
    len = length(xPositions{i});
    xMatrix(i, 1:len) = xPositions{i};
    yMatrix(i, 1:len) = yPositions{i};
  %  distanceMatrix(i, 1:len) = distanceToROI{i};
end

Dfin=xMatrix';
Dfiny=yMatrix';
%Distance_traj=distanceMatrix';
% % Output the X and Y position matrices
% disp('X Positions:');
% disp(xMatrix);
% 
% disp('Y Positions:');
% disp(yMatrix);