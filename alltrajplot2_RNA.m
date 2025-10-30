

% Description: Plots2D particle trajectories for all cell in ascending 
% order of length with optimized sorting or of linearity 
% Author: Yanyu Zhu, Date: 5-1-2025
% after running alltraj
% Dfin = x-coordinate
% Dfiny = y-coordinate

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
tic

  % reading the Dfin and Dfiny Matrices from the workspace   
        Dfin=evalin('base','Dfin');
        Dfiny=evalin('base','Dfiny');
        Dfiny=-Dfiny;
        % Dfin=evalin('base','Dfin_new(:,821:845)');   % use this when you
        % % plot first or last 50 trajectories 
        % Dfiny=evalin('base','Dfiny_new(:,821:845)');
        
        lr = size(Dfin,1);
     
        lc = size(Dfin,2);
        
        cc=hsv(size(Dfin,2));
        
        for i = 1:lc
            
            x = Dfin(1,i);
            y = Dfiny(1,i);

            x = nanmean(Dfin(:,i));
            y = nanmean(Dfiny(:,i));
            
          for j = 1:lr

            Dfin(j,i) = Dfin(j,i)-x;           % shifting all points by first coordinate
            Dfiny(j,i) = Dfiny(j,i)-y;         % staring from zero
            
            
            
          end
        end
        
        assignin('base','Dfinshifted1',Dfin);
        assignin('base','Dfinyshifted1',Dfiny);
 
        % 
        % for i = 1:lc                             % rsquare excurtion
        % 
        %  % r(i)=sqrt((max(Dfin(1:lr,i))-min(Dfin(1:lr,i)))^2+(max(Dfiny(1:lr,i))...
        %  %     -min(Dfiny(1:lr,i)))^2);
        % 
        %  r(i)=(sqrt((max(Dfin(1:lr,i))-min(Dfin(1:lr,i)))^2+(max(Dfiny(1:lr,i))...
        %      -min(Dfiny(1:lr,i)))^2))/sqrt(sum(~isnan(Dfin(1:lr,i))));
        % 
        %     % for j= 1: lr
        %     %    Rg(j,i)=(Dfin(j,i)-nanmean(Dfin(1:lr,i)))^2+(Dfiny(j,i)-nanmean(Dfiny(1:lr,i)))^2
        %     %       % sorting by radius of gyration;
        %     % end 
        %     % r(i)=sqrt(nansum(Rg(1:lr,i)));   %  /sqrt(sum(~isnan(Dfin(1:lr,i))));
        % end
        
%%
    num_trajectories = size(Dfin, 2);
    mean_linearities = zeros(1, num_trajectories);
    
   % Iterate through each trajectory
    for j = 1:num_trajectories
        x = Dfin(:, j);
        y = Dfiny(:, j);
        
        % Find valid (non-NaN) indices
        valid_indices = ~isnan(x) & ~isnan(y);
        
        if sum(valid_indices) < 2
            % If fewer than 2 valid points, mean linearity is undefined
            mean_linearities(j) = NaN; 
            fprintf('Mean Linearity for Trajectory %d: Undefined (not enough points)\n', j);
            continue;
        end
        
        % Find the start and end points (first and last valid points)
        start_point = [x(find(valid_indices, 1, 'first')), y(find(valid_indices, 1, 'first'))];
        end_point = [x(find(valid_indices, 1, 'last')), y(find(valid_indices, 1, 'last'))];
        
        % Calculate straight line length from start to end
        straight_line_length = norm(end_point - start_point);
        
        % Calculate actual trajectory length using valid points
        trajectory_length = 0;
        for i = 1:length(x)-1
            if valid_indices(i) && valid_indices(i+1)
                segment_length = norm([x(i+1) - x(i), y(i+1) - y(i)]);
                trajectory_length = trajectory_length + segment_length;
            end
        end
        
        % Calculate mean linearity
        if trajectory_length > 0
            mean_linearities(j) = straight_line_length / trajectory_length;
        else
            mean_linearities(j) = 0;  % If no trajectory length, set to zero to avoid division by zero
        end
       
    end

      r=mean_linearities;
%%
        assignin('base','r_excurtion',r);
        
        [rsort,rindex] = sort(r);
        
        assignin('base','r_ex_sorted',rsort);
        assignin('base','sorted_index',rindex);
        
       for i = 1:lc
           
           Dfin_new(1:lr,i) = Dfin(1:lr,rindex(i));
           Dfiny_new(1:lr,i) = Dfiny(1:lr,rindex(i));
         %  distance_new(1:lr,i) = Distance_traj(1:lr,rindex(i));
         %  distance_change_new(i)=distanceChange(rindex(i));
       end
        
        % assignin('base','Dfinnew_excur',Dfin_new);
        %  assignin('base','Dfinynew_excur',Dfiny_new);
    
   
     % for i=1:lc
     %       step(i)=sum(~isnan(Dfin_new(:,i)));
     %       velo(i)=sqrt((Dfin_new(step(i),i)-(Dfin_new(1,i)))^2+(Dfiny_new(step(i),i)-(Dfiny_new(1,i)))^2)/(3*(step(i)-1));
     % end
     % 
     % assignin('base','step',step);
     %    assignin('base','velo',velo);    

%%
     k=10;
     c=0;
     figure(8)
        for i = 1:lc
            
       %    plot(Dfin(1:lr,i)+c,Dfiny(1:lr,i)-k);
           plot(Dfin_new(1:lr,i)+c,Dfiny_new(1:lr,i)-k);
           
           c=c+5;
           if mod(i,5)==0
           k=k+5;
           c=0;
           end
           hold on
           axis image
          
           
        end
      % colorarray = rand(lc,3);
      % colororder(colorarray);
       
%daspect([1 1]);
xlabel('X (um)');
ylabel('Y (um)');
% xlim([40*pixl 86*pixl])
% ylim([-85*pixl -40*pixl ])
fontsize(20,"points");
%title('Locus 2');

toc 
%%
% for kk=1:size(distance_new,2)
%     nanlendistance_new(:,kk)
% 
% end
% 
% 
%     [numSteps, numTrajectories] = size(distance_new);
%     stepLengths = zeros(1, numTrajectories);
% 
%     % Loop through each trajectory
%     for i = 1:numTrajectories
%         % Find the indices of non-NaN values in the current trajectory
%         nonNaNIndices = find(~isnan(distance_new(:, i)));
% 
%         % Calculate the step length as the number of non-NaN values
%         stepLengths(i) = length(nonNaNIndices);
%     end

%%
       
     

            
            
           
        
