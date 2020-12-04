function [instances,threshold] = ZeroVelocityFinder(TotalAcc,TotalAngVel,RangeStand,RangeSit)
% The points where the total acceleration and the total rotational velocity
% of the sensors are close to zero are found. It is assumed that this 
% occurs sometime during stance when the leg is close to vertical. These
% points will serve as an initial guess at when the joint angles of the leg
% are zero.

% Inputs
% TotalAcc: A vector containing the magnitude of the acceleration signal
% TotalAngVel: A vector containing the magnitude of the rotational
%   velocity signal
% RangeStand: A range of indices during which the subject is standing still
% RangeSit: A range of indices during which the subject is seated in the
%   calibration pose

% Outputs
% Instances: Indices of minimal motion where rotational velocity is close
%   to zero. These are treated as instances of flat foot during mid-stance.
% threshold: The number of frames separating indices of minimal motion
%   (Instances). Detected points will be treated as coincident if they are
%   within +/- the threshold to avoid duplicate minimal motion points for a
%   single flat foot detection.
correctVel = 0;
while correctVel ==0;
    %The total magnitude of rotational velocity is found by adding the
    %magnitude of rotational velocity of each sensor. The result of this may not
    %be particularly meaningful but is a convenient way to boost the signal
    %as the sensors which should be in agreement when there is no velocity
    %on the body.
    
    % The user is prompted to select a point on the graph above which there
    % are no valleys they want to scope out. Then the user is asked to enter
    % the minimum distance apart adjacent valleys should be. These inputs
    % are fed into Matlab's findpeaks function to get points corresponding
    % to the times when rotational velocity is minimal.
    figure
    hold on
    plot(1:length(TotalAngVel),TotalAngVel)
    title('Total Net Rotational Velocity of All Sensors')
    clc
    disp('This section uses the Matlab "findpeaks" function to determine')
    disp('the points corresponding to minimal rotational velocity')
    disp('Click on the plot: Maximum vally height')
    maxHeightVel = ginput(1);
    minDistVel = input('Enter in the command prompt: Minimum Valley Distance? (enter a number)');
    clc
    
    [pksR,locsR] = findpeaks(-1*TotalAngVel,'MinPeakDistance',minDistVel,'MinPeakHeight',-1*maxHeightVel(2));
    plot(locsR,-1*pksR,'ro')
    hold off
 
    correctVel = input('Continue? (click zero to try again, otherwise just hit enter)');
end
close all

correctAcc = 0;
while correctAcc ==0;
    %The points of minimal total acceleration are found in the same manner
    %as those for rotational velocity.

    figure
    hold on
    plot(1:length(TotalAcc),TotalAcc)
    title('Total Net Acceleration of All Sensors')
    clc
    disp('This section uses the Matlab "findpeaks" function to determine')
    disp('the points corresponding to minimal acceleration')
    disp('Click on the plot: Maximum valley height')
    maxHeightAcc = ginput(1);
    minDistAcc = input('Enter in the command prompt: Minimum Valley Distance? (enter a number)');
    clc
    
    [pksA,locsA] = findpeaks(-1*TotalAcc,'MinPeakDistance',minDistAcc,'MinPeakHeight',-1*maxHeightAcc(2));
    plot(locsA,-1*pksA,'ro')
    hold off
      
    correctAcc = input('Continue? (click zero to try again, otherwise just hit enter)');
   clc
end

correctTol = 0;
while correctTol ==0;
%The points found for minimal acceleration and minimal rotational velocity
%don't overlap perfectly. In this section the user is prompted to enter a
%threshold value indicating how close the points should be together to
%treat as coincidental.
    
figure
hold on
plot(locsA,zeros(length(locsA),1),'ro')
plot(locsR,zeros(length(locsR),1),'bo')
legend('Acceleration Minima','Rotation Minima')
plot(0,50)

disp('This plot shows the indexes where the minimal acceleration')
disp('and minimal rotational velocity points occur. The goal here')
disp('is to choose a threshold within which to treat these points')
disp('as occurring simultaneously. The goal is to look for a disperse')
disp('and evenly spread set of points. Multiple values may be tested')
disp('so experiment untill you get satisfactory results')
threshold = input('Enter in the command prompt: Threshold to treat minima points as a zero point? (enter a number)');
clc
instInd = 1;
for p = 1:length(locsA)
    for q = 1:length(locsR)
    if locsA(p) < locsR(q) + threshold && locsA(p) > locsR(q)-threshold
        instances(instInd) = mean(locsA(p),locsR(q));
        instInd = instInd+1;
    end
    end
end

for p = 1:length(instances)
    for k = 1:length(instances)
        if instances(p) == instances(k) && p~=k
            instances(k) = NaN;
        end  
    end
end
instances(find(isnan(instances))) = [];

if length(RangeSit) == 2
    deletions = instances < max(max(RangeStand),max(RangeSit));
else
    deletions = instances < max(max(RangeStand));
end
instances(deletions) = [];
%ensures that the minimal motion points only exist after the initial
%stationary calibration procedure

clear time driftX
     
plot(instances,ones(length(instances),1),'ko')
hold off

correctTol = input('Continue? (click zero to try again, otherwise just hit enter)');
clc
end