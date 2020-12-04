function [RotAng,instances] = Angle_Corrected(HeelRot,instances,slopeSpan)
%% Creates a corrected angle plot for the heel sensor
%This function is used to calculate drift corrected angle from the
%rotational velocity output of the heel sensor's gyroscope (note that any
%sensor can be used with this function, the heel sensor just produces the
%best waveform). This is an approximation which will be used primarily
%to determine ranges of steady state gait.
 
% The leg is approximately vertical during the times when the acceleration
% and rotational velocity is minimal. From this we can assume that all of
% the joint angles at these times are approximately zero and apply a linear
% correction for drift between these zero points. The results of this
% integration are valid during steady state gait where the subject isn't
% turning at the end of the mat (The large spikes that occur at each turn
% are due to the fact that the angle is increasing by 360 degrees around
% vertical and are completely invalid after drift correction. The
% assumptions made for the drift correction are only valid between
% consecutive minimum points during steady state gait.)

% Inputs
% HeelRot: The three axis rotational velocity signal
% Instances: Indices where the rotational velocity is approximately zero
%   corresponding to flat foot.
% slopeSpan: The number of frames above and below the flat foot points
%   (Instances) for which a refinement will be performed. The flat foot
%   points (Instances) are used as an initial guess, then an algorithm
%   looks for a point of minimal slope in the primary component of angle.
%   This minimum is the new zero velocity update point and the angle is 
%   recalculated.

% Outputs
% RotAng: The drift corrected angle of the heel sensor. This is found by
%   integrating the rotational velocity signal, then
%   using zero velocity update points to determine where the angle should
%   be zero. Drift is assumed linear and subtracted from subsequent zero
%   velocity update points.
% instances: The zero velocity update points. These are refined from the
%   initial guess given by the input instances.
%% Initial integration using approximated minimum motion points
%The rotational velocity between each minimum point (instances) is
%integrated to yield uncorrected angle. The drift (assumed linear) is then
%subtracted out so that the first and last points are zero. The drift is
%represented as a linear equation connecting the previous corrected point
%(who's value is zero) to the next uncorrected point (who's value is not
%zero due to drift). When this drift equation is subtracted out the
%uncorrected point becomes zero (as it should be) and all of the points
%between these zero points are shifted proportionally.

close all
InitialInstances = instances;
% while correct == 0;
instances = InitialInstances;
%Angle
RotAng = zeros(instances(1),3);
for i = 1:length(instances)-1
    vel = HeelRot(instances(i):instances(i+1),:);
    span = instances(i+1)-instances(i);
    time = (instances(i):instances(i+1))';

    angX = cumtrapz(vel(:,1));
    angY = cumtrapz(vel(:,2));
    angZ = cumtrapz(vel(:,3));
    
    %Linear Drift is Assumed
    slopeX = (angX(end)-angX(1))/span;
    driftX = slopeX*(time-time(1));
    angX = angX-driftX;
    
    slopeY = (angY(end)-angY(1))/span;
    driftY = slopeY*(time-time(1));
    angY = angY-driftY;
    
    slopeZ = (angZ(end)-angZ(1))/span;
    driftZ = slopeZ*(time-time(1));
    angZ = angZ-driftZ;
    
    tempVel = [angX,angY,angZ];
    RotAng = [RotAng;tempVel(2:end,:);];
    

end
RotAng = [RotAng;zeros(length(HeelRot)-length(RotAng),3);];


figure
hold on
plot(1:length(RotAng(:,1)),RotAng(:,1),'r')
plot(1:length(RotAng(:,1)),RotAng(:,2),'g')
plot(1:length(RotAng(:,1)),RotAng(:,3),'b')
plot(instances,zeros(length(instances),1),'bo')
text(instances,zeros(length(instances),1),num2str([1:length(instances)]'))
legend('component 1','component 2','component 3','zero angle instances')
hold off

%% Refinement of the initial minimal motion points.
% The user is prompted to select the primary component of angle based on the
% three sensor axes. For this to work, one of the sensor axes should be
% fairly in line with the primary axis of rotation, however this is only an
% approximation so it does not have to be perfect.

% The zero angle approximation is then further refined. The variable
% 'slopeSpan' is used to define a range of indexes above and below the
% initial minimal motion points to look for instances of minimal slope for
% the primary component of angle. The waveform should be flattest for the
% primary component of rotation during the middle of flatfoot.

comp = input('What is the primary component of rotation(1, 2 , or 3)?');
big = abs(RotAng) > ones(length(RotAng(:,1)),3)*10;
smallRotAng = RotAng;
smallRotAng(big) = NaN;
%smallRotAng contains all of the values of angle less than 10 degrees. This
%variable is used when looking at slope to eliminate the possibility of
%mistakenly finding instances of minimal slope in the peaks and valleys of
%the angle plot.
RotAng_Slope = zeros(length(smallRotAng(:,1)),3);

for i = 1:length(smallRotAng(:,1))-1
    RotAng_Slope(i,:) = smallRotAng(i+1,:)-smallRotAng(i,:);
end

for j = 1:length(instances)-1

        instances_temp = instances(j)-slopeSpan + ...
            (find(abs(RotAng_Slope(instances(j)-slopeSpan:instances(j)+slopeSpan,comp)) ...
            == min(abs(RotAng_Slope(instances(j)-slopeSpan:instances(j)+slopeSpan,comp)))));
        instances(j) = round(instances_temp(1));
end

%% Re-integration of rotational velocity with updated minima points
RotAng = zeros(instances(1),3);
for i = 1:length(instances)-1
    vel = HeelRot(instances(i):instances(i+1),:);
    span = instances(i+1)-instances(i);
    time = (instances(i):instances(i+1))';

    angX = cumtrapz(vel(:,1));
    angY = cumtrapz(vel(:,2));
    angZ = cumtrapz(vel(:,3));
    
    %Linear Drift is Assumed
    slopeX = (angX(end)-angX(1))/span;
    driftX = slopeX*(time-time(1));
    angX = angX-driftX;
    
    slopeY = (angY(end)-angY(1))/span;
    driftY = slopeY*(time-time(1));
    angY = angY-driftY;
    
    slopeZ = (angZ(end)-angZ(1))/span;
    driftZ = slopeZ*(time-time(1));
    angZ = angZ-driftZ;
    
    tempVel = [angX,angY,angZ];
    RotAng = [RotAng;tempVel(2:end,:);];
    

end
RotAng = [RotAng;zeros(length(HeelRot)-length(RotAng),3);];

figure
hold on
plot(1:length(RotAng(:,1)),RotAng(:,1),'r')
plot(1:length(RotAng(:,1)),RotAng(:,2),'g')
plot(1:length(RotAng(:,1)),RotAng(:,3),'b')
plot(instances,zeros(length(instances),1),'ro')
hold off
end






