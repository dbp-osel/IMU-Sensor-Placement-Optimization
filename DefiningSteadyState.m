function [Range] = DefiningSteadyState(Ang,instances)

%This function prompts the user to select ranges of steady state walking
%spanning between minima points. The user will first be asked to define how
%many regions of steady state walking there are. For instance, if a subject
%walks across a room, turns around, and walks back, there would be two
%instances of steady state gait corresponding to each pass of the room. The
%user will then be asked to define the ranges of each steady state region.

% Inputs
% Ang: The drift corrected angle calculated by the Angle_Corrected
%   function.
% instances: The zero velocity update points found by the Angle_Corrected
%   function.

% Outputs
% Range: A range of indices when steady state walking occurs.

close all
correct = 0;

figure
hold on
plot(1:length(Ang(:,1)),(Ang))
plot(instances,zeros(length(instances),1),'bo')
text(instances,zeros(length(instances),1),num2str([1:length(instances)]'))
hold off
title('Rotational Angle')
legend('rotx','roty','rotz')

while correct == 0;
numRegions = input('How many steady state regions are there? (enter an integer)');
Range = [];
for qq = 1:numRegions
    disp(strcat('Define region',32,num2str(qq),32,'as a 1x2 vector corresponding to the zero angle instances (blue dots)'))
Ranges = input(strcat('Region',32,num2str(qq)));
Range = [Range,instances(Ranges(1)):instances(Ranges(2))];
end

figure
plot(1:length(Ang(Range,1)),(Ang(Range,:)))
title('Rotational Angle')
legend('rotx','roty','rotz')

correct = input('Is this correct? (1/0)');
end
end



