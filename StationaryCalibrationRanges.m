function [RangeStand,RangeSit] = StationaryCalibrationRanges(Acc)
%The user is prompted to select a range on the plot in during which the
%subject is standing upright and still. This range is used for a few things
%including calculations of the SI vectors based off of the direction of
%gravity.

%Input
%Acc: The acceleration signal from a sensor (typically on the foot)
%   type: nx3 double

%Outputs
%RangeStand: A range of indeces when the subject stands
%   type: 1x2 double
%RangeSit: A range of indeces when the subject sits
%   type: 1x2 double

figure
plot(1:length(Acc(:,1)),Acc)
title('Acceleration')
legend('accx','accy','accz')

clc
disp('Click on the plot: Choose Standing Calibration Start Point')
[Gravity_start,~] = ginput(1);
clc
disp('Click on the plot: Choose Standing Calibration End Point')
[Gravity_end,~] = ginput(1);
clc

RangeStand = [round(Gravity_start) round(Gravity_end)];

clc
disp('Click on the plot: Choose Sitting Calibration Start Point')
disp('If only performing a functional calibration just hit enter')
[Gravity_start,~] = ginput(1);
clc
disp('Click on the plot: Choose Sitting Calibration End Point')
disp('If only performing a functional calibration just hit enter')
[Gravity_end,~] = ginput(1);
clc

RangeSit = [round(Gravity_start) round(Gravity_end)];
end



