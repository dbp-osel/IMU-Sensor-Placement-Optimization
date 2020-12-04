function [SI] = Axis_Estimation_SI_Only(Acc,GravSpan)
%% Finding Superior Inferior Axis Through Gravity
%At the start when the subject is standing still the main contributor to
%acceleration should be gravity. The acceleration vectors from each sensor
%are averaged during this time period to find the gravity vector associated
%with each sensor.

% Inputs
% Acc: Acceleration signal
% GravSpan: A range of indices when the subject stands still

% Output
% SI: A unit vector aligned with the SI direction (Superior is positive)

Gravity_start = GravSpan(1);
Gravity_end = GravSpan(2);

Gravity = mean(Acc(round(Gravity_start):round(Gravity_end),:),1);
SI = Gravity/norm(Gravity);
end

