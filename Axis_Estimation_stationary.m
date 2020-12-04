function [SI,ML,AP] = Axis_Estimation_stationary(Acc,RangeStand,RangeSit)
% This function uses the gravity vector during two static poses to define
% the coordinate systems of the body segments. During this calibration the
% subject begins standing stationary, then sits with their legs extended
% while leaning back. During the seated pose care is taken so that the long
% axes of the body segments lie in the sagittal plane. This function takes
% the acceleration signal (Acc) as well as ranges of indices of still
% standing (RangeStand) and of sitting (RangeSit). The components of the
% gravity vector are averaged for the standing and sitting time ranges. The
% standing gravity vector is normalized and taken as the SI vector. The
% sitting vector is normalized and a cross product between this vector and
% the SI vector gives the ML vector. The cross of SI and ML gives AP. The
% cross of AP and SI then redefines ML to ensure orthogonality. 

% Inputs
% Acc: The acceleration signal.
% RangeStand: A range of indices when the subject stands still.
% RangeSit: A range of indices when the subject sits in the specific
%   posture mentioned above.

% Outputs
% SI: A unit vector aligned with the SI direction (superior is positive)
% ML: A unit vector aligned with the ML direction (lateral is positive)
% AP: A unit vector aligned with the AP direction (anterior is positive)

spanG = length(Acc(RangeStand(1):RangeStand(2),1));
GravVect = [sum(Acc(RangeStand(1):RangeStand(2),1))/spanG,sum(Acc(RangeStand(1):RangeStand(2),2))/spanG,sum(Acc(RangeStand(1):RangeStand(2),3))/spanG];
SI = GravVect./norm(GravVect);

spanC = length(Acc(RangeSit(1):RangeSit(2),1));
CoplanerVect = [sum(Acc(RangeSit(1):RangeSit(2),1))/spanC,sum(Acc(RangeSit(1):RangeSit(2),2))/spanC,sum(Acc(RangeSit(1):RangeSit(2),3))/spanC];
NormCoplanerVect = CoplanerVect./norm(CoplanerVect);

ML = cross(NormCoplanerVect,SI);
ML = ML./norm(ML);

AP = cross(SI,ML);
AP = AP./norm(AP);

ML = cross(AP,SI);
ML = ML./norm(ML);
end

