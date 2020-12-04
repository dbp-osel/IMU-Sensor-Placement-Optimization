function [components] = CoordinateTransformQuatPreSpecifiedAxes(Quats,axis1,axis2,axis3)
% This function rotates a unit vector lying along a given axis of the sensor
% as in order to determine its motion in space at every given time point.
% This unit vector is defined by its' components 'axis1', 'axis2', and
% 'axis3' and is rotated using quaternions inputted as the variable 'Quats'.
% 'Quats' in an nx4 matrix containing all of the quaternions so that each
% row corresponds to a single quaternion. The first component of each
% quaternion is the scaler component while the next three components are the
% i,j, and k components.

span = length(Quats(:,1));
zero_ = zeros(span,1);
one_ = ones(span,1);

comp1 = QuaternionConjugation(Quats,[zero_,one_*axis1(1),one_*axis1(2),one_*axis1(3)]);
comp2 = QuaternionConjugation(Quats,[zero_,one_*axis2(1),one_*axis2(2),one_*axis2(3)]);
comp3 = QuaternionConjugation(Quats,[zero_,one_*axis3(1),one_*axis3(2),one_*axis3(3)]);

components.Axis1 = comp1(:,[2 3 4])';
components.Axis2 = comp2(:,[2 3 4])';
components.Axis3 = comp3(:,[2 3 4])';


%All three components (comp1,comp2, and comp3) are put into a structure
%called components. Since these are unit quaternions the scaler portion
%corresponding to the first column was omitted (these are all zero). Only
%the vector portion is kept (these are essentially unit vectors).
end