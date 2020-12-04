function [proj] = DCMProjection(vector,Axis1,Axis2,Axis3)
%This function projects a vector onto another coordinate system using
%DCMs. The inputs into this function must be unit vectors.

%The input 'vector' consists of unit vectors described by a 3xn matrix where
%each column represents a single unit vector. The inputs 'Axis1', 'Axis2', and
%'Axis3' are also unit vectors described by 3xn matrices. These inputs
%describe a coordinate system the vectors described by the input 'vector' are
%projected onto.

%This function is used to project a moving axis onto a coordinate system
%that is also assumed to be moving. This is used to describe the motion of
%moving body segments relative to other moving body segments. For instance,
%if a vector describing the long axis of the shin is projected into a
%coordinate system describing the thigh, then the transformed vector in
%the thigh coordinate system can be used to calculate knee angles.

proj = zeros(3,length(Axis1(1,:)));
for i = 1:length(Axis1(1,:))
     proj(:,i) = [Axis1(:,i)';Axis2(:,i)';Axis3(:,i)';]*vector(:,i); %DCM projection
end
norm = sqrt(sum(proj.^2,1));
proj = [proj(1,:)./norm;proj(2,:)./norm;proj(3,:)./norm;]; %re-normalization to play it safe
end

