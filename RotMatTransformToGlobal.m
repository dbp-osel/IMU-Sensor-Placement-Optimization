function [newXax,newYax,newZax] = RotMatTransformToGlobal(u,v,w,Xax,Yax,Zax)
% This function constructs a direction cosine matrix (DCM) constructed from
% the unit vectors u, v, and w and uses this DCM (A) to rotate three
% vectors corresponding to one basis to the "global coordinate system" they
% are represented in terms of.

% The unit vectors u, v, and w are taken from the beginning of each trial
% when the subject is standing still. These vectors are used to create a
% rotation matrix which rotates the coordinate system of a given body
% segment at all time points such that these axes initially line up with a
% global coordinate system (ie [1 0 0],[0 1 0],[0 0 1]). Since the
% coordinate systems of each body segment are more or less parallel at the
% start of each trial, we can assume that all of the sensors are in terms
% of the same global coordinate system after this transformation. In other
% words, since all of the axes of each sensor should be lined up at the
% start of each trial, we can simply rotate the coordinate systems so they
% do indeed all overlap at the beginning of the trial.

A = [u(1) u(2) u(3);
     v(1) v(2) v(3);
     w(1) w(2) w(3);];

span = length(Xax(1,:));
    newXax = zeros(3,span);
    newYax = newXax;
    newZax = newXax;
    
for ii = 1:span
    newXax(:,ii) = A*Xax(:,ii);
    newYax(:,ii) = A*Yax(:,ii);
    newZax(:,ii) = A*Zax(:,ii);
end

end

