function [ rotated ] = QuaternionConjugation(q,p)
% Conjugation of p with q (p => qpq^-1). 
% This function rotates the quaternion p using the rotational quaternion q
% through conjugation.

% Inputs
% p: quaternion to be rotated
% q: the quaternion that rotates p

% Outputs
% rotated: The rotated quaternion
%% Basic notes on quaternion math
 
% ***What is a quaternion and why is it useful?***
% A quaternion is like a 4D vector with one scaler, and three imaginary
% parts.
%       q = [w + i + j + k]
% Recall that the square root of -1 is i and that i*i = -1, i*i*i = -i,
% i*i*i*i = 1, i*i*i*i*i = i, etc making complex numbers of the form a + bi
% useful for representing rotation. cos(x) + i*sin(x) represents rotation
% by the angle x around a unit circle in the real/imaginary plane (recall
% Euler's Formula). As far as I can tell this is essentially that in 3D.
% Hamilton (the guy who figured this out) just introduced two more
% imaginary numbers (j and k) to represent two more axes of rotation.
 
% ***Inverse quaternions***
% The inverse of a quaternion can be found by taking the negative of the
% imaginary components while keeping the scaler portion the same.
%       say q = [w+i+j+k] where w is scaler and i,j, and k are imaginary
%       The inverse of q is q^-1 = [w-i-j-k]
 
% ***Pure quaternions***
% A pure quaternion is a quaternion that has no scaler portion. Pure
% quaternions can be treated like vectors. In this instance p is a pure
% quaternion and represents one of the three unit vectors corresponding to
% the coordinate axes of the sensor where i, j, and k lie on the axes of
% the sensor (ie p=[0w+1i+0j+0k], p=[0w+0i+1j+0k], or p=[0w+0i+0j+1k]). In
% this form the imaginary portions (ijk) can be treated as a unit vectors.
% The output rotated is also a pure quaternion who's ijk components
% represent unit vectors pointing in the direction of the sensor axis.
 
% ***Quaternion Conjugation***
% p => pqp^-1 The conjugate of a quaternion can be calculated by
% multiplying q into p and then by multiplying the result (qp) into the
% inverse of q (q^-1). To rotate a unit vector (p = [ai + bj + ck]) around
% an arbitrary axis ([xi + yj + zk]) by an angle (a) we first calculate q as
%       q = cos(a/2) + [xi + yi + zk]*sin(a/2)
% and q^-1 as
%       q^-1 = cos(a/2) - [xi + yj + zk]*sin(a/2)
% Note that this formula is very similar to Euler's formula 
% e^ia = cos(a) + i*sin(a). This can sort of be thought as using Euler's 
% function in 3D space. Also notice that the angle in our trig function is
% divided by two. This is because multiplying the unit vector (p) by q and
% q^-1 actually rotates the vector twice so it is necessary to halve the
% desired angle in q and q^-1 to get the desired result. The reason the
% rotation is performed this way is because it cancels out the change in
% the scaler portion so that the resultant quaternion remains a unit vector
% (ie has a zero scaler value). The final calculation looks something like
% this
%   p => [cos(a/2)+[xi+yi+zk]*sin(a/2)]*[ai+bj+ck]*[cos(a/2)-[xi+yj+zk]*sin(a/2)]
% The Xsens IMU sensors automatically calculate the four components of the
% rotational quaternion q, so all we have to do is first find the inverse
% by changing the sign of the last three components (imaginary components)
% then perform the calculation qpq^-1 using basic quaternion
% multiplications which can be done with a quaternion multiplication table
% (see QuaternionMultiplication.m for notes on that). 

%% Calculations
q_inv = [q(:,1),-1*q(:,2),-1*q(:,3),-1*q(:,4)]; %Xsens does [i j k scaler], not [scaler i j k]
%The inverse of q is found by taking the negative of the last three imaginary components
rotated = QuaternionMultiplication(QuaternionMultiplication(q,p),q_inv);
%The rotation is performed by taking qpq^-1 which is done first by
%multiplying q into p, then the result into q^-1. 

end

