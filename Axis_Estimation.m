function [SI,ML,AP] = Axis_Estimation(Acc,Rot,GravSpan,SteadySpan,guess)
%% This function uses the gravitational acceleration vector while the 
% subject is stationary to determine the SI direction and principal
% component analysis (PCA) of the rotational velocity to estimate the ML
% direction as the subject walks at steady state. The AP direction is the
% cross product between SI and ML vectors. After AP is found the AP and SI
% vector cross product re-defines ML to ensure the three vectors are
% mutually orthogonal.

% Inputs
% Acc: The acceleration signal
% Rot: The rotational velocity signal
% GravSpan: A span of indices when the subject is standing still in a
%   neutral posture.
% SteadySpan: A span of indices when the subject is walking straight at a
%   steady pace.
% guess: An initial guess of the direction of the ML vector. Depending on
%   the alignment of the sensor the ML estimate will either be the 1st pca
%   component, or the negative of the first pca component. Comparison with
%   an initial guess prevents sign errors.

% Outputs
% SI: A unit vector aligned with the SI direction (superior is positive)
% ML: A unit vector aligned with the ML direction (lateral is positive)
% AP: A unit vector aligned with the AP direction (anterior is positive)

%% Finding Superior Inferior Axis Through Gravity
%At the start when the subject is standing still the main contributor to
%acceleration should be gravity. The acceleration vectors from each sensor
%are averaged during this time period to find the gravity vector associated
%with each sensor. This time period is defined by the variable GravSpan
%which contains a 2 element vector with the start and stop indexes. Note
%that the force caused by gravity is actually upward since the force the
%sensor feels is really the resultant force of the body it is attached to
%on it (Remember Newton's third law?).

Gravity_start = GravSpan(1);
Gravity_end = GravSpan(2);

Gravity = mean(Acc(round(Gravity_start):round(Gravity_end),:),1);
SI = Gravity/norm(Gravity);

%% Finding Medial Lateral through Principal Component Analysis (PCA)
% The principal axis of rotation is assumed to be about the medial-lateral
% axis. The first component of the PCA algorithm is in the direction that
% describes the greatest possible variance of the data. This is done by
% fitting a line to the data that minimizes the distance between each data
% point and its' orthogonal projection onto this line. This is similar to
% linear regression, the difference being that linear regression minimizes
% the distance between the fitted line and the data points along the axis
% of some dependent variable while PCA finds a line that minimizes the
% distance in all possible directions.

% Before PCA analysis the data is run through a low pass moving average
% filter to eliminate possible outliers in the data.
Rotfilt = Rot;
filterWeights = [1 2 3 2 1]';
filterSpan = (length(filterWeights)-1)/2;
filterMean = sum(filterWeights);
for zz = filterSpan+1:length(Rot(:,1))-filterSpan
    Rotfilt(zz) = (Rot(zz-filterSpan:zz+filterSpan)*filterWeights)/filterMean;
end
Rot = Rotfilt;

% PCA is performed using Matlab's pca function
pcaCoef = pca(Rot(SteadySpan,:));

% Often the first PCA component will be aligned along the ML axis in the
% wrong direction. The first PCA component is compared to the initial guess
% given by the user to determine whether or not this vector points in the
% right direction (laterally rather than medially). For this to work the
% guess provided by the user must be within +/- 0.5 of the actual vector
% for each component. If the first PCA component is within +/- 0.5 for each
% component of the guess, then the ML vector is set equal to it. If the
% negative of the first PCA component is within +/- 0.5 for each component
% of the guess, then the ML vector is set equal to the negative of the
% first PCA component. If neither case is true, the ML vector is set to
% equal the first PCA component and a warning message is displayed. Sign
% errors with Euler angles (ie when the plot is upside down) later on can
% usually be traced back here. The code below generally works, but is not
% always the perfect fix for this particular problem. If you are running
% into issues, you can try changing the tolerance or the initial guess. You
% can also uncomment the lines of code concerned with switching the axis
% manually if you are confident you know which direction the vector should
% point.

tolerance = 0.5;
if abs(guess(1)-pcaCoef(1,1)') < tolerance && abs(guess(2)-pcaCoef(2,1)') < tolerance && abs(guess(3)-pcaCoef(3,1)') < tolerance
    ML_f = pcaCoef(:,1)';
elseif abs(-1*guess(1)-pcaCoef(1,1)') < tolerance && abs(-1*guess(2)-pcaCoef(2,1)') < tolerance && abs(-1*guess(3)-pcaCoef(3,1)') < tolerance
    ML_f = -1*pcaCoef(:,1)';
else
    ML_f = pcaCoef(:,1)';
    warning('PCA axis does not align well with guess, assuming correct direction!')
end
% % manualSwitch = input('Apply a negative sign to the first PCA component manually (1 for yes, otherwise just hit enter)');
% % if manualSwitch == 1
% %     ML_f = -1*ML_f;
% % end

% The AP vector is defined as the cross product between the ML and SI
% vectors. To ensure orthogonality of the coordinate system, the ML axis is
% redefined as the cross between the SI and AP vectors, then the AP vector
% is redefined as the cross between the new ML vector and the SI vector.

ML_f = ML_f./norm(ML_f);
AP = cross(SI,ML_f);
AP = AP./norm(AP);

ML = ML_f;


end

