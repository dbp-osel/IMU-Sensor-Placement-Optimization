function [reflex reabd rerot reOpt] = resampling(flexang,abdang,rotang,comp,FR1,FR2)
%% Re-sampling data based on system frame rates
% This section calculates the times in seconds when the task ends according
% to both systems based on the frame rates provided. This is used to find
% which system runs for a shorter period of time (if there is any
% discrepency). This shorter time is multiplied by the frame rates of each
% system to find the the final time point when both systems are running in
% terms of the frames of each system. The data from each system is then
% re-sampled. The sample times are rounded to the nearest frame introducing
% some error due to rounding. This error can be minimized by making sure
% the larger frame rate is a multiple of the smaller frame rate (for
% instance, if Vicon is running at 120Hz and the IMUs at 60Hz then the IMUs
% will be sampled at every frame and Vicon at every other frame).

% Inputs
% flexang: The flexion/extension angle vector calculated from the IMU data
% abdang: The abduction/adduction angle vector calculated from the IMU data
% rotang: The internal/external rotation angle vector calculated from the IMU data
% comp: Comparison data from the optical system as an nx3 matrix. Column 1
% 	is flexion/extension, column 2 is abduction/adduction, and column 3 is
% 	internal/external rotation.
% FR1: Frame rate of the IMU system
% FR2: Frame rate of the optical system

% Outputs
% reflex: The resampled flexion/extension angle from the IMU data
% reabd: The resampled abduction/adduction angle from the IMU data
% rerot: The resampled internal/external rotation from the IMU data
% reOpt: The resampled optical motion capture angles
%   sampleInds: 

spanIMU = length(rotang{1,1}(1,:));
spanOpt = length(comp(:,1));

minTime = min([FR1^-1*spanIMU,FR2^-1*spanOpt]);

endIMU = FR1*minTime;
endOpt = FR2*minTime;

minEnd = min([endIMU,endOpt]);
sampleIMU = round(1:endIMU/minEnd:endIMU);
sampleOpt = round(1:endOpt/minEnd:endOpt);

logicInds = false(spanOpt,1);

if length(sampleIMU) ~= length(sampleOpt)
    % The sampling vectors sampleIMU and sampleOpt may be different lengths
    % due to rounding error. This if statement handels this case by
    % truncating the longer of the two so that it is the same length as the
    % shorter. A warning message is also displayed conveying the lengths of
    % each vector before truncation. A difference of a frame or two here
    % and there is to be expected, but large differences may indicate
    % error.
    lengths = [length(sampleIMU) length(sampleOpt)];
    [minSpan,whichMin] = min(lengths);
    warning(strcat('Length of IMU sample vector:',32,num2str(lengths(1)),10, ...
        'End of Vicon sample vector:',32,num2str(lengths(2)),10, ...
        'The larger vector was truncated to the length of the shorter vector.'))      
    if whichMin == 1
        sampleOpt = sampleOpt(1:minSpan);
    elseif whichMin == 2
        sampleIMU = sampleIMU(1:minSpan);
    end    
end
%% Handeling cases of NaNs
% This section identifies cases of nonexistant numbers (NaNs) introduced in
% the Comparison_Tool script. These NaNs are applied to regions of the
% Vicon data where the angle values are set to zero due to a loss of
% visable markers. 
notNaNs = ~isnan(comp(sampleOpt,1));
%%
realSampleOpt = sampleOpt(notNaNs);
reOpt = comp(realSampleOpt,:);

[ispan jspan] = size(rotang);
for ii = 1:ispan
    for jj = 1:jspan
        reflex{ii,jj} = flexang{ii,jj}(sampleIMU(notNaNs));
        reabd{ii,jj} = abdang{ii,jj}(sampleIMU(notNaNs));
        rerot{ii,jj} = rotang{ii,jj}(sampleIMU(notNaNs));
    end
end