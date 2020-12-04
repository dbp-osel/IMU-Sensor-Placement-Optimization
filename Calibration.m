%% Set up workspace
%This section clears every variable except for the structure 'data' and
%loads the structure 'data' from 'Trial Data Imported.mat' if it does not
%exist. This allows for this function to be run multiple times much faster
%since the large vatriable data does not need to be re-loaded every time.
%This can be helpful if generating multiple plots but can run into errors
%if the user happens to be using another variable named data.

% The cell array 'data' contains data from all of the files with the data
% from each file corresponding to a single row. The files were not re-named
% after the trial was completed and are named based on the default naming
% scheme of the Xsens Mtw Workstation which names each file using the base
% station serial number, followed by a time stamp, followed by the
% individual sensor's serial number. Matlab orders files in both numeric
% and alphabetical order prioritizing numbers over letters (ie 1 before 2,
% a before b, and any number before any letter). Because of this the data
% is organized in the following ways. Every 11 files correspond to the 11
% sensors in a given trial because they all share the same time stamp (ie
% the first 11 files are trial 1, the next 11 are trial 2). For each trial
% every sensor is listed in order based on sensor serial number.

clearvars -except data TaskList RVelColumns RowStart QuatColumnOne NumberSensors num AccColumns guess
dir = pwd;
if exist('data')
    loadNew = input('Do you want to use the same trial data? Yes or no (1/0)?');
else
    loadNew = 0;
end
if loadNew == 0
    disp('Load trial data ("P_xxx.mat")')
    cd(strcat(dir,'/Imported Subject Data'))
    [Datfile Datdir] = uigetfile;
    load(strcat(Datdir,Datfile))
    cd(dir)
    clc
end
%% Select Test
%Below is a list of every task performed in the order. The user is prompted
%to select a trial. The user's input is stored as the variable j which is
%fed into a linear equation which calculates the index (ind). The index
%corresponds to the first row corresponding to that given trial and is used
%later when referencing data.

if exist('TaskList')
    loadSub = input('Do you want to use the same subject information? Yes or no (1/0)?');
else
    loadSub = 0;
end
if loadSub == 0
    cd(strcat(dir,'/Subject Information'))
    disp('Load subject information ("P_xxx_Info.mat")')
    [SOfile SOdir] = uigetfile;
    load(strcat(SOdir,SOfile))
    cd(dir)
    clc
end
disp(TaskList)

j = input('Select Task');
       
ind = NumberSensors*(j-1)+1;
clc
%% Create acceleration and rotational velocity variables
%Acceleration and rotational velocity data are taken from the data matrices
%in row 2 of data. These data are indexed from these matices using the
%AccColumns and RVelColumns variables defined in the subject information
%files. Individual sensors are assigned to a segment based on ind, which
%defines where the data for this specific task starts in the file folder,
%and Num which is defined in the 'Subject Information' file and indicates
%the order the sensors appear in the folder based on MATLAB's alpha-numeric
%ordering scheme. Gravity free acceleration is found by subtracting gravity
%at a still point defined at frame 1. Net acceleration is found from the
%gravity free acceleration. Net rotational velocity is also found.
span = min([length(data{ind+num.LowThighAnt,2}(:,1)),length(data{ind+num.ShinBone,2}(:,1)), ...
    length(data{ind+num.LowShankLat,2}(:,1)),length(data{ind+num.Sacrum,2}(:,1)), ...
    length(data{ind+num.LowThighPos,2}(:,1)),length(data{ind+num.MidThighLat,2}(:,1)), ...
    length(data{ind+num.Heel,2}(:,1)), ...
    length(data{ind+num.LowThighLat,2}(:,1)),length(data{ind+num.L4L5,2}(:,1)), ...
    length(data{ind+num.MidShankLat,2}(:,1)),length(data{ind+num.DFoot,2}(:,1))]);
    %indexing to span eliminates the possibility of "index exceeds matrix
    %dimensions" type errors since sensors occasionally don't capture the
    %exact same number of data points.
% Acceleration
AccLowThighAnt = data{ind+num.LowThighAnt,2}(1:span,AccColumns); % Low Thigh Anterior
AccShinBone = data{ind+num.ShinBone,2}(1:span,AccColumns); % Shin Flat Bone
AccLowShankLat = data{ind+num.LowShankLat,2}(1:span,AccColumns); % Low Shank Lateral
AccSacrum = data{ind+num.Sacrum,2}(1:span,AccColumns); % Sacrum
AccLowThighPos = data{ind+num.LowThighPos,2}(1:span,AccColumns); % Low Thigh Posterior
AccMidThighLat = data{ind+num.MidThighLat,2}(1:span,AccColumns); % Middle Thigh Lateral
AccHeel = data{ind+num.Heel,2}(1:span,AccColumns); % Heel
AccLowThighLat = data{ind+num.LowThighLat,2}(1:span,AccColumns); % Low Thigh Lateral
AccL4L5 = data{ind+num.L4L5,2}(1:span,AccColumns); % L4-L5 Lumbar Spine
AccMidShankLat = data{ind+num.MidShankLat,2}(1:span,AccColumns); % Middle Shank Lateral
AccDFoot = data{ind+num.DFoot,2}(1:span,AccColumns); % Dorsal Foot

%Acceleration without gravity
stillpoint = 1;
Acc_noGLowThighAnt = [AccLowThighAnt(1:span,1)-AccLowThighAnt(stillpoint,1),AccLowThighAnt(1:span,2)-AccLowThighAnt(stillpoint,2),AccLowThighAnt(1:span,3)-AccLowThighAnt(stillpoint,3)]; % Low Shank, Lateral
Acc_noGShinBone = [AccShinBone(1:span,1)-AccShinBone(stillpoint,1),AccShinBone(1:span,2)-AccShinBone(stillpoint,2),AccShinBone(1:span,3)-AccShinBone(stillpoint,3)]; % Low Thigh, Lateral
Acc_noGLowShankLat = [AccLowShankLat(1:span,1)-AccLowShankLat(stillpoint,1),AccLowShankLat(1:span,2)-AccLowShankLat(stillpoint,2),AccLowShankLat(1:span,3)-AccLowShankLat(stillpoint,3)]; % 783 Mid Shank, Lateral
Acc_noGSacrum = [AccSacrum(1:span,1)-AccSacrum(stillpoint,1),AccSacrum(1:span,2)-AccSacrum(stillpoint,2),AccSacrum(1:span,3)-AccSacrum(stillpoint,3)]; % Low Thigh, Posterior
Acc_noGLowThighPos = [AccLowThighPos(1:span,1)-AccLowThighPos(stillpoint,1),AccLowThighPos(1:span,2)-AccLowThighPos(stillpoint,2),AccLowThighPos(1:span,3)-AccLowThighPos(stillpoint,3)]; % Shank, flat bone
Acc_noGMidThighLat = [AccMidThighLat(1:span,1)-AccMidThighLat(stillpoint,1),AccMidThighLat(1:span,2)-AccMidThighLat(stillpoint,2),AccMidThighLat(1:span,3)-AccMidThighLat(stillpoint,3)]; % Low Thigh, Anterior
Acc_noGHeel = [AccHeel(1:span,1)-AccHeel(stillpoint,1),AccHeel(1:span,2)-AccHeel(stillpoint,2),AccHeel(1:span,3)-AccHeel(stillpoint,3)]; % Dorsal Foot
Acc_noGLowThighLat = [AccLowThighLat(1:span,1)-AccLowThighLat(stillpoint,1),AccLowThighLat(1:span,2)-AccLowThighLat(stillpoint,2),AccLowThighLat(1:span,3)-AccLowThighLat(stillpoint,3)]; % Mid Thigh, Lateral
Acc_noGL4L5 = [AccL4L5(1:span,1)-AccL4L5(stillpoint,1),AccL4L5(1:span,2)-AccL4L5(stillpoint,2),AccL4L5(1:span,3)-AccL4L5(stillpoint,3)]; % L4-L5
Acc_noGMidShankLat = [AccMidShankLat(1:span,1)-AccMidShankLat(stillpoint,1),AccMidShankLat(1:span,2)-AccMidShankLat(stillpoint,2),AccMidShankLat(1:span,3)-AccMidShankLat(stillpoint,3)]; % Sacrum
Acc_noGDFoot = [AccDFoot(1:span,1)-AccDFoot(stillpoint,1),AccDFoot(1:span,2)-AccDFoot(stillpoint,2),AccDFoot(1:span,3)-AccDFoot(stillpoint,3)]; % Heel

AccNet = sqrt(sum(Acc_noGLowThighAnt.^2,2))+sqrt(sum(Acc_noGShinBone.^2,2))+sqrt(sum(Acc_noGLowShankLat.^2,2))+sqrt(sum(Acc_noGSacrum.^2,2))+sqrt(sum(Acc_noGLowThighPos.^2,2))+ ...
    sqrt(sum(Acc_noGMidThighLat.^2,2))+sqrt(sum(Acc_noGHeel.^2,2))+sqrt(sum(Acc_noGLowThighLat.^2,2))+sqrt(sum(Acc_noGL4L5.^2,2))+sqrt(sum(Acc_noGMidShankLat.^2,2))+sqrt(sum(Acc_noGDFoot.^2,2));

% Rotational Velocity
RotLowThighAnt = data{ind+num.LowThighAnt,2}(1:span,RVelColumns); % Low Thigh Anterior
RotShinBone = data{ind+num.ShinBone,2}(1:span,RVelColumns); % Shin Flat Bone
RotLowShankLat = data{ind+num.LowShankLat,2}(1:span,RVelColumns); % Low Shank Lateral
RotSacrum = data{ind+num.Sacrum,2}(1:span,RVelColumns); % Sacrum
RotLowThighPos = data{ind+num.LowThighPos,2}(1:span,RVelColumns); % Low Thigh Posterior
RotMidThighLat = data{ind+num.MidThighLat,2}(1:span,RVelColumns); % Middle Thigh Lateral
RotHeel = data{ind+num.Heel,2}(1:span,RVelColumns); % Heel
RotLowThighLat = data{ind+num.LowThighLat,2}(1:span,RVelColumns); % Low Thigh Lateral
RotL4L5 = data{ind+num.L4L5,2}(1:span,RVelColumns); % L4-L5 Lumbar Spine
RotMidShankLat = data{ind+num.MidShankLat,2}(1:span,RVelColumns); % Middle Shank Lateral
RotDFoot = data{ind+num.DFoot,2}(1:span,RVelColumns); % Dorsal Foot

RotNet = sqrt(sum(RotLowThighAnt.^2,2))+sqrt(sum(RotShinBone.^2,2))+sqrt(sum(RotLowShankLat.^2,2))+sqrt(sum(RotSacrum.^2,2))+sqrt(sum(RotLowThighPos.^2,2))+ ...
    sqrt(sum(RotMidThighLat.^2,2))+sqrt(sum(RotHeel.^2,2))+sqrt(sum(RotLowThighLat.^2,2))+sqrt(sum(RotL4L5.^2,2))+sqrt(sum(RotMidShankLat.^2,2))+sqrt(sum(RotDFoot.^2,2));
%% Selecting Calibration approaches
%Two calibration approaches can be implemented with this script. The
%functional calibration uses an outwalk protocol while the static
%calibration uses static poses.
clc
disp('Select which calibration(s) you would like to perform')
disp('The calibration task accomodates both a static approach and a functional')
disp('approach. Both approaches use gravitational acceleration when the subject')
disp('stands still to determine the SI axes of each sensor. The static approach')
disp('then uses a seated pose to define the ML axes. During this pose the')
disp('subject leans back and extends their legs. Care is taken so that the long')
disp('axes of each body segment lie parallel to the sagittal plane. Cross products')
disp('between the SI axes and the gravitational acceleration vectors during this')
disp('pose yeilds the ML axes since they both lie parallel to the sagittal plane.')
disp('The functional approach uses principal component analysis (PCA) as the')
disp('subject walks in order to find the 1st component of rotational velocity.')
disp('The main component is assumed to be the component about ML. Cross products')
disp('between SI and ML give AP.')

functionalTru = input('Perform functional Calibration? (1 for yes otherwise just hit enter)');
staticTru = input('Perform Static Calibration? (1 for yes otherwise just hit enter)');
%% Defining ranges of still standing, still sitting, and steady state gait
%The following interactive functions use input from the user to determine
%regions of still standing and steady state gait. These functions look
%exclusively at data from the heel sensor which provides (in my opinion)
%more readily interpretable waveforms. This procedure could be done with
%other sensors. These functions are easier to run through section by
%section in case the user makes a mistake and must re-run one of these
%functions.

% Defining regions of still standing
%The user will be prompted to click two points on the heel acceleration
%plot corresponding to the range of time when the subject is standing
%still and when the subject is seated.
[RangeStand,RangeSit] = StationaryCalibrationRanges(AccHeel);

if functionalTru == 1
    % Zero angle instances initial estimate
    %The user will be prompted to scope out minima points on acceleration and
    %rotational velocity plots. The user will then provide a threshold value by
    %which to treat these points as coincidental. These points are used to
    %estimate instances of flat foot.
    [instances_initial,threshold] = ZeroVelocityFinder(AccNet,RotNet,RangeStand,RangeSit);
    
    % Zero angle instances refinement and heel waveform
    %The user will be asked to identify the primary component of angle. This
    %component will be used to further refine the flat foot instances and
    %calculate an estimate of drift corrected angle.
    [Heel,instances] = Angle_Corrected(RotHeel,instances_initial,threshold);
    
    % Defining regions of steady state gait
    % The estimate of drift corrected angle as well as the corrected flat foot
    % instances are used to identify regions of steady state gait. The user
    % will be prompted to identify the number of times steady state gait occurs
    % as well as the span of time which each occurs for.
    [RangeSteady] = DefiningSteadyState(Heel,instances);
end
%% Perform estimations of anatomical axes: Method 1, functional calibration

% This section uses the raw acceleration and raw rotational velocity data
% from each sensor to estimate the superior-inferior, medial-lateral, and
% anterior-posterior anatomical directions of each body segment in terms of
% the sensor coordinate systems. 

% The superior-inferior axes are assumed to be parallel with gravity when
% the subject is standing still. As such these axes are defined in terms of
% unit vectors parallel to the average acceleration vector (acceleration
% should be due to gravity alone) when the subject is standing still. The
% medial-lateral axis is assumed to be the primary axis of rotation during
% steady state gait for all sensors except the torso sensors. The
% medial-lateral axis is approximated as the first component found through
% Principle Component Analysis (PCA) using Matlab's 'pca' function for the
% data points within the steady state ranges. The medial-lateral axes of
% the torso sensors cannot be readily found using this method because the
% torso does not rotate about the medial-lateral axis significantly during
% steady state gait. Instead the medial lateral axes of each torso sensor
% is assumed to be parallel to the sensor's [0 1 0] axis. The user should
% take care to align this axis of these sensors when attaching these
% sensors to each subject.

clc
if functionalTru == 1

    % Low Thigh Anterior
    [SI_LowThighAnt,ML_LowThighAnt,AP_LowThighAnt] = Axis_Estimation(AccLowThighAnt,RotLowThighAnt,RangeStand,RangeSteady,guess.LowThighAnt);
    Zero.LowThighAnt = [SI_LowThighAnt;ML_LowThighAnt;AP_LowThighAnt;];
    % Shin Flat Bone
    [SI_ShinBone,ML_ShinBone,AP_ShinBone] = Axis_Estimation(AccShinBone,RotShinBone,RangeStand,RangeSteady,guess.ShinBone);
    Zero.ShinBone = [SI_ShinBone;ML_ShinBone;AP_ShinBone;];
    % Low Shank Lateral
    [SI_LowShankLat,ML_LowShankLat,AP_LowShankLat] = Axis_Estimation(AccLowShankLat,RotLowShankLat,RangeStand,RangeSteady,guess.LowShankLat);
    Zero.LowShankLat = [SI_LowShankLat;ML_LowShankLat;AP_LowShankLat;];
    % Sacrum
    [SI_Sacrum] = Axis_Estimation_SI_Only(AccSacrum,RangeStand);
    AP_Sacrum = cross(SI_Sacrum,guess.Sacrum);
    ML_Sacrum = guess.Sacrum;
    
%     ML_Sacrum = cross(AP_Sacrum,SI_Sacrum);
%     AP_Sacrum = cross(SI_Sacrum,ML_Sacrum);
%     AP_Sacrum = AP_Sacrum./norm(AP_Sacrum);
%     ML_Sacrum = ML_Sacrum./norm(ML_Sacrum);
%     SI_Sacrum = SI_Sacrum./norm(SI_Sacrum);
    Zero.Sacrum = [SI_Sacrum;ML_Sacrum;AP_Sacrum;];
    % Low Thigh Posterior
    [SI_LowThighPos,ML_LowThighPos,AP_LowThighPos] = Axis_Estimation(AccLowThighPos,RotLowThighPos,RangeStand,RangeSteady,guess.LowThighPos);
    Zero.LowThighPos = [SI_LowThighPos;ML_LowThighPos;AP_LowThighPos;];
    % Middle Thigh Lateral
    [SI_MidThighLat,ML_MidThighLat,AP_MidThighLat] = Axis_Estimation(AccMidThighLat,RotMidThighLat,RangeStand,RangeSteady,guess.MidThighLat);
    Zero.MidThighLat = [SI_MidThighLat;ML_MidThighLat;AP_MidThighLat;];
    % Heel
    [SI_Heel,ML_Heel,AP_Heel] = Axis_Estimation(AccHeel,RotHeel,RangeStand,RangeSteady,guess.Heel);
    Zero.Heel = [SI_Heel;ML_Heel;AP_Heel;];
    % Low Thigh Lateral
    [SI_LowThighLat,ML_LowThighLat,AP_LowThighLat] = Axis_Estimation(AccLowThighLat,RotLowThighLat,RangeStand,RangeSteady,guess.LowThighLat);
    Zero.LowThighLat = [SI_LowThighLat;ML_LowThighLat;AP_LowThighLat;];
    % L4-L5 Lumbar Spine
    [SI_L4L5] = Axis_Estimation_SI_Only(AccL4L5,RangeStand);
    AP_L4L5 = cross(SI_L4L5,guess.L4L5);
    ML_L4L5 = guess.L4L5;
    
%     ML_L4L5 = cross(AP_L4L5,SI_L4L5);
%     AP_L4L5 = cross(SI_L4L5,ML_L4L5);
%     AP_L4L5 = AP_L4L5./norm(AP_L4L5);
%     ML_L4L5 = ML_L4L5./norm(ML_L4L5);
%     SI_L4L5 = SI_L4L5./norm(SI_L4L5);
    Zero.L4L5 = [SI_L4L5;ML_L4L5;AP_L4L5;];
    % Middle Shank Lateral
    [SI_MidShankLat,ML_MidShankLat,AP_MidShankLat] = Axis_Estimation(AccMidShankLat,RotMidShankLat,RangeStand,RangeSteady,guess.MidShankLat);
    Zero.MidShankLat = [SI_MidShankLat;ML_MidShankLat;AP_MidShankLat;];
    % Dorsal Foot
    [SI_DFoot,ML_DFoot,AP_DFoot] = Axis_Estimation(AccDFoot,RotDFoot,RangeStand,RangeSteady,guess.DFoot);
    Zero.DFoot = [SI_DFoot;ML_DFoot;AP_DFoot;];
    %% Save functional calibration file
    disp('Save the functional calibration file')  
    cd(strcat(dir,'\Calibrations'))
    uisave('Zero')
    cd(dir)
    clc
end
%% Perform estimations of anatomical axes: Method 2, stationary calibration
%This calibration method uses two vectors obtained from gravitational
%acceleration while the subject is stationary. The first vector is obtained
%when the subject is standing still and, like with the previous method,
%this vector is assumed to lie on the superior-inferior axis. The second
%vector is obtained when the subject is sitting whith their legs stretched
%out and their torso leaned back. The vectors obtained in this manner are
%co-planer to the vectors obtained when the subject is standing. A cross
%product is taken between the seated and standing vectors to define the
%medial-lateral axis. Another cross product is taken between the
%medial-lateral and superior-inferior vectors to get the anterior-posterior
%vector. Further cross products between these three vectors are performed
%to make sure the coordinate system being defined is orthogonal.


clc
if staticTru == 1
    clear Zero
    % Low Thigh Anterior
    [SI_LowThighAnt,ML_LowThighAnt,AP_LowThighAnt] = Axis_Estimation_stationary(AccLowThighAnt,RangeStand,RangeSit);
    Zero.LowThighAnt = [SI_LowThighAnt;ML_LowThighAnt;AP_LowThighAnt;];
    % Shin Flat Bone
    [SI_ShinBone,ML_ShinBone,AP_ShinBone] = Axis_Estimation_stationary(AccShinBone,RangeStand,RangeSit);
    Zero.ShinBone = [SI_ShinBone;ML_ShinBone;AP_ShinBone;];
    % Low Shank Lateral
    [SI_LowShankLat,ML_LowShankLat,AP_LowShankLat] = Axis_Estimation_stationary(AccLowShankLat,RangeStand,RangeSit);
    Zero.LowShankLat = [SI_LowShankLat;ML_LowShankLat;AP_LowShankLat;];
    % Sacrum
    [SI_Sacrum,ML_Sacrum,AP_Sacrum] = Axis_Estimation_stationary(AccSacrum,RangeStand,RangeSit);
    Zero.Sacrum = [SI_Sacrum;ML_Sacrum;AP_Sacrum;];
    % Low Thigh Posterior
    [SI_8C6,ML_8C6,AP_8C6] = Axis_Estimation_stationary(AccLowThighPos,RangeStand,RangeSit);
    Zero.LowThighPos = [SI_8C6;ML_8C6;AP_8C6;];
    % Middle Thigh Lateral
    [SI_MidThighLat,ML_MidThighLat,AP_MidThighLat] = Axis_Estimation_stationary(AccMidThighLat,RangeStand,RangeSit);
    Zero.MidThighLat = [SI_MidThighLat;ML_MidThighLat;AP_MidThighLat;];
    % Heel
    [SI_Heel,ML_Heel,AP_Heel] = Axis_Estimation_stationary(AccHeel,RangeStand,RangeSit);
    Zero.Heel = [SI_Heel;ML_Heel;AP_Heel;];
    % Low Thigh Lateral
    [SI_LowThighLat,ML_LowThighLat,AP_LowThighLat] = Axis_Estimation_stationary(AccLowThighLat,RangeStand,RangeSit);
    Zero.LowThighLat = [SI_LowThighLat;ML_LowThighLat;AP_LowThighLat;];
    % L4-L5 Lumbar Spine
    [SI_L4L5,ML_L4L5,AP_L4L5] = Axis_Estimation_stationary(AccL4L5,RangeStand,RangeSit);
    Zero.L4L5 = [SI_L4L5;ML_L4L5;AP_L4L5;];
    % Middle Shank Lateral
    [SI_MidShankLat,ML_MidShankLat,AP_MidShankLat] = Axis_Estimation_stationary(AccMidShankLat,RangeStand,RangeSit);
    Zero.MidShankLat = [SI_MidShankLat;ML_MidShankLat;AP_MidShankLat;];
    % Dorsal Foot
    [SI_DFoot,ML_DFoot,AP_DFoot] = Axis_Estimation_stationary(AccDFoot,RangeStand,RangeSit);
    Zero.DFoot = [SI_DFoot;ML_DFoot;AP_DFoot;];
    %% Save static calibration file
    disp('Save the static calibration file')
    cd(strcat(dir,'\Calibrations'))
    uisave('Zero')
    cd(dir)
    clc
end


