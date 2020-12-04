%% Load Optical Motion Capture Data
% This checks the existence of the Optical_Angles variable. If it exists,
% the user is asked whether or not they wish to re-use the optical motion
% capture data already loaded. If the answer is no or if the variable
% 'Optical_Angles' is not present, then the user will be prompted to load
% the appropriate file using the uigetfile function. These .mat files can
% be found in the 'Imported Optical Motion Capture Data' folder.
clearvars -except Optical_Angles Optical_tasks flexangHip flexangKnee flexangAnkle abdangHip ...
    abdangKnee abdangAnkle rotangHip rotangKnee rotangAnkle Hip_Names Knee_Names Ankle_Names ...
    Optical_Events Optical_inds
dir = pwd;
try
    amIreal = Optical_tasks{1,1};
    loadOpt = input('Use the same optical motion capture data? (1 for yes/0 for no)');
    if loadOpt == 0
        cd(strcat(dir,'/Optical Motion Capture Angles'))
        disp('Load optical motion capture data')
        [Datfile Datdir] = uigetfile;
        cd(dir)
        load(strcat(Datdir,Datfile))
    end
catch
    disp('Load optical motion capture data')
    cd(strcat(dir,'/Optical Motion Capture Angles'))
    [Datfile Datdir] = uigetfile;
    cd(dir)
    load(strcat(Datdir,Datfile))     
end
clear amIreal
clc
%% Select the appropriate task and extract angles
for pp = 1:length(Optical_tasks)
    disp(strcat(num2str(pp),':',32,Optical_tasks{pp},10))
end
OptTaskNum = input('Select task');
%% Storing Vicon angles into variables and inserting NaNs into empty data regions
% This section uses the number of the task selected by the user in the
% previous section to identify and store the appropraite joint angle into
% variables. Vicon inserts zero angle values during the times when it
% cannot calculate angles due to a loss of marker visibility. These values
% are replaced with NaNs to make it easier to distinguish between real and
% non-real data.

OGRightHip = Optical_Angles{OptTaskNum,1}.RHipAngles;
OGRightKnee = Optical_Angles{OptTaskNum,1}.RKneeAngles;
OGRightAnkle = Optical_Angles{OptTaskNum,1}.RAnkleAngles;

%These next few lines test for two conditions, 1) if the value of angle is
%zero and 2) if the slope of the graph is zero. The slope is checked simply
%by subtracting adjacent points (if they are the same the slope is zero).
%This is done in both directions (so that the points n-1 and n+1 are both
%subjtracted from n) and the results of these are given as two seperate
%logical vectors. These vectors are added ("OR" boolean logic) producing
%another logical vector that is true when the slope is zero going away from
%the point in question in either direction. The points where the angle is
%zero are found as a logical vector. The angle logical vector and slope
%logical vector are multiplied ("AND" boolean logic). The find function is
%used to convert the resulting logicals to indices and the angles
%corresponding to these indeces are set equal to NaN. This should fix the
%issue where Vicon sets empty gaps in the data to equal zero. It is not a
%perfect fix as it is possible to have zero angle and zero slope
%simultaneously, however this scenario is practically impossible given that
%there are 13 decimal places and random noise. Furthermore even if a frame
%or two was misidentified it would have a negligible effect on the
%subsequent calculations.

OGRightHip(find((OGRightHip(:,1) == 0).*[0;((OGRightHip(2:end,1)-(OGRightHip(1:end-1,1))) == 0)] ...
    +[((OGRightHip(1:end-1,1)-(OGRightHip(2:end,1))) == 0);0]),:) = NaN;
OGRightKnee(find((OGRightKnee(:,1) == 0).*[0;((OGRightKnee(2:end,1)-(OGRightKnee(1:end-1,1))) == 0)] ...
    +[((OGRightKnee(1:end-1,1)-(OGRightKnee(2:end,1))) == 0);0]),:) = NaN;
OGRightAnkle(find((OGRightAnkle(:,1) == 0).*[0;((OGRightAnkle(2:end,1)-(OGRightAnkle(1:end-1,1))) == 0)] ...
    +[((OGRightAnkle(1:end-1,1)-(OGRightAnkle(2:end,1))) == 0);0]),:) = NaN;
% These angles are zeroed in the same manner as the IMU data

% OGRightHip = [OGRightHip(:,1)-OGRightHip(startHip,1),OGRightHip(:,2)-OGRightHip(startHip,2),OGRightHip(:,3)-OGRightHip(startHip,3)];
% OGRightKnee = [OGRightKnee(:,1)-OGRightKnee(startKnee,1),OGRightKnee(:,2)-OGRightKnee(startKnee,2),OGRightKnee(:,3)-OGRightKnee(startKnee,3)];
% OGRightAnkle = [OGRightAnkle(:,1)-OGRightAnkle(startAnkle,1),OGRightAnkle(:,2)-OGRightAnkle(startAnkle,2),OGRightAnkle(:,3)-OGRightAnkle(startAnkle,3)];

%% Load IMU Motion Capture Data
% This checks the existence of the Hip_Names, Knee_Names, and Ankle_Names
% variables. If any exist, the user is asked whether or not they wish to
% re-use the IMU motion capture data already loaded. If the answer is no or
% if these variable are not present, then the user will be prompted to load
% the appropriate file using the uigetfile function. These .mat files can
% be found in the 'IMU Comparison Data' folder.

if exist('Hip_Names','var')||exist('Knee_Names','var')||exist('Ankle_Names','var')
    loadOpt = input('Use the same IMU motion capture data? (1 for yes/0 for no)');
    if loadOpt == 0
        cd(strcat(dir,'/IMU Angles'))
        disp('Load IMU motion capture data')
        [Datfile Datdir] = uigetfile;
        cd(dir)
        load(strcat(Datdir,Datfile))
    end
else
    cd(strcat(dir,'/IMU Angles'))
    disp('Load IMU motion capture data')
    [Datfile Datdir] = uigetfile;
    cd(dir)
    load(strcat(Datdir,Datfile))
end
clc

%% Input Frame Rates of IMUs and Optical Motion Capture Systems
FR1 = input('IMU motion capture frame rate');
FR2 = input('Optical motion capture frame rate');
%% Resample and Crop data
[flexangHip_resamp,abdangHip_resamp,rotangHip_resamp,RightHip] = resampling(flexangHip,abdangHip,rotangHip,OGRightHip,FR1,FR2);
[flexangKnee_resamp,abdangKnee_resamp,rotangKnee_resamp,RightKnee] = resampling(flexangKnee,abdangKnee,rotangKnee,OGRightKnee,FR1,FR2);
[flexangAnkle_resamp,abdangAnkle_resamp,rotangAnkle_resamp,RightAnkle] = resampling(flexangAnkle,abdangAnkle,rotangAnkle,OGRightAnkle,FR1,FR2);
%% Offset IMU data to match Optical motion capture data at the start of the trial
disp('Offset IMU motion capture data so that it equals')
disp('the optical motion capture data on the first frame?')
zerod = input('(1 for yes otherwise just hit enter)');
if zerod == 1
    startHip = 1;
    if isnan(RightHip(startHip,1))
        for jj = 1:length(RightHip(:,1))
            if ~isnan(RightHip(jj,1))
                startHip = jj;
                warning(strcat('Offsetting hip angles based on the',32,num2str(startHip),32,'frame.',10, ...
                    'Check plots to verify correct offset has occured.'))
                break
            end
        end
    end
    startKnee = 1;
    if isnan(RightKnee(startKnee,1))
        for jj = 1:length(RightKnee(:,1))
            if ~isnan(RightKnee(jj,1))
                startKnee = jj;
                warning(strcat('Offsetting knee angles based on the',32,num2str(startKnee),32,'frame.',10, ...
                    'Check plots to verify correct offset has occured.'))
                break
            end
        end
    end
    startAnkle = 1;
    if isnan(RightAnkle(startAnkle,1))
        for jj = 1:length(RightAnkle(:,1))
            if ~isnan(RightAnkle(jj,1))
                startAnkle = jj;
                warning(strcat('Offsetting ankle angles based on the',32,num2str(startAnkle),32,'frame.',10, ...
                    'Check plots to verify correct offset has occured.'))
                break
            end
        end
    end
    [ispan jspan] = size(rotangHip_resamp);
    for ii = 1:ispan
        for jj = 1:jspan
            flexangHip_resamp{ii,jj} = flexangHip_resamp{ii,jj} - RightHip(startHip,1);  
            abdangHip_resamp{ii,jj} = abdangHip_resamp{ii,jj} + RightHip(startHip,2); 
            rotangHip_resamp{ii,jj} = rotangHip_resamp{ii,jj} + RightHip(startHip,3); 
        end
    end
    [ispan jspan] = size(rotangKnee_resamp);
    for ii = 1:ispan
        for jj = 1:jspan
            flexangKnee_resamp{ii,jj} = flexangKnee_resamp{ii,jj} + RightKnee(startKnee,1); 
            abdangKnee_resamp{ii,jj} = abdangKnee_resamp{ii,jj} + RightKnee(startKnee,2); 
            rotangKnee_resamp{ii,jj} = rotangKnee_resamp{ii,jj} + RightKnee(startKnee,3); 
        end
    end
    [ispan jspan] = size(rotangAnkle_resamp);
    for ii = 1:ispan
        for jj = 1:jspan
            flexangAnkle_resamp{ii,jj} = flexangAnkle_resamp{ii,jj} - RightAnkle(startAnkle,1); 
            abdangAnkle_resamp{ii,jj} = abdangAnkle_resamp{ii,jj} + RightAnkle(startAnkle,2); 
            rotangAnkle_resamp{ii,jj} = rotangAnkle_resamp{ii,jj} + RightAnkle(startAnkle,3); 
        end
    end
end
%% Generate plots
% These variables are used later in plotting to generate unique plots for
% every sensor combination. Each sensor above the joint is given a unique
% color and every sensor below the joint is given a unique line style.
colors = {'r','b','g','k'};
lines = {'-','--','-.',':'};
%% Hip Angles
figure
[espan fspan] = size(flexangHip_resamp);
legInd = 1;
leg = {};
for ee = 1:espan
    for ff = 1:fspan
        %   Flexion
        subplot(3,1,1)
        hold on
        plot(1:length(flexangHip_resamp{ee,ff}),-1*flexangHip_resamp{ee,ff},strcat(colors{ee},lines{ff}))
        hold off
        %   Abduction
        subplot(3,1,2)
        hold on
        plot(1:length(abdangHip_resamp{ee,ff}),abdangHip_resamp{ee,ff},strcat(colors{ee},lines{ff}))
        hold off
        %   Rotation
        subplot(3,1,3)
        hold on
        plot(1:length(rotangHip_resamp{ee,ff}),rotangHip_resamp{ee,ff},strcat(colors{ee},lines{ff}))
        hold off
        
        leg{legInd} = Hip_Names{ee,ff};
        legInd = legInd + 1;
    end
end

%   Flexion
subplot(3,1,1)
hold on
OptflexplotHip = plot(1:length(RightHip(:,1)),RightHip(:,1),'k');
title(strcat('Hip Angles',10,'Flexion'))
ylabel('Degrees')
hold off
%   Abduction
subplot(3,1,2)
hold on
OptabdplotHip = plot(1:length(RightHip(:,2)),RightHip(:,2),'k');
title('Abduction')
ylabel('Degrees')
hold off
%   Rotation
subplot(3,1,3)
hold on
OptrotplotHip = plot(1:length(RightHip(:,3)),RightHip(:,3),'k');
title('Rotation')
ylabel('Degrees')
hold off

leg{end+1} = 'Optical';
hold on
subplot(3,1,1)
legend(leg)
hold off

%% Knee Angles
figure
[espan fspan] = size(flexangKnee_resamp);
legInd = 1;
leg = {};
for ee = 1:espan
    for ff = 1:fspan
        %   Flexion
        subplot(3,1,1)
        hold on
        plot(1:length(flexangKnee_resamp{ee,ff}),flexangKnee_resamp{ee,ff},strcat(colors{ee},lines{ff}))
        hold off
        %   Abduction
        subplot(3,1,2)
        hold on
        plot(1:length(abdangKnee_resamp{ee,ff}),abdangKnee_resamp{ee,ff},strcat(colors{ee},lines{ff}))
        hold off
        %   Rotation
        subplot(3,1,3)
        hold on
        plot(1:length(rotangKnee_resamp{ee,ff}),rotangKnee_resamp{ee,ff},strcat(colors{ee},lines{ff}))
        hold off
        
        leg{legInd} = Knee_Names{ee,ff};
        legInd = legInd + 1;
    end
end
%   Flexion
subplot(3,1,1)
hold on
plot(1:length(RightKnee(:,1)),RightKnee(:,1),'k')
title(strcat('Knee Angles',10,'Flexion'))
ylabel('Degrees')
hold off
%   Abduction
subplot(3,1,2)
hold on
plot(1:length(RightKnee(:,2)),RightKnee(:,2),'k')
title('Abduction')
ylabel('Degrees')
hold off
%   Rotation
subplot(3,1,3)
hold on
plot(1:length(RightKnee(:,3)),RightKnee(:,3),'k')
title('Rotation')
ylabel('Degrees')
hold off

leg{end+1} = 'Optical';
hold on
subplot(3,1,1)
legend(leg)
hold off
%% Ankle Angles
figure
[espan fspan] = size(flexangAnkle_resamp);
legInd = 1;
leg = {};
for ee = 1:espan
    for ff = 1:fspan
        %   Flexion
        subplot(3,1,1)
        hold on
        plot(1:length(flexangAnkle_resamp{ee,ff}),-1*flexangAnkle_resamp{ee,ff},strcat(colors{ee},lines{ff}))
        hold off
        %   Abduction
        subplot(3,1,2)
        hold on
        plot(1:length(abdangAnkle_resamp{ee,ff}),abdangAnkle_resamp{ee,ff},strcat(colors{ee},lines{ff}))
        hold off
        %   Rotation
        subplot(3,1,3)
        hold on
        plot(1:length(rotangAnkle_resamp{ee,ff}),rotangAnkle_resamp{ee,ff},strcat(colors{ee},lines{ff}))
        hold off
        
        leg{legInd} = Ankle_Names{ee,ff};
        legInd = legInd + 1;
    end
end
%   Flexion
subplot(3,1,1)
hold on
plot(1:length(RightAnkle(:,1)),RightAnkle(:,1),'k')
title(strcat('Ankle Angles',10,'Flexion'))
ylabel('Degrees')
hold off
%   Abduction
subplot(3,1,2)
hold on
plot(1:length(RightAnkle(:,2)),RightAnkle(:,2),'k')
title('Abduction')
ylabel('Degrees')
hold off
%   Rotation
subplot(3,1,3)
hold on
plot(1:length(RightAnkle(:,3)),RightAnkle(:,3),'k')
title('Rotation')
ylabel('Degrees')
hold off

leg{end+1} = 'Optical';
hold on
subplot(3,1,1)
legend(leg)
hold off
