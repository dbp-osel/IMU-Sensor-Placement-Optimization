%% Set up workspace
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

clearvars -except data TaskList RVelColumns RowStart QuatColumnOne NumberSensors num AccColumns guess FR
dir = pwd;
if exist('data','var')
    loadNew = input('Do you want to use the same trial data? Yes or no (1/0)?');
else
    loadNew = 0;
end
if loadNew == 0
    cd(strcat(dir,'/Imported Subject Data'))
    disp('Load trial data ("P_xxx.mat")')
    [Datfile Datdir] = uigetfile;
    cd(dir)
    load(strcat(Datdir,Datfile))
    clc
end
%% Load Calibration File
% In this section the user is prompted to load a calibration file. This
% file will contain a structure called 'Zero' with unit vectors aligned
% along all of the corrected anatomical axes for each sensor.
if exist('TaskList','var')
    loadSub = input('Do you want to use the same subject information? Yes or no (1/0)?');
else
    loadSub = 0;
end
if loadSub == 0
    cd(strcat(dir,'\Subject Information'))
    disp('Load subject information ("P_xxx_Info.mat")')
    [SOfile SOdir] = uigetfile;
    load(strcat(SOdir,SOfile))
    clc
    cd(dir)
end

disp(TaskList)
cd(strcat(dir,'\Calibrations'))
disp('Select calibration file')
[Cfile Cdir] = uigetfile;
load(strcat(Cdir,Cfile))
cd(dir)
clc

%% Select Test
%Below is a list of every task performed in the order. The user is prompted
%to select a trial. The user's input is stored as the variable j which is
%fed into a linear equation which calculates the index (ind). The index
%corresponds to the first row corresponding to that given trial and is used
%later when referencing data.
disp(TaskList)

j = input('Select Task');
       
ind = NumberSensors*(j-1)+1;
clc
%% Select Joints
% The user is prompted to select which joints to calculate Euler angles
% about.
disp(strcat('1: Hip Angles',10,'2: Knee Angles',10,'3: Ankle Angles',10,'4: All Joints, All Angles'))
jointOpt = input('Select Joint (enter all that apply as a vector)');
HipOpt = 0;
KneeOpt = 0;
AnkleOpt = 0;
skip = 0;
for qq = 1:length(jointOpt)
    if jointOpt(qq) == 1
        HipOpt = 1;
    elseif jointOpt(qq) == 2
        KneeOpt = 1;
    elseif jointOpt(qq) == 3
        AnkleOpt = 1;
    elseif jointOpt(qq) == 4
        HipOpt = 1;
        KneeOpt = 1;
        AnkleOpt = 1;
        FootOpt = [1 2];
        ShankOpt = [1 2 3];
        ThighOpt = [1 2 3 4];
        TorsoOpt = [1 2];
        skip = 1;  
    end
end
clc
%% Select Sensors
% The user is prompted to select the applicable sensors to use given the
% joints they selected. Any pair of sensors from the segments above and
% below the joint in question can be used to calculate the Euler angles
% about that joint. After the user selects the sensor they want to look at,
% Euler angles will be calculated for every possible combination of the
% sensors selected.
if skip == 0
    
    FootOpt = [];
    ShankOpt = [];
    ThighOpt = [];
    TorsoOpt = [];
    
    if AnkleOpt == 1
        disp(strcat('1: Heel',10,'2: Dorsal Foot'))
        FootOpt = input('Select Foot Sensors (enter all that apply as a vector)');
        %1 => Heel
        %2 => Dorsal Foot
        clc
    end
    
    if AnkleOpt == 1 || KneeOpt == 1
        disp(strcat('1: Shin Flat Bone',10,'2: MidShank Lateral',10,'3: LowShank Lateral'))
        ShankOpt = input('Select Shank Sensors (enter all that apply as a vector)');
        %1 => Shin Flat Bone
        %2 => MidShank Lateral
        %3 => LowShank Lateral
        clc
    end
    
    if KneeOpt == 1 || HipOpt == 1
        disp(strcat('1: LowThigh Anterior',10,'2: MidThigh Lateral',10, ...
            '3: LowThigh Lateral',10,'4: LowThigh Posterior'))
        ThighOpt = input('Select Thigh Sensors (enter all that apply as a vector)');
        %1 => LowThigh Anterior
        %2 => MidThigh Lateral
        %3 => LowThigh Lateral
        %4 => LowThigh Posterior
        clc
    end
    
    if HipOpt == 1
        disp(strcat('1: Sacrum',10,'2: L4-L5'))
        TorsoOpt = input('Select Torso Sensors (enter all that apply as a vector)');
        %1 => Sacrum
        %2 => L4-L5
        clc
    end
end
%% Perform Rotations on Body Segment Coordinate Systems
% This section performs calculations to rotate initial unit vectors
% provided by the calibration file using the quaternion outputs of the
% sensors. The result is a series of coordinate systems defined by
% orthogonal unit vectors that are rotated so that at each time point they
% describe the orientation of the given sensor in space. This is performed
% by the 'CoordinateTransformQuatPreSpecifiedAxes' function which takes the
% quaternions and the three orthogonal unit vectors of a given sensor as
% input. The quaternions are referenced in the cell array 'data' using the
% index variable 'ind' (which corresponds to the first row of a given task)
% plus some number corresponding to the order of the sensors based on
% Matlab's numeric/alphabetical ordering scheme. After the sensor
% coordinate systems are rotated another function,
% 'RotMatTransformToGlobal' rotates the coordinate system of a given sensor
% at every time point so that it is in line with the global coordinate
% system at the beginning of the trial when all of the coordinate systems
% of all of the sensors are more or less in line with each other.
span = min([length(data{ind+num.LowThighAnt,2}(:,1)),length(data{ind+num.ShinBone,2}(:,1)), ...
    length(data{ind+num.LowShankLat,2}(:,1)),length(data{ind+num.Sacrum,2}(:,1)), ...
    length(data{ind+num.LowThighPos,2}(:,1)),length(data{ind+num.MidThighLat,2}(:,1)), ...
    length(data{ind+num.Heel,2}(:,1)), ...
    length(data{ind+num.LowThighLat,2}(:,1)),length(data{ind+num.L4L5,2}(:,1)), ...
    length(data{ind+num.MidShankLat,2}(:,1)),length(data{ind+num.DFoot,2}(:,1))]);


OffsetPoint = 10;

for aa = 1:length(FootOpt)
    
    if FootOpt(aa) == 1;
        footSupinf.Heel = data{ind+num.Heel,2}(OffsetPoint,18:20)./norm(data{ind+num.Heel,2}(OffsetPoint,18:20));
        footMedlat.Heel = Zero.Heel(2,:);
        footAntpos.Heel = cross(footSupinf.Heel,footMedlat.Heel);
        footMedlat.Heel = cross(footAntpos.Heel,footSupinf.Heel);
        
        FootCoord.Heel = CoordinateTransformQuatPreSpecifiedAxes(data{ind+num.Heel,3}(1:span,:),footSupinf.Heel,footMedlat.Heel,footAntpos.Heel);
        
        [FootCoord.Heel.Axis1,FootCoord.Heel.Axis2,FootCoord.Heel.Axis3] = RotMatTransformToGlobal( ...
            FootCoord.Heel.Axis1(:,OffsetPoint),FootCoord.Heel.Axis2(:,OffsetPoint),FootCoord.Heel.Axis3(:,OffsetPoint), ...
            FootCoord.Heel.Axis1,FootCoord.Heel.Axis2,FootCoord.Heel.Axis3);
        
    elseif FootOpt(aa) == 2;
        footSupinf.DFoot = data{ind+num.DFoot,2}(OffsetPoint,18:20)./norm(data{ind+num.DFoot,2}(OffsetPoint,18:20));
        footMedlat.DFoot = Zero.DFoot(2,:);
        footAntpos.DFoot = cross(footSupinf.DFoot,footMedlat.DFoot);
        footMedlat.DFoot = cross(footAntpos.DFoot,footSupinf.DFoot);
        
        FootCoord.DFoot = CoordinateTransformQuatPreSpecifiedAxes(data{ind+num.DFoot,3}(1:span,:),footSupinf.DFoot,footMedlat.DFoot,footAntpos.DFoot);
        
        [FootCoord.DFoot.Axis1,FootCoord.DFoot.Axis2,FootCoord.DFoot.Axis3] = RotMatTransformToGlobal( ...
            FootCoord.DFoot.Axis1(:,OffsetPoint),FootCoord.DFoot.Axis2(:,OffsetPoint),FootCoord.DFoot.Axis3(:,OffsetPoint), ...
            FootCoord.DFoot.Axis1,FootCoord.DFoot.Axis2,FootCoord.DFoot.Axis3);
    end
    
    FootSensors = fieldnames(FootCoord);
    
end

for bb = 1:length(ShankOpt)
    if ShankOpt(bb) == 1;
        shinSupinf.ShinBone=data{ind+num.ShinBone,2}(OffsetPoint,18:20)./norm(data{ind+num.ShinBone,2}(OffsetPoint,18:20));
        shinMedlat.ShinBone=Zero.ShinBone(2,:);
        shinAntpos.ShinBone = cross(shinSupinf.ShinBone,shinMedlat.ShinBone);
        shinMedlat.ShinBone = cross(shinAntpos.ShinBone,shinSupinf.ShinBone);
        
        ShankCoord.ShinBone = CoordinateTransformQuatPreSpecifiedAxes(data{ind+num.ShinBone,3}(1:span,:),shinSupinf.ShinBone,shinMedlat.ShinBone,shinAntpos.ShinBone);
        
        [ShankCoord.ShinBone.Axis1,ShankCoord.ShinBone.Axis2,ShankCoord.ShinBone.Axis3] = RotMatTransformToGlobal( ...
            ShankCoord.ShinBone.Axis1(:,OffsetPoint),ShankCoord.ShinBone.Axis2(:,OffsetPoint),ShankCoord.ShinBone.Axis3(:,OffsetPoint), ...
            ShankCoord.ShinBone.Axis1,ShankCoord.ShinBone.Axis2,ShankCoord.ShinBone.Axis3);
        
    elseif ShankOpt(bb) == 2;
        shinSupinf.MidShankLat=data{ind+num.MidShankLat,2}(OffsetPoint,18:20)./norm(data{ind+num.MidShankLat,2}(OffsetPoint,18:20));
        shinMedlat.MidShankLat=Zero.MidShankLat(2,:);
        shinAntpos.MidShankLat = cross(shinSupinf.MidShankLat,shinMedlat.MidShankLat);
        shinMedlat.MidShankLat = cross(shinAntpos.MidShankLat,shinSupinf.MidShankLat);
        
        ShankCoord.MidShankLat = CoordinateTransformQuatPreSpecifiedAxes(data{ind+num.MidShankLat,3}(1:span,:),shinSupinf.MidShankLat,shinMedlat.MidShankLat,shinAntpos.MidShankLat);
        
        [ShankCoord.MidShankLat.Axis1,ShankCoord.MidShankLat.Axis2,ShankCoord.MidShankLat.Axis3] = RotMatTransformToGlobal( ...
            ShankCoord.MidShankLat.Axis1(:,OffsetPoint),ShankCoord.MidShankLat.Axis2(:,OffsetPoint),ShankCoord.MidShankLat.Axis3(:,OffsetPoint), ...
            ShankCoord.MidShankLat.Axis1,ShankCoord.MidShankLat.Axis2,ShankCoord.MidShankLat.Axis3);
        
    elseif ShankOpt(bb) == 3;
        shinSupinf.LowShankLat=data{ind+num.LowShankLat,2}(OffsetPoint,18:20)./norm(data{ind+num.LowShankLat,2}(OffsetPoint,18:20));
        shinMedlat.LowShankLat=Zero.LowShankLat(2,:);
        shinAntpos.LowShankLat = cross(shinSupinf.LowShankLat,shinMedlat.LowShankLat);
        shinMedlat.LowShankLat = cross(shinAntpos.LowShankLat,shinSupinf.LowShankLat);
        
        ShankCoord.LowShankLat = CoordinateTransformQuatPreSpecifiedAxes(data{ind+num.LowShankLat,3}(1:span,:),shinSupinf.LowShankLat,shinMedlat.LowShankLat,shinAntpos.LowShankLat);
        
        [ShankCoord.LowShankLat.Axis1,ShankCoord.LowShankLat.Axis2,ShankCoord.LowShankLat.Axis3] = RotMatTransformToGlobal( ...
            ShankCoord.LowShankLat.Axis1(:,OffsetPoint),ShankCoord.LowShankLat.Axis2(:,OffsetPoint),ShankCoord.LowShankLat.Axis3(:,OffsetPoint), ...
            ShankCoord.LowShankLat.Axis1,ShankCoord.LowShankLat.Axis2,ShankCoord.LowShankLat.Axis3);
    end
    ShankSensors = fieldnames(ShankCoord);
    
end

for cc = 1:length(ThighOpt)
    if ThighOpt(cc) == 1;
        thighSupinf.LowThighAnt=data{ind+num.LowThighAnt,2}(OffsetPoint,18:20)./norm(data{ind+num.LowThighAnt,2}(OffsetPoint,18:20));
        thighMedlat.LowThighAnt=Zero.LowThighAnt(2,:);
        thighAntpos.LowThighAnt = cross(thighSupinf.LowThighAnt,thighMedlat.LowThighAnt);
        thighMedlat.LowThighAnt = cross(thighAntpos.LowThighAnt,thighSupinf.LowThighAnt);
        
        ThighCoord.LowThighAnt = CoordinateTransformQuatPreSpecifiedAxes(data{ind+num.LowThighAnt,3}(1:span,:),thighSupinf.LowThighAnt,thighMedlat.LowThighAnt,thighAntpos.LowThighAnt);
        
        [ThighCoord.LowThighAnt.Axis1,ThighCoord.LowThighAnt.Axis2,ThighCoord.LowThighAnt.Axis3] = RotMatTransformToGlobal( ...
            ThighCoord.LowThighAnt.Axis1(:,OffsetPoint),ThighCoord.LowThighAnt.Axis2(:,OffsetPoint),ThighCoord.LowThighAnt.Axis3(:,OffsetPoint), ...
            ThighCoord.LowThighAnt.Axis1,ThighCoord.LowThighAnt.Axis2,ThighCoord.LowThighAnt.Axis3);
        
    elseif ThighOpt(cc) == 2;
        thighSupinf.MidThighLat=data{ind+num.MidThighLat,2}(OffsetPoint,18:20)./norm(data{ind+num.MidThighLat,2}(OffsetPoint,18:20));
        thighMedlat.MidThighLat=Zero.MidThighLat(2,:);
        thighAntpos.MidThighLat = cross(thighSupinf.MidThighLat,thighMedlat.MidThighLat);
        thighMedlat.MidThighLat = cross(thighAntpos.MidThighLat,thighSupinf.MidThighLat);
        
        ThighCoord.MidThighLat = CoordinateTransformQuatPreSpecifiedAxes(data{ind+num.MidThighLat,3}(1:span,:),thighSupinf.MidThighLat,thighMedlat.MidThighLat,thighAntpos.MidThighLat);
        
        [ThighCoord.MidThighLat.Axis1,ThighCoord.MidThighLat.Axis2,ThighCoord.MidThighLat.Axis3] = RotMatTransformToGlobal( ...
            ThighCoord.MidThighLat.Axis1(:,OffsetPoint),ThighCoord.MidThighLat.Axis2(:,OffsetPoint),ThighCoord.MidThighLat.Axis3(:,OffsetPoint), ...
            ThighCoord.MidThighLat.Axis1,ThighCoord.MidThighLat.Axis2,ThighCoord.MidThighLat.Axis3);
        
    elseif ThighOpt(cc) == 3;
        thighSupinf.LowThighLat=data{ind+num.LowThighLat,2}(OffsetPoint,18:20)./norm(data{ind+num.LowThighLat,2}(OffsetPoint,18:20));
        thighMedlat.LowThighLat=Zero.LowThighLat(2,:);
        thighAntpos.LowThighLat = cross(thighSupinf.LowThighLat,thighMedlat.LowThighLat);
        thighMedlat.LowThighLat = cross(thighAntpos.LowThighLat,thighSupinf.LowThighLat);
        
        ThighCoord.LowThighLat = CoordinateTransformQuatPreSpecifiedAxes(data{ind+num.LowThighLat,3}(1:span,:),thighSupinf.LowThighLat,thighMedlat.LowThighLat,thighAntpos.LowThighLat);
        
        [ThighCoord.LowThighLat.Axis1,ThighCoord.LowThighLat.Axis2,ThighCoord.LowThighLat.Axis3] = RotMatTransformToGlobal( ...
            ThighCoord.LowThighLat.Axis1(:,OffsetPoint),ThighCoord.LowThighLat.Axis2(:,OffsetPoint),ThighCoord.LowThighLat.Axis3(:,OffsetPoint), ...
            ThighCoord.LowThighLat.Axis1,ThighCoord.LowThighLat.Axis2,ThighCoord.LowThighLat.Axis3);
        
    elseif ThighOpt(cc) == 4;
        thighSupinf.LowThighPos=data{ind+num.LowThighPos,2}(OffsetPoint,18:20)./norm(data{ind+num.LowThighPos,2}(OffsetPoint,18:20));
        thighMedlat.LowThighPos=Zero.LowThighPos(2,:);
        thighAntpos.LowThighPos = cross(thighSupinf.LowThighPos,thighMedlat.LowThighPos);
        thighMedlat.LowThighPos = cross(thighAntpos.LowThighPos,thighSupinf.LowThighPos);
        
        ThighCoord.LowThighPos = CoordinateTransformQuatPreSpecifiedAxes(data{ind+num.LowThighPos,3}(1:span,:),thighSupinf.LowThighPos,thighMedlat.LowThighPos,thighAntpos.LowThighPos);
        
        [ThighCoord.LowThighPos.Axis1,ThighCoord.LowThighPos.Axis2,ThighCoord.LowThighPos.Axis3] = RotMatTransformToGlobal( ...
            ThighCoord.LowThighPos.Axis1(:,OffsetPoint),ThighCoord.LowThighPos.Axis2(:,OffsetPoint),ThighCoord.LowThighPos.Axis3(:,OffsetPoint), ...
            ThighCoord.LowThighPos.Axis1,ThighCoord.LowThighPos.Axis2,ThighCoord.LowThighPos.Axis3);
        
    end
    
    ThighSensors = fieldnames(ThighCoord);
    
end

for dd = 1:length(TorsoOpt)
    if TorsoOpt(dd) == 1;
        torsoSupinf.Sacrum=data{ind+num.Sacrum,2}(OffsetPoint,18:20)./norm(data{ind+num.Sacrum,2}(OffsetPoint,18:20));
        torsoMedlat.Sacrum=Zero.Sacrum(2,:);
        torsoAntpos.Sacrum = cross(torsoSupinf.Sacrum,torsoMedlat.Sacrum);
        torsoMedlat.Sacrum = cross(torsoAntpos.Sacrum,torsoSupinf.Sacrum);
        
        TorsoCoord.Sacrum = CoordinateTransformQuatPreSpecifiedAxes(data{ind+num.Sacrum,3}(1:span,:),torsoSupinf.Sacrum,torsoMedlat.Sacrum,torsoAntpos.Sacrum);
        
        [TorsoCoord.Sacrum.Axis1,TorsoCoord.Sacrum.Axis2,TorsoCoord.Sacrum.Axis3] = RotMatTransformToGlobal( ...
            TorsoCoord.Sacrum.Axis1(:,OffsetPoint),TorsoCoord.Sacrum.Axis2(:,OffsetPoint),TorsoCoord.Sacrum.Axis3(:,OffsetPoint), ...
            TorsoCoord.Sacrum.Axis1,TorsoCoord.Sacrum.Axis2,TorsoCoord.Sacrum.Axis3);
        
    elseif TorsoOpt(dd) == 2;
        torsoSupinf.L4L5 = data{ind+num.L4L5,2}(OffsetPoint,18:20)./norm(data{ind+num.L4L5,2}(OffsetPoint,18:20));
        torsoMedlat.L4L5 = Zero.L4L5(2,:);
        torsoAntpos.L4L5 = cross(torsoSupinf.L4L5,torsoMedlat.L4L5);
        torsoMedlat.L4L5 = cross(torsoAntpos.L4L5,torsoSupinf.L4L5);
        
        TorsoCoord.L4L5 = CoordinateTransformQuatPreSpecifiedAxes(data{ind+num.L4L5,3}(1:span,:),torsoSupinf.L4L5,torsoMedlat.L4L5,torsoAntpos.L4L5);
        
        [TorsoCoord.L4L5.Axis1,TorsoCoord.L4L5.Axis2,TorsoCoord.L4L5.Axis3] = RotMatTransformToGlobal( ...
            TorsoCoord.L4L5.Axis1(:,OffsetPoint),TorsoCoord.L4L5.Axis2(:,OffsetPoint),TorsoCoord.L4L5.Axis3(:,OffsetPoint), ...
            TorsoCoord.L4L5.Axis1,TorsoCoord.L4L5.Axis2,TorsoCoord.L4L5.Axis3);
        
    end
    
    TorsoSensors = fieldnames(TorsoCoord);
    
end
%% Set up plotting colors and line styles
% These variables are used later in plotting to generate unique plots for
% every sensor combination. Each sensor above the joint is given a unique
% color and every sensor below the joint is given a unique line style.
colors = {'r','b','g','k'};
lines = {'-','--','-.',':'};

%% Hip Angles
% The hip Euler angles are calculated using the 'HipEulerAngle' function
% and plotted. This is done for every possible sensor combination on the
% same figure. A legend is provided listing every combination.
if HipOpt == 1
    figure
    HipLegend = cell(length(TorsoSensors)*length(ThighSensors),1);
    legInd = 1;
    
    
    for ee = 1:length(TorsoSensors)
        for ff = 1:length(ThighSensors)
            
            ToC = TorsoCoord.(TorsoSensors{ee});
            ThC = ThighCoord.(ThighSensors{ff});
            
            HipLegend{legInd} = strcat('Torso',32,TorsoSensors{ee},32,'/',32,'Thigh',32,ThighSensors{ff});            
            
            [flexangHip{ee,ff},abdangHip{ee,ff},rotangHip{ee,ff}] = AnglesYXZ(ToC.Axis1,ToC.Axis2,ToC.Axis3,ThC.Axis1, ...
                ThC.Axis2,ThC.Axis3,[0 -1 0],[0 0 1]);
            
            Hip_Names{ee,ff} = HipLegend{legInd};            
            legInd = legInd + 1;
            
      
            %   Flexion
            subplot(3,1,1)
            hold on
            plot(linspace(0,(1/FR)*length(flexangHip{ee,ff}),length(flexangHip{ee,ff})),flexangHip{ee,ff},strcat(colors{ee},lines{ff}))
            title(strcat('Hip Angles',10,'Flexion'))
            ylabel('Degrees')
            hold off
            %   Abduction
            subplot(3,1,2)
            hold on
            plot(linspace(0,(1/FR)*length(abdangHip{ee,ff}),length(abdangHip{ee,ff})),abdangHip{ee,ff},strcat(colors{ee},lines{ff}))
            title('Abduction')
            ylabel('Degrees')
            hold off
            %   Rotation
            subplot(3,1,3)
            hold on
            plot(linspace(0,(1/FR)*length(rotangHip{ee,ff}),length(rotangHip{ee,ff})),rotangHip{ee,ff},strcat(colors{ee},lines{ff}))
            title('Rotation')
            ylabel('Degrees')
            hold off


        end
    end

    hold on
    subplot(3,1,1)
    legend(HipLegend)
    hold off
 
    
end
%% Knee Angles
% The knee Euler angles are calculated using the 'KneeEulerAngle' function
% and plotted. This is done for every possible sensor combination on the
% same figure. A legend is provided listing every combination.
if KneeOpt == 1
    figure
    KneeLegend = cell(length(ThighSensors)*length(ShankSensors),1);
    legInd = 1;
    for ee = 1:length(ShankSensors)
        for ff = 1:length(ThighSensors)

            ShC = ShankCoord.(ShankSensors{ee});
            ThC = ThighCoord.(ThighSensors{ff});
            
            KneeLegend{legInd} = strcat('Shank',32,ShankSensors{ee},32,'/',32,'Thigh',32,ThighSensors{ff});
         
             [flexangKnee{ee,ff},abdangKnee{ee,ff},rotangKnee{ee,ff}] = AnglesYXZ(ThC.Axis1,ThC.Axis2,ThC.Axis3, ...
                ShC.Axis1,ShC.Axis2,ShC.Axis3,[0 -1 0],[0 0 1]);
                 
            Knee_Names{ee,ff} = KneeLegend{legInd};
            legInd = legInd + 1;
            %   Flexion
            subplot(3,1,1)
            hold on
            plot(linspace(0,(1/FR)*length(flexangKnee{ee,ff}),length(flexangKnee{ee,ff})),flexangKnee{ee,ff},strcat(colors{ee},lines{ff}))
            title(strcat('Knee Angles',10,'Flexion'))
            ylabel('Degrees')
            hold off
            %   Abduction
            subplot(3,1,2)
            hold on
            plot(linspace(0,(1/FR)*length(abdangKnee{ee,ff}),length(abdangKnee{ee,ff})),abdangKnee{ee,ff},strcat(colors{ee},lines{ff}))
            title('Abduction')
            ylabel('Degrees')
            hold off
            %   Rotation
            subplot(3,1,3)
            hold on
            plot(linspace(0,(1/FR)*length(rotangKnee{ee,ff}),length(rotangKnee{ee,ff})),rotangKnee{ee,ff},strcat(colors{ee},lines{ff}))
            title('Rotation')
            ylabel('Degrees')
            hold off
        end
    end
    
    hold on
    subplot(3,1,1)
    legend(KneeLegend)
    hold off
    
end
%% Ankle Angles
% The ankle Euler angles are calculated using the 'AnkleEulerAngle'
% function and plotted. This is done for every possible sensor combination
% on the same figure. A legend is provided listing every combination.
if AnkleOpt == 1
    figure
    AnkleLegend = cell(length(ShankSensors)*length(FootSensors),1);
    legInd = 1;
    for ee = 1:length(ShankSensors)
        for ff = 1:length(FootSensors)

            ShC = ShankCoord.(ShankSensors{ee});
            FoC = FootCoord.(FootSensors{ff});
            
            AnkleLegend{legInd} = strcat('Shank',32,ShankSensors{ee},32,'/',32,'Foot',32,FootSensors{ff});
             
            [flexangAnkle{ee,ff},abdangAnkle{ee,ff},rotangAnkle{ee,ff}] = AnglesYXZ(ShC.Axis1,ShC.Axis2,ShC.Axis3, ...
                FoC.Axis1,FoC.Axis2,FoC.Axis3,[0 -1 0],[0 0 1]);
                   
            Ankle_Names{ee,ff} = AnkleLegend{legInd};
            legInd = legInd + 1;
            %   Flexion
            subplot(3,1,1)
            hold on
            plot(linspace(0,(1/FR)*length(flexangAnkle{ee,ff}),length(flexangAnkle{ee,ff})),flexangAnkle{ee,ff},strcat(colors{ee},lines{ff}))
            title(strcat('Ankle Angles',10,'Flexion'))
            ylabel('Degrees')
            hold off
            %   Abduction
            subplot(3,1,2)
            hold on
            plot(linspace(0,(1/FR)*length(abdangAnkle{ee,ff}),length(abdangAnkle{ee,ff})),abdangAnkle{ee,ff},strcat(colors{ee},lines{ff}))
            title('Abduction')
            ylabel('Degrees')
            hold off
            %   Rotation
            subplot(3,1,3)
            hold on
            plot(linspace(0,(1/FR)*length(rotangAnkle{ee,ff}),length(rotangAnkle{ee,ff})),rotangAnkle{ee,ff},strcat(colors{ee},lines{ff}))
            title('Rotation')
            ylabel('Degrees')
            hold off
        end
    end
    
    hold on
    subplot(3,1,1)
    legend(AnkleLegend)
    hold off
    
end
%% Save Angle Data
saveit = input('Do you want to save the calculated angles to a .mat file? (1 yes/ 0 no)');
if saveit == 1
    saveInd = 1;
    vars2save = {};
    vars2check = {'abdangAnkle','abdangKnee','abdangHip','rotangAnkle','rotangKnee', ...
        'rotangHip','flexangAnkle','flexangKnee','flexangHip','Ankle_Names', ...
        'Knee_Names','Hip_Names'};
    for qq = 1:12
        if exist(vars2check{qq},'var')
            vars2save{saveInd} = vars2check{qq};
            saveInd = saveInd + 1;
        end
    end
    cd(strcat(dir,'/IMU Angles'))
    uisave(vars2save)
    cd(dir)
end