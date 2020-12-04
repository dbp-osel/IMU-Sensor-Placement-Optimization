function [files] = IMUImport(location,column1,row1,TaskList,num,NumberSensors)
%Imports the IMU files and puts them in a structure containing file name,
%corresponding data, and quaternions. The input (location) is the
%directory where the IMU files are stored. Note that this function will read
%any text file so be sure only IMU files are in this folder. The output
%(files) is a cell array containing the file name in the first column. The
%second column contains all of the data in each file as is in the form of a
%matrix. The third column of the cell array contains quaternions. The
%starting column for the quaternions is given by the column1 function
%input. The starting row for all the data is given by row1.

%Inputs
%location: The directory of the folder containing the text files from Xsens
%   type: char
%column1: The first column in the text file data containing quaternions
%   type: 1x1 double
%row1: The first row in the text file data containing numeric data
%   type: 1x1 double
%TaskList: The list of tasks contained in the folder
%   type: char
%Num: Structure containing sensor names and the orders they appear in the
%folder
%   type: struct
%NumberSensors: The number of sensors used.
%   type: 1x1 double

%Output
%files: Contains text file data, quaternions, sensor names, task names, and
%file names. Each row corresponds to a single file.
%   type: cell

direct = pwd;
directData = strcat(location);
cd(directData);
filenames = dir('*.txt');
files = cell(length(filenames),4);

for i = 1:length(filenames(:,1))
    fopen(filenames(i).name);
    
    try
        data_temp = dlmread(filenames(i).name,'	',row1,0);
        
        Quat_temp = [data_temp(:,column1),data_temp(:,column1+1),data_temp(:,column1+2),data_temp(:,column1+3)];
        
        files(i,1:3) = [{filenames(i).name},{data_temp},Quat_temp, []];
        clear data_temp Quat_temp
    catch
        disp(strcat('Unable to read',filenames(i).name,' with dlmread',10,'Check file contents'))
        %Some of the files may be empty if the recording didn't start right
        %or if there is some other error in the file. This try catch
        %statement makes sure that the entire code doesn't bug out if one
        %trial is faulty.
    end
    fclose('all');
    %If fclose isn't in the loop then too many files will open causing
    %Matlab to bug out. Probably a memory issue.

end
fclose('all');
cd(direct);

newlines = find(TaskList == 10);
newlines = [1 newlines length(TaskList)];

taskNum = 1;
taskCount = 0;

SensorNames = fieldnames(num);
SensorNums = zeros(length(SensorNames),1);

for z = 1:length(SensorNames)
    SensorNums(z) = eval(strcat('num.',SensorNames{z}));
end
for q = 1:length(filenames(:,1))
    try
        files(q,4:5) = [{TaskList(newlines(taskNum):newlines(taskNum+1))} ...
            {SensorNames{find(SensorNums == taskCount)}}];
        taskCount = taskCount + 1;
        if taskCount == NumberSensors
            taskCount = 0;
            taskNum = taskNum + 1;
        end
    catch
        
    end
end
end
