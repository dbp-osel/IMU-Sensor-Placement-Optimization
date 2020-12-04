function [flexang,abdang,rotang] = AnglesYXZ(UpperCoord1,UpperCoord2,UpperCoord3,LowerCoordSI,LowerCoordML,LowerCoordAP,ML,AP,varargin)
%Inputs
%UpperCoord 1,2, and 3 specify 3xn matrices where columns represent unit
%vectors pointing along the three axes of the upper sensor. LowerCoord SI,
%ML and AP specify 3xn matrices of unit vectors along the three axes of the
%lower sensor. ML and AP specify the directions of the ML and AP axes
%during calibration. varargin has 2 possible inputs. If 'wrap' is specified
%then angle wrapping correction above or below +/-180 degrees will occur.
%If 'threshold' is specified as a varargin input then the next varargin
%input should be a numeric value. If the slope of the angle wrapping
%"jumps" is greater than this value then the correction will be applied.
%The default value of this threshold is 50.
%Outputs
%flexion extension angle, abduction/adduction angle, internal/external
%rotation angle in that order. A YXZ rotation order is used.
       
slopeThreshold = 50;
Wrapit = false;
for qq = 1:length(varargin)
if strcmp(varargin{qq},'wrap')
    Wrapit = true;
elseif strcmp(varargin{qq},'threshold')
    slopeThreshold = varargin{qq+1};
end
end
Vect1 = DCMProjection(LowerCoordSI,UpperCoord1,UpperCoord2,UpperCoord3);
Vect2 = DCMProjection(LowerCoordML,UpperCoord1,UpperCoord2,UpperCoord3);
Vect3 = DCMProjection(LowerCoordAP,UpperCoord1,UpperCoord2,UpperCoord3);


%% First Euler rotation about Y
% The rotation about Y is described as the planer rotation in the pelvis
% sensor Z-X plane. As such the angle described by the components of the
% thigh sensor's z-axis (down) on the pelvis sensor's x and z axes is found.

flexang = atan2d((Vect1(3,:)),(Vect1(1,:)));

%% Second Euler Rotation about X' A second rotation about the X' axis can be
% found in the z'-y plane (same as z'-y' plane as y didn't move). The
% component along z' can be found from the components in the original z-x
% plane using the pythagorean theorem. The second Euler angle is found
% between the z' component and the y component (determined earlier from the
% projection). Essentially we are finding the rotation in the original x-z
% plane, then finding a second rotation in the plane defined by the axis
% created by the first rotation and the original y axis.

abdang = atan2d(Vect1(2,:),sqrt(sum(Vect1([1,3],:).^2,1)));

% abdang = atand(sqrt(sum(RupaVect1([1,3],:).^2,1))./RupaVect1(2,:));
%using pythagorean theorem is valid for the X' component since it
%defines the positive direction and cannot be negative in this reference
%frame

if Wrapit == true
    [flexang] = WrapOpt(flexang,slopeThreshold);
    [abdang] = WrapOpt(abdang,slopeThreshold);
end

  %% Third Euler Angle about Z''
  %Quaternion describing rotation about Y
  onez = ones(1,length(Vect1(1,:)));
  zeroz = zeros(1,length(Vect1(1,:)));
  
  quatFlex = [cosd(flexang./2);sind(flexang./2).*ML(1);sind(flexang./2).*ML(2);sind(flexang./2).*ML(3);];
  %Rotation about Y (ML axis)
  
  %Use quatFlex to rotate X axis to find X'
  Xprime = QuaternionConjugation(quatFlex',[zeroz;onez*AP(1);onez*AP(2);onez*AP(3);]');
  Xprime(:,1) = [];
  
  XprimeProj = DCMProjection(Xprime',Vect1,Vect2,Vect3);

  rotang = atan2d(XprimeProj(2,:),XprimeProj(3,:));

  if Wrapit == true
      [rotang] = WrapOpt(rotang,slopeThreshold);
  end
%% Wrap correction
    function [ang] = WrapOpt(ang,slopeThreshold)
        
        % Sometimes angle wrapping occurs since Matlab's atan2d function only
        % produces angles between +/- 180 degrees. A slope threshold of
        % between 50 and 100 seems to work but you'll probably have to play
        % with it. 50 is the default if you don't specify an option.

        angSlope = [0 ang(2:end)-ang(1:end-1)];
        
        
        if max(angSlope >= slopeThreshold)
            %The derivitive is calculated
            
            [~,poslocs] = findpeaks(angSlope,'MinPeakHeight',slopeThreshold);
            [~,neglocs] = findpeaks(-1*angSlope,'MinPeakHeight',slopeThreshold);
            %findpeaks locates spikes in the slope function indicitive of angle wrapping
            
            if ~isempty(poslocs)
                for ii = 1:length(poslocs)
                    ang(poslocs(ii):end) =  ang(poslocs(ii):end)-360;
                end
            end
            if ~isempty(neglocs)
                for ii = 1:length(neglocs)
                    ang(neglocs(ii):end) =  ang(neglocs(ii):end)+360;
                end
            end
            % 360 degrees is added or subtracted for subsequent points every time angle
            % wrapping occurs depending on the direction of wrapping. If the points
            % shift down, 360 degrees is added. If the points shift up, 360 degrees is
            % subtracted
        end
    end
end
