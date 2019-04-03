function [visAngleX, visAngleY, visAngPerPixX, visAngPerPixY] = transform_eyedata(data)
% pixel coords for 1920x1080 screen, upper left is (0,0)
    % these dimensions can be found in PHYSICAL.INI or your presentation/PTB
    % settings.
    screenRight=1919;
    screenLeft=0;
    screenTop=0;
    screenBottom=1079;
    
    % range of voltages and recording range as defined in FINAL.INI
    minVoltage=-5;
    maxVoltage=5;
    minRange=0;
    maxRange=1;
    
    totDist=780; % distance eye-screen in mm
    
    % if you want to calculate the visual degree per pixel, physical dimensions
    % of the pixels have to be known (so projected screen size).
    screenWidth_y=270;
    screenWidth_x=480; % has to be 'projected width' in mm (not actual width)
    
    % convert ADC-voltref back to GAZE position
    xGaze=[];
    yGaze=[];
    for iTrial=1:size(data.trial,2);% these are trials with right-left voluntary saccades
        
        
        idxChanHor=find(strcmp(data.label,'UADC005'));
        idxChanVer=find(strcmp(data.label,'UADC006'));
        
        voltageHor=data.trial{iTrial}(idxChanHor,:);
        voltageVer=data.trial{iTrial}(idxChanVer,:);
        
        % see eyelink1000_UserManual for formulas:
        RangeHor = (voltageHor-minVoltage)./(maxVoltage-minVoltage); % voltage range proportion
        S_h = RangeHor.*(maxRange-minRange)+minRange; % proportion of screen width or height
        
        R_v = (voltageVer-minVoltage)./(maxVoltage-minVoltage);
        S_v = R_v.*(maxRange-minRange)+minRange;
        
        xGaze{iTrial} = S_h.*(screenRight-screenLeft+1)+screenLeft;
        yGaze{iTrial} = S_v.*(screenBottom-screenTop+1)+screenTop;
        
        
        % gaze position on center of screen (where fixation dot is)
        Xgaze_origin=xGaze;
        Ygaze_origin=yGaze;
        
        xGaze{iTrial} = xGaze{iTrial} - round(screenRight/2);%%presentation/PTB assumes 1 more pixel than there actually is.
        %presenation x-position = screen position + 1 pixel to the right:
        %presentation x-pos(0) = 1919/2
        
        yGaze{iTrial} = yGaze{iTrial} - round(screenBottom/2);
    end
    % Xgaze and Ygaze now are the eye positions on the screen (relative to
    % fixation point). We can now convert them to visual angles (relative to
    % fixation point).
    
    % if you want to know the visual angle per pixel, calculate this
    %
    % X
    screenResX=screenRight+1;
    visAngX_rad = 2 * atan(screenWidth_x/2/totDist);
    visAngX_deg = visAngX_rad * (180/pi);
    visAngPerPixX = visAngX_deg / screenResX;
    
    % Y
    screenResY=screenBottom+1;
    visAngY_rad = 2 * atan(screenWidth_y/2/totDist);
    visAngY_deg = visAngY_rad * (180/pi);
    visAngPerPixY = visAngY_deg / screenResY;
    
    % visual angle in reference to middle of the screen
    for iTrial = 1:size(data.trial,2)
        visAngleX{iTrial} = rad2deg(atan((xGaze{iTrial}/2)/(totDist*screenResX/screenWidth_x)))*2;
        visAngleY{iTrial} = rad2deg(atan((yGaze{iTrial}/2)/(totDist*screenResY/screenWidth_y)))*2;
    end