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
%{
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
%} 

% visual angle in reference to middle of the screen
for iTrial = 1:size(data.trial,2)
visAngleX{iTrial} = rad2deg(atan((xGaze{iTrial}/2)/(totDist*screenResX/screenWidth_x)))*2;
visAngleY{iTrial} = rad2deg(atan((yGaze{iTrial}/2)/(totDist*screenResY/screenWidth_y)))*2;
end

% find onset and offset latencies of sacades with ft_detect_movement
cfg=[];
cfg.channel = {'UADC005', 'UADC006'};
dataEye = ft_selectdata(cfg, data);
cfg                     = [];
cfg.method              = 'velocity2D';
cfg.velocity2D.mindur   = 14; % minimum microsaccade durantion in samples (default = 3) (=12ms Engbert & Kliegl)
cfg.velocity2D.velthres = 6; % threshold for velocity outlier detection (default = 6);
[~, movement]           = ft_detect_movement(cfg,dataEye);

% movement contains saccade on- and offset samples. find the corresponding
% visual angles and compute the difference (for x and y). save these and
% select trial as containing eyemovement if visualangle>1
trlIdx = zeros(1,size(movement, 1));
saccadeAngle = zeros(2,size(movement, 1));
artfctdef.eye.saccade.artifact = [];

for iSaccade = 1:size(movement, 1)
    % find the trial where the saccade belongs to
    trlIdx(iSaccade) = find(dataEye.sampleinfo(:,1)<=movement(iSaccade,1) & dataEye.sampleinfo(:,2)>=movement(iSaccade,2));
    % how many samples after the start of the trial does the saccade take
    % place and when does it end?
    saccadeStart = movement(iSaccade,1) - dataEye.sampleinfo(trlIdx(iSaccade),1) + 1;
    saccadeEnd = movement(iSaccade,2) - dataEye.sampleinfo(trlIdx(iSaccade),1) + 1;
    % The difference in the visual angle w.r.t. fixation at the start and 
    % end of the saccade is the saccade visual angle (column 1 horizontal,
    % column 2 vertical).
    saccadeAngle(1, iSaccade) = abs(visAngleX{trlIdx(iSaccade)}(saccadeEnd) - visAngleX{trlIdx(iSaccade)}(saccadeStart));
    saccadeAngle(2, iSaccade) = abs(visAngleY{trlIdx(iSaccade)}(saccadeEnd) - visAngleY{trlIdx(iSaccade)}(saccadeStart));
    if saccadeAngle(:, iSaccade)>1 % otherwise it is a microsaccade.
        artfctdef.eye.saccade.artifact(end+1,:) = movement(iSaccade, 1:2);
    end
end

