function concentric_grating_experiment(fileName, isLive)
% This function is copied from DriftDemo and modified to display a (masked)
% animated concentric grating moving inward. At random time points, the
% grating will make a 180 degree phase shift, which has to be reported by
% the subject by a button press. 10% of the trials have no phase shift.

% NOTE: costs a lot of memory because all images are made individually

% _________________________________________________________________________
%
% This is a very simple, bare bones demo on how to do frame animation. For
% much more efficient ways to draw gratings and gabors, have a look at
% DriftDemo2, DriftDemo3, DriftDemo4, ProceduralGaborDemo, GarboriumDemo,
% ProceduralGarboriumDemo and DriftWaitDemo.

% NOTE: DriftDemo2 is more efficient because only one texture is
% generated, which is moved every frame to create the illusion of movement.
% However, this is not possible for a concentric grating, so every frame a
% new texture has to be generated.
% _________________________________________________________________________

% HISTORY
% 7/6/2016: change grating to concentric grating
% 8/6/2016: add circular and gaussian mask, introduce phase-shift
% 9/6/2016: make trials, introduce baseline period
% 10/6/2016: add no shift trials
% 15/5/2016: change background color
% 16/6/2016: edit Gaussian mask, add (rectangular) hanning mask
% 16/6/2016: add circular hanning mask
% 17/6/2016: debugged circular hanning mask and response buttons
% 20/6/2016: change location to lower visual field, lef t and right.
% 21/6/2016: edit registereing of user input, record timing, edit locations
%   of gratings (central, 2 peripheral) and introduce trial
%   counterbalancing (semirandom instead of random locations)
% 22/6/2016: solve bug in hanning mask. Before, only [gray, white] values
%   were dampened, not [black, gray] values. Because the grating
%   fluctuated around gray when the mask was applied. Now it fluctuates
%   around 0 (it is symmetric, [-gray,gray]), and afterwards the grating is
%   moved to the [black, white] interval.
%   solved bug in saving response
% 22/6/2016: introduce Bitsi responses, solve bugs in saving responses.
% 23/6/2016: changed naming convention
% 28/6/2016: deleted offset trigger, edited bitsi button response

%% Function settings


if nargin<1
    fileName = 'test';
end

if isempty(fileName)
    fileName = 'test';
end

if nargin<2
    isLive = 0;
end

if isempty(isLive)
    isLive = 0;
end

try
    %% Standard PTB settings
    % This script calls Psychtoolbox commands available only in OpenGL-
    % based versions of the Psychtoolbox. (So far, the OS X Psychtoolbox is
    % the only OpenGL-base Psychtoolbox.)  The Psychtoolbox command 
    % AssertPsychOpenGL will issue an error message if someone tries to 
    % execute this script on a computer without an OpenGL Psychtoolbox
    AssertOpenGL;
    
    % Get the list of screens and choose the one with the highest screen number.
    % Screen 0 is, by definition, the display with the menu bar. Often when
    % two monitors are connected the one without the menu bar is used as
    % the stimulus display.  Chosing the display with the highest dislay number is
    % a best guess about where you want the stimulus displayed.
    screens=Screen('Screens');
    if isLive
        screenNo=max(screens);
    else
        screenNo=1;
    end
    
    % Find the color values which correspond to white and black: Usually
    % black is always 0 and white 255, but this rule is not true if one of
    % the high precision framebuffer modes is enabled via the
    % PsychImaging() commmand, so we query the true values via the
    % functions WhiteIndex and BlackIndex:
    white=WhiteIndex(screenNo);
    black=BlackIndex(screenNo);
    
    % Round gray to integral number, to avoid roundoff artifacts with some
    % graphics cards:
    gray=round((white+black)/2);
    
    % This makes sure that on floating point framebuffers we still get a
    % well defined gray. It isn't strictly neccessary in this demo:
    if gray == white
        gray=white / 2;
    end
    
    % Contrast 'inc'rement range for given white and gray values:
    inc=white-gray;
    
    % Open a double buffered fullscreen window and select a black background
    % color:
    [window, windowRect]=Screen('OpenWindow',screenNo, gray);
    
    % Get the centre coordinate of the window
    [xCenter, yCenter] = RectCenter(windowRect);
    
    % Query the frame duration
    ifi = Screen('GetFlipInterval', window);
    frameRate=Screen('FrameRate',screenNo); %1/ifi
    
    % Set up alpha-blending for smooth (anti-aliased) lines
    Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
    
    
    HideCursor;
    
    
    %% Trigger settings
    if isLive
        btsi = Bitsi('com1');% when there is no com connected, use keyboard as input
    else
        btsi = Bitsi('');
    end
    
    xp = expconstants();
    
    %% Input settings
    
    KbName('UnifyKeyNames'); % converts operating specific keynames to universal keyname
    escape = KbName('ESCAPE');
    space = KbName('space');
    
    log=[]; %save trial latencies, response and response time
    
    %% Experimental settings
    % calculate the size of the stimulus according to the used setup
    if isLive
        screen.width = 52; %MEG screenwidth in cm
        screen.viewingDistance = 76; % in cm
    else
        screen.width = 52;%desktop screenwidth in cm
        screen.viewingDistance = 59; % in cm
    end
    screen.resolutionX = 1920; %desktop and projector screen in pixels
    screen.resolutionY = 1080;
    
    screen.totVisDeg = 2*atan(screen.width/(2*screen.viewingDistance))*(180/pi); % formula for calculating visual degrees
    screen.pixPerDeg = screen.resolutionX/screen.totVisDeg;
    
    visualAngleGrating = 12.5;
    visualAngleLocation = 15;
    gratingSize = visualAngleGrating*screen.pixPerDeg; % calculate how big 
    % the grating should be (can not be rounded off, might result in odd integer)
    gratingRadius = round(gratingSize/2); % the grating can only exist of integers, so round
    gratingSize = 2*gratingRadius; % to prevent consistency errors, redifine gratingSize
    rLocation = round(visualAngleLocation*screen.pixPerDeg/2);
    
    instructions = 'Instructions \n\n Fixate at the middle of the screen. \n\n Press the button if you see a flash in the animation';
    nTrialsPerCondition = 150;
    nConditions=3;
    nTrials=nConditions*nTrialsPerCondition;
    conditions = [-ones(nTrialsPerCondition,1); zeros(nTrialsPerCondition,1);...
        ones(nTrialsPerCondition,1)];
    conditions = conditions(randperm(size(conditions,1)),:); % shuffle conditions
    
    chanceNoShift = 0.1; % in 20% of the trials, have no shift.
    preStimTime=1.5; % fixation, baseline period
    preShiftTime=1; % cue can come 1 seconds after start trial
    shiftRange = 2; % present the shift at random in 2 seconds
    postShiftTime = 0.75; % keep animation going 1 seconds after shift
    driftFreq = 1; % in one second, every pixel of the grating completes one cylce (black-white-black)
    maxRespTime = 0.7; % maximum response time after shift
    nFramesInCycle=round(driftFreq/ifi); % temporal period, in frames, of the drifting grating
    
    % settings grating
    [x,y]=meshgrid(-gratingRadius:gratingRadius,-gratingRadius:gratingRadius);
    f=1.5*2*pi; % period of the grating.
    
    
    %% Generate stimulus
    
    % make circular hanning mask
    L = 2*gratingRadius+1;
    w1D = hann(L); % 1D hann window
    xx = linspace(-gratingRadius,gratingRadius,L);
    [X,Y] = meshgrid(xx);
    r = sqrt( X.^2 + Y.^2 );
    w2D = zeros(L);
    w2D(r<=gratingRadius) = interp1(xx,w1D,r(r<=gratingRadius)); % create
    % 2D hanning window.
    
    % Generate grating texture
    % Compute each frame of the movie and convert those frames stored in
    % MATLAB matrices, into Psychtoolbox OpenGL textures using 'MakeTexture';
    % NOTE: very inefficient, costs a lot of memory.. alternative?
    tex=zeros(nFramesInCycle,1);
    
    for jFrame=1:nFramesInCycle
        phase=(jFrame/nFramesInCycle)*2*pi; %change the phase of the grating according to the framenumber
        m=sin(sqrt(x.^2+y.^2) / f + phase);
        grating = (w2D.*(inc*m)+gray);
        % inc*m fluctuates from [-gray, gray]. Multiply this with the
        % hanning mask to let the grating die off at 0. Now add gray to let
        % the grating fluctuate from [black, white], converging at gray.
        tex(jFrame)=Screen('MakeTexture', window, grating); % multiply mask with grating
    end
    
    % set location gratings
    gratingDim = [0 0 2*gratingRadius 2*gratingRadius];
    gratingYpos = zeros(nTrials,1);
    for iTrl=1:nTrials % gratings in periphery are presented in lower visual
        % field, central stimuli are central.
        if conditions(iTrl,1)==1 || conditions(iTrl,1)==-1
            gratingYpos(iTrl,1) = yCenter+0.25*screen.resolutionY;
        else
            gratingYpos(iTrl,1) = yCenter;
        end
    end
    gratingXpos = xCenter + conditions*rLocation;
    log.Xpos = conditions;
    
    shiftLatency = ((preShiftTime+shiftRange)-preShiftTime).*rand(nTrials,1)...
        + preShiftTime; % generate latencies at which a shift occurs for every trial
    
    trlNoShift = sort(randperm(nTrials, nTrials*chanceNoShift)); % generate a number of random (w/o replacement) trials that won't have a phase shift
    shiftFrame = round(shiftLatency*frameRate);
    log.trlNoShift = trlNoShift;
    
    KbName('UnifyKeyNames'); % convert operating specific keynames to general keynames
    escape = KbName('ESCAPE');
    space = KbName('space');
    
    %% Start Experiment
    
    % Run the movie animation for a variable period: end it 1.5 seconds
    % after phase shift.
    movieDurationSecs=shiftLatency+postShiftTime;
    % Convert movieDuration in seconds to duration in frames to draw:
    nFramesTotal=round(movieDurationSecs * frameRate);
    % assign the right texture indices to frames
    for iTrl=1:nTrials
        if iTrl==1
            % Provide instructions at the beginning of the experiment
            Screen('FillRect', window, gray);
            DrawFormattedText(window, instructions,...
                'center', 'center', white, 55, [], [], 1.3);
            Screen('Flip', window);
            KbWait; % wait until key is pressed
            pause(0.5)
        end
        frameTexId=mod(0:(nFramesTotal(iTrl)-1), nFramesInCycle) + 1; %assign the right texture index to each frame
        
        position = CenterRectOnPointd(gratingDim, gratingXpos(iTrl), ...
                    gratingYpos(iTrl)); %move the object of size gratingDim to those coordinates
        
        
        % Use realtime priority for better timing precision:
        priorityLevel=MaxPriority(window);
        Priority(priorityLevel);
        
        % Baseline period, gray screen with fixation dot
        Screen('DrawDots', window, [xCenter yCenter], 20, [255 0 0], [], 2);
        btsi.sendTrigger(xp.TRIG_ONSET_TRIAL);
        Screen('Flip', window)
        pause(preStimTime)
        
        [~, ~, keyCode] = KbCheck(-3); % all keyCodes should be zero at start of trial.
        flag.prevResp=0; % reset flag
        btsi.clearResponses(); % reset Bitsi response
        btsi.sendTrigger(xp.TRIG_ONSET_GRATING);
        t0= GetSecs; % time onset grating
        
        % Animation loop:
        for jFrame=1:nFramesTotal(iTrl)
            
            isShiftFrame = (jFrame==shiftFrame(iTrl));
            if isShiftFrame
                t1=GetSecs; % get time stamp at the moment of phase shift in order to calculate reaction time
            end
            
            % at the moment the phase-shift should occur, jump to the texture
            % with a 180 degree phase shift, and add this shift for the
            % remainder of the video. Otherwise, don't shift.
            isPostShift = (jFrame>=shiftFrame(iTrl));
            isPostHalfCycle = (frameTexId(jFrame)>nFramesInCycle/2);
            isPreHalfCycle = (frameTexId(jFrame)<=nFramesInCycle/2);
            isNoShiftTrl = any(trlNoShift==iTrl);
            if isPostShift && isPostHalfCycle && ~isNoShiftTrl
                shift = -nFramesInCycle/2; % if more than half of the period already passed, subtract the shift
                % enable abort using esc
            elseif isPostShift && isPreHalfCycle && ~isNoShiftTrl
                shift = nFramesInCycle/2; % otherwise add the shift (because the period loops)
            else
                shift=0;
            end
            
            % Draw the appropriate image
            Screen('DrawTexture', window, tex(frameTexId(jFrame) + shift),...
                [], position);
            Screen('DrawDots', window, [xCenter yCenter], 20, [255 0 0], [], 2);
            
            % NOTE
            % Show it at next display vertical retrace. Please check DriftDemo2
            % and later, as well as DriftWaitDemo for much better approaches to
            % guarantee a robust and constant animation display timing! This is
            % very basic and not best practice!
            
            if isShiftFrame;
                btsi.sendTrigger(xp.TRIG_SHIFT);
            end
            Screen('Flip', window); %% Flipping the screen takes up to 30ms every frame.. 1.8s every cycle..
            
            
            %% Save response
            
            % check whether there already was a response in a previous
            % frame. If so, flag it so new responses won't be recorded
            % again.
            if keyCode(space)
                flag.prevResp=1;
            end
            [~, ~, keyCode] = KbCheck(-3);
            
            % check whether experiment has to be aborted
            if keyCode(escape) % end the experiment
                save(fileName, 'log');
                sca;
                Screen('CloseAll')
            end
            
            isWithinResponseTime = (jFrame/frameRate<(shiftLatency(iTrl)...
                                   + maxRespTime));
            if isLive
                [resp, secs] = btsi.getResponse(0.0001, true);
                % ...getResponse(x,[]), x should be as small as possible,
                % because this step will wait maximal x seconds or until
                % there is a response.
                if resp == 'd' ...% right red button down event (left button down event is 'h')
                        && ~flag.prevResp && isWithinResponseTime % no previous response and within response time
                    if jFrame>shiftFrame(iTrl) % response only correct after shift
                        log.response(iTrl) = 1;
                        log.responseTime(iTrl) = secs-t1;
                    else % response before shift: incorrect trial
                        log.response(iTrl) = 2;
                        log.responseTime(iTrl) = 0;
                    end
                elseif jFrame==nFramesTotal(iTrl) && ~flag.prevResp
                    log.response(iTrl) = 0;
                    log.responseTime(iTrl)=0;
                end
            else
                [~, secs, keyCode] = KbCheck(-3);
                if keyCode(space) && ~flag.prevResp &&  isWithinResponseTime % no previous response and within response time
                    if jFrame>shiftFrame(iTrl) % response only correct after shift
                        log.response(iTrl) = 1;
                        log.responseTime(iTrl) = secs-t1;
                    else % response before shift: incorrect trial
                        log.response(iTrl) = 2;
                        log.responseTime(iTrl) = 0;
                    end
                elseif jFrame==nFramesTotal(iTrl) && ~flag.prevResp
                    log.response(iTrl) = 0;
                    log.responseTime(iTrl)=0;
                end
            end % stop looking for user input
        end % end frame presentation
        log.setDuration(iTrl) = shiftLatency(iTrl); % set duration till shift
        log.realDuration(iTrl) = t1-t0; % real (measured) duration till shift
        log.completeDuration(iTrl) = shiftLatency(iTrl) + postShiftTime;
        Priority(0);
    end % end trial
    save(fileName, 'log')
    
    
    % Close all textures. This is not strictly needed, as
    % Screen('CloseAll') would do it anyway. However, it avoids warnings by
    % Psychtoolbox about unclosed textures. The warnings trigger if more
    % than 10 textures are open at invocation of Screen('CloseAll') and we
    % have 12 textues here:
    Screen('Close');
    
    % Close window:
    Screen('CloseAll');
    
catch
    %this "catch" section executes in case of an error in the "try" section
    %above.  Importantly, it closes the onscreen window if its open.
    Screen('CloseAll');
    Priority(0);
    psychrethrow(psychlasterror);
end %try..catch..
