function concentric_grating(isLive, drawMask)
% This function is copied from DriftDemo and modified to display an
% animated concentric grating moving inward.

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
%
% _________________________________________________________________________

% HISTORY
% 7/6/2016: change grating to concentric grating
% 8/6/2016: add circular and gaussian mask, introduce phase-shift

%% Function settings
if nargin<1
    isLive = 0;
end

if isempty(isLive)
    isLive = 0;
end

if nargin<2
    drawMask = 0;
end

if isempty(drawMask)
    drawMask = 0;
end
try
    
    %% Standard PTB settings
    % This script calls Psychtoolbox commands available only in OpenGL-based
    % versions of the Psychtoolbox. (So far, the OS X Psychtoolbox is the
    % only OpenGL-base Psychtoolbox.)  The Psychtoolbox command AssertPsychOpenGL will issue
    % an error message if someone tries to execute this script on a computer without
    % an OpenGL Psychtoolbox
    AssertOpenGL;
    
    % Get the list of screens and choose the one with the highest screen number.
    % Screen 0 is, by definition, the display with the menu bar. Often when
    % two monitors are connected the one without the menu bar is used as
    % the stimulus display.  Chosing the display with the highest dislay number is
    % a best guess about where you want the stimulus displayed.
    screens=Screen('Screens');
    screenNumber=max(screens);
    
    % Find the color values which correspond to white and black: Usually
    % black is always 0 and white 255, but this rule is not true if one of
    % the high precision framebuffer modes is enabled via the
    % PsychImaging() commmand, so we query the true values via the
    % functions WhiteIndex and BlackIndex:
    white=WhiteIndex(screenNumber);
    black=BlackIndex(screenNumber);
    
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
    [window, windowRect]=Screen('OpenWindow',screenNumber, black);
    
    % Get the centre coordinate of the window
    [xCenter, yCenter] = RectCenter(windowRect);
    
    % Query the frame duration
    ifi = Screen('GetFlipInterval', window);
    frameRate=Screen('FrameRate',screenNumber); %1/ifi
    
    % Set up alpha-blending for smooth (anti-aliased) lines
    Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
    
    HideCursor;
    
    %% Experimental settings
    if isLive
        screen.width = 49; %MEG screenwidth in cm
        screen.viewingDistance = 76; % in cm
        screen.resolutionX = 1024; %MEG screen in pixels
    else
        screen.width = 52;%desktop screenwidth in cm
        screen.viewingDistance = 59; % in cm
        screen.resolutionX = 1920; %desktop screen in pixels
    end
    
    screen.totVisDeg = 2*atan(screen.width/(2*screen.viewingDistance))*(180/pi);
    screen.pixPerDeg = screen.resolutionX/screen.totVisDeg;
    
    visualAngleGrating = 25;
    gratingSize = visualAngleGrating*screen.pixPerDeg; % calculate how big the grating should be (can not be rounded off, might result in odd integer)
    gratingRadius = round(gratingSize/2); % the grating can only exist of integers, so round
    gratingSize = 2*gratingRadius; % to prevent consistency errors, redifine gratingSize
    
    % Run the movie animation for a fixed period.
    movieDurationSecs=8;
    % Convert movieDuration in seconds to duration in frames to draw:
    movieDurationFrames=round(movieDurationSecs * frameRate);
    
    shiftPoint = 2+4*rand(1,1); % generate a random number between 2 and 6 seconds for the timing of the phase shift
    shiftFrame = round(shiftPoint*frameRate);
    driftFreq = 1; % in one second, every pixel of the grating completes one cylce (black-white-black)
    numFrames=round(driftFreq/ifi); % temporal period, in frames, of the drifting grating
        
    % settings grating
    [x,y]=meshgrid(-gratingRadius:gratingRadius,-gratingRadius:gratingRadius);
    f=1.5*2*pi; % cycles/pixel
    
    %% Generate stimulus
    % generate a mask to cover the edges of the grating, making it a circle
    R = gratingRadius ; % radius
    M = zeros(2*gratingRadius+1, 2*gratingRadius+1, 1); %zero matrix same size as the grid
    M(:,:,1) = sqrt(x.^2 + y.^2) < R ; %disk formula: all masked parts of the grid must be zero.
    
    
    % Generate grating texture
    % Compute each frame of the movie and convert those frames stored in
    % MATLAB matrices, into Psychtoolbox OpenGL textures using 'MakeTexture';
    
    %________________________________________________
    % very inefficient, costs a lot of memory.. alternative?
    
    
    for i=1:numFrames
        phase=(i/numFrames)*2*pi; %change the phase of the grating according to the framenumber
        m=sin(sqrt(x.^2+y.^2) / f + phase);
        tex(i)=Screen('MakeTexture', window, M.*(gray+inc*m)); % multiply mask with grating
    end
    
    %________________________________________________

    mask=ones(2*gratingRadius+1, 2*gratingRadius+1, 2) * black;
    mask(:, :, 2)=white * (1 - exp(-((x/300).^2)-((y/300).^2)));
    masktex=Screen('MakeTexture', window, mask);
    %% Present stimuli 
    % assign the right texture indices to frames
    movieFrameIndices=mod(0:(movieDurationFrames-1), numFrames) + 1;
    
    % Use realtime priority for better timing precision:
    priorityLevel=MaxPriority(window);
    Priority(priorityLevel);
    
    
    
    
    % Animation loop:
    tic
    for i=1:movieDurationFrames
        % at the moment the phase-shift should occur, jump to the texture
        % with a 180 degree phase shift, and add this shift for the
        % remainder of the video. Otherwise, don't shift.
        if i>=shiftFrame && movieFrameIndices(i)>numFrames/2
            shift = -numFrames/2; % if more than half of the period already passed, subtract the shift
        elseif i>=shiftFrame && movieFrameIndices(i)<=numFrames/2
            shift = numFrames/2; % otherwise add the shift (because the period loops)
        else
            shift=0;
        end

        Screen('DrawTexture', window, tex(movieFrameIndices(i) + shift)); 
        if drawMask
        Screen('DrawTexture', window, masktex);
        end
        % Show it at next display vertical retrace. Please check DriftDemo2
        % and later, as well as DriftWaitDemo for much better approaches to
        % guarantee a robust and constant animation display timing! This is
        % very basic and not best practice!
        Screen('Flip', window);
    end
    toc
    Priority(0);
    
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
