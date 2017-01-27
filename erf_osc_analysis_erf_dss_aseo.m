function erf_osc_analysis_erf_dss_aseo(subj, isPilot)
if nargin<1
    subj = 4;
end
if isempty(subj)
    subj = 4;
end
if nargin<2
    isPilot = true;
end
if isempty(isPilot);
    isPilot = true;
end


%% initiate diary
workSpace = whos;
diary('tmpDiary') % save command window output
fname = mfilename('fullpath')
datetime

fid = fopen(fullfile([fname '.m']));
tline = fgets(fid); % returns first line of fid
while ischar(tline) % at the end of the script tline=-1
    disp(tline) % display tline
    tline = fgets(fid); % returns the next line of fid
end
fclose(fid);

for i = 1:numel(workSpace) % list all workspace variables
    workSpace(i).name % list the variable name
    printstruct(eval(workSpace(i).name)) % show its value(s)
end

%% load data
erf_osc_datainfo;
if isPilot
    data = load(sprintf('/project/3011085.02/processed/pilot-%03d/ses-meg01/cleandata.mat', subj), 'dataClean');
    load(pilotsubjects(subj).logfile);% load log file
else
    data = load(sprintf('/project/3011085.02/processed/subj-%03d/ses-meg01/cleandata.mat', subj), 'dataClean');
    load(subjects(subj).logfile);% load log file
end
data = data.dataClean;
fs = data.fsample;

% select only shift trials, with valid response
idxM        = find(data.trialinfo(:,5)>0 & data.trialinfo(:,6)>0);
cfg         = [];
cfg.trials  = idxM;
cfg.channel = 'MEG';
data        = ft_selectdata(cfg, data);

% redifine zero point to shift onset
cfg         = [];
cfg.offset  = -(data.trialinfo(:,5)-data.trialinfo(:,4));
dataShift   = ft_redefinetrial(cfg, data);

cfg         = [];
cfg.latency = [-0.05 0.75];
dataShift   = ft_selectdata(cfg, dataShift);

%% DSS (denoising source seperation)
% use DSS component analysis to select a handful of components that should
% represent the ERF (data is timelocked to stimilus reversal)

s.X = 1; % 1 source signal estimate?
nComponent = 2;

% run a dss decomposition
params          = [];
params.time     = dataShift.time;
params.demean   = 'prezero';
params.pre      = 0.05*fs; % 50 ms pre stim reversal as baseline
params.pst      = 0.75*fs; % 750 ms after stim reversal

[~,~,avgorig] = denoise_avg2(params,dataShift.trial,s);

cfg                     = [];
cfg.method              = 'dss';
cfg.dss.denf.function   = 'denoise_avg2';
cfg.dss.denf.params     = params;
cfg.dss.wdim            = 100;
cfg.numcomponent        = nComponent;
cfg.cellmode            = 'yes';

% cfg for plotting
cfgp            = [];
cfgp.layout     = 'CTF275_helmet.mat';
cfgp.component  = 1:nComponent;
comp            = ft_componentanalysis(cfg, dataShift);
figure;ft_topoplotIC(cfgp, comp);drawnow;

% get the 'average'
[~,~,avgcomp] = denoise_avg2(params,comp.trial,s);

npos = sum(avgcomp>0,2);
if npos<100,
    fprintf('converting polarity\n');
    avgcomp     = -avgcomp;
    comp.topo   = -comp.topo;
    comp.trial  = -comp.trial;
end
N = size(avgcomp,2);

%% Estimate number of peaks and estimate single-trial variation of those peaks
% run the peakfit algorithm to estimate how many peaks have to be used to
% approach the real (ERF) signal. Then use the 
% estimate how many peaks (max 5) have to be used to approach real signal

% get the parameters for the peaks
% taper       = tukeywin(N,0.2)';

% peak fit algorithm
% center      = params.pre+round(params.pst/2)+1;
% window      = N-1-params.pre;
% maxnumpeaks = 3;

% f1 = cell(1,maxnumpeaks);
% f2 = cell(1,maxnumpeaks);
% for each component, fit 1-5 peaks and
for k = 1%1:nComponent;
    % peakfit algorithm
    %{
    for numpeaks = 1:maxnumpeaks
        [f1{numpeaks},f2{numpeaks},tmp1,tmp2,tmp3,xi,m{numpeaks}]=peakfit_jm(avgcomp(k,:).*taper,center,window,numpeaks,1,[],15,0,0,[],1,1);
    end
    qc          = cat(1,f2{:});
    [~,bestfit] = max(qc(:,2));
    
    initcomp = m{bestfit};
    width    = f1{bestfit}(:,4);
    area     = f1{bestfit}(:,end);
    sel      = width>20 & abs(area)>0.0001.*median(abs(area));
    
    initcomp    = initcomp(sel,:);
    width       = width(sel);
    area        = area(sel);
    jitter      = (width./2)*[-1 1];
    jitter(:,1) = max(jitter(:,1),-45);
    jitter(:,2) = min(jitter(:,2), 45);
    
    [srt, ix] = sort(sum(initcomp,2),'descend');
    initcomp  = initcomp(ix,:);
    jitter    = jitter(ix,:);
    
    
    figure;plot(initcomp','b');hold on;plot(avgcomp(k,:),'r');drawnow;
    %} 
    figure; plot(comp.time{1}, avgcomp(k,:),'r')
%     input('Please enter the peak latencies in erf_osc_datainfo. press any key to continue');
    erf_osc_datainfo;
    
    
    cfgb         = [];
    cfgb.channel = comp.label(k);
    comp_sel     = ft_selectdata(cfgb, comp);
    
    % estimate the single-trial variables for the peaks estimated with the
    % peak fitting algorithm
%     [r1(k),r2(k)] = doASEO(comp_sel,'initcomp',repmat(taper',[1 size(initcomp,1)]).*initcomp','jitter',jitter, 'numiteration', 1);
    
    % use the ASEO algorithm with visually selected latencies for the (2)
    % peaks.
    waveformInitSet = pilotsubjects(subj).aseo;
    waveformInitSet = waveformInitSet';
% waveformInitSet = [119 161; 162 226; 230 459]';
    waveformInitSet = waveformInitSet(:);
    jitter = [-50 50; -50 50; -50 50];
    ASEOiteration = 1;
    [q1(k), q2(k)] = doASEO(comp_sel, 'waveformInitSet', waveformInitSet, 'jitter', jitter, 'numiteration', ASEOiteration);
end

%% save
if isPilot
    filename = sprintf('/project/3011085.02/results/erf/pilot-%03d/dss_ASEO', subj);
else
    filename = sprintf('/project/3011085.02/results/erf/subj-%03d/dss_ASEO', subj);
end
save(fullfile([filename '.mat']),'comp', 'q1', 'q2', 'avgorig', 'avgcomp')
diary off
movefile('tmpDiary', fullfile([filename, '.txt']));


end
