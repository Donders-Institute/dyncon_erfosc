function [data, avgorig, avgnew, avgcomp, comp, r1, r2] = erf_osc_analysis_erf_dss_aseo(subj, isPilot, varargin)
if nargin<1
    subj = 1;
end
if isempty(subj)
    subj = 1;
end
if nargin<2
    isPilot = true;
end
if isempty(isPilot);
    isPilot = true;
end
if nargin<3
    varargin=[];
end
if isempty(varargin)
    varargin=[];
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
    data = load(sprintf('/home/electromag/matves/Data/ERF_oscillation/clean_data/pilot/%02d/cleandata.mat', subj), 'dataClean');
    load(pilotsubjects(subj).logfile);% load log file
else
    data = load(sprintf('/home/electromag/matves/Data/ERF_oscillation/clean_data/experiment/%02d/cleandata.mat', subj), 'dataClean');
    load(subjects(subj).logfile);% load log file
end
data = data.dataClean;

% select only shift trials, with valid response
idxM = find(data.trialinfo(:,5)>0 & data.trialinfo(:,6)>0);
cfg=[];
cfg.trials = idxM;
cfg.channel = 'MEG';
dataM = ft_selectdata(cfg, data);

% redifine zero point to shift onset
cfg=[];
cfg.offset = -(dataM.trialinfo(:,5)-dataM.trialinfo(:,4));
dataShiftM = ft_redefinetrial(cfg, dataM);

cfg=[];
cfg.latency = [-0.05 0.75];
dataShiftM = ft_selectdata(cfg, dataShiftM);

%%
Ncomp = ft_getopt(varargin, 'Ncomp', 5);

addpath('/home/electromag/matves/MATLAB/dss_aseo/peakfit');
addpath('/home/electromag/matves/MATLAB/dss_aseo/ASEO');
addpath('/home/electromag/matves/MATLAB/dss_aseo/dss2_1-0');

s.X = 1; % 1 source signal estimate?

% run a dss decomposition
params      = [];
params.time = dataShiftM.time;
params.demean = 'prezero';
params.pre = 60; % 50 ms pre shift as baseline
params.pst = 900; % 750 ms after shift

% [~,~,avgorig] = denoise_avg2(params,dataOnsetMBl.trial,s);

cfg          = [];
cfg.method   = 'dss';
cfg.dss.denf.function = 'denoise_avg2';
cfg.dss.denf.params = params;
cfg.dss.wdim = 100;
cfg.numcomponent = 1;
cfg.cellmode = 'yes';

% cfg for plotting
cfgp = [];
cfgp.layout = 'CTF275_helmet.mat';
cfgp.component = 1;

for k = 1:Ncomp;
  comp(k) = ft_componentanalysis(cfg, dataShiftM);
  figure;ft_topoplotIC(cfgp, comp(k));drawnow;
  
  % get the 'average'
  [~,~,avgcomp(k,:)] = denoise_avg2(params,comp(k).trial,s);
  
  npos = sum(avgcomp(k,:)>0);
  if npos<100,
    fprintf('converting polarity\n');
    avgcomp(k,:) = -avgcomp(k,:);
    comp(k).topo = -comp(k).topo;
    comp(k).trial = -comp(k).trial;
  end
  N = size(avgcomp,2);
  
  % run the peakfit algorithm
  
  % get the parameters for the peaks
  taper  = tukeywin(N,0.2)';
  center = params.pre+round(params.pst/2)+1;
  window = N-1-params.pre;
  maxnumpeaks = 5;
  
  f1 = cell(1,maxnumpeaks);
  f2 = cell(1,maxnumpeaks);
  for numpeaks = 1:maxnumpeaks
    [f1{numpeaks},f2{numpeaks},tmp1,tmp2,tmp3,xi,m{numpeaks}]=peakfit_jm(avgcomp(k,:).*taper,center,window,numpeaks,1,[],15,0,0,[],1,1);
  end
  qc          = cat(1,f2{:});
  [~,bestfit] = max(qc(:,2));
  
  initcomp = m{bestfit};
  width = f1{bestfit}(:,4);
  area  = f1{bestfit}(:,end);
  sel   = width>20 & abs(area)>0.0001.*median(abs(area));
  
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
  
  [r1(k),r2(k)] = doASEO(comp(k),'initcomp',repmat(taper',[1 size(initcomp,1)]).*initcomp','jitter',jitter);
  
  % % create a projection matrix to remove the estimated component
  % P = eye(numel(data.label))-comp.topo*comp.unmixing;
  % data.trial = P*data.trial + comp.topo*r2.trial; % remove the component and add the residuals
  
%   data.trial   = data.trial - comp(k).topo*r1(k).trial; % remove the modelled data
%   [~,~,avgnew(:,:,k)] = denoise_avg2(params,data.trial,s);
end
