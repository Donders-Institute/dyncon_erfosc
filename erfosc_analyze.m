clear
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global ft_default
ft_default = [];
ft_default.checksize = inf;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% the analysis follows the strategy outlined on https://doi.org/10.1016/j.neuroimage.2018.11.029
reproducescript_dir = '/project_ext/3010029/reproducescript/';
subject_list = [1 2 3 4 5 6 7 8 9 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33];

for subj=1
  ft_default.reproducescript = [reproducescript_dir, 'reproduce_', sprintf('%03d', subj)];
  erfosc_execute_pipeline('erfosc_script', subj, {'getdata', true}, {'dofreq', 'true'}, ...
    {'dodics', true}, {'dohighfreq', true}, {'dofreq_short', true}, {'dolcmv_parc', true},...
    {'dolcmv_norm', true}, {'docorrpow_lcmv', true}, {'dopartialcorr', false}, {'dosave', false});
  ft_default.reproducescript = []; % disable
end

% for subj=1
%   ft_default.reproducescript = [reproducescript_dir, 'reproduce_Group', {'dolcmv_corr', false});
%   ft_default.reproducescript = []; % disable
% end



function erfosc_execute_pipeline(pipelinename, subj, varargin)

if numel(varargin)>0
  for k = 1:numel(varargin)
    eval([varargin{k}{1},'=varargin{k}{2}']);
  end
end
eval(pipelinename);
end