function erfosc_execute_pipeline(pipelinename, subj, varargin)

% ERFOSC_EXECUTE_PIPELINE serves the purpose to make a script executable by qsub.
% supply it with the name of the script that has to be run, and the
% subject identifier. Other variables need to be passed as key-value pairs
%
% example use in combination with qsubfeval:
%
% qsubfeval('erfosc_superscript', 1, {'someothervariable' value}, 'memreq', memreq, 'timreq', timreq);
% 
% example use (standalone; not so useful):
%
% erfosc_execute_pipeline('erfosc_superscript', 1);



if numel(varargin)>0
  for k = 1:numel(varargin)
    eval([varargin{k}{1},'=varargin{k}{2}']);
  end
end
eval(pipelinename);
