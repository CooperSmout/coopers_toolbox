function rmvar(filename,varargin)

% RMVAR FILENAME VAR1 VAR2... removes the variables VAR1, VAR2... from the Mat-File FILENAME. If FILENAME has no extension RMVAR looks for FILENAME.mat    
% RMVAR('Filename','VAR1','VAR2' ...) does the same. 

narginchk(2,inf); 
m = matfile(filename,'writable',true);
vars = whos('-file',filename); 
removeThese = {}; 
for ii=1:numel(varargin); 
    if ~any(strcmp({vars(:).name},varargin{ii})); 
        warning([ varargin{ii} ' isn''t saved in ' filename ]); 
    else 
        removeThese = [ removeThese , varargin{ii} ]; 
        m.(varargin{ii}) = []; % clear variable in matfile
    end; 
end; 
newfile = rmfield(load(filename),removeThese);
save(filename,'-struct','newfile'); 