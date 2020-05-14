

function removeMatfileVar(m,varargin)
%

% create temporary matfile
tmpname = [m.Properties.Source(1:end-4) '_TEMP.mat'];
newm = matfile(tmpname,'writable',true);

% copy all fields to new matfile
fields = who(m);
for F = 1:length(fields)
    if ~any(strcmp(fields{F},varargin))
        newm.(fields{F}) = m.(fields{F});
    end
end

% rename original matfile
movefile(tmpname,m.Properties.Source)
    
end

    