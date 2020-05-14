function saveasCS(h, name, varargin)

% SAVEASCS Saves matlab figure after making the following changes:
%   - adds title to figure before saving
%   - title is derived from the filename (excluding the path, if one exists)
%
%   see SAVEAS for full documentation
%
%   created by Cooper Smout (ORCID: 0000-0003-1144-3272)


% default format
if nargin<3
    formats = {'fig','png'};
else
    formats = varargin;
end

% change default paper size to match onscreen size
if any(strcmp(formats,'png')) && isequal(get(h,'paperposition'),[0.6350    6.3500   20.3200   15.2400])
    screenSize = get(h,'position');
    screenX = screenSize(3);
    screenY = screenSize(4);
    set(h,'paperpositionmode','auto')
    set(h,'papersize',[screenX/10 screenY/10])
end

% title
figure(h)
[~,tit] = fileparts(name);
suptitle(strrep(tit,'_',' '))

% save
for f = 1:length(formats)
    saveas(h,name,formats{f});
end


