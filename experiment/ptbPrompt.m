
% screen prompt
% set up for text size ~20
% input is any number of text strings

function ptbPrompt(key, varargin)

nLines = numel(varargin) + 2;

mon = Screen('Windows');
mon = mon(1);
[~, y] = RectCenter(opt.ptb.mainRect)

for line = 1:numel(varargin)
    
    DrawFormattedText(mon, varargin{line}, 'center', y + 30*((nLines+1)/2 - line)); % cgtext('Failed to detect heartbeat, restarting trial. ',0,0)
    % cgtext(varargin{line}, 0, 30*((nLines+1)/2 - line))
    
end

DrawFormattedText(mon, '(press SPACEBAR to continue)', 'center', y + 30*((nLines+1)/2 - nLines)); 
Screen('Flip', mon); 
KbWait4keys(key)