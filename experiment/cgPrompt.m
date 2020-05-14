
% screen prompt
% set up for text size ~20
% input is any number of text strings

function cgPrompt(varargin)

nLines = numel(varargin) + 2;

for line = 1:numel(varargin)
    cgtext(varargin{line}, 0, 30*((nLines+1)/2 - line))
end
cgtext('(press SPACEBAR to continue)', 0, 30*((nLines+1)/2 - nLines))
cgflip(0,0,0)
waitkeydown( inf,71 );