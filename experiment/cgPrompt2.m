
% screen prompt
% set up for text size ~20
% input is any number of text strings

function cgPrompt2(varargin)



% concatenate strings
strings{1} = '';
for arg = 1:nargin
    if isempty(varargin{arg})
        strings{end+1} = ''; % make next line blank
        strings{end+1} = ''; % start new line after blank
        
    else
        strings{end} = [strings{end} ' ' varargin{arg}];
    end
end


% cycle through string groups
blockCount = 0;
for S = 1:length(strings)

    % break into lines
    limit = 90;
    chars = length(strings{S});
    nLines = ceil(chars/limit);
    charsPerLine = chars/nLines;
    lineCount = 1;
    lines{blockCount+lineCount} = '';
    
    % split string group into lines
    for C = 1:chars
        if strcmp(strings{S}(C),' ') && (C > lineCount*limit)
            lineCount = lineCount + 1;
            lines{blockCount+lineCount} = '';
        else
            lines{blockCount+lineCount}(end+1) = strings{S}(C);
        end
    end
    
    blockCount = blockCount + lineCount;
    
end
    
    
% display text on screen
for L = 1:length(lines)
    cgtext(lines{L}, 0, 30*((length(lines)+1)/2 - L))
end
cgtext('(press SPACEBAR to continue)', 0, -30*length(lines)/2 - 50)%((length(lines)+1)/2 - nLines))
cgflip(0,0,0)
waitkeydown( inf,71 );