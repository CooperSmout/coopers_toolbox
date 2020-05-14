function [fmt_str, C] = struct2string(opt)

%// Extract field data
fields = repmat(fieldnames(opt), numel(opt), 1);
values = struct2cell(opt);

% remove structure fields
idDel = cellfun(@isstruct,values); 
fields = fields(~idDel);
values = values(~idDel);

%// Convert all numerical values to strings
idx = cellfun(@isnumeric, values); 
values(idx) = cellfun(@num2str, values(idx), 'UniformOutput', 0);

%// Combine field names and values in the same array
C = {fields{:}; values{:}};
fmt_str = repmat('%s,', 1, size(C, 2));
end