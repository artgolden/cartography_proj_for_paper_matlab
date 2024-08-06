function path = custom_fullfile(varargin)
    custom_filesep = '/';
    path = strrep(fullfile(varargin{:}), filesep, custom_filesep);
end
