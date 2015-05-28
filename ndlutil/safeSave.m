function safeSave(fname, varargin)

% SAFESAVE Safe save
% FORMAT
% DESC Saves given variables in such a way that existing files are
%   only overwritten with a complete save, thus making sure that a
%   complete copy exists even if the call is interrupted
% ARG fname : name of the file to save to (should include extension)
% ARG ... : the variables to save (either as names (strings) or as
%   variables
%
% SEEALSO save
%
% COPYRIGHT : Antti Honkela, 2010

% NDLUTIL

  saveargs = cell(size(varargin));
  for k=1:length(varargin),
    if isempty(inputname(k+1)),
      saveargs{k} = varargin{k};
    else
      saveargs{k} = inputname(k+1);
    end
  end
  argstring = sprintf(', ''%s''', saveargs{:});
  
  % First time saving to this file
  if ~exist(fname, 'file'),
    callstr = sprintf('save(''%s''%s)', fname, argstring);
    evalin('caller', callstr);
    return;
  end
  
  % If file exists, first write to a different name
  [pathstr, name, ext] = fileparts(fname);
  tmpfname = fullfile(pathstr,[name '_tmp' ext]);
  callstr = sprintf('save(''%s''%s)', tmpfname, argstring);
  evalin('caller', callstr);
  
  % Then move to overwrite the original
  movefile(tmpfname, fname);
end
