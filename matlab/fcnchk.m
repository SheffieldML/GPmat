function [f,msg] = fcnchk(fun,varargin)
%FCNCHK Check FUNFUN function argument.
%   FCNCHK(FUN,...) returns an inline object based on FUN if FUN
%   is a string containing parentheses, variables, and math
%   operators.  FCNCHK simply returns FUN if FUN is a function handle, 
%   or a matlab object with an feval method (such as an inline object). 
%   If FUN is a string name of a function (e.g. 'sin'), FCNCHK returns a
%   function handle to that function.
%
%   FCNCHK is a helper function for FMINBND, FMINSEARCH, FZERO, etc. so they
%   can compute with string expressions in addition to m-file functions.
%
%   FCNCHK(FUN,...,'vectorized') processes the string (e.g., replacing
%   '*' with '.*') to produce a vectorized function.
%
%   When FUN contains an expression then FCNCHK(FUN,...) is the same as
%   INLINE(FUN,...) except that the optional trailing argument 'vectorized'
%   can be used to produce a vectorized function.
%
%   [F,MSG] = FCNCHK(...) returns an empty structure in MSG if successful
%   or an error message structure if not.
%
%   See also ERROR, INLINE, @.

%   Copyright 1984-2004 The MathWorks, Inc.
%   $Revision: 1.29.4.8 $  $Date: 2004/07/28 04:25:42 $

message = '';
msgident = '';

nin = nargin;
if (nin>1) && strcmp(varargin{end},'vectorized')
    vectorizing = 1;
    nin = nin-1;
else
    vectorizing = 0;
end

if ischar(fun)
    fun = strtrim_local_function(fun);
    % Check for non-alphanumeric characters that must be part of an
    % expression.
    if isempty(fun),
        f = inline('[]');
    elseif ~vectorizing && isidentifier_local_function(fun)
        f = str2func(fun); % Must be a function name only
        % Note that we avoid collision of f = str2func(fun) with any local
        % function named fun, by uglifying the local function's name
        if isequal('x',fun)
            warning('MATLAB:fcnchk:AmbiguousX', ...
                ['Ambiguous expression or function input.\n The string ''x'' will be ',...
                'interpreted as the name of a ',...
                'function called ''x'' \n (e.g., x.m) and not as the mathematical expression ''x'' (i.e., f(x)=x). \n ',...
                'Use the anonymous function:  @(x)x  ',...
                'if you meant the mathematical expression ''x''.']);
        end
    else
        if vectorizing
            f = inline(vectorize(fun),varargin{1:nin-1});
            var = argnames(f);
            f = inline([formula(f) '.*ones(size(' var{1} '))'],var{1:end});
        else
            f = inline(fun,varargin{1:nin-1});
        end 
    end
elseif isa(fun,'function_handle') 
    f = fun; 
    % is it a matlab object with a feval method?
elseif isobject(fun)
    % delay the methods call unless we know it is an object to avoid runtime error for compiler
    meths = methods(class(fun));
    if any(strmatch('feval',meths,'exact'))
       if vectorizing && any(strmatch('vectorize',meths,'exact'))
          f = vectorize(fun);
       else
          f = fun;
       end
    else % no feval method
        f = '';
        message = 'If FUN is a MATLAB object, it must have an feval method.';
        msgident = 'MATLAB:fcnchk:objectMissingFevalMethod';
    end
else
    f = '';
    message = ['FUN must be a function, a valid string expression, ', ...
            sprintf('\n'),'or an inline function object.'];
    msgident = 'MATLAB:fcnchk:invalidFunctionSpecifier';
end

% If no errors and nothing to report then we are done.
if nargout < 2 && isempty(message)
    return
end

% compute MSG
if isempty(message)
    msg.message = '';
    msg.identifier = '';
    msg = msg(zeros(0,1)); % make sure msg is the right dimension
else
    msg.message = message;
    msg.identifier = msgident;
end

if nargout < 2, error(msg); end


%------------------------------------------
function s1 = strtrim_local_function(s)
%STRTRIM_LOCAL_FUNCTION Trim spaces from string.
% Note that we avoid collision with line 45: f = str2func('strtrim')
% by uglifying the local function's name

if isempty(s)
    s1 = s;
else
    % remove leading and trailing blanks (including nulls)
    c = find(s ~= ' ' & s ~= 0);
    s1 = s(min(c):max(c));
end

%-------------------------------------------
function tf = isidentifier_local_function(str)
% Note that we avoid collision with line 45: f = str2func('isidentifier')
% by uglifying the local function's name


tf = false;

if ~isempty(str)
    first = str(1);
    if (isletter(first))
        letters = isletter(str);
        numerals = (48 <= str) & (str <= 57);
        underscore = (95 == str);
        tf = all(letters | numerals | underscore);
    end
end

