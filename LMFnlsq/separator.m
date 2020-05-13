%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   SEPARATOR.M     2005-02-24  %   Separating line with text in the middle
%       (c)         M. Balda    %           modified on 2006-03-30
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The function creates a separating text line on the screen and in records 
% (diaries). User defined text can be placed in the middle of the line.
% The length of the separator is text-length independent.
%
% Forms of calls:
% ~~~~~~~~~~~~~~~
%   separator                   %   default values of params will be used
%   separator(txt)
%   separator(txt,chr)  
%   separator(txt,chr,len)
%       txt     inserted text; '' default
%       chr     separating character:
%               '-' default in the first call,
%               old c
%       len     half-length of the line,    default 25 characters
% Examples:
% ~~~~~~~~
%   separator('EXAMPLE')
%   separator('2nd EXAMPLE','=')
%   separator('3rd EXAMPLE')
%   separator('Another example','*')
%   %   Outputs:
%   %   -------------------- EXAMPLE ---------------------
%   %   ================== 2nd EXAMPLE ===================
%   %   ================== 3rd EXAMPLE ===================
%   %   **************** Another example *****************

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function separator(txt,chr,len)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent len_ chr_

if isempty(len_), len_=25; end
if isempty(chr_), chr_='-'; end

if nargin<1, txt = ''; end
if nargin>1, chr_=chr; end
if nargin>2, len_=len; end 

if length(txt)>0, txt=[' ' txt ' ']; end
n    = fix((2*len_-length(txt))/2);
Line = char(ones(1,n)*chr_);
Line = [Line txt Line];
len  = length(Line);
if len-fix(len/2)*2~=0
    Line = [Line chr_];
end

fprintf('\n%s\n\n',Line);
