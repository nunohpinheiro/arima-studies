%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%   fig.m           2001-06-27  %       updated  2005-02-17
%     (c)           M. Balda    %                2006-08-15
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function fig.m is prepared for fast positioning of a figure windows on 
% screen by means of the window position code. The general positioning is 
% also possible. Position codes belong to parts of the screen as follows:
%       -------------   -------------   -------------   -------------
%       |  1  |  2  |   |           |   |     |     |   |     8     |
%       |-----+-----|   |     5     |   |  6  |  7  |   |-----------|
%       |  3  |  4  |   |           |   |     |     |   |     9     |
%       -------------   -------------   -------------   -------------
% The other way to define general position of a figure window is via the
% position vector of normalized measures [left bottom width height].
% It obtainss the position code equal 10.
%
%  Function calling:
%  ~~~~~~~~~~~~~~~~~
%   fih = fig(pos);     %  Creates a figure window and its handle
%       pos     integer 1 - 9; = position code of the figure, or
%               vector of four real elements <= 1 = position measures
%       fih     handle of the figure
%
%  Examples:
%  ~~~~~~~~~
%   hf = fig([.25, .25, .5, .5]); % figure in the middle of the screen
%   fig(2);  %  figure window in the right upper quarter of the screen
%   fig(10); %  last user's figure position (here in the screen center)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fih = fig(pos)
%~~~~~~~~~~~~~~~~~~~~~~
persistent fig_

if isempty(fig_)
    fig_= zeros(10,4);
    y = .48;                        %   screen dependent constant
    fig_(1,:) = [0,  y, .5,   y];   %   upper left
    fig_(2,:) = [.5, y, .5,   y];   %   upper right
    fig_(3,:) = [0,  0, .5,   y];   %   lower left
    fig_(4,:) = [.5, 0, .5,   y];   %   lower right
    fig_(5,:) = [0,  0,  1, 2*y];   %   full screen
    fig_(6,:) = [0,  0, .5, 2*y];   %   left  half screen
    fig_(7,:) = [.5, 0, .5, 2*y];   %   right half screen
    fig_(8,:) = [0,  y,  1,   y];   %   upper half screen
    fig_(9,:) = [0,  0,  1,   y];   %   lower half screen
end
if nargin<1, return, end            %   good for startup.m

if length(pos)==1                   %   Integer position code
    pos = fig_(pos,:);
else                                %   fig window in normalized coords
    fig_(10,:) = pos;               %   user's figure window position
end

units = get(0,'Units')
fih = figure('Units','normalized','Position',pos,'toolbar','none');
set(gcf,'Units',units)