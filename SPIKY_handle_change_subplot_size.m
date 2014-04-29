% This function which is called by ‘SPIKY_handle_subplot.m’ allows changing the position (x and y) and the size (width and height)
% of the matrix subplots in the figure 

function [] = SPIKY_handle_change_subplot_size(varargin)

stri={'xpos';'ypos';'width';'height'};
pos = find(strcmp(varargin{4},stri));

if nargin<=5
    posi = get(gco,'Position');
    if strncmp(varargin{5},'+',1) || strncmp(varargin{5},'-',1)
        posi(pos) = posi(pos) + str2num(varargin{5});
    else
        posi(pos) = str2num(varargin{5});
    end
    set(gco,'Position',posi)
else
    if strncmp(varargin{5},'+',1) || strncmp(varargin{5},'-',1)
        for lhc1=1:size(varargin{3},1)
            for lhc2=1:size(varargin{3},2)
                for lhc3=1:size(varargin{3},3)
                    if varargin{3}(lhc1,lhc2,lhc3)>0
                        posi = get(varargin{3}(lhc1,lhc2,lhc3),'Position');
                        posi(pos) = posi(pos) + str2num(varargin{5});
                        set(varargin{3}(lhc1,lhc2,lhc3),'Position',posi)
                    end
                end
            end
        end
    else
        for lhc1=1:size(varargin{3},1)
            for lhc2=1:size(varargin{3},2)
                for lhc3=1:size(varargin{3},3)
                    if varargin{3}(lhc1,lhc2,lhc3)>0
                        posi=get(varargin{3}(lhc1,lhc2,lhc3),'Position');
                        posi(pos) = str2num(varargin{5});
                        set(varargin{3}(lhc1,lhc2,lhc3),'Position',posi)
                    end
                end
            end
        end
    end
end
end

