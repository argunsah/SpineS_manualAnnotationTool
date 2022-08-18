function [X,Y] = ginputzoom;

button = 0;
zlvl   = 1;
xl     = get(gca,'xlim');
xlen   = 1600;
yl     = get(gca,'ylim');
ylen   = 1600;

X = [];
Y = [];

while button~=3
    % Get the mouse position on the axes (needed for binary image editing) and button number
    [x,y,button] = myginput(1,'arrow');
    if x<xlen & y<ylen
        % Determine if it is a zoom-in or zoom-out
        if button==122
            zlvl = zlvl*2;
            zoom(2);
        elseif button==120
            zlvl = zlvl/2;
            if zlvl<1, zlvl=1; end % No zoom level smaller than 1
            zoom(0.5);
        elseif button==1
               X = [X x];
               Y = [Y y];
               hold on, plot(x,y,'.y');
        end
        % Change the axis limits to where the mouse click has occurred
        % and make sure that the display window is within the image dimensions
        xlimit = [x-xlen/zlvl/2+0.5 x+xlen/zlvl/2+0.5];
        if xlimit(1)<0.5, xlimit=[0.5 xlen/zlvl+0.5]; end
        if xlimit(2)>0.5+xlen, xlimit=[xlen-xlen/zlvl+0.5 xlen+0.5]; end
        xlim(xlimit);

        ylimit = [y-ylen/zlvl/2+0.5 y+ylen/zlvl/2+0.5];
        if ylimit(1)<=0.5, ylimit=[0.5 ylen/zlvl+0.5]; end
        if ylimit(2)>=0.5+ylen, ylimit=[ylen-ylen/zlvl+0.5 ylen+0.5]; end
        ylim(ylimit);
    end
end

hA = gca;
resetplotview(hA,'InitializeCurrentView');
set(hA,'xlim',[1 1600]);
set(hA,'ylim',[1 1600]);