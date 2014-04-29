% This function can be used to obtain from piecewise constant measure profiles
% (that can be extracted via the red mouse button while hovering over a profile)
% the average value as well as the x- and y-vectors for plotting.

function [ave,plot_x_values,plot_y_values]=SPIKY_f_pico(isi,pico,tmin)

cum_isi=cumsum([0 isi]);
ave=sum(pico.*repmat(isi,size(pico,1),1),2)/cum_isi(end);
plot_x_values=tmin+sort([cum_isi(1:end) cum_isi(2:end-1)]);
plot_y_values=reshape([pico; pico],1,length(pico)*2);

