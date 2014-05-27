% SPIKY_loop_surro --- Copyright Thomas Kreuz, Nebojsa Bozanic; March 2014
%
% This is similar to "SPIKY_loop.m" only that you do not look at different datasets;
% rather you can compare the results obtained for one dataset against the results
% obtained for spike train surrogates generated from that dataset.  
% [This program is not yet completely finished, there will be some refinements soon.]
%
% More information on the program and the spike train distances can be found under
% "http://www.fi.isc.cnr.it/users/thomas.kreuz/Source-Code/SPIKY.html" and/or in
%
% Kreuz T, Chicharro D, Houghton C, Andrzejak RG, Mormann F: Monitoring spike train synchrony. J Neurophysiol 109, 1457 (2013)
% Kreuz T: SPIKE-distance. Scholarpedia 7(12):30652 (2012).
% Kreuz T: Measures of spike train synchrony. Scholarpedia, 6(10):11934 (2011).
% Kreuz T, Chicharro D, Greschner M, Andrzejak RG: Time-resolved and time-scale adaptive measures of spike train synchrony. J Neurosci Methods 195, 92 (2011).
% Kreuz T, Chicharro D, Andrzejak RG, Haas JS, Abarbanel HDI: Measuring multiple spike train synchrony. J Neurosci Methods 183, 287 (2009).
% Kreuz T, Haas JS, Morelli A, Abarbanel HDI, Politi A: Measuring spike train synchrony. J Neurosci Methods 165, 151 (2007)
%
% For questions and comments please contact us at "thomaskreuz (at) cnr.it".
%
%
%
% Input:
% ======
%
% Cell 'spikes' with two or more spike trains (each cell array contains the spike times of one spike train)
%
% Parameter structure 'para' that describe the data (see below)
%
% tmin:            Beginning of recording
% tmax:            End of recording
% dts:             Sampling interval, precision of spike times
%                  [!!! Please take care that this value is not larger than the actual sampling size,
%                   otherwise two spikes can occur at the same time instant and this can lead to problems in the algorithm !!!]
% select_measures: Vector with measure selection (for order see below)
%
%
% Output (Structure 'SPIKY_results'):
% =============================
%
%    SPIKY_results.<Measure>.name:     Name of selected measures (helps to identify the order within all other variables)
%    SPIKY_results.<Measure>.distance: Level of dissimilarity over all spike trains and the whole interval
%                                just one value, obtained by averaging over both spike trains and time
%    SPIKY_results.<Measure>.matrix:   Pairwise distance matrices, obtained by averaging over time
%    SPIKY_results.<Measure>.x:        Time-values of overall dissimilarity profile
%    SPIKY_results.<Measure>.y:        Overall dissimilarity profile obtained by averaging over spike train pairs
%
% Note: For the ISI-distance the function 'SPIKY_f_pico' can be used to obtain the average value as well as
% x- and y-vectors for plotting (see example below):
%
% [overall_dissimilarity,plot_x_values,plot_y_values] = SPIKY_f_pico(SPIKY_results.isi,SPIKY_results.dissimilarity_profiles{1},para.tmin);
%

clear all
close all
clc

para=struct('tmin',[],'tmax',[],'dts',[],'select_measures',[]);            % Initialization of parameter structure



dataset=2;            % 1-Frequency mismatch,2-Spiking events,3-Random


m_para.all_measures_string={'ISI';'SPIKE';'SPIKE_realtime';'SPIKE_future';'PSTH';};  % order of select_measures

para.select_measures=[0 1 0 0 1];  % Select measures (0-calculate,1-do not calculate)

num_surros=2;
para.choice=4;   % 1-keep spike train numbers,2-keep interval distribution,3-keep pooled spike train,4-keep PSTH

plotting=3;           % +1:spikes,+2:dissimilarity profile,+4:dissimilarity matrix


% ################################################### Example spike trains

if dataset==1                    % Frequency mismatch (from Fig. 2a of 2013 paper)
    para.tmin=0; para.tmax=1300; para.dts=1;
    num_trains=2;
    spikes=cell(1,num_trains);
    spikes{1}=(100:100:1200);
    spikes{2}=(100:110:1200);
elseif dataset==2              % Spiking events (from Fig. 2b of 2013 paper)
    para.tmin=0; para.tmax=4000; para.dts=1;
    num_trains=50; num_spikes=40;
    noise=[1:-1/(num_spikes/4-1):0 0];
    num_noises=length(noise);
    num_events=5;
    spikes=cell(1,num_trains);
    for trc=1:num_trains
        spikes{trc}=sort(rand(1,num_spikes/2),2)*para.tmax/2;
        for nc=1:num_events
            spikes{trc}=[spikes{trc} nc*para.tmax/2/num_events+50*noise(ceil(num_noises-(nc-1)*num_noises/num_events)).*randn];
        end
        spikes{trc}=[spikes{trc} 100*(num_spikes/2-1)+200*(1:num_spikes/4+1)+50*noise.*randn(1,num_spikes/4+1)];
    end
    
elseif dataset==3
    para.tmin=0; para.tmax=4000; para.dts=1;
    num_trains=50; num_spikes=40;
    spikes=cell(1,num_trains);
    for trc=1:num_trains
        spikes{trc}=sort(rand(1,num_spikes),2)*para.tmax;
    end
end

d_para=para;
SPIKY_check_spikes
if ret==1
    return
end
para=d_para;

ori_spikes=spikes;
surro_results=cell(1,num_surros);
surro_distances=zeros(1,num_surros);
for sc=1:num_surros+1

    % ################################################### Actual call of function !!!!!

    if sc<num_surros+1
        spikes=SPIKY_f_spike_train_surrogates(ori_spikes,para);
        SPIKY_results = SPIKY_loop_f_distances(spikes,para);
        surro_results{sc}=SPIKY_results;
        surro_distances(sc)=SPIKY_results.SPIKE.distance;
        if sc==num_surros
            mean_surro_distance=mean(surro_distances);
            std_surro_distance=std(surro_distances);
        end
    else
        spikes=ori_spikes;
        SPIKY_results = SPIKY_loop_f_distances(spikes,para);
        ori_results=SPIKY_results;
    end
    % ################################################### Example plotting (just as a demonstration)

    num_plots=(mod(plotting,2)>0)+(mod(plotting,4)>1)+(mod(plotting,8)>3);
    if num_plots>0
        measures=find(para.select_measures);
        for mc=1:length(measures)
            measure=measures(mc);
            measure_var=m_para.all_measures_string{measure};
            measure_name=regexprep(measure_var,'_','-');

            figure(sc*100+mc); clf
            set(gcf,'Name',measure_name)
            set(gcf,'Units','normalized','Position',[0.0525 0.0342 0.8854 0.8867])
            subplotc=0;

            if mod(plotting,2)>0
                subplotc=subplotc+1;
                subplot(num_plots,1,subplotc)                                      % Spikes
                for trc=1:length(spikes)
                    for spc=1:length(spikes{trc})
                        line(spikes{trc}(spc)*ones(1,2),[trc-1 trc])
                    end
                end
                xlim([para.tmin para.tmax])
                ylim([0 length(spikes)])
                title ('Spike trains','FontWeight','bold','FontSize',14)
            end

            if mod(plotting,4)>1
                subplotc=subplotc+1;
                subplot(num_plots,1,subplotc)                                      % Dissimilarity profile
                if measure==1                                % piecewise constant profiles have first to be transformed
                    isi_x=SPIKY_results.(measure_var).x;
                    isi_y=SPIKY_results.(measure_var).y;
                    plot_y_values=zeros(size(isi_y,1),length(isi_x)*2);
                    for pc=1:size(isi_y,1)
                        [overall_dissimilarity,plot_x_values,plot_y_values(pc,:)] = SPIKY_f_pico(isi_x,isi_y(pc,:),para.tmin);
                    end
                    hold on
                    plot(plot_x_values,plot_y_values)
                    plot(plot_x_values,plot_y_values(1,:),'LineWidth',1.5)
                elseif ismember(measure,[2 3 4])             % piecewise linear profiles can be plotted right away
                    x=SPIKY_results.(measure_var).x;
                    y=SPIKY_results.(measure_var).y;
                    hold on
                    plot(x,y)
                    plot(x,y(1,:),'LineWidth',1.5)
                elseif measure==5                            % PSTH
                    x=SPIKY_results.(measure_var).x;
                    y=SPIKY_results.(measure_var).y;
                    hold on
                    plot(x,y)
                    plot(x,y(1,:),'LineWidth',1.5)
                end
                xlim([para.tmin para.tmax])
                if measure~=5
                    line([para.tmin para.tmax],SPIKY_results.(measure_var).distance*ones(1,2),'LineWidth',2)
                    if sc==num_surros+1
                        line([para.tmin para.tmax],(mean_surro_distance+std_surro_distance)*ones(1,2),'Color','g','LineWidth',1)
                        line([para.tmin para.tmax],(mean_surro_distance-std_surro_distance)*ones(1,2),'Color','g','LineWidth',1)
                        line([para.tmin para.tmax],mean_surro_distance*ones(1,2),'Color','r','LineWidth',2)
                    end
                end
                title ([measure_name,'   ---   Dissimilarity profile'],'FontWeight','bold','FontSize',14)
            end

            if mod(plotting,8)>3 && measure<5
                subplotc=subplotc+1;
                subplot(num_plots,1,subplotc)                                      % Dissimilarity matrix
                mat=SPIKY_results.(measure_var).matrix;
                imagesc(mat)
                axis square
                if size(mat,1)<10
                    set(gca,'XTick',1:size(mat,1),'YTick',1:size(mat,1))
                end
                colorbar
                title ([measure_name,'   ---   Dissimilarity matrix'],'FontWeight','bold','FontSize',14)
            end
        end
    end
end
