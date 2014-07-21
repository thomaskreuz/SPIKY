% SPIKY_loop_examples --- Copyright Thomas Kreuz, Nebojsa Bozanic; March 2014
%
% 'SPIKY_loop_examples' demonstrates a few ways of how to use 'SPIKY_loop'.
%
% 'SPIKY_loop' is complementary to the graphical user interface 'SPIKY'.
% Both programs can be used to calculate time-resolved spike train distances (ISI and SPIKE) between two (or more) spike trains.
% However, whereas SPIKY was mainly designed to facilitate the detailed analysis of one dataset,
% 'SPIKY_loop' is meant to be used in order to compare the SPIKY_results for many different datasets (e.g. in some kind of loop).
% The source codes use a minimum number of input and output variables (described below).
% 'SPIKY_loop' is the main program where the variables are set and from where the funtion 'SPIKY_loop_f_distances' is called.
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


dataset=1;            % 1-Frequency mismatch,2-Spiking events,3-Splay state vs. identical,4-Clustering


m_para.all_measures_string={'ISI';'SPIKE';'SPIKE_realtime';'SPIKE_future';'PSTH';};  % order of select_measures

para.select_measures=[0 1 0 0 0];  % Select measures (0-calculate,1-do not calculate)


plotting=7;           % +1:spikes,+2:dissimilarity profile,+4:dissimilarity matrix


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
elseif dataset==3                % Splay state vs. identical  (from Fig. 1 of 2011 paper)
    para.tmin=0; para.tmax=para.tmin+800; para.dts=1;
    num_trains=20; num_spikes=10;
    spikes=cell(1,num_trains);
    spikes{1}(1:num_spikes)=(0:100:(num_spikes-1)*100);
    for trac=2:num_trains
        spikes{trac}(1)=para.tmin+spikes{1}(1);
        spikes{trac}(2:5)=(spikes{1}(1:4)+(trac-1)*100/num_trains);
        spikes{trac}(6:num_spikes)=para.tmin+spikes{1}(5:num_spikes-1);
        spikes{trac}=spikes{trac}(spikes{trac}<=para.tmax);
    end
elseif dataset==4
    para.tmin=0; para.tmax=4000; para.dts=1;
    para.all_train_group_names={'G1';'G2';'G3';'G4'};
    para.all_train_group_sizes=[10 10 10 10];

    num_trains=40; num_spikes=16;
    noise=[0.1 0.15 0.2 0.25 0.2 0.15 0.1 0.1];
    spikes=zeros(num_trains,num_spikes);
    for nc=1:num_spikes/8
        spikes(1:num_trains/2,nc)=(nc-0.5)/num_spikes*para.tmax+50*noise(1).*randn(1,num_trains/2)';
        spikes(num_trains/2+(1:num_trains/2),nc)=nc/num_spikes*para.tmax+50*noise(1).*randn(1,num_trains/2)';
    end
    for nc=num_spikes/8+(1:num_spikes/8)
        spikes(num_trains/4+(1:num_trains/2),nc)=(nc-0.5)/num_spikes*para.tmax+50*noise(2).*randn(1,num_trains/2)';
        spikes([1:num_trains/4 num_trains*3/4+(1:num_trains/4)],nc)=nc/num_spikes*para.tmax+50*noise(2).*randn(1,num_trains/2)';
    end
    for nc=num_spikes/4+(1:num_spikes/8)
        spikes([1:num_trains/4 num_trains/2+(1:num_trains/4)],nc)=(nc-0.5)/num_spikes*para.tmax+50*noise(3).*randn(1,num_trains/2)';
        spikes([num_trains/4+(1:num_trains/4) num_trains*3/4+(1:num_trains/4)],nc)=nc/num_spikes*para.tmax+50*noise(3).*randn(1,num_trains/2)';
    end
    rand_st=randperm(num_trains);
    for nc=num_spikes*3/8+(1:num_spikes/8)
        spikes(rand_st(1:num_trains/2),nc)=(nc-0.5)/num_spikes*para.tmax+50*noise(4).*randn(1,num_trains/2)';
        spikes(rand_st(num_trains/2+(1:num_trains/2)),nc)=nc/num_spikes*para.tmax+50*noise(4).*randn(1,num_trains/2)';
    end
    for nc=num_spikes/2+(1:num_spikes/8)
        spikes(1:num_trains/4,nc)=(nc-0.25)/num_spikes*para.tmax+50*noise(5).*randn(1,num_trains/4)';
        spikes(num_trains/4+(1:num_trains/2),nc)=(nc-0.5)/num_spikes*para.tmax+50*noise(5).*randn(1,num_trains/2)';
        spikes(num_trains*3/4+(1:num_trains/4),nc)=nc/num_spikes*para.tmax+50*noise(5).*randn(1,num_trains/4)';
    end
    spikes(spikes>0)=spikes(spikes>0)-60;
    for nc=num_spikes*5/8+(1:num_spikes/8)
        spikes(1:num_trains/4,nc)=nc/num_spikes*para.tmax+50*noise(6).*randn(1,num_trains/4)'-30;
        spikes(num_trains/4+(1:num_trains/4),nc)=(nc-0.25)/num_spikes*para.tmax+50*noise(6).*randn(1,num_trains/4)'-30;
        spikes(num_trains/2+(1:num_trains/4),nc)=(nc-0.5)/num_spikes*para.tmax+50*noise(6).*randn(1,num_trains/4)'-30;
        spikes(num_trains*3/4+(1:num_trains/4),nc)=(nc-0.75)/num_spikes*para.tmax+50*noise(6).*randn(1,num_trains/4)'-30;
    end
    for nc=num_spikes*6/8+(1:num_spikes/8)
        spikes(1:num_trains/8,nc)=(nc-0.11)/num_spikes*para.tmax+50*noise(7)/2.*randn(1,num_trains/8)';
        spikes(num_trains/8+(1:num_trains/8),nc)=(nc-0.22)/num_spikes*para.tmax+50*noise(7)/2.*randn(1,num_trains/8)';
        spikes(num_trains/4+(1:num_trains/8),nc)=(nc-0.33)/num_spikes*para.tmax+50*noise(7)/2.*randn(1,num_trains/8)';
        spikes(num_trains*3/8+(1:num_trains/8),nc)=(nc-0.44)/num_spikes*para.tmax+50*noise(7)/2.*randn(1,num_trains/8)';
        spikes(num_trains/2+(1:num_trains/8),nc)=(nc-0.55)/num_spikes*para.tmax+50*noise(7)/2.*randn(1,num_trains/8)';
        spikes(num_trains*5/8+(1:num_trains/8),nc)=(nc-0.66)/num_spikes*para.tmax+50*noise(7)/2.*randn(1,num_trains/8)';
        spikes(num_trains*3/4+(1:num_trains/8),nc)=(nc-0.77)/num_spikes*para.tmax+50*noise(7)/2.*randn(1,num_trains/8)';
        spikes(num_trains*7/8+(1:num_trains/8),nc)=(nc-0.88)/num_spikes*para.tmax+50*noise(7)/2.*randn(1,num_trains/8)';
    end
    for nc=num_spikes*7/8+(1:num_spikes/8)
        spikes(1:num_trains,nc)=nc/num_spikes*para.tmax-para.tmax/num_spikes.*rand(1,num_trains)';
    end

    para.instants=250:500:3750;
    para.selective_averages={[0 4000];...
        [0 500]; [500 1000]; [1000 1500]; [1500 2000]; [2000 2500]; [2500 3000]; [3000 3500]; [3500 4000];...
        [0 1000]; [500 1500]; [1000 2000]; [1500 2500]; [2000 3000]; [2500 3500]; [3000 4000];...
        [0 500 1000 1500]; [500 1000 1500 2000]; [1000 1500 2000 2500]; [1500 2000 2500 3000]; [2000 2500 3000 3500]; [2500 3000 3500 4000];...
        [0 500 1000 1500 2000 2500 3000 3500]; [500 1000 1500 2000 2500 3000 3500 4000];...     % Selected average over different intervals
        [0 4000]};
    tracs=1:num_trains;
    para.triggered_averages=cell(1,length(tracs));
    for trac=1:length(tracs)
        num_spikes=find(spikes(tracs(trac),:),1,'last');
        para.triggered_averages{trac}=round(spikes(tracs(trac),1:num_spikes)/para.dts)*para.dts;       % Triggered averaging over all time instants when a certain neuron fires
    end

    para.thin_separators=[];
    para.thick_separators=[];
    para.thin_markers= (500:500:para.tmax-500);
    para.thick_markers=[];

    %para.interval_names={'2 Cluster - AABB';'2 Cluster - ABBA';'2 Cluster - ABAB';'2 Cluster - Random association';...
    %    '3 Cluster - ABBC';'4 Cluster - ABCD';'8 Cluster - ABCDEFGH';'Random Spiking'};
    %para.interval_divisions=500:500:para.tmax-500; % Edges of subsections
elseif dataset==5
    MatInputStruct = load('IM334_IM335_spike_and_box_ts.mat');   %# pass file path as string
    StructName = fieldnames(MatInputStruct);

    % Joshua exception
    VariableNames = {['IM334_unit']};
    % also we should consider substructers
    VariableNameSubStruct = ['ts'];

    spikes = [];
    for cntr_vars = 1 : numel(VariableNames)

        spikesLocal = cell(size(MatInputStruct.(VariableNames{cntr_vars})));

        for cntr = 1 : size(MatInputStruct.(VariableNames{cntr_vars}), 2)
            spikesLocal{cntr} = single(MatInputStruct.(VariableNames{cntr_vars})(cntr).(VariableNameSubStruct));
        end

        spikes = [spikes spikesLocal];
    end

    clear MatInputStruct StructName FileName PathName VariableNames VariableNameSubStruct spikesLocal
end

d_para=para;
SPIKY_check_spikes
if ret==1
    return
end
para=d_para;

% ################################################### Actual call of function !!!!!

SPIKY_results = SPIKY_loop_f_distances(spikes,para) %#ok<NOPTS>

% ################################################### Example plotting (just as a demonstration)

num_plots=(mod(plotting,2)>0)+(mod(plotting,4)>1)+(mod(plotting,8)>3);
if num_plots>0
    measures=find(para.select_measures);
    for mc=1:length(measures)
        measure=measures(mc);
        measure_var=m_para.all_measures_string{measure};
        measure_name=regexprep(measure_var,'_','-');

        figure(mc); clf
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
            title ([measure_name,'   ---   Dissimilarity matrix'],'FontWeight','bold','FontSize',14)
        end
    end
end

