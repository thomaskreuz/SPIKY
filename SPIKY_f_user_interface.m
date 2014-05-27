% This function is the only file you might need to deal with. It consists of two parts. 
% The upper part can be used to define some of the standard parameters of SPIKY. The lower part
% can be used to generate predefined spike train datasets which can then be called via the listbox in the ?Selection: Data? panel.
% Always make sure that the variable 'listbox_str' labels correctly the datasets defined in the subsequent indices.

function [spikes,d_para,f_para,d_para_default,f_para_default,s_para_default,p_para_default,listbox_str]=SPIKY_f_user_interface(d_para,f_para,calc)

% structure 'd_para': parameters that describe the data, for a description see comments at the end of this file
d_para_default=struct('tmin',[],'tmax',[],'dts',[],'max_total_spikes',100000,'max_memo_init',10000000,'num_trains',[],...
    'select_train_mode',1,'select_train_groups',[],'preselect_trains',[],'latency_onset',[],...
    'num_all_train_groups',1,'all_train_group_names',[],'all_train_group_sizes',[],...
    'thick_separators',[],'thin_separators',[],'thick_markers',[],'thin_markers',[],'interval_divisions',[],'interval_names','',...
    'instants_str','','selective_averages_str','','triggered_averages_str','','spikes_variable','',...
    'comment_string','SPIKY','example',3);

% structure 'f_para': parameters that determine the appearance of the figures (and the movie), for a description see comments at the end of this file
f_para_default=struct('imagespath',['.',filesep],'moviespath',['.',filesep],...    % Default values
    'matpath',['.',filesep],'print_mode',0,'movie_mode',0,'publication',0,'comment_string','','num_fig',123,...
    'pos_fig',[0.5 0.01 0.5 0.8867],'supo1',[0.1 0.1 0.8 0.8],'hints',0,'show_title',1,'edge_correction',1,...
    'time_unit_string','','x_offset',0,'x_scale',1,'x_realtime_mode',0,'extreme_spikes',1,'ma_mode',1,'mao',20,...
    'psth_window',0,'psth_num_bins',1000,'frames_per_second',1,'num_average_frames',1,'profile_mode',1,...
    'profile_norm_mode',1,'profile_average_line',0,'color_norm_mode',1,'colorbar',1,'group_matrices',0,'dendrograms',0,...
    'histogram',0,'spike_train_color_coding_mode',1,'spike_col',2,'subplot_size',[],'plot_mode',5,...
    'subplot_posi',[0 1  0  0    2 0 0]);
% subplot_posi: Stim Spikes    PSTH   ISI   SPIKE SPIKE-realtime SPIKE-future

% SPIKE-pico SPIKE-rt-pico SPIKE-fut-pico   Vic vR
f_para_default.subplot_posi(8:12)=[0 0 0   0 0];  % please ignore, for testing purposes only
f_para_default.rel_subplot_size=ones(1,length(f_para_default.subplot_posi(f_para_default.subplot_posi>0)));
f_para_default.run_test=0;

f_para_default.regexp_str_scalar_positive_integer='[^1234567890]'; f_para_default.regexp_str_scalar_integer='[^-1234567890]';
f_para_default.regexp_str_scalar_positive_float='[^1234567890-e\.]'; f_para_default.regexp_str_scalar_float='[^-1234567890e\.]';
f_para_default.regexp_str_vector_positive_integers='[^1234567890: ]'; f_para_default.regexp_str_vector_integers='[^-1234567890: ]';
f_para_default.regexp_str_vector_positive_floats='[^1234567890:e \.]'; f_para_default.regexp_str_vector_floats='[^-1234567890:e \.]';
f_para_default.regexp_str_cell_floats='[^-1234567890:e[]{},; \.]';
% structure 's_para': parameters that describe the appearance of the individual subplots (measure time profiles)
% for a description see comments at the end of this file
% (this structure is less relevant for you since most of these parameters are set automatically)
s_para_default=struct('window_mode',1,'nma',1,'causal',1,'itmin',[],'itmax',[],...
    'num_subplots',[],'xl',[],'yl',[],'plot_mode',1);

% structure 'p_para': parameters that describe the appearance of the individual plot elements (described at the end of each line)
p_para_default=struct(...
    'stim_vis','on','stim_col','k','stim_ls','-','stim_lw',1,...                                % lines: stimulus
    'spike_vis','on','spike_col','k','spike_ls','-','spike_lw',1,'spike_marked_col','r',...     % lines: spikes (raster plot)
    'tbounds_vis','on','tbounds_col','k','tbounds_ls','-.','tbounds_lw',1,...                   % lines: overall time bounds
    'sp_seps_vis','on','sp_seps_col','k','sp_seps_ls','-','sp_seps_lw',2,...                    % lines: subplot separators
    'stim_bounds_vis','on','stim_bounds_col','k','stim_bounds_ls',':','stim_bounds_lw',1,...    % lines: stimulus time bounds
    'sp_bounds_vis','on','sp_bounds_col','k','sp_bounds_ls',':','sp_bounds_lw',1,...            % lines: subplot bounds
    'st_sep_vis','on','st_sep_col','k','st_sep_ls',':','st_sep_lw',1,...                        % lines: spike train separators
    'onset_vis','on','onset_col','b','onset_ls','-','onset_lw',2,...                            % lines: onset / offset
    'thick_sep_vis','on','thick_sep_col','r','thick_sep_marked_col','b','thick_sep_ls','-','thick_sep_lw',2,... % lines: thick spike train separators
    'thin_sep_vis','on','thin_sep_col','r','thin_sep_marked_col','b','thin_sep_ls','--','thin_sep_lw',1,...     % lines: thin spike train separators
    'thick_mar_vis','on','thick_mar_col','r','thick_mar_marked_col','b','thick_mar_ls','-','thick_mar_lw',2,... % lines: thick time markers
    'thin_mar_vis','on','thin_mar_col','r','thin_mar_marked_col','b','thin_mar_ls','--','thin_mar_lw',1.5,...   % lines: thin time markers
    'sp_thick_mar_vis','on','sp_thick_mar_col','r','sp_thick_mar_ls','-','sp_thick_mar_lw',2,...% lines: thick subplot time markers
    'sp_thin_mar_vis','on','sp_thin_mar_col','r','sp_thin_mar_ls','--','sp_thin_mar_lw',1.5,... % lines: thin subplot time markers
    'sgs_vis','on','sgs_col','k','sgs_marked_col','r','sgs_ls',':','sgs_lw',2,...               % lines: spike train group separators
    'mat_sgs_vis','on','mat_sgs_col','w','mat_sgs_ls',':','mat_sgs_lw',1,...                    % lines: matrix spike train group separators
    'mat_thick_sep_vis','on','mat_thick_sep_col','k','mat_thick_sep_ls','-','mat_thick_sep_lw',1,...% lines: thick matrix spike train separators
    'mat_thin_sep_vis','on','mat_thin_sep_col','k','mat_thin_sep_ls',':','mat_thin_sep_lw',1,...% lines: thin matrix spike train separators
    'prof_vis','on','prof_col','k','prof_ls','-','prof_lw',1,...                                % lines: measure profiles
    'ma_prof_vis','on','ma_prof_col','c','ma_prof_ls','-','ma_prof_lw',2,...                    % lines: measure profiles (analysis window)
    'ave_vis','on','ave_col','b','ave_ls','--','ave_lw',2,...                                   % lines: average of measure profiles
    'com_vis','on','com_col','k','com_ls',':','com_lw',1,...                                    % lines: common spikes
    'extreme_vis','on','extreme_col','k','extreme_ls',':','extreme_lw',1,...                    % lines: extrems spikes
    'dendrol_vis','on','dendrol_ls','-','dendrol_lw',2,...                                      % lines: dendrogram lines
    'mov_vis','on','mov_col','g','mov_ls','-','mov_lw',2,...                                    % lines: moving line
    'instants_vis','on','instants_col','g','instants_marked_col','r','instants_ls','-','instants_lw',2,...% lines: instants
    'selave_vis','on','selave_col','g','selave_active_col','b','selave_marked_col','r','selave_overlap_col','m','selave_ls','-','selave_lw',3,...% lines: selective averaging lines
    'trigave_vis','on','trigave_col','g','trigave_active_col','b','trigave_marked_col','r','trigave_ls','none','trigave_lw',3,...
    'trigave_symb_top','v','trigave_symb_bot','^','trigave_symb','+',...                        % plot: triggered averaging symbols
    'title_vis','on','title_col','k','title_fs',18,'title_fw','bold','title_fa','normal',...                              % texts: title
    'xlab_vis','on','xlab_col','k','xlab_fs',15,'xlab_fw','bold','xlab_fa','normal',...                                   % texts: profile x-label
    'prof_title_vis','on','prof_title_col','k','prof_title_fs',14,'prof_title_fw','bold','prof_title_fa','normal',...     % texts: profile title
    'prof_ave_vis','on','prof_ave_col','b','prof_ave_fs',16,'prof_ave_fw','bold','prof_ave_fa','normal',...               % texts: profile average
    'prof_tick_vis','on','prof_tick_col','k','prof_tick_fs',13,'prof_tick_fw','normal','prof_tick_fa','normal',...        % texts: profile ticks
    'mat_title_vis','on','mat_title_col','k','mat_title_fs',14,'mat_title_fw','bold','mat_title_fa','normal',...          % texts: matrix title
    'mat_label_vis','on','mat_label_col','k','mat_label_fs',14,'mat_label_fw','bold','mat_label_fa','normal',...          % texts: matrix label
    'mat_tick_vis','on','mat_tick_col','k','mat_tick_fs',12,'mat_tick_fw','normal','mat_tick_fa','normal',...             % texts: matrix ticks
    'measure_vis','on','measure_col','k','measure_fs',16,'measure_fw','bold','measure_fa','normal',...                    % texts: measure names
    'stimf_vis','on','stimf_col','k','stimf_fs',13,'stimf_fw','bold','stimf_fa','normal',...                              % texts: stimulus
    'group_names_vis','on','group_names_col','k','group_names_fs',15,'group_names_fw','bold','group_names_fa','normal',...% texts: group names
    'hist_max_vis','on','hist_max_col','k','hist_max_fs',12,'hist_max_fw','bold','hist_max_fa','normal',...               % texts: histogram maximum
    'mat_height',0.2,'mat_width',0.2,... % subplots: images
    'image_vis','on',...                 % subplots: images
    'colpat_vis','on'...                 % spike train group color patches
);


spikes=[];
listbox_str={'Frequency mismatch';...   % 1
    'Spiking events';...                % 2
    'Clustering';...                    % 3
    'Non-spurious events';...           % 4
    'Poisson Divergence';...            % 5
    'Splay state vs. identical';...     % 6
    'Edge-Test';...                     % 7
    'Testfile-Txt';...                  % 8
    'Testfile-Mat-ca';...               % 9
    'Testfile-Mat-zp';...               % 10
    'Testfile-Mat-01';...               % 11
};
%    'Ladder';...                        % 12
%    'paolo_spikes';...                  % 12

if calc
    switch d_para.sel_index

        case 1                                % Bi: Frequency mismatch

            %d_para.tmin=1000;
            d_para.tmin=0;
            d_para.tmax=d_para.tmin+1300;
            %d_para.tmax=d_para.tmin+600;
            d_para.dts=1;
            
            num_trains=2;
            spikes=cell(1,num_trains);
            spikes{1}=d_para.tmin+(100:100:1200); %+rand(1,num_spikes)*50;
            spikes{2}=d_para.tmin+(100:110:1200);
            %spikes{1}=d_para.tmin+(100:100:800); %+rand(1,num_spikes)*50;
            %spikes{1}=[100 200 310 400 510 600 720 800 930 1000 1110 1200];
            %spikes{2}=[100 210 330 410 530 640 760 870 970 1080 1200];
            %spikes{1}=[100 200 310 400 510 600 720 800];
            %spikes{2}=[100 210 330 410 530 640 760 870 970 1080 1200];
            %spikes{1}=[100 200 310];
            %spikes{2}=[100 210 330 410 530];
            %spikes{1}=[310 400 510];
            %spikes{2}=[100 210 330 410 530];

            %d_para.instants=d_para.tmin+[220 350 1070];

            %d_para.selective_averages{1}=d_para.tmin+[250 350];
            %d_para.selective_averages{1}=[d_para.tmin d_para.tmax];
            
            %d_para.triggered_averages=cell(1);
            %d_para.triggered_averages{1}=d_para.tmin+[280 730 900 1070];
            
            d_para.comment_string='SPIKY_Bi-Frequency-Mismatch';
            
            %f_para.plot_mode=1;

        case 2                                % Multi: decreasing noise + 5 events

            d_para.tmin=0;
            d_para.tmax=4000;
            d_para.dts=1;

            num_trains=50; num_spikes=40;               % Set spikes
            noise=[1:-1/(num_spikes/4-1):0 0];
            num_noises=length(noise);
            num_events=5;
            spikes=cell(1,num_trains);
            for trac=1:num_trains
                spikes{trac}=sort(rand(1,num_spikes/2),2)*d_para.tmax/2;
                for nc=1:num_events
                    spikes{trac}=[spikes{trac} nc*d_para.tmax/2/num_events+50*noise(ceil(num_noises-(nc-1)*num_noises/num_events)).*randn];
                end
                spikes{trac}=[spikes{trac} 100*(num_spikes/2-1)+200*(1:num_spikes/4+1)+50*noise.*randn(1,num_spikes/4+1)];
            end

            d_para.comment_string='SPIKY_Multi-Events';

        case 3                                      % Clustering
            
            parts=3;   % 1-first half,2-second half,3-all

            d_para.tmin=0;
            d_para.tmax=4000;
            d_para.dts=1;

            d_para.all_train_group_names={'G1';'G2';'G3';'G4'};
            d_para.all_train_group_sizes=[10 10 10 10];
            %d_para.num_all_train_groups=length(d_para.all_train_group_sizes);

            num_trains=40; num_spikes=16;
            noise=[0.1 0.15 0.2 0.25 0.2 0.15 0.1 0.1];

            matspikes=zeros(num_trains,num_spikes);
            %d_para.interval_names{1}='2 Cluster - AABB';
            for nc=1:num_spikes/8
                matspikes(1:num_trains/2,nc)=(nc-0.5)/num_spikes*d_para.tmax+50*noise(1).*randn(1,num_trains/2)';
                matspikes(num_trains/2+(1:num_trains/2),nc)=nc/num_spikes*d_para.tmax+50*noise(1).*randn(1,num_trains/2)';
            end

            %d_para.interval_names{2}='2 Cluster - ABBA';
            for nc=num_spikes/8+(1:num_spikes/8)
                matspikes(num_trains/4+(1:num_trains/2),nc)=(nc-0.5)/num_spikes*d_para.tmax+50*noise(2).*randn(1,num_trains/2)';
                matspikes([1:num_trains/4 num_trains*3/4+(1:num_trains/4)],nc)=nc/num_spikes*d_para.tmax+50*noise(2).*randn(1,num_trains/2)';
            end

            %d_para.interval_names{3}='2 Cluster - ABAB';
            for nc=num_spikes/4+(1:num_spikes/8)
                matspikes([1:num_trains/4 num_trains/2+(1:num_trains/4)],nc)=(nc-0.5)/num_spikes*d_para.tmax+50*noise(3).*randn(1,num_trains/2)';
                matspikes([num_trains/4+(1:num_trains/4) num_trains*3/4+(1:num_trains/4)],nc)=nc/num_spikes*d_para.tmax+50*noise(3).*randn(1,num_trains/2)';
            end

            %d_para.interval_names{4}='2 Cluster - Random association';
            rand_st=randperm(num_trains);
            for nc=num_spikes*3/8+(1:num_spikes/8)
                matspikes(rand_st(1:num_trains/2),nc)=(nc-0.5)/num_spikes*d_para.tmax+50*noise(4).*randn(1,num_trains/2)';
                matspikes(rand_st(num_trains/2+(1:num_trains/2)),nc)=nc/num_spikes*d_para.tmax+50*noise(4).*randn(1,num_trains/2)';
            end

            %d_para.interval_names{5}='3 Cluster - ABBC';
            for nc=num_spikes/2+(1:num_spikes/8)
                matspikes(1:num_trains/4,nc)=(nc-0.25)/num_spikes*d_para.tmax+50*noise(5).*randn(1,num_trains/4)';
                matspikes(num_trains/4+(1:num_trains/2),nc)=(nc-0.5)/num_spikes*d_para.tmax+50*noise(5).*randn(1,num_trains/2)';
                matspikes(num_trains*3/4+(1:num_trains/4),nc)=nc/num_spikes*d_para.tmax+50*noise(5).*randn(1,num_trains/4)';
            end
            matspikes(matspikes>0)=matspikes(matspikes>0)-60;

            %d_para.interval_names{6}='4 Cluster - ABCD';
            for nc=num_spikes*5/8+(1:num_spikes/8)
                matspikes(1:num_trains/4,nc)=nc/num_spikes*d_para.tmax+50*noise(6).*randn(1,num_trains/4)'-30;
                matspikes(num_trains/4+(1:num_trains/4),nc)=(nc-0.25)/num_spikes*d_para.tmax+50*noise(6).*randn(1,num_trains/4)'-30;
                matspikes(num_trains/2+(1:num_trains/4),nc)=(nc-0.5)/num_spikes*d_para.tmax+50*noise(6).*randn(1,num_trains/4)'-30;
                matspikes(num_trains*3/4+(1:num_trains/4),nc)=(nc-0.75)/num_spikes*d_para.tmax+50*noise(6).*randn(1,num_trains/4)'-30;
            end

            %d_para.interval_names{7}='8 Cluster - ABCDEFGH';
            for nc=num_spikes*6/8+(1:num_spikes/8)
                matspikes(1:num_trains/8,nc)=(nc-0.11)/num_spikes*d_para.tmax+50*noise(7)/2.*randn(1,num_trains/8)';
                matspikes(num_trains/8+(1:num_trains/8),nc)=(nc-0.22)/num_spikes*d_para.tmax+50*noise(7)/2.*randn(1,num_trains/8)';
                matspikes(num_trains/4+(1:num_trains/8),nc)=(nc-0.33)/num_spikes*d_para.tmax+50*noise(7)/2.*randn(1,num_trains/8)';
                matspikes(num_trains*3/8+(1:num_trains/8),nc)=(nc-0.44)/num_spikes*d_para.tmax+50*noise(7)/2.*randn(1,num_trains/8)';
                matspikes(num_trains/2+(1:num_trains/8),nc)=(nc-0.55)/num_spikes*d_para.tmax+50*noise(7)/2.*randn(1,num_trains/8)';
                matspikes(num_trains*5/8+(1:num_trains/8),nc)=(nc-0.66)/num_spikes*d_para.tmax+50*noise(7)/2.*randn(1,num_trains/8)';
                matspikes(num_trains*3/4+(1:num_trains/8),nc)=(nc-0.77)/num_spikes*d_para.tmax+50*noise(7)/2.*randn(1,num_trains/8)';
                matspikes(num_trains*7/8+(1:num_trains/8),nc)=(nc-0.88)/num_spikes*d_para.tmax+50*noise(7)/2.*randn(1,num_trains/8)';
            end

            %d_para.interval_names{8}='Random Spiking';
            for nc=num_spikes*7/8+(1:num_spikes/8)
                matspikes(1:num_trains,nc)=nc/num_spikes*d_para.tmax-d_para.tmax/num_spikes.*rand(1,num_trains)';
            end
            
            if parts==1
                d_para.tmin=0;
                d_para.tmax=2000;
            elseif parts==2
                d_para.tmin=2000;
                d_para.tmax=4000;
            end
            spikes=cell(1,num_trains);
            for trac=1:num_trains
                spikes{trac}=matspikes(trac,(matspikes(trac,:)>d_para.tmin & matspikes(trac,:)<d_para.tmax))-d_para.tmin;
            end
            d_para.tmax=d_para.tmax-d_para.tmin;
            d_para.tmin=0;
            clear matspikes
            %spikes{14}=[];
            %spikes{14}=spikes{14}(1:7);

            %d_para.instants= [250 1250 2750];
            %d_para.instants= (1:4000);            
            d_para.instants= (250:500:3750);

            d_para.selective_averages={ [0 4000];...
                 [0 500]; [500 1000]; [1000 1500]; [1500 2000]; [2000 2500]; [2500 3000]; [3000 3500]; [3500 4000];...
                 [0 1000]; [500 1500]; [1000 2000]; [1500 2500]; [2000 3000]; [2500 3500]; [3000 4000];...
                 [0 500 1000 1500]; [500 1000 1500 2000]; [1000 1500 2000 2500]; [1500 2000 2500 3000]; [2000 2500 3000 3500]; [2500 3000 3500 4000];...
                 [0 500 1000 1500 2000 2500 3000 3500]; [500 1000 1500 2000 2500 3000 3500 4000];...     % Selected average over different intervals
                 [0 4000]};
            %d_para.selective_averages={ [0 d_para.tmax-d_para.tmin]; [0 500]; [0 500 1000 1500]};
            %d_para.selective_averages={ [0 500]};
            %d_para.selective_averages=[];

            %d_para.selective_averages_str='{[610 740]};';     % Selected average over different intervals
            %d_para.selective_averages{1}=[610 740];
            %d_para.selective_averages{2}=[0 4000];

            % tracs=1+[0 1 2 3]*num_trains/4;
            %tracs=1;
            tracs=1:num_trains;
            d_para.triggered_averages=cell(1,length(tracs));
            for trac=1:length(tracs)
                %num_spikes=length(spikes{tracs(trac)});
                d_para.triggered_averages{trac}=round(spikes{tracs(trac)}/d_para.dts)*d_para.dts;       % Triggered averaging over all time instants when a certain neuron fires
            end
            %d_para.triggered_averages{trac}=d_para.triggered_averages{trac}(1:2);
            %d_para.triggered_averages{7}=round(spikes{tracs(7)}(1:2)/d_para.dts)*d_para.dts;       % Triggered averaging over all time instants when a certain neuron fires
            %d_para.triggered_averages{2}= [1200 1400];
            %d_para.triggered_averages=[];
            %d_para.triggered_averages_str='{[1300 1800 2300]; [1200 2200 3200]; [1750 2050 2750]; [1234 2345]}';
            %d_para.triggered_averages_str='SPIKY_get_trig_ave';
            %d_para.triggered_averages_str='';

            %d_para.instants=[2740];
            %d_para.selective_averages={[0 4000]; 2740+[0 5]};
            %d_para.triggered_averages=cell(1);
            %d_para.triggered_averages{1}=2740;

            d_para.thin_separators=[];
            d_para.thick_separators=[];
            d_para.thin_markers=500:500:d_para.tmax-500;
            d_para.thick_markers=[];
            %d_para.thick_markers=[];
            %d_para.thin_separators=[5 15 25 35];
            %d_para.thick_separators=[20];

            d_para.interval_names={'2 Cluster - AABB';'2 Cluster - ABBA';'2 Cluster - ABAB';'2 Cluster - Random association';...
                '3 Cluster - ABBC';'4 Cluster - ABCD';'8 Cluster - ABCDEFGH';'Random Spiking'};
            d_para.interval_divisions=500:500:d_para.tmax-500; % Edges of subsections
            d_para.comment_string='SPIKY_Clustering';
            f_para.comment_string='';
            f_para.spike_train_color_coding_mode=2;
            %f_para.group_matrices=0;
            %f_para.dendrograms=0;

        case 4                                % Non-spurious events

            d_para.tmin=0;
            d_para.tmax=400;
            d_para.dts=1;

            num_trains=50;
            num_spikes=7;

            thr1=100; % thr2=thr1+200;

            balance=0;

            while balance~=1

                randy=rand(1,num_trains);
                matspikes=zeros(num_trains,num_spikes);
                matspikes(1:num_trains,1)=0.001;
                matspikes(1:num_trains,2)=1+randy*80;
                matspikes(1:num_trains,3)=repmat(thr1,num_trains,1)+10*(rand(num_trains,1)-0.5);
                matspikes(1:num_trains,4)=200;
                matspikes(1:num_trains,5)=201+randy*80;
                matspikes(1:num_trains,num_spikes)=d_para.tmax; %-0.001;

                matspikes(:,num_spikes-1)=matspikes(:,3)+200;
                matspikes((matspikes(:,3)>thr1),num_spikes-1)=0;

                larger=sum(matspikes(:,num_spikes-1)>0);
                balance=larger*2/num_trains;

                for trac=1:num_trains
                    dummy=sort(matspikes(trac,1:num_spikes));
                    matspikes(trac,1:num_spikes)=[dummy(dummy>0) dummy(dummy==0)];
                end
                if balance==1
                    break
                end
            end
            spikes=cell(1,num_trains);
            for trac=1:num_trains
                spikes{trac}=matspikes(trac,:);
            end
            clear matspikes
            d_para.comment_string='SPIKY_Non-spurious';

        case 5                                % Poisson (Divergence)

            num_all_spikes=200;
            num_trig_trac1_spikes=12;
            %num_all_spikes=10; num_trig_trac1_spikes=3;

            num_trains=20;
            trig_trac1=1;
            trig_tracs1=[4 8 11 16 19];
            d_para.tmin=0;
            d_para.tmax=100; %min(matspikes(:,num_all_spikes))*1.0001;
            d_para.dts=0.001;

            matspikes=zeros(num_trains,num_all_spikes);
            for trac=setdiff(1:num_trains,trig_trac1)
                dummy=SPIKY_f_poisson(num_all_spikes,1,0);
                dummy(dummy<d_para.dts)=d_para.dts;
                matspikes(trac,1:num_all_spikes)=cumsum(dummy);
            end
            dummy=SPIKY_f_poisson(num_trig_trac1_spikes,1,0);
            dummy(dummy<d_para.dts)=d_para.dts;
            matspikes(trig_trac1,1:num_trig_trac1_spikes)=cumsum(dummy);

            matspikes(matspikes<d_para.tmin | matspikes>d_para.tmax)=0;

            num_spikes=zeros(1,num_trains);
            for trac=1:num_trains
                num_spikes(trac)=find(matspikes(trac,1:num_all_spikes),1,'last');
            end

            matspikes(trig_trac1,1:num_trig_trac1_spikes)=matspikes(trig_trac1,1:num_trig_trac1_spikes)/max(matspikes(trig_trac1,1:num_trig_trac1_spikes))*d_para.tmax*0.97;
            for trac=setdiff(1:num_trains,trig_trac1)
                matspikes(trac,1:num_spikes(trac))=matspikes(trac,1:num_spikes(trac))/max(matspikes(trac,1:num_spikes(trac)))*d_para.tmax*0.995;
                matspikes(trac,num_spikes(trac))=matspikes(trac,num_spikes(trac))-rand*5;
                matspikes(trac,1:num_spikes(trac))=sort(matspikes(trac,1:num_spikes(trac)));
            end

            max_num_spikes=max(num_spikes);
            matspikes=matspikes(:,1:max_num_spikes);

            for trac=trig_tracs1
                indy=[];
                for spic=1:num_spikes(trig_trac1)
                    rem_indy=setdiff(1:max_num_spikes,indy);
                    [dummy,index]=min(abs(matspikes(trac,rem_indy)-matspikes(trig_trac1,spic)));
                    indy=[indy rem_indy(index)];
                    matspikes(trac,rem_indy(index))=matspikes(trig_trac1,spic)+0.05*rand;
                end
                matspikes(trac,1:num_spikes(trac))=sort(matspikes(trac,1:num_spikes(trac)));
            end
            matspikes(matspikes<d_para.tmin | matspikes>d_para.tmax)=0;

            spikes=cell(1,num_trains);
            for trac=1:num_trains
                spikes{trac}=matspikes(trac,:);
            end
            clear matspikes
            
            %num_spikes2=cellfun('length',spikes);
            d_para.comment_string='SPIKY_Poisson-Driver';

        case 6                                % 2011 paper, Fig. 2a   (splay state vs. identical)

            d_para.tmin=0;
            d_para.tmax=800;
            d_para.dts=1;

            num_trains=20;
            num_spikes=10;
            spikes=cell(1,num_trains);
            spikes{1}(1:num_spikes)= (0:100:(num_spikes-1)*100);
            for trac=2:num_trains
                spikes{trac}(1)= spikes{1}(1);
                spikes{trac}(2:5)= (spikes{1}(1:4)+(trac-1)*100/num_trains);
                spikes{trac}(6:num_spikes)= spikes{1}(5:num_spikes-1);
                spikes{trac}=spikes{trac}(spikes{trac}<=d_para.tmax);
            end
            d_para.comment_string='SPIKY_Splay';

        case 7                                       % Edge-Test

            d_para.tmin=0;
            d_para.tmax=1;
            d_para.dts=0.001;
            num_trains=3;
            spikes=cell(1,num_trains);
            spikes{1}= [0.2 0.3 0.4 0.6 0.7 0.9];
            spikes{2}= [0.1 0.25 0.5 0.75 0.8];
            %spikes{1}= [0.2 0.5 0.9];
            %spikes{2}= [0.1 0.75];
            spikes{3}= [0.05 0.08 0.1 0.25 0.3 0.75 0.95 0.98];
            %spikes{4}= 0.75;

            d_para.instants= [0.6];

            %d_para.selective_averages={[d_para.tmin d_para.tmax]; [0 0.5]};
            d_para.selective_averages={ [0.15 0.6]};

            % d_para.triggered_averages=cell(1,3);
            % d_para.triggered_averages{1}= [0.1 0.4 0.7];
            % d_para.triggered_averages{2}= [0.1:0.1:0.8];
            % d_para.triggered_averages{3}= [0.6 0.7];

            d_para.comment_string='SPIKY_Edge-Test';

        case 8                                        % Testfile-Txt

            textfile='SPIKY_testdata.txt';
            spikes=dlmread(textfile);
            d_para.tmin=0;
            d_para.tmax=4000;
            d_para.dts=1;
            d_para.comment_string='Testfile-Txt';

        case 9                                        % Testfile-Mat (cell arrays)

            d_para.matfile='SPIKY_testdata_ca.mat';
            d_para.tmin=0;
            d_para.tmax=4000;
            d_para.dts=1;
            d_para.comment_string='Testfile-Mat-ca';

        case 10                                       % Testfile-Mat (matrix with zero padding)

            d_para.matfile='SPIKY_testdata_zp.mat';
            d_para.tmin=0;
            d_para.tmax=4000;
            d_para.dts=1;
            d_para.comment_string='Testfile-Mat-zp';

        case 11                                       % Testfile-Mat (matrix with 0 and 1)

            d_para.matfile='SPIKY_testdata_01.mat';
            d_para.tmin=0;
            d_para.tmax=4000;
            d_para.dts=1;
            d_para.comment_string='Testfile-Mat-01';


        case 12                                       % Ladder

            d_para.tmin=0;
            d_para.tmax=10000;
            d_para.dts=1;

            num_trains=30;               % Set spikes
            
            spikes=cell(1,num_trains);
            for trac=1:num_trains
                spikes{trac}=d_para.tmin+(0:num_trains)/num_trains*(d_para.tmax-d_para.tmin);
                spikes{trac}(num_trains+1+(1:trac))=d_para.tmin+(-0.5+(1:trac))/num_trains*(d_para.tmax-d_para.tmin);
                spikes{trac}=sort(spikes{trac});
            end

            d_para.comment_string='Ladder';

%         case 12                                       % Paolo data (matrix with 0 and 1)
% 
%             d_para.matfile='paolo_spikes.mat';
%             d_para.dts=0.025;
%             d_para.tmin=0;
%             d_para.tmax=70559*d_para.dts;
%             d_para.comment_string='Paolo';

        otherwise

            entries = d_para.entries;
            entry=char(entries(d_para.sel_index));
            if stracmp(entry(length(entry)-3:length(entry)),'.mat')
                d_para.matfile=entry;
            elseif stracmp(entry(length(entry)-3:length(entry)),'.txt')
                txtfile=entry;
                spikes=dlmread(txtfile);
            end
            d_para.comment_string=entry;
    end

    SPIKY_check_spikes
    if ret==1
        spikes=[];
        return
    end

end
end

%
% Here are short descriptions for some of the parameters that can be set in the beginning of this file and which are loaded
% when SPIKY is started (but not when it is reset, thus in order for changes to become active you have to restart SPIKY !!!)
%
%
% structure 'd_para': parameters that describe the data (some of these will automatically be determined from the data)
% ====================================================================================================================
%
% tmin: Beginning of recording, will automatically be determined from the data 
% tmax: End of recording, will automatically be determined from the data 
% dts: Sampling interval, precision of spike times, will automatically be determined from the data 
% max_total_spikes: Maximum number of spikes that will be plotted
% max_memo_init: Memory management, should be big enough to hold the basic matrices but small enough to not run out of memory
% num_trains: Number of spike trains, will automatically be determined from the data 
% select_train_mode: Selection of spike trains (1-all,2-selected groups,3-selected trains)
% select_train_groups: Selected spike train groups (if 'select_train_mode==2')
% preselect_trains: Selected spike trains (if 'select_train_mode==3')
% latency_onset: latency of first spike after this time instant, can be used as a criterion for sorting spike trains 
% num_all_train_groups: Number of spike train groups
% all_train_group_names: Names of spike train groups
% all_train_group_sizes: Sizes of respective spike train groups
% thick_separators: Very relevant seperations between groups of spike trains
% thin_separators: Relevant seperations between groups of spike trains
% thick_markers: Even more relevant time instants
% thin_markers: Relevant time instants
% interval_divisions: Edges of subsections, can be useful in titles of the movie
% interval_names: Captions for subsections, can be shown in titles of the movie (e.g., for time instants within these intervals)
% instants_str: One or more continuous intervals for selective temporal averaging
% select_averages_str: Time instants of interest
% trigger_averages_str: Non-continuous time-instants for triggered temporal averaging, external (e.g. certain stimulus properties) but also internal (e.g. certain event times)
% spikes_variable: Default name of variable / field that cpontains the spike times
% comment_string: Additional comment on the example, will be used in file and figure names
% example: Position of preselected item in the ListBox in the 'Selection: Data' panel
%
%
% structure 'f_para': parameters that determine the appearance of the figures (and the movie)
% ===========================================================================================
%
% imagespath: Path where images (postscript) will be stored
% moviespath: Path where movies (avi) will be stored
% matpath: Path where Matlab files (mat) will be stored
% print_mode: Saving to postscript file? (0-no,1-yes)
% movie_mode: Record a movie? (0-no,1-yes)
% publication: Omits otherwise helpful information, such as additional comments (0-no,1-yes)
% comment_string: Additional comment on the example, will be used in file and figure names
% num_fig: Number of figure
% pos_fig: Position of figure
% supo1: Position (in relative units) of the first subplot which contains spikes and dissimilarity profiles
% hints: Determines whether short hints will be shown when hovering with the mouse cursor above the SPIKY elements of interest
% show_title: Determines whether a figure title will be shown (above the first subplot containing the spikes) 
% edge_correction: Determines whether the edge effect (spurious drop to zero dissimilarity at the edges caused by the auxiliary spikes) is corrected or not
% time_unit_string: Time unit, used in labels
% x_offset: Offset for time axis (useful for example if you select an intermediate interval but want to have the time scale start at t=0) 
% x_scale: Scale factor of time unit
% x_realtime_mode: Sets reference system for the movie, if set during the movie the time position is kept constant while the spikes move
%       (simulating an online recorded stream of incoming spike times)
% extreme_spikes: Determines whether dotted lines are shown at the position of the extreme spikes
%       (the last first spike and the first last spike which indicate the onset of the edge effects)
% filename: Name under which images, movies and Matlab files will be stored
% ma_mode: Moving average mode: (0-no,1-only,2-both)
% mao: Order of moving average (piecewise constant and piecewise linear dissimilarity profiles)
% psth_window: Kernel width of Gaussian filter for the Peristimulus time histogram. Set to 0 for no filtering.
% psth_num_bins: Number of bins for the Peristimulus time histogram.
% frames_per_second: Well, frames per second for the movie
% num_average_frames: Number of frames the averages are shown at the end of the movie (if this is too small they are hardly visible)
% plot_mode: +1:profiles,+2:frame comparison,+4:frame sequence (the latter two are mutually exclusive, otherwise binary addition allows combinations)
% profile_mode: Allows to additionally/exclusively show dissimilarity profiles averaged over certain spike train groups only
% profile_norm_mode: Normalization of averaged bivariate dissimilarity profiles (1-Absolute maximum value 'one',2-Overall maximum value,3-Individual maximum value)
% profile_average_line: Determines whether a line at the mean value is shown for each dissimilarity profile
% color_norm_mode: normalization of pairwise color matrices (1-Absolute maximum value,2-Overall maximum value,3-Each one individually)
% colorbar: Determines whether a colorbar is shown next to the dissimilarity matrices
% group_matrices: Allows tracing the overall synchronization among groups of spike trains (0-no,1-yes)
% dendrograms: Show cluster trees from distance matrices (0-no,1-yes)
% histogram: Show histogram with number of spikes for each spike train (on the right of the spike trains)
% spike_train_color_coding_mode: Labeling of spike trains used also in the dendrograms (0-no,1-groups,2-trains)
% spike_col: color spikes according to group of train, depending on spike_train_color_coding_mode (0-no,+1-spike trains,+2-bars at the edges)
% subplot_size: Vector with relative size of subplots
% subplot_posi: Vector with order of subplots
%
%
% structure 's_para': parameters that describe the appearance of the individual subplots (measure time profiles)
% ==============================================================================================================
% (this structure is less relevant for you since most of these parameters are set automatically)
%
% window_mode: time interval selection (1-all (recording limits),2:outer spikes,3-inner spikes,4-smaller analysis window)
% nma: Selected moving averages (depends on f_para.ma_mode)
% causal: determines kind of moving average, set automatically for each measure (0-no,1-yes)
% itmin: Beginning of analysis window, set automatically
% itmax: End of analysis window, set automatically
% num_subplots: depends on f_para.subplot_posi, set automatically
% xl: x-limits for plotting, set automatically
% yl: y-limits for plotting, set automatically
% plot_mode: equals f_para.plot_mode
%

