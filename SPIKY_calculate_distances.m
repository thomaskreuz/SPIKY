% This calculates the selected time-resolved spike train distances once the button 'Calculate' is pressed.
% In case of very long datasets it also initiates the memory management.

disp(' '); disp(' '); disp(' ');

m_para.all_measures_str={'Stimulus';'Spikes';'PS';'I'; 'S';'S_r';'S_f'; 'Ss';'Ss_r';'SCs'; 'Si';'Si_r';'Si_f';'PT';'PN';'Comb';'V';'R'};
m_para.all_measures_string={'Stimulus';'Spikes';'PSTH';'ISI'; 'SPIKE';'SPIKE_realtime';'SPIKE_future';...
    'SPIKEsamp';'SPIKEsamp_realtime';'SPIKEsamp_future'; 'SPIKEi';'SPIKEi_realtime';'SPIKEi_future';'PT';'PN';'Comb';'Victor';'vanRossum'};

m_para.num_diff_measures=2;
m_para.num_non_bi_measures=3;

m_para.stimulus=1;
m_para.spikes=2;
m_para.psth=3;
m_para.isi_pico=4;
m_para.spike_pili=5;
m_para.realtime_spike_pili=6;
m_para.future_spike_pili=7;
m_para.spike_samp=8;
m_para.realtime_spike_samp=9;
m_para.future_spike_samp=10;
m_para.spike_pico=11;
m_para.realtime_spike_pico=12;
m_para.future_spike_pico=13;
m_para.test_pili_1=14;
m_para.test_pili_2=15;
m_para.test_pili_3=16;
m_para.victor=17;
m_para.van_rossum=18;

m_para.pili_measures=[m_para.spike_pili m_para.realtime_spike_pili m_para.future_spike_pili ...
    m_para.test_pili_1 m_para.test_pili_2 m_para.test_pili_3];
m_para.samp_measures=[m_para.spike_samp m_para.realtime_spike_samp m_para.future_spike_samp];
m_para.pico_measures=[m_para.isi_pico m_para.spike_pico m_para.realtime_spike_pico m_para.future_spike_pico];
m_para.realtime_measures=[m_para.realtime_spike_pico m_para.realtime_spike_samp m_para.realtime_spike_pili];
m_para.future_measures=[m_para.future_spike_pico m_para.future_spike_samp m_para.future_spike_pili];
m_para.bi_measures=[m_para.pico_measures m_para.samp_measures m_para.pili_measures];

f_para.subplot_posi(m_para.stimulus)=str2num(get(handles.subplot_stimulus_posi_edit,'String')); %#ok<*ST2NM>
f_para.subplot_posi(m_para.spikes)=str2num(get(handles.subplot_spikes_posi_edit,'String'));
f_para.subplot_posi(m_para.psth)=str2num(get(handles.subplot_psth_posi_edit,'String'));
f_para.subplot_posi(m_para.isi_pico)=str2num(get(handles.subplot_isi_posi_edit,'String'));
f_para.subplot_posi(m_para.spike_pili)=str2num(get(handles.subplot_spike_posi_edit,'String'));
f_para.subplot_posi(m_para.realtime_spike_pili)=str2num(get(handles.subplot_spike_realtime_posi_edit,'String'));
f_para.subplot_posi(m_para.future_spike_pili)=str2num(get(handles.subplot_spike_future_posi_edit,'String'));


m_para.as_vect=SPIKY_f_all_sorted(f_para.subplot_posi);    % position number, not just 0 and 1
if ~all(m_para.as_vect==f_para.subplot_posi)
    f_para.subplot_posi=m_para.as_vect;
    set(handles.subplot_stimulus_posi_edit,'String',num2str(f_para.subplot_posi(m_para.stimulus)))
    set(handles.subplot_spikes_posi_edit,'String',num2str(f_para.subplot_posi(m_para.spikes)))
    set(handles.subplot_psth_posi_edit,'String',num2str(f_para.subplot_posi(m_para.psth)))
    set(handles.subplot_isi_posi_edit,'String',num2str(f_para.subplot_posi(m_para.isi_pico)))
    set(handles.subplot_spike_posi_edit,'String',num2str(f_para.subplot_posi(m_para.spike_pili)))
    set(handles.subplot_spike_realtime_posi_edit,'String',num2str(f_para.subplot_posi(m_para.realtime_spike_pili)))
    set(handles.subplot_spike_future_posi_edit,'String',num2str(f_para.subplot_posi(m_para.future_spike_pili)))
end

m_para.num_all_measures=length(f_para.subplot_posi);
select_measures=intersect(m_para.num_diff_measures+1:m_para.num_all_measures,find(f_para.subplot_posi));
[dummy,ms_indy]=sort(f_para.subplot_posi(select_measures));
m_para.select_measures=select_measures(ms_indy);
m_para.select_bi_measures=m_para.select_measures(ismember(m_para.select_measures,m_para.bi_measures));
m_para.num_sel_measures=length(m_para.select_measures);
m_para.num_sel_bi_measures=length(m_para.select_bi_measures);

m_para.select_pili_measures=intersect(m_para.pili_measures,m_para.select_measures);
m_para.num_pili_measures=length(m_para.select_pili_measures);
m_para.select_samp_measures=intersect(m_para.samp_measures,m_para.select_measures);
m_para.num_samp_measures=length(m_para.select_samp_measures);
m_para.select_pico_measures=intersect(m_para.pico_measures,m_para.select_measures);
m_para.num_pico_measures=length(m_para.select_pico_measures);

m_para.measure_indy=zeros(1,m_para.num_all_measures);
m_para.measure_indy(m_para.pili_measures)=SPIKY_f_all_sorted(f_para.subplot_posi(m_para.pili_measures));
m_para.measure_indy(m_para.samp_measures)=SPIKY_f_all_sorted(f_para.subplot_posi(m_para.samp_measures));
m_para.measure_indy(m_para.pico_measures)=SPIKY_f_all_sorted(f_para.subplot_posi(m_para.pico_measures));

aaa=f_para.subplot_posi(m_para.num_diff_measures+1:m_para.num_all_measures);
bbb=aaa+(1:(m_para.num_all_measures-m_para.num_diff_measures))/(m_para.num_all_measures-m_para.num_diff_measures+1);
bbb(bbb<1)=0;
m_para.measure_all_indy=[zeros(1,m_para.num_diff_measures) SPIKY_f_all_sorted(bbb)];

ccc=aaa(aaa>0)-min(aaa(aaa>0))+1;
[ddd,eee]=sort(ccc);
fff=find(aaa);
m_res.mat_str=m_para.all_measures_str(m_para.num_diff_measures+fff(eee));

aaa=f_para.subplot_posi(m_para.num_non_bi_measures+1:m_para.num_all_measures);   % assumes 1 non-bi measure (psth) directly after num_diff
bbb=aaa+(1:(m_para.num_all_measures-m_para.num_non_bi_measures))/(m_para.num_all_measures-m_para.num_non_bi_measures+1);
bbb(bbb<1)=0;
m_para.measure_bi_indy=[zeros(1,m_para.num_non_bi_measures) SPIKY_f_all_sorted(bbb)];
%m_res.mat_str=m_para.all_measures_str(m_para.num_diff_measures+m_para.measure_all_indy(m_para.measure_all_indy>0));

ccc=aaa(aaa>0)-min(aaa(aaa>0))+1;
[ddd,eee]=sort(ccc);
fff=find(aaa);
m_res.bi_mat_str=m_para.all_measures_str(m_para.num_non_bi_measures+fff(eee));

m_para.sel_measures_str=m_para.all_measures_str(m_para.select_measures);
m_para.sel_bi_measures_str=m_para.all_measures_str(m_para.select_bi_measures);
m_para.sel_pili_measures_str=m_para.all_measures_str(m_para.select_pili_measures);
m_para.sel_samp_measures_str=m_para.all_measures_str(m_para.select_samp_measures);
m_para.sel_pico_measures_str=m_para.all_measures_str(m_para.select_pico_measures);

m_para.pili_measures_indy=find(ismember(m_para.select_measures,m_para.pili_measures));
m_para.samp_measures_indy=find(ismember(m_para.select_measures,m_para.samp_measures));
m_para.pico_measures_indy=find(ismember(m_para.select_measures,m_para.pico_measures));

m_para.pili_bi_measures_indy=find(ismember(m_para.select_bi_measures,m_para.pili_measures));
m_para.samp_bi_measures_indy=find(ismember(m_para.select_bi_measures,m_para.samp_measures));
m_para.pico_bi_measures_indy=find(ismember(m_para.select_bi_measures,m_para.pico_measures));

if str2double(get(handles.dpara_tmin_edit,'String'))>=str2double(get(handles.dpara_tmax_edit,'String'))
    set(0,'DefaultUIControlFontSize',16);
    mbh=msgbox(sprintf('The beginning of the recording can not be later\nthan the end of the recording!'),'Warning','warn','modal');
    htxt = findobj(mbh,'Type','text');
    set(htxt,'FontSize',12,'FontWeight','bold')
    mb_pos=get(mbh,'Position');
    set(mbh,'Position',[mb_pos(1:2) mb_pos(3)*1.5 mb_pos(4)])
    uiwait(mbh);
    set(handles.dpara_tmin_edit,'String',num2str(d_para.tmin))
    set(handles.dpara_tmax_edit,'String',num2str(d_para.tmax))
    ret=1;
    return;
end

%num_coins=cell(1,d_para.num_trains);
uspikes=cell(1,d_para.num_trains);
for trac=1:d_para.num_trains
    uspikes{trac}=spikes{trac}(spikes{trac}>=d_para.tmin & spikes{trac}<=d_para.tmax);
    uspikes{trac}=unique([d_para.tmin uspikes{trac} d_para.tmax]);
%     dummy=unique(uspikes{trac});
%     for uic=1:length(dummy)-1
%         num_coins{trac}(uic)=sum(uspikes{trac}==dummy(uic));
%     end
%     uspikes{trac}=dummy;
end
num_uspikes=cellfun('length',uspikes);
max_num_uspikes=max(num_uspikes);

% #########################################################################
% #########################################################################
% #########################################################################

if f_para.subplot_posi(m_para.victor)
    q=[0 5 10];
    m_res.vic=zeros(length(q),d_para.num_trains,d_para.num_trains);
    for trac1=1:d_para.num_trains-1
        for trac2=trac1+1:d_para.num_trains
            m_res.vic(:,trac1,trac2)=SPIKY_Victor_MEX(uspikes{trac1},uspikes{trac2},single(q));
            m_res.vic(1:length(q),trac2,trac1)=m_res.vic(1:length(q),trac1,trac2);
        end
    end
end
if f_para.subplot_posi(m_para.van_rossum)   % still problematic !!!!!!!
    tau=1;
    m_res.vr=zeros(length(tau),d_para.num_trains,d_para.num_trains);
    m_res.vr2=zeros(length(tau),d_para.num_trains,d_para.num_trains);
    m_res.vr3=zeros(length(tau),d_para.num_trains,d_para.num_trains);
    for tauc=1:length(tau)
        for trac1=1:d_para.num_trains-1
            for trac2=trac1+1:d_para.num_trains
                m_res.vr(tauc,trac1,trac2)=SPIKY_vanRossum(uspikes{trac1}/d_para.tmax*2,uspikes{trac2}/d_para.tmax*2,tau*2);
                m_res.vr2(tauc,trac1,trac2)=SPIKY_vanRossum(uspikes{trac1}/d_para.tmax,uspikes{trac2}/d_para.tmax,tau);
                m_res.vr3(tauc,trac1,trac2)=SPIKY_vanRossum(uspikes{trac1},uspikes{trac2},tau);
            end
        end
        m_res.vr(tauc,:,:)=shiftdim(m_res.vr(tauc,:,:),1)+shiftdim(m_res.vr(tauc,:,:),1)';
        m_res.vr2(tauc,:,:)=shiftdim(m_res.vr2(tauc,:,:),1)+shiftdim(m_res.vr2(tauc,:,:),1)';
        m_res.vr3(tauc,:,:)=shiftdim(m_res.vr3(tauc,:,:),1)+shiftdim(m_res.vr3(tauc,:,:),1)';   % still problematic !!!!!!!
    end
end

% #########################################################################
% #########################################################################
% #########################################################################

if any(f_para.subplot_posi(m_para.isi_pico:m_para.num_all_measures))                                                        % ISI-SPIKE

    dummy=[0 num_uspikes];
    all_indy=zeros(1,sum(num_uspikes));
    for trac=1:d_para.num_trains
        all_indy(sum(dummy(1:trac))+(1:num_uspikes(trac)))=trac*ones(1,num_uspikes(trac));
    end
    [all_spikes,indy]=sort([uspikes{:}]);
    all_trains=all_indy(indy);
    all_trains(1:d_para.num_trains)=0;
    all_trains(end-d_para.num_trains+1:end)=0;

    all_spikes=all_spikes(d_para.num_trains:end-d_para.num_trains+1);
    all_trains=all_trains(d_para.num_trains:end-d_para.num_trains+1);
    all_isi=zeros(1,length(all_spikes)-1,'single');
    all_isi(1:length(all_spikes)-1)=diff(all_spikes);
    num_all_isi=length(all_isi);
    m_res.isi=all_isi(all_isi>0);
    m_res.num_isi=length(m_res.isi);

    m_res.isi_pos=cumsum([d_para.tmin m_res.isi(1:end-1)])+m_res.isi/2;
    m_res.cum_isi=all_spikes(1)+[0 cumsum(m_res.isi)];
    clear all_indy indy all_spikes

    % We define the ISI as going from right after the last spike to the
    % exact time of the next spike. So the previous spike is not part of the ISI
    % while the following is. So when you are on a spike this is the
    % following spike and the one before is the previous spike.

    dummy=[spikes{:}];
    dummy=double(unique(dummy(dummy>d_para.tmin & dummy<d_para.tmax)));
    m_res.pili_supi=sort([d_para.tmin dummy dummy d_para.tmax]);
    m_res.pili_len=length(m_res.pili_supi);
    if m_res.pili_len~=2*m_res.num_isi
        disp(' '); disp(' ')
        error(['m_res.pili_len 2*m_res.num_isi =',num2str([m_res.pili_len 2*m_res.num_isi])])
    end
    m_res.pili_supi_indy=round((m_res.pili_supi-d_para.tmin)/d_para.dts);   % dtm
    isi_pili=reshape(repmat(m_res.isi,2,1),1,2*m_res.num_isi);
    isi_indy_pili=reshape(repmat(1:m_res.num_isi,2,1),1,2*m_res.num_isi);


    isis=cell(1,d_para.num_trains);
    ints=zeros(d_para.num_trains,num_all_isi,'single');
    for trac=1:d_para.num_trains
        isis{trac}=diff(uspikes{trac});
        %isis{trac}=isis{trac}(isis{trac}~=0)./double(num_coins{trac}(isis{trac}~=0));
        if f_para.edge_correction && num_uspikes(trac)>4                                                                 % $$$$$$$$$$$$
            isis{trac}(1)=max(isis{trac}(1:2));
            isis{trac}(end)=max(isis{trac}(end-1:end));
        end
        ivs=[1 find(all_trains==trac)];
        ive=[ivs(2:length(ivs))-1 num_all_isi];
        for ic=1:num_uspikes(trac)-1
            ints(trac,ivs(ic):ive(ic))=isis{trac}(ic);
        end
    end
    ints=ints(:,all_isi>0);
    clear isis % num_coins


    if f_para.subplot_posi(m_para.test_pili_2)
        %num_spikes=cellfun('length',spikes);
        first_spike=cellfun(@(x) x(1), spikes, 'UniformOutput',false);
        last_spike=cellfun(@(x) x(end), spikes, 'UniformOutput',false);

        train_spike_status=zeros(d_para.num_trains,m_res.num_isi);
        for trac=1:d_para.num_trains
            train_spike_status(trac,1:find(m_res.cum_isi==first_spike{trac})-1)=1;
            train_spike_status(trac,find(m_res.cum_isi==last_spike{trac}):m_res.num_isi)=3;
        end

        mat_indy=nchoosek(1:d_para.num_trains,2);

        pair_spike_status=zeros(d_para.num_pairs,m_res.num_isi);
        for pac=1:d_para.num_pairs
            pair_spike_status(pac,1:m_res.num_isi)=train_spike_status(mat_indy(pac,1),1:m_res.num_isi)+...
                train_spike_status(mat_indy(pac,2),1:m_res.num_isi);
        end
    end


    % ###########################################################################################################################################
    % ############################################################## Memory management ##########################################################
    % ###########################################################################################################################################

    %d_para.max_memo_init=100000;      % set in SPIKY_f_user_interface: memory management, should be big enough to hold the basic matrices but small enough to not run out of memory

    m_res.num_samples=fix((d_para.tmax-d_para.tmin)/d_para.dts);
    m_res.num_measure_samples=fix((d_para.tmax-d_para.tmin)/d_para.dtm);

    m_para.memo_num_pi_measures=2*m_para.num_pili_measures+m_para.num_pico_measures;

    memo_fact=(m_para.memo_num_pi_measures*m_res.num_isi+m_para.num_samp_measures*m_res.num_measure_samples);
    memo_pi_fact=m_para.memo_num_pi_measures*d_para.num_pairs;
    memo=memo_pi_fact*m_res.num_isi;

    run_test=0;
    if run_test==0 || ~ismember(get(handles.Data_listbox,'Value'),[1 2 3])
        max_memo=d_para.max_memo_init;
    elseif get(handles.Data_listbox,'Value')==1
        max_memo=30; % $$$$$$$$$$$
    elseif get(handles.Data_listbox,'Value')==2
        max_memo=800; % $$$$$$$$$$$
    elseif get(handles.Data_listbox,'Value')==3
        max_memo=800; % $$$$$$$$$$$
    end
    

    r_para.num_samp_runs=1;
    r_para.num_pi_runs=1;
    if memo>max_memo
        num_init_runs=ceil(memo/max_memo);
        if run_test==0 || ~ismember(get(handles.Data_listbox,'Value'),[1 2 3])
            max_samp_len=fix(m_res.num_measure_samples/num_init_runs); % in samples
            max_pi_len=fix(m_res.num_isi/num_init_runs);
        elseif get(handles.Data_listbox,'Value')==1
            max_samp_len=870; % in samples      % $$$$$$$$$$
            max_pi_len=15; % in isi             % $$$$$$$$$$
        elseif get(handles.Data_listbox,'Value')==2
            max_samp_len=2000000; % in samples  % $$$$$$$$$$
            max_pi_len=500; % in isi            % $$$$$$$$$$
        elseif get(handles.Data_listbox,'Value')==3
            max_samp_len=2000000; % in samples  % $$$$$$$$$$
            max_pi_len=100; % in isi            % $$$$$$$$$$
        end

        if m_para.memo_num_pi_measures>0
            r_para.num_pi_runs=ceil(m_res.num_isi/max_pi_len);
            if m_para.num_samp_measures>0 && r_para.num_samp_runs==0
                set(0,'DefaultUIControlFontSize',16);
                mbh=msgbox(sprintf('Dataset might be too large.\nPlease increase the value of the variable ''max_memo'' !!!'),'Warning','warn','modal');
                htxt = findobj(mbh,'Type','text');
                set(htxt,'FontSize',12,'FontWeight','bold')
                mb_pos=get(mbh,'Position');
                set(mbh,'Position',[mb_pos(1:2) mb_pos(3)*1.5 mb_pos(4)])
                uiwait(mbh);
                ret=1;
                return
            end
        end

        r_para.run_pico_ends=cumsum(fix([max_pi_len*ones(1,r_para.num_pi_runs-1) m_res.num_isi-max_pi_len*(r_para.num_pi_runs-1)]));
        r_para.run_pico_starts=[1 r_para.run_pico_ends(1:end-1)+1];
        r_para.run_pico_lengths=r_para.run_pico_ends-r_para.run_pico_starts+1;
        if m_para.num_pili_measures>0
            r_para.run_pili_ends=2*r_para.run_pico_ends;
            r_para.run_pili_starts=[1 r_para.run_pili_ends(1:end-1)+1];
            r_para.run_pili_lengths=r_para.run_pili_ends-r_para.run_pili_starts+1;
        end

        if m_para.num_samp_measures>0
            r_para.num_samp_runs=ceil(m_res.num_measure_samples/max_samp_len);
            if r_para.num_samp_runs==0
                set(0,'DefaultUIControlFontSize',16);
                mbh=msgbox(sprintf('Dataset might be too large.\nPlease increase the value of the variable ''max_memo'' !!!'),'Warning','warn','modal');
                htxt = findobj(mbh,'Type','text');
                set(htxt,'FontSize',12,'FontWeight','bold')
                mb_pos=get(mbh,'Position');
                set(mbh,'Position',[mb_pos(1:2) mb_pos(3)*1.5 mb_pos(4)])
                uiwait(mbh);
                ret=1;
                return
            end

            run_samp_ends_init=cumsum(fix([max_samp_len*ones(1,r_para.num_samp_runs-1) m_res.num_measure_samples-max_samp_len*(r_para.num_samp_runs-1)]));
            run_samp_starts_init=[1 run_samp_ends_init(1:end-1)+1];

            r_para.run_samp_ends=zeros(1,r_para.num_samp_runs);
            for samp_ruc=1:r_para.num_samp_runs
                r_para.run_samp_ends(samp_ruc)=round((m_res.cum_isi(r_para.run_pico_ends(samp_ruc)+1)-d_para.tmin)/d_para.dtm);
            end
            r_para.run_samp_starts=[1 r_para.run_samp_ends(1:end-1)+1];
            r_para.run_samp_lengths=r_para.run_samp_ends-r_para.run_samp_starts+1;
        end
    end
    if m_para.num_samp_measures>0 && r_para.num_samp_runs==1
        r_para.run_samp_lengths=m_res.num_measure_samples;
        r_para.run_samp_starts=1;
        r_para.run_samp_ends=m_res.num_measure_samples;
    end
    if r_para.num_pi_runs==1
        r_para.run_pico_lengths=m_res.num_isi;
        r_para.run_pico_starts=1;
        r_para.run_pico_ends=m_res.num_isi;
        if m_para.num_pili_measures>0
            r_para.run_pili_lengths=2*m_res.num_isi;
            r_para.run_pili_starts=1;
            r_para.run_pili_ends=2*m_res.num_isi;
        end
    end

    if m_para.num_samp_measures>0                                                         % TIME
        start=cell(1,d_para.num_trains);                                                % integers relative to data sampling
        ende=cell(1,d_para.num_trains);
        for trac=1:d_para.num_trains
            start{trac}=fix((uspikes{trac}(1:num_uspikes(trac)-1)+d_para.dts/10-d_para.tmin)/d_para.dts)+1;
            ende{trac}=fix((uspikes{trac}(2:num_uspikes(trac))+d_para.dts/10-d_para.tmin)/d_para.dts);
        end
    end

    ave_bi_vect=zeros(m_para.num_sel_bi_measures,d_para.num_pairs);
    if r_para.num_samp_runs>1 || r_para.num_pi_runs>1
        disp(' ');
        disp('Large data set. Please be patient.')
        disp(' ');
        if r_para.num_pi_runs>1
            disp(['Number of calculation loop runs: ',num2str(r_para.num_pi_runs)])
            pwbh = waitbar(0,'Large data set. Please be patient.');
        end
        if r_para.num_samp_runs>1
            disp(['Number of calculation loop runs: ',num2str(r_para.num_samp_runs)])
            swbh = waitbar(0,'Large data set. Please be patient.');
        end
    end


    % #####################################################################################################################################
    % ################################################################# Pi-Measures #######################################################
    % #####################################################################################################################################
    
    for pi_ruc=r_para.num_pi_runs:-1:1
        
        % ###########################################################################################################################################
        % ################################################################# Pico-Pili-Quantities #####################################################
        % ###########################################################################################################################################

        if any(f_para.subplot_posi([m_para.spike_pico m_para.realtime_spike_pico m_para.spike_pili m_para.realtime_spike_pili m_para.test_pili_1]))  % SPIKE-Pre-Pico
            if max_num_uspikes<256                                                % integers relative to data sampling
                previs_indy=zeros(d_para.num_trains,max_num_uspikes-1,'uint8');
                prev_spikes_indy=zeros(d_para.num_trains,num_all_isi,'uint8');
            elseif max_num_uspikes<65536
                previs_indy=zeros(d_para.num_trains,max_num_uspikes-1,'uint16');
                prev_spikes_indy=zeros(d_para.num_trains,num_all_isi,'uint16');
            elseif max_num_uspikes<2^32
                previs_indy=zeros(d_para.num_trains,max_num_uspikes-1,'uint32');
                prev_spikes_indy=zeros(d_para.num_trains,num_all_isi,'uint32');
            else
                previs_indy=zeros(d_para.num_trains,max_num_uspikes-1,'uint64');
                prev_spikes_indy=zeros(d_para.num_trains,num_all_isi,'uint64');
            end
            for trac=1:d_para.num_trains
                previs_indy(trac,1:num_uspikes(trac)-1)=1:num_uspikes(trac)-1;
                ivs=[1 find(all_trains==trac)];
                ive=[ivs(2:length(ivs))-1 num_all_isi];
                for ic=1:num_uspikes(trac)-1
                    prev_spikes_indy(trac,ivs(ic):ive(ic))=previs_indy(trac,ic);
                end
                if f_para.edge_correction && num_uspikes(trac)>3           % $$$$$$$$$$$$
                    prev_spikes_indy(trac,prev_spikes_indy(trac,:)==1)=2;
                end
            end
            prev_spikes_indy=prev_spikes_indy(:,all_isi>0);
            clear previs_indy

            previs=zeros(d_para.num_trains,max_num_uspikes-1,'single');
            prev_spikes=zeros(d_para.num_trains,num_all_isi,'single');
            for trac=1:d_para.num_trains
                previs(trac,1:num_uspikes(trac)-1)=uspikes{trac}(1:num_uspikes(trac)-1);
                ivs=[1 find(all_trains==trac)];
                ive=[ivs(2:length(ivs))-1 num_all_isi];
                for ic=1:num_uspikes(trac)-1
                    prev_spikes(trac,ivs(ic):ive(ic))=previs(trac,ic);
                end
            end
            prev_spikes=prev_spikes(:,all_isi>0);
            clear previs

            if any(f_para.subplot_posi([m_para.spike_pili m_para.realtime_spike_pili m_para.test_pili_1]))  % SPIKE-Pre-Pili ***************
                prev_spikes_indy_pili=zeros(d_para.num_trains,m_res.pili_len,'uint64');
                prev_spikes_pili=zeros(d_para.num_trains,m_res.pili_len,'single');
                for trac=1:d_para.num_trains
                    prev_spikes_indy_pili(trac,:)=reshape(repmat(prev_spikes_indy(trac,1:m_res.num_isi),2,1),1,m_res.pili_len);
                    if f_para.edge_correction && num_uspikes(trac)>3                                                                 % $$$$$$$$$$$$
                        prev_spikes_indy_pili(trac,prev_spikes_indy_pili(trac,:)==1)=2;
                    end
                    prev_spikes_pili(trac,:)=m_res.pili_supi-reshape([prev_spikes(trac,1:m_res.num_isi); prev_spikes(trac,1:m_res.num_isi)],1,m_res.pili_len);
                end
            end
        end

        if any(f_para.subplot_posi([m_para.spike_pico m_para.future_spike_pico m_para.spike_pili m_para.future_spike_pili m_para.test_pili_1]))  % SPIKE-Future-Pico
            if max_num_uspikes<256                                                % integers relative to data sampling
                follis_indy=zeros(d_para.num_trains,max_num_uspikes-1,'uint8');
                foll_spikes_indy=zeros(d_para.num_trains,num_all_isi,'uint8');
            elseif max_num_uspikes<65536
                follis_indy=zeros(d_para.num_trains,max_num_uspikes-1,'uint16');
                foll_spikes_indy=zeros(d_para.num_trains,num_all_isi,'uint16');
            elseif max_num_uspikes<2^32
                follis_indy=zeros(d_para.num_trains,max_num_uspikes-1,'uint32');
                foll_spikes_indy=zeros(d_para.num_trains,num_all_isi,'uint32');
            else
                follis_indy=zeros(d_para.num_trains,max_num_uspikes-1,'uint64');
                foll_spikes_indy=zeros(d_para.num_trains,num_all_isi,'uint64');
            end
            for trac=1:d_para.num_trains
                follis_indy(trac,1:num_uspikes(trac)-1)=2:num_uspikes(trac);
                ivs=[1 find(all_trains==trac)];
                ive=[ivs(2:length(ivs))-1 num_all_isi];
                for ic=1:num_uspikes(trac)-1
                    foll_spikes_indy(trac,ivs(ic):ive(ic))=follis_indy(trac,ic);
                end
                if f_para.edge_correction && num_uspikes(trac)>3    % $$$$$$$$$$$$
                    foll_spikes_indy(trac,foll_spikes_indy(trac,:)==num_uspikes(trac))=num_uspikes(trac)-1;
                end
            end
            foll_spikes_indy=foll_spikes_indy(:,all_isi>0);
            clear follis_indy

            follis=zeros(d_para.num_trains,max_num_uspikes-1,'single');
            foll_spikes=zeros(d_para.num_trains,num_all_isi,'single');
            for trac=1:d_para.num_trains
                follis(trac,1:num_uspikes(trac)-1)=uspikes{trac}(2:num_uspikes(trac));
                ivs=[1 find(all_trains==trac)];
                ive=[ivs(2:length(ivs))-1 num_all_isi];
                for ic=1:num_uspikes(trac)-1
                    foll_spikes(trac,ivs(ic):ive(ic))=follis(trac,ic);
                end
            end
            foll_spikes=foll_spikes(:,all_isi>0);
            clear follis

            if any(f_para.subplot_posi([m_para.spike_pili m_para.future_spike_pili m_para.test_pili_1]))   % SPIKE-Future-Pili **************
                foll_spikes_indy_pili=zeros(d_para.num_trains,m_res.pili_len,'uint64');
                foll_spikes_pili=zeros(d_para.num_trains,m_res.pili_len,'single');
                for trac=1:d_para.num_trains
                    foll_spikes_indy_pili(trac,:)=reshape(repmat(foll_spikes_indy(trac,1:m_res.num_isi),2,1),1,m_res.pili_len);
                    if f_para.edge_correction && num_uspikes(trac)>3                                                                 % $$$$$$$$$$$$
                        foll_spikes_indy_pili(trac,foll_spikes_indy_pili(trac,:)==num_uspikes(trac))=num_uspikes(trac)-1;
                    end
                    foll_spikes_pili(trac,:)=reshape([foll_spikes(trac,1:m_res.num_isi); foll_spikes(trac,1:m_res.num_isi)],1,m_res.pili_len)-m_res.pili_supi;
                end
            end
        end
        
        % #####################################################################################################################################
        % ################################################################# Pili-Measures #####################################################
        % #####################################################################################################################################

        localuspikes=cell(1, d_para.num_trains);
        
        temp_prev=zeros(d_para.num_trains, r_para.run_pili_lengths(pi_ruc),'int32');
        temp_foll=zeros(d_para.num_trains, r_para.run_pili_lengths(pi_ruc),'int32');
        index=cell(1,d_para.num_trains);
        for trac=1:d_para.num_trains
            index{trac}=[prev_spikes_indy_pili(trac,r_para.run_pili_starts(pi_ruc)) foll_spikes_indy_pili(trac,r_para.run_pili_ends(pi_ruc))];
            temp_prev(trac,:)=int32(prev_spikes_indy_pili(trac,r_para.run_pili_starts(pi_ruc):r_para.run_pili_ends(pi_ruc)))-...
                int32(prev_spikes_indy_pili(trac,r_para.run_pili_starts(pi_ruc)))+1;
            temp_foll(trac,:)=int32(foll_spikes_indy_pili(trac,r_para.run_pili_starts(pi_ruc):r_para.run_pili_ends(pi_ruc)))-...
                int32(foll_spikes_indy_pili(trac,r_para.run_pili_starts(pi_ruc)))+1;
            localuspikes{trac}=uspikes{trac}(int32(index{trac}(1)):int32(index{trac}(2)));
        end
        local_num_uspikes=cellfun('length',localuspikes);
        if exist('SPIKY_udists_MEX.mexw32','file')
            udists=SPIKY_udists_MEX(int32(d_para.num_trains),int32(local_num_uspikes),localuspikes);
        else
            udists=cell(d_para.num_trains);
            for trac1=1:d_para.num_trains
                if d_para.num_trains>=100 && max(local_num_uspikes)>1000
                    disp(['udistc = ',num2str([1 trac1])])
                end
                for trac2=setdiff(1:d_para.num_trains,trac1)
                    udists{trac1,trac2}=zeros(1,local_num_uspikes(trac1),'single');
                end
            end
            for trac1=1:d_para.num_trains
                if d_para.num_trains>=100 && max(local_num_uspikes)>1000
                    disp(['udistc = ',num2str([2 trac1])])
                end
                for trac2=setdiff(1:d_para.num_trains,trac1)
                    for spc=1:local_num_uspikes(trac1)
                        udists{trac1,trac2}(spc)=min(abs(localuspikes{trac1}(spc)-localuspikes{trac2}(1:local_num_uspikes(trac2))));
                    end
                end
            end
        end
        
        for trac=1:d_para.num_trains
            for tric=1:d_para.num_trains
                if (trac~=tric)
                    udists{trac,tric}(1)=min(abs(uspikes{trac}(index{trac}(1))-uspikes{tric}));
                    udists{trac,tric}(end)=min(abs(uspikes{trac}(index{trac}(2))-uspikes{tric}));
                end
            end
        end        
        
        if m_para.num_pili_measures>0
            odds=1:2:r_para.run_pili_lengths(pi_ruc);
            evens=odds+1;
            m_res.pili_measures_mat=zeros(m_para.num_pili_measures,d_para.num_pairs,r_para.run_pili_lengths(pi_ruc),'single');
            if f_para.subplot_posi(m_para.spike_pili)                                                   % SPIKE-Pili
                if exist('SPIKY_SPIKE_MEX.mexw32','file')
                    m_res.pili_measures_mat(m_para.measure_indy(m_para.spike_pili),1:d_para.num_pairs,1:r_para.run_pili_lengths(pi_ruc)) = ...
                        SPIKY_SPIKE_MEX(int32(d_para.num_pairs),int32(r_para.run_pili_lengths(pi_ruc)),int32(d_para.num_trains),...
                        foll_spikes_pili(:,r_para.run_pili_starts(pi_ruc):r_para.run_pili_ends(pi_ruc)),...
                        prev_spikes_pili(:,r_para.run_pili_starts(pi_ruc):r_para.run_pili_ends(pi_ruc)),...
                        int32(isi_indy_pili(:,r_para.run_pili_starts(pi_ruc):r_para.run_pili_ends(pi_ruc))),...
                        temp_prev,temp_foll,...
                        int32(r_para.run_pico_starts(pi_ruc)),ints(:,r_para.run_pico_starts(pi_ruc):r_para.run_pico_ends(pi_ruc)),...
                        udists);
%                     new = squeeze(m_res.pili_measures_mat(m_para.measure_indy(m_para.spike_pili),1:d_para.num_pairs,1:r_para.run_pili_lengths(pi_ruc)));
%                     old = ...
%                         SPIKY_SPIKE_MEX(int32(d_para.num_pairs),int32(r_para.run_pili_lengths(pi_ruc)),int32(d_para.num_trains),...
%                         foll_spikes_pili(:,r_para.run_pili_starts(pi_ruc):r_para.run_pili_ends(pi_ruc)),...
%                         prev_spikes_pili(:,r_para.run_pili_starts(pi_ruc):r_para.run_pili_ends(pi_ruc)),...
%                         int32(isi_indy_pili(:,r_para.run_pili_starts(pi_ruc):r_para.run_pili_ends(pi_ruc))),...
%                         int32(prev_spikes_indy_pili(:,r_para.run_pili_starts(pi_ruc):r_para.run_pili_ends(pi_ruc))),...
%                         int32(foll_spikes_indy_pili(:,r_para.run_pili_starts(pi_ruc):r_para.run_pili_ends(pi_ruc))),...
%                         int32(r_para.run_pico_starts(pi_ruc)),ints(:,r_para.run_pico_starts(pi_ruc):r_para.run_pico_ends(pi_ruc)),...
%                         udists);
                else
                    run_pico_range=r_para.run_pico_starts(pi_ruc):r_para.run_pico_ends(pi_ruc);
                    pac=0;
                    for trac1=1:d_para.num_trains-1
                        for trac2=trac1+1:d_para.num_trains
                            pac=pac+1;
                            dummy=double(isi_indy_pili(r_para.run_pili_starts(pi_ruc):r_para.run_pili_ends(pi_ruc)))-r_para.run_pico_starts(pi_ruc)+1;
                            m_res.pili_measures_mat(m_para.measure_indy(m_para.spike_pili),pac,1:r_para.run_pili_lengths(pi_ruc)) = ...
                                ((udists{trac1,trac2}(temp_prev(trac1,:)).*...
                                foll_spikes_pili(trac1,r_para.run_pili_starts(pi_ruc):r_para.run_pili_ends(pi_ruc))+...
                                udists{trac1,trac2}(temp_foll(trac1,:)).*...
                                prev_spikes_pili(trac1,r_para.run_pili_starts(pi_ruc):r_para.run_pili_ends(pi_ruc)))./...
                                ints(trac1,isi_indy_pili(r_para.run_pili_starts(pi_ruc):r_para.run_pili_ends(pi_ruc))).*...
                                ints(trac2,isi_indy_pili(r_para.run_pili_starts(pi_ruc):r_para.run_pili_ends(pi_ruc)))+...
                                (udists{trac2,trac1}(temp_prev(trac2,:)).*...
                                foll_spikes_pili(trac2,r_para.run_pili_starts(pi_ruc):r_para.run_pili_ends(pi_ruc))+...
                                udists{trac2,trac1}(temp_foll(trac2,:)).*...
                                prev_spikes_pili(trac2,r_para.run_pili_starts(pi_ruc):r_para.run_pili_ends(pi_ruc)))./...
                                ints(trac2,isi_indy_pili(r_para.run_pili_starts(pi_ruc):r_para.run_pili_ends(pi_ruc))).*...
                                ints(trac1,isi_indy_pili(r_para.run_pili_starts(pi_ruc):r_para.run_pili_ends(pi_ruc))))./...
                                ((ints(trac1,run_pico_range(dummy))+ints(trac2,run_pico_range(dummy))).^2/2);
%                             m_res.pili_measures_mat(m_para.measure_indy(m_para.spike_pili),pac,1:r_para.run_pili_lengths(pi_ruc)) = ...
%                                 ((udists{trac1,trac2}(prev_spikes_indy_pili(trac1,r_para.run_pili_starts(pi_ruc):r_para.run_pili_ends(pi_ruc))).*...
%                                 foll_spikes_pili(trac1,r_para.run_pili_starts(pi_ruc):r_para.run_pili_ends(pi_ruc))+...
%                                 udists{trac1,trac2}(foll_spikes_indy_pili(trac1,r_para.run_pili_starts(pi_ruc):r_para.run_pili_ends(pi_ruc))).*...
%                                 prev_spikes_pili(trac1,r_para.run_pili_starts(pi_ruc):r_para.run_pili_ends(pi_ruc)))./...
%                                 ints(trac1,isi_indy_pili(r_para.run_pili_starts(pi_ruc):r_para.run_pili_ends(pi_ruc))).*...
%                                 ints(trac2,isi_indy_pili(r_para.run_pili_starts(pi_ruc):r_para.run_pili_ends(pi_ruc)))+...
%                                 (udists{trac2,trac1}(prev_spikes_indy_pili(trac2,r_para.run_pili_starts(pi_ruc):r_para.run_pili_ends(pi_ruc))).*...
%                                 foll_spikes_pili(trac2,r_para.run_pili_starts(pi_ruc):r_para.run_pili_ends(pi_ruc))+...
%                                 udists{trac2,trac1}(foll_spikes_indy_pili(trac2,r_para.run_pili_starts(pi_ruc):r_para.run_pili_ends(pi_ruc))).*...
%                                 prev_spikes_pili(trac2,r_para.run_pili_starts(pi_ruc):r_para.run_pili_ends(pi_ruc)))./...
%                                 ints(trac2,isi_indy_pili(r_para.run_pili_starts(pi_ruc):r_para.run_pili_ends(pi_ruc))).*...
%                                 ints(trac1,isi_indy_pili(r_para.run_pili_starts(pi_ruc):r_para.run_pili_ends(pi_ruc))))./...
%                                 ((ints(trac1,run_pico_range(dummy))+ints(trac2,run_pico_range(dummy))).^2/2);
                        end
                    end
                end
                ave_bi_vect(logical(m_para.select_bi_measures==m_para.spike_pili),:)=...
                    sum((m_res.pili_measures_mat(m_para.measure_indy(m_para.spike_pili),:,odds)+...
                    m_res.pili_measures_mat(m_para.measure_indy(m_para.spike_pili),:,evens))/2.*...
                    repmat(shiftdim(m_res.isi(r_para.run_pico_starts(pi_ruc):r_para.run_pico_ends(pi_ruc)),-1),[1,d_para.num_pairs]),3)/...
                    sum(m_res.isi(r_para.run_pico_starts(pi_ruc):r_para.run_pico_ends(pi_ruc)));
            end

            if f_para.subplot_posi(m_para.realtime_spike_pili)                                           % REALTIME-Pili
                if exist('SPIKY_realtimeSPIKE_MEX.mexw32','file')
                    m_res.pili_measures_mat(m_para.measure_indy(m_para.realtime_spike_pili),1:d_para.num_pairs,1:r_para.run_pili_lengths(pi_ruc)) = ...
                        SPIKY_realtimeSPIKE_MEX(int32(d_para.num_pairs),int32(r_para.run_pili_lengths(pi_ruc)),int32(d_para.num_trains),...
                        prev_spikes_pili(:,r_para.run_pili_starts(pi_ruc):r_para.run_pili_ends(pi_ruc)),...
                        temp_prev,udists);
%                     m_res.pili_measures_mat(m_para.measure_indy(m_para.realtime_spike_pili),1:d_para.num_pairs,1:r_para.run_pili_lengths(pi_ruc)) = ...
%                         SPIKY_realtimeSPIKE_MEX(int32(d_para.num_pairs),int32(r_para.run_pili_lengths(pi_ruc)),int32(d_para.num_trains),...
%                         prev_spikes_pili(:,r_para.run_pili_starts(pi_ruc):r_para.run_pili_ends(pi_ruc)),...
%                         int32(prev_spikes_indy_pili(:,r_para.run_pili_starts(pi_ruc):r_para.run_pili_ends(pi_ruc))),udists);
                else
                    run_pili_range=r_para.run_pili_starts(pi_ruc):r_para.run_pili_ends(pi_ruc);
                    pac=0;
                    for trac1=1:d_para.num_trains-1
                        for trac2=trac1+1:d_para.num_trains
                            pac=pac+1;
                            normy=(prev_spikes_pili(trac1,run_pili_range)+prev_spikes_pili(trac2,run_pili_range));
                            dummy=(prev_spikes_pili(trac1,run_pili_range)<prev_spikes_pili(trac2,run_pili_range)-0.00000001);
                            dummy1=trac1*(1-dummy)+trac2*dummy;   % index of spike train with earlier spike
                            dummy2=trac2*(1-dummy)+trac1*dummy;   % index of spike train with later spike
                            for sc=1:r_para.run_pili_lengths(pi_ruc)
                                m_res.pili_measures_mat(m_para.measure_indy(m_para.realtime_spike_pili),pac,sc)= ...
                                    (abs(prev_spikes_pili(trac1,run_pili_range(sc))-...
                                    prev_spikes_pili(trac2,run_pili_range(sc)))+...   % later spike
                                    udists{dummy1(sc),dummy2(sc)}(temp_prev(dummy1(sc),sc)))...      % earlier spike
                                    /(2*normy(sc)+(normy(sc)==0));
%                                 m_res.pili_measures_mat(m_para.measure_indy(m_para.realtime_spike_pili),pac,sc)= ...
%                                     (abs(prev_spikes_pili(trac1,run_pili_range(sc))-...
%                                     prev_spikes_pili(trac2,run_pili_range(sc)))+...   % later spike
%                                     udists{dummy1(sc),dummy2(sc)}(prev_spikes_indy_pili(dummy1(sc),run_pili_range(sc))))...      % earlier spike
%                                     /(2*normy(sc)+(normy(sc)==0));
                            end
                        end
                    end
                    clear dummy; clear dummy1; clear dummy2
                end
                aves=(log(1./m_res.pili_measures_mat(m_para.measure_indy(m_para.realtime_spike_pili),:,evens))-...
                    log(1./m_res.pili_measures_mat(m_para.measure_indy(m_para.realtime_spike_pili),:,odds)))./...
                    (1./m_res.pili_measures_mat(m_para.measure_indy(m_para.realtime_spike_pili),:,evens)-...
                    1./m_res.pili_measures_mat(m_para.measure_indy(m_para.realtime_spike_pili),:,odds));
                aves(isnan(aves))=0;
                ave_bi_vect(logical(m_para.select_bi_measures==m_para.realtime_spike_pili),:)=...
                    sum(aves.*repmat(shiftdim(m_res.isi(r_para.run_pico_starts(pi_ruc):r_para.run_pico_ends(pi_ruc)),-1),...
                    [1,d_para.num_pairs]),3)/sum(m_res.isi(r_para.run_pico_starts(pi_ruc):r_para.run_pico_ends(pi_ruc)));
            end

            if f_para.subplot_posi(m_para.future_spike_pili)                                             % FUTURE-Pili
                if exist('SPIKY_futureSPIKE_MEX.mexw32','file')
                    m_res.pili_measures_mat(m_para.measure_indy(m_para.future_spike_pili),1:d_para.num_pairs,1:r_para.run_pili_lengths(pi_ruc)) = ...
                        SPIKY_futureSPIKE_MEX(int32(d_para.num_pairs),int32(r_para.run_pili_lengths(pi_ruc)),int32(d_para.num_trains),...
                        foll_spikes_pili(:,r_para.run_pili_starts(pi_ruc):r_para.run_pili_ends(pi_ruc)),...
                        temp_foll,udists);
%                     m_res.pili_measures_mat(m_para.measure_indy(m_para.future_spike_pili),1:d_para.num_pairs,1:r_para.run_pili_lengths(pi_ruc)) = ...
%                         SPIKY_futureSPIKE_MEX(int32(d_para.num_pairs),int32(r_para.run_pili_lengths(pi_ruc)),int32(d_para.num_trains),...
%                         foll_spikes_pili(:,r_para.run_pili_starts(pi_ruc):r_para.run_pili_ends(pi_ruc)),...
%                         int32(foll_spikes_pili_indy(:,r_para.run_pili_starts(pi_ruc):r_para.run_pili_ends(pi_ruc))),udists);
                else
                    run_pili_range=r_para.run_pili_starts(pi_ruc):r_para.run_pili_ends(pi_ruc);
                    pac=0;
                    for trac1=1:d_para.num_trains-1
                        for trac2=trac1+1:d_para.num_trains
                            pac=pac+1;
                            normy=(foll_spikes_pili(trac1,run_pili_range)+foll_spikes_pili(trac2,run_pili_range));
                            dummy=(foll_spikes_pili(trac1,run_pili_range)<foll_spikes_pili(trac2,run_pili_range)-0.00000001);
                            dummy1=trac1*(1-dummy)+trac2*dummy;   % index of spike train with earlier spike
                            dummy2=trac2*(1-dummy)+trac1*dummy;   % index of spike train with later spike
                            for sc=1:r_para.run_pili_lengths(pi_ruc)
                                m_res.pili_measures_mat(m_para.measure_indy(m_para.future_spike_pili),pac,sc)= ...
                                    (abs(foll_spikes_pili(trac1,run_pili_range(sc))-...
                                    foll_spikes_pili(trac2,run_pili_range(sc)))+...   % later spike
                                    udists{dummy1(sc),dummy2(sc)}(temp_foll(dummy1(sc),sc)))...      % earlier spike
                                    /(2*normy(sc)+(normy(sc)==0));
%                                 m_res.pili_measures_mat(m_para.measure_indy(m_para.future_spike_pili),pac,sc)= ...
%                                     (abs(foll_spikes_pili(trac1,run_pili_range(sc))-...
%                                     foll_spikes_pili(trac2,run_pili_range(sc)))+...   % later spike
%                                     udists{dummy1(sc),dummy2(sc)}(foll_spikes_indy_pili(dummy1(sc),run_pili_range(sc))))...      % earlier spike
%                                     /(2*normy(sc)+(normy(sc)==0));
                            end
                        end
                    end
                    clear dummy; clear dummy1; clear dummy2
                end
                aves=(log(1./m_res.pili_measures_mat(m_para.measure_indy(m_para.future_spike_pili),:,evens))-...
                    log(1./m_res.pili_measures_mat(m_para.measure_indy(m_para.future_spike_pili),:,odds)))./...
                    (1./m_res.pili_measures_mat(m_para.measure_indy(m_para.future_spike_pili),:,evens)-...
                    1./m_res.pili_measures_mat(m_para.measure_indy(m_para.future_spike_pili),:,odds));
                aves(isnan(aves))=0;
                ave_bi_vect(logical(m_para.select_bi_measures==m_para.future_spike_pili),:)=...
                    sum(aves.*repmat(shiftdim(m_res.isi(r_para.run_pico_starts(pi_ruc):r_para.run_pico_ends(pi_ruc)),-1),...
                    [1,d_para.num_pairs]),3)/sum(m_res.isi(r_para.run_pico_starts(pi_ruc):r_para.run_pico_ends(pi_ruc)));
            end

            if f_para.subplot_posi(m_para.test_pili_1)   % $$$$$$                                                 % Test-1
                m_res.pili_measures_mat(m_para.measure_indy(m_para.test_pili_1),1:d_para.num_pairs,1:r_para.run_pili_lengths(pi_ruc)) = ...
                    SPIKY_SPIKE_MEX(int32(d_para.num_pairs),int32(r_para.run_pili_lengths(pi_ruc)),int32(d_para.num_trains),...
                    foll_spikes_pili(:,r_para.run_pili_starts(pi_ruc):r_para.run_pili_ends(pi_ruc)),...
                    prev_spikes_pili(:,r_para.run_pili_starts(pi_ruc):r_para.run_pili_ends(pi_ruc)),...
                    int32(isi_indy_pili(:,r_para.run_pili_starts(pi_ruc):r_para.run_pili_ends(pi_ruc))),...
                    int32(prev_spikes_indy_pili2(:,r_para.run_pili_starts(pi_ruc):r_para.run_pili_ends(pi_ruc))),...
                    int32(foll_spikes_indy_pili2(:,r_para.run_pili_starts(pi_ruc):r_para.run_pili_ends(pi_ruc))),...
                    int32(r_para.run_pico_starts(pi_ruc)),ints2(:,r_para.run_pico_starts(pi_ruc):r_para.run_pico_ends(pi_ruc)),udists);
                ave_bi_vect(logical(m_para.select_bi_measures==m_para.test_pili_1),:)=...
                    sum((m_res.pili_measures_mat(m_para.measure_indy(m_para.test_pili_1),:,odds)+...
                    m_res.pili_measures_mat(m_para.measure_indy(m_para.test_pili_1),:,evens))/2.*...
                    repmat(shiftdim(m_res.isi,-1),[1,d_para.num_pairs]),3)/sum(m_res.isi);
            end

            if f_para.subplot_posi(m_para.test_pili_2)                                                    % Test-2
                b=zeros(d_para.num_pairs,m_res.num_isi*2);
                b(:,1:2:end)=pair_spike_status;
                b(:,2:2:end)=pair_spike_status;
                m_res.pili_measures_mat(m_para.measure_indy(m_para.test_pili_2),1:d_para.num_pairs,1:r_para.run_pili_lengths(pi_ruc)) = ...
                    (b==0).*shiftdim(m_res.pili_measures_mat(1,:,:),1)+...
                    (b==2).*shiftdim(m_res.pili_measures_mat(2,:,:),1)+...
                    (b==6).*shiftdim(m_res.pili_measures_mat(3,:,:),1);
                ave_bi_vect(logical(m_para.select_bi_measures==m_para.test_pili_2),:)=...
                    sum((m_res.pili_measures_mat(m_para.measure_indy(m_para.test_pili_2),:,odds)+...
                    m_res.pili_measures_mat(m_para.measure_indy(m_para.test_pili_2),:,evens))/2.*...
                    repmat(shiftdim(m_res.isi,-1),[1,d_para.num_pairs]),3)/sum(m_res.isi);
            end

            if f_para.subplot_posi(m_para.test_pili_3)                                                    % Test-3
                b=zeros(d_para.num_pairs,m_res.num_isi*2);
                b(:,1:2:end)=pair_spike_status;
                b(:,2:2:end)=pair_spike_status;
                m_res.pili_measures_mat(m_para.measure_indy(m_para.test_pili_3),1:d_para.num_pairs,1:r_para.run_pili_lengths(pi_ruc)) = ...
                    (b==0).*shiftdim(m_res.pili_measures_mat(1,:,:),1)+...
                    (b==2).*shiftdim(m_res.pili_measures_mat(2,:,:),1)+...
                    (b==6).*shiftdim(m_res.pili_measures_mat(3,:,:),1);
                ave_bi_vect(logical(m_para.select_bi_measures==m_para.test_pili_3),:)=...
                    sum((m_res.pili_measures_mat(m_para.measure_indy(m_para.test_pili_3),:,odds)+...
                    m_res.pili_measures_mat(m_para.measure_indy(m_para.test_pili_3),:,evens))/2.*...
                    repmat(shiftdim(m_res.isi,-1),[1,d_para.num_pairs]),3)/sum(m_res.isi);
            end
        end
        
        % #####################################################################################################################################
        % ################################################################# Pico-Measures #####################################################
        % #####################################################################################################################################

        if m_para.num_pico_measures>0                                                         % Pico
            m_res.pico_measures_mat=zeros(m_para.num_pico_measures,d_para.num_pairs,r_para.run_pico_lengths(pi_ruc),'single');

            if f_para.subplot_posi(m_para.isi_pico)                                                 % ISI
                if exist('SPIKY_ISI_MEX.mexw32','file')
                    m_res.pico_measures_mat(m_para.measure_indy(m_para.isi_pico),1:d_para.num_pairs,1:r_para.run_pico_lengths(pi_ruc)) = ...
                        abs(SPIKY_ISI_MEX(int32(d_para.num_pairs),int32(r_para.run_pico_lengths(pi_ruc)),int32(d_para.num_trains),...
                        ints(:,r_para.run_pico_starts(pi_ruc):r_para.run_pico_ends(pi_ruc))));

                else
                    pac=0;
                    run_pico_range=r_para.run_pico_starts(pi_ruc):r_para.run_pico_ends(pi_ruc);
                    for trac1=1:d_para.num_trains-1
                        for trac2=trac1+1:d_para.num_trains
                            pac=pac+1;
                            dummy1=find(ints(trac1,run_pico_range)<ints(trac2,run_pico_range));
                            m_res.pico_measures_mat(m_para.measure_indy(m_para.isi_pico),pac,dummy1)=abs(ints(trac1,run_pico_range(dummy1))./ints(trac2,run_pico_range(dummy1))-1);
                            dummy2=find(ints(trac1,run_pico_range)>=ints(trac2,run_pico_range) & ints(trac1,run_pico_range)~=0);
                            m_res.pico_measures_mat(m_para.measure_indy(m_para.isi_pico),pac,dummy2)=abs(ints(trac2,run_pico_range(dummy2))./ints(trac1,run_pico_range(dummy2))-1);
                        end
                    end
                end
            end

            if f_para.subplot_posi(m_para.spike_pico)                                         % SPIKE-Pico
                m_res.pico_measures_mat(m_para.measure_indy(m_para.spike_pico),1:d_para.num_pairs,1:r_para.run_pico_lengths(pi_ruc)) = ...
                    SPIKY_SPIKEpico_MEX(int32(d_para.num_pairs),int32(r_para.run_pico_lengths(pi_ruc)),int32(d_para.num_trains),...
                    foll_spikes(:,r_para.run_pico_starts(pi_ruc):r_para.run_pico_ends(pi_ruc)),m_res.isi_pos(:,r_para.run_pico_starts(pi_ruc):...
                    r_para.run_pico_ends(pi_ruc)),int32(prev_spikes_indy(:,r_para.run_pico_starts(pi_ruc):r_para.run_pico_ends(pi_ruc))),...
                    ints(:,r_para.run_pico_starts(pi_ruc):r_para.run_pico_ends(pi_ruc)),int32(foll_spikes_indy(:,r_para.run_pico_starts(pi_ruc):...
                    r_para.run_pico_ends(pi_ruc))),prev_spikes(:,r_para.run_pico_starts(pi_ruc):r_para.run_pico_ends(pi_ruc)),udists);
            end

            if f_para.subplot_posi(m_para.realtime_spike_pico)                                  % REALTIME-PICO
                m_res.pico_measures_mat(m_para.measure_indy(m_para.realtime_spike_pico),1:d_para.num_pairs,1:r_para.run_pico_lengths(pi_ruc)) = ...
                    SPIKY_realtimeSPIKEpico_MEX(int32(d_para.num_pairs),int32(r_para.run_pico_lengths(pi_ruc)),int32(d_para.num_trains),...
                    prev_spikes(:,r_para.run_pico_starts(pi_ruc):r_para.run_pico_ends(pi_ruc)),m_res.cum_isi(:,[r_para.run_pico_starts(pi_ruc):...
                    r_para.run_pico_ends(pi_ruc) r_para.run_pico_ends(pi_ruc)+1]),m_res.isi(:,r_para.run_pico_starts(pi_ruc):r_para.run_pico_ends(pi_ruc)),...
                    int32(prev_spikes_indy(:,r_para.run_pico_starts(pi_ruc):r_para.run_pico_ends(pi_ruc))),udists);
            end

            if f_para.subplot_posi(m_para.future_spike_pico)                                    % FUTURE-PICO
                m_res.pico_measures_mat(m_para.measure_indy(m_para.future_spike_pico),1:d_para.num_pairs,1:r_para.run_pico_lengths(pi_ruc)) = ...
                    SPIKY_futureSPIKEpico_MEX(int32(d_para.num_pairs),int32(r_para.run_pico_lengths(pi_ruc)),int32(d_para.num_trains),...
                    foll_spikes(:,r_para.run_pico_starts(pi_ruc):r_para.run_pico_ends(pi_ruc)),m_res.cum_isi(:,[r_para.run_pico_starts(pi_ruc):...
                    r_para.run_pico_ends(pi_ruc) r_para.run_pico_ends(pi_ruc)+1]),m_res.isi(:,r_para.run_pico_starts(pi_ruc):r_para.run_pico_ends(pi_ruc)),...
                    int32(foll_spikes_indy(:,r_para.run_pico_starts(pi_ruc):r_para.run_pico_ends(pi_ruc))),udists);
            end

            ave_bi_vect(m_para.pico_bi_measures_indy,:)=ave_bi_vect(m_para.pico_bi_measures_indy,:)+sum(m_res.pico_measures_mat.*...
                repmat(shiftdim(m_res.isi(r_para.run_pico_starts(pi_ruc):r_para.run_pico_ends(pi_ruc)),-1),[m_para.num_pico_measures,...
                d_para.num_pairs]),3);
        end        
        
        if r_para.num_pi_runs>1
            disp(['Calculation-Loop-Info: = ',num2str(r_para.num_pi_runs+1-pi_ruc),'  (',num2str(r_para.num_pi_runs),')'])
            save (['SPIKY_pi_AZBYCX_',num2str(pi_ruc)],'-struct','m_res','pi*_measures_mat')
            waitbar((r_para.num_pi_runs+1-pi_ruc)/r_para.num_pi_runs,pwbh,['Calculation-Loop-Info: ',num2str(r_para.num_pi_runs+1-pi_ruc),'  (',num2str(r_para.num_pi_runs),')'])
            if pi_ruc==1
                delete(pwbh)
            end
        end
    end
    r_para.pi_ruc=pi_ruc;
    
    % #####################################################################################################################################
    % ################################################################# Samp-Measures #####################################################
    % #####################################################################################################################################

    if m_para.num_samp_measures>0                                                         % TIME

        for samp_ruc=r_para.num_samp_runs:-1:1

            if any(f_para.subplot_posi([m_para.realtime_spike_samp m_para.spike_samp]))
                if r_para.num_samp_runs==1
                    previ=zeros(d_para.num_trains,m_res.num_samples,'single');                          % distance to previous spike
                    if max_num_uspikes<256
                        prev_indy=zeros(d_para.num_trains,m_res.num_samples,'uint8');
                    elseif max_num_uspikes<65536
                        prev_indy=zeros(d_para.num_trains,m_res.num_samples,'uint16');
                    elseif max_num_uspikes<2^32
                        prev_indy=zeros(d_para.num_trains,m_res.num_samples,'uint32');
                    else
                        prev_indy=zeros(d_para.num_trains,m_res.num_samples,'uint64');
                    end
                    prev_indy2=prev_indy;
                    for trac=1:d_para.num_trains
                        for sc=1:num_uspikes(trac)-1
                            previ(trac,start{trac}(sc):ende{trac}(sc))=(1:single(ende{trac}(sc)-start{trac}(sc)+1))*d_para.dts;
                            prev_indy(trac,start{trac}(sc):ende{trac}(sc))=sc;                         % index of preceding spike
                        end
                        prev_indy2(trac,1:m_res.num_samples)=prev_indy(trac,1:m_res.num_samples);
                        if num_uspikes(trac)>3                                                         % $$$$$$$$$$$$
                            prev_indy2(trac,prev_indy2(trac,:)==1)=2;
                        end
                    end
                    previ=previ(:,d_para.dsf:d_para.dsf:m_res.num_samples);
                    prev_indy=prev_indy(:,d_para.dsf:d_para.dsf:end);
                    prev_indy2=prev_indy2(:,d_para.dsf:d_para.dsf:end);
                else
                    previ=zeros(d_para.num_trains,r_para.run_samp_lengths(samp_ruc),'single');            % distance to previous spike
                    if max_num_uspikes<256
                        prev_indy=zeros(d_para.num_trains,r_para.run_samp_lengths(samp_ruc),'uint8');
                    elseif max_num_uspikes<65536
                        prev_indy=zeros(d_para.num_trains,r_para.run_samp_lengths(samp_ruc),'uint16');
                    elseif max_num_uspikes<2^32
                        prev_indy=zeros(d_para.num_trains,r_para.run_samp_lengths(samp_ruc),'uint32');
                    else
                        prev_indy=zeros(d_para.num_trains,r_para.run_samp_lengths(samp_ruc),'uint64');
                    end
                    for trac=1:d_para.num_trains
                        first=find(start{trac}(1:num_uspikes(trac)-1)<=r_para.run_samp_starts(samp_ruc)*d_para.dsf,1,'last');
                        last=find(ende{trac}(1:num_uspikes(trac)-1)>=r_para.run_samp_ends(samp_ruc)*d_para.dsf,1,'first');
                        if first<last
                            vect1=(r_para.run_samp_starts(samp_ruc):fix(ende{trac}(first)/d_para.dsf));
                            previ(trac,1:length(vect1))=(single(vect1*d_para.dsf)-single(start{trac}(first)))*d_para.dts;
                            prev_indy(trac,1:length(vect1))=first;
                            for sc=first+1:last-1
                                vect1=(ceil(start{trac}(sc)/d_para.dsf):fix(ende{trac}(sc)/d_para.dsf));
                                previ(trac,vect1-r_para.run_samp_starts(samp_ruc)+1)=(single(vect1*d_para.dsf)-single(start{trac}(sc)))*d_para.dts;
                                prev_indy(trac,vect1-r_para.run_samp_starts(samp_ruc)+1)=sc;
                            end
                            vect1=(ceil(start{trac}(last)/d_para.dsf):r_para.run_samp_ends(samp_ruc));
                            previ(trac,vect1-r_para.run_samp_starts(samp_ruc)+1)=(single(vect1*d_para.dsf)-single(start{trac}(last)))*d_para.dts;
                            prev_indy(trac,vect1-r_para.run_samp_starts(samp_ruc)+1)=last;
                        else
                            previ(trac,1:r_para.run_samp_lengths(samp_ruc))=(single((r_para.run_samp_starts(samp_ruc):...
                                r_para.run_samp_ends(samp_ruc))*d_para.dsf)-single(start{trac}(first)))*d_para.dts;
                            prev_indy(trac,1:r_para.run_samp_lengths(samp_ruc))=first;
                        end
                    end
                end
            end

            if any(f_para.subplot_posi([m_para.future_spike_samp m_para.spike_samp]))
                if r_para.num_samp_runs==1
                    if max_num_uspikes<256
                        foll_indy=zeros(d_para.num_trains,m_res.num_samples,'uint8');
                    elseif max_num_uspikes<65536
                        foll_indy=zeros(d_para.num_trains,m_res.num_samples,'uint16');
                    elseif max_num_uspikes<2^32
                        foll_indy=zeros(d_para.num_trains,m_res.num_samples,'uint32');
                    else
                        foll_indy=zeros(d_para.num_trains,m_res.num_samples,'uint64');
                    end
                    foll_indy2=foll_indy;
                    folli=zeros(d_para.num_trains,m_res.num_samples,'single');
                    for trac=1:d_para.num_trains
                        for sc=1:num_uspikes(trac)-1
                            folli(trac,start{trac}(sc):ende{trac}(sc))=((double(ende{trac}(sc)-start{trac}(sc)))*...
                                d_para.dts:-d_para.dts:0);           % distance to following spike
                            foll_indy(trac,start{trac}(sc):ende{trac}(sc))=sc+1;                     % index of following spike
                        end
                        foll_indy2(trac,1:m_res.num_samples)=foll_indy(trac,1:m_res.num_samples);
                        if num_uspikes(trac)>3                                                         % $$$$$$$$$$$$
                            foll_indy2(trac,foll_indy2(trac,:)==num_uspikes(trac))=num_uspikes(trac)-1;
                        end
                    end
                    folli=folli(:,d_para.dsf:d_para.dsf:end);
                    foll_indy=foll_indy(:,d_para.dsf:d_para.dsf:end);
                    foll_indy2=foll_indy2(:,d_para.dsf:d_para.dsf:end);
                else
                    folli=zeros(d_para.num_trains,r_para.run_samp_lengths(samp_ruc),'single');            % distance to following spike
                    if max_num_uspikes<256
                        foll_indy=zeros(d_para.num_trains,r_para.run_samp_lengths(samp_ruc),'uint8');
                    elseif max_num_uspikes<65536
                        foll_indy=zeros(d_para.num_trains,r_para.run_samp_lengths(samp_ruc),'uint16');
                    elseif max_num_uspikes<2^32
                        foll_indy=zeros(d_para.num_trains,r_para.run_samp_lengths(samp_ruc),'uint32');
                    else
                        foll_indy=zeros(d_para.num_trains,r_para.run_samp_lengths(samp_ruc),'uint64');
                    end
                    for trac=1:d_para.num_trains
                        first=find(start{trac}(1:num_uspikes(trac)-1)<=r_para.run_samp_starts(samp_ruc)*d_para.dsf,1,'last');
                        last=find(ende{trac}(1:num_uspikes(trac)-1)>=r_para.run_samp_ends(samp_ruc)*d_para.dsf,1,'first');
                        if first<last
                            vect1=(r_para.run_samp_starts(samp_ruc):fix(ende{trac}(first)/d_para.dsf));
                            folli(trac,1:length(vect1))=(single(ende{trac}(first))-single(vect1*d_para.dsf)+1)*d_para.dts;
                            foll_indy(trac,1:length(vect1))=first+1;
                            for sc=first+1:last-1
                                vect1=(ceil(start{trac}(sc)/d_para.dsf):fix(ende{trac}(sc)/d_para.dsf));
                                folli(trac,vect1-r_para.run_samp_starts(samp_ruc)+1)=(single(ende{trac}(sc))-single(vect1*d_para.dsf)+1)*d_para.dts;
                                foll_indy(trac,vect1-r_para.run_samp_starts(samp_ruc)+1)=sc+1;
                            end
                            vect1=(ceil(start{trac}(last)/d_para.dsf):r_para.run_samp_ends(samp_ruc));
                            folli(trac,vect1-r_para.run_samp_starts(samp_ruc)+1)=(single(ende{trac}(last))-single(vect1*d_para.dsf)+1)*d_para.dts;
                            foll_indy(trac,vect1-r_para.run_samp_starts(samp_ruc)+1)=last+1;
                        else
                            folli(trac,1:r_para.run_samp_lengths(samp_ruc))=(single(ende{trac}(last))-...
                                single((r_para.run_samp_starts(samp_ruc):r_para.run_samp_ends(samp_ruc))*d_para.dsf)+1)*d_para.dts;
                            foll_indy(trac,1:r_para.run_samp_lengths(samp_ruc))=last+1;
                        end
                    end
                end
            end

            if ~exist('ustart','var')
                ustart=unique([start{:}]);
                uende=unique([ende{:}]);
            end
            if r_para.num_samp_runs==1
                isi_indy=zeros(1,m_res.num_measure_samples,'uint64');
                for ic=1:length(ustart)
                    isi_indy(ustart(ic):uende(ic))=ic;
                end
                isi_indy=isi_indy(d_para.dsf:d_para.dsf:end);
            else
                isi_indy=zeros(1,r_para.run_samp_lengths(samp_ruc),'uint64');
                first=find(ustart<=r_para.run_samp_starts(samp_ruc)*d_para.dsf,1,'last');
                last=find(uende>=r_para.run_samp_ends(samp_ruc)*d_para.dsf,1,'first');
                if first<last
                    isi_indy(1:length(r_para.run_samp_starts(samp_ruc)*d_para.dsf:d_para.dsf:uende(first)))=first;
                    for sc=first+1:last-1
                        isi_indy((ceil(ustart(sc)/d_para.dsf):fix(uende(sc)/d_para.dsf))-r_para.run_samp_starts(samp_ruc)+1)=sc;
                    end
                    isi_indy((ceil(ustart(last)/d_para.dsf):r_para.run_samp_ends(samp_ruc))-r_para.run_samp_starts(samp_ruc)+1)=last;
                else
                    isi_indy(1:r_para.run_samp_lengths(samp_ruc))=last;
                end
            end

            % ######################################################################################

            m_res.samp_measures_mat=zeros(m_para.num_samp_measures,d_para.num_pairs,r_para.run_samp_lengths(samp_ruc),'single');

            if f_para.subplot_posi(m_para.spike_samp)                                         % SPIKE-SAMP
                m_res.samp_measures_mat(m_para.measure_indy(m_para.spike_samp),1:d_para.num_pairs,1:r_para.run_samp_lengths(samp_ruc)) = ...
                    SPIKY_SPIKE_MEX(int32(d_para.num_pairs),int32(r_para.run_samp_lengths(samp_ruc)),int32(d_para.num_trains),...
                    folli,previ,int32(isi_indy),int32(prev_indy),int32(foll_indy),int32(r_para.run_pico_starts(samp_ruc)),...
                    ints(:,r_para.run_pico_starts(samp_ruc):r_para.run_pico_ends(samp_ruc)),udists);
            end

            if f_para.subplot_posi(m_para.realtime_spike_samp)                                       % REALTIME-SAMP
                m_res.samp_measures_mat(m_para.measure_indy(m_para.realtime_spike_samp),1:d_para.num_pairs,1:r_para.run_samp_lengths(samp_ruc)) = ...
                    SPIKY_realtimeSPIKE_MEX(int32(d_para.num_pairs),int32(r_para.run_samp_lengths(samp_ruc)),int32(d_para.num_trains),...
                    previ,int32(prev_indy),udists);
            end

            if f_para.subplot_posi(m_para.future_spike_samp)                                         % FUTURE-SAMP
                m_res.samp_measures_mat(m_para.measure_indy(m_para.future_spike_samp),1:d_para.num_pairs,1:r_para.run_samp_lengths(samp_ruc)) = ...
                    SPIKY_futureSPIKE_MEX(int32(d_para.num_pairs),int32(r_para.run_samp_lengths(samp_ruc)),int32(d_para.num_trains),...
                    folli,int32(foll_indy),udists);
            end
            ave_bi_vect(m_para.samp_bi_measures_indy,:)=ave_bi_vect(m_para.samp_bi_measures_indy,:)+mean(m_res.samp_measures_mat,3)*r_para.run_samp_lengths(samp_ruc);

            % ######################################################################################

            if r_para.num_samp_runs>1
                disp(['Calculation-Loop-info = ',num2str(r_para.num_samp_runs+1-samp_ruc),'  (',num2str(r_para.num_samp_runs),')'])
                save (['SPIKY_samp_AZBYCX_',num2str(samp_ruc)],'-struct','m_res','samp_measures_mat')
                waitbar((r_para.num_samp_runs+1-samp_ruc)/r_para.num_samp_runs,swbh,['Calculation-Loop-Info: ',num2str(r_para.num_samp_runs+1-samp_ruc),'  (',num2str(r_para.num_samp_runs),')'])
                if samp_ruc==1
                    delete(swbh)
                end
            end
        end
        clear start; clear ende;
        clear ustart; clear uende;
    else
        samp_ruc=1;
    end
    r_para.samp_ruc=samp_ruc;

    clear uspikes all_isi all_trains

    if m_para.num_pico_measures>0
        ave_bi_vect(m_para.pico_bi_measures_indy,:)=ave_bi_vect(m_para.pico_bi_measures_indy,:)/sum(m_res.isi);
    end
    if m_para.num_samp_measures>0
        ave_bi_vect(m_para.samp_bi_measures_indy,:)=ave_bi_vect(m_para.samp_bi_measures_indy,:)/m_res.num_measure_samples;
    end
    all_distances=mean(ave_bi_vect,2)';
    mat_indy=nchoosek(1:d_para.num_trains,2);

    results=[];
    results.spikes=spikes;
    for mac=1:m_para.num_sel_bi_measures
        eval(['results.',char(m_para.all_measures_string(m_para.select_bi_measures(mac))),...
            '.name=''',char(m_para.all_measures_string(m_para.select_bi_measures(mac))),''';'])
        eval(['results.',char(m_para.all_measures_string(m_para.select_bi_measures(mac))),...
            '.distance=all_distances(mac);'])
        eval(['results.',char(m_para.all_measures_string(m_para.select_bi_measures(mac))),...
            '.matrix=zeros(d_para.num_trains,d_para.num_trains,''single'');'])
        eval(['results.',char(m_para.all_measures_string(m_para.select_bi_measures(mac))),...
            '.matrix(sub2ind([d_para.num_trains d_para.num_trains],mat_indy(:,1),mat_indy(:,2)))=ave_bi_vect(mac,:);'])
        eval(['results.',char(m_para.all_measures_string(m_para.select_bi_measures(mac))),...
            '.matrix(sub2ind([d_para.num_trains d_para.num_trains],mat_indy(:,2),mat_indy(:,1)))=ave_bi_vect(mac,:);'])
        if ismember(m_para.select_bi_measures(mac),m_para.pico_measures)
            eval(['results.',char(m_para.all_measures_string(m_para.select_bi_measures(mac))),...
                '.x=m_res.isi;'])
            eval(['results.',char(m_para.all_measures_string(m_para.select_bi_measures(mac))),...
                '.y=shiftdim(permute(mean(m_res.pico_measures_mat(m_para.measure_indy(m_para.select_bi_measures(mac)),:,:),2),[2 1 3]),1);'])
        elseif ismember(m_para.select_bi_measures(mac),m_para.pili_measures)
            eval(['results.',char(m_para.all_measures_string(m_para.select_bi_measures(mac))),...
                '.x=m_res.pili_supi;'])
            eval(['results.',char(m_para.all_measures_string(m_para.select_bi_measures(mac))),...
                '.y=shiftdim(permute(mean(m_res.pili_measures_mat(m_para.measure_indy(m_para.select_bi_measures(mac)),:,:),2),[2 1 3]),1);'])
        elseif ismember(m_para.select_bi_measures(mac),m_para.samp_measures)
            eval(['results.',char(m_para.all_measures_string(m_para.select_bi_measures(mac))),...
                '.x=single(r_para.run_samp_starts:d_para.dts:r_para.run_samp_ends);'])
            eval(['results.',char(m_para.all_measures_string(m_para.select_bi_measures(mac))),...
                '.y=shiftdim(permute(mean(m_res.samp_measures_mat(m_para.measure_indy(m_para.select_bi_measures(mac)),:,:),2),[2 1 3]),1);'])
        end
    end

    if isfield(d_para,'instants') && ~isempty(d_para.instants)
        d_para.instants=unique(d_para.instants);
        d_para.num_instants=length(d_para.instants);
    else
        d_para.instants=[];
        d_para.num_instants=0;
    end

    if isfield(d_para,'selective_averages') && ~isempty(d_para.selective_averages)
        d_para.num_selective_averages=length(d_para.selective_averages);
        d_para.num_sel_ave=zeros(1,d_para.num_selective_averages);
        for sac=1:d_para.num_selective_averages
            d_para.num_sel_ave(sac)=length(d_para.selective_averages{sac})/2;
        end
    else
        d_para.num_selective_averages=0;
    end

    if isfield(d_para,'triggered_averages') && ~isempty(d_para.triggered_averages)
        d_para.num_triggered_averages=length(d_para.triggered_averages);
    else
        d_para.num_triggered_averages=0;
    end
    f_para.select_trains=1:d_para.num_trains;
    f_para.group_vect=d_para.group_vect;
    f_para.max_total_spikes=d_para.max_total_spikes;
end

