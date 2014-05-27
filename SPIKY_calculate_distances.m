% This calculates the selected time-resolved spike train distances once the button 'Calculate' is pressed.
% In case of very long datasets it also initiates the memory management.

disp(' '); disp(' ')

%
% ==================================
%  To do - List:
% ===============
% Stimulus
% reset in STG
% Movie (new Matlab-function)
% Poisson line
% effect of edge correction on ISI-distance
% More than 2 measures and only matrices in frame comparison
% ==================================
% 

m_para.all_measures_str={'Stimulus';'Spikes';'PS';'I'; 'S';'S_r';'S_f';'Si';'Si_r';'Si_f';'V';'R'};
m_para.all_measures_string={'Stimulus';'Spikes';'PSTH';'ISI'; 'SPIKE';'SPIKE_realtime';'SPIKE_future';...
    'SPIKEi';'SPIKEi_realtime';'SPIKEi_future';'Victor';'vanRossum'};

m_para.num_diff_measures=2;
m_para.num_non_bi_measures=3;

m_para.stimulus=1;
m_para.spikes=2;
m_para.psth=3;
m_para.isi_pico=4;
m_para.spike_pili=5;
m_para.realtime_spike_pili=6;
m_para.future_spike_pili=7;
m_para.spike_pico=8;
m_para.realtime_spike_pico=9;
m_para.future_spike_pico=10;
m_para.victor=11;
m_para.van_rossum=12;

m_para.pili_measures=[m_para.spike_pili m_para.realtime_spike_pili m_para.future_spike_pili];
m_para.pico_measures=[m_para.isi_pico m_para.spike_pico m_para.realtime_spike_pico m_para.future_spike_pico];
m_para.realtime_measures=[m_para.realtime_spike_pico m_para.realtime_spike_pili];
m_para.future_measures=[m_para.future_spike_pico m_para.future_spike_pili];
m_para.bi_measures=[m_para.pico_measures m_para.pili_measures];

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
m_para.select_pico_measures=intersect(m_para.pico_measures,m_para.select_measures);
m_para.num_pico_measures=length(m_para.select_pico_measures);

m_para.measure_indy=zeros(1,m_para.num_all_measures);
m_para.measure_indy(m_para.pili_measures)=SPIKY_f_all_sorted(f_para.subplot_posi(m_para.pili_measures));
m_para.measure_indy(m_para.pico_measures)=SPIKY_f_all_sorted(f_para.subplot_posi(m_para.pico_measures));

aaa=f_para.subplot_posi(m_para.num_diff_measures+1:m_para.num_all_measures);
bbb=aaa+(1:(m_para.num_all_measures-m_para.num_diff_measures))/(m_para.num_all_measures-m_para.num_diff_measures+1);
bbb(bbb<1)=0;
m_para.measure_all_indy=[zeros(1,m_para.num_diff_measures) SPIKY_f_all_sorted(bbb)];

ccc=aaa(aaa>0)-min(aaa(aaa>0))+1;
[dummy,eee]=sort(ccc);
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
m_para.sel_pico_measures_str=m_para.all_measures_str(m_para.select_pico_measures);

m_para.pili_measures_indy=find(ismember(m_para.select_measures,m_para.pili_measures));
m_para.pico_measures_indy=find(ismember(m_para.select_measures,m_para.pico_measures));

m_para.pili_bi_measures_indy=find(ismember(m_para.select_bi_measures,m_para.pili_measures));
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

uspikes=cell(1,d_para.num_trains);
for trac=1:d_para.num_trains
    uspikes{trac}=spikes{trac}(spikes{trac}>=d_para.tmin & spikes{trac}<=d_para.tmax);
    uspikes{trac}=unique([d_para.tmin uspikes{trac} d_para.tmax]);
end
num_uspikes=cellfun('length',uspikes);
max_num_uspikes=max(num_uspikes);

% #########################################################################
% ######################################################################### Victor, van Rossum
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
    %m_res.vr2=zeros(length(tau),d_para.num_trains,d_para.num_trains);
    %m_res.vr3=zeros(length(tau),d_para.num_trains,d_para.num_trains);
    for tauc=1:length(tau)
        for trac1=1:d_para.num_trains-1
            for trac2=trac1+1:d_para.num_trains
                m_res.vr(tauc,trac1,trac2)=SPIKY_vanRossum(uspikes{trac1}/d_para.tmax*2,uspikes{trac2}/d_para.tmax*2,tau*2);
                %m_res.vr2(tauc,trac1,trac2)=SPIKY_vanRossum(uspikes{trac1}/d_para.tmax,uspikes{trac2}/d_para.tmax,tau);
                %m_res.vr3(tauc,trac1,trac2)=SPIKY_vanRossum(uspikes{trac1},uspikes{trac2},tau);
            end
        end
        m_res.vr(tauc,:,:)=shiftdim(m_res.vr(tauc,:,:),1)+shiftdim(m_res.vr(tauc,:,:),1)';
        %m_res.vr2(tauc,:,:)=shiftdim(m_res.vr2(tauc,:,:),1)+shiftdim(m_res.vr2(tauc,:,:),1)';
        %m_res.vr3(tauc,:,:)=shiftdim(m_res.vr3(tauc,:,:),1)+shiftdim(m_res.vr3(tauc,:,:),1)';   % still problematic !!!!!!!
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
    m_res.pili_supi_indy=round((m_res.pili_supi-d_para.tmin)/d_para.dts);
    clear dummy
            
    % ###########################################################################################################################################
    % ############################################################## Memory management ##########################################################
    % ###########################################################################################################################################
    
    m_para.memo_num_measures=2*m_para.num_pili_measures+m_para.num_pico_measures;
    memo_fact=m_para.memo_num_measures*d_para.num_pairs;
    memo=memo_fact*m_res.num_isi;
    
    if f_para.run_test==0 || ~ismember(get(handles.Data_listbox,'Value'),[1 2 3])
        max_memo=d_para.max_memo_init;
    elseif get(handles.Data_listbox,'Value')==1
        max_memo=30; % $$$$$$$$$$$
        %max_memo=10; % $$$$$$$$$$$
    elseif get(handles.Data_listbox,'Value')==2
        max_memo=800; % $$$$$$$$$$$
    elseif get(handles.Data_listbox,'Value')==3
        max_memo=800; % $$$$$$$$$$$
    end
    
    
    r_para.num_runs=1;
    if memo>max_memo
        num_init_runs=ceil(memo/max_memo);
        if f_para.run_test==0 || ~ismember(get(handles.Data_listbox,'Value'),[1 2 3])
            max_pico_len=fix(m_res.num_isi/num_init_runs);
        elseif get(handles.Data_listbox,'Value')==1
            max_pico_len=15; % in isi             % $$$$$$$$$$
            %max_pico_len=5; % in isi             % $$$$$$$$$$
        elseif get(handles.Data_listbox,'Value')==2
            max_pico_len=500; % in isi            % $$$$$$$$$$
        elseif get(handles.Data_listbox,'Value')==3
            max_pico_len=100; % in isi            % $$$$$$$$$$
        end
        
        if m_para.memo_num_measures>0
            r_para.num_runs=ceil(m_res.num_isi/max_pico_len);
            if r_para.num_runs==0
                set(0,'DefaultUIControlFontSize',16);
                mbh=msgbox(sprintf('Dataset might be too large.\nPlease increase the value of the variable\n ''d_para.max_memo_init'' in ''SPIKY_f_user_interface'' !!!'),'Warning','warn','modal');
                htxt = findobj(mbh,'Type','text');
                set(htxt,'FontSize',12,'FontWeight','bold')
                mb_pos=get(mbh,'Position');
                set(mbh,'Position',[mb_pos(1:2) mb_pos(3)*1.5 mb_pos(4)])
                uiwait(mbh);
                ret=1;
                return
            end
        end
        
        r_para.run_pico_ends=cumsum(fix([max_pico_len*ones(1,r_para.num_runs-1) m_res.num_isi-max_pico_len*(r_para.num_runs-1)]));
        r_para.run_pico_starts=[1 r_para.run_pico_ends(1:end-1)+1];
        r_para.run_pico_lengths=r_para.run_pico_ends-r_para.run_pico_starts+1;
        if m_para.num_pili_measures>0
            r_para.run_pili_ends=2*r_para.run_pico_ends;
            r_para.run_pili_starts=[1 r_para.run_pili_ends(1:end-1)+1];
            r_para.run_pili_lengths=r_para.run_pili_ends-r_para.run_pili_starts+1;
        end
    end
    
    if r_para.num_runs==1
        r_para.run_pico_lengths=m_res.num_isi;
        r_para.run_pico_starts=1;
        r_para.run_pico_ends=m_res.num_isi;
        if m_para.num_pili_measures>0
            r_para.run_pili_lengths=2*m_res.num_isi;
            r_para.run_pili_starts=1;
            r_para.run_pili_ends=2*m_res.num_isi;
        end
    end
    
    empties=find(all_isi==0);
    ivs=cell(1,d_para.num_trains);
    ive=cell(1,d_para.num_trains);
    for trac=1:d_para.num_trains
        dummy1=[1 find(all_trains==trac)];
        dummy2=[dummy1(2:length(dummy1))-1 num_all_isi];
        len=dummy2-dummy1+1-histc(empties,dummy1);
        ive{trac}=cumsum(len);
        ivs{trac}=[1 ive{trac}(1:end-1)+1];
    end
    clear all_isi all_trains empties
    
    ave_bi_vect=zeros(m_para.num_sel_bi_measures,d_para.num_pairs);
    if r_para.num_runs>1
        disp(' ');
        disp('Large data set. Please be patient.')
        disp(' ');
        disp(['Number of calculation loop runs: ',num2str(r_para.num_runs)])
        pwbh = waitbar(0,'Large data set. Please be patient.','CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
        setappdata(pwbh,'canceling',0)
    end
    
    if f_para.edge_correction 
        if r_para.num_runs>1    % Edge-correction for many runs
            fl_isis=repmat({zeros(1,4)},1,d_para.num_trains);
            for trac=1:d_para.num_trains
                if num_uspikes(trac)>3
                    fl_isis{trac}(1:2)=diff(uspikes{trac}(1:3));
                    fl_isis{trac}(3:4)=diff(uspikes{trac}(end-2:end));
                end
            end
        end
        if any(f_para.subplot_posi([m_para.spike_pico m_para.realtime_spike_pico m_para.spike_pili m_para.realtime_spike_pili]))
            prev_edge_cor_indy=cell(1,d_para.num_trains);
        end
        if any(f_para.subplot_posi([m_para.spike_pico m_para.future_spike_pico m_para.spike_pili m_para.future_spike_pili]))
            foll_edge_cor_indy=cell(1,d_para.num_trains);
        end
    end
    
    run_ivs=cell(1,d_para.num_trains);
    run_ive=cell(1,d_para.num_trains);
    for ruc=r_para.num_runs:-1:1
        if r_para.num_runs>1 && getappdata(pwbh,'canceling')
            delete(pwbh)
            ret=1;
            return
        end
        
        % ###########################################################################################################################################
        % ################################################################# Pico-Pili-Quantities #####################################################
        % ###########################################################################################################################################
        
        if any(f_para.subplot_posi(m_para.bi_measures))
            firsts=zeros(1,d_para.num_trains);
            lasts=zeros(1,d_para.num_trains);
            for trac=1:d_para.num_trains
                firsts(trac)=find(uspikes{trac}(1:num_uspikes(trac))<=m_res.cum_isi(r_para.run_pico_starts(ruc)),1,'last');
                lasts(trac)=find(uspikes{trac}(1:num_uspikes(trac))<m_res.cum_isi(r_para.run_pico_ends(ruc)+1),1,'last');
                run_ivs{trac}=ivs{trac}(firsts(trac):lasts(trac))-r_para.run_pico_starts(ruc)+1;
                run_ive{trac}=ive{trac}(firsts(trac):lasts(trac))-r_para.run_pico_starts(ruc)+1;
                run_ivs{trac}(run_ivs{trac}<1)=1;
                run_ive{trac}(run_ive{trac}>r_para.run_pico_lengths(ruc))=r_para.run_pico_lengths(ruc);
            end
            run_num_ints=lasts-firsts+1;
        end
        
        if any(f_para.subplot_posi([m_para.spike_pico m_para.realtime_spike_pico m_para.spike_pili m_para.realtime_spike_pili]))  % SPIKE-Pre-Pico
            previs_indy=zeros(d_para.num_trains,max(run_num_ints),'uint32');
            previs=zeros(d_para.num_trains,max(run_num_ints),'single');
            prev_spikes_indy=zeros(d_para.num_trains,r_para.run_pico_lengths(ruc),'uint32');
            prev_spikes=zeros(d_para.num_trains,r_para.run_pico_lengths(ruc),'single');
            for trac=1:d_para.num_trains
                previs_indy(trac,1:run_num_ints(trac))=firsts(trac):lasts(trac);
                previs(trac,1:run_num_ints(trac))=uspikes{trac}(previs_indy(trac,1:run_num_ints(trac)));
                for ic=1:run_num_ints(trac)
                    prev_spikes_indy(trac,run_ivs{trac}(ic):run_ive{trac}(ic))=previs_indy(trac,ic);
                    prev_spikes(trac,run_ivs{trac}(ic):run_ive{trac}(ic))=previs(trac,ic);
                end
                if f_para.edge_correction && num_uspikes(trac)>3 && firsts(trac)==1
                    prev_edge_cor_indy{trac}=find(prev_spikes_indy(trac,:)==1);
                end
                prev_spikes_indy(trac,:)=prev_spikes_indy(trac,:)-uint32(firsts(trac)-1);
            end
            clear previs_indy previs
            
            if any(f_para.subplot_posi([m_para.spike_pili m_para.realtime_spike_pili]))                     % SPIKE-Pre-Pili
                prev_spikes_indy_pili=zeros(d_para.num_trains,r_para.run_pili_lengths(ruc),'uint32');
                prev_spikes_pili=zeros(d_para.num_trains,r_para.run_pili_lengths(ruc),'single');
                for trac=1:d_para.num_trains
                    prev_spikes_indy_pili(trac,:)=reshape(repmat(prev_spikes_indy(trac,1:r_para.run_pico_lengths(ruc)),2,1),1,r_para.run_pili_lengths(ruc));
                    prev_spikes_pili(trac,:)=m_res.pili_supi(r_para.run_pili_starts(ruc):r_para.run_pili_ends(ruc))-...
                        reshape([prev_spikes(trac,1:r_para.run_pico_lengths(ruc)); prev_spikes(trac,1:r_para.run_pico_lengths(ruc))],...
                        1,r_para.run_pili_lengths(ruc));
                end
            end
        end
        
        if any(f_para.subplot_posi([m_para.spike_pico m_para.future_spike_pico m_para.spike_pili m_para.future_spike_pili]))  % SPIKE-Future-Pico
            follis_indy=zeros(d_para.num_trains,max(run_num_ints),'uint32');
            follis=zeros(d_para.num_trains,max(run_num_ints),'single');
            foll_spikes_indy=zeros(d_para.num_trains,r_para.run_pico_lengths(ruc),'uint32');
            foll_spikes=zeros(d_para.num_trains,r_para.run_pico_lengths(ruc),'single');
            for trac=1:d_para.num_trains
                follis_indy(trac,1:run_num_ints(trac))=(firsts(trac):lasts(trac))+1;
                follis(trac,1:run_num_ints(trac))=uspikes{trac}(follis_indy(trac,1:run_num_ints(trac)));
                for ic=1:run_num_ints(trac)
                    foll_spikes_indy(trac,run_ivs{trac}(ic):run_ive{trac}(ic))=follis_indy(trac,ic);
                    foll_spikes(trac,run_ivs{trac}(ic):run_ive{trac}(ic))=follis(trac,ic);
                end
                if f_para.edge_correction && num_uspikes(trac)>3 && lasts(trac)==num_uspikes(trac)-1
                    foll_edge_cor_indy{trac}=find(foll_spikes_indy(trac,:)==num_uspikes(trac));
                end
                foll_spikes_indy(trac,:)=foll_spikes_indy(trac,:)-uint32(firsts(trac)-1);
            end
            clear follis_indy follis
            
            if any(f_para.subplot_posi([m_para.spike_pili m_para.future_spike_pili]))                       % SPIKE-Future-Pili
                foll_spikes_indy_pili=zeros(d_para.num_trains,r_para.run_pili_lengths(ruc),'uint32');
                foll_spikes_pili=zeros(d_para.num_trains,r_para.run_pili_lengths(ruc),'single');
                for trac=1:d_para.num_trains
                    foll_spikes_indy_pili(trac,:)=reshape(repmat(foll_spikes_indy(trac,1:r_para.run_pico_lengths(ruc)),2,1),1,r_para.run_pili_lengths(ruc));
                    foll_spikes_pili(trac,:)=reshape([foll_spikes(trac,1:r_para.run_pico_lengths(ruc)); foll_spikes(trac,1:r_para.run_pico_lengths(ruc))],...
                        1,r_para.run_pili_lengths(ruc))-m_res.pili_supi(r_para.run_pili_starts(ruc):r_para.run_pili_ends(ruc));
                end
            end
        end
        
        run_uspikes=cell(1,d_para.num_trains);
        for trac=1:d_para.num_trains
            run_uspikes{trac}=uspikes{trac}(int32(firsts(trac)):int32(lasts(trac)+1));
        end
        run_num_uspikes=run_num_ints+1;
        if exist(['SPIKY_udists_MEX.',mexext],'file')
            run_udists=SPIKY_udists_MEX(int32(d_para.num_trains),int32(run_num_uspikes),run_uspikes);
        else
            run_udists=cell(d_para.num_trains);
            for trac1=1:d_para.num_trains
                for trac2=setdiff(1:d_para.num_trains,trac1)
                    run_udists{trac1,trac2}=zeros(1,run_num_uspikes(trac1),'single');
                    for spc=1:run_num_uspikes(trac1)
                        run_udists{trac1,trac2}(spc)=min(abs(run_uspikes{trac1}(spc)-run_uspikes{trac2}));
                    end
                end
            end
        end
        for trac=1:d_para.num_trains
            for trac2=1:d_para.num_trains
                if (trac~=trac2)
                    run_udists{trac,trac2}(1)=min(abs(uspikes{trac}(firsts(trac))-uspikes{trac2}));
                    run_udists{trac,trac2}(end)=min(abs(uspikes{trac}(lasts(trac)+1)-uspikes{trac2}));
                end
            end
        end
        
        if any(f_para.subplot_posi([m_para.isi_pico m_para.spike_pili]))
            isis=cell(1,d_para.num_trains);
            ints=zeros(d_para.num_trains,r_para.run_pico_lengths(ruc),'single');
            for trac=1:d_para.num_trains
                isis{trac}=diff(uspikes{trac}(firsts(trac):lasts(trac)+1));
                for ic=1:run_num_ints(trac)
                    ints(trac,run_ivs{trac}(ic):run_ive{trac}(ic))=isis{trac}(ic);
                end
            end
        end
        
        % #####################################################################################################################################
        % ########################################################### Pico-Measures I (ISI) ###################################################
        % #####################################################################################################################################

        if m_para.num_pico_measures>0                                                                                                   % Pico
            m_res.pico_measures_mat=zeros(m_para.num_pico_measures,d_para.num_pairs,r_para.run_pico_lengths(ruc),'single');
            
            if f_para.subplot_posi(m_para.isi_pico)                        % ISI (calculated first, then edge-correction of ints for SPIKE)
                if exist(['SPIKY_ISI_MEX.',mexext],'file')
                    m_res.pico_measures_mat(m_para.measure_indy(m_para.isi_pico),1:d_para.num_pairs,1:r_para.run_pico_lengths(ruc)) = ...
                        abs(SPIKY_ISI_MEX(int32(d_para.num_pairs),int32(r_para.run_pico_lengths(ruc)),int32(d_para.num_trains),ints));
                else
                    pac=0;
                    for trac1=1:d_para.num_trains-1
                        for trac2=trac1+1:d_para.num_trains
                            pac=pac+1;
                            dummy1=find(ints(trac1,:)<ints(trac2,:));
                            m_res.pico_measures_mat(m_para.measure_indy(m_para.isi_pico),pac,dummy1)=abs(ints(trac1,dummy1)./ints(trac2,dummy1)-1);
                            dummy2=find(ints(trac1,:)>=ints(trac2,:) & ints(trac1,:)~=0);
                            m_res.pico_measures_mat(m_para.measure_indy(m_para.isi_pico),pac,dummy2)=abs(ints(trac2,dummy2)./ints(trac1,dummy2)-1);
                        end
                    end
                end
            end
        end
        
        % #####################################################################################################################################
        % ################################################################# Pili-Measures #####################################################
        % #####################################################################################################################################

        if m_para.num_pili_measures>0
            odds=1:2:r_para.run_pili_lengths(ruc);
            evens=odds+1;
            m_res.pili_measures_mat=zeros(m_para.num_pili_measures,d_para.num_pairs,r_para.run_pili_lengths(ruc),'single');

            if f_para.subplot_posi(m_para.realtime_spike_pili)                                           % REALTIME-Pili
                if exist(['SPIKY_realtimeSPIKE_MEX.',mexext],'file')
                   m_res.pili_measures_mat(m_para.measure_indy(m_para.realtime_spike_pili),1:d_para.num_pairs,1:r_para.run_pili_lengths(ruc)) = ...
                        SPIKY_realtimeSPIKE_MEX(int32(d_para.num_pairs),int32(r_para.run_pili_lengths(ruc)),int32(d_para.num_trains),...
                        prev_spikes_pili,int32(prev_spikes_indy_pili),run_udists);
                else
                    pac=0;
                    for trac1=1:d_para.num_trains-1
                        for trac2=trac1+1:d_para.num_trains
                            pac=pac+1;
                            normy=(prev_spikes_pili(trac1,:)+prev_spikes_pili(trac2,:));
                            dummy=(prev_spikes_pili(trac1,:)<prev_spikes_pili(trac2,:)-0.00000001);
                            dummy1=trac1*(1-dummy)+trac2*dummy;   % index of spike train with earlier spike
                            dummy2=trac2*(1-dummy)+trac1*dummy;   % index of spike train with later spike
                            for sc=1:r_para.run_pili_lengths(ruc)
                                m_res.pili_measures_mat(m_para.measure_indy(m_para.realtime_spike_pili),pac,sc)= ...
                                    (abs(prev_spikes_pili(trac1,sc)-prev_spikes_pili(trac2,sc))+...   % later spike
                                    run_udists{dummy1(sc),dummy2(sc)}(prev_spikes_indy_pili(dummy1(sc),sc)))...      % earlier spike
                                    /(2*normy(sc)+(normy(sc)==0));
                            end
                        end
                    end
                    clear dummy dummy1 dummy2
                end
                aves=(log(1./m_res.pili_measures_mat(m_para.measure_indy(m_para.realtime_spike_pili),:,evens))-...
                    log(1./m_res.pili_measures_mat(m_para.measure_indy(m_para.realtime_spike_pili),:,odds)))./...
                    (1./m_res.pili_measures_mat(m_para.measure_indy(m_para.realtime_spike_pili),:,evens)-...
                    1./m_res.pili_measures_mat(m_para.measure_indy(m_para.realtime_spike_pili),:,odds));
                aves(isnan(aves))=0;
                ave_bi_vect(logical(m_para.select_bi_measures==m_para.realtime_spike_pili),:)=...
                    sum(aves.*repmat(shiftdim(m_res.isi(r_para.run_pico_starts(ruc):r_para.run_pico_ends(ruc)),-1),...
                    [1,d_para.num_pairs]),3)/sum(m_res.isi(r_para.run_pico_starts(ruc):r_para.run_pico_ends(ruc)));
            end
            
            if f_para.subplot_posi(m_para.future_spike_pili)                                             % FUTURE-Pili
                if exist(['SPIKY_futureSPIKE_MEX.',mexext],'file')
                    m_res.pili_measures_mat(m_para.measure_indy(m_para.future_spike_pili),1:d_para.num_pairs,1:r_para.run_pili_lengths(ruc)) = ...
                        SPIKY_futureSPIKE_MEX(int32(d_para.num_pairs),int32(r_para.run_pili_lengths(ruc)),int32(d_para.num_trains),...
                        foll_spikes_pili,int32(foll_spikes_indy_pili),run_udists);
                else
                    pac=0;
                    for trac1=1:d_para.num_trains-1
                        for trac2=trac1+1:d_para.num_trains
                            pac=pac+1;
                            normy=(foll_spikes_pili(trac1,:)+foll_spikes_pili(trac2,:));
                            dummy=(foll_spikes_pili(trac1,:)<foll_spikes_pili(trac2,:)-0.00000001);
                            dummy1=trac1*(1-dummy)+trac2*dummy;   % index of spike train with earlier spike
                            dummy2=trac2*(1-dummy)+trac1*dummy;   % index of spike train with later spike
                            for sc=1:r_para.run_pili_lengths(ruc)
                                m_res.pili_measures_mat(m_para.measure_indy(m_para.future_spike_pili),pac,sc)= ...
                                    (abs(foll_spikes_pili(trac1,sc)-foll_spikes_pili(trac2,sc))+...   % later spike
                                    run_udists{dummy1(sc),dummy2(sc)}(foll_spikes_indy_pili(dummy1(sc),sc)))...      % earlier spike
                                    /(2*normy(sc)+(normy(sc)==0));
                            end
                        end
                    end
                    clear dummy dummy1 dummy2
                end
                aves=(log(1./m_res.pili_measures_mat(m_para.measure_indy(m_para.future_spike_pili),:,evens))-...
                    log(1./m_res.pili_measures_mat(m_para.measure_indy(m_para.future_spike_pili),:,odds)))./...
                    (1./m_res.pili_measures_mat(m_para.measure_indy(m_para.future_spike_pili),:,evens)-...
                    1./m_res.pili_measures_mat(m_para.measure_indy(m_para.future_spike_pili),:,odds));
                aves(isnan(aves))=0;
                ave_bi_vect(logical(m_para.select_bi_measures==m_para.future_spike_pili),:)=...
                    sum(aves.*repmat(shiftdim(m_res.isi(r_para.run_pico_starts(ruc):r_para.run_pico_ends(ruc)),-1),...
                    [1,d_para.num_pairs]),3)/sum(m_res.isi(r_para.run_pico_starts(ruc):r_para.run_pico_ends(ruc)));
            end
            
            if f_para.subplot_posi(m_para.spike_pili)                                                   % SPIKE-Pili
                if f_para.edge_correction
                    for trac=1:d_para.num_trains
                        if num_uspikes(trac)>3
                            if r_para.num_runs==1
                                isis{trac}(1)=max(isis{trac}(1:2));
                                isis{trac}(end)=max(isis{trac}(end-1:end));
                            else
                                if firsts(trac)==1
                                    isis{trac}(1)=max(fl_isis{trac}(1:2));
                                end
                                if lasts(trac)==num_uspikes(trac)-1
                                    isis{trac}(end)=max(fl_isis{trac}(3:4));
                                end
                            end
                            for ic=1:run_num_ints(trac)
                                ints(trac,run_ivs{trac}(ic):run_ive{trac}(ic))=isis{trac}(ic);
                            end
                        end
                        if firsts(trac)==1
                            %prev_spikes_indy(prev_edge_cor_indy{trac})=prev_spikes_indy(prev_edge_cor_indy{trac})+1;
                            prev_spikes_indy_pili(trac,[prev_edge_cor_indy{trac}*2-1 prev_edge_cor_indy{trac}*2])=...
                                prev_spikes_indy_pili(trac,[prev_edge_cor_indy{trac}*2-1 prev_edge_cor_indy{trac}*2])+1;
                        end
                        if lasts(trac)==num_uspikes(trac)-1
                            %foll_spikes_indy(foll_edge_cor_indy{trac})=foll_spikes_indy(foll_edge_cor_indy{trac})-1;
                            foll_spikes_indy_pili(trac,[foll_edge_cor_indy{trac}*2-1 foll_edge_cor_indy{trac}*2])=...
                                foll_spikes_indy_pili(trac,[foll_edge_cor_indy{trac}*2-1 foll_edge_cor_indy{trac}*2])-1;
                        end
                    end                    
                end
                clear isis
                isi_indy_pili=reshape(repmat(r_para.run_pico_starts(ruc):r_para.run_pico_ends(ruc),2,1),1,2*r_para.run_pico_lengths(ruc))-...
                    r_para.run_pico_starts(ruc)+1;

                if exist(['SPIKY_SPIKE_MEX.',mexext],'file')
                    m_res.pili_measures_mat(m_para.measure_indy(m_para.spike_pili),1:d_para.num_pairs,1:r_para.run_pili_lengths(ruc)) = ...
                        SPIKY_SPIKE_MEX(int32(d_para.num_pairs),int32(r_para.run_pili_lengths(ruc)),int32(d_para.num_trains),...
                        foll_spikes_pili,prev_spikes_pili,int32(isi_indy_pili-1),int32(prev_spikes_indy_pili),int32(foll_spikes_indy_pili),...
                        ints,run_udists);
                else
                    pac=0;
                    for trac1=1:d_para.num_trains-1
                        for trac2=trac1+1:d_para.num_trains
                            pac=pac+1;
                            m_res.pili_measures_mat(m_para.measure_indy(m_para.spike_pili),pac,1:r_para.run_pili_lengths(ruc)) = ...
                                ((run_udists{trac1,trac2}(prev_spikes_indy_pili(trac1,:)).*foll_spikes_pili(trac1,:)+...
                                run_udists{trac1,trac2}(foll_spikes_indy_pili(trac1,:)).*prev_spikes_pili(trac1,:))./...
                                ints(trac1,isi_indy_pili).*ints(trac2,isi_indy_pili)+...
                                (run_udists{trac2,trac1}(prev_spikes_indy_pili(trac2,:)).*foll_spikes_pili(trac2,:)+...
                                run_udists{trac2,trac1}(foll_spikes_indy_pili(trac2,:)).*prev_spikes_pili(trac2,:))./...
                                ints(trac2,isi_indy_pili).*ints(trac1,isi_indy_pili))./...
                                ((ints(trac1,isi_indy_pili)+ints(trac2,isi_indy_pili)).^2/2);
                        end
                    end
                end
                ave_bi_vect(logical(m_para.select_bi_measures==m_para.spike_pili),:)=...
                    sum((m_res.pili_measures_mat(m_para.measure_indy(m_para.spike_pili),:,odds)+...
                    m_res.pili_measures_mat(m_para.measure_indy(m_para.spike_pili),:,evens))/2.*...
                    repmat(shiftdim(m_res.isi(r_para.run_pico_starts(ruc):r_para.run_pico_ends(ruc)),-1),[1,d_para.num_pairs]),3)/...
                    sum(m_res.isi(r_para.run_pico_starts(ruc):r_para.run_pico_ends(ruc)));
            end
        end
        
        % #####################################################################################################################################
        % ########################################################### [ Pico-Measures II ] ####################################################
        % #####################################################################################################################################
        
        if m_para.num_pico_measures>0                                                                                   % Pico
            if f_para.subplot_posi(m_para.spike_pico)                                         % SPIKE-Pico
                if r_para.ruc==r_para.num_runs
                    m_res.isi_pos=cumsum([d_para.tmin m_res.isi(1:end-1)])+m_res.isi/2;
                end
                m_res.pico_measures_mat(m_para.measure_indy(m_para.spike_pico),1:d_para.num_pairs,1:r_para.run_pico_lengths(ruc)) = ...
                    SPIKY_SPIKEpico_MEX(int32(d_para.num_pairs),int32(r_para.run_pico_lengths(ruc)),int32(d_para.num_trains),...
                    foll_spikes,m_res.isi_pos(:,r_para.run_pico_starts(ruc):r_para.run_pico_ends(ruc)),int32(prev_spikes_indy),...
                    ints,int32(foll_spikes_indy),prev_spikes,run_udists);
            end
            
            if f_para.subplot_posi(m_para.realtime_spike_pico)                                  % REALTIME-PICO
                m_res.pico_measures_mat(m_para.measure_indy(m_para.realtime_spike_pico),1:d_para.num_pairs,1:r_para.run_pico_lengths(ruc)) = ...
                    SPIKY_realtimeSPIKEpico_MEX(int32(d_para.num_pairs),int32(r_para.run_pico_lengths(ruc)),int32(d_para.num_trains),...
                    prev_spikes,m_res.cum_isi(:,[r_para.run_pico_starts(ruc):r_para.run_pico_ends(ruc) r_para.run_pico_ends(ruc)+1]),...
                    m_res.isi(:,r_para.run_pico_starts(ruc):r_para.run_pico_ends(ruc)),int32(prev_spikes_indy),run_udists);
            end
            
            if f_para.subplot_posi(m_para.future_spike_pico)                                    % FUTURE-PICO
                m_res.pico_measures_mat(m_para.measure_indy(m_para.future_spike_pico),1:d_para.num_pairs,1:r_para.run_pico_lengths(ruc)) = ...
                    SPIKY_futureSPIKEpico_MEX(int32(d_para.num_pairs),int32(r_para.run_pico_lengths(ruc)),int32(d_para.num_trains),...
                    foll_spikes,m_res.cum_isi(:,[r_para.run_pico_starts(ruc):r_para.run_pico_ends(ruc) r_para.run_pico_ends(ruc)+1]),...
                    m_res.isi(:,r_para.run_pico_starts(ruc):r_para.run_pico_ends(ruc)),int32(foll_spikes_indy),run_udists);
            end
            
            ave_bi_vect(m_para.pico_bi_measures_indy,:)=ave_bi_vect(m_para.pico_bi_measures_indy,:)+sum(m_res.pico_measures_mat.*...
                repmat(shiftdim(m_res.isi(r_para.run_pico_starts(ruc):r_para.run_pico_ends(ruc)),-1),[m_para.num_pico_measures,...
                d_para.num_pairs]),3);
        end
        
        
        if r_para.num_runs>1
            disp(['Calculation-Loop-Info: = ',num2str(r_para.num_runs+1-ruc),'  (',num2str(r_para.num_runs),')'])
            if m_para.num_pico_measures>0
                eval(['pico' num2str(ruc) '= m_res.pico_measures_mat;']);
            end
            if m_para.num_pili_measures>0
                eval(['pili' num2str(ruc) '= m_res.pili_measures_mat;']);
            end
            if ruc==r_para.num_runs
                save('SPIKY_AZBYCX',['pi*' num2str(ruc)])
            else
                save('SPIKY_AZBYCX',['pi*' num2str(ruc)],'-append')
            end
            waitbar((r_para.num_runs+1-ruc)/r_para.num_runs,pwbh,['Calculation-Loop-Info: ',...
                num2str(r_para.num_runs+1-ruc),'  (',num2str(r_para.num_runs),')'])
            if ruc==1
                delete(pwbh)
            end
        end
    end
    r_para.ruc=ruc;
    clear uspikes
    
    if m_para.num_pico_measures>0
        ave_bi_vect(m_para.pico_bi_measures_indy,:)=ave_bi_vect(m_para.pico_bi_measures_indy,:)/sum(m_res.isi);
    end
    all_distances=mean(ave_bi_vect,2)';
    mat_indy=nchoosek(1:d_para.num_trains,2);
    
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
else
    results=[];
    r_para=[];
    m_para.memo_num_measures=0;
end
f_para.max_total_spikes=d_para.max_total_spikes;

