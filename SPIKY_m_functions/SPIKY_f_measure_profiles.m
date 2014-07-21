%samp This function plots the individual time profiles of the selected spike train distances and calculates their average.

function [fave,yt,ytl,p_para]=SPIKY_f_measure_profiles(fave,fyt,fytl,fx,fy,fheadline,...
    fsp_paras,fminval,fmaxval,profi_ave,fdata_type,s_para,p_para,mac)

fsp_posi=fsp_paras(1);
fsp_size=fsp_paras(2);
fsp_start=fsp_paras(3);

for nmac=1:length(s_para.nma)
    mova=0;
    if s_para.nma(nmac)==2 && fdata_type==3
        continue
    end
    if s_para.nma(nmac)==2 && size(fy,1)==1                                                       % moving average
        if s_para.causal==0                                                % regular (average over both past and future values)
            if fdata_type==1          % piecewise constant
                for trac=1:size(fy,1)
                    fy(trac,:)=SPIKY_f_moving_average_weighted(fy(trac,:),fx,s_para.mao);
                end
            elseif fdata_type==2      % piecewise linear
                mova=1;
                for trac=1:size(fy,1)
                    %fx=fx(2:2:end)-fx([1 2:2:end-2]);
                    fy(trac,1:length(fx))=mean([fy(trac,2:2:end); fy(trac,1:2:end)]);
                    fy(trac,1:length(fx))=SPIKY_f_moving_average_weighted(fy(trac,1:length(fx)),fx,s_para.mao);
                end
                fy=fy(:,1:length(fx));
            elseif fdata_type==3      % sampled
                if size(fy,1)==1
                    fy=SPIKY_f_moving_average(fy,s_para.samp_mao);
                else
                    fy=SPIKY_f_moving_average_para(fy,s_para.samp_mao);
                end
            end
        elseif s_para.causal==1                                            % realtime (average over past values only)
            if fdata_type==1          % piecewise constant
                for trac=1:size(fy,1)
                    fy(trac,:)=SPIKY_f_moving_average_weighted_p(fy(trac,:),fx,s_para.mao);
                end
            elseif fdata_type==2      % piecewise linear
                mova=1;
                for trac=1:size(fy,1)
                    %fx=fx(2:2:end)-fx([1 2:2:end-2]);
                    fy(trac,1:length(fx))=mean([fy(trac,2:2:end); fy(trac,1:2:end)]);
                    fy(trac,1:length(fx))=SPIKY_f_moving_average_weighted_p(fy(trac,1:length(fx)),fx,s_para.mao);
                end
                fy=fy(:,1:length(fx));
            elseif fdata_type==3      % sampled
                if size(fy,1)==1
                    fy=SPIKY_f_moving_average_p(fy,s_para.samp_mao);
                else
                    fy=SPIKY_f_moving_average_para_p(fy,s_para.samp_mao);
                end
            end
        elseif s_para.causal==2                                            % future (average over future values only)
            if fdata_type==1          % piecewise constant
                for trac=1:size(fy,1)
                    fy(trac,:)=SPIKY_f_moving_average_weighted_f(fy(trac,:),fx,s_para.mao);
                end
            elseif fdata_type==2      % piecewise linear
                mova=1;
                for trac=1:size(fy,1)
                    %fx=fx(2:2:end)-fx([1 2:2:end-2]);
                    fy(trac,1:length(fx))=mean([fy(trac,2:2:end); fy(trac,1:2:end)]);
                    fy(trac,1:length(fx))=SPIKY_f_moving_average_weighted_f(fy(trac,1:length(fx)),fx,s_para.mao);
                end
                fy=fy(:,1:length(fx));
            elseif fdata_type==3      % sampled
                if size(fy,1)==1
                    fy=SPIKY_f_moving_average_f(fy,s_para.samp_mao);
                else
                    fy=SPIKY_f_moving_average_para_f(fy,s_para.samp_mao);
                end
            end
        end
        %fheadline=[fheadline,'*'];
    end

    if fdata_type==1 || mova==1                                            % using vectors of length "num_isi" (piecewise constant)
        ctfx=cumsum([0 fx]);
        pfx=s_para.itmin+sort([ctfx(1:end) ctfx(2:end-1)]);
        if size(fy,1)==1
            pfy=reshape([fy; fy],1,length(fy)*2);
        end
        ave=sum(fy.*repmat(fx,size(fy,1),1),2)/sum(fx);
    elseif fdata_type==2                                                   % using vectors of length "2*num_isi" (piecewise linear)
        ctfx=cumsum([0 fx]);
        pfx=s_para.itmin+sort([ctfx(1:end) ctfx(2:end-1)]);
        pfy=fy;
        if s_para.causal==0
            odds=1:2:length(fx)*2;
            evens=odds+1;
            ave=sum((fy(:,odds)+fy(:,evens))/2.*repmat(fx,size(fy,1),1),2)/sum(fx);
        else                                             % 1-realtime,2-future
            ave=profi_ave';
        end
    elseif fdata_type==3                                                   % using vectors of length "len" (sampled)
        pfx=fx;
        pfy=fy;
        ave=mean(fy,2); %*(s_para.psth+(s_para.psth==0));
        ave=ave*(s_para.psth+(s_para.psth==0));
    end
    
    if mod(s_para.plot_mode,2)>0 && fsp_posi>0                                                              % plotting
        if s_para.profile_norm_mode<4
            num_vals=2;
            if s_para.num_subplots<4
                intvals=[0.5 1];
            elseif s_para.num_subplots<8
                intvals=[0.4 0.8];
            else
                intvals=[0.3 0.6];
            end
            intlab=unique([0 SPIKY_f_lab(intvals*fmaxval,num_vals,1,1)]);
        else
            num_vals=3;
            if s_para.num_subplots<4
                intvals=[fminval fmaxval];
            elseif s_para.num_subplots<8
                intvals=[fminval+0.1*(fmaxval-fminval) fmaxval-0.1*(fmaxval-fminval)];
            else
                intvals=[fminval+0.25*(fmaxval-fminval) fmaxval-0.25*(fmaxval-fminval)];
            end
            intlab=unique(SPIKY_f_lab(intvals,num_vals,0,1));
            intlab=unique(intlab([1 end]));
        end
        if s_para.psth>0
            intlab(1)=0;
        end
        %[fminval fmaxval]
        yt=[fyt s_para.yl(2)-fsp_start+(0.05+(intlab-fminval)/(fmaxval-fminval))/1.1*fsp_size];
        intlab(end)=intlab(end)*(s_para.psth+(s_para.psth==0));
        ytl=[fytl intlab];
        
        text(s_para.xl(1)-0.08*(s_para.xl(2)-s_para.xl(1)),s_para.yl(2)-fsp_start+0.6/1.1*fsp_size,fheadline,...
            'Visible',p_para.measure_vis,'Color',p_para.measure_col,'FontSize',p_para.measure_fs,'FontWeight',p_para.measure_fw,...
            'FontAngle',p_para.measure_fa,'UIContextMenu',p_para.measure_cmenu);

        line(s_para.xl,s_para.yl(2)-fsp_start+0.05/1.1*fsp_size*ones(1,2),...
            'Visible',p_para.sp_bounds_vis,'Color',p_para.sp_bounds_col,'LineStyle',p_para.sp_bounds_ls,...
            'LineWidth',p_para.sp_bounds_lw,'UIContextMenu',p_para.sp_bounds_cmenu);
        line(s_para.xl,s_para.yl(2)-fsp_start+1.05/1.1*fsp_size*ones(1,2),...
            'Visible',p_para.sp_bounds_vis,'Color',p_para.sp_bounds_col,'LineStyle',p_para.sp_bounds_ls,...
            'LineWidth',p_para.sp_bounds_lw,'UIContextMenu',p_para.sp_bounds_cmenu);
        
        if size(fy,1)==1                                                                % one single profile
            
            if s_para.nma(nmac)==1
                if length(pfx)~=length(pfy)
                    p_para.prof_lh(p_para.supc,1)=plot(pfx,s_para.yl(2)-fsp_start+((0.05+pfy-fminval)/(fmaxval-fminval))/1.1*fsp_size,...  % '.-',
                        'Visible',p_para.prof_vis,'Color',p_para.prof_col,'LineStyle',p_para.prof_ls,...
                        'LineWidth',p_para.prof_lw,'UIContextMenu',p_para.prof_cmenu);
                else
                    p_para.prof_lh(p_para.supc,1)=plot(pfx,s_para.yl(2)-fsp_start+(0.05+(pfy-fminval)/(fmaxval-fminval))/1.1*fsp_size,...
                        'Visible',p_para.prof_vis,'Color',p_para.prof_col,'LineStyle',p_para.prof_ls,...
                        'LineWidth',p_para.prof_lw,'UIContextMenu',p_para.prof_cmenu);
                end
                set(p_para.prof_lh(p_para.supc),'UserData',[pfx mac; pfy 0],'Tag',fheadline)
            else
                p_para.ma_prof_lh(p_para.supc,1)=plot(pfx,s_para.yl(2)-fsp_start+(0.05+(pfy-fminval)/(fmaxval-fminval))/1.1*fsp_size,...
                    'Visible',p_para.ma_prof_vis,'Color',p_para.ma_prof_col,'LineStyle',p_para.ma_prof_ls,...
                    'LineWidth',p_para.ma_prof_lw,'UIContextMenu',p_para.ma_prof_cmenu);
                set(p_para.ma_prof_lh(p_para.supc),'UserData',[pfx mac; pfy 0],'Tag',fheadline)
            end
            if s_para.profile_average_line && nmac==1
                p_para.ave_lh(p_para.supc)=line([s_para.itmin s_para.itmax],s_para.yl(2)-fsp_start+(0.05+(ave-fminval)/(fmaxval-fminval))/1.1*fsp_size*ones(1,2),...
                    'Visible',p_para.ave_vis,'Color',p_para.ave_col,'LineStyle',p_para.ave_ls,...
                    'LineWidth',p_para.ave_lw,'UIContextMenu',p_para.ave_cmenu);
                p_para.prof_ave_fh(p_para.supc)=text(s_para.xl(2)-0.1*(s_para.xl(2)-s_para.xl(1)),s_para.yl(2)-fsp_start+0.93/1.1*fsp_size,num2str(ave,3),...
                    'Color',p_para.prof_ave_col,'FontSize',p_para.prof_ave_fs,'FontWeight',p_para.prof_ave_fw,'UIContextMenu',p_para.prof_ave_cmenu);
            end

        elseif nmac<2                                                                   % ISIs and rates for individual trials

            if fdata_type==1                                               % using vectors of length "num_isi" (piecewise constant)
                pfyy=zeros(size(fy));
                for trac=1:size(fy,1)
                    fyy=fy(trac,:);
                    pfyy(trac,1:length(fyy)*2)=reshape([fyy; fyy],1,length(fyy)*2);
                    if isfield(s_para,'dcols') && size(s_para.dcols,1)==size(fy,1)
                        if s_para.profile_mode==3 && trac==1
                            p_para.prof_lh(p_para.supc,trac)=plot(pfx,s_para.yl(2)-fsp_start+(0.05+(pfyy(trac,1:length(fyy)*2)-fminval)/(fmaxval-fminval))/1.1*fsp_size,...
                                'Visible',p_para.prof_vis,'Color',s_para.dcols(trac,:),'LineStyle',p_para.prof_ls,...
                                'LineWidth',p_para.prof_lw+1,'UIContextMenu',p_para.prof_cmenu);
                        else
                            p_para.prof_lh(p_para.supc,trac)=plot(pfx,s_para.yl(2)-fsp_start+(0.05+(pfyy(trac,1:length(fyy)*2)-fminval)/(fmaxval-fminval))/1.1*fsp_size,...
                                'Visible',p_para.prof_vis,'Color',s_para.dcols(trac,:),'LineStyle',p_para.prof_ls,...
                                'LineWidth',p_para.prof_lw,'UIContextMenu',p_para.prof_cmenu);
                        end
                    else
                        p_para.prof_lh(p_para.supc,trac)=plot(pfx,s_para.yl(2)-fsp_start+(0.05+(pfyy(trac,1:length(fyy)*2)-fminval)/(fmaxval-fminval))/1.1*fsp_size,...
                            'Visible',p_para.prof_vis,'Color',p_para.prof_col,'LineStyle',p_para.prof_ls,...
                            'LineWidth',p_para.prof_lw,'UIContextMenu',p_para.prof_cmenu);
                    end
                end
                set(p_para.prof_lh(p_para.supc),'UserData',[pfx mac; pfyy zeros(1,size(fy,1))'],'Tag',fheadline)
            else                                    % using vectors of length length "2*num_isi" (piecewise linear) or "len" (sampled) 
                for trac=1:size(fy,1)
                    pfyy=fy(trac,:);
                    if isfield(s_para,'dcols') && size(s_para.dcols,1)==size(fy,1)
                        if s_para.profile_mode==3 && trac==1
                            p_para.prof_lh(p_para.supc,trac)=plot(pfx,s_para.yl(2)-fsp_start+(0.05+(pfyy-fminval)/(fmaxval-fminval))/1.1*fsp_size,...
                                'Visible',p_para.prof_vis,'Color',s_para.dcols(trac,:),'LineStyle',p_para.prof_ls,...
                                'LineWidth',p_para.prof_lw+1,'UIContextMenu',p_para.prof_cmenu);
                        else
                            p_para.prof_lh(p_para.supc,trac)=plot(pfx,s_para.yl(2)-fsp_start+(0.05+(pfyy-fminval)/(fmaxval-fminval))/1.1*fsp_size,...
                                'Visible',p_para.prof_vis,'Color',s_para.dcols(trac,:),'LineStyle',p_para.prof_ls,...
                                'LineWidth',p_para.prof_lw,'UIContextMenu',p_para.prof_cmenu);
                        end
                    else
                        p_para.prof_lh(p_para.supc,trac)=plot(pfx,s_para.yl(2)-fsp_start+(0.05+(pfyy-fminval)/(fmaxval-fminval))/1.1*fsp_size,...
                            'Visible',p_para.prof_vis,'Color',p_para.prof_col,'LineStyle',p_para.prof_ls,...
                            'LineWidth',p_para.prof_lw,'UIContextMenu',p_para.prof_cmenu);
                    end
                    if s_para.profile_average_line && trac==1
                        p_para.ave_lh(p_para.supc)=line([s_para.itmin s_para.itmax],s_para.yl(2)-fsp_start+(0.05+(ave(1)-fminval)/(fmaxval-fminval))/1.1*fsp_size*ones(1,2),...
                            'Visible',p_para.ave_vis,'Color',p_para.ave_col,'LineStyle',p_para.ave_ls,...
                            'LineWidth',p_para.ave_lw,'UIContextMenu',p_para.ave_cmenu);
                        p_para.prof_ave_fh(p_para.supc)=text(s_para.xl(2)-0.1*(s_para.xl(2)-s_para.xl(1)),s_para.yl(2)-fsp_start+0.93/1.1*fsp_size,num2str(ave(1),3),...
                            'Color',p_para.prof_ave_col,'FontSize',p_para.prof_ave_fs,'FontWeight',p_para.prof_ave_fw,'UIContextMenu',p_para.prof_ave_cmenu);
                    end
                end
                set(p_para.prof_lh(p_para.supc),'UserData',[pfx mac; fy zeros(1,size(fy,1))'],'Tag',fheadline)
            end
        end
    else
        yt=[];
        ytl=[];
    end
end
%if size(fy,1)==1
    fave=[fave ave];
%else
end

