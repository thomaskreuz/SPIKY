% This function provides an input mask (for keyboard and mouse) for selecting temporal intervals for the movie (their names and their divisions).

function SPIKY_select_intervals(hObject, eventdata, handles)

set(handles.figure1,'Visible','off')

d_para=getappdata(handles.figure1,'data_parameters');
f_para=getappdata(handles.figure1,'figure_parameters');
p_para=getappdata(handles.figure1,'plot_parameters');
s_para=getappdata(handles.figure1,'subplot_parameters');

SG_fig=figure('units','normalized','menubar','none','position',[0.05 0.17 0.4 0.66],'Name','Select group separators',...
    'NumberTitle','off','Color',[0.9294 0.9176 0.851],'DeleteFcn',{@SG_Close_callback}); % ,'WindowStyle','modal'

SID_panel=uipanel('units','normalized','position',[0.03 0.66 0.94 0.3],'Title','Group separators','FontSize',15,...
    'FontWeight','bold','HighlightColor','k');
SID_edit=uicontrol('style','edit','units','normalized','position',[0.06 0.845 0.88 0.05],...
    'FontSize',15,'FontUnits','normalized','BackgroundColor','w');
SID_Cancel_pushbutton=uicontrol('style','pushbutton','units','normalized','position',[0.1 0.765 0.22 0.05],'string','Cancel',...
    'FontSize',15,'FontUnits','normalized','BackgroundColor',[0.8353 0.8235 0.7922],'CallBack',{@SID_callback});
SID_Reset_pushbutton=uicontrol('style','pushbutton','units','normalized','position',[0.39 0.765 0.22 0.05],'string','Reset',...
    'FontSize',15,'FontUnits','normalized','BackgroundColor',[0.8353 0.8235 0.7922],'CallBack',{@SID_callback});
SID_Apply_pushbutton=uicontrol('style','pushbutton','units','normalized','position',[0.68 0.765 0.22 0.05],'string','Apply',...
    'FontSize',15,'FontUnits','normalized','BackgroundColor',[0.8353 0.8235 0.7922],'CallBack',{@SID_callback});
SID_Load_pushbutton=uicontrol('style','pushbutton','units','normalized','position',[0.19 0.69 0.28 0.05],'string','Load from file',...
    'FontSize',15,'FontUnits','normalized','FontWeight','bold','BackgroundColor',[0.8353 0.8235 0.7922],'CallBack',{@SID_Load_interval_divisions});
SID_OK_pushbutton=uicontrol('style','pushbutton','units','normalized','position',[0.53 0.69 0.28 0.05],'string','OK',...
    'FontSize',15,'FontUnits','normalized','FontWeight','bold','BackgroundColor',[0.8353 0.8235 0.7922],'CallBack',{@SID_callback});
uicontrol(SID_OK_pushbutton)

SIN_panel=uipanel('units','normalized','position',[0.03 0.05 0.94 0.58],'Title','Group names','FontSize',15,...
    'FontWeight','bold','Visible','off');
SIN_edit=uicontrol('style','edit','units','normalized','position',[0.06 0.51 0.88 0.05],...
    'FontSize',15,'FontUnits','normalized','BackgroundColor','w','Visible','off');
SIN_Down_pushbutton=uicontrol('style','pushbutton','units','normalized','position',[0.1 0.425 0.22 0.05],...
    'string','Down','FontSize',15,'FontUnits','normalized','BackgroundColor',[0.8353 0.8235 0.7922],...
    'CallBack',{@SIN_ListBox_callback},'Visible','off');
SIN_Up_pushbutton=uicontrol('style','pushbutton','units','normalized','position',[0.39 0.425 0.22 0.05],...
    'string','Up','FontSize',15,'FontUnits','normalized','BackgroundColor',[0.8353 0.8235 0.7922],...
    'CallBack',{@SIN_ListBox_callback},'Visible','off');
SIN_Replace_pushbutton=uicontrol('style','pushbutton','units','normalized','position',[0.68 0.425 0.22 0.05],...
    'string','Replace','FontSize',15,'FontUnits','normalized','BackgroundColor',[0.8353 0.8235 0.7922],...
    'CallBack',{@SIN_ListBox_callback},'Visible','off');
SIN_groups_listbox=uicontrol('style','listbox','units','normalized','position',[0.08 0.17 0.4 0.22],...
    'FontSize',15,'FontUnits','normalized','BackgroundColor','w','min',0,'max',1,'Visible','off','Enable','off');
SIN_names_listbox=uicontrol('style','listbox','units','normalized','position',[0.52 0.17 0.4 0.22],...
    'FontSize',15,'FontUnits','normalized','BackgroundColor','w','CallBack',{@SIN_ListBox_callback},'min',0,'max',1,'Visible','off');
SIN_Load_pushbutton=uicontrol('style','pushbutton','units','normalized','position',[0.06 0.085 0.22 0.05],...
    'string','Load from file','FontSize',15,'FontUnits','normalized','BackgroundColor',[0.8353 0.8235 0.7922],...
    'CallBack',{@SIN_Load_group_names},'Visible','off');
SIN_Cancel_pushbutton=uicontrol('style','pushbutton','units','normalized','position',[0.3 0.085 0.2 0.05],...
    'string','Cancel','FontSize',15,'FontUnits','normalized','BackgroundColor',[0.8353 0.8235 0.7922],...
    'CallBack',{@SIN_callback},'Visible','off');
SIN_Reset_pushbutton=uicontrol('style','pushbutton','units','normalized','position',[0.52 0.085 0.2 0.05],...
    'string','Reset','FontSize',15,'FontUnits','normalized','BackgroundColor',[0.8353 0.8235 0.7922],...
    'CallBack',{@SIN_callback},'Visible','off');
SIN_OK_pushbutton=uicontrol('style','pushbutton','units','normalized','position',[0.74 0.085 0.2 0.05],...
    'string','OK','FontSize',15,'FontUnits','normalized','FontWeight','bold','BackgroundColor',[0.8353 0.8235 0.7922],...
    'CallBack',{@SIN_callback},'Visible','off');
uicontrol(SIN_OK_pushbutton)

figure(f_para.num_fig);
xlim([s_para.itmin-0.04*s_para.itrange s_para.itmax+0.04*s_para.itrange])
SID_UserData.fh=gcf;
SID_UserData.ah=gca;
set(SID_UserData.fh,'Units','Normalized')
SID_UserData.num_trains=d_para.num_trains;
SID_UserData.xlim=xlim;
SID_UserData.tmin=d_para.tmin;
SID_UserData.tmax=d_para.tmax;
SID_UserData.dts=d_para.dts;
SID_UserData.flag=0;
SID_UserData.marked_indy=[];

if isfield(d_para,'all_train_group_sizes')
    SID_UserData.group_sizes=d_para.all_train_group_sizes;
    SID_UserData.interval_divisions=cumsum(SID_UserData.group_sizes);
    SID_UserData.interval_divisions=unique(round(SID_UserData.interval_divisions(SID_UserData.interval_divisions>0 & SID_UserData.interval_divisions<SID_UserData.num_trains)));
    set(SID_edit,'String',regexprep(num2str(SID_UserData.interval_divisions),'\s+',' '))
    SID_UserData.lh=zeros(1,length(SID_UserData.interval_divisions));
    for lhc=1:length(SID_UserData.interval_divisions)
        SID_UserData.lh(lhc)=line([SID_UserData.tmin SID_UserData.tmax],(0.05+((SID_UserData.num_trains-SID_UserData.interval_divisions(lhc))/SID_UserData.num_trains))*ones(1,2),...
            'Color',p_para.SID_col,'LineStyle',p_para.SID_ls,'LineWidth',p_para.SID_lw);
    end
else
    SID_UserData.interval_divisions=[];
end
SID_UserData.tx=uicontrol('style','tex','String','','unit','normalized','backg',get(SID_UserData.fh,'Color'),...
    'position',[0.36 0.907 0.4 0.04],'FontSize',18,'FontWeight','bold','HorizontalAlignment','left');
SID_UserData.spike_lh=getappdata(handles.figure1,'spike_lh');

set(SID_UserData.fh,'Userdata',SID_UserData)
SID_UserData.cm=uicontextmenu;
SID_UserData.um=uimenu(SID_UserData.cm,'label','Delete group separator(s)','CallBack',{@SID_delete_interval_divisions,SID_UserData});
set(SID_UserData.fh,'WindowButtonMotionFcn',{@SID_get_coordinates,SID_UserData},'KeyPressFcn',{@SID_keyboard,SID_UserData});
set(SID_UserData.lh,'ButtonDownFcn',{@SID_start_move_interval_divisions,SID_UserData},'UIContextMenu',SID_UserData.cm)
set(SID_UserData.ah,'ButtonDownFcn',{@SID_pick_interval_divisions,SID_UserData},'UIContextMenu',[])
set(SID_UserData.spike_lh,'ButtonDownFcn',{@SID_pick_interval_divisions,SID_UserData},'UIContextMenu',[])
set(SID_UserData.fh,'Userdata',SID_UserData)

    function SID_get_coordinates(varargin)
        SID_UserData=varargin{3};
        SID_UserData=get(SID_UserData.fh,'UserData');

        ax_pos=get(SID_UserData.ah,'CurrentPoint');
        ax_x_ok=SID_UserData.tmin<=ax_pos(1,1)&& ax_pos(1,1)<=SID_UserData.tmax;
        ax_y_ok=0.05<=ax_pos(1,2) && ax_pos(1,2)<=1.05;

        if ax_x_ok && ax_y_ok
            ax_y=SID_UserData.num_trains-round((ax_pos(1,2)-0.05)*SID_UserData.num_trains);
            ax_y(ax_y==0 && ax_y_ok)=1;
            ax_y(ax_y==SID_UserData.num_trains && ax_y_ok)=SID_UserData.num_trains-1;            
            set(SID_UserData.tx,'str',['After spike train: ',num2str(ax_y)]);
        else
            set(SID_UserData.tx,'str','Out of range');
        end
    end


    function SID_pick_interval_divisions(varargin)                                   % SID_marked_indy changes
        SID_UserData=varargin{3};
        SID_UserData=get(SID_UserData.fh,'UserData');

        ax_pos=get(SID_UserData.ah,'CurrentPoint');
        ax_x_ok=SID_UserData.tmin<=ax_pos(1,1)&& ax_pos(1,1)<=SID_UserData.tmax;
        ax_y_ok=0.05<=ax_pos(1,2) && ax_pos(1,2)<=1.05;
        ax_y=SID_UserData.num_trains-round((ax_pos(1,2)-0.05)*SID_UserData.num_trains);
        ax_y(ax_y==0 && ax_y_ok)=1;
        ax_y(ax_y==SID_UserData.num_trains && ax_y_ok)=SID_UserData.num_trains-1;

        modifiers=get(SID_UserData.fh,'CurrentModifier');
        if ax_x_ok && ax_y_ok
            seltype=get(SID_UserData.fh,'SelectionType');            % Right-or-left click?
            ctrlIsPressed=ismember('control',modifiers);
            if strcmp(seltype,'normal') || ctrlIsPressed                      % right or middle
                if ~isempty(SID_UserData.interval_divisions)
                    num_lh=length(SID_UserData.interval_divisions);
                    for lhc=1:num_lh
                        if ishandle(SID_UserData.lh(lhc))
                            delete(SID_UserData.lh(lhc))
                        end
                    end
                    SID_UserData.interval_divisions=unique([SID_UserData.interval_divisions ax_y]);
                else
                    SID_UserData.interval_divisions=ax_y;
                end
                num_lh=length(SID_UserData.interval_divisions);
                SID_UserData.lh=zeros(1,num_lh);
                for selc=1:num_lh
                    if ismember(selc,SID_UserData.marked_indy)
                        SID_UserData.lh(selc)=line([SID_UserData.tmin SID_UserData.tmax],(0.05+((SID_UserData.num_trains-SID_UserData.interval_divisions(selc))/SID_UserData.num_trains))*ones(1,2),...
                            'Visible',p_para.SID_vis,'Color',p_para.SID_marked_col,'LineStyle',p_para.SID_ls,'LineWidth',p_para.SID_lw);
                    else
                        SID_UserData.lh(selc)=line([SID_UserData.tmin SID_UserData.tmax],(0.05+((SID_UserData.num_trains-SID_UserData.interval_divisions(selc))/SID_UserData.num_trains))*ones(1,2),...
                            'Visible',p_para.SID_vis,'Color',p_para.SID_col,'LineStyle',p_para.SID_ls,'LineWidth',p_para.SID_lw);
                    end
                end

                SID_str=num2str(SID_UserData.interval_divisions);
                if length(SID_UserData.interval_divisions)>1
                    SID_str=regexprep(SID_str,'\s+',' ');
                end
                set(SID_edit,'String',SID_str)
            end

            shftIsPressed=ismember('shift',modifiers);
            if ctrlIsPressed
                if isfield(SID_UserData,'marked')
                    SID_UserData.marked=unique([SID_UserData.marked ax_y]);
                else
                    SID_UserData.marked=ax_y;
                end
                [dummy,this,dummy2]=intersect(SID_UserData.interval_divisions,SID_UserData.marked);
                SID_UserData.marked_indy=this;
                set(SID_UserData.lh,'Color',p_para.SID_col)
                set(SID_UserData.lh(SID_UserData.marked_indy),'Color',p_para.SID_marked_col)
                SID_UserData.flag=1;
            elseif ~shftIsPressed
                set(SID_UserData.lh,'Color',p_para.SID_col)
                SID_UserData.marked=[];
                SID_UserData.flag=0;
                SID_UserData.marked_indy=[];
            end
        end

        shftIsPressed=ismember('shift',modifiers);
        dummy_marked=SID_UserData.marked_indy;
        if shftIsPressed
            set(SID_UserData.fh,'WindowButtonUpFcn',[])
            ax_pos=get(SID_UserData.ah,'CurrentPoint');
            first_corner_x=ax_pos(1,1);
            first_corner_y=ax_pos(1,2);
            window=rectangle('Position',[first_corner_x, first_corner_y, 0.01, 0.01]);
            group_seps=(0.05+((SID_UserData.num_trains-SID_UserData.interval_divisions)/SID_UserData.num_trains));
            while shftIsPressed
                ax_pos=get(SID_UserData.ah,'CurrentPoint');
                second_corner_x=ax_pos(1,1);
                second_corner_y=ax_pos(1,2);
                if second_corner_x ~= first_corner_x && second_corner_y ~= first_corner_y
                    set(window,'Position',[min(first_corner_x, second_corner_x), min(first_corner_y, second_corner_y), abs(second_corner_x-first_corner_x), abs(second_corner_y-first_corner_y)])
                    drawnow
                    bottom_mark=min(first_corner_y, second_corner_y);
                    top_mark=max(first_corner_y, second_corner_y);
                    SID_UserData.marked_indy=unique([dummy_marked find(group_seps>=bottom_mark & group_seps<=top_mark)]);
                    SID_UserData.flag=(~isempty(SID_UserData.marked_indy));
                    set(SID_UserData.lh,'Color',p_para.SID_col)
                    set(SID_UserData.lh(SID_UserData.marked_indy),'Color',p_para.SID_marked_col)
                end
                pause(0.001);
                modifiers=get(SID_UserData.fh,'CurrentModifier');
                shftIsPressed=ismember('shift',modifiers);
            end
            delete(window)
            SID_UserData.marked=SID_UserData.interval_divisions(SID_UserData.marked_indy);
        end
        
        %pick_marky=SID_UserData.marked
        %pick_marky_indy=SID_UserData.marked_indy

        set(SID_UserData.fh,'Userdata',SID_UserData)
        set(SID_UserData.um,'CallBack',{@SID_delete_interval_divisions,SID_UserData})
        set(SID_UserData.fh,'WindowButtonMotionFcn',{@SID_get_coordinates,SID_UserData},'KeyPressFcn',{@SID_keyboard,SID_UserData})
        set(SID_UserData.lh,'ButtonDownFcn',{@SID_start_move_interval_divisions,SID_UserData},'UIContextMenu',SID_UserData.cm)
        set(SID_UserData.ah,'ButtonDownFcn',{@SID_pick_interval_divisions,SID_UserData},'UIContextMenu',[])
        set(SID_UserData.spike_lh,'ButtonDownFcn',{@SID_pick_interval_divisions,SID_UserData},'UIContextMenu',[])
        set(SID_UserData.fh,'Userdata',SID_UserData)
    end


    function SID_delete_interval_divisions(varargin)
        SID_UserData=varargin{3};
        SID_UserData=get(SID_UserData.fh,'UserData');
        
        %delete_marky=SID_UserData.marked
        %delete_marky_indy=SID_UserData.marked_indy

        if (SID_UserData.flag)
            SID_UserData.interval_divisions=setdiff(SID_UserData.interval_divisions,SID_UserData.interval_divisions(SID_UserData.marked_indy));
            SID_UserData.marked_indy=[];
            SID_UserData.flag=0;
        else
            ax_y=get(gco,'YData');
            SID_UserData.interval_divisions=setdiff(SID_UserData.interval_divisions,SID_UserData.num_trains-ceil((ax_y(1)-0.05)*SID_UserData.num_trains));
        end
        for lhc=1:length(SID_UserData.lh)
            if ishandle(SID_UserData.lh(lhc))
                delete(SID_UserData.lh(lhc))
            end
        end
        num_lh=length(SID_UserData.interval_divisions);
        SID_UserData.lh=zeros(1,num_lh);
        for selc=1:num_lh
            SID_UserData.lh(selc)=line([SID_UserData.tmin SID_UserData.tmax],(0.05+((SID_UserData.num_trains-SID_UserData.interval_divisions(selc))/SID_UserData.num_trains))*ones(1,2),...
                'Visible',p_para.SID_vis,'Color',p_para.SID_col,'LineStyle',p_para.SID_ls,'LineWidth',p_para.SID_lw);
        end
        SID_str=num2str(SID_UserData.interval_divisions);
        if num_lh>1
            SID_str=regexprep(SID_str,'\s+',' ');
        end
        set(SID_edit,'String',SID_str)

        set(SID_UserData.fh,'Userdata',SID_UserData)
        set(SID_UserData.um,'CallBack',{@SID_delete_interval_divisions,SID_UserData})
        set(SID_UserData.fh,'WindowButtonMotionFcn',{@SID_get_coordinates,SID_UserData},'KeyPressFcn',{@SID_keyboard,SID_UserData})
        set(SID_UserData.lh,'ButtonDownFcn',{@SID_start_move_interval_divisions,SID_UserData},'UIContextMenu',SID_UserData.cm)
        set(SID_UserData.ah,'ButtonDownFcn',{@SID_pick_interval_divisions,SID_UserData},'UIContextMenu',[])
        set(SID_UserData.spike_lh,'ButtonDownFcn',{@SID_pick_interval_divisions,SID_UserData},'UIContextMenu',[])
        set(SID_UserData.fh,'Userdata',SID_UserData)
    end


    function SID_start_move_interval_divisions(varargin)
        SID_UserData=varargin{3};
        SID_UserData=get(SID_UserData.fh,'UserData');

        seltype=get(SID_UserData.fh,'SelectionType'); % Right-or-left click?
        %ax_pos=get(SID_UserData.ah,'CurrentPoint');
        if strcmp(seltype,'alt')
            modifiers=get(SID_UserData.fh,'CurrentModifier');
            ctrlIsPressed=ismember('control',modifiers);
            if ctrlIsPressed
                ax_y=get(gco,'YData');

                if isfield(SID_UserData,'marked')
                    SID_UserData.marked=unique([SID_UserData.marked ax_y]);
                else
                    SID_UserData.marked=ax_y;
                end
                [dummy,this,dummy2]=intersect(SID_UserData.interval_divisions,SID_UserData.marked);
                SID_UserData.marked_indy=this;
                if ~SID_UserData.flag
                    SID_UserData.flag=1;
                end
                set(gco,'Color',p_para.SID_marked_col);

                set(SID_UserData.fh,'Userdata',SID_UserData)
                set(SID_UserData.um,'CallBack',{@SID_delete_interval_divisions,SID_UserData})
                set(SID_UserData.fh,'WindowButtonMotionFcn',{@SID_get_coordinates,SID_UserData},'KeyPressFcn',{@SID_keyboard,SID_UserData})
                set(SID_UserData.lh,'ButtonDownFcn',{@SID_start_move_interval_divisions,SID_UserData},'UIContextMenu',SID_UserData.cm)
                set(SID_UserData.ah,'ButtonDownFcn',{@SID_pick_interval_divisions,SID_UserData},'UIContextMenu',[])
                set(SID_UserData.spike_lh,'ButtonDownFcn',{@SID_pick_interval_divisions,SID_UserData},'UIContextMenu',[])
                set(SID_UserData.fh,'Userdata',SID_UserData)
            end
        else
            if ~SID_UserData.flag
                num_lh=length(SID_UserData.interval_divisions);
                for selc=1:num_lh
                    set(SID_UserData.lh(selc),'Color',p_para.SID_col);
                end
                SID_UserData.marked_indy=find(SID_UserData.lh(:) == gco);
                dummy=get(gco,'YData');
                SID_UserData.initial_ST=SID_UserData.num_trains-floor((dummy(1)-0.05)*SID_UserData.num_trains);
                SID_UserData.initial_YPos=SID_UserData.initial_ST;
            else
                %SID_UserData.initial_YPos=SID_UserData.num_trains-floor((ax_pos(1,2)-0.05)*SID_UserData.num_trains);
                SID_UserData.initial_ST=SID_UserData.interval_divisions(SID_UserData.marked_indy);
                SID_UserData.initial_YPos=SID_UserData.initial_ST(1);
            end
            set(SID_UserData.fh,'WindowButtonMotionFcn',{@SID_move_interval_divisions,SID_UserData})
            set(SID_UserData.ah,'ButtonDownFcn','')
            set(SID_UserData.spike_lh,'ButtonDownFcn','')
            set(SID_UserData.fh,'Userdata',SID_UserData)
        end        
    end


    function SID_move_interval_divisions(varargin)
        SID_UserData=varargin{3};
        SID_UserData=get(SID_UserData.fh,'UserData');

        ax_pos=get(SID_UserData.ah,'CurrentPoint');
        ax_x_ok=SID_UserData.tmin<=ax_pos(1,1)&& ax_pos(1,1)<=SID_UserData.tmax;
        ax_y_ok=0.05<=ax_pos(1,2) && ax_pos(1,2)<=1.05;
        
        ax_y=SID_UserData.num_trains-round((ax_pos(1,2)-0.05)*SID_UserData.num_trains);
        ax_y(ax_y==0 && ax_y_ok)=1;
        ax_y(ax_y==SID_UserData.num_trains && ax_y_ok)=SID_UserData.num_trains-1;
        
        for idx=1:length(SID_UserData.marked_indy)
            if SID_UserData.initial_ST(idx)+ax_y-SID_UserData.initial_YPos>0 && ...
                    SID_UserData.initial_ST(idx)+ax_y-SID_UserData.initial_YPos<SID_UserData.num_trains
                set(SID_UserData.lh(SID_UserData.marked_indy(idx)),'Color',p_para.SID_marked_col,...
                    'YData',(0.05+((SID_UserData.num_trains-(SID_UserData.initial_ST(idx)+ax_y-SID_UserData.initial_YPos))/...
                    SID_UserData.num_trains))*ones(1,2),'LineWidth',2)
            else
                set(SID_UserData.lh(SID_UserData.marked_indy(idx)),'YData',(0.05+((SID_UserData.num_trains-(SID_UserData.initial_ST(idx)+ax_y-SID_UserData.initial_YPos))/...
                    SID_UserData.num_trains))*ones(1,2),'Color','w')
            end
        end
        drawnow;

        if ax_x_ok && ax_y_ok
            %set(SID_UserData.tx,'str',['After spike train: ',num2str(ax_y)]);
            ax_y(ax_y==0 && ax_y_ok)=1;
            ax_y(ax_y==SID_UserData.num_trains && ax_y_ok)=SID_UserData.num_trains-1;
            set(SID_UserData.tx,'str',['After spike train: ',num2str(ax_y)]);
        else
            set(SID_UserData.tx,'str','Out of range');
        end
        set(SID_UserData.fh,'WindowButtonUpFcn',{@SID_stop_move_interval_divisions,SID_UserData})
        set(SID_UserData.fh,'Userdata',SID_UserData)
    end


    function SID_stop_move_interval_divisions(varargin)
        SID_UserData=varargin{3};
        SID_UserData=get(SID_UserData.fh,'UserData');

        ax_pos=get(SID_UserData.ah,'CurrentPoint');
        ax_y_ok=0.05<=ax_pos(1,2) && ax_pos(1,2)<=1.05;

        ax_y=SID_UserData.num_trains-round((ax_pos(1,2)-0.05)*SID_UserData.num_trains);
        ax_y(ax_y==0 && ax_y_ok)=1;
        ax_y(ax_y==SID_UserData.num_trains && ax_y_ok)=SID_UserData.num_trains-1;

        for idx=1:length(SID_UserData.marked_indy)
            if SID_UserData.initial_ST(idx)+ax_y-SID_UserData.initial_YPos>0 && ...
                    SID_UserData.initial_ST(idx)+ax_y-SID_UserData.initial_YPos<SID_UserData.num_trains
                SID_UserData.interval_divisions(SID_UserData.marked_indy(idx))=SID_UserData.initial_ST(idx)+ax_y-SID_UserData.initial_YPos;
            else
                SID_UserData.interval_divisions=SID_UserData.interval_divisions(setdiff(1:length(SID_UserData.interval_divisions),SID_UserData.marked_indy(idx)));
                SID_UserData.marked_indy(idx+1:end)=SID_UserData.marked_indy(idx+1:end)-1;
            end
        end
        SID_UserData.interval_divisions=unique(SID_UserData.interval_divisions);
        SID_UserData.marked_indy=[];
        SID_UserData.flag=0;
        for lhc=1:size(SID_UserData.lh,2)
            if ishandle(SID_UserData.lh(lhc))
                delete(SID_UserData.lh(lhc))
            end
        end
        num_lh=length(SID_UserData.interval_divisions);
        SID_UserData.lh=zeros(1,num_lh);
        for selc=1:num_lh
            SID_UserData.lh(selc)=line([SID_UserData.tmin SID_UserData.tmax],(0.05+((SID_UserData.num_trains-SID_UserData.interval_divisions(selc))/SID_UserData.num_trains))*ones(1,2),...
                'Visible',p_para.SID_vis,'Color',p_para.SID_col,'LineStyle',p_para.SID_ls,'LineWidth',p_para.SID_lw);
        end
        SID_str=num2str(SID_UserData.interval_divisions);
        if num_lh>1
            SID_str=regexprep(SID_str,'\s+',' ');
        end
        set(SID_edit,'String',SID_str)
        set(SID_UserData.fh,'Pointer','arrow')
        drawnow;

        set(SID_UserData.fh,'Userdata',SID_UserData)
        set(SID_UserData.um,'CallBack',{@SID_delete_interval_divisions,SID_UserData})
        set(SID_UserData.fh,'WindowButtonMotionFcn',{@SID_get_coordinates,SID_UserData},'KeyPressFcn',{@SID_keyboard,SID_UserData},...
            'WindowButtonUpFcn',[])
        set(SID_UserData.lh,'ButtonDownFcn',{@SID_start_move_interval_divisions,SID_UserData},'UIContextMenu',SID_UserData.cm)
        set(SID_UserData.ah,'ButtonDownFcn',{@SID_pick_interval_divisions,SID_UserData},'UIContextMenu',[])
        set(SID_UserData.spike_lh,'ButtonDownFcn',{@SID_pick_interval_divisions,SID_UserData},'UIContextMenu',[])
        set(SID_UserData.fh,'Userdata',SID_UserData)
    end


    function SID_keyboard(varargin)
        SID_UserData=varargin{3};
        SID_UserData=get(SID_UserData.fh,'UserData');

        if strcmp(varargin{2}.Key,'delete')
            if (SID_UserData.flag)
                SID_UserData.interval_divisions=setdiff(SID_UserData.interval_divisions,SID_UserData.interval_divisions(SID_UserData.marked_indy));
                SID_UserData.marked_indy=[];
                SID_UserData.flag=0;
                for lhc=1:length(SID_UserData.lh)
                    if ishandle(SID_UserData.lh(lhc))
                        delete(SID_UserData.lh(lhc))
                    end
                end
                num_lh=length(SID_UserData.interval_divisions);
                SID_UserData.lh=zeros(1,num_lh);
                for selc=1:num_lh
                    SID_UserData.lh(selc)=line([SID_UserData.tmin SID_UserData.tmax],(0.05+((SID_UserData.num_trains-SID_UserData.interval_divisions(selc))/SID_UserData.num_trains))*ones(1,2),...
                        'Visible',p_para.SID_vis,'Color',p_para.SID_col,'LineStyle',p_para.SID_ls,'LineWidth',p_para.SID_lw);
                end
                SID_str=num2str(SID_UserData.interval_divisions);
                if num_lh>1
                    SID_str=regexprep(SID_str,'\s+',' ');
                end
                set(SID_edit,'String',SID_str)

                set(SID_UserData.fh,'Userdata',SID_UserData)
                set(SID_UserData.um,'CallBack',{@SID_delete_interval_divisions,SID_UserData})
                set(SID_UserData.fh,'WindowButtonMotionFcn',{@SID_get_coordinates,SID_UserData},'KeyPressFcn',{@SID_keyboard,SID_UserData})
                set(SID_UserData.lh,'ButtonDownFcn',{@SID_start_move_interval_divisions,SID_UserData},'UIContextMenu',SID_UserData.cm)
                set(SID_UserData.ah,'ButtonDownFcn',{@SID_pick_interval_divisions,SID_UserData},'UIContextMenu',[])
                set(SID_UserData.spike_lh,'ButtonDownFcn',{@SID_pick_interval_divisions,SID_UserData},'UIContextMenu',[])
                set(SID_UserData.fh,'Userdata',SID_UserData)
            end
        elseif strcmp(varargin{2}.Key,'a')
            modifiers=get(SID_UserData.fh,'CurrentModifier');
            ctrlIsPressed=ismember('control',modifiers);
            if ctrlIsPressed
                num_lh=length(SID_UserData.interval_divisions);
                SID_UserData.marked=SID_UserData.interval_divisions;
                SID_UserData.marked_indy=1:num_lh;
                for selc=1:num_lh
                    set(SID_UserData.lh(selc),'Color',p_para.SID_marked_col);
                end
                SID_UserData.flag=1;

                set(SID_UserData.fh,'Userdata',SID_UserData)
                set(SID_UserData.um,'CallBack',{@SID_delete_interval_divisions,SID_UserData})
                set(SID_UserData.fh,'WindowButtonMotionFcn',{@SID_get_coordinates,SID_UserData},'KeyPressFcn',{@SID_keyboard,SID_UserData})
                set(SID_UserData.lh,'ButtonDownFcn',{@SID_start_move_interval_divisions,SID_UserData},'UIContextMenu',SID_UserData.cm)
                set(SID_UserData.ah,'ButtonDownFcn',{@SID_start_pick,SID_UserData})
                set(SID_UserData.spike_lh,'ButtonDownFcn',{@SID_pick_interval_divisions,SID_UserData},'UIContextMenu',[])
                set(SID_UserData.fh,'Userdata',SID_UserData)
            end
        elseif strcmp(varargin{2}.Key,'c')
            modifiers=get(SID_UserData.fh,'CurrentModifier');
            ctrlIsPressed=ismember('control',modifiers);
            if ctrlIsPressed
                num_lh=length(SID_UserData.interval_divisions);
                ax_y=SID_UserData.interval_divisions(SID_UserData.marked_indy);
                for idx=1:length(SID_UserData.marked_indy)
                    SID_UserData.lh(num_lh+idx)=line([SID_UserData.tmin SID_UserData.tmax],(0.05+((SID_UserData.num_trains-SID_UserData.interval_divisions(SID_UserData.marked_indy(idx)))/SID_UserData.num_trains))*ones(1,2),...
                        'Visible',p_para.SID_vis,'Color',p_para.SID_marked_col,'LineStyle',p_para.SID_ls,'LineWidth',p_para.SID_lw);
                end
                SID_UserData.marked=SID_UserData.interval_divisions(SID_UserData.marked_indy);
                SID_UserData.marked_indy=num_lh+(1:length(SID_UserData.marked_indy));
                SID_UserData.interval_divisions(SID_UserData.marked_indy)=SID_UserData.marked;
                SID_UserData.flag=1;
                ax_pos=get(SID_UserData.ah,'CurrentPoint');
                SID_UserData.initial_YPos=SID_UserData.num_trains-floor((ax_pos(1,2)-0.05)*SID_UserData.num_trains);
                SID_UserData.initial_ST=ax_y;
                set(SID_UserData.ah,'ButtonDownFcn','')
                set(SID_UserData.spike_lh,'ButtonDownFcn','')
                set(SID_UserData.fh,'WindowButtonMotionFcn',{@SID_move_interval_divisions,SID_UserData})
                set(SID_UserData.fh,'Userdata',SID_UserData)
            end
        end
    end


    function SID_Load_interval_divisions(varargin)
        SID_UserData=get(f_para.num_fig,'UserData');
        if ~isfield(d_para,'matfile')
            d_para.matfile=uigetfile('*.mat','Pick a .mat-file');
        end
        data=load(d_para.matfile);
        if isstruct(data)
            data.default_variable='';
            data.content='group separators';
            SPIKY('SPIKY_select_variable',gcbo,data,handles)
            variable=getappdata(handles.figure1,'variable');
            if ~isempty(variable) && ~strcmp(variable,'A9ZB1YC8X')
                SID_UserData.interval_divisions=squeeze(data.(variable));
            end
        elseif (iscell(data) && length(data)>1 && isnumeric(data{1})) || ...
                (isnumeric(data) && ndims(data)==2 && size(data,1)>1)
            SID_UserData.interval_divisions=squeeze(data);
        end
        if ~isnumeric(SID_UserData.interval_divisions) || ndims(SID_UserData.interval_divisions)~=2 || ~any(size(SID_UserData.interval_divisions)==1)
            if ~strcmp(variable,'A9ZB1YC8X')
                set(0,'DefaultUIControlFontSize',16);
                mbh=msgbox(sprintf('No vector has been chosen. Please try again!'),'Warning','warn','modal');
                htxt = findobj(mbh,'Type','text');
                set(htxt,'FontSize',12,'FontWeight','bold')
                mb_pos=get(mbh,'Position');
                set(mbh,'Position',[mb_pos(1:2) mb_pos(3)*1.5 mb_pos(4)])
                uiwait(mbh);
            end
            ret=1;
            return
        end
        if sum(SID_UserData.interval_divisions)==d_para.num_trains
            SID_UserData.interval_divisions=cumsum(SID_UserData.interval_divisions);
            SID_UserData.interval_divisions=SID_UserData.interval_divisions(1:end-1);
        end
        if size(SID_UserData.interval_divisions,2)<size(SID_UserData.interval_divisions,1)
            SID_UserData.interval_divisions=SID_UserData.interval_divisions';
        end
        group_sep_str=num2str(SID_UserData.interval_divisions);
        if length(SID_UserData.interval_divisions)>1
            group_sep_str=regexprep(group_sep_str,'\s+',' ');
        end
        set(SID_edit,'String',group_sep_str)
        set(SID_UserData.fh,'Userdata',SID_UserData)
    end


    function SID_callback(varargin)
        d_para=getappdata(handles.figure1,'data_parameters');
        f_para=getappdata(handles.figure1,'figure_parameters');
        figure(f_para.num_fig);
        SID_UserData=get(gcf,'Userdata');

        if gcbo==SID_Reset_pushbutton || gcbo==SID_Cancel_pushbutton || gcbo==SG_fig
            if isfield(SID_UserData,'lh')
                for rc=1:size(SID_UserData.lh,1)
                    for lhc=1:size(SID_UserData.lh,2)
                        if ishandle(SID_UserData.lh(lhc))
                            delete(SID_UserData.lh(lhc))
                        end
                    end
                end
                SID_UserData.lh=[];
            end
            if isfield(SID_UserData,'interval_divisions')
                SID_UserData.interval_divisions=[];
            end
            d_para.interval_divisions=[];
            set(SID_edit,'string',[])
            if gcbo==SID_Reset_pushbutton
                set(SID_UserData.fh,'Userdata',SID_UserData)
                set(SID_UserData.um,'CallBack',{@SID_delete_interval_divisions,SID_UserData})
                set(SID_UserData.fh,'WindowButtonMotionFcn',{@SID_get_coordinates,SID_UserData},'KeyPressFcn',{@SID_keyboard,SID_UserData})
                set(SID_UserData.lh,'ButtonDownFcn',{@SID_start_move_interval_divisions,SID_UserData},'UIContextMenu',SID_UserData.cm)
                set(SID_UserData.ah,'ButtonDownFcn',{@SID_pick_interval_divisions,SID_UserData},'UIContextMenu',[])
                set(SID_UserData.spike_lh,'ButtonDownFcn',{@SID_pick_interval_divisions,SID_UserData},'UIContextMenu',[])
                set(SID_UserData.fh,'Userdata',SID_UserData)
            elseif gcbo==SID_Cancel_pushbutton || gcbo==SG_fig
                SID_UserData.cm=[];
                SID_UserData.um=[];
                set(SID_UserData.fh,'WindowButtonMotionFcn',[],'WindowButtonUpFcn',[],'KeyPressFcn',[])
                set(SID_UserData.lh,'ButtonDownFcn',[],'UIContextMenu',[])
                set(SID_UserData.ah,'ButtonDownFcn',[],'UIContextMenu',[])
                set(SID_UserData.spike_lh,'ButtonDownFcn',[],'UIContextMenu',[])
                set(SID_UserData.fh,'Userdata',SID_UserData)
                delete(SG_fig)
                set(handles.figure1,'Visible','on')
            end
        elseif gcbo==SID_OK_pushbutton || gcbo==SID_Apply_pushbutton
            [new_values,conv_ok]=str2num(get(SID_edit,'String'));
            if conv_ok==0 && ~isempty(new_values)
                set(0,'DefaultUIControlFontSize',16);
                mbh=msgbox(sprintf('The values entered are not numeric!'),'Warning','warn','modal');
                htxt = findobj(mbh,'Type','text');
                set(htxt,'FontSize',12,'FontWeight','bold')
                mb_pos=get(mbh,'Position');
                set(mbh,'Position',[mb_pos(1:2) mb_pos(3)*1.5 mb_pos(4)])
                uiwait(mbh);
                return
            end
            SID_UserData.interval_divisions=unique([new_values SID_UserData.interval_divisions]);
            set(SID_edit,'String',num2str(SID_UserData.interval_divisions))
            delete(SID_UserData.lh)
            SID_UserData=rmfield(SID_UserData,'lh');

            figure(f_para.num_fig);
            SID_cmenu=uicontextmenu;
            SID_UserData.lh=zeros(1,length(SID_UserData.interval_divisions));
            for lhc=1:length(SID_UserData.interval_divisions)
                SID_UserData.lh(lhc)=line([SID_UserData.tmin SID_UserData.tmax],(0.05+((SID_UserData.num_trains-SID_UserData.interval_divisions(lhc))/SID_UserData.num_trains))*ones(1,2),...
                    'Color',p_para.SID_col,'LineStyle',p_para.SID_ls,'LineWidth',p_para.SID_lw,'UIContextMenu',SID_cmenu);
            end
            set(SID_UserData.fh,'Userdata',SID_UserData)

            if gcbo==SID_Apply_pushbutton
                set(SID_UserData.fh,'Userdata',SID_UserData)
                set(SID_UserData.um,'CallBack',{@SID_delete_interval_divisions,SID_UserData})
                set(SID_UserData.fh,'WindowButtonMotionFcn',{@SID_get_coordinates,SID_UserData},'KeyPressFcn',{@SID_keyboard,SID_UserData})
                set(SID_UserData.lh,'ButtonDownFcn',{@SID_start_move_interval_divisions,SID_UserData},'UIContextMenu',SID_UserData.cm)
                set(SID_UserData.ah,'ButtonDownFcn',{@SID_pick_interval_divisions,SID_UserData},'UIContextMenu',[])
                set(SID_UserData.spike_lh,'ButtonDownFcn',{@SID_pick_interval_divisions,SID_UserData},'UIContextMenu',[])
                set(SID_UserData.fh,'Userdata',SID_UserData)
            else
                set(SID_panel,'HighlightColor','w')
                set(SID_edit,'Enable','off')
                set(SID_Cancel_pushbutton,'Enable','off')
                set(SID_Reset_pushbutton,'Enable','off')
                set(SID_Apply_pushbutton,'Enable','off')
                set(SID_Load_pushbutton,'Enable','off')
                set(SID_OK_pushbutton,'Enable','off')
                set(SID_UserData.tx,'str','');
                SID_UserData.marked_indy=[];
                SID_UserData.flag=0;
                SID_UserData.cm=[];
                SID_UserData.um=[];
                set(SID_UserData.fh,'WindowButtonMotionFcn',[],'WindowButtonUpFcn',[],'KeyPressFcn',[])
                set(SID_UserData.lh,'ButtonDownFcn',[],'UIContextMenu',[])
                set(SID_UserData.ah,'ButtonDownFcn',[],'UIContextMenu',[])
                set(SID_UserData.spike_lh,'ButtonDownFcn',[],'UIContextMenu',[])

                % ########################################################

                SIN_UserData=SID_UserData;

                set(SIN_panel,'Visible','on','HighlightColor','k')
                set(SIN_edit,'Visible','on')
                set(SIN_Down_pushbutton,'Visible','on')
                set(SIN_Up_pushbutton,'Visible','on')
                set(SIN_Replace_pushbutton,'Visible','on')
                set(SIN_groups_listbox,'Visible','on')
                set(SIN_names_listbox,'Visible','on')
                set(SIN_Load_pushbutton,'Visible','on')
                set(SIN_Cancel_pushbutton,'Visible','on')
                set(SIN_Reset_pushbutton,'Visible','on')
                set(SIN_OK_pushbutton,'Visible','on','FontWeight','bold')
                uicontrol(SIN_OK_pushbutton)

                set(SID_edit,'String',num2str(SIN_UserData.interval_divisions))
                SIN_UserData.interval_divisions=unique(SIN_UserData.interval_divisions);
                sep_vect=[0 SIN_UserData.interval_divisions d_para.num_trains];
                SIN_UserData.all_train_group_sizes=diff(sep_vect);
                
                SIN_UserData.num_all_train_groups=length(SIN_UserData.all_train_group_sizes);
                SIN_UserData.group_first=sep_vect(1:SIN_UserData.num_all_train_groups)+1;
                SIN_UserData.group_last=sep_vect(2:SIN_UserData.num_all_train_groups+1);
                SIN_UserData.groups_str=cell(1,SIN_UserData.num_all_train_groups);
                for stgc=1:SIN_UserData.num_all_train_groups
                    if SIN_UserData.group_first(stgc)==SIN_UserData.group_last(stgc)
                        SIN_UserData.groups_str{stgc}=num2str(SIN_UserData.group_first(stgc));
                    else
                        SIN_UserData.groups_str{stgc}=[num2str(SIN_UserData.group_first(stgc)),' - ',num2str(SIN_UserData.group_last(stgc))];
                    end
                end
                set(SIN_groups_listbox,'string',SIN_UserData.groups_str)
                if length(d_para.all_train_group_names)==SIN_UserData.num_all_train_groups
                    set(SIN_names_listbox,'string',d_para.all_train_group_names)
                    set(SIN_edit,'String',d_para.all_train_group_names{1});
                end
                set(SIN_names_listbox,'value',1)
                set(SIN_UserData.fh,'Userdata',SIN_UserData)
            end
        end
    end


    function SIN_ListBox_callback(varargin)
        figure(f_para.num_fig);
        SIN_UserData=get(gcf,'Userdata');
        SIN_UserData=get(SIN_UserData.fh,'UserData');

        lb_pos=get(SIN_names_listbox,'Value');
        lb_num_strings=length(get(SIN_groups_listbox,'String'));
        lb_strings=get(SIN_names_listbox,'String');
        if isempty(lb_strings)
            lb_strings=cell(1,lb_num_strings);
        end
        if gcbo==SIN_names_listbox
            set(SIN_edit,'String',lb_strings{lb_pos});            
            set(SIN_groups_listbox,'Value',lb_pos);            
        elseif gcbo==SIN_Down_pushbutton
            if lb_pos<lb_num_strings
                dummy=lb_strings{lb_pos+1};
                lb_strings{lb_pos+1}=lb_strings{lb_pos};
                lb_strings{lb_pos}=dummy;
                set(SIN_groups_listbox,'Value',lb_pos+1)
                set(SIN_names_listbox,'Value',lb_pos+1)
            end
        elseif gcbo==SIN_Up_pushbutton
            if lb_pos>1 && ~isempty(lb_strings{lb_pos})
                dummy=lb_strings{lb_pos-1};
                lb_strings{lb_pos-1}=lb_strings{lb_pos};
                lb_strings{lb_pos}=dummy;
                set(SIN_groups_listbox,'Value',lb_pos-1)
                set(SIN_names_listbox,'Value',lb_pos-1)
            end
        elseif gcbo==SIN_Replace_pushbutton
            SIN_str=get(SIN_edit,'String');
            if ~isempty(SIN_str)
                lb_strings{lb_pos}=SIN_str;
            end
        end
        set(SIN_names_listbox,'String',lb_strings)
        lb_pos=get(SIN_names_listbox,'Value');
        if lb_pos>0
            set(SIN_edit,'String',lb_strings{lb_pos})
        else
            lb_pos=1;
            set(SIN_groups_listbox,'Value',lb_pos);
            set(SIN_names_listbox,'Value',lb_pos);
            set(SIN_edit,'String','')
        end
        set(SIN_UserData.fh,'Userdata',SIN_UserData)
    end


    function SIN_Load_group_names(varargin)
        SIN_UserData=get(f_para.num_fig,'UserData');
        if ~isfield(d_para,'matfile')
            d_para.matfile=uigetfile('*.mat','Pick a .mat-file');
        end
        data=load(d_para.matfile);
        if isstruct(data)
            data.default_variable='';
            data.content='group names';
            SPIKY('SPIKY_select_variable',gcbo,data,handles)
            variable=getappdata(handles.figure1,'variable');
            if ~isempty(variable) && ~strcmp(variable,'A9ZB1YC8X')
                SIN_UserData.group_names=data.(variable);
            end
        else
            SIN_UserData.group_names=data;
        end
        if ~iscell(SIN_UserData.group_names) || length(SIN_UserData.group_names)~=SIN_UserData.num_all_train_groups
            if ~strcmp(variable,'A9ZB1YC8X')
                set(0,'DefaultUIControlFontSize',16);
                if ~iscell(SIN_UserData.group_names)
                    mbh=msgbox(sprintf('No string array has been chosen. Please try again!'),'Warning','warn','modal');
                else
                    mbh=msgbox(sprintf('The selected string array does not have the correct size. Please try again!'),'Warning','warn','modal');
                end
                htxt = findobj(mbh,'Type','text');
                set(htxt,'FontSize',12,'FontWeight','bold')
                mb_pos=get(mbh,'Position');
                set(mbh,'Position',[mb_pos(1:2) mb_pos(3)*1.5 mb_pos(4)])
                uiwait(mbh);
            end
            ret=1;
            return
        end
        set(SIN_names_listbox,'String',SIN_UserData.group_names,'Value',1)
        set(SIN_edit,'String',SIN_UserData.group_names{1})
        set(SIN_UserData.fh,'Userdata',SIN_UserData)
    end


    function SIN_callback(varargin)
        d_para=getappdata(handles.figure1,'data_parameters');
        f_para=getappdata(handles.figure1,'figure_parameters');
        figure(f_para.num_fig);
        SIN_UserData=get(gcf,'Userdata');
        if gcbo==SIN_Reset_pushbutton || gcbo==SIN_Cancel_pushbutton || gcbo==SG_fig
            set(SIN_edit,'string',[])
            set(SIN_names_listbox,'String',[])
            if gcbo==SIN_Reset_pushbutton
                set(SIN_UserData.fh,'Userdata',SIN_UserData)
            elseif gcbo==SIN_Cancel_pushbutton || gcbo==SG_fig
                set(SIN_UserData.fh,'Userdata',SIN_UserData)
                delete(SG_fig)
                set(handles.figure1,'Visible','on')
            end
        elseif gcbo==SIN_OK_pushbutton
            figure(f_para.num_fig);
            set(SIN_UserData.fh,'Userdata',SIN_UserData)
            
            d_para.interval_divisions=SIN_UserData.interval_divisions;
            d_para.all_train_group_sizes=diff([0 d_para.interval_divisions d_para.num_trains]);
            setappdata(handles.figure1,'data_parameters',d_para)
            set(handles.dpara_group_sizes_edit,'String',regexprep(num2str(d_para.all_train_group_sizes),'\s+',' '))
            d_para_all_train_group_names=get(SIN_names_listbox,'String');
            d_para_all_train_group_names_str='';
            for strc=1:numel(d_para.all_train_group_sizes)
                d_para_all_train_group_names_str=[d_para_all_train_group_names_str,char(d_para_all_train_group_names{strc}),'; '];
            end
            set(handles.dpara_group_names_edit,'String',d_para_all_train_group_names_str)

            set(SIN_OK_pushbutton,'UserData',1)
            set(handles.figure1,'Visible','on')
            delete(SG_fig)
        end
    end


    function SG_Close_callback(varargin)
        figure(f_para.num_fig);
        xlim([s_para.itmin-0.02*s_para.itrange s_para.itmax+0.02*s_para.itrange])
        SG_UserData=get(gcf,'Userdata');
        if isfield(SG_UserData,'lh')
            for cc=1:length(SG_UserData.lh)
                if ishandle(SG_UserData.lh(cc))
                    delete(SG_UserData.lh(cc))
                end
            end
            SG_UserData.lh=[];
        end
        if isfield(SG_UserData,'lh2')
            for rc=1:2
                for cc=1:size(SG_UserData.lh2,2)
                    if ishandle(SG_UserData.lh2(rc,cc))
                        delete(SG_UserData.lh2(rc,cc))
                    end
                end
            end
            SG_UserData.lh2=[];
        end
        SG_UserData.cm=[];
        SG_UserData.um=[];
        set(SG_UserData.fh,'WindowButtonMotionFcn',[],'WindowButtonUpFcn',[],'KeyPressFcn',[])
        set(SG_UserData.ah,'ButtonDownFcn',[],'UIContextMenu',[])
        set(SG_UserData.spike_lh,'ButtonDownFcn',[],'UIContextMenu',[])
        set(SG_UserData.tx,'str','');
        set(SG_UserData.fh,'Userdata',[])
        delete(SG_fig)
        set(handles.figure1,'Visible','on')
    end
end
