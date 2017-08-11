function varargout = am_PDFcard(varargin)
    gui_Singleton = 1;
    gui_State = struct('gui_Name',          mfilename, ...
                       'gui_Singleton',     gui_Singleton, ...
                       'gui_OpeningFcn',    @am_PDFcard_OpeningFcn, ...
                       'gui_OutputFcn',     @am_PDFcard_OutputFcn, ...
                       'gui_LayoutFcn',     [], ...
                       'gui_Callback',      []);
    if nargin && ischar(varargin{1})
       gui_State.gui_Callback = str2func(varargin{1});
    end
    
    if nargout
        [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
    else
        gui_mainfcn(gui_State, varargin{:});
    end
end

% --- Executes just before am_PDFcard is made visible.
function am_PDFcard_OpeningFcn(hObject, eventdata, handles, varargin)
    handles.output = hObject;
    guidata(hObject,handles);
    if nargin == 3
        initial_dir = pwd;
    elseif nargin > 4
        if strcmpi(varargin{1},'dir')
            if exist(varargin{2},'dir')
                initial_dir = varargin{2};
            else
                errordlg('Input argument must be a valid directory','Input Argument Error!')
                return
            end
        else
            errordlg('Unrecognized input argument','Input Argument Error!');
            return;
        end
    end
    load_listbox(initial_dir,handles);
end

% --- Outputs from this function are returned to the command line.
function varargout = am_PDFcard_OutputFcn(hObject, eventdata, handles)
    varargout{1} = handles.output;
end

function listbox1_CreateFcn(hObject, eventdata, handles) 
    set(hObject,'BackgroundColor','white');
end

function figure1_CreateFcn(hObject, eventdata, handles)
    set(hObject,'color','w'); setappdata(hObject, 'StartPath', pwd); addpath(pwd);
end

function figure1_DeleteFcn(hObject, eventdata, handles)
    if isappdata(hObject, 'StartPath'); rmpath(getappdata(hObject, 'StartPath')); end
end

function axes1_CreateFcn(hObject, eventdata, handles)
    set(hObject,'color','w');
    handles.axes1.YScale='log';
    handles.axes1.YLabel.String='Intensity';
    handles.axes1.XLabel.String='2\theta [deg]';
end

function edit1_CreateFcn(hObject, eventdata, handles)
    set(hObject,'BackgroundColor','white'); hObject.String={''};
end

function listbox1_Callback(hObject, eventdata, handles)
    import am_dft.*
    index_selected = get(handles.listbox1,'Value');
    file_list = get(handles.listbox1,'String');
    filename = file_list{index_selected};
    if     strcmp(get(handles.figure1,'SelectionType'),'normal')
        if contains(filename,'poscar')
            uc = load_poscar(filename);
            plot_millers(handles.axes1,uc); 
    %         handles.axes1.YScale='log';
    %         handles.axes1.YLabel.String='Intensity';
        end
    elseif strcmp(get(handles.figure1,'SelectionType'),'open')
        if  handles.is_dir(index_selected)
            cd(filename);
            load_listbox(pwd,handles);
        elseif contains(filename,'poscar')
            fprintf('%s\n',repmat('-',100));
            fprintf('%s\n',filename);
            type(filename)
        end
    end
end

function [h]     = plot_millers(h,uc)
    %
    import am_lib.* am_mbe.*
    
    % generate hkl list
    hv = get_atomic_emission_line_energy(get_atomic_number('Cu'),'kalpha1'); lambda = get_photon_energy(hv);
    % get miller indicies [hkl]
    N=6; hkl=permn_([0:N],3).'; k=inv(uc.bas*0.1).'*hkl; 
    % % exclude values which cannot be reached by diffractometer
    th2_max=110; ex_=lambda*normc_(k)/2<sind(th2_max/2); hkl=hkl(:,ex_); k=k(:,ex_);
    % % find unique hkls (proper way to do it would be to use symmetry relations, but whatever)
    [~,i] = unique(rnd_(normc_(k))); k=k(:,i(2:end)); hkl=hkl(:,i(2:end));
    % compute angles
    th2 = 2*asind( lambda * normc_(k) / 2 );

    % get structure factor
    [Z,~,i] = unique(get_atomic_number({uc.symb{uc.species}}));
    f0 = permute(get_atomic_xray_form_factor(Z,hv,th2),[1,3,2]);
    Fhkl = sum(f0(i,:).*exp(2i*pi*uc.tau.'*hkl),1); Fhkl2= abs(Fhkl).^2; 
    Fhkl2=Fhkl2./max(Fhkl2);

    % plot Bragg peaks
    ex_=Fhkl2>0.01; label_threshold = 0; plot(h,0,0);
    for i = 1:numel(th2); if ex_(i)
        line(h,[th2(i),th2(i)],[0,Fhkl2(i)]);
        if Fhkl2(i) > label_threshold
            text(h,th2(i),Fhkl2(i),sprintf('  [%i%i%i] %.3f^\\circ',hkl(:,i),th2(i)),'rotation',90);
        end
    end; end
    set(h,'XLim',[0 115],'YLim',[0 2],'XTick',[0:10:130],'XMinorTick','on','YTick',[]);
    h.XLabel.String='2\theta [deg]'; h.YLabel.String='Intensity [a.u.]'; 
%     box on; set(gca,'YTick',[]); grid on; set(gca,'XMinorGrid','on');

end

% ------------------------------------------------------------
% Read the current directory and sort the names
% ------------------------------------------------------------
function load_listbox(dir_path,handles)
    cd(dir_path);
    dir_struct = dir(dir_path); file_names = {dir_struct.name}; isdir = [dir_struct.isdir];
    element_txt = strtrim(handles.edit1.String);
    if ~strcmp(element_txt,'')
        % if not empty
        element_list = strsplit(element_txt{:},' ');
        % get files names
        if handles.togglebutton2.Value==1
            % any
            ex_ = contains(file_names,element_list,'IgnoreCase',true);
        else
            % all
            ex_ = true(1,numel(file_names));
            for i = 1:numel(element_list)
                ex_ = and(ex_,contains(file_names,element_list{i},'IgnoreCase',true));
            end
        end
    else
        % no filter
        ex_ = true(1,numel(file_names));
    end
    % append directories
    ex_ = or( ex_, isdir );
    % remove hidden files 
    ex_ = and( ex_, ~startsWith(file_names,'.') );
    % append . and ..
    ex_ = or( ex_, strcmp(file_names,'.') );
    ex_ = or( ex_, strcmp(file_names,'..'));
    % sort by name and put directories first
    [~,handles.sorted_index] = sortrows( strcat(num2cell((~isdir)*1+48).',{dir_struct.name}.') );
    % filter based on ex_
    ex_ = ex_(handles.sorted_index);
    handles.file_names = file_names(handles.sorted_index); 
    handles.file_names = handles.file_names(ex_);
    handles.is_dir     = [dir_struct.isdir];
    handles.is_dir     = handles.is_dir(handles.sorted_index);
    handles.is_dir     = handles.is_dir(ex_);
    guidata(handles.figure1,handles)
    set(handles.listbox1,'String',handles.file_names,'Value',1)
end

function edit1_Callback(hObject, eventdata, handles)
    load_listbox(pwd,handles)
end

function togglebutton2_Callback(hObject, eventdata, handles)
    if hObject.Value == 1
        hObject.String={'Any'};
    else
        hObject.String={'All'};
    end
    load_listbox(pwd,handles)
end
