function [fileName,path, Find] = Datafilecreation(fileName,path, windowF, data)
    if ~exist('data', 'var')
        for i=1:3
            data.X.Sys1.g(i) = 2.0023;
            data.X.Sys2.g(i) = 2.0023;
            data.X.Sys1.gStrain(i) = 0.0023;
            data.X.Sys2.gStrain(i) = 0.0023;
        end
        data.X.Sys1.S = 0.5;
        data.X.Sys2.S = 0.5;
        data.X.Sys1.weight = 1;
        data.X.Sys2.weight = 1;
        data.X.Sys1.lw = 0.1;
        data.X.Sys2.lw = 0.1;
        data.Q.Sys1.lw = 0.1;
        data.Q.Sys2.lw = 0.1;
        data.comment = '';
        Find=0;
        fileName='0';
        path='0';
    end
    if ~exist('windowF', 'var')
        windowF = figure('Color', [0.9255 0.9137 0.8471], 'Name', 'Parameters finder',...
                    'DockControl', 'off', 'Units', 'Pixels', 'Position', [100 50 800 600]);
    end
    Sys1.S.edit = uicontrol('Parent', windowF,'Style', 'edit','FontSize', 15,...
                    'Units', 'normalized','String',data.X.Sys1.S,...
                    'Position', [0.12 0.80 0.1 0.05], 'Callback', @updateSys);
    Sys2.S.edit = uicontrol('Parent', windowF,'Style', 'edit','FontSize', 15,...
                    'Units', 'normalized','String',data.X.Sys2.S,...
                    'Position', [0.62 0.80 0.1 0.05], 'Callback', @updateSys);
    Sys1.weight.edit = uicontrol('Parent', windowF,'Style', 'edit','FontSize', 15,...
                    'Units', 'normalized','String',data.X.Sys1.weight,...
                    'Position', [0.38 0.80 0.1 0.05], 'Callback', @updateSys);
    Sys2.weight.edit = uicontrol('Parent', windowF,'Style', 'edit','FontSize', 15,...
                    'Units', 'normalized','String',data.X.Sys2.weight,...
                    'Position', [0.88 0.80 0.1 0.05], 'Callback', @updateSys);
                
                
    for i=1:3
        Sys1.g(i).edit = uicontrol('Parent', windowF,'Style', 'edit','FontSize', 15,...
                    'Units', 'normalized','String',data.X.Sys1.g(i),...
                    'Position', [i*0.115+0.02 0.65 0.11 0.05], 'Callback', @updateSys);
        Sys2.g(i).edit = uicontrol('Parent', windowF,'Style', 'edit','FontSize', 15,...
                    'Units', 'normalized','String',data.X.Sys2.g(i),...
                    'Position', [0.53+i*0.115 0.65 0.11 0.05], 'Callback', @updateSys);
        Sys1.gStrain(i).edit = uicontrol('Parent', windowF,'Style', 'edit','FontSize', 15,...
                    'Units', 'normalized','String',data.X.Sys1.gStrain(i),...
                    'Position', [i*0.115+0.02 0.60 0.11 0.05], 'Callback', @updateSys);
        Sys2.gStrain(i).edit = uicontrol('Parent', windowF,'Style', 'edit','FontSize', 15,...
                    'Units', 'normalized','String',data.X.Sys2.gStrain(i),...
                    'Position', [0.53+i*0.115 0.60 0.11 0.05], 'Callback', @updateSys);
    end
                
    Sys1.lwX.edit = uicontrol('Parent', windowF,'Style', 'edit','FontSize', 15,...
                    'Units', 'normalized','String',data.X.Sys1.lw,...
                    'Position', [0.3 0.50 0.1 0.05], 'Callback', @updateSys);
    Sys1.lwQ.edit = uicontrol('Parent', windowF,'Style', 'edit','FontSize', 15,...
                    'Units', 'normalized','String',data.Q.Sys1.lw,...
                    'Position', [0.3 0.45 0.1 0.05], 'Callback', @updateSys);
    Sys2.lwX.edit = uicontrol('Parent', windowF,'Style', 'edit','FontSize', 15,...
                    'Units', 'normalized','String',data.X.Sys2.lw,...
                    'Position', [0.8 0.50 0.1 0.05], 'Callback', @updateSys);
    Sys2.lwQ.edit = uicontrol('Parent', windowF,'Style', 'edit','FontSize', 15,...
                    'Units', 'normalized','String',data.Q.Sys2.lw,...
                    'Position', [0.8 0.45 0.1 0.05], 'Callback', @updateSys);                         
                
                
    Sys1.S.label = uicontrol('Parent', windowF,'Style', 'text','FontSize', 15,...
                    'Units', 'normalized','String','Spin',...
                    'Position', [0.02 0.80 0.1 0.05]);
    Sys2.S.label = uicontrol('Parent', windowF,'Style', 'text','FontSize', 15,...
                    'Units', 'normalized','String','Spin',...
                    'Position', [0.52 0.80 0.1 0.05]);
    Sys1.weight.label = uicontrol('Parent', windowF,'Style', 'text','FontSize', 15,...
                    'Units', 'normalized','String','Weight',...
                    'Position', [0.28 0.80 0.1 0.05]);
    Sys2.weight.label = uicontrol('Parent', windowF,'Style', 'text','FontSize', 15,...
                    'Units', 'normalized','String','Weight',...
                    'Position', [0.78 0.80 0.1 0.05]);
                          
       Sys1.glabel = uicontrol('Parent', windowF,'Style', 'text','FontSize', 15,...
                    'Units', 'normalized','String','g-factor',...
                    'Position', [0.02 0.65 0.11 0.05]);
       str= {'x','y','z'};
        for i=1:3
        Sys1.glabel(i) = uicontrol('Parent', windowF,'Style', 'text','FontSize', 15,...
                    'Units', 'normalized','String',str(i),...
                    'Position', [i*0.115+0.02 0.7 0.11 0.05]);
        Sys2.glabel(i) = uicontrol('Parent', windowF,'Style', 'text','FontSize', 15,...
                    'Units', 'normalized','String',str(i),...
                    'Position', [0.53+i*0.115 0.7 0.11 0.05]);
    end
             
       Sys2.glabel = uicontrol('Parent', windowF,'Style', 'text','FontSize', 15,...
                    'Units', 'normalized','String','g-factor',...
                    'Position', [0.53 0.65 0.11 0.05]);
        Sys1.gStrainlabel = uicontrol('Parent', windowF,'Style', 'text','FontSize', 15,...
                    'Units', 'normalized','String','gStrain',...
                    'Position', [0.02 0.60 0.11 0.05]);
        Sys2.gStrainlabel = uicontrol('Parent', windowF,'Style', 'text','FontSize', 15,...
                    'Units', 'normalized','String','gStrain',...
                    'Position', [0.53 0.60 0.11 0.05]);
                
    Sys1.lwX.label = uicontrol('Parent', windowF,'Style', 'text','FontSize', 15,...
                    'Units', 'normalized','String','lw Xband(mT)',...
                    'Position', [0.02 0.50 0.28 0.05]);
    Sys1.lwQ.label = uicontrol('Parent', windowF,'Style', 'text','FontSize', 15,...
                    'Units', 'normalized','String','lw Qband(mT)',...
                    'Position', [0.02 0.45 0.28 0.05]);
    Sys2.lwX.label = uicontrol('Parent', windowF,'Style', 'text','FontSize', 15,...
                    'Units', 'normalized','String','lw Xband(mT)',...
                    'Position', [0.52 0.50 0.28 0.05]);
    Sys2.lwQ.label = uicontrol('Parent', windowF,'Style', 'text','FontSize', 15,...
                    'Units', 'normalized','String','lw Qband(mT)',...
                    'Position', [0.52 0.45 0.28 0.05]);                    
                                   
                
    Sys1.label = uicontrol('Parent', windowF,'Style', 'text','FontSize', 20,...
                    'Units', 'normalized','String','System1',...
                    'Position', [0.12 0.90 0.3 0.05]);
    Sys2.label = uicontrol('Parent', windowF,'Style', 'text','FontSize', 20,...
                    'Units', 'normalized','String','System2',...
                    'Position', [0.63 0.90 0.3 0.05]);
    
    comment.label = uicontrol('Parent', windowF,'Style', 'text','FontSize', 15,...
                    'Units', 'normalized','String','Comment',...
                    'Position', [0.3 0.35 0.4 0.05]);
    comment.edit = uicontrol('Parent', windowF,'Style', 'edit','FontSize', 10,...
                    'Units', 'normalized',...
                    'Position', [0.1 0.25 0.8 0.1], 'Callback', @updateSys);
    
    if ~strcmp(fileName,'0') && ~strcmp(path,'0')
        savebutton = uicontrol('Parent', windowF, 'Style', 'pushbutton','Units', 'normalized', ...
            'Position', [0.66 0.1 0.3 0.1], 'String', 'Save', 'Callback', @saveSys); 
        openbutton = uicontrol('Parent', windowF, 'Style', 'pushbutton','Units', 'normalized', ...
            'Position', [0.04 0.1 0.3 0.1], 'String', 'Open', 'Callback', {@openSys, 0,0,0});  
        openSys(0,0,0);
    else
        savebutton = uicontrol('Parent', windowF, 'Style', 'pushbutton','Units', 'normalized', ...
            'Position', [0.55 0.1 0.35 0.1], 'String', 'Save', 'Callback', @saveSys); 
        openbutton = uicontrol('Parent', windowF, 'Style', 'pushbutton','Units', 'normalized', ...
            'Position', [0.1 0.1 0.35 0.1], 'String', 'Open', 'Callback', {@openSys, 1});  
        waitfor(windowF);
    end
    
    function updateSys(hObj, event)
        for i=1:3
                data.X.Sys1.g(i) = str2double(Sys1.g(i).edit.String);
                data.X.Sys2.g(i) = str2double(Sys2.g(i).edit.String);
                data.X.Sys1.gStrain(i) = str2double(Sys1.gStrain(i).edit.String);
                data.X.Sys2.gStrain(i) = str2double(Sys2.gStrain(i).edit.String);
        end
        data.X.Sys1.S = str2double(Sys1.S.edit.String);
        data.X.Sys2.S = str2double(Sys2.S.edit.String);
        data.X.Sys1.weight = str2double(Sys1.weight.edit.String);
        data.X.Sys2.weight = str2double(Sys2.weight.edit.String);
        data.Q = data.X;
        data.X.Sys1.lw = str2double(Sys1.lwX.edit.String);
        data.X.Sys2.lw = str2double(Sys2.lwX.edit.String);
        data.Q.Sys1.lw = str2double(Sys1.lwQ.edit.String);
        data.Q.Sys2.lw = str2double(Sys2.lwQ.edit.String);
        data.comment = comment.edit.String;
        assignin('base', 'data', data);
    end

    function saveSys(hObj, event)
        choice = questdlg('Choose', ' ', 'Create new file','Add to existing file', 'Create new file');
            switch choice
            case 'Create new file'
                 [fileName,path] = uiputfile('*.mat','Name of the file, where data will be stored'); 
            case 'Add to existing file'
                if exist('fileName', 'var') && exist('path', 'var') 
                    [fileName,path] = uigetfile('*.mat','Name of the file, where data will be stored', strcat(path, fileName)); 
                else
                    [fileName,path] = uigetfile('*.mat','Name of the file, where data will be stored'); 
                end
            end
 
        updateSys();

            switch choice
            case 'Create new file'
               save(strcat(path,fileName),'data');
               msgbox('Done');
               Find=1;
            case 'Add to existing file'
               save(strcat(path,fileName),'data', '-v6', '-append');
               msgbox('Done'); 
               Find=1;
            end    
    end

    function openSys(~,~, fl)
         if fl == 1
            [fileName,path, fileind] = uigetfile('*.mat','Name of the file, where data will be stored'); 
            if fileind
                Find=1;
            end
         end
         D = load([path,fileName],'-mat');
         data = D.data;
         Sys1.S.edit.String = data.Q.Sys1.S;
         Sys2.S.edit.String = data.Q.Sys2.S;
         Sys1.weight.edit.String = data.Q.Sys1.weight;
         Sys2.weight.edit.String = data.Q.Sys2.weight;
         for i=1:3
            Sys1.g(i).edit.String = data.Q.Sys1.g(i);
            Sys2.g(i).edit.String = data.Q.Sys2.g(i);
            Sys1.gStrain(i).edit.String = data.Q.Sys1.gStrain(i);
            Sys2.gStrain(i).edit.String = data.Q.Sys2.gStrain(i);
        end

        Sys1.lwX.edit.String = data.X.Sys1.lw;
        Sys1.lwQ.edit.String = data.Q.Sys1.lw;
        Sys2.lwX.edit.String = data.X.Sys2.lw;
        Sys2.lwQ.edit.String = data.Q.Sys2.lw;
        
        comment.edit.String = data.comment; 
    end
end