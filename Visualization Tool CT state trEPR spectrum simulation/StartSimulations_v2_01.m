% The visualization toolbox makes trEPR spectra simulations of CT states 
% easier. The tool can be used for powders - all orientations in the 
% magnetic field are possible and equally probable. Two models can be used  
% for relative orientations between radicals - fixed and random relative 
% orientation.


function StartSimulations_v2(data)
    clear all;
    clc;close all;
    global Data Dat Xband Qband Q X DS...
           xdataX  ydataX xdataQ ydataQ ydataXtot ydataQtot...
           xdataXfp xdataQfp ydataXfp ydataQfp ydataXtotfp ydataQtotfp ...
           ydataXtotnorm ydataQtotnorm ydataXtotfpnorm ydataQtotfpnorm;
    % DS. - save into data file
    % Q., X. - downloaded from Xband and Qband file 
    % Data - download from file with system data
    % Dat - parameters which need to be global and related to Data
    % Xband, Qband - parameters related to Xband and Qband
   addpath('private');
    Flagfinish = 0;
    maxJ=60;
    minJ=-60;                       
    maxD=60;
    minD=-60;
    maxJFP=60;
    minJFP=-60;                       
    maxDFP=60;
    minDFP=-60;
    maxa=pi;
    mina=0;
    DS.addbrX = 0.0;
    DS.addbrQ = 0.0;
    DS.JFP = 2;
    DS.DFP = 3;
    DS.J = 2;
    DS.D = 3;
    DS.E = 0;
    DS.normXFP = 1;
    DS.normQFP = 1;
    DS.normX = 1;
    DS.normQ = 1;
    DS.JFPn = 1;
    DS.DFPn = 1;
    DS.JFPw = 0;
    DS.DFPw = 0;
    DS.npFP = 8;
    DS.np = 20;
    DS.smoothX = 5;
    DS.smoothQ = 5;
    DS.fitlevelX = 1;
    DS.fitlevelQ = 1;
    Dat.regime = 'bg';
    h=NaN;
    %% windows creation, labels
    
    window2 = figure('Color', [0.92 0.92 0.84], 'Name', '',...
                    'Units', 'Pixels', 'Position', [100 50 800 600]);
    Xband2 = axes('Parent', window2,'Units', 'normalized', 'Position', [0.07 0.55 0.40 0.40]);
    Qband2 = axes('Parent', window2,'Units', 'normalized','Position', [0.55 0.55 0.40 0.40]);
    Xband2FP = axes('Parent', window2,'Units', 'normalized', 'Position', [0.07 0.05 0.40 0.40]);
    Qband2FP = axes('Parent', window2,'Units', 'normalized','Position', [0.55 0.05 0.40 0.40]);
    
    window = figure('Color', [0.92 0.92 0.84], 'Name', 'Parameters finder',...
                    'Units', 'Pixels', 'Position', [100 50 800 600]);
    Xband.plot = axes('Parent', window,'Units', 'normalized', 'Position', [0.03 0.39 0.415 0.58]);
    Qband.plot = axes('Parent', window,'Units', 'normalized','Position', [0.53 0.39 0.415 0.58]);
    fullpowderrect = annotation('rectangle','Units', 'normalized', 'Position', [0.00 0.178 0.5 0.085], 'FaceColor', [0.8, 0.9, 1]);
    fixdorientrect = annotation('rectangle','Units', 'normalized', 'Position', [0.00 0.095 0.5 0.085], 'FaceColor', [1, 0.9, 0.9]);
    systemrect = annotation('rectangle','Units', 'normalized', 'Position', [0.00 0.00 0.5 0.095], 'FaceColor', [0.9 0.9 0.9]);
    fullpowdertext = uicontrol('Parent', window,'Style', 'text','String', 'FULL POWDER (FP)','FontSize', 9,'FontWeight', 'bold',...
                    'Units', 'normalized','Position', [0.005 0.238 0.07 0.021],'Backgroundcolor',[0.8, 0.9, 1]);
    fixedorienttext = uicontrol('Parent', window,'Style', 'text','String', 'FRO IN POWER','FontSize', 9,'FontWeight', 'bold',...
                    'Units', 'normalized','Position', [0.005 0.153 0.09 0.021],'Backgroundcolor',[1, 0.9, 0.9]);
    systemtext = uicontrol('Parent', window,'Style', 'text','String', 'SYSTEM PARAMETERS','FontSize', 9,'FontWeight', 'bold',...
                    'Units', 'normalized','Position', [0.005 0.068 0.08 0.021],'Backgroundcolor',[0.9 0.9 0.9]);
   
            Xband.B1.edit = uicontrol('Parent', window, 'Style', 'edit', 'FontSize', 9, 'Units', 'normalized',...
                'Position', [0.05 0.27 0.06 0.025], 'Callback', @plotwindow);
            Xband.B2.edit = uicontrol('Parent', window, 'Style', 'edit', 'FontSize', 9, 'Units', 'normalized',...
                'Position', [0.15 0.27 0.06 0.025], 'Callback', @plotwindow);
            Xband.time.edit = uicontrol('Parent', window, 'Style', 'edit', 'FontSize', 9, 'Units', 'normalized',...
                'Position', [0.45 0.30 0.05 0.025], 'Callback', @updatetime2X);
            Xband.filebutton = uicontrol('Parent', window, 'Style', 'pushbutton','Units', 'normalized', ...
                'Position', [0.51 0.001 0.15 0.025],'FontSize', 8, 'Callback', @openXband);
            Xband.time.slider = uicontrol('Parent', window, 'Style', 'slider', 'Units', 'normalized',...
                'Position', [0.05 0.30 0.4 0.03], 'Callback', @updatetimeX, 'SliderStep', [0.001 0.001]);
            Qband.time.slider = uicontrol('Parent', window, 'Style', 'slider', 'Units', 'normalized',...
                 'Position', [0.55 0.30 0.4 0.03],'Callback', @updatetimeQ, 'SliderStep', [0.001 0.001]);
            Qband.B1.edit = uicontrol('Parent', window, 'Style', 'edit', 'FontSize', 9, 'Units', 'normalized',...
                 'Position', [0.84 0.27 0.06 0.025], 'Callback', @plotwindow);
            Qband.B2.edit = uicontrol('Parent', window, 'Style', 'edit', 'FontSize', 9, 'Units', 'normalized',...
                 'Position', [0.94 0.27 0.06 0.025], 'Callback', @plotwindow);  
            Qband.time.edit = uicontrol('Parent', window, 'Style', 'edit', 'FontSize', 9, 'Units', 'normalized',...
                'Position', [0.95 0.30 0.05 0.025], 'Callback', @updatetime2Q);
            Qband.filebutton = uicontrol('Parent', window, 'Style', 'pushbutton','Units', 'normalized',...
                'Position', [0.84 0.001 0.15 0.025],'FontSize', 8, 'Callback', @openQband);

            JFP.slider = uicontrol('Parent', window, 'Style', 'slider','Units', 'normalized', ...
                'Position', [0.025 0.215 0.08 0.025], 'Min', minJFP,'Max', maxJFP,...
                'Callback', @updateGraphFP, 'SliderStep', [0.001 0.001],'Backgroundcolor',[0.5 0.5 0.5]);
            DFP.slider = uicontrol('Parent', window, 'Style', 'slider','Units', 'normalized',...
                'Position', [0.025 0.185 0.08 0.025],'Min', minDFP,'Max', maxDFP,...
                'Callback', @updateGraphFP, 'SliderStep', [0.001 0.001],'Backgroundcolor',[0.5 0.5 0.5]);
            
            J.slider = uicontrol('Parent', window,'Style', 'slider','Units', 'normalized',...
                'Position', [0.025 0.13 0.08 0.025],'Min', minJ, 'Max', maxJ,...
                'Callback', @updateGraph1, 'SliderStep', [0.001 0.001], 'Backgroundcolor',[1 0.8 0.8]);
            D.slider = uicontrol('Parent', window, 'Style', 'slider', 'Units', 'normalized',...
                 'Position', [0.025 0.10 0.08 0.025],'Min', minD,'Max', maxD,...
                 'Callback', @updateGraph1, 'SliderStep', [0.001 0.001], 'Backgroundcolor',[1 0.8 0.8]);

             for i=1:6
                Dat.aslider(i) = uicontrol('Parent', window,'Style', 'slider','Units', 'normalized',...
                     'Position', [0.223+mod(i-1,3)*0.094, 0.13-round(i/8)*0.03, 0.052, 0.025],'Min', mina,'Max', maxa,...
                     'Callback', @updateGraph1, 'SliderStep', [0.005 0.005], 'Backgroundcolor',[1 0.8 0.8]);
            end

            Qband.normQ.slider = uicontrol('Parent', window,'Style', 'slider','Units', 'normalized',...
                  'Position', [0.95 0.39 0.02 0.58],'Min', -100,'Max', 100, 'Value', 1,...
                  'background',[1.0 0.5 0.5],'Callback', @updateGraph1Norm, 'SliderStep', [0.001 0.001]);
            Xband.normX.slider = uicontrol('Parent', window,'Style', 'slider','Units', 'normalized',...
                  'Position', [0.45 0.39 0.02 0.58],'Min', -100,'Max', 100, 'Value', 1,...
                  'background',[1.0 0.5 0.5],'Callback', @updateGraph1Norm, 'SliderStep', [0.001 0.001]);      
            Xband.normXFP.slider = uicontrol('Parent', window,'Style', 'slider','Units', 'normalized',...
                  'Position', [0.475 0.39 0.02 0.58],'Min', -100, 'Max', 100, 'Value', 1,...
                  'background',[0.4 0.5 1],'Callback', @updateGraphFPNorm, 'SliderStep', [0.001 0.001]);
            Qband.normQFP.slider = uicontrol('Parent', window,'Style', 'slider','Units', 'normalized',...
                  'Position', [0.975 0.39 0.02 0.58],'Min', -100,'Max', 100, 'Value', 1,...
                  'background',[0.4 0.5 1],'Callback', @updateGraphFPNorm, 'SliderStep', [0.0001 0.001]);       

            JFP.edit = uicontrol('Parent', window, 'Style', 'edit', 'FontSize', 9, 'Units', 'normalized',...
                            'Position', [0.106 0.215 0.04 0.025], 'Callback', @updateGraphFP2, 'Backgroundcolor',[0.9 0.9 0.9]);
            DFP.edit = uicontrol('Parent', window, 'Style', 'edit', 'FontSize', 9, 'Units', 'normalized',...
                            'Position', [0.106 0.185 0.04 0.025], 'Callback', @updateGraphFP2, 'Backgroundcolor',[0.9 0.9 0.9]);
            EFP.edit = uicontrol('Parent', window, 'Style', 'edit', 'FontSize', 9, 'Units', 'normalized',...
                            'Position', [0.147 0.185 0.04 0.025], 'Callback', @updateGraphFP2, 'Backgroundcolor',[0.9 0.9 0.9]);

            Dat.JFPwidth.edit = uicontrol('Parent', window, 'Style', 'edit', 'FontSize', 9, 'Units', 'normalized',...
                            'Position', [0.215 0.214 0.04 0.025], 'Callback', @updateGraphFP2, 'Backgroundcolor',[0.7 0.7 0.7]);
            Dat.DFPwidth.edit = uicontrol('Parent', window, 'Style', 'edit', 'FontSize', 9, 'Units', 'normalized',...
                            'Position', [0.215 0.187 0.04 0.025], 'Callback', @updateGraphFP2, 'Backgroundcolor',[0.7 0.7 0.7]);
            Dat.JFPnumber.edit = uicontrol('Parent', window, 'Style', 'edit', 'FontSize', 9, 'Units', 'normalized',...
                            'Position', [0.395 0.214 0.04 0.025], 'Callback', @updateGraphFP2, 'Backgroundcolor',[0.7 0.7 0.7]);
            Dat.DFPnumber.edit = uicontrol('Parent', window, 'Style', 'edit', 'FontSize', 9, 'Units', 'normalized',...
                            'Position', [0.395 0.187 0.04 0.025], 'Callback', @updateGraphFP2, 'Backgroundcolor',[0.7 0.7 0.7]);
            Dat.JFPrange = uicontrol('Parent', window, 'Style', 'text', 'FontSize', 9, 'Units', 'normalized',...
                            'Position', [0.265 0.214 0.12 0.025], 'Backgroundcolor',[0.7 0.7 0.7]);
            Dat.DFPrange = uicontrol('Parent', window, 'Style', 'text', 'FontSize', 9, 'Units', 'normalized',...
                            'Position', [0.265 0.187 0.12 0.025], 'Backgroundcolor',[0.7 0.7 0.7]);

            J.edit = uicontrol('Parent', window, 'Style', 'edit', 'FontSize', 9, 'Units', 'normalized',...
                            'Position', [0.106 0.13 0.04 0.025], 'Callback', @updateGraph2, 'Backgroundcolor',[1 0.95 0.95]);
            D.edit = uicontrol('Parent', window, 'Style', 'edit', 'FontSize', 9, 'Units', 'normalized',...
                            'Position', [0.106 0.10 0.04 0.025], 'Callback', @updateGraph2, 'Backgroundcolor',[1 0.95 0.95]);
            E.edit = uicontrol('Parent', window, 'Style', 'edit', 'FontSize', 9, 'Units', 'normalized',...
                            'Position', [0.147 0.10 0.04 0.025], 'Callback', @updateGraph2, 'Backgroundcolor',[1 0.95 0.95]);
            for i=1:6
                Dat.aedit(i) = uicontrol('Parent', window, 'Style', 'edit', 'FontSize', 9, 'Units', 'normalized',...
                            'Position', [0.273+mod(i-1,3)*0.094, 0.13-round(i/8)*0.03, 0.038, 0.025], 'Callback', @updateGraph2, 'Backgroundcolor',[1 0.95 0.95], 'String', i);
            end            
            Xband.normX.edit = uicontrol('Parent', window, 'Style', 'edit', 'FontSize', 9, 'Units', 'normalized',...
                            'Position', [0.445 0.97 0.0275 0.03], 'Callback', @updateGraph2Norm); 
            Qband.normQ.edit = uicontrol('Parent', window, 'Style', 'edit', 'FontSize', 9, 'Units', 'normalized',...
                            'Position', [0.945 0.97 0.0275 0.03], 'Callback', @updateGraph2Norm);    
            Xband.normXFP.edit = uicontrol('Parent', window, 'Style', 'edit', 'FontSize', 9, 'Units', 'normalized',...
                            'Position', [0.4725 0.97 0.0275 0.03], 'Callback', @updateGraphFP2Norm);  
            Qband.normQFP.edit = uicontrol('Parent', window, 'Style', 'edit', 'FontSize', 9, 'Units', 'normalized',...
                            'Position', [0.9725 0.97 0.0275 0.03], 'Callback', @updateGraphFP2Norm);   
            Dat.addbrX.edit = uicontrol('Parent', window, 'Style', 'edit', 'FontSize', 9, 'Visible', 'off', ...
                            'Units', 'normalized','Position', [0.5425 0.138 0.06 0.025], 'Callback', @addbroadX);  
            Dat.addbrQ.edit = uicontrol('Parent', window, 'Style', 'edit', 'FontSize', 9, 'Visible', 'off', ...
                            'Units', 'normalized','Position', [0.6675 0.138 0.06 0.025], 'Callback', @addbroadQ);     
            Dat.npFP.edit = uicontrol('Parent', window, 'Style', 'edit', 'FontSize', 9,  'Visible', 'off', ...
                            'Units', 'normalized','Position', [0.67 0.088 0.06 0.025], 'Callback', {@npFP,1});   
            Dat.np.edit = uicontrol('Parent', window, 'Style', 'edit', 'FontSize', 9, 'Visible', 'off', ...
                            'Units', 'normalized','Position', [0.67 0.062 0.06 0.025], 'Callback', {@np,1});   

            Xband.fitlevel.edit = uicontrol('Parent', window, 'Style', 'edit', 'FontSize', 9, 'Visible', 'off',...
                        'Units', 'normalized','Position', [0.5425 0.136 0.06 0.025], 'Callback', @Xbgcorr2);  
            Qband.fitlevel.edit = uicontrol('Parent', window, 'Style', 'edit', 'FontSize', 9, 'Visible', 'off',...
                         'Units', 'normalized','Position', [0.6675 0.136 0.06 0.025], 'Callback', @Qbgcorr2);    
            Xband.smooth.edit = uicontrol('Parent', window, 'Style', 'edit', 'FontSize', 9, 'Visible', 'off',...
                         'String', 5,'Units', 'normalized','Position', [0.5425 0.076 0.06 0.025], 'Callback', @Xbgcorr2);  
            Qband.smooth.edit = uicontrol('Parent', window, 'Style', 'edit', 'FontSize', 9,'Visible', 'off',...
                         'Value', 5,'Units', 'normalized','Position', [0.6675 0.076 0.06 0.025], 'Callback', @Qbgcorr2);    
                     
    for i=1:6
       DS.a(i).a=0;
    end
    
    R2X_edit = uicontrol('Parent', window, 'Style', 'edit', 'FontSize', 9, 'Units', 'normalized',...
                    'Position', [0.30 0.27 0.06 0.025]); 
    R2Q_edit = uicontrol('Parent', window, 'Style', 'edit', 'FontSize', 9, 'Units', 'normalized',...
                    'Position', [0.73 0.27 0.06 0.025]);             
    R2fo_edit = uicontrol('Parent', window, 'Style', 'edit', 'FontSize', 9, 'Units', 'normalized',...
                    'Position', [0.44 0.27 0.06 0.025]); 
    R2fp_edit = uicontrol('Parent', window, 'Style', 'edit', 'FontSize', 9, 'Units', 'normalized',...
                    'Position', [0.58 0.27 0.06 0.025]); 
    widthfp.label = uicontrol('Parent', window, 'Style', 'text', 'FontSize', 9, 'Units', 'normalized',...
                    'String', 'Width of range R(MHz)','Position', [0.168 0.238 0.14 0.021], 'Backgroundcolor',[0.8, 0.9, 1]);
    valuefp.label = uicontrol('Parent', window, 'Style', 'text', 'FontSize', 9, 'Units', 'normalized',...
                    'String', 'Value(MHz)','Position', [0.089 0.238 0.089 0.021], 'Backgroundcolor',[0.8, 0.9, 1]);
    valuefix.label = uicontrol('Parent', window, 'Style', 'text', 'FontSize', 9, 'Units', 'normalized',...
                    'String', 'Value(MHz)','Position', [0.089 0.153 0.089 0.021], 'Backgroundcolor',[1, 0.9, 0.9]);
    rangefp.label = uicontrol('Parent', window, 'Style', 'text', 'FontSize', 9, 'Units', 'normalized',...
                    'String', 'Range(MHz)','Position', [0.278 0.238 0.10 0.021], 'Backgroundcolor',[0.8, 0.9, 1]);
    numberfp.label = uicontrol('Parent', window, 'Style', 'text', 'FontSize', 9, 'Units', 'normalized',...
                    'String', 'N of points','Position', [0.368 0.238 0.10 0.021], 'Backgroundcolor',[0.8, 0.9, 1]);
    JFP.label = uicontrol('Parent', window,'Style', 'text','String', 'J','FontSize', 9,'FontWeight', 'bold',...
                    'Units', 'normalized','Position', [0.003 0.215 0.015 0.021],'Backgroundcolor',[0.8, 0.9, 1]);
    DFP.label = uicontrol('Parent', window,'Style', 'text','String', 'D','FontSize', 9,'FontWeight', 'bold',...
                    'Units', 'normalized','Position', [0.003 0.185 0.015 0.021],'Backgroundcolor',[0.8, 0.9, 1]);
    J.label = uicontrol('Parent', window,'Style', 'text','String', 'J','FontSize', 9,'FontWeight', 'bold',...
                    'Units', 'normalized','Position', [0.003 0.13 0.015 0.021], 'Backgroundcolor',[1 0.9 0.9]);
    D.label = uicontrol('Parent', window,'Style', 'text','String', 'D','FontSize', 9,'FontWeight', 'bold',...
                    'Units', 'normalized','Position', [0.003 0.10 0.015 0.021],'Backgroundcolor',[1 0.9 0.9]);
    angle.label = uicontrol('Parent', window,'Style', 'text','String', 'E u l e r    a n g l e s    i n    r a d i a n s',...
                    'FontSize', 9,'Units', 'normalized','Position', [0.22 0.153 0.27 0.021], 'Backgroundcolor',[1, 0.9, 0.9]);
    anSys1.label = uicontrol('Parent', window,'Style', 'text','Units', 'normalized', 'FontSize', 9,...
                     'Position', [0.183, 0.13, 0.04, 0.025], 'Backgroundcolor',[1 0.8 0.8],'String', 'Sys1');  
    anSys2.label = uicontrol('Parent', window,'Style', 'text','Units', 'normalized', 'FontSize', 9,...
                     'Position', [0.183, 0.10, 0.04, 0.025], 'Backgroundcolor',[1 0.8 0.8],'String', 'Sys2'); 
    drawasnglesbutton = uicontrol('Parent', window, 'Style', 'pushbutton','Units', 'normalized', ...
                'Position', [0.180 0.153 0.04 0.025],'FontSize', 8, 'String', 'Draw','Callback', @draworient);         
                        
    Dat.BGrect = annotation('rectangle','Units', 'normalized', 'Position', [0.51 0.03 0.25 0.18], 'FaceColor', [0.92 0.92 0.84],'visible','off');    
    Dat.Xbandrect = annotation('rectangle','Units', 'normalized', 'Position', [0.51 0.06 0.125 0.15], 'FaceColor', [0.8 0.8 0.8],'visible','off'); %[0.95 0.95 1]
    Dat.Qbandrect = annotation('rectangle','Units', 'normalized', 'Position', [0.635 0.06 0.125 0.15], 'FaceColor', [0.7 0.7 0.7],'visible','off');%[0.88 0.88 1]
    Dat.Xbandlabel = uicontrol('Parent', window, 'Style', 'text', 'FontSize', 9, 'Units', 'normalized','FontWeight', 'bold',...
        'String','Xband','Position', [0.52 0.185 0.1 0.020], 'Backgroundcolor',[0.8 0.8 0.8],'visible','off');
    Dat.Qbandlabel = uicontrol('Parent', window, 'Style', 'text', 'FontSize', 9, 'Units', 'normalized','FontWeight', 'bold',...
        'String','Qband','Position', [0.645 0.185 0.1 0.020], 'Backgroundcolor',[0.7 0.7 0.7],'visible','off');
    Dat.addbrlabel = uicontrol('Parent', window, 'Style', 'text', 'FontSize', 9, 'Units', 'normalized',...
        'String','Additional broadening(mT)','Position', [0.52 0.166 0.23 0.019], 'Backgroundcolor',[0.95, 0.95, 0.95],'visible','off');
    Dat.numberofpointslabel = uicontrol('Parent', window, 'Style', 'text', 'FontSize', 9, 'Units', 'normalized',...
        'String','Number of points','Position', [0.56 0.115 0.15 0.019], 'Backgroundcolor',[0.95, 0.95, 0.95],'visible','off');
    Dat.npfulllabel = uicontrol('Parent', window, 'Style', 'text', 'FontSize', 9, 'Units', 'normalized',...
        'String','For full powder','Position', [0.515 0.09 0.155 0.019], 'Backgroundcolor',[0.95, 0.95, 0.95],'visible','off');
    Dat.npfixlabel = uicontrol('Parent', window, 'Style', 'text', 'FontSize', 9, 'Units', 'normalized',...
        'String','For fixed orientation','Position', [0.515 0.065 0.155 0.019], 'Backgroundcolor',[0.95, 0.95, 0.95],'visible','off');
    Dat.systembutton = uicontrol('Parent', window, 'Style', 'pushbutton','Units', 'normalized',...
        'Position', [0.09 0.066 0.40 0.027],'FontSize', 9,'String', 'Change system parameters', 'Callback', @systemchangeregime);
     
    Dat.fitlevellabel = uicontrol('Parent', window,'Style', 'text','String', 'Fitlevel(polynom power):',...
                    'FontSize', 9,'Units', 'normalized','Position', [0.52 0.166 0.23 0.019], 'Backgroundcolor',[0.95, 0.95, 0.95],'visible','off');              
    Dat.smoothlabel = uicontrol('Parent', window,'Style', 'text','String', 'Smoothing power:',...
                    'FontSize', 9,'Units', 'normalized','Position', [0.52 0.106 0.23 0.019],'visible','off');              
    Dat.bgwithoutsmoothlabel = uicontrol('Parent', window,'Style', 'text','String', 'Bg corrected without smoothing:',...
                    'FontSize', 9,'Units', 'normalized','Position', [0.51 0.035 0.23 0.025],'visible','off');              
    Dat.bgwithoutsmooth_checkbox = uicontrol('Parent', window,'Style', 'checkbox','FontSize', 30,'visible','off',...
                    'Backgroundcolor',[0.92 0.92 0.84],'Value',1,'Units', 'normalized','Position', [0.74 0.035 0.02 0.02], 'Callback', @plotwindow); 
                
    BX1_label = uicontrol('Parent', window,'Style', 'text','String', 'B1:','FontSize', 9,...
                    'Units', 'normalized','Position', [0.02 0.27 0.03 0.025]);       
    BX2_label = uicontrol('Parent', window,'Style', 'text','String', 'B2:','FontSize', 9,...
                    'Units', 'normalized','Position', [0.12 0.27 0.03 0.025]);
    BQ1_label = uicontrol('Parent', window,'Style', 'text','String', 'B1:','FontSize', 9,...
                    'Units', 'normalized','Position', [0.81 0.27 0.03 0.025]);       
    BQ2_label = uicontrol('Parent', window,'Style', 'text','String', 'B2:','FontSize', 9,...
                    'Units', 'normalized','Position', [0.91 0.27 0.03 0.025]); 
    R2X_label = uicontrol('Parent', window,'Style', 'text','String', 'R^2 X:','FontSize', 9,...
                    'Units', 'normalized','Position', [0.23 0.27 0.07 0.025]);
    R2Q_label = uicontrol('Parent', window,'Style', 'text','String', 'R^2 Q:','FontSize', 9,...
                    'Units', 'normalized','Position', [0.65 0.27 0.07 0.025]);            
    R2fo_label = uicontrol('Parent', window,'Style', 'text','String', 'R^2 FRO:','FontSize', 9,...
                    'Units', 'normalized','Position',  [0.37 0.27 0.07 0.025]);
    R2fp_label = uicontrol('Parent', window,'Style', 'text','String', 'R^2 FP:','FontSize', 9,...
                    'Units', 'normalized','Position',  [0.51 0.27 0.07 0.025]);
    
    fullpowder_label = uicontrol('Parent', window,'Style', 'text','String', 'FP:','FontSize', 9,...
                    'Units', 'normalized','Position', [0.583 0.24 0.043 0.021],'Backgroundcolor',[0.7 0.8 1]);
    fullpowder_checkbox = uicontrol('Parent', window,'Style', 'checkbox','FontSize', 30,...
                    'Backgroundcolor',[0.8 0.8 0.8],'Value',1,'Units', 'normalized','Position', [0.626 0.24 0.02 0.021], 'Callback', @plotwindow);
    fixorient_label = uicontrol('Parent', window,'Style', 'text','String', 'FRO:','FontSize', 9,...
                    'Units', 'normalized','Position', [0.51 0.24 0.043 0.021],'Backgroundcolor',[1, 0.8, 0.8]);
    fixorient_checkbox = uicontrol('Parent', window,'Style', 'checkbox','FontSize', 30,...
                    'Backgroundcolor',[1, 0.8, 0.8],'Value',1,'Units', 'normalized','Position', [0.553 0.24 0.02 0.021], 'Callback', @plotwindow); 
    bgcorr_label = uicontrol('Parent', window,'Style', 'text','String', 'Backgr:','FontSize', 9,...
                    'Units', 'normalized','Position', [0.656 0.24 0.043 0.021],'Backgroundcolor',[0.95 0.95 0.9]);
    bgcorr_checkbox = uicontrol('Parent', window,'Style', 'checkbox','FontSize', 30,...
                    'Backgroundcolor',[0.92 0.92 0.84],'Value',1,'Units', 'normalized','Position', [0.699 0.24 0.02 0.021], 'Callback', @bgcorr); 
    initSpectra_label = uicontrol('Parent', window,'Style', 'text','String', 'Initial:','FontSize', 9,...
                    'Units', 'normalized','Position', [0.729 0.24 0.043 0.021],'Backgroundcolor',[0.95 0.95 0.9]);
    initSpectra_checkbox = uicontrol('Parent', window,'Style', 'checkbox','FontSize', 30,...
                    'Backgroundcolor',[0.92 0.92 0.84],'Value',1,'Units', 'normalized','Position', [0.772 0.24 0.02 0.021], 'Callback', @plotwindow); 
    savedSpectra_label = uicontrol('Parent', window,'Style', 'text','String', 'Saved:','FontSize', 9,...
                    'Units', 'normalized','Position', [0.802 0.24 0.043 0.021],'Backgroundcolor',[0.95 0.95 0.9]);
    savedSpectra_checkbox = uicontrol('Parent', window,'Style', 'checkbox','FontSize', 30,...
                    'Backgroundcolor',[0.92 0.92 0.84],'Value',0,'Units', 'normalized','Position', [0.845 0.24 0.02 0.021], 'Callback', @plotwindow); 
    expSpectra_label = uicontrol('Parent', window,'Style', 'text','String', 'Experim:','FontSize', 9,...
                    'Units', 'normalized','Position', [0.875 0.24 0.043 0.021],'Backgroundcolor',[0.95 0.95 0.9]);
    Freq_button = uicontrol('Parent', window,'Style', 'pushbutton','String', 'Freq','FontSize', 9,...
                    'Units', 'normalized','Position', [0.948 0.24 0.043 0.025],'Backgroundcolor',[0.95 0.95 0.9], 'Callback', @frequency_change);
    expSpectra_checkbox = uicontrol('Parent', window,'Style', 'checkbox','FontSize', 30,...
                    'Backgroundcolor',[0.92 0.92 0.84],'Value',1,'Units', 'normalized','Position', [0.918 0.24 0.02 0.021], 'Callback', @plotwindow); 
                
    timeX.label = uicontrol('Parent', window,'Style', 'text','String', 't(mus):','FontSize', 8,...
                    'Units', 'normalized','Position', [0.0 0.30 0.05 0.03]);   
    timeQ.label = uicontrol('Parent', window,'Style', 'text','String', 't(mus):','FontSize', 8,...
                    'Units', 'normalized','Position', [0.50 0.30 0.05 0.03]);       
    Parameters_label = uicontrol('Parent', window,'Style', 'listbox',...
                    'FontSize', 8,'Units', 'normalized','Position', [0.002 0.001 0.498 0.065], 'Backgroundcolor', [0.92 0.92 0.84],'Value', [], 'max',2,...
     'min',0);                         
                         
    Dat.loadlist = uicontrol('Style','listbox','Units','normalized','Position',[0.77 0.06 0.23 0.15],'Callback',@plotwindow);
    Loadtmp.filebutton = uicontrol('Parent', window, 'Style', 'pushbutton','Units', 'normalized',...
        'Position', [0.77 0.03 0.11 0.025],'FontSize', 9,'String', 'Load selected parameters', 'Callback', @loadtmpData);
    Deletetmp.filebutton = uicontrol('Parent', window, 'Style', 'pushbutton','Units', 'normalized',...
        'Position', [0.89 0.03 0.11 0.025],'FontSize', 9,'String', 'Delete selected parameters', 'Callback', @deletetmpData);
    Save.filebutton = uicontrol('Parent', window, 'Style', 'pushbutton','Units', 'normalized',...
        'Position', [0.9 0.212 0.1 0.025],'Backgroundcolor', 'yellow','FontSize', 9,'String', 'Save', 'Callback', @saveData);
    Savetmp.filebutton = uicontrol('Parent', window, 'Style', 'pushbutton','Units', 'normalized',...
        'Position', [0.795 0.212 0.1 0.025],'FontSize', 9,'String', 'Save temporaly', 'Callback', @savetmpData);
    Find.button = uicontrol('Parent', window, 'Style', 'pushbutton','Units', 'normalized',...
        'Position', [0.45 0.212 0.04 0.025],'FontSize', 9,'String', 'fit D and J', 'Callback', @findfpDJ);
    StartFP.button = uicontrol('Parent', window, 'Style', 'pushbutton','Units', 'normalized',...
        'Position', [0.45 0.185 0.04 0.025],'FontSize', 9,'String', 'Start simulation', 'Callback', @updateGraphFP2);
    Bgcorrbutton.button = uicontrol('Parent', window, 'Style', 'pushbutton','Units', 'normalized',...
        'Position', [0.51 0.212 0.09 0.025],'FontSize', 9,'String', 'Background corrections', 'Callback', @bgcorrregime);
    Advanced.button = uicontrol('Parent', window, 'Style', 'pushbutton','Units', 'normalized',...
        'Position', [0.61 0.212 0.09 0.025],'FontSize', 9,'String', 'Advanced settings', 'Callback', @advancedregime);
    Hide.button = uicontrol('Parent', window, 'Style', 'pushbutton','Units', 'normalized',...
        'Position', [0.71 0.212 0.05 0.025],'FontSize', 9,'String', 'Hide', 'Callback', @regularregime);

    %% initializing
    openData();
    if Flagfinish
        return
    end
    openXband();
    if Flagfinish
        return
    end
    openQband();      
    if Flagfinish
        return
    end
    
    np(0,0,0);
    npFP(0,0,0);
    simulationXfp(); 
    simulationQfp(); 
    simulationX(); 
    simulationQ(); 
    updatetimeQ();
    updatetimeX();
    Dat.list = {'last'};
    savetmpData(0,0,'last');
    bgcorrregime();
    set(Parameters_label, 'String',...
            sprintf('Sys1:g=%1.4f %1.4f %1.4f gStrain=%1.4f %1.4f %1.4f\n lwX=%1.4f lwQ=%1.4f \nSys2:g=%1.4f %1.4f %1.4f gStrain= %1.4f %1.4f %1.4f \nlwX=%1.4f  lwQ=%1.4f \nmwFreqX=%3.4f  mwFreqQ=%3.4f\n', ...
            data.X.Sys1.g, data.X.Sys1.gStrain, data.X.Sys1.lw, data.Q.Sys1.lw, data.X.Sys2.g, data.X.Sys2.gStrain,data.X.Sys2.lw,data.Q.Sys2.lw, X.Exp.mwFreq, Q.Exp.mwFreq));    

%     saveas(window,'window.bmp');    
        
    %% FUNCTIONS
    
    %% change of J,D, normalization
    function updateGraph1(~, ~)
        savetmpData(0,0,'last');
        DS.J = J.slider.Value;
        J.edit.String = DS.J;
        if DS.J > maxJ
            maxJ = DS.J;
            set(J.slider, 'Max', DS.J);
        elseif DS.J < minJ
            minJ = DS.J;
            set(J.slider, 'Min', DS.J);
        end
        DS.D = D.slider.Value;
        D.edit.String = DS.D;
        if DS.D > maxD
            maxD = DS.D;
            set(D.slider, 'Max', DS.D);
        elseif DS.D < minD
            minD = DS.D;
            set(D.slider, 'Min', DS.D);
        end
        
        for i=1:6
            DS.a(i).a = Dat.aslider(i).Value;
            Dat.aedit(i).String = DS.a(i).a;
        end
        
        Q.euler1 = [DS.a(1).a DS.a(2).a DS.a(3).a];  
        Q.euler2 = [DS.a(4).a DS.a(5).a DS.a(6).a];
        
        if ~isempty(findobj('type','figure','Name','Relative orientation'))
            draworient();
        end
        
        simulationX();
        simulationQ();

        if fixorient_checkbox.Value
            plotwindow();
        end
        plotwindow2fixed();
    end

    function updateGraph1Norm(~, ~)
        DS.normX = Xband.normX.slider.Value; 
        DS.normQ = Qband.normQ.slider.Value;
        Xband.normX.edit.String = DS.normX; 
        Qband.normQ.edit.String = DS.normQ;

        ydataXtotnorm = ydataXtot/max(ydataXtot) * DS.normX;
        ydataQtotnorm = ydataQtot/max(ydataQtot) * DS.normQ;
        
        if fixorient_checkbox.Value
            plotwindow();
        end
    end

    function updateGraphFP(~, ~)
        savetmpData(0,0,'last');
        DS.JFP = JFP.slider.Value;
        DS.DFP = DFP.slider.Value;
        JFP.edit.String=DS.JFP;
        DFP.edit.String=DS.DFP;
        if DS.JFP > maxJFP
            maxJFP = DS.JFP;
            set(JFP.slider, 'Max', DS.JFP);
        elseif DS.JFP < minJFP
            minJFP = DS.JFP;
            set(JFP.slider, 'Min', DS.JFP);
        end
        
        if DS.DFP > maxDFP
            maxDFP = DS.DFP;
            set(DFP.slider, 'Max', DS.DFP);
        elseif DS.DFP < minDFP
            minDFP = DS.DFP;
            set(DFP.slider, 'Min', DS.DFP);
        end
        
        if DS.DFP>0
            Dat.DFPrange.String = sprintf('%2.3f:%2.3f', DS.DFP/(2-(DS.DFP/(DS.DFP+DS.DFPw/2))^(1/3))^3, DS.DFP+DS.DFPw/2);
        else
            Dat.DFPrange.String = sprintf('%2.3f:%2.3f', DS.DFP-DS.DFPw/2, DS.DFP/(2-(DS.DFP/(DS.DFP-DS.DFPw/2))^(1/3))^3);
        end
        if DS.JFP>0
            Dat.JFPrange.String = sprintf('%2.3f:%2.3f', DS.JFP^2/(DS.JFP+DS.JFPw/2), DS.JFP+DS.JFPw/2);
        else
            Dat.JFPrange.String = sprintf('%2.3f:%2.3f', DS.JFP-DS.JFPw/2, DS.JFP^2/(DS.JFP-DS.JFPw/2));
        end
            
        if (DS.JFPw~=0 || DS.DFPw~=0) && StartFP.button.Value == 1
               if ishandle(h)
                  delete (h);
               end
               h = msgbox('Wait'); 
               SiJ = 1/DS.JFPn; 
               SiD = 1/DS.DFPn; 

               if JFP.slider.Value >0
                    normJ = (JFP.slider.Value/(JFP.slider.Value+DS.JFPw/2));
               else
                    normJ = (JFP.slider.Value/(JFP.slider.Value-DS.JFPw/2));
               end
               
               if normJ<0
                   error('J cannot change sign!');
               end
               if DFP.slider.Value >0
                   normD = (DFP.slider.Value/(DFP.slider.Value+DS.DFPw/2))^(1/3);
               else
                   normD = (DFP.slider.Value/(DFP.slider.Value-DS.DFPw/2))^(1/3);
               end
               if normD<0
                   error('D cannot change sign or be zero!')
               end
               aX = zeros(1, length(ydataXtotfpnorm));
               aQ = zeros(1, length(ydataQtotfpnorm));
               k=1;
               k2=1;
               for i=0:(DS.JFPn-1)/2
                 if  i~=0
                     DS.JFP = JFP.slider.Value*normJ^(-norminv(i*SiJ + 0.5)/norminv(1-SiJ/2));
                 else
                     DS.JFP = JFP.slider.Value;
                 end
                 JFPslide(k) = DS.JFP; k=k+1;
                 for j=0:(DS.DFPn-1)/2
                     if  j~=0
                         DS.DFP = DFP.slider.Value/(1+norminv(j*SiD + 0.5)/norminv(1-SiD/2)*(normD-1))^3;
                         DFPslide(k2) = DS.DFP; k2=k2+1;
                         simulationXfp();
                         simulationQfp();
                         aX = aX + ydataXtotfpnorm;
                         aQ = aQ + ydataQtotfpnorm;
                         DS.DFP = DFP.slider.Value/(1+norminv(j*SiD + 0.5)/norminv(1-SiD/2)*(1-normD))^3;
                     else
                         DS.DFP = DFP.slider.Value;
                     end
                     DFPslide(k2) = DS.DFP; k2=k2+1;
                     simulationXfp();
                     simulationQfp();
                     aX = aX + ydataXtotfpnorm;
                     aQ = aQ + ydataQtotfpnorm;
                 end
                 if  i~=0
                     DS.JFP = JFP.slider.Value*normJ^(norminv(i*SiJ + 0.5)/norminv(1-SiJ/2));
                     JFPslide(k) = DS.JFP; k=k+1;
                     for j=0:(DS.DFPn-1)/2
                         if  j~=0
                             DS.DFP = DFP.slider.Value/(1+norminv(j*SiD + 0.5)/norminv(1-SiD/2)*(normD-1))^3;
                             DFPslide(k2) = DS.DFP; k2=k2+1;
                             simulationXfp();
                             simulationQfp();
                             aX = aX + ydataXtotfpnorm;
                             aQ = aQ + ydataQtotfpnorm;
                             DS.DFP = DFP.slider.Value/(1+norminv(j*SiD + 0.5)/norminv(1-SiD/2)*(1-normD))^3;
                         else
                             DS.DFP = DFP.slider.Value;
                         end
                         DFPslide(k2) = DS.DFP; k2=k2+1;
                         simulationXfp();
                         simulationQfp();
                         aX = aX + ydataXtotfpnorm;
                         aQ = aQ + ydataQtotfpnorm;
                     end
                 end
               end
               figure(4)
               subplot(2,1,1)
               plot(DFPslide, zeros(1, length(DFPslide)),'.')
               legend('D')
               subplot(2,1,2)
               plot(JFPslide, zeros(1, length(JFPslide)),'.')
               legend('J')
               ydataXtotfp = aX;
               ydataXtotfpnorm = ydataXtotfp/max(ydataXtotfp) * DS.normXFP;
               ydataQtotfp = aQ;
               ydataQtotfpnorm = ydataQtotfp/max(ydataQtotfp) * DS.normQFP;
               DS.JFP = JFP.slider.Value;
               DS.DFP = DFP.slider.Value;
               delete(h);
               h = msgbox('Done');
        elseif (DS.JFPw~=0 || DS.DFPw~=0)
           if ishandle(h)
                delete (h);
           end
           h = msgbox('Press "start simulation" to start simulation');
        else
           simulationXfp();
           simulationQfp();  
        end


        if fullpowder_checkbox.Value    
            plotwindow();
        end
        plotwindow2fp();
    end

    function updateGraphFPNorm(~, ~)
        DS.normXFP = Xband.normXFP.slider.Value;
        DS.normQFP = Qband.normQFP.slider.Value;
        Xband.normXFP.edit.String=DS.normXFP;
        Qband.normQFP.edit.String=DS.normQFP;
        
        ydataXtotfpnorm = ydataXtotfp/max(ydataXtotfp) * DS.normXFP;
        ydataQtotfpnorm = ydataQtotfp/max(ydataQtotfp) * DS.normQFP;
        
        if fullpowder_checkbox.Value    
            plotwindow();
        end
    end

    function updateGraph2(~, ~)
        DS.J= str2double(J.edit.String);
        J.slider.Value=DS.J;
        DS.D= str2double(D.edit.String);
        D.slider.Value=DS.D;
        DS.E= str2double(E.edit.String);
        for i=1:6
            DS.a(i).a = str2double(Dat.aedit(i).String);
            Dat.aslider(i).Value = DS.a(i).a;
        end
        
        updateGraph1();
    end

    function updateGraph2Norm(~, ~)
        DS.normX = str2double(Xband.normX.edit.String);
        DS.normQ = str2double(Qband.normQ.edit.String);
        Xband.normX.slider.Value = DS.normX;
        Qband.normQ.slider.Value = DS.normQ;

        updateGraph1Norm();
    end

    function updateGraphFP2Norm(~, ~)
        DS.normXFP = str2double(Xband.normXFP.edit.String);
        DS.normQFP = str2double(Qband.normQFP.edit.String);
        Xband.normXFP.slider.Value = DS.normXFP;
        Qband.normQFP.slider.Value = DS.normQFP;
        
        updateGraphFPNorm();
    end

    function updateGraphFP2(~, ~)
        DS.JFP = str2double(JFP.edit.String);
        DS.DFP = str2double(DFP.edit.String);
        JFP.slider.Value =DS.JFP;
        DFP.slider.Value = DS.DFP;
        if DS.EFP == 0 && str2double(EFP.edit.String)~=0
            [Dat.DOrientationFP.phidip,Dat.DOrientationFP.thetadip,Dat.DOrientationFP.Weight] = sphgrid('Th', DS.npFP);
        end
        DS.EFP = str2double(EFP.edit.String);
        DS.JFPn = round((str2double(Dat.JFPnumber.edit.String)-1)/2)*2+1;
        DS.DFPn = round((str2double(Dat.DFPnumber.edit.String)-1)/2)*2+1;
        Dat.JFPnumber.edit.String = DS.JFPn;
        Dat.DFPnumber.edit.String = DS.DFPn;
        DS.JFPw = str2double(Dat.JFPwidth.edit.String);
        DS.DFPw = str2double(Dat.DFPwidth.edit.String);
        
        updateGraphFP();
    end
    
%% residual
    function residual(~, ~)
            if bgcorr_checkbox.Value
                Xbgcorr();
                Qbgcorr();
            end

            ind = [find(X.B0<=str2double(Xband.B1.edit.String)), find(X.B0>=str2double(Xband.B2.edit.String))];
            smXtemp = X.smooth(setdiff(1:length(X.B0),ind),DS.sliceX);
            ind = [find(Q.B0<=str2double(Qband.B1.edit.String)), find(Q.B0>=str2double(Qband.B2.edit.String))];
            smQtemp = Q.smooth(setdiff(1:length(Q.B0),ind),DS.sliceQ);
        
            ind = [find(xdataXfp<=str2double(Xband.B1.edit.String)), find(xdataXfp>=str2double(Xband.B2.edit.String))];
            ydataXtemp = ydataXtotfpnorm(setdiff(1:length(ydataXtotfp),ind));
            ind = [find(xdataQfp<=str2double(Qband.B1.edit.String)), find(xdataQfp>=str2double(Qband.B2.edit.String))];
            ydataQtemp = ydataQtotfpnorm(setdiff(1:length(ydataQtotfp),ind));
            resnormX = sum((ydataXtemp-smXtemp').^2);
            resnormQ = sum((ydataQtemp-smQtemp').^2);
            if fullpowder_checkbox.Value
                R2X = 1-resnormX/sum((smXtemp-mean(smXtemp)).^2);  
                R2Q = 1-resnormQ/sum((smQtemp-mean(smQtemp)).^2);
            end
            R2fp = 1-(resnormX+resnormQ)/(sum((smXtemp-mean(smXtemp)).^2)+sum((smQtemp-mean(smQtemp)).^2)); 
            ind = [find(xdataX<=str2double(Xband.B1.edit.String)), find(xdataX>=str2double(Xband.B2.edit.String))];
            ydataXtemp = ydataXtotnorm(setdiff(1:length(ydataXtot),ind));
            ind = [find(xdataQ<=str2double(Qband.B1.edit.String)), find(xdataQ>=str2double(Qband.B2.edit.String))];
            ydataQtemp = ydataQtotnorm(setdiff(1:length(ydataQtot),ind));
            resnormX = sum((ydataXtemp-smXtemp').^2);
            resnormQ = sum((ydataQtemp-smQtemp').^2);
            if fixorient_checkbox.Value
                R2X = 1-resnormX/sum((smXtemp-mean(smXtemp)).^2);  
                R2Q = 1-resnormQ/sum((smQtemp-mean(smQtemp)).^2);
            end
            R2fo = 1-(resnormX+resnormQ)/(sum((smXtemp-mean(smXtemp)).^2)+sum((smQtemp-mean(smQtemp)).^2)); 
            if ~fullpowder_checkbox.Value && ~fixorient_checkbox.Value
                R2X = 0; R2Q = 0;
            end
            set(R2fp_edit, 'String', R2fp);
            set(R2fo_edit, 'String', R2fo);
            set(R2X_edit, 'String', R2X);
            set(R2Q_edit, 'String', R2Q);
    end

%% addbr, time is changed
    function addbroadX (~, ~)
        savetmpData(0,0,'last');
        DS.addbrX = str2double(Dat.addbrX.edit.String);
        simulationXfp();
        simulationX();
        plotwindow();
    end

    function addbroadQ (~, ~)
        savetmpData(0,0,'last');
        DS.addbrQ = str2double(Dat.addbrQ.edit.String);
        simulationQfp();
        simulationQ();
        plotwindow();
    end

    function updatetimeX (~, ~)
        savetmpData(0,0,'last');
%         index = find(Xband.time.slider.Value >= X.t);
%         DS.sliceX= index(end); %,'last'
        [~, DS.sliceX] = min(abs(X.t-Xband.time.slider.Value));
        DS.tX = X.t(DS.sliceX);
        Xband.time.edit.String = DS.tX;
        Xband.time.slider.Value = DS.tX;
        
        if isfield(Dat, 'I')
            if Data(Dat.I).tX == DS.tX
                DS.fitlevelX = Data(Dat.I).fitlevelX;
                Xband.fitlevel.edit.String =  num2str(DS.fitlevelX);
                DS.smoothX = Data(Dat.I).smoothX;
                Xband.smooth.edit.String =  num2str(DS.smoothX);
                Xband.B1.edit.String = Data(Dat.I).B1X;
                Xband.B2.edit.String = Data(Dat.I).B2X;
            end
        end
        
        Xbgcorr();
        plotwindow();
    end

    function updatetime2X (~, ~)
        Xband.time.slider.Value = str2double(Xband.time.edit.String);
        updatetimeX();
    end 
 
    function updatetimeQ (~, ~)
        savetmpData(0,0,'last');
%         index = find(Qband.time.slider.Value >= Q.t);
%         DS.sliceQ= index(end); %,'last'
        [~, DS.sliceQ] = min(abs(Q.t-Qband.time.slider.Value));
        DS.tQ = Q.t(DS.sliceQ);
        Qband.time.edit.String = Q.t(DS.sliceQ);
        Qband.time.slider.Value = Q.t(DS.sliceQ);
            
         if isfield(Dat, 'I')
            if Data(Dat.I).tQ == DS.tQ
                DS.fitlevelQ =  Data(Dat.I).fitlevelQ;
                Qband.fitlevel.edit.String =  num2str(DS.fitlevelQ);
                DS.smoothQ =  Data(Dat.I).smoothQ;
                Qband.smooth.edit.String =  num2str(DS.smoothQ);
                Qband.B1.edit.String =  Data(Dat.I).B1Q;
                Qband.B2.edit.String =  Data(Dat.I).B2Q;
            end
        end
        
        Qbgcorr();
        plotwindow();
    end

    function updatetime2Q (~, ~)
        Qband.time.slider.Value = str2double(Qband.time.edit.String);
        updatetimeQ();
    end

    function np(~,~, fl)
        savetmpData(0,0,'last');
        DS.np = str2double(Dat.np.edit.String);
        [Dat.gOrientation.S1.phi,Dat.gOrientation.S1.theta,Dat.gOrientation.S1.Weight] = sphgridpart_Ci(pi/2, 0, 0);
        [Dat.gOrientation.S2.phi,Dat.gOrientation.S2.theta,Dat.gOrientation.S2.Weight] = sphgridpart_Ci(pi/2, 0, 0);
        [Dat.DOrientation.phidip,Dat.DOrientation.thetadip,Dat.DOrientation.Weight] = sphgrid('Ci', DS.np);
         
        if fl~=0
            simulationX();
            simulationQ();
            plotwindow();
        end
    end

    function npFP(~,~, fl)
        savetmpData(0,0,'last');
        DS.npFP = str2double(Dat.npFP.edit.String);  
        [Dat.gOrientationFP.S1.phi,Dat.gOrientationFP.S1.theta,Dat.gOrientationFP.S1.Weight] =  sphgrid('Th', DS.npFP);
        if DS.EFP ==0
            [Dat.DOrientationFP.phidip,Dat.DOrientationFP.thetadip,Dat.DOrientationFP.Weight] = sphgrid('Dinfh', DS.npFP);
        else
            [Dat.DOrientationFP.phidip,Dat.DOrientationFP.thetadip,Dat.DOrientationFP.Weight] = sphgrid('Th', DS.npFP);
        end
        if fl~=0 
            simulationXfp();
            simulationQfp();
            plotwindow();
        end
    end


%% bachground corrections
    function bgcorr(~, ~)
        if bgcorr_checkbox.Value
            bgcorrregime();
            Xbgcorr();
            Qbgcorr();
        else
            regularregime();
        end
        plotwindow();
    end

    function Xbgcorr(~, ~)
            DS.fitlevelX = round(str2double(Xband.fitlevel.edit.String));
            DS.smoothX = round(str2double(Xband.smooth.edit.String));

            Xband.fitlevel.edit.String = DS.fitlevelX;
            Xband.smooth.edit.String = DS.smoothX;
            
            fit_begin = 1;
            fit_end = length(X.B0);    %% Fit Bereich bis zu den letzten Werten
            [~, Xband.fit_range] = min(abs(X.B0 - str2double(Xband.B1.edit.String)));
            [~, Xband.fit_range2] = min(abs(X.B0 - str2double(Xband.B2.edit.String)));
            
            X.s(:,DS.sliceX) = X.s(:,DS.sliceX)/max(X.s(:,DS.sliceX));
            B_exp_fit = [X.B0(fit_begin:Xband.fit_range),X.B0(Xband.fit_range2:fit_end)];
            S_exp_fit = [X.s(fit_begin:Xband.fit_range,DS.sliceX);X.s(Xband.fit_range2:fit_end,DS.sliceX)]';
            fit1 = polyfit(B_exp_fit,S_exp_fit, DS.fitlevelX);
            
            Xband.S_exp_background = polyval(fit1,X.B0);
            
            X.sFit(:,DS.sliceX) = X.s(:,DS.sliceX)-Xband.S_exp_background';
            X.smooth(:,DS.sliceX) = datasmooth(X.sFit(:,DS.sliceX),DS.smoothX,'savgol');
     end

    function Qbgcorr(~, ~)
            DS.fitlevelQ = round(str2double(Qband.fitlevel.edit.String));
            DS.smoothQ = round(str2double(Qband.smooth.edit.String));

            Qband.fitlevel.edit.String = DS.fitlevelQ;
            Qband.smooth.edit.String = DS.smoothQ;
            
            fit_begin = 1;
            fit_end = length(Q.B0);    %% Fit Bereich bis zu den letzten Werten
            [~, Qband.fit_range] = min(abs(Q.B0 - str2double(Qband.B1.edit.String)));
            [~, Qband.fit_range2] = min(abs(Q.B0 - str2double(Qband.B2.edit.String)));
            
            Q.s(:,DS.sliceQ) = Q.s(:,DS.sliceQ)/ max(Q.s(:,DS.sliceQ));
            B_exp_fit = [Q.B0(fit_begin:Qband.fit_range),Q.B0(Qband.fit_range2:fit_end)];
            S_exp_fit = [Q.s(fit_begin:Qband.fit_range,DS.sliceQ);Q.s(Qband.fit_range2:fit_end,DS.sliceQ)]';
            fit1 = polyfit(B_exp_fit,S_exp_fit,DS.fitlevelQ);
            
            Qband.S_exp_background = polyval(fit1,Q.B0);
            
            Q.sFit(:,DS.sliceQ) = Q.s(:,DS.sliceQ)-Qband.S_exp_background';
            Q.smooth(:,DS.sliceQ) = datasmooth(Q.s(:,DS.sliceQ)-Qband.S_exp_background',DS.smoothQ,'savgol');
    end

    function Xbgcorr2(~, ~)
        Xbgcorr();
        plotwindow();
    end

    function Qbgcorr2(~, ~)
        Qbgcorr();
        plotwindow();
    end

%% PLOT
    function plotwindow(~, ~)
        cla(Qband.plot);  hold (Qband.plot,'on');
        grid(Qband.plot,'on')
        cla(Xband.plot);  hold (Xband.plot,'on');
        grid(Xband.plot,'on')
        residual();
        lstr = 1;
        legendstr = string(6);
        if initSpectra_checkbox.Value
            plot(Qband.plot,Q.B0, Q.s(:,DS.sliceQ), 'Color',[0.5 0.5 0.5],'LineWidth', 0.3); 
            plot(Xband.plot,X.B0, X.s(:,DS.sliceX),  'Color',[0.5 0.5 0.5], 'LineWidth', 0.3); 
            legendstr(lstr)= 'experiment(without corrections)'; lstr = lstr+1;
        end
        
        if bgcorr_checkbox.Value
            [~, Qband.fit_range] = min(abs(Q.B0 - str2double(Qband.B1.edit.String)));
            [~, Qband.fit_range2] = min(abs(Q.B0 - str2double(Qband.B2.edit.String)));
            [~, Xband.fit_range] = min(abs(X.B0 - str2double(Xband.B1.edit.String)));
            [~, Xband.fit_range2] = min(abs(X.B0 - str2double(Xband.B2.edit.String)));
           
            plot(Qband.plot,Q.B0, Qband.S_exp_background,'-', 'Color', [0.7  0.7  0.7],'LineWidth', 1.0); 
            plot(Xband.plot,X.B0, Xband.S_exp_background,'-', 'Color', [0.7  0.7  0.7],'LineWidth', 1.0); 
            legendstr(lstr)= 'background'; lstr = lstr+1;
            if Dat.bgwithoutsmooth_checkbox.Value
                plot(Qband.plot,Q.B0, Q.sFit(:,DS.sliceQ),'-', 'Color',[0.5 0.5 0.5 0.5],'LineWidth', 0.1); 
                plot(Xband.plot,X.B0, X.sFit(:,DS.sliceX),'-',  'Color',[0.5 0.5 0.5 0.5],'LineWidth', 0.1);  
                legendstr(lstr)= 'experiment(bg corrected)'; lstr = lstr+1;
            end
        end
        
        if expSpectra_checkbox.Value
            plot(Qband.plot, Q.B0, Q.smooth(:,DS.sliceQ),'k','LineWidth', 1.5);
            plot(Xband.plot, X.B0, X.smooth(:,DS.sliceX),'k','LineWidth', 1.5);
        end
        
        if  bgcorr_checkbox.Value
            legendstr(lstr)= 'experiment(bg corrected and smoothed)';lstr = lstr+1;
        else
            legendstr(lstr)= 'experiment(corrected)';lstr = lstr+1;
        end
        
        ind = Dat.loadlist.Value;
        if savedSpectra_checkbox.Value
            if fixorient_checkbox.Value
                plot(Qband.plot, Q.B0, (Dat.listdata{Dat.loadlist.Value}.DStmpQ.fix*DS.normQ),'--','Color','r','LineWidth', 1.5);
                plot(Xband.plot, X.B0, (Dat.listdata{Dat.loadlist.Value}.DStmpX.fix*DS.normX),'--','Color','r','LineWidth', 1.5);
                legendstr(lstr)= [Dat.loadlist.String{ind},' FRO'];lstr = lstr+1;
            end
            if fullpowder_checkbox.Value
                plot(Qband.plot, Q.B0, (Dat.listdata{Dat.loadlist.Value}.DStmpQ.fp*DS.normQFP),'--','Color','b','LineWidth', 1.5);
                plot(Xband.plot, X.B0, (Dat.listdata{Dat.loadlist.Value}.DStmpX.fp*DS.normXFP),'--','Color','b','LineWidth', 1.5);
                legendstr(lstr)= [Dat.loadlist.String{ind},' FP'];lstr = lstr+1;
            end
        end
        
        if fixorient_checkbox.Value
            plot(Qband.plot, xdataQ, ydataQtotnorm,'r','LineWidth', 1.5);
            plot(Xband.plot, xdataX, ydataXtotnorm,'r','LineWidth', 1.5);
            legendstr(lstr)= 'FRO simulation';lstr = lstr+1;
        end
        
        if fullpowder_checkbox.Value
            plot(Qband.plot, xdataQfp, ydataQtotfpnorm,'b', 'LineWidth', 1.5);
            plot(Xband.plot, xdataXfp, ydataXtotfpnorm,'b', 'LineWidth', 1.5);
            legendstr(lstr)= 'FP simulation'; lstr = lstr+1;
        end
       
        xlabel(Qband.plot, 'Magnetic field (mT)');
        ylabel(Qband.plot, 'Intensity (a. u.)');
        title(Qband.plot, 'Qband');
 
        xlabel(Xband.plot, 'Magnetic field (mT)');
        ylabel(Xband.plot, 'Intensity (a. u.)');
        title(Xband.plot, 'Xband');
        
        axes(Xband.plot)        
        legend(Xband.plot, legendstr); legend boxoff;
        axes(Qband.plot)
        legend(Qband.plot, legendstr); legend boxoff;
        
        
        if bgcorr_checkbox.Value
            axes(Qband.plot)
            ylQ = ylim;
            p1Q=fill([Q.B0(1) Q.B0(Qband.fit_range) Q.B0(Qband.fit_range) Q.B0(1)],[ylQ(1) ylQ(1) ylQ(2) ylQ(2)],[0.9 0.9 0.9],'LineStyle','none'); hold on
            p2Q=fill([Q.B0(Qband.fit_range2) Q.B0(end) Q.B0(end) Q.B0(Qband.fit_range2)],[ylQ(1) ylQ(1) ylQ(2) ylQ(2)],[0.9 0.9 0.9],'LineStyle','none'); hold on
            set(p1Q,'FaceAlpha',0.4);  set(p2Q,'FaceAlpha',0.4); 

            axes(Xband.plot)
            ylX = ylim;
            p1X=fill([X.B0(1) X.B0(Xband.fit_range) X.B0(Xband.fit_range) X.B0(1)],[ylX(1) ylX(1) ylX(2) ylX(2)],[0.9 0.9 0.9],'LineStyle','none'); hold on
            p2X=fill([X.B0(Xband.fit_range2) X.B0(end) X.B0(end) X.B0(Xband.fit_range2)],[ylX(1) ylX(1) ylX(2) ylX(2)],[0.9 0.9 0.9],'LineStyle','none'); hold on
            set(p1X,'FaceAlpha',0.4);  set(p2X,'FaceAlpha',0.4);
        end
        plotwindow2fixed();
        plotwindow2fp();
    end

    function plotwindow2fixed()
        cla(Qband2);  hold (Qband2,'on');
        cla(Xband2);  hold (Xband2,'on');
        for l = 1:4
%             plot(Qband2,Q.B0,ydataQ{l}/max(ydataQ{1}),'b','LineStyle','-');
%             plot(Xband2,X.B0,ydataX{l}/max(ydataX{1}),'b','LineStyle','-');
            plot(Qband2,Q.B0,ydataQ{l},'k','LineStyle','-');
            plot(Xband2,X.B0,ydataX{l},'k','LineStyle','-');
        end 
%         plot(Qband2,Q.B0,ydataQtot/max(ydataQtot),'LineWidth',1.5,'Color','k');        
%         plot(Xband2,X.B0,ydataXtot/max(ydataXtot),'LineWidth',1.5,'Color','k');  
                plot(Qband2,Q.B0,ydataQtot,'LineWidth',1.5,'Color','r');        
        plot(Xband2,X.B0,ydataXtot,'LineWidth',1.5,'Color','r');  
        
        if savedSpectra_checkbox.Value
            for l = 1:4
%                 plot(Qband2,Q.B0,Dat.listdata{Dat.loadlist.Value}.DStmpQ.fixtemp{l}/max(Dat.listdata{Dat.loadlist.Value}.DStmpQ.fixtemp{1}),'b','LineStyle','--');
%                 plot(Xband2,X.B0,Dat.listdata{Dat.loadlist.Value}.DStmpX.fixtemp{l}/max(Dat.listdata{Dat.loadlist.Value}.DStmpX.fixtemp{1}),'b','LineStyle','--');
                plot(Qband2,Q.B0,Dat.listdata{Dat.loadlist.Value}.DStmpQ.fixtemp{l},'k','LineStyle','--');
                plot(Xband2,X.B0,Dat.listdata{Dat.loadlist.Value}.DStmpX.fixtemp{l},'k','LineStyle','--');
            end 
%             plot(Qband2, Q.B0, Dat.listdata{Dat.loadlist.Value}.DStmpQ.fix,'--','Color','k');
%             plot(Xband2, X.B0, Dat.listdata{Dat.loadlist.Value}.DStmpX.fix,'--','Color','k');
             plot(Qband2, Q.B0, Dat.listdata{Dat.loadlist.Value}.DStmpQ.fixtemp{1}+Dat.listdata{Dat.loadlist.Value}.DStmpQ.fixtemp{2}+Dat.listdata{Dat.loadlist.Value}.DStmpQ.fixtemp{3}+Dat.listdata{Dat.loadlist.Value}.DStmpQ.fixtemp{4},'--','Color','r');
            plot(Xband2, X.B0, Dat.listdata{Dat.loadlist.Value}.DStmpX.fixtemp{1}+Dat.listdata{Dat.loadlist.Value}.DStmpX.fixtemp{2}+Dat.listdata{Dat.loadlist.Value}.DStmpX.fixtemp{3}+Dat.listdata{Dat.loadlist.Value}.DStmpX.fixtemp{4},'--','Color','r');
        end
        
        title(Qband2, 'Qband fixed relative orientation in power');
        title(Xband2, 'Xband fixed relative orientation in power');
        xlabel(Qband2, 'Magnetic field (mT)');
        ylabel(Qband2, 'Intensity (a. u.)');
        xlabel(Xband2, 'Magnetic field (mT)');
        ylabel(Xband2, 'Intensity (a. u.)');
        
    end

    function plotwindow2fp() 
        cla(Qband2FP);  hold (Qband2FP,'on');
        cla(Xband2FP);  hold (Xband2FP,'on');
        for l = 1:4
%             plot(Qband2FP,Q.B0,ydataQfp{l}/max(ydataQfp{1}),'b');
%             plot(Xband2FP,X.B0,ydataXfp{l}/max(ydataXfp{1}),'b');
            plot(Qband2FP,Q.B0,ydataQfp{l},'k');
            plot(Xband2FP,X.B0,ydataXfp{l},'k');
        end 
%         plot(Qband2FP,Q.B0,ydataQtotfp/max(ydataQtotfp),'LineWidth',1.5,'Color','k');        
%         plot(Xband2FP,X.B0,ydataXtotfp/max(ydataXtotfp),'LineWidth',1.5,'Color','k');        
     plot(Qband2FP,Q.B0,ydataQtotfp,'LineWidth',1.5,'Color','b');        
        plot(Xband2FP,X.B0,ydataXtotfp,'LineWidth',1.5,'Color','b');        
        if savedSpectra_checkbox.Value   
            for l = 1:4
%                 plot(Qband2FP,Q.B0,Dat.listdata{Dat.loadlist.Value}.DStmpQ.fptemp{l}/max(Dat.listdata{Dat.loadlist.Value}.DStmpQ.fptemp{1}),'b','LineStyle','--');
%                 plot(Xband2FP,X.B0,Dat.listdata{Dat.loadlist.Value}.DStmpX.fptemp{l}/max(Dat.listdata{Dat.loadlist.Value}.DStmpX.fptemp{1}),'b','LineStyle','--');
                plot(Qband2FP,Q.B0,Dat.listdata{Dat.loadlist.Value}.DStmpQ.fptemp{l},'k','LineStyle','--');
                plot(Xband2FP,X.B0,Dat.listdata{Dat.loadlist.Value}.DStmpX.fptemp{l},'k','LineStyle','--');

            end 
%             plot(Qband2FP,Q.B0,Dat.listdata{Dat.loadlist.Value}.DStmpQ.fp,'k','LineStyle','--');
%             plot(Xband2FP,X.B0,Dat.listdata{Dat.loadlist.Value}.DStmpX.fp,'k','LineStyle','--');
            plot(Qband2FP,Q.B0,Dat.listdata{Dat.loadlist.Value}.DStmpQ.fptemp{1}+Dat.listdata{Dat.loadlist.Value}.DStmpQ.fptemp{2}+Dat.listdata{Dat.loadlist.Value}.DStmpQ.fptemp{3}+Dat.listdata{Dat.loadlist.Value}.DStmpQ.fptemp{4},'b','LineStyle','--');
            plot(Xband2FP,X.B0,Dat.listdata{Dat.loadlist.Value}.DStmpX.fptemp{1}+Dat.listdata{Dat.loadlist.Value}.DStmpX.fptemp{2}+Dat.listdata{Dat.loadlist.Value}.DStmpX.fptemp{3}+Dat.listdata{Dat.loadlist.Value}.DStmpX.fptemp{4},'b','LineStyle','--');
        end
        
        title(Qband2FP, 'Qband full powder');
        title(Xband2FP, 'Xband full powder');
        xlabel(Qband2FP, 'Magnetic field (mT)');
        ylabel(Qband2FP, 'Intensity (a. u.)');
        xlabel(Xband2FP, 'Magnetic field (mT)');
        ylabel(Xband2FP, 'Intensity (a. u.)');
        
    end

%% simulation functions
    function simulationX()
        if DS.E == 0
            [xdataX,ydataX] = ctssimaniso_git_v3_6_1(data.X.Sys1,data.X.Sys2,X.Exp,Dat.gOrientation,Dat.DOrientation,DS.J,DS.D,DS.addbrX, Q.euler1, Q.euler2, length(X.B0)); 
        else
            [xdataX,ydataX] = ctssimaniso_git_v4(data.X.Sys1,data.X.Sys2,X.Exp,Dat.gOrientation,Dat.DOrientation,DS.J,DS.D,DS.E,DS.addbrX, Q.euler1, Q.euler2, length(X.B0));
        end
        ydataXtot = ydataX{1} + ydataX{2} + ydataX{3} + ydataX{4};
        ydataXtotnorm = ydataXtot/max(ydataXtot) * DS.normX;
    end
    function simulationQ()
        if DS.E == 0
            [xdataQ,ydataQ] = ctssimaniso_git_v3_6_1(data.Q.Sys1,data.Q.Sys2,Q.Exp,Dat.gOrientation,Dat.DOrientation,DS.J,DS.D,DS.addbrQ, Q.euler1, Q.euler2, length(Q.B0)); 
        else
            [xdataQ,ydataQ] = ctssimaniso_git_v4(data.Q.Sys1,data.Q.Sys2,Q.Exp,Dat.gOrientation,Dat.DOrientation,DS.J,DS.D,DS.E,DS.addbrQ, Q.euler1, Q.euler2, length(Q.B0)); 
        end
        ydataQtot = ydataQ{1} + ydataQ{2} + ydataQ{3} + ydataQ{4};   
        ydataQtotnorm = ydataQtot/max(ydataQtot) * DS.normQ;
    end
    function simulationXfp()
        if DS.EFP == 0
            [xdataXfp,ydataXfp] = ctssimaniso_git_initial_v2_1(data.X.Sys1,data.X.Sys2,X.Exp,Dat.gOrientationFP,Dat.DOrientationFP,DS.JFP,DS.DFP,DS.addbrX, length(X.B0)); 
        else
            [xdataXfp,ydataXfp] = ctssimaniso_git_initial_v2_4(data.X.Sys1,data.X.Sys2,X.Exp,Dat.gOrientationFP,Dat.DOrientationFP,DS.JFP,DS.DFP,DS.EFP,DS.addbrX, length(X.B0)); 
        end
        ydataXtotfp = ydataXfp{1} + ydataXfp{2} + ydataXfp{3} + ydataXfp{4};
        ydataXtotfpnorm = ydataXtotfp/max(ydataXtotfp) * DS.normXFP;
    end
    function simulationQfp()
        if DS.EFP == 0
            [xdataQfp,ydataQfp] = ctssimaniso_git_initial_v2_1(data.Q.Sys1,data.Q.Sys2,Q.Exp,Dat.gOrientationFP,Dat.DOrientationFP,DS.JFP,DS.DFP,DS.addbrQ, length(Q.B0)); 
        else
            [xdataQfp,ydataQfp] = ctssimaniso_git_initial_v2_4(data.Q.Sys1,data.Q.Sys2,Q.Exp,Dat.gOrientationFP,Dat.DOrientationFP,DS.JFP,DS.DFP,DS.EFP,DS.addbrQ, length(Q.B0)); 
        end
        ydataQtotfp = ydataQfp{1} + ydataQfp{2} + ydataQfp{3} + ydataQfp{4};
        ydataQtotfpnorm = ydataQtotfp/max(ydataQtotfp) * DS.normQFP;     
    end

    function frequency_change(~, ~)
        windowFreqX = figure('Color', [0.9255 0.9137 0.8471], 'Name', 'Parameters finder',...
                    'DockControl', 'off', 'Units', 'Pixels', 'Position', [100 50 800 600]);
        windowFreqQ = figure('Color', [0.9255 0.9137 0.8471], 'Name', 'Parameters finder',...
                    'DockControl', 'off', 'Units', 'Pixels', 'Position', [100 50 800 600]);
        if ~isfield(Xband,'freq')
            Xband.freq = X.Exp.mwFreq;
        end
        if ~isfield(Qband,'freq')
            Qband.freq = Q.Exp.mwFreq;
        end
        Xband.freq = FrequencyChange(windowFreqX, data.X, Dat, X.Exp, DS, DS.addbrX, DS.normX, DS.normXFP, X.B0, X.smooth(:,DS.sliceX), Xband.freq);
        Qband.freq = FrequencyChange(windowFreqQ, data.Q, Dat, Q.Exp, DS, DS.addbrQ, DS.normQ, DS.normQFP, Q.B0, Q.smooth(:,DS.sliceQ), Qband.freq);
    end

%% open files
    function openXband(~, ~)
        [DS.fileNameX,Xband.path,FilterIndex] = uigetfile('Xband.mat','Select the MATLAB code file for Xband');
        if FilterIndex~=0
            X = load([Xband.path,DS.fileNameX],'-mat');

            if isfield(X, 'file')
                X.B0 = X.file.B0;
                X.s = X.file.s_trcorr-(X.file.s_trcorr(1,:)+X.file.s_trcorr(length(X.file.s_trcorr(:,1)),:))/2;
                X.t = X.file.t;
                X.Exp.mwFreq = X.file.mw_freq_set;
                if X.Exp.mwFreq >10^9
                   X.Exp.mwFreq = X.Exp.mwFreq/10^9;
                end
                DS.mwFreqX = X.Exp.mwFreq;
                X.Exp.Range = [min(X.B0), max(X.B0)];
            else
               error('Cannot load the file for Xband. Wrong structure. Should have field "file" that contains B0, s,t,mw_freq_set')
           end

            X.t = round(X.t*50)/50;
            DS.sliceX = 1;
            DS.tX = (min(X.t) + max(X.t))/2;
            DS.B1X = min(X.B0);
            DS.B2X = max(X.B0);
            loadparameters(1);

            if isfield (DS, 'fileNameQ')
                Xbgcorr();
                simulationXfp();
                updateGraph1();
                updateGraphFP();
                Dat.list = {'last'};
                rmfield(Dat, 'listdata');
                savetmpData(0,0,'last');
                Dat.loadlist.Value=1;
                if ~(fullpowder_checkbox.Value || fixorient_checkbox.Value)
                    plotwindow();
                end
            end
        else
            Flagfinish = 1;
        end
    end

    function openQband(~, ~)
       [DS.fileNameQ,Qband.path,FilterIndex] = uigetfile('Qband.mat','Select the MATLAB code file for Qband');
       if FilterIndex~=0
           Q = load([Qband.path,DS.fileNameQ],'-mat');

           if isfield(Q, 'file')
                Q.B0 = Q.file.B0;
                Q.s = Q.file.s_trcorr-(Q.file.s_trcorr(1,:)+Q.file.s_trcorr(length(Q.file.s_trcorr(:,1)),:))/2;
                Q.t = Q.file.t;
                Q.Exp.mwFreq = Q.file.mw_freq_set;
                if Q.Exp.mwFreq >10^9
                    Q.Exp.mwFreq = Q.Exp.mwFreq/10^9;
                end
                DS.mwFreqQ = Q.Exp.mwFreq;
                Q.Exp.Range = [min(Q.B0), max(Q.B0)];
          else
               error('Cannot load the file for Qband. Wrong structure. Should have field "file" that contains B0, s,t,mw_freq_set')
          end

          Q.t = round(Q.t*50)/50;
          DS.sliceQ = 1;

          DS.tQ = (min(Q.t) + max(Q.t))/2;
          DS.B1Q = min(Q.B0);
          DS.B2Q = max(Q.B0);

          loadparameters(1);

          if ~isempty(xdataXfp)
            Qbgcorr();
            simulationQfp();
            updateGraph1();
            updateGraphFP();
            Dat.list = {'last'};
            rmfield(Dat, 'listdata');
            savetmpData(0,0,'last');
            Dat.loadlist.Value=1;
            if ~(fullpowder_checkbox.Value || fixorient_checkbox.Value)
                    plotwindow();
            end
          end
       else
          Flagfinish = 1;
       end
    end

    function openData(~, ~)
        choice = questdlg('Choose file for system data', ' ', 'Create/change file','Open file', 'Open file');
            switch choice
            case 'Create/change file'
                 [Dat.fileName,Dat.path, FilterIndex] = Datafilecreation(); 
            case 'Open file'
                 [Dat.fileName,Dat.path, FilterIndex] = uigetfile('System.mat','Select the MATLAB code file for system data'); 
            end
        if ~isempty(choice)
            if FilterIndex
                Data = load([Dat.path,Dat.fileName],'-mat');
                data = Data.data;
                if isfield (Data, 'simulations')
                    Data = Data.simulations;
                else
                   Data = struct([]);
                end
                Dat.filebutton = uicontrol('Parent', window, 'Style', 'pushbutton','Units', 'normalized',...
                     'Position', [0.675 0.001 0.15 0.025],'FontSize', 9,'String', strcat('Sys:', Dat.fileName), 'Callback', @openData);

                loadparameters(1);
                if isfield(DS, 'fileNameQ') && isfield(DS, 'fileNameX') 
                    simulationQfp();
                    simulationXfp();
                    simulationQ();
                    simulationX();
                    plotwindow();
                    Dat.list = {'last'};
                    rmfield(Dat, 'listdata');
                    savetmpData(0,0,'last');
                    Dat.loadlist.Value=1;
                end
            else
                openData();  
            end
        else
            Flagfinish = 1;
        end

    end

    function systemchangeregime(~,~)
        windowF = figure('Color', [0.9255 0.9137 0.8471], 'Name', 'Parameters finder',...
                    'DockControl', 'off', 'Units', 'Pixels', 'Position', [100 50 800 600]);
        updatebutton = uicontrol('Parent', windowF, 'Style', 'pushbutton','Units', 'normalized', ...
            'Position', [0.35 0.1 0.3 0.1], 'String', 'Update without saving','Callback', @updateSys); 
        Datafilecreation(Dat.fileName,Dat.path, windowF, data);
    end

    function updateSys(~,~)
        savetmpData(0,0,'last');
        data = evalin('base', 'data');
    
        simulationQfp();
        simulationXfp();
        simulationQ();
        simulationX();
        plotwindow();
        set(Parameters_label, 'String',...
            sprintf('Sys1:g=%1.4f %1.4f %1.4f gStrain=%1.4f %1.4f %1.4f\n lwX=%1.4f lwQ=%1.4f \nSys2:g=%1.4f %1.4f %1.4f gStrain= %1.4f %1.4f %1.4f \nlwX=%1.4f  lwQ=%1.4f \nmwFreqX=%3.4f  mwFreqQ=%3.4f\n', ...
            data.X.Sys1.g, data.X.Sys1.gStrain, data.X.Sys1.lw, data.Q.Sys1.lw, data.X.Sys2.g, data.X.Sys2.gStrain,data.X.Sys2.lw,data.Q.Sys2.lw, X.Exp.mwFreq, Q.Exp.mwFreq));    
    end

    function loadparameters(fl)
         if ~isempty(Data) && isfield(DS, 'fileNameQ') && isfield(DS, 'fileNameX')  && fl
            Dat.I = length(Data);
            j=1;
            for k = 1:Dat.I 
               if strcmp(Data(k).fileNameX, DS.fileNameX) && strcmp(Data(k).fileNameQ, DS.fileNameQ)
                   str{j} = Data(k).Name;
                   str2(j) = k;j=j+1;
               end
            end
            if exist('str', 'var')
                [s] = listdlg('PromptString','Select a file:','SelectionMode','single','ListString',str,'ListSize', [300 160]);
                if ~isempty(s)
                    Dat.I = str2(s);
                    DS = Data(Dat.I);
                end
            end
         end
         if isfield(DS, 'fileNameQ') && isfield(DS, 'fileNameX') && isfield(Dat, 'fileName')
            Xband.B1.edit.String=num2str(DS.B1X);
            Xband.B2.edit.String=num2str(DS.B2X);
            Qband.B1.edit.String=num2str(DS.B1Q);
            Qband.B2.edit.String=num2str(DS.B2Q);
            Xband.time.edit.String=DS.tX;
            Qband.time.edit.String=DS.tQ;
            Xband.time.slider.Min=min(X.t); Xband.time.slider.Max=max(X.t); Xband.time.slider.Value=DS.tX;
            Qband.time.slider.Min=min(Q.t); Qband.time.slider.Max=max(Q.t); Qband.time.slider.Value=DS.tQ;
            Xband.filebutton.String=strcat('X:', DS.fileNameX);
            Qband.filebutton.String=strcat('Q:', DS.fileNameQ);

            JFP.slider.Value=DS.JFP;
            DFP.slider.Value=DS.DFP;
            if ~isfield(DS, 'EFP')
                DS.EFP = 0;
                DS.E = 0;
                EFP.edit.String=DS.EFP;
                E.edit.String=DS.E;
            end
            J.slider.Value=DS.J;
            D.slider.Value=DS.D;
            for i=1:6
                Dat.aslider(i).Value=DS.a(i).a;
                Dat.aedit(i).String=DS.a(i).a;
            end
            Q.euler1 = [DS.a(1).a DS.a(2).a DS.a(3).a];
            Q.euler2 = [DS.a(4).a DS.a(5).a DS.a(6).a];
            Q.normQ.slider.Value=DS.normQ;
            Q.normQFP.slider.Value=DS.normQFP;
            Q.normX.slider.Value=DS.normX;
            Q.normXFP.slider.Value=DS.normXFP;
            
            JFP.edit.String=DS.JFP;
            DFP.edit.String=DS.DFP;
            Dat.JFPwidth.edit.String=DS.JFPw;
            Dat.DFPwidth.edit.String=DS.DFPw;
            Dat.JFPnumber.edit.String=DS.JFPn;
            Dat.DFPnumber.edit.String=DS.DFPn;
            if DS.D>0
                Dat.DFPrange.String = sprintf('%2.3f:%2.3f', DS.DFP/(2-(DS.DFP/(DS.DFP+DS.DFPw/2))^(1/3))^3, DS.DFP+DS.DFPw/2);
            else
                Dat.DFPrange.String = sprintf('%2.3f:%2.3f', DS.DFP-DS.DFPw/2, DS.DFP/(2-(DS.DFP/(DS.DFP-DS.DFPw/2))^(1/3))^3);
            end
            if DS.J>0
                Dat.JFPrange.String = sprintf('%2.3f:%2.3f', DS.JFP^2/(DS.JFP+DS.JFPw/2), DS.JFP+DS.JFPw/2);
            else
                Dat.JFPrange.String = sprintf('%2.3f:%2.3f', DS.JFP-DS.JFPw/2, DS.JFP^2/(DS.JFP-DS.JFPw/2));
            end
            
            J.edit.String=DS.J;
            D.edit.String=DS.D;

            Xband.normX.edit.String=DS.normX;
            Xband.normXFP.edit.String=DS.normXFP;
            Qband.normQ.edit.String=DS.normQ;
            Qband.normQFP.edit.String=DS.normQFP;
            Dat.addbrX.edit.String=DS.addbrX;
            Dat.addbrQ.edit.String=DS.addbrQ;
            Dat.npFP.edit.String=DS.npFP;
            Dat.np.edit.String=DS.np;
            Xband.fitlevel.edit.String=DS.fitlevelX;
            Xband.smooth.edit.String=DS.smoothX;
            Qband.fitlevel.edit.String=DS.fitlevelQ;
            Qband.smooth.edit.String=DS.smoothQ;
            set(Parameters_label, 'String',...
                sprintf('Sys1:g=%1.4f %1.4f %1.4f gStrain=%1.4f %1.4f %1.4f\n lwX=%1.4f lwQ=%1.4f \nSys2:g=%1.4f %1.4f %1.4f gStrain= %1.4f %1.4f %1.4f \nlwX=%1.4f  lwQ=%1.4f \nmwFreqX=%3.4f  mwFreqQ=%3.4f\n', ...
                data.X.Sys1.g, data.X.Sys1.gStrain, data.X.Sys1.lw, data.Q.Sys1.lw, data.X.Sys2.g, data.X.Sys2.gStrain,data.X.Sys2.lw,data.Q.Sys2.lw, X.Exp.mwFreq, Q.Exp.mwFreq));    
         end
end

%% save files
    function saveData(~, ~)
            DS.B1X = str2double(Xband.B1.edit.String); 
            DS.B2X = str2double(Xband.B2.edit.String); 
            DS.B1Q = str2double(Qband.B1.edit.String); 
            DS.B2Q = str2double(Qband.B2.edit.String); 
            choice = bttnChoiseDialog({'System data file','Picture','ASCII','mat-file'}, 'In which format to save? Save to', 'mat-file', 'Save');
            switch choice
            case 1
                prompt = {'Name in System data file and file:'};
                dlg_title = 'Input';
                num_lines = 1;
                defaultans = {strcat(datestr(now,'ddmmyy'),',J=',num2str(round(DS.J)),',D=[',num2str(round(DS.D)),num2str(round(DS.E)),...
                    '],JFP=',num2str(round(DS.JFP)),',DFP=[',num2str(round(DS.DFP)),num2str(round(DS.EFP)), ']'), ''};
                answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
                        
                rmfield(DS,'R2X');
                rmfield(DS,'R2Q');
                DS.R2fo = str2double(R2fo_edit.String);
                DS.R2fp = str2double(R2fp_edit.String);
                DS.Name = answer{1};
                if length(Data) == 0
                    Data = DS;  
                else
                    Data(length(Data)+1) = DS;  
                end
                simulations = Data;
                save(strcat(Dat.path,Dat.fileName),'simulations','-v6','-append');
                h = msgbox('Saved');
                
            case 4
                [fileName,path] = uiputfile('*.mat','Name of the file, where data will be stored');                
                fid = fopen(strcat(path, fileName),'wt');
                if fullpowder_checkbox.Value
                    Q_band.FullPowder.B0_simulation = xdataQfp;
                    Q_band.FullPowder.spectrum_simulation = ydataQfp;
                    Q_band.FullPowder.spectra_sum_scaled_simulation = ydataQtotfpnorm;
                    X_band.FullPowder.B0_simulation = xdataXfp;
                    X_band.FullPowder.spectrum_simulation = ydataXfp;
                    X_band.FullPowder.spectra_sum_scaled_simulation = ydataXtotfpnorm;
                end
                if fixorient_checkbox.Value
                    Q_band.FixedRElativeOrient.B0_simulation = xdataQ;
                    Q_band.FixedRElativeOrient.spectrum_simulation = ydataQ;
                    Q_band.FixedRElativeOrient.spectra_sum_scaled_simulation = ydataQtotnorm;
                    X_band.FixedRElativeOrient.B0_simulation = xdataX;
                    X_band.FixedRElativeOrient.spectrum_simulation = ydataX;
                    X_band.FixedRElativeOrient.spectra_sum_scaled_simulation = ydataXtotnorm;
                end
                Q_band.B0 = Q.B0;
                Q_band.spectrum_initial = Q.s(:,DS.sliceQ);
                Q_band.spectrum_bgcorrected = Q.sFit(:,DS.sliceQ);
                Q_band.spectrum_bgcorrected_smoothed = Q.smooth(:,DS.sliceQ);
                X_band.B0 = X.B0;
                X_band.spectrum_initial = X.s(:,DS.sliceX);
                X_band.spectrum_bgcorrected = X.sFit(:,DS.sliceX);
                X_band.spectrum_bgcorrected_smoothed = X.smooth(:,DS.sliceX);
                
                info=string();
                info=strcat(info, sprintf('Files used for simulation:\nXband: %s\nQband: %s\nSystem data: %s', DS.fileNameX, DS.fileNameQ, Dat.fileName));
                info=strcat(info, sprintf('\n%30s %15s %15s\n\n', 'Parametrs', 'Xband', 'Qband')); 
                info=strcat(info, sprintf('\n%30s %15.5f %15.5f\n', 'Frequency(GHz):', DS.mwFreqX, DS.mwFreqQ)); 
                info=strcat(info, sprintf('\n%30s %15.5f %15.5f\n', 'Time slice(mus):', DS.tX, DS.tQ));  
                info=strcat(info, sprintf('\n%30s %15u %15u\n', 'Slice:', DS.sliceX, DS.sliceQ));  
                
                info=strcat(info, sprintf('\n\n                       Parameters used for background corrections: \n')); 
                info=strcat(info, sprintf('\n%30s %15u %15u\n', 'Smooth power:', DS.smoothX, DS.smoothQ));  
                info=strcat(info, sprintf('\n%30s %15u %15u\n', 'Polynom power used:', DS.fitlevelX, DS.fitlevelQ));  
                info=strcat(info, sprintf('\n%30s %15.5f %15.5f\n', 'Lower bound(mT):', DS.B1X, DS.B1Q)); 
                info=strcat(info, sprintf('\n%30s %15.5f %15.5f\n', 'Upper bound(mT):', DS.B2X, DS.B2Q)); 
                
                info=strcat(info, sprintf('\n\n                          Parameters used for simulations: \n')); 
                info=strcat(info, sprintf('\n%30s %15.5f %15.5f\n', 'Additional broadening(mT)', DS.addbrX, DS.addbrQ)); 
                
                if fullpowder_checkbox.Value
                    info=strcat(info, sprintf('\n\n                                       Full powder : \n')); 
                    info=strcat(info, sprintf('\n%30s %20u\n', 'Number of points used', DS.npFP)); 
                    info=strcat(info, sprintf('\n%30s %23.5f\n', 'Exchange coupling(MHz)', DS.JFP)); 
                    info=strcat(info, sprintf('\n%30s %23.5f\n', 'Exchange coupling width(MHz)', DS.JFPw)); 
                    info=strcat(info, sprintf('\n%30s %20s\n', 'Range of exchange coupling(MHz)',  Dat.JFPrange.String)); 
                    info=strcat(info, sprintf('\n%30s %20u\n', 'Exchange coupling number of poin', DS.JFPn)); 
                    info=strcat(info, sprintf('\n%30s %23.5f\n', 'Dipolar coupling D(MHz)', DS.DFP)); 
                    info=strcat(info, sprintf('\n%30s %23.5f\n', 'Dipolar coupling E(MHz)', DS.EFP));
                    info=strcat(info, sprintf('\n%30s %23.5f\n', 'Dipolar coupling width(MHz)', DS.DFPw)); 
                    info=strcat(info, sprintf('\n%30s %20s\n', 'Range of dipolar coupling(MHz)',  Dat.DFPrange.String)); 
                    info=strcat(info, sprintf('\n%30s %20u\n', 'Dipolar coupling number of poin', DS.DFPn)); 
                    info=strcat(info, sprintf('\n%30s %20u\n', 'R2', DS.R2fp)); 
                    info=strcat(info, sprintf('\n%30s %15.5f %15.5f\n', 'Scaling', DS.normXFP, DS.normQFP)); 
                end

                if fixorient_checkbox.Value
                    info=strcat(info, sprintf('\n\n                                    Fixed orientation: \n')); 
                    info=strcat(info, sprintf('\n%30s %20u\n', 'Number of points used', DS.np)); 
                    info=strcat(info, sprintf('\n%30s %23.5f\n', 'Exchange coupling(MHz)', DS.J)); 
                    info=strcat(info, sprintf('\n%30s %23.5f\n', 'Dipolar coupling D(MHz)', DS.D)); 
                    info=strcat(info, sprintf('\n%30s %23.5f\n', 'Dipolar coupling E(MHz)', DS.E)); 
                    info=strcat(info, sprintf('\n%30s %20u\n', 'R2', DS.R2fo)); 
                    info=strcat(info, sprintf('\n%30s %15.5f %15.5f\n', 'Scaling', DS.normX, DS.normQ)); 
                    info=strcat(info, sprintf('\nEuler angles(radians): System1(%f %f %f) System2(%f %f %f)\n', DS.a(1).a,DS.a(2).a,DS.a(3).a,DS.a(4).a,DS.a(5).a,DS.a(6).a)); 
                end
                
                info=strcat(info, sprintf('\n\n%15s %25s %25s\n', 'Parametrs', 'System1', 'System2'));  
                info=strcat(info, sprintf('\n%15s %25.1f %25f\n', 'Spin', data.X.Sys1.S, data.X.Sys2.S)); 
                info=strcat(info, sprintf('\n%15s %25.5f %25.5f\n', 'Weight', data.X.Sys1.weight, data.X.Sys2.weight)); 
                info=strcat(info, sprintf('\n%15s %25.5f %25.5f\n', 'Linewidth Xband', data.X.Sys1.lw, data.X.Sys2.lw)); 
                info=strcat(info, sprintf('\n%15s %25.5f %25.5f\n', 'Linewidth Qband', data.Q.Sys1.lw, data.Q.Sys2.lw)); 
                info=strcat(info, sprintf('\n%15s (%f %f %f) (%f %f %f)\n', 'g', data.X.Sys1.g, data.X.Sys2.g)); 
                info=strcat(info, sprintf('\n%15s (%f %f %f) (%f %f %f)\n', 'gStrain', data.X.Sys1.gStrain, data.X.Sys2.gStrain)); 
                
                
                save(strcat(path,fileName),'info', 'X_band', 'Q_band');
                
            case 3
                [fileName,path] = uiputfile('*.txt','Name of the file, where data will be stored');
                fid = fopen(strcat(path, fileName),'wt');
                fprintf(fid, 'Files used for simulation:\nXband: %s\nQband: %s\nSystem data: %s\n\n\n', DS.fileNameX, DS.fileNameQ, Dat.fileName);
                fprintf(fid, '%30s %15s %15s\n\n', 'Parametrs', 'Xband', 'Qband'); 
                fprintf(fid, '%30s %15.5f %15.5f\n', 'Frequency(GHz):', DS.mwFreqX, DS.mwFreqQ);
                fprintf(fid, '%30s %15.5f %15.5f\n', 'Time slice(mus):', DS.tX, DS.tQ); 
                fprintf(fid, '%30s %15u %15u\n', 'Slice:', DS.sliceX, DS.sliceQ); 
                
                fprintf(fid, '                       Parameters used for background corrections: \n');
                fprintf(fid, '%30s %15u %15u\n', 'Smooth power:', DS.smoothX, DS.smoothQ); 
                fprintf(fid, '%30s %15u %15u\n', 'Polynom power used:', DS.fitlevelX, DS.fitlevelQ); 
                fprintf(fid, '%30s %15.5f %15.5f\n', 'Lower bound(mT):', DS.B1X, DS.B1Q);
                fprintf(fid, '%30s %15.5f %15.5f\n', 'Upper bound(mT):', DS.B2X, DS.B2Q);
                
                fprintf(fid, '                          Parameters used for simulations: \n');
                fprintf(fid, '%30s %15.5f %15.5f\n', 'Additional broadening(mT)', DS.addbrX, DS.addbrQ);
                
                fprintf(fid, '                                       Full powder : \n');
                fprintf(fid, '%30s %20u\n', 'Number of points used', DS.npFP);
                fprintf(fid, '%30s %23.5f\n', 'Exchange coupling(MHz)', DS.JFP);
                fprintf(fid, '%30s %23.5f\n', 'Exchange coupling width(MHz)', DS.JFPw);
                fprintf(fid, '%30s %20s\n', 'Range of exchange coupling(MHz)',  Dat.JFPrange.String);
                fprintf(fid, '%30s %20u\n', 'Exchange coupling number of poin', DS.JFPn);
                fprintf(fid, '%30s %23.5f\n', 'Dipolar coupling D(MHz)', DS.DFP);
                fprintf(fid, '%30s %23.5f\n', 'Dipolar coupling E(MHz)', DS.EFP);
                fprintf(fid, '%30s %23.5f\n', 'Dipolar coupling width(MHz)', DS.DFPw);
                fprintf(fid, '%30s %20s\n', 'Range of dipolar coupling(MHz)',  Dat.DFPrange.String);
                fprintf(fid, '%30s %20u\n', 'Dipolar coupling number of poin', DS.DFPn);
                fprintf(fid, '%30s %20u\n', 'R2', DS.R2fp);
                fprintf(fid, '%30s %15.5f %15.5f\n', 'Scaling', DS.normXFP, DS.normQFP);
                
                fprintf(fid, '                                    Fixed orientation: \n');
                fprintf(fid, '%30s %20u\n', 'Number of points used', DS.np);
                fprintf(fid, '%30s %23.5f\n', 'Exchange coupling(MHz)', DS.J);
                fprintf(fid, '%30s %23.5f\n', 'Dipolar coupling D(MHz)', DS.D);
                fprintf(fid, '%30s %23.5f\n', 'Dipolar coupling E(MHz)', DS.E);
                fprintf(fid, '%30s %20u\n', 'R2', DS.R2fo);
                fprintf(fid, '%30s %15.5f %15.5f\n', 'Scaling', DS.normX, DS.normQ);
                fprintf(fid, 'Euler angles(radians): System1(%f %f %f) System2(%f %f %f)\n', DS.a(1).a,DS.a(2).a,DS.a(3).a,DS.a(4).a,DS.a(5).a,DS.a(6).a);
                
                fprintf(fid, '\n\n%15s %25s %25s\n', 'Parametrs', 'System1', 'System2'); 
                fprintf(fid, '%15s %25.1f %25f\n', 'Spin', data.X.Sys1.S, data.X.Sys2.S);
                fprintf(fid, '%15s %25.5f %25.5f\n', 'Weight', data.X.Sys1.weight, data.X.Sys2.weight);
                fprintf(fid, '%15s %25.5f %25.5f\n', 'Linewidth Xband', data.X.Sys1.lw, data.X.Sys2.lw);
                fprintf(fid, '%15s %25.5f %25.5f\n', 'Linewidth Qband', data.Q.Sys1.lw, data.Q.Sys2.lw);
                fprintf(fid, '%15s (%f %f %f) (%f %f %f)\n', 'g', data.X.Sys1.g, data.X.Sys2.g);
                fprintf(fid, '%15s (%f %f %f) (%f %f %f)\n', 'gStrain', data.X.Sys1.gStrain, data.X.Sys2.gStrain);
                fprintf(fid, '\n %15s  %15s  %15s  %15s  %15s  %15s  %15s  %15s  %15s  %15s  %15s  %15s   %15s  %15s  %15s  %15s   %15s  %15s  %15s  %15s   %15s  %15s  %15s  %15s   %15s  %15s  %15s  %15s\n',...
                    'XbandB0(mT)', 'QbandB0(mT)', 'initialX','initialQ','bgcorrX','bgcorrQ','bgcorrsmoothX','bgcorrsmoothQ',...
                    'FullPowderSimX','FullPowderSimX1','FullPowderSimX2','FullPowderSimX3','FullPowderSimX4','FullPowderSimQ','FullPowderSimQ1','FullPowderSimQ2','FullPowderSimQ3','FullPowderSimQ4',...
                    'FixedOrSimX', 'FixedOrSimX1','FixedOrSimX2','FixedOrSimX3','FixedOrSimX4','FixedOrSimQ','FixedOrSimQ1','FixedOrSimQ2','FixedOrSimQ3','FixedOrSimQ4');
                ad0X = zeros(1,length(Q.B0)-length(X.B0));
                ad0Q = zeros(1,length(X.B0)-length(Q.B0));
                fprintf(fid, ' %15.5f  %15.5f  %15.5f  %15.5f  %15.5f  %15.5f  %15.5f  %15.5f  %15.5f  %15.5f  %15.5f  %15.5f  %15.5f  %15.5f  %15.5f  %15.5f  %15.5f  %15.5f  %15.5f  %15.5f  %15.5f  %15.5f  %15.5f  %15.5f  %15.5f  %15.5f  %15.5f  %15.5f\n',...
                    [X.B0,ad0X;Q.B0,ad0Q;X.s(:,DS.sliceX)',ad0X;Q.s(:,DS.sliceQ)',ad0Q;X.sFit(:,DS.sliceX)',ad0X;Q.sFit(:,DS.sliceQ)',...
                    ad0Q;X.smooth(:,DS.sliceX)',ad0X;Q.smooth(:,DS.sliceQ)',ad0Q;...
                    ydataXtotfp,ad0X;ydataXfp{1},ad0X;ydataXfp{2},ad0X;ydataXfp{3},ad0X;ydataXfp{4},ad0X;...
                    ydataQtotfp,ad0Q;ydataQfp{1},ad0Q;ydataQfp{2},ad0Q;ydataQfp{3},ad0Q;ydataQfp{4},ad0Q;...
                    ydataXtot,ad0X;ydataX{1},ad0X;ydataX{2},ad0X;ydataX{3},ad0X;ydataX{4},ad0X;...
                    ydataQtot,ad0Q;ydataQ{1},ad0Q;ydataQ{2},ad0Q;ydataQ{3},ad0Q;ydataQ{4},ad0Q]);
                fclose(fid);
                h = msgbox('Saved in ASCII file'); 

                case 2
                 choice = questdlg('In which format to save? Save to', 'Format', '.bmp', '.eps','.eps');
                 switch choice
                 case '.bmp'  
                    [fileName1,path1,FilterIndex] = uiputfile('*.bmp','Name of the file, where data will be stored');
                    format = '-dbmp';
                 case '.eps'
                     [fileName1,path1,FilterIndex] = uiputfile('*.eps','Name of the file, where data will be stored');
                     format = '-depsc';
                 end
                 
                 if ~isempty(choice)
                 if FilterIndex~=0
                    h2 = figure; 
                    set(h2,'Units','normalized','Position',[1 1 1 1]);                    
                    Xband.plot.Position = [0.0750    0.3500    0.4150    0.5800];
                    Qband.plot.Position = [0.5600    0.3500    0.4150    0.5800];
                    X1=copyobj([Xband.plot, legend(Xband.plot)],h2);
                    Q1=copyobj([Qband.plot, legend(Qband.plot)],h2);
                    str=string();
                    if fullpowder_checkbox.Value == 1 || fixorient_checkbox.Value == 1
                            str2 = sprintf('%sTime for Xband %.3fmus; for Qband %.3fmus', str, DS.tX, DS.tQ);
                            str = sprintf('%sAdditional broadening for Xband %2.3fmT; for Qband %2.3fmT', str, DS.addbrX, DS.addbrQ);
                            if fullpowder_checkbox.Value == 1 
                                if DS.JFPw==0
                                    str = sprintf('%s\n\nFULL POWDER\nJ = %2.3fMHz', str, DS.JFP);
                                else
                                    str = sprintf('%s\n\nFULL POWDER\nJ = %sMHz with mean value %2.3fMHz and total number of points = %u', str, Dat.JFPrange.String, DS.JFP, DS.JFPn);
                                end
                                if DS.DFPw==0
                                    str = sprintf('%s; D = [%2.3f %2.3f]MHz', str,DS.DFP,DS.EFP);
                                else
                                    str = sprintf('%s; \nD = %sMHz with mean value %2.3fMHz and total number of points = %u', str, Dat.DFPrange.String, DS.DFP, DS.DFPn);
                                end
                                str = sprintf('%s\nNumber of points used for simulation %u', str, DS.npFP);
                            end
                            if fixorient_checkbox.Value == 1 
                                str = sprintf('%s\n\nFIXED ORIENTATION POWDER\nJ = %2.3fMHz', str,DS.J);
                                str = sprintf('%s; D = [%2.3f %2.3f]MHz', str,DS.D,DS.E);
                                str = sprintf('%s\nNumber of points used for simulation %u', str, DS.np);
                                str = sprintf('%s\nEuler angles in radians. System1:%1.3f %1.3f %1.3f; System2:%1.3f %1.3f %1.3f', str, DS.a(1).a,DS.a(2).a,DS.a(3).a,DS.a(4).a,DS.a(5).a,DS.a(6).a);
                            end
                            
                            if savedSpectra_checkbox.Value == 1
                                ind = Dat.loadlist.Value;
                                str = sprintf('%s\n\nSaved spectra(dashed lines) "%s"', str, Dat.loadlist.String{ind});
                            end
                    end
                    str2 = sprintf('%s\n\nSYSTEM PARAMETRS\nSystem1:\ng=%1.4f %1.4f %1.4f gStrain=%1.4f %1.4f %1.4f\nlwX=%1.4f lwQ=%1.4f \nSystem2:\ng=%1.4f %1.4f %1.4f gStrain=%1.4f %1.4f %1.4f \nlwX=%1.4f  lwQ=%1.4f \nmwFreqX=%3.4fGHz mwFreqQ=%3.4fGHz\n', ...
                          str2,data.X.Sys1.g, data.X.Sys1.gStrain, data.X.Sys1.lw, data.Q.Sys1.lw, data.X.Sys2.g, data.X.Sys2.gStrain,data.X.Sys2.lw,data.Q.Sys2.lw, X.Exp.mwFreq, Q.Exp.mwFreq);    
                    uicontrol('Parent', h2, 'Style', 'text', 'fontunits', 'normalized', 'fontsize', 0.055, 'Units', 'normalized',...
                                'String', str,'Position', [0.03 0.01 0.46 0.258], 'Backgroundcolor',[1, 1, 1],'HorizontalAlignment','left');
                    uicontrol('Parent', h2, 'Style', 'text', 'fontunits', 'normalized', 'fontsize', 0.055, 'Units', 'normalized',...
                                'String', str2,'Position', [0.51 0.01 0.46 0.258], 'Backgroundcolor',[1, 1, 1],'HorizontalAlignment','left');
                    print(h2,format,[strcat(path1, fileName1)])
                    Xband.plot.Position = [0.0300    0.3900    0.4150    0.5800];
                    Qband.plot.Position = [0.5300    0.3900    0.4150    0.5800];
                    clf 
                    close(h2)
                    h = msgbox('Saved as picture'); 
                 end
                 end
            end
    end

%% D and J fitting
    function findfpDJ(~, ~)
           h = msgbox('Wait'); 
           params = [DS.JFP, DS.DFP];  

           ydata = [Q.smooth(:,DS.sliceQ); X.smooth(:,DS.sliceX)];         
           
           figure('Name', 'Before fitting');
           str = {['J=', num2str(DS.JFP), ' D=', num2str(DS.DFP), ' NormX=', num2str(DS.normXFP), ' NormQ=',  num2str(DS.normQFP)]};
           uicontrol( 'Units','normalized','Position', [0.0 0.95 0.3 0.05], 'style','text','String', str, 'FontSize', 12, 'FontWeight', 'bold');
           subplot (2, 1, 1);
           plot(X.B0, X.smooth(:,DS.sliceX), 'b', xdataXfp, ydataXtotfpnorm, 'k') 
           legend('background corrected smoothed', 'simulation','Location','Best'); legend boxoff
           title('Xband');
           xlabel('B(mT)');
           ylabel('I');
           subplot (2, 1, 2);
           plot(Q.B0, Q.smooth(:,DS.sliceQ), 'b', xdataQfp, ydataQtotfpnorm, 'k') 
           legend('background corrected smoothed', 'simulation','Location','Best'); legend boxoff
           title('Qband');
           xlabel('B(mT)');
           ylabel('I');
            
           opts = optimoptions(@lsqcurvefit,'Algorithm','levenberg-marquardt','Display','iter-detailed');
           lb = [0, 0];
           ub = [maxJ, maxD];
           [parameter] = lsqcurvefit(@func,params,[],ydata,lb,ub,opts);
           JFP.edit.String = parameter(1);
           DFP.edit.String = parameter(2);
%            normXFP_edit.String = parameter(4);
%            normQFP_edit.String = parameter(3);
%            updateGraphFP2Norm();
           updateGraphFP2();
           delete(h);
           h = msgbox('Done'); 
    end

    function ysim = func(x,~)
       [xsim1,ysimtemp1] = ctssimaniso_git_initial_v2_1(data.Q.Sys1,data.Q.Sys2,Q.Exp,Dat.gOrientation,Dat.DOrientation,x(1),x(2),DS.addbrQ,length(Q.B0)); 
       ysimtot1 = ysimtemp1{1} +  ysimtemp1{2} + ysimtemp1{3} + ysimtemp1{4};

       [xsim2,ysimtemp2] = ctssimaniso_git_initial_v2_1(data.X.Sys1,data.X.Sys2,X.Exp,Dat.gOrientation,Dat.DOrientation,x(1),x(2),DS.addbrX,length(X.B0)); 
       ysimtot2 = ysimtemp2{1} +  ysimtemp2{2} + ysimtemp2{3} + ysimtemp2{4};
       
       ys = [ysimtot1'/ abs(max(ysimtot1)*DS.normQFP); ysimtot2' / abs(max(ysimtot2) * DS.normXFP)];
%        figure(5)
%        plot(xsim1,Q.smooth(:,DS.sliceQ), xsim1, ysimtot1'/ abs(max(ysimtot1)*DS.normQFP))
%        sum((ys-[Q.smooth(:,DS.sliceQ); X.smooth(:,DS.sliceX)]).^2)
       ysim = ys;
    end

%% different regimes
    function advancedregime(~, ~)
        set(Dat.fitlevellabel, 'Visible', 'off');
        set(Dat.smoothlabel, 'Visible', 'off');
        set(Dat.bgwithoutsmoothlabel, 'Visible', 'off');
        set(Dat.bgwithoutsmooth_checkbox, 'Visible', 'off');
        set(Xband.fitlevel.edit, 'Visible', 'off');
        set(Xband.smooth.edit, 'Visible', 'off');
        set(Qband.fitlevel.edit, 'Visible', 'off');
        set(Qband.smooth.edit, 'Visible', 'off');
        set(Dat.BGrect, 'Visible','off');
        set(Dat.Xbandrect, 'Visible','on');
        Dat.Xbandrect.Position = [0.51 0.138 0.125 0.072];
        set(Dat.Qbandrect, 'Visible','on');
        Dat.Qbandrect.Position = [0.635 0.138 0.125 0.072];
        set(Dat.Xbandlabel, 'Visible','on');
        set(Dat.Qbandlabel, 'Visible','on');
        set(Dat.addbrlabel, 'Visible','on');
        set(Dat.numberofpointslabel, 'Visible','on');
        set(Dat.npfulllabel, 'Visible','on');
        set(Dat.npfixlabel, 'Visible','on');
        set(Dat.addbrX.edit, 'Visible','on');
        set(Dat.addbrQ.edit, 'Visible','on');
        set(Dat.npFP.edit, 'Visible','on');
        set(Dat.np.edit, 'Visible','on');
        Dat.regime = 'adv';
    end

    function bgcorrregime(~, ~)
        set(Dat.BGrect, 'Visible','off');
        set(Dat.Xbandrect, 'Visible','on');
        Dat.Xbandrect.Position = [0.51 0.065 0.125 0.145];
        set(Dat.Qbandrect, 'Visible','on');
        Dat.Qbandrect.Position = [0.635 0.065 0.125 0.145];
        set(Dat.Xbandlabel, 'Visible','on');
        set(Dat.Qbandlabel, 'Visible','on');
        set(Dat.addbrlabel, 'Visible','off');
        set(Dat.numberofpointslabel, 'Visible','off');
        set(Dat.npfulllabel, 'Visible','off');
        set(Dat.npfixlabel, 'Visible','off');
        set(Dat.addbrX.edit, 'Visible','off');
        set(Dat.addbrQ.edit, 'Visible','off');
        set(Dat.npFP.edit, 'Visible','off');
        set(Dat.np.edit, 'Visible','off');  
        set(Dat.fitlevellabel, 'Visible', 'on');
        set(Dat.smoothlabel, 'Visible', 'on');
        set(Dat.bgwithoutsmoothlabel, 'Visible', 'on');
        set(Dat.bgwithoutsmooth_checkbox, 'Visible', 'on');
        set(Xband.fitlevel.edit, 'Visible', 'on');
        set(Xband.smooth.edit, 'Visible', 'on');
        set(Qband.fitlevel.edit, 'Visible', 'on');
        set(Qband.smooth.edit, 'Visible', 'on');
        Dat.regime = 'bg';
    end

    function regularregime(~, ~)
        set(Dat.BGrect, 'Visible','off');
        set(Dat.Xbandrect, 'Visible','off');
        Dat.Xbandrect.Position = [0.51 0.065 0.125 0.145];
        set(Dat.Qbandrect, 'Visible','off');
        Dat.Qbandrect.Position = [0.635 0.065 0.125 0.145];
        set(Dat.Xbandlabel, 'Visible','off');
        set(Dat.Qbandlabel, 'Visible','off');
        set(Dat.addbrlabel, 'Visible','off');
        set(Dat.numberofpointslabel, 'Visible','off');
        set(Dat.npfulllabel, 'Visible','off');
        set(Dat.npfixlabel, 'Visible','off');
        set(Dat.addbrX.edit, 'Visible','off');
        set(Dat.addbrQ.edit, 'Visible','off');
        set(Dat.npFP.edit, 'Visible','off');
        set(Dat.np.edit, 'Visible','off');  
        set(Dat.fitlevellabel, 'Visible', 'off');
        set(Dat.smoothlabel, 'Visible', 'off');
        set(Dat.bgwithoutsmoothlabel, 'Visible', 'off');
        set(Dat.bgwithoutsmooth_checkbox, 'Visible', 'off');
        set(Xband.fitlevel.edit, 'Visible', 'off');
        set(Xband.smooth.edit, 'Visible', 'off');
        set(Qband.fitlevel.edit, 'Visible', 'off');
        set(Qband.smooth.edit, 'Visible', 'off');
        Dat.regime = 'hide';
    end

%% temporal saving
    function savetmpData(~, ~, str)
        if exist('str', 'var')
            answer{1} = 'last';
            ind = 1;
        else
            prompt = {'Name:'};
            dlg_title = 'Save temporaly';
            num_lines = 1;
            defaultans = {strcat('J=',num2str(round(DS.J)),',D=',num2str(round(DS.D)),...
                ',JFP=',num2str(round(DS.JFP)),',DFP=',num2str(round(DS.DFP)), ''), ''};
            defaultans = strrep(defaultans, '-', '_');
            answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
            ind = length(Dat.list) +1;
        end

        DStmpX.B1 = str2double(Xband.B1.edit.String); 
        DStmpX.B2 = str2double(Xband.B2.edit.String); 
        DStmpX.fp = ydataXtotfp / max(ydataXtotfp);
        DStmpX.fptemp = ydataXfp;
        DStmpX.fix = ydataXtot / max(ydataXtot);
        DStmpX.fixtemp = ydataX;
        
        DStmpQ.B1 = str2double(Qband.B1.edit.String); 
        DStmpQ.B2 = str2double(Qband.B2.edit.String); 
        DStmpQ.fp = ydataQtotfp / max(ydataQtotfp);
        DStmpQ.fptemp = ydataQfp;
        DStmpQ.fix = ydataQtot / max(ydataQtot);
        DStmpQ.fixtemp = ydataQ;
        
        DS.R2X = str2double(R2X_edit.String);
        DS.R2Q = str2double(R2Q_edit.String);
        DS.R2fp = str2double(R2fp_edit.String);
        DS.R2fo = str2double(R2fo_edit.String);
        
        Dat.list(ind) = answer(1); 
        Dat.listdata{ind}.DS = DS; 
        Dat.listdata{ind}.DStmpX = DStmpX; 
        Dat.listdata{ind}.DStmpQ = DStmpQ; 
        Dat.listdata{ind}.data = data; 
        Dat.loadlist.String = Dat.list;     
    end

    function loadtmpData(~, ~)
        ind = Dat.loadlist.Value; 
        DS = Dat.listdata{ind}.DS; 
        DStmpX = Dat.listdata{ind}.DStmpX; 
        DStmpQ = Dat.listdata{ind}.DStmpQ; 
        data = Dat.listdata{ind}.data;
         
        Xband.B1.edit.String = DStmpX.B1; 
        Xband.B2.edit.String = DStmpX.B2; 
        
        Qband.B1.edit.String = DStmpQ.B1; 
        Qband.B2.edit.String = DStmpQ.B2; 
                           
        R2X_edit.String = DS.R2X;
        R2Q_edit.String = DS.R2Q;
        R2fp_edit.String = DS.R2fp; 
        R2fo_edit.String = DS.R2fo; 
        
        Qband.time.slider.Value = Q.t(DS.sliceQ);
        Xband.time.slider.Value = X.t(DS.sliceX);
        Qband.time.edit.String = Q.t(DS.sliceQ);
        Xband.time.edit.String = X.t(DS.sliceX);
        ydataQtotfp = DStmpQ.fp;
        ydataQtotfpnorm = ydataQtotfp*DS.normQFP;
        ydataQfp = DStmpQ.fptemp;
        ydataQtot = DStmpQ.fix;
        ydataQtotnorm = ydataQtot*DS.normQ;
        ydataQ = DStmpQ.fixtemp;
        ydataXtotfp = DStmpX.fp;
        ydataXtotfpnorm = ydataXtotfp*DS.normXFP;
        ydataXfp = DStmpX.fptemp;
        ydataXtot = DStmpX.fix;
        ydataXtotnorm = ydataXtot*DS.normX;
        ydataX = DStmpX.fixtemp;
                   
        loadparameters(0);
        
        Xbgcorr();
        Qbgcorr();
        np(0,0,0);
        npFP(0,0,0);
        plotwindow();
        
        switch  Dat.regime
            case 'hide'
                regularregime();
            case 'adv'
                advancedregime();
            case 'bg'
                bgcorrregime();
        end
                
    end

    function deletetmpData(~, ~)
         if Dat.loadlist.Value ~= 1  
            Dat.list(Dat.loadlist.Value) = []; 
            Dat.listdata(Dat.loadlist.Value) = [];
            Dat.loadlist.Value = Dat.loadlist.Value - 1;
            Dat.loadlist.String = Dat.list;  
         end
    end

%% draw orientation
    function draworient(~, ~)
        if isempty(findobj('type','figure','Name','Initial relative orientation'))
             f=figure('Name','Relative orientation');
             hh1=axes('Parent', f);
        else            
            f=findobj('type','figure','Name','Relative orientation');
            hh1=gca(f);
            vw = hh1.View;
            clf(f);  
            hh1=axes('Parent', f);
        end
        R1=eulerR(Q.euler1);
        R2=eulerR(Q.euler2);
        hold(hh1, 'on')
        plot3(hh1,[0, 0],[0, 0], [-3, 3],'b','LineStyle', '-', 'LineWidth', 2)
        plot3(hh1,[0, R1(1,1), 0, R1(1,2), 0,R1(1,3)],[0,R1(2,1), 0,R1(2,2), 0,R1(2,3)], [0,R1(3,1),0,R1(3,2),0, R1(3,3)],'r')
        plot3(hh1,[0, R2(1,1), 0,R2(1,2), 0,R2(1,3)],[0,R2(2,1), 0,R2(2,2), 0,R2(2,3)], [-1.105,-1.105+R2(3,1),-1.105,-1.105+R2(3,2),-1.105,-1.105+R2(3,3)],'k')
        plot3(hh1,[0, R2(1,1), 0,R2(1,2), 0,R2(1,3)],[0,R2(2,1), 0,R2(2,2), 0,R2(2,3)], [0.005,0.005+R2(3,1),0.005,0.005+R2(3,2),0.005,0.005+R2(3,3)],'k','LineStyle', ':')
        legend(hh1,'Magnetic field','System1', 'System2', 'Location','Best')
        axes(hh1)
        text([R1(1,1), R1(1,2),R1(1,3)],[R1(2,1),R1(2,2),R1(2,3)], [R1(3,1),R1(3,2),R1(3,3)],['gx'; 'gy'; 'gz'],'Color', 'r');
        text([R2(1,1), R2(1,2),R2(1,3)],[R2(2,1),R2(2,2),R2(2,3)], [-1.105+R2(3,1)+0.1,-1.105+R2(3,2)+0.1,-1.105+R2(3,3)+0.1],['gx'; 'gy'; 'gz']);
        text([R2(1,1), R2(1,2),R2(1,3)],[R2(2,1),R2(2,2),R2(2,3)], [R2(3,1)+0.1,R2(3,2)+0.1,R2(3,3)+0.1],['gx'; 'gy'; 'gz']);
        text (0.05, 0.05, 1.5,'B', 'Color', 'b');
        grid(hh1, 'on')
        grid(hh1, 'minor')
        xlim([-1.7;1.7]);
        ylim([-1.7;1.7]);
        zlim([-1.7;1.7]);
        if exist('vw', 'var')
            view(vw);
        end
        
        function [r] = roty(a)
            r = [cos(a) 0 -sin(a); 0 1 0; sin(a) 0 cos(a)];
        end
        function [r] = rotz(a)
            r = [cos(a) sin(a) 0; -sin(a) cos(a) 0; 0 0 1];
        end
        function [r] = eulerR(a)
            r = rotz(a(3))*roty(a(2))*rotz(a(1));
        end

    end

end


