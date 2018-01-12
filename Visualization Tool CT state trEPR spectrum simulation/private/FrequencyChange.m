function [mwFreq] = FrequencyChange(windowFreq, data1, Dat, Exp, DS, addbr, norm, normFP, B0, spec, Freq0)
%         euler1 = [DS.a(1).a DS.a(2).a DS.a(3).a];  
%         euler2 = [DS.a(4).a DS.a(5).a DS.a(6).a];
%         
%         if DS.E == 0
%             [xdata,ydata1] = ctssimaniso_git_v3_6_1(data1.Sys1,data1.Sys2,Exp,Dat.gOrientation,Dat.DOrientation,DS.J,DS.D,addbr, euler1, euler2, length(B0)); 
%         else
%             [xdata,ydata1] = ctssimaniso_git_v4(data1.Sys1,data1.Sys2,Exp,Dat.gOrientation,Dat.DOrientation,DS.J,DS.D,DS.E,addbr, euler1, euler2, length(B0));
%         end
%         ydatatot1 = ydata1{1} + ydata1{2} + ydata1{3} + ydata1{4};
%         ydatatotnorm1 = ydatatot1/max(ydatatot1) * norm;
%    
%         if DS.EFP == 0
%             [xdatafp,ydatafp1] = ctssimaniso_git_initial_v2_1(data1.Sys1,data1.Sys2,Exp,Dat.gOrientationFP,Dat.DOrientationFP,DS.JFP,DS.DFP,DS.addbr, length(B0)); 
%         else
%             [xdatafp,ydatafp1] = ctssimaniso_git_initial_v2_4(data1.Sys1,data1.Sys2,Exp,Dat.gOrientationFP,Dat.DOrientationFP,DS.JFP,DS.DFP,DS.EFP,DS.addbr, length(B0)); 
%         end
%         ydatatotfp1 = ydatafp1{1} + ydatafp1{2} + ydatafp1{3} + ydatafp1{4};
%         ydatatotfpnorm1 = ydatatotfp1/max(ydatatotfp1) * normFP;
%     
        
    Exp.mwFreq = Freq0;
    mwFreq = Exp.mwFreq;
    
    Freq.plot = axes('Parent', windowFreq,'Units', 'normalized', 'Position', [0.03 0.19 0.94 0.78]);
    
    Freq.edit = uicontrol('Parent', windowFreq,'Style', 'edit','FontSize', 15,...
                    'Units', 'normalized','String', Exp.mwFreq,...
                    'Position', [0.82 0.05 0.09 0.05], 'Callback', @freq_change);
    Freq.label = uicontrol('Parent', windowFreq,'Style', 'text','FontSize', 15,...
                    'Units', 'normalized','String', 'GHz',...
                    'Position', [0.92 0.05 0.07 0.05], 'Callback', @freq_change);
    Freq.slider = uicontrol('Parent', windowFreq,'Style', 'slider','Units', 'normalized',...
                  'Position', [0.01 0.05 0.8 0.03],'Min', 0,'Max', 100, 'Value', Exp.mwFreq,...
                  'background',[1.0 0.5 0.5],'Callback', @freq_change2, 'SliderStep', [0.01 0.01]);
              
    plotf();
    
    function freq_change(~, ~)
        Freq.slider.Value = str2double(Freq.edit.String);
        end1 = length(B0);
        B0 = [B0(1):((B0(end1)-B0(1))/(end1-1)):B0(end1)]'*Freq.slider.Value/Exp.mwFreq;
        Exp.mwFreq = Freq.slider.Value;
        mwFreq = Exp.mwFreq;
        Exp.Range = [B0(1) B0(end1)];
        plotf();
    end

    function freq_change2(~, ~)
        Freq.edit.String = Freq.slider.Value;
        freq_change();
    end

    function plotf()
        euler1 = [DS.a(1).a DS.a(2).a DS.a(3).a];  
        euler2 = [DS.a(4).a DS.a(5).a DS.a(6).a];
        
        if DS.E == 0
            [xdata,ydata] = ctssimaniso_git_v3_6_1(data1.Sys1,data1.Sys2,Exp,Dat.gOrientation,Dat.DOrientation,DS.J,DS.D,addbr, euler1, euler2, length(B0)); 
        else
            [xdata,ydata] = ctssimaniso_git_v4(data1.Sys1,data1.Sys2,Exp,Dat.gOrientation,Dat.DOrientation,DS.J,DS.D,DS.E,addbr, euler1, euler2, length(B0));
        end
        ydatatot = ydata{1} + ydata{2} + ydata{3} + ydata{4};
        ydatatotnorm = ydatatot/max(ydatatot) * norm;
   
        if DS.EFP == 0
            [xdatafp,ydatafp] = ctssimaniso_git_initial_v2_1(data1.Sys1,data1.Sys2,Exp,Dat.gOrientationFP,Dat.DOrientationFP,DS.JFP,DS.DFP, addbr, length(B0)); 
        else
            [xdatafp,ydatafp] = ctssimaniso_git_initial_v2_4(data1.Sys1,data1.Sys2,Exp,Dat.gOrientationFP,Dat.DOrientationFP,DS.JFP,DS.DFP,DS.EFP, addbr, length(B0)); 
        end
        ydatatotfp = ydatafp{1} + ydatafp{2} + ydatafp{3} + ydatafp{4};
        ydatatotfpnorm = ydatatotfp/max(ydatatotfp) * normFP;
    
        
        
        cla(Freq.plot);  hold (Freq.plot,'on');
        grid(Freq.plot,'on')

        plot(Freq.plot, B0, spec,'k','LineWidth', 1.5); 
        plot(Freq.plot, xdata, ydatatotnorm,'r','LineWidth', 1.5);
%         plot(Freq.plot, xdata, ydatatotnorm1,'r','LineWidth', 1.5);
        plot(Freq.plot, xdatafp, ydatatotfpnorm,'b', 'LineWidth', 1.5);
%         plot(Freq.plot, xdatafp, ydatatotfpnorm1,'b', 'LineWidth', 1.5);
        
        xlabel(Freq.plot, 'Magnetic field (mT)');
        ylabel(Freq.plot, 'Intensity (a. u.)');
        title(Freq.plot, 'Qband');

        legend(Freq.plot, 'Exp.', 'Fixed relative orientation', 'Full powder'); legend boxoff;
    end
  
end
