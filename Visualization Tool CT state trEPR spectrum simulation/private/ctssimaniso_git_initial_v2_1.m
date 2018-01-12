function [sim_field,simtot] =  ctssimaniso_git_initial_v2_1(Sys1,Sys2,Exp,gOrientation,DOrientation,J,D,additionalbroadening, nB0)
%  ctssimaniso(Sys1,Sys2,Exp,gOrientation,DOrientation,J,D,additionalbroadening)
%  v2.0
%  Last Changes: 01.11.2016
%  Author: Felix Kraffert
 
   if ~isfield(gOrientation,'S2')
       gOrientation.S2 = gOrientation.S1;
   end
   
   if ~exist('nB0', 'var')
       nB0 = 1/2e+3*1e+6*4;
   end
   
   Jiso = J; % MHz
   Jdip = D; % MHz
   simstep = (abs(Exp.Range(2)-Exp.Range(1))*2e+3)/2e+3/(nB0-1);
    sim_field = [Exp.Range(1):simstep:Exp.Range(2)];      % Range for stick spectra for all possible res-fields 
    lw = 1* simstep;               %0.005; % line width for broadening of stick spectra
 
    Popul = [0 1 1 0];                % vector of populations, Singluet precursor populates only state 2 and 3, is only used for the sign of transisitons
    trans =  [1 2;1 3;2 4;3 4];
    alpha = 1.0;                     % 0.4 Gaussian character of line for convolution (only manual simulation), 1 = Gaussian, 0 = Lorentzian
    
     % phi,theta and Weight Combination for all orientations (euler angles)
    phi1= gOrientation.S1.phi;
    theta1 = gOrientation.S1.theta;
    Weights1 = gOrientation.S1.Weight;   
    
    phi2= gOrientation.S2.phi;
    theta2 = gOrientation.S2.theta;
    Weights2 = gOrientation.S2.Weight;  
    
    %theta for summation over all dipolar axes
    phidip = DOrientation.phidip;
    thetadip = DOrientation.thetadip;
    thetaw = DOrientation.Weight; 
    
    

    npoints = length(sim_field);    % calculating lengthes for loops in parfor loop to make it more efficient
    lthetadip = length(thetadip);
    ltheata1 = length(theta1);
    ltheata2 = length(theta2);
       
    clear sim
    
    % define the linewidth of the four transistions from the Spin System 1
    % and 2
    lwbroaden = [Sys1.lw  Sys2.lw Sys2.lw Sys1.lw];

    I = [];
    x = [];

    simtot1 = zeros(1,npoints);
    simtot2 = zeros(1,npoints);
    simtot3 = zeros(1,npoints);
    simtot4 = zeros(1,npoints);

    geff_1 =  zeros(1,ltheata1); %new v2.0
    geff_strain_1= zeros(1,ltheata1);%new v2.0
    dBstrain_1 = zeros(1,ltheata1);%new v2.0


%%    %%
parfor ii = 1:ltheata1 %loop over geff for Sys1   %parfor
    geff_1= geff(Sys1.g,theta1(ii),phi1(ii));
    geff_strain_1 = geff(Sys1.gStrain,theta1(ii),phi1(ii)); %new v2.0  % calculation of the gstrain for the effective g-value
    dBstrain_1 = ((planck*Exp.mwFreq*1e+9/((geff_1)*bmagn))-(planck*Exp.mwFreq*1e+9/((geff_1+geff_strain_1)*bmagn)))*1e3;%new v2.0 %mT calculation of the gstrain in mT
    
    simj12 = zeros(ltheata1,npoints);    %definig Matrices for parfor loop
    simj13 = zeros(ltheata1,npoints);
    simj24 = zeros(ltheata1,npoints);
    simj34 = zeros(ltheata1,npoints);
    
    geff_2 =  zeros(1,ltheata2);%new v2.0
    geff_strain_2= zeros(1,ltheata2);%new v2.0
    dBstrain_2 = zeros(1,ltheata2);%new v2.0
    simsuml = zeros(1,npoints);%new v2.0
    
       
    for Jl = 1:ltheata2 %loop over geff for Sys2
        geff_2= geff(Sys2.g,theta2(Jl),phi2(Jl));  % analytical calculating off geff
        geff_strain_2 = geff(Sys2.gStrain,theta2(Jl),phi2(Jl)); %new v2.0 % calculation of the gstrain for the effective g-value
        dBstrain_2 = ((planck*Exp.mwFreq*1e+9/((geff_2)*bmagn))-(planck*Exp.mwFreq*1e+9/((geff_2+geff_strain_2)*bmagn)))*1e3;%new v2.0 %mT calculation of the gstrain in mT
        
        xsum12 = zeros(1,npoints);
        xsum13 = zeros(1,npoints);
        xsum24 = zeros(1,npoints);
        xsum34 = zeros(1,npoints);
        
        a = 0.5*(geff_1+geff_2)*bmagn/planck*1e-9;   % analytical calculation off alpha and beta for both g_eff values for each combination
        b = 0.5*abs(geff_1-geff_2)*bmagn/planck*1e-9;
        Jiso15 = Jiso/2;  %% Jiso is divided by 2 due to definition of Hamiltonian
        
        for it = 1:lthetadip % loop over dipolar coupling angle theta
            I = Itotal(thetadip(it),geff_1,geff_2,Jiso15,1.5*Jdip,Exp.mwFreq)*thetaw(it);   % Calculation of transistion Intensity in dependece of dipolar and exchange coupling
            
            d = 1.5*Jdip*( (cos(thetadip(it)))^2 - 1/3 );   %% !! d is multiplied by 1.5 due to definition of dzz = 2D/3
            c = (Jiso15+0.5*d);            % used formula1,  from Swaneburg and Hore paper
            
            % Analytical calculation of the four transitions, the
            % selection uses the four not complex solutions of the
            % eight possibillities from the Energy levels
            B12 = a*(Exp.mwFreq*1e3+Jiso15-d)/(a^2-b^2) + sqrt( a^2*(Exp.mwFreq*1e3+Jiso15-d)^2/(a^2-b^2)^2 - ((Exp.mwFreq*1e3+Jiso15-d)^2-c^2)/(a^2-b^2) );
            B13 = a*(Exp.mwFreq*1e3+Jiso15-d)/(a^2-b^2) - sqrt( a^2*(Exp.mwFreq*1e3+Jiso15-d)^2/(a^2-b^2)^2 - ((Exp.mwFreq*1e3+Jiso15-d)^2-c^2)/(a^2-b^2) );
            B24 = a*(Exp.mwFreq*1e3-Jiso15+d)/(a^2-b^2) - sqrt( a^2*(Exp.mwFreq*1e3-Jiso15+d)^2/(a^2-b^2)^2 - ((Exp.mwFreq*1e3-Jiso15+d)^2-c^2)/(a^2-b^2) );
            B34 = a*(Exp.mwFreq*1e3-Jiso15+d)/(a^2-b^2) + sqrt( a^2*(Exp.mwFreq*1e3-Jiso15+d)^2/(a^2-b^2)^2 - ((Exp.mwFreq*1e3-Jiso15+d)^2-c^2)/(a^2-b^2) );
%             B12new=1/(2*geff_1*geff_2*bmagn/planck*1e-9)*(Exp.mwFreq*1e3+Jiso15-d)*((geff_1+geff_2)+sqrt(geff_1*geff_2)*(d+2*Jiso15));
%             [B12, B12new]
%             (geff_1*geff_2)*(d+2*Jiso15)^2/(geff_1-geff_2)^2/(Exp.mwFreq*1e3+Jiso15-d)^2  -  (geff_1*geff_2)*(d+2*Jiso15)^2/(geff_1-geff_2)^2/(Exp.mwFreq*1e3)^2;

            for ik = 1:4 %loop over all resonances for different dipolar
                % couplingto generate the stick spectra with
                % four transitions and the gaussian broadening with lw
                % and the additionalbroadening originate in Heisenberg
                % Broadening, also gstrain is include as dBstrain
                % Popul is used to generate the sign of the
                % transition
                switch ik
                    case 1  %'12'
                        if isreal(B12)
                            
                            if geff_1 <= geff_2  % either flip s1 or s2 depending which g is smaller to avoid wrong
                                
                                xsum12 = xsum12 + gaussian(sim_field,(B12), sqrt((lwbroaden(ik)+dBstrain_1)^2+additionalbroadening^2+lw^2)) *I*(Popul(2)-Popul(1))*Weights2(Jl)*Weights1(ii);
                                
                            else
                                xsum13 = xsum13 + gaussian(sim_field,(B12),sqrt((lwbroaden(2)+dBstrain_2)^2+additionalbroadening^2+lw^2))*I*(Popul(3)-Popul(1))*Weights2(Jl)*Weights1(ii);
                            end
                            
                        end
                    case 2 %'13'
                        if isreal(B13)
                            if geff_1 <= geff_2  % either flip s1 or s2 depending which g is smaller to avoid wrong
                                xsum13 = xsum13 + gaussian(sim_field,(B13) ,sqrt((lwbroaden(ik)+dBstrain_2)^2+additionalbroadening^2+lw^2))*I*(Popul(3)-Popul(1))*Weights2(Jl)*Weights1(ii);
                            else
                                xsum12 = xsum12 + gaussian(sim_field,(B13) ,sqrt((lwbroaden(1)+dBstrain_1)^2+additionalbroadening^2+lw^2))*I*(Popul(2)-Popul(1))*Weights2(Jl)*Weights1(ii);
                            end
                        end
                    case 3 % '24'
                        if isreal(B24)
                            if geff_1 <= geff_2  % either flip s1 or s2 depending which g is smaller to avoid wrong
                                xsum24 = xsum24 + gaussian(sim_field,(B24),sqrt((lwbroaden(ik)+dBstrain_2)^2+additionalbroadening^2+lw^2))*I*(Popul(4)-Popul(2))*Weights2(Jl)*Weights1(ii);
                            else
                                xsum34 = xsum34 + gaussian(sim_field,(B24),sqrt((lwbroaden(4)+dBstrain_1)^2+additionalbroadening^2+lw^2))*I*(Popul(4)-Popul(3))*Weights2(Jl)*Weights1(ii);
                            end
                        end
                    case 4 % '34'
                        if isreal(B34)
                            if geff_1 <= geff_2  % either flip s1 or s2 depending which g is smaller to avoid wrong
                                xsum34 = xsum34 + gaussian(sim_field,(B34),sqrt((lwbroaden(ik)+dBstrain_1)^2+additionalbroadening^2+lw^2))*I*(Popul(4)-Popul(3))*Weights2(Jl)*Weights1(ii);
                            else
                                xsum24 = xsum24 + gaussian(sim_field,(B34),sqrt((lwbroaden(3)+dBstrain_2)^2+additionalbroadening^2+lw^2))*I*(Popul(4)-Popul(2))*Weights2(Jl)*Weights1(ii);
                            end
                        end
                        
                end
                
            end
        end
        simj12(Jl,:) = xsum12;
        simj13(Jl,:) = xsum13;
        simj24(Jl,:) = xsum24;
        simj34(Jl,:) = xsum34;
    end
    sim12{ii} = simj12;
    sim13{ii} = simj13;
    sim24{ii} = simj24;
    sim34{ii} = simj34;
end
sim{1} = sim12;
sim{2} = sim13;
sim{3} = sim24;
sim{4} = sim34;


   
simtot1 = zeros(1,npoints);
simtot2 = zeros(1,npoints);
simtot3 = zeros(1,npoints);
simtot4 = zeros(1,npoints);

% add up all different orientations for all 4 transitions
for l = 1:4
    for li = 1:ltheata1  %loop over geff_1
        simsuml = zeros(1,npoints);
        
        for Jl = 1:ltheata2  % loop over geff_2
            simsuml = simsuml + sim{l}{li}(Jl,:);
        end
        switch l
            
            case 1
                simtot1 = simtot1 + simsuml;
            case 2
                simtot2 = simtot2 + simsuml;
            case 3
                simtot3 = simtot3 + simsuml;
            case 4
                simtot4 = simtot4 + simsuml;
        end
    end
end

simtot{1} = simtot1;
simtot{2} = simtot2;
simtot{3} = simtot3;
simtot{4} = simtot4;
end
 
%% geff %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [geffout] = geff(g,theta,phi,varargin)
% give theta and phi in 'radian' or 'degree'
% the default value is radian
% geff(theta,phi,'radian')

geffout = sqrt(g(1)^2*sin(theta)^2*cos(phi)^2+g(2)^2*sin(theta)^2*sin(phi)^2+g(3)^2*cos(theta)^2);  % g(theta, phi)
end
%% I_total %%%%%%%%%%%%%%%%%%%%%%%%%%

function [I_total] = Itotal(theta,g1_iso,g2_iso,J,D,freq,varargin)
% give theta 'radian' or 'degree'
% the default value is radian
% phi_dip(theta,g1_iso,g2_iso,J,D,'degree')
% give J,D i MHz
% freq in GHz

dw =  abs(g1_iso-g2_iso)*freq/(g1_iso+g2_iso)*2*1000;

% if strcmp(p.Results.mode, 'degree')
%     theta = theta / 180 * pi;
% end
phi = 0.5*atan( dw / ( 2*J + 1/3*D*(3*(cos(theta))^2-1) ) );

I_total = sin(phi)^2*cos(phi)^2;
end