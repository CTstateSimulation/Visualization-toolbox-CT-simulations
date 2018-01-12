%Version the same as V3_6 but with addition: Number of points for B0

%Initial order of loops. The idea: rotation of 1st spin (theta1), than of the 2nd one (theta2), rotation of the system all together (thetadip)
%geff(theta + thetadip, phi + phidip)


%  Last Changes: 15.12.2016 (the only right working version on this date)
%                12.10.2017 E is added
%  Author: Felix Kraffert, Daria Dymnikova

function [sim_field,simtot] =  ctssimaniso_git_v4(Sys1,Sys2,Exp,gOrientation,DOrientation,J,D,E,additionalbroadening,euler1, euler2, nB0)
   if ~isfield(gOrientation,'S2')
       gOrientation.S2 = gOrientation.S1;
   end
   
   if ~exist('nB0', 'var')
       nB0 = 1/2e+3*1e+6*4;
   end
   
    if ~exist('distr', 'var')
       distr = 'Equally Distant';
    end
   
   if ~strcmp(distr,'Equally Distant') && ~strcmp(distr,'Gaussian')
      error('distriburion should be "Equally Distant" or "Gaussian"')    
   end
   
   Sys1.g = eulerR(euler1)*[Sys1.g(1) 0 0; 0 Sys1.g(2) 0; 0 0 Sys1.g(3)]*(eulerR(euler1)');
   Sys2.g = eulerR(euler2)*[Sys2.g(1) 0 0; 0 Sys2.g(2) 0; 0 0 Sys2.g(3)]*(eulerR(euler2)');
   Sys1.gStrain = eulerR(euler1)*[Sys1.gStrain(1) 0 0; 0 Sys1.gStrain(2) 0; 0 0 Sys1.gStrain(3)]*(eulerR(euler1)');
   Sys2.gStrain = eulerR(euler2)*[Sys2.gStrain(1) 0 0; 0 Sys2.gStrain(2) 0; 0 0 Sys2.gStrain(3)]*(eulerR(euler2)');
   
%    X1 = Sys1.g
%    X2 = Sys2.g
%    v1=X1(:,3);
%    v2=X2(:,3);
%    v1 = v1/norm(v1);
%    v2 = v2/norm(v2);
%    angle =  atan2(norm(cross(v1,v2)),dot(v1,v2));
%    axis = cross(v1,v2)/norm(cross(v1,v2));
%    K = [ 0 -axis(3) axis(2) ; axis(3) 0 -axis(1) ; -axis(2) axis(1) 0];
%    R = eye(3) + sin(angle)*K + (1-cos(angle))*K*K;
%    v2-R*v1
%    X1=R * X1;
%    X1
%    X2
%    phiD = subspace(X1(:,1), X2(:,1))
%    phiPROV = subspace(X1(:,2), X2(:,2))
   
   Jiso = J; % MHz
   Jdip = D; % MHz
   JdipE = E; % MHz
   simstep =  (abs(Exp.Range(2)-Exp.Range(1))*2e+3)/2e+3/(nB0-1);   % resolution of sim_field
   
    sim_field = [Exp.Range(1):simstep:Exp.Range(2)];      % Range for stick spectra for all possible res-fields 
    lw = 5* simstep;               %0.005; % line width for broadening of stick spectra
 
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

%     PointsConsid2 = zeros(lthetadip, 4);
%     FieldG = zeros(lthetadip,6);
%     lg = 1;

    normalizCoeff = norminv(1-10^-2)^2 / max(theta1) / max(theta2);

    l = 0;
   for i = 1:ltheata1
      if theta1(i) == 0  || theta1(i) == pi/2
        l = l + 1;   
      else  
        Weights1(i) = Weights1(i) / 2;
        Weights1(ltheata1 + i - l) = Weights1(i);
        theta1(ltheata1 + i - l) = -theta1(i); 
        phi1(ltheata1 + i - l) = phi1(i);
      end
   end


   l = 0;
   for i = 1:ltheata2
      if theta2(i) == 0 || theta2(i) == pi/2
        l = l + 1;  
      else
        Weights2(i) = Weights2(i) / 2;
        Weights2(ltheata1 + i - l) = Weights2(i);
        theta2(ltheata2 + i - l) = -theta2(i); 
        phi2(ltheata2 + i - l) = phi2(i);
      end
   end
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
%    lpoints=1;   

parfor it = 1:lthetadip  %loop over geff for Sys1   %parfor
    
    simj12 = zeros(ltheata1,npoints);    %definig Matrices for parfor loop
    simj13 = zeros(ltheata1,npoints);
    simj24 = zeros(ltheata1,npoints);
    simj34 = zeros(ltheata1,npoints);
    
    geff_2 =  zeros(1,ltheata2);%new v2.0
    geff_strain_2= zeros(1,ltheata2);%new v2.0
    dBstrain_2 = zeros(1,ltheata2);%new v2.0
    simsuml = zeros(1,npoints);%new v2.0
     
    for ii = 1:ltheata1%loop over geff for Sys2
      
        xsum12 = zeros(1,npoints);
        xsum13 = zeros(1,npoints);
        xsum24 = zeros(1,npoints);
        xsum34 = zeros(1,npoints);

        for Jl = 1:ltheata2 % loop over dipolar coupling angle theta
          
            X1 = [theta1(ii) + thetadip(it), phi1(ii) + phidip(it)];
            X2 = [theta2(Jl) + thetadip(it), phi2(Jl) + phidip(it)];

            if X1(1) < 0
                X1(1) = -X1(1);
            end
            if X1(1) >= pi
                X1(1) = X1(1) - pi;
            end
            if X1(1) > pi/2
                X1(1) = pi - X1(1);
            end
            
            if X2(1) < 0
                X2(1) = -X2(1);
            end
            if X2(1) >= pi
                X2(1) = X2(1) - pi;
            end
            if X2(1) > pi/2
                X2(1) = pi - X2(1);
            end

            if length(theta1(1,:))== 1
                k1 = dsearchn(theta1,X1(1));
                k2 = dsearchn(theta2,X2(1));
            else
                k1 = dsearchn(theta1',X1(1));
                k2 = dsearchn(theta2',X2(1));
            end
%             
%             if length(thetadip(1,:))== 1
%                 kd1 = dsearchn(thetadip,X1(1));
%                 kd2 = dsearchn(thetadip,X2(1));
%             else
%                 kd1 = dsearchn(thetadip',X1(1));
%                 kd2 = dsearchn(thetadip',X2(1));
%             end
%             
%             
%             PointsConsid2(kd1, 1) = PointsConsid2(kd1, 1) + 1;
%             PointsConsid2(kd2, 2) = PointsConsid2(kd2, 2) + 1;
%             PointsConsid2(it, 3) = PointsConsid2(it, 3) + 1;
%             PointsConsid2(it, 4) = PointsConsid2(it, 4) + Weights1(k1)*Weights2(k2)*thetaw(it);
            
%             pointsused2(lpoints,:) = [thetadip(it),theta1(ii) + thetadip(it), phi1(ii) + phidip(it), theta2(Jl) + thetadip(it), phi2(Jl) + phidip(it)];
%             lpoints = lpoints+1;
            
            geff_1 = geff(Sys1.g,X1(1),X1(2));
            geff_strain_1 = geff(Sys1.gStrain,X1(1),X1(2)); %new v2.0  % calculation of the gstrain for the effective g-value
            dBstrain_1 = ((planck*Exp.mwFreq*1e+9/((geff_1)*bmagn))-(planck*Exp.mwFreq*1e+9/((geff_1+geff_strain_1)*bmagn)))*1e3;%new v2.0 %mT calculation of the gstrain in mT
    
            geff_2 = geff(Sys2.g,X2(1),X2(2));  % analytical calculating off geff
            geff_strain_2 = geff(Sys2.gStrain,X2(1),X2(2)); %new v2.0 % calculation of the gstrain for the effective g-value
            dBstrain_2 = ((planck*Exp.mwFreq*1e+9/((geff_2)*bmagn))-(planck*Exp.mwFreq*1e+9/((geff_2+geff_strain_2)*bmagn)))*1e3;%new v2.0 %mT calculation of the gstrain in mT
         
            a = 0.5*(geff_1+geff_2)*bmagn/planck*1e-9;   % analytical calculation off alpha and beta for both g_eff values for each combination
            b = 0.5*abs(geff_1-geff_2)*bmagn/planck*1e-9;
            Jiso15 = Jiso/2;  %% Jiso is divided by 2 due to definition of Hamiltonian
            d = 1.5* Jdip*( (cos(thetadip(it)))^2 - 1/3 ) + JdipE*(sin(thetadip(it)))^2*(2*(cos(phidip(it)))^2-1); %% !! d is multiplied by 1.5 due to definition of dzz = 2D/3
            I = Itotal(thetadip(it),geff_1,geff_2,Jiso15,1.5*D,Exp.mwFreq)*thetaw(it);   % Calculation of transistion Intensity in dependece of dipolar and exchange coupling
            
            if strcmp(distr,'Gaussian')
                I = I * normcdf(sqrt(theta1(ii)^2 + phi1(ii)^2)) * normcdf(sqrt(theta2(Jl)^2 + phi2(Jl)^2)) * normalizCoeff;
                %norminv(1-Error)^2 / max(theta1) / max(theta2) is needed
                %for normalization Error determines the border of Gaussian
                %distribution
            end
            
            c = (Jiso15+0.5*d);            % used formula1,  from Swaneburg and Hore paper
            
            % Analytical calculation of the four transitions, the
            % selection uses the four not complex solutions of the
            % eight possibillities from the Energy levels
            B12 = a*(Exp.mwFreq*1e3+Jiso15-d)/(a^2-b^2) + sqrt( a^2*(Exp.mwFreq*1e3+Jiso15-d)^2/(a^2-b^2)^2 - ((Exp.mwFreq*1e3+Jiso15-d)^2-c^2)/(a^2-b^2) );
            B13 = a*(Exp.mwFreq*1e3+Jiso15-d)/(a^2-b^2) - sqrt( a^2*(Exp.mwFreq*1e3+Jiso15-d)^2/(a^2-b^2)^2 - ((Exp.mwFreq*1e3+Jiso15-d)^2-c^2)/(a^2-b^2) );
            B24 = a*(Exp.mwFreq*1e3-Jiso15+d)/(a^2-b^2) - sqrt( a^2*(Exp.mwFreq*1e3-Jiso15+d)^2/(a^2-b^2)^2 - ((Exp.mwFreq*1e3-Jiso15+d)^2-c^2)/(a^2-b^2) );
            B34 = a*(Exp.mwFreq*1e3-Jiso15+d)/(a^2-b^2) + sqrt( a^2*(Exp.mwFreq*1e3-Jiso15+d)^2/(a^2-b^2)^2 - ((Exp.mwFreq*1e3-Jiso15+d)^2-c^2)/(a^2-b^2) );
%             
%             FieldG(lg, :) = [geff_1, geff_2, B12, B13, B24, B34];
%             lg = lg+1;
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
                                xsum12 = xsum12 + gaussian(sim_field,(B12),sqrt((lwbroaden(ik)+dBstrain_1)^2+additionalbroadening^2+lw^2))*I*(Popul(2)-Popul(1))*Weights2(k2)*Weights1(k1);
                            else
                                xsum13 = xsum13 + gaussian(sim_field,(B12),sqrt((lwbroaden(2)+dBstrain_2)^2+additionalbroadening^2+lw^2))*I*(Popul(3)-Popul(1))*Weights2(k2)*Weights1(k1);
                            end
                        end
                    case 2 %'13'
                        if isreal(B13)
                            if geff_1 <= geff_2  % either flip s1 or s2 depending which g is smaller to avoid wrong
                                xsum13 = xsum13 + gaussian(sim_field,(B13),sqrt((lwbroaden(ik)+dBstrain_2)^2+additionalbroadening^2+lw^2))*I*(Popul(3)-Popul(1))*Weights2(k2)*Weights1(k1);
                            else
                                xsum12 = xsum12 + gaussian(sim_field,(B13),sqrt((lwbroaden(1)+dBstrain_1)^2+additionalbroadening^2+lw^2))*I*(Popul(2)-Popul(1))*Weights2(k2)*Weights1(k1);
                            end
                        end
                    case 3 % '24'
                        if isreal(B24)
                            if geff_1 <= geff_2  % either flip s1 or s2 depending which g is smaller to avoid wrong
                                xsum24 = xsum24 + gaussian(sim_field,(B24),sqrt((lwbroaden(ik)+dBstrain_2)^2+additionalbroadening^2+lw^2))*I*(Popul(4)-Popul(2))*Weights2(k2)*Weights1(k1);
                            else
                                xsum34 = xsum34 + gaussian(sim_field,(B24),sqrt((lwbroaden(4)+dBstrain_1)^2+additionalbroadening^2+lw^2))*I*(Popul(4)-Popul(3))*Weights2(k2)*Weights1(k1);
                            end
                        end
                    case 4 % '34'
                        if isreal(B34)
                            if geff_1 <= geff_2  % either flip s1 or s2 depending which g is smaller to avoid wrong
                                xsum34 = xsum34 + gaussian(sim_field,(B34),sqrt((lwbroaden(ik)+dBstrain_1)^2+additionalbroadening^2+lw^2))*I*(Popul(4)-Popul(3))*Weights2(k2)*Weights1(k1);
                            else
                                xsum24 = xsum24 + gaussian(sim_field,(B34),sqrt((lwbroaden(3)+dBstrain_2)^2+additionalbroadening^2+lw^2))*I*(Popul(4)-Popul(2))*Weights2(k2)*Weights1(k1);
                            end
                        end
                end  
            end
%  end
        end
            simj12(ii,:) = xsum12;
            simj13(ii,:) = xsum13;
            simj24(ii,:) = xsum24;
            simj34(ii,:) = xsum34;
    end
        sim12{it} = simj12;
        sim13{it} = simj13;
        sim24{it} = simj24;
        sim34{it} = simj34;
end
sim{1} = sim12;
sim{2} = sim13;
sim{3} = sim24;
sim{4} = sim34;
% 
% 
% save ('PointsConsid2.mat', 'PointsConsid2');
% save('pointsused2.mat', 'pointsused2');
% save('FieldG.mat', 'FieldG');

simtot1 = zeros(1,npoints);
simtot2 = zeros(1,npoints);
simtot3 = zeros(1,npoints);
simtot4 = zeros(1,npoints);

% add up all different orientations for all 4 transitions
for l = 1:4
    for it = 1:lthetadip  %loop over geff_1
        simsuml = zeros(1,npoints);
        
        for ii = 1:ltheata1  % loop over geff_2
            simsuml = simsuml + sim{l}{it}(ii,:);
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

function [geffout] = geff(g,theta,phi)
% give theta and phi in 'radian' or 'degree'
% the default value is radian
% geff(theta,phi,'radian')
B0 = ang2vec(phi,theta);
%geffout = sqrt(g(1)^2*sin(theta)^2*cos(phi)^2+g(2)^2*sin(theta)^2*sin(phi)^2+g(3)^2*cos(theta)^2);  % g(theta, phi)
geffout = sqrt(B0'*(g*g')*B0);  %should be B0*(g*g')*B0', but Bo - colomn, so it is right
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

%rotz(phi) rotz - rotation around z axis
function [r] = rotx(a)
r = [1 0 0; 0 cos(a) sin(a); 0 -sin(a) cos(a)];
end
function [r] = roty(a)
r = [cos(a) 0 -sin(a); 0 1 0; sin(a) 0 cos(a)];
end
function [r] = rotz(a)
r = [cos(a) sin(a) 0; -sin(a) cos(a) 0; 0 0 1];
end
function [r] = eulerR(a)
r = rotz(a(3))*roty(a(2))*rotz(a(1));
% r2 = rotz(a(1))*roty(a(2))*rotz(a(3));
% if ~isequal(r,r2)
%     error([num2str(a(1)), num2str(a(2)), num2str(a(3))])
% end
end