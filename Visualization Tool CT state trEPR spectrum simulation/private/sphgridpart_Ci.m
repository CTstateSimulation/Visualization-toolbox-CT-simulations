function [phi, theta, weights] = sphgridpart(DeviationAngle, NpointsTheta, NpointsPhi, ErrorValue, Distr)
% [phi, theta, weights] = sphgridpart(DeviationAngle, NpointsTheta, NpointsPhi, ErrorValue, Distr)
%     [phi1, theta1, weights1] = sphgrid('Ci',8);
%     DeviationAngle = pi/2; 
%     NpointsTheta = 10;
%     NpointsPhi = 20;
%     ErrorValue = 10^-4;
%     Distr = 'Equally Distant';
%     phi0 = 0;
%     theta0 = pi/3;
    % Gaussian distribution, get points on line with Gaussian distribution
    %Gaussian1Ddestribution (NpointsTheta, ErrorValue). The integral of gaussian
    %function is 1. And the border is inf. ErrorValue of 10^-6  means that 10^-6 part of points won't be
    %included. From this border is calculated. For error = 10^-6 border
    %is 4,7534. Max border which is possible - 10.
  if ~exist('ErrorValue', 'var')
     ErrorValue = 10^-2;
  end  
  if ~exist('distr', 'var')
     Distr = 'Equally Distant';
  end
  if ~strcmp(Distr,'Equally Distant') && ~strcmp(Distr,'Gaussian')
     error('distriburion should be "Equally Distant" or "Gaussian"')    
  end 
   
  if (NpointsTheta == 0) || (NpointsPhi == 0)
    phi = 0';
    theta = 0';
    weights = 1;
 
  elseif NpointsTheta == 1
    weights = zeros (NpointsPhi, 1); 
    phi = zeros (NpointsPhi, 1); 
    theta = zeros (NpointsPhi, 1); 
    for i = 1:NpointsPhi
        weights(i) = 4*pi/NpointsPhi;
        phi(i) = (i-1) * 2*pi / NpointsPhi;
        theta(i) = DeviationAngle;
    end
    
  else    
    if strcmp (Distr,'Gaussian')
        [ThetaI, ThetaMax]  = Gaussian1Ddestribution (NpointsTheta, ErrorValue); %has values from 0 to some max
    elseif strcmp (Distr,'Equally Distant')
        [ThetaI, ThetaMax]  = onedestribution ( NpointsTheta); 
    end
    
    ThetaI = ThetaI .* DeviationAngle ./ ThetaMax; %now it has values from 0 to deviationangle 
    %y =zeros (NpointsTheta, 1);
    %plot (ThetaI, y, '.')
    
    % Angle for rotation of 1D distribution
    alpha = 2*pi / NpointsPhi;
  
    % Massive of all points and massive of weights
    points = zeros (NpointsPhi * NpointsTheta, 2);
    
    %weight it is sphere segment S=pi(h^2+r^2) - S(of privious segment)
    weights = zeros (NpointsPhi * NpointsTheta, 1); 
    
    ThetaV = zeros (NpointsTheta, 1); %points of voronoi cell to calculate weights
    for j = 1:NpointsTheta-1 
        ThetaV(j) = (ThetaI(j) + ThetaI(j+1))/2;
    end
    ThetaV(NpointsTheta) = DeviationAngle;
    
    [x,~,z] = sph2cart (0, pi/2 - ThetaV, 1);
%     plot(x,z,'.');
    h = 1 - z; 
    r = x;
    
    %area of segment S=pi (h^2 + r^2)
    weights(1) = pi() * (h(1)^2 + r(1)^2) / NpointsPhi * 2;
    
    for j = 2:NpointsTheta 
        weights(j) = abs(pi() * (h(j)^2 + r(j)^2 - h(j-1)^2 - r(j-1)^2) / NpointsPhi) * 2;
    end
        
    points(1, :) = 0;
    for i = 0 : (NpointsPhi - 1)
        for j = 1:NpointsTheta 
            points(i * NpointsTheta + j, 2) = ThetaI(j);
            points(i * NpointsTheta + j, 1) = alpha * i;
            weights(i * NpointsTheta + j) = weights(j);
        end
    end
    
    phi = points(:,1);
    theta = points(:,2);

%       [x,y,z] = sph2cart (phi, pi/2 -theta, 1);
%        x = sin(theta) .* cos(phi);
%        y = sin(theta) .* sin(phi);
%        z = cos(theta);
       
%      [x1,y1,z1] = sph2cart (phi1,pi/2 - theta1, 1);
%      hold on
%         plot3(x,y,z,'+');
%         ylim([-2 2]); xlim([-2 2]); zlim([-2 2]); xlabel('x'); ylabel('y'); zlabel('z')
%      plot3 (x1,y1,z1,'.');
%      hold off
  end
end    
