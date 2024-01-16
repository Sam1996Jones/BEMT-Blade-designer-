%%function to calculate coefficient of power, axial induction factor,
%angular induction factor, Beta, angle of incidence, thrust and torque
%inputs are:
% number of blades, maximum radius, tip speed ratio, wind velocity, segment
% radii, gamma, number of aerofoil segments and if tip loss should be
% included (Y/N)


%%
function [Cp,a,a_,BetaDeg,i,Thrust,Torque] = DesignCalcs(NumBlades,Radius,lambda2,Vel,r2,c2,Gam2,NumAero,TipLoss,AeroType)

B = NumBlades;              % Number of blades;
R = Radius;              % Rotor radius
lambda = lambda2;         % Tip speed ratio
V = Vel;              % Expected wind velocity
Omega = (lambda*V)/R;       % Rotating blade speed
Q = 1;                  %No tip loss initially
%% First / Initial Guess
r = r2;              % Segment radius
c = c2;           % Chord Length
Gamma = Gam2;       % Gamma

lambda_r = (Omega .* r) ./ V;
Sigma_ = (B .* c) ./ (2 .* pi .* r);            % Local solidity
Beta =(pi/2) - ( (2/3) * atan(1./lambda_r) ) ;      %rads
i = Gamma - rad2deg(Beta);       % Angle of incidence in degrees

% if TipLoss == 'Y'            %calculate tip loss correction
%     Q = (2/pi).*acos(exp(-( (B./(2.*(1-r/R)))./((r/R).*cos(Beta)))));
% end

CL = 1; %initial guess
a = (1 + (4 .* cos(Beta).^2)./(Sigma_ .* CL .* sin(Beta))).^-1;     % Axial induction factor
a_ = (1 - (3.*a))./((4.*a)-1);       % Angular induction factor

%% iteration
olda(5) = 0;
olda_(5) = 0;
counter = 1;
while ( abs((round(a_(1),5))-(round(olda_(1),5) )) ~=0 )
    olda = a;
    olda_ = a_;
    Beta = atan( (lambda_r .* (1 + a_))./(1-a));
    BetaDeg = rad2deg(Beta);
    i = Gamma - rad2deg(Beta);
    CL = CalculateCL(AeroType,i);
    a = (1 + (4 .* cos(Beta).^2)./(Sigma_ .* CL .* sin(Beta))).^-1; 
    a_ = ( (Sigma_ .* CL)./(4.*lambda_r.*cos(Beta))) .* (1-a);
    if counter > 500
        break
    end
    counter = counter + 1;
end

if TipLoss == 'Y'            %calculate tip loss correction
    Q = (2/pi).*acos(exp(-( (B./2.*((1-r/R)))./((r/R).*cos(Beta)))));
end

%% output
for z = 1:NumAero
    integral(z) = ( (lambda_r(z+1)-lambda_r(z)) ./2 ) * (   ((lambda_r(z).^3).*a_(z).*(1 - a(z))) + ((lambda_r(z+1).^3).*a_(z+1).*(1 - a(z+1))) ); 
end
area= sum(integral);
Cp = 8./lambda.^2 .*area;

for z = 1: NumAero
    ThrustIntegral(z) = ( (r(z+1)-r(z))./2) * ( (4.*r(z).*a(z).*(1 - a(z))) + (4.*r(z+1).*a(z+1).*(1 - a(z+1))) ); 
end
ThrustArea = sum(ThrustIntegral);
%Thrust = (pi .* 1000 .* V.^2 .* ThrustArea)./(1000.*B); %kN per blade
Thrust = (pi .* 1000 .* V.^2 .* ThrustArea)./1000; %kN Total

for z = 1: NumAero
    TorqueIntegral(z) = ( (r(z+1)-r(z))./2) * ( (4.* r(z).^3 .*a_(z).*(1 - a(z))) + (4.*r(z+1).^3.*a_(z+1).*(1 - a(z+1))) ); 
end
TorqueArea = sum(TorqueIntegral);
%Torque = (pi .* 1000 .* V .* Omega .*TorqueArea)./(1000.*B); %kNm per blade
Torque = (pi .* 1000 .* V .* Omega .*TorqueArea)./1000; %kNm Total