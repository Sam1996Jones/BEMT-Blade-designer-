%% calculates the CL value depending on the aerofoil selected
% inputs - aerofoil, angle of incidence

%%
function [CL] = CalculateCL(aerofoil,i)

if aerofoil == 1
%     for p = 1:6
%         if i(p) < 14
%             CL(p) = (-0.0028 .* i(p).^2) + (0.1342 .* i(p));       % Derived equation for NACA 0012 profile on excel
%         elseif i(p) == 15                                      % Different equations for sector profiles
%             CL(p) = 0.8;
%         elseif i(p) > 15 && i(p) < 20
%             CL(p) = (0.0064.*i(p).^2) - (0.2449.*i(p)) + 3.0231;
%         else
%             CL(p) = (0.0005.*i(p).^2) - (0.0005.*i(p)) + 0.5075;
%         end
%     end
      CL = 0.084 .* i;
elseif aerofoil == 2
    CL = 0.084 .* i;
end