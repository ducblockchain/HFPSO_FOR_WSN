function [xunit,yunit] = circle(x,y,r)
hold on
th = 0:pi/50:2*pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
% h = plot(xunit, yunit);
% hold off
end
function [xunit, yunit, zunit] = sphere3D(x, y, z, r)
    % Generate points on the surface of a sphere with center (x, y, z) and radius r
    [theta, phi] = meshgrid(linspace(0, 2*pi, 50), linspace(0, pi, 25));
    xunit = r * sin(phi) .* cos(theta) + x;
    yunit = r * sin(phi) .* sin(theta) + y;
    zunit = r * cos(phi) + z;
end
% function [xunit, yunit, zunit] = sphere3D(x, y, z, r)
%     Generate points on the surface of a sphere with center (x, y, z) and radius r
%     [theta, phi] = meshgrid(linspace(0, 2*pi, 50), linspace(0, pi, 25));
%     xunit = r * sin(phi) .* cos(theta) + x;
%     yunit = r * sin(phi) .* sin(theta) + y;
%     zunit = r * cos(phi) + z;
% end