function [lat, lon] = ecef_to_latlon(r)

% r: Nx3 matrix or 3x1 vector

x = r(:,1);
y = r(:,2);
z = r(:,3);

R = sqrt(x.^2 + y.^2 + z.^2);

lat = asin(z ./ R);        % radians
lon = atan2(y, x);         % radians

end