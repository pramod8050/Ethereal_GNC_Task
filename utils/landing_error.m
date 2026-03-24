function error = landing_error(lat, lon, lat_t, lon_t, R_e)

% Haversine formula

dlat = lat - lat_t;
dlon = lon - lon_t;

a = sin(dlat/2).^2 + cos(lat).*cos(lat_t).*sin(dlon/2).^2;
c = 2 * atan2(sqrt(a), sqrt(1 - a));

error = R_e * c;   % distance in meters

end