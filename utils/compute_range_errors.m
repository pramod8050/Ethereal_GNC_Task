function [downrange, crossrange] = compute_range_errors(r_current, lat_t, lon_t, R_e)

% Target position in ECEF
r_target = R_e * [cos(lat_t)*cos(lon_t);
                  cos(lat_t)*sin(lon_t);
                  sin(lat_t)];

% Local east direction at target
east = [-sin(lon_t);
         cos(lon_t);
         0];

% Local north direction at target
north = [-sin(lat_t)*cos(lon_t);
         -sin(lat_t)*sin(lon_t);
          cos(lat_t)];

% Position error in ECEF
delta_r = r_current - r_target;

% Project onto tangent-plane axes
crossrange = dot(delta_r, east);
downrange  = dot(delta_r, north);

end