function [value, isterminal, direction] = event_ground(t, x, params)

r = x(1:3);
alt = norm(r) - params.R_e;

value = alt;        % stop when altitude = 0
isterminal = 1;     % stop integration
direction = -1;     % detect decreasing altitude

end