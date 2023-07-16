function UV_normalized = NormalizeUV(UV,varargin)
if nargin <2
    map = UV;
else
    map=varargin{1};

end
delta_x = max(map(:,1))-min(map(:,1));
delta_y = max(map(:,2))-min(map(:,2));
if delta_x>delta_y
    s = delta_x;
else
    s = delta_y;
end
UV_normalized = (UV - [min(map(:,1)),min(map(:,2))])/s;

% UV_normalized = (map - [min(map(:,1)),min(map(:,2))])/s;
end

