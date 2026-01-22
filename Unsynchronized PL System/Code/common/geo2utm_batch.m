function utm_coords = geo2utm_batch(geo_data, i)
% 批量转换：geo_data 是 Nx3 矩阵，每行为 [lat, lon, h]

n = size(geo_data, 1);
utm_coords = zeros(n, 3);

for j = 1:n
    [E, N, U] = geo2utm(geo_data(j,1), geo_data(j,2), geo_data(j,3), i);
    utm_coords(j,:) = [E, N, U];
end

end