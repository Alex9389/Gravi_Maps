function arr_map = load_map(b_grad, l_grad, size_map_b, size_map_l, step_map_sec)
b_rad = b_grad*pi/180;
l_rad = l_grad*pi/180;
step_map = pi/180/60/60*step_map_sec;
arr_map = zeros(size_map_b, size_map_l);
% i0 = floor(size_map_b/2)+1;
% j0 = floor(size_map_l/2)+1;
for i=1:size_map_b
%   bb = b_rad + step_map*(i-i0);
  bb = b_rad + step_map*(i);
  for j=1:size_map_l
    ll = l_rad + step_map*(j);
%     ll = l_rad + step_map*(j-j0);
    arr_map(i,j) = hmap(bb,ll);
  end
end

end