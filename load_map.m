%% Создание матрицы необходимой размерности и координатами с информацией о высоте рельефа 
function arr_map = load_map(b_grad, l_grad, size_map_b, size_map_l, step_map_sec)
%Перевод градусов в радианы
b_rad = b_grad*pi/180;
l_rad = l_grad*pi/180;
%Шаг карты в радианах
step_map = pi/180/60/60*step_map_sec;
%Создание массива необходимой размерности
arr_map = zeros(size_map_b, size_map_l);
%Заполнение массива (получение данных из карты)
for i=1:size_map_b
  bb = b_rad + step_map*(i);
  for j=1:size_map_l
    ll = l_rad + step_map*(j);
    arr_map(i,j) = hmap(bb,ll);
  end
end

end