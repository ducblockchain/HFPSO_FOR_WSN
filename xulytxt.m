% Đọc nội dung file hiện tại
fileID = fopen('fitness_value.txt', 'r');
data = textscan(fileID, 'Iteration %d: Best Cost = %f');
fclose(fileID);

% Trích xuất dữ liệu từ textscan
iterations = data{1};
costs = data{2};

% Mở file để ghi dữ liệu mới
fileID = fopen('fitnessvalue.txt', 'w');

% Ghi dữ liệu mới vào file
for i = 1:length(iterations)
    fprintf(fileID, '%d,%.4f\n', iterations(i), costs(i) );  % Nhân với 100 như yêu cầu
end

% Đóng file
fclose(fileID);
