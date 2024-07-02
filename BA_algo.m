function CH_array = BA_algo(SN, centroids)
    meas2 = [vertcat(SN(:).x), vertcat(SN(:).y)];
    dataset_size = size(meas2);
    meas = [];
    for i = 1:dataset_size
        if SN(i).E > 0
            meas = [meas; meas2(i, :), SN(i).id];
        end
    end
    dataset_size = size(meas);

    % INIT BAT ALGORITHM PARAMETERS
    centroids = max(1, centroids);    % số cụm (centroids)
    dimensions = 2;                   % số chiều của mỗi centroid
    numBats = 20;                     % số dơi
    maxGenerations = 50;              % số lần lặp tối đa
    A = 0.5;                          % biên độ điều chỉnh
    r = 0.5;                          % tỷ lệ phát xung
    Qmin = 0;                         % tần số tối thiểu
    Qmax = 2;                         % tần số tối đa

    % Khởi tạo các dơi
    bats = repmat(min(meas(:, [1, 2])), numBats, 1) + (repmat(max(meas(:, [1, 2])) - min(meas(:, [1, 2])), numBats, 1)) .* rand(numBats, 2, centroids);
    velocities = zeros(numBats, 2, centroids);
    fitness = inf(numBats, 1);

    % Điểm tốt nhất toàn cục
    globalBestPosition = zeros(2, centroids);
    globalBestScore = inf;

    for gen = 1:maxGenerations
        Q = Qmin + (Qmax - Qmin) .* rand(numBats, 1);
        for i = 1:numBats
            % Cập nhật vị trí và vận tốc của các dơi
            velocities(i, :, :) = velocities(i, :, :) + (bats(i, :, :) - globalBestPosition) .* Q(i);
            newPosition = bats(i, :, :) + velocities(i, :, :);

            % Áp dụng việc dò tìm local
            if rand > r
                newPosition = globalBestPosition + 0.001 * randn(1, 2, centroids);
            end

            % Đảm bảo các hạt trong phạm vi giới hạn
            newPosition = max(newPosition, min(meas(:, [1, 2])));
            newPosition = min(newPosition, max(meas(:, [1, 2])));

            % Tính toán hàm fitness
            distances = zeros(dataset_size(1), centroids, 1);
            for centroid = 1:centroids
                for data_vector = 1:dataset_size(1)
                    distances(data_vector, centroid, 1) = norm(newPosition(centroid, :) - meas(data_vector, [1, 2]));
                end
            end
            c = zeros(dataset_size(1), 1);
            [~, index] = min(distances(:, :, 1), [], 2);
            c(:, 1) = index;
            newScore = fitness_function_for_CH(SN, distances, floor(meas(:, 3)), c, 1, centroid);

            % Chấp nhận các giải pháp mới với xác suất A
            if (newScore < fitness(i)) && (rand < A)
                bats(i, :, :) = newPosition;
                fitness(i) = newScore;
            end

            % Cập nhật điểm tốt nhất toàn cục
            if newScore < globalBestScore
                globalBestScore = newScore;
                globalBestPosition = newPosition;
            end
        end

        % Hiển thị thông tin quá trình
        disp(['Generation ', num2str(gen), ': Best Score = ', num2str(globalBestScore)]);
    end

    % Lấy vị trí tối ưu cuối cùng
    finalPos = globalBestPosition;

    % Gán Cluster Head (CH) dựa trên vị trí tối ưu cuối cùng
    distances = zeros(dataset_size(1), centroids, 1);
    for centroid = 1:centroids
        for data_vector = 1:dataset_size(1)
            distances(data_vector, centroid, 1) = norm(finalPos(centroid, :) - meas(data_vector, [1, 2]));
        end
    end
    [~, ind] = min(distances(:, :, 1));
    CH_array = floor(meas(ind, 3));
end
