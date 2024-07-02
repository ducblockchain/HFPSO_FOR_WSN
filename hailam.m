N = 100;                       % number of nodes
area = [100, 100, 100];         % nodes deployment area in meter (3D area)
Trange = 7;                    % transmission range of sensor node in meter
nodes.pos = [area(1) .* rand(N, 1), area(2) .* rand(N, 1), area(3) .* rand(N, 1)]; % nodes geographical locations in 3D
lambda = 0.125;                % signal wavelength in meter
nodes.major = Trange;          % major axis for elliptical range in meter
nodes.minor = lambda * Trange; % minor axis for elliptical range in meter
redundantNo = round(10 * N / 100); % number of healing nodes

% Calculate distances and check transmission range
cnt = 1;
for ii = 1:N
    for jj = 1:N
        if ii ~= jj
            nodes.distance(ii, jj) = pdist([nodes.pos(ii, :); nodes.pos(jj, :)]);
            if nodes.distance(ii, jj) < Trange || nodes.distance(ii, jj) == Trange
                nodes.inrange(ii, jj) = 1;
            else
                nodes.inrange(ii, jj) = 0;
            end
        end
    end
end

% Plot the nodes deployment in 3D
figure
F5 = plot3(nodes.pos(:, 1), nodes.pos(:, 2), nodes.pos(:, 3), '.', 'color', 'r');
hold on
for ii = 1:N % plot the spherical transmission range
    [nodes.sphere.x{ii}, nodes.sphere.y{ii}, nodes.sphere.z{ii}] = sphere3D(nodes.pos(ii, 1), nodes.pos(ii, 2), nodes.pos(ii, 3), Trange);
    surf(nodes.sphere.x{ii}, nodes.sphere.y{ii}, nodes.sphere.z{ii}, 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'FaceColor', [0.25, 0.25, 0.25]);
    hold on
end
axis on
xlabel('x(m)')
ylabel('y(m)')
zlabel('z(m)')
title('Initial Placement of Nodes with spherical transmission range')
grid on
%PSO
nvars = 3 * N; 
numParticles = 30; % Số lượng hạt
maxIterations = 1000; % Số lần lặp tối đa
w = 0.72; % Trọng số quán tính
c1 = 1.49; % Hệ số nhận thức
c2 = 1.49; % Hệ số xã hội
% Giới hạn dưới và trên của các biến
lb = zeros(nvars, 1);
ub = [area(1).*ones(N,1); area(2).*ones(N,1); area(3).*ones(N,1)];

% Khởi tạo các hạt
particlePositions = repmat(lb', numParticles, 1) + (repmat(ub', numParticles, 1) - repmat(lb', numParticles, 1)) .* rand(numParticles, nvars);
particleVelocities = zeros(numParticles, nvars);
personalBestPositions = particlePositions;
personalBestScores = inf(numParticles, 1);

% Khởi tạo điểm tốt nhất toàn cục
globalBestPosition = zeros(1, nvars);
globalBestScore = inf;

for iter = 1:maxIterations
    for i = 1:numParticles
        % Tính toán hàm fitness
        score = fitness_function_for_position_test(particlePositions(i, :), Trange, area);

        % Cập nhật điểm tốt nhất cá nhân
        if score < personalBestScores(i)
            personalBestScores(i) = score;
            personalBestPositions(i, :) = particlePositions(i, :);
        end

        % Cập nhật điểm tốt nhất toàn cục
        if score < globalBestScore
            globalBestScore = score;
            globalBestPosition = particlePositions(i, :);
        end
    end

    % Cập nhật vận tốc và vị trí của các hạt
    for i = 1:numParticles
        r1 = rand(1, nvars);
        r2 = rand(1, nvars);
        cognitiveComponent = c1 * r1 .* (personalBestPositions(i, :) - particlePositions(i, :));
        socialComponent = c2 * r2 .* (globalBestPosition - particlePositions(i, :));
        particleVelocities(i, :) = w * particleVelocities(i, :) + cognitiveComponent + socialComponent;
        particlePositions(i, :) = particlePositions(i, :) + particleVelocities(i, :);

        % Giữ các hạt trong phạm vi giới hạn
        particlePositions(i, :) = max(particlePositions(i, :), lb');
        particlePositions(i, :) = min(particlePositions(i, :), ub');
    end

    % Hiển thị thông tin quá trình
    disp(['Iteration ', num2str(iter), ': Best Score = ', num2str(globalBestScore)]);
end

% Lấy vị trí tối ưu cuối cùng
finalPos = reshape(globalBestPosition, [numel(globalBestPosition) / 3, 3]);

T = table(finalPos(:,1), finalPos(:,2), finalPos(:,3));
writetable(T, 'somefile.txt')

figure
plot3(finalPos(:,1), finalPos(:,2), finalPos(:,3), '.', 'color', 'r');
hold on
for ii = 1:N % plot the spherical transmission range
    [sphereX, sphereY, sphereZ] = sphere3D(finalPos(ii, 1), finalPos(ii, 2), finalPos(ii, 3), Trange);
    surf(sphereX, sphereY, sphereZ, 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'FaceColor', [0.25, 0.25, 0.25]);
    hold on
end
axis on
xlabel('x(m)')
ylabel('y(m)')
zlabel('z(m)')
title('Optimized location of Nodes with spherical transmission range')
grid on

