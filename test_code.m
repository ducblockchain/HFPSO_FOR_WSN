%%

N=100;                       % number of nodes
area=[100,100];              % nodes deployment area in meter
Trange=7;                   % transmission range of sensor node in meter
nodes.pos=area(1).*rand(N,2);% nodes geographical locations
lambda=0.125;                % signal wavelength in meter
nodes.major = Trange;        % major axis for ellpitical range in meter
nodes.minor = lambda*Trange;  % minro axis for ellipitical range in meter
% redundantNo=9;               % number of healing nodes   
redundantNo=round(10*N/100);
%% plot the nodes deployment
cnt=1;
for ii=1:N      
    for jj=1:N
        if ii~=jj
            nodes.distance(ii,jj)=pdist([nodes.pos(ii,:);nodes.pos(jj,:)]);
            if nodes.distance(ii,jj)<Trange || nodes.distance(ii,jj)==Trange
                nodes.inrange(ii,jj)=1;
            else
                nodes.inrange(ii,jj)=0;
            end
        end
    end
end
figure
F5=plot(nodes.pos(:,1),nodes.pos(:,2),'.','color','r');
hold on
for ii=1:N                   % plot the circular transmission range
    [nodes.circle.x(ii,:),nodes.circle.y(ii,:)]=circle(nodes.pos(ii,1),nodes.pos(ii,2),Trange);
    F6=fill(nodes.circle.x(ii,:),nodes.circle.y(ii,:),[0.25,0.25,0.25]);
    alpha 0.3
    hold on
end
axis on
xlabel('x(m)')
ylabel('y(m)')
title('Initial Placement of Nodes with circular transmission range')
%% plot delauny triangle
TRI = delaunay(nodes.pos(:,1),nodes.pos(:,2));
figure(2)
F5 = plot(nodes.pos(:,1),nodes.pos(:,2),'.','color','r');
hold on
for ii=1:N                   % plot the circular transmission range
    [nodes.circle.x(ii,:),nodes.circle.y(ii,:)]=circle(nodes.pos(ii,1),nodes.pos(ii,2),Trange);
    F6=fill(nodes.circle.x(ii,:),nodes.circle.y(ii,:),[0.25,0.25,0.25]);
    alpha 0.3
    hold on
end
axis on
xlabel('x(m)')
ylabel('y(m)')
title('Coverage hole in initila position of Nodes')
hold on
triplot(TRI,nodes.pos(:,1),nodes.pos(:,2))
%% Hole detection
[holeDetected.circle,Circmcenter.circle,circumradius.circle]=holeDetection(TRI,nodes,F5,F6,Trange,area,2,1);
display(['--> No of detected Holes for Circular = ',num2str(numel(find(holeDetected.circle)))])
%% PSO optimize position of rest wsn nodes to cover the hole
nvars = 2 * 100; 
numParticles = 30; % Số lượng hạt
maxIterations = 1000; % Số lần lặp tối đa
w = 0.72; % Trọng số quán tính
c1 = 1.49; % Hệ số nhận thức
c2 = 1.49; % Hệ số xã hội

% Giới hạn dưới và trên của các biến
lb = zeros(nvars, 1);
ub = area(1) .* ones(nvars, 1);

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
finalPos = reshape(globalBestPosition, [numel(globalBestPosition) / 2, 2]);

T = table(finalPos(:,1), finalPos(:,2));
writetable(T, 'somefile.txt')
figure
plot(finalPos(:,1),finalPos(:,2),'.','color','r');
hold on
for ii=1:100                 % plot the circular transmission range
    [finalcircle.x(ii,:),finalcircle.y(ii,:)]=circle(finalPos(ii,1),finalPos(ii,2),Trange);
    fill(finalcircle.x(ii,:),finalcircle.y(ii,:),[0.25,0.25,0.25]);
    alpha 0.3
    hold on
end
axis on
xlabel('x(m)')
ylabel('y(m)')
title('Optimized location of Nodes with circular transmission range')

tp1;