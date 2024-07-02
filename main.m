%% Khởi tạo các tham số
N = 100;                       % số lượng nút
area = [120, 120];             % khu vực triển khai các nút tính bằng mét
Trange = 7;                    % phạm vi truyền của nút cảm biến tính bằng mét
nodes.pos = area(1) * rand(N, 2); % vị trí địa lý của các nút
lambda = 0.125;                % bước sóng tín hiệu tính bằng mét
nodes.major = Trange;          % trục chính cho phạm vi elip tính bằng mét
nodes.minor = lambda * Trange; % trục nhỏ cho phạm vi elip tính bằng mét
redundantNo = round(10 * N / 100); % số lượng nút dự phòng
data = readtable('random.txt');
nodes.pos = [data.Var1, data.Var2];
%% Vẽ sơ đồ triển khai các nút
cnt = 1;
for ii = 1:N
    for jj = 1:N
        if ii ~= jj
            nodes.distance(ii, jj) = pdist([nodes.pos(ii, :); nodes.pos(jj, :)]);
            if nodes.distance(ii, jj) <= Trange
                nodes.inrange(ii, jj) = 1;
            else
                nodes.inrange(ii, jj) = 0;
            end
        end
    end
end
figure;
plot(nodes.pos(:, 1), nodes.pos(:, 2), '.', 'color', 'r');
hold on;
for ii = 1:N % vẽ phạm vi truyền hình tròn
    [nodes.circle.x(ii, :), nodes.circle.y(ii, :)] = circle(nodes.pos(ii, 1), nodes.pos(ii, 2), Trange);
    fill(nodes.circle.x(ii, :), nodes.circle.y(ii, :), [0.25, 0.25, 0.25]);
    alpha(0.3);
    hold on;
end
axis on;
xlabel('x(m)');
ylabel('y(m)');
title('Vị trí ban đầu ngẫu nhiên của các nút cảm biến với phạm vi đo hình tròn');

%% Sử dụng thuật toán HFPSO để tối ưu hóa vị trí các nút
nVars = 2 * N; % Số lượng biến (x và y của mỗi nút)
fitness_function = @(x) fitness_function_for_position(x, Trange, area);
LB = zeros(1, nVars);
UB = ones(1, nVars) * max(area);
iter = 4000;
swarm_size = 30;
c1 = 1.49445;
c2 = 1.49445;
vmax_coef = 0.1;
fitness_over_time = zeros(iter, 1);
[f_valc, xc] = hfpso_v4(iter, swarm_size, c1, c2, LB, UB, nVars, vmax_coef, fitness_function);

% Hiển thị kết quả
finalPos = reshape(xc, [numel(xc) / 2, 2]);
T = table(finalPos(:,1), finalPos(:,2));
writetable(T, 'somefile.txt')
figure;
plot(finalPos(:, 1), finalPos(:, 2), '.', 'color', 'r');
hold on;
for ii = 1:N % vẽ phạm vi truyền hình tròn
    [finalcircle.x(ii, :), finalcircle.y(ii, :)] = circle(finalPos(ii, 1), finalPos(ii, 2), Trange);
    fill(finalcircle.x(ii, :), finalcircle.y(ii, :), [0.25, 0.25, 0.25]);
    alpha(0.3);
    hold on;
end
axis on;
xlabel('x(m)');
ylabel('y(m)');
title('Vị trí ban đầu tối ưu của các nút cảm biến với phạm vi đo hình tròn');
test_codeLASTcolap;
%% Các hàm phụ trợ

function [f_valc, xc] = hfpso_v4(iter, swarm_size, c1, c2, LB, UB, D, vmax_coef, fhd)
    rand('state', sum(100*clock));
    v_max = vmax_coef * (UB - LB);
    v_min = -v_max;
    
    particles_x = zeros(swarm_size, D);
    particles_v = zeros(swarm_size, D);
    
    for sw = 1:swarm_size
        for ds = 1:D
            particles_x(sw, ds) = LB(ds) + rand * (UB(ds) - LB(ds));
            particles_v(sw, ds) = v_min(ds) + rand * (v_max(ds) - v_min(ds));
        end
    end
    
    f_val = zeros(swarm_size, 1);
    for piiz = 1:swarm_size
        f_val(piiz, 1) = feval(fhd, particles_x(piiz, :));
    end
    
    p_best = particles_x;
    p_best_val = f_val;
    [~, index] = min(f_val(:, 1));
    g_best = particles_x(index, :);
    g_best_val = f_val(index, 1);
    dmax = (UB - LB) * sqrt(D);
    
    g_best_t = zeros(iter, D);
    g_best_val_t = zeros(iter, 1);
    variabless = zeros(iter, D);
    valuess = zeros(iter, 1);
    
    for i = 1:iter
        w_linear = 0.9 - ((0.9 - 0.5) / iter) * i; % Linear Decreasing Inertia Weight
        w = w_linear;
        
        for j = 1:swarm_size
            if (i > 2) && (f_val(j, 1) <= g_best_val_t(i-2, :))
                rij = norm(particles_x(j, :) - g_best_t(i-2, :)) ./ dmax;
                alpha = 0.2;
                beta0 = 2; 
                m = 2;
                gamma = 1;
                beta = beta0 * exp(-gamma * rij.^m);  
                e = rand(1, D) - 0.5;
                prev_pos = particles_x(j, :);
                particles_x(j, :) = particles_x(j, :) + beta .* (particles_x(j, :) - g_best_t(i-2, :)) + alpha .* e;
                
                for k = 1:D
                    if particles_x(j, k) > UB(k)
                        particles_x(j, k) = UB(k); 
                    end
                    if particles_x(j, k) < LB(k)
                        particles_x(j, k) = LB(k); 
                    end
                end
                
                particles_v(j, :) = particles_x(j, :) - prev_pos;
                
                for k = 1:D
                    if particles_v(j, k) > v_max(k)
                        particles_v(j, k) = v_max(k); 
                    end
                    if particles_v(j, k) < v_min(k)
                        particles_v(j, k) = v_min(k);
                    end 
                end
            else
                for k = 1:D
                    r1 = rand();
                    r2 = rand();
                    particles_v(j, k) = w * particles_v(j, k) + c1 * r1 * (p_best(j, k) - particles_x(j, k)) + c2 * r2 * (g_best(1, k) - particles_x(j, k)); 
                end
                
                for k = 1:D
                    if particles_v(j, k) > v_max(k)
                        particles_v(j, k) = v_max(k); 
                    end
                    if particles_v(j, k) < v_min(k)
                        particles_v(j, k) = v_min(k);
                    end
                end
                
                particles_x(j, :) = particles_x(j, :) + particles_v(j, :);
                
                for k = 1:D
                    if particles_x(j, k) > UB(k)
                        particles_x(j, k) = UB(k); 
                    end
                    if particles_x(j, k) < LB(k)
                        particles_x(j, k) = LB(k); 
                    end
                end
            end
        end
        
        for piiz = 1:swarm_size
            f_val(piiz, 1) = feval(fhd, particles_x(piiz, :));
        end
        
        for j = 1:swarm_size
            if f_val(j, 1) < p_best_val(j, 1)
                p_best(j, :) = particles_x(j, :);
                p_best_val(j, 1) = f_val(j, 1);
            end
            if p_best_val(j, 1) < g_best_val
                g_best = particles_x(j, :);
                g_best_val = p_best_val(j, 1);
            end
        end
        
        g_best_t(i, :) = g_best;
        g_best_val_t(i, :) = g_best_val;
        variabless(i, :) = g_best;
        valuess(i, 1) = g_best_val;
        disp(['Iteration ' num2str(i) ': Best Cost = ' num2str(valuess(i, 1))]);
        
    end
    figure;
    plot(1:iter, valuess, 'o', 'MarkerSize', 2, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b'); % 'o' để vẽ các điểm, 4 là kích thước của điểm
    xlabel('Iteration');
    ylabel('Best Cost');
    title('Best Cost vs. Iteration HFPSO ');
    xc = variabless(iter, :);
    f_valc = valuess(iter, :);
end



function [x, y] = circle(x_center, y_center, radius)
    theta = linspace(0, 2 * pi, 100);
    x = x_center + radius * cos(theta);
    y = y_center + radius * sin(theta);
end



