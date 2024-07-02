% function CH_array = HFPSO_algo(SN, centroids)
%     meas2 = [vertcat(SN(:).x), vertcat(SN(:).y)];
%     dataset_size = size(meas2);
%     meas = [];
%     for i = 1:dataset_size(1)
%         if SN(i).E > 0
%             meas = [meas; meas2(i, :), SN(i).id];
%         end
%     end
%     dataset_size = size(meas);
% 
%     % INIT PARTICLE SWARM
%     centroids = max(1, centroids);    % number of centroids
%     dimensions = 2;                   % number of dimensions
%     particles = 30;                   % number of particles in swarm
%     iterations = 300;                  % number of iterations
% 
%     % EXECUTE K-MEANS
%     hybrid_pso = true;                % activate hybrid PSO
%     if hybrid_pso
%         [idx, KMEANS_CENTROIDS] = kmeans(meas(:, [1, 2]), centroids, 'dist', 'sqEuclidean', 'display', 'off', 'start', 'uniform', 'onlinephase', 'off');
%     end
% 
%     % GLOBAL PARAMETERS (the paper reports this values 0.72;1.49;1.49)
%     w  = 0.72; % INERTIA
%     c1 = 1.49; % COGNITIVE
%     c2 = 1.49; % SOCIAL
% 
%     % SETTING UP PSO DATA STRUCTURES
%     swarm_vel = rand(centroids, dimensions, particles) * 0.1;
%     swarm_pos = rand(centroids, dimensions, particles);
%     ranges = max(meas(:, [1, 2])) - min(meas(:, [1, 2])); % scale
%     swarm_pos = swarm_pos .* repmat(ranges, centroids, 1, particles) + repmat(min(meas(:, [1, 2])), centroids, 1, particles);
%     swarm_fitness = Inf(particles, 1);
%     if hybrid_pso
%         swarm_pos(:,:,1) = KMEANS_CENTROIDS;
%     end
% 
%     % Initialize Firefly algorithm parameters
%     alpha = 0.2;        % Attraction coefficient base value
%     beta0 = 2;          % Attraction coefficient base value
%     gamma = 1;          % Light absorption coefficient
%     m = 2;              % Exponent in the attraction formula
% 
%     global_fitness_values = zeros(iterations, 1);
% 
%     for iteration = 1:iterations
%         % CALCULATE EUCLIDEAN DISTANCES TO ALL CENTROIDS
%         distances = zeros(dataset_size(1), centroids, particles);
%         for particle = 1:particles
%             for centroid = 1:centroids
%                 for data_vector = 1:dataset_size(1)
%                     distances(data_vector, centroid, particle) = norm(swarm_pos(centroid, :, particle) - meas(data_vector, [1, 2]));
%                 end
%             end
%         end
% 
%         % ASSIGN MEASURES with CLUSTERS    
%         c = zeros(dataset_size(1), particles);
%         for particle = 1:particles
%             [~, index] = min(distances(:, :, particle), [], 2);
%             c(:, particle) = index;
%         end
% 
%         % CALCULATE GLOBAL FITNESS and LOCAL FITNESS
%         average_fitness = zeros(particles, 1);
%         for particle = 1:particles
%             for centroid = 1:centroids
%                 if any(c(:, particle) == centroid)
%                     local_fitness = fitness_function_for_CH(SN, distances, floor(meas(:, 3)), c, particle, centroid);
%                     average_fitness(particle, 1) = average_fitness(particle, 1) + local_fitness;
%                 end
%             end
%             average_fitness(particle, 1) = average_fitness(particle, 1) / centroids;
%             if average_fitness(particle, 1) < swarm_fitness(particle)
%                 swarm_fitness(particle) = average_fitness(particle, 1);
%                 swarm_best(:, :, particle) = swarm_pos(:, :, particle); % LOCAL BEST FITNESS
%             end
%         end    
%         [global_fitness, index] = min(swarm_fitness); % GLOBAL BEST FITNESS
%         global_fitness_values(iteration) = global_fitness; % Store global fitness value
%         swarm_overall_pose = swarm_pos(:, :, index);  % GLOBAL BEST POSITION
% 
%         % SAMPLE r1 AND r2 FROM UNIFORM DISTRIBUTION [0..1]
%         r1 = rand;
%         r2 = rand;
% 
%         % UPDATE CLUSTER CENTROIDS USING HFPSO
%         for particle = 1:particles
%             % Apply PSO updates
%             inertia = w * swarm_vel(:, :, particle);
%             cognitive = c1 * r1 * (swarm_best(:, :, particle) - swarm_pos(:, :, particle));
%             social = c2 * r2 * (swarm_overall_pose - swarm_pos(:, :, particle));
%             vel = inertia + cognitive + social;
% 
%             % Apply Firefly updates if necessary
%             if iteration > 2 && average_fitness(particle) <= global_fitness
%                 rij = norm(swarm_pos(:, :, particle) - swarm_overall_pose) ./ max(ranges);
%                 beta = beta0 * exp(-gamma * rij.^m);  
%                 e = rand(1, dimensions) - 0.5;
%                 prev_pos = swarm_pos(:, :, particle);
%                 swarm_pos(:, :, particle) = swarm_pos(:, :, particle) + beta .* (swarm_pos(:, :, particle) - swarm_overall_pose) + alpha .* e;
% 
%                 % Ensure boundaries are respected
%                 for k = 1:dimensions
%                     if swarm_pos(:, k, particle) > max(meas(:, k))
%                         swarm_pos(:, k, particle) = max(meas(:, k));
%                     elseif swarm_pos(:, k, particle) < min(meas(:, k))
%                         swarm_pos(:, k, particle) = min(meas(:, k));
%                     end
%                 end
% 
%                 % Update velocity
%                 swarm_vel(:, :, particle) = swarm_pos(:, :, particle) - prev_pos;
%             else
%                 % Ensure boundaries are respected
%                 for k = 1:dimensions
%                     if swarm_pos(:, k, particle) > max(meas(:, k))
%                         swarm_pos(:, k, particle) = max(meas(:, k));
%                     elseif swarm_pos(:, k, particle) < min(meas(:, k))
%                         swarm_pos(:, k, particle) = min(meas(:, k));
%                     end
%                 end
% 
%                 % Update position and velocity
%                 swarm_pos(:, :, particle) = swarm_pos(:, :, particle) + vel;
%                 swarm_vel(:, :, particle) = vel;
%             end
%         end
%     end
% 
%     % PLOT THE ASSOCIATIONS WITH RESPECT TO THE CLUSTER
%     particle = index; 
%     [~, ind] = min(distances(:, :, particle));
%     CH_array = floor(meas(ind, 3));
%      % Plot the fitness values
%     figure;
%     plot(1:iteration, global_fitness_values, 'o', 'MarkerSize', 2, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b'); % 'o' để vẽ các điểm, 4 là kích thước của điểm
%     title('Global Fitness Values Over Iterations HFPSO Cluster head selection');
%     xlabel('Iteration');
%     ylabel('Fitness Value');
% end



function CH_array = HFPSO_algo(SN, centroids)
    meas2 = [vertcat(SN(:).x), vertcat(SN(:).y)];
    dataset_size = size(meas2);
    meas = [];
    for i = 1:dataset_size(1)
        if SN(i).E > 0
            meas = [meas; meas2(i, :), SN(i).id];
        end
    end
    dataset_size = size(meas);

    % INIT PARTICLE SWARM
    centroids = max(1, centroids);    % number of centroids
    dimensions = 2;                   % number of dimensions
    particles = 30;                   % number of particles in swarm
    iterations = 2000;                  % number of iterations

    % EXECUTE K-MEANS
    hybrid_pso = true;                % activate hybrid PSO
    if hybrid_pso
        [idx, KMEANS_CENTROIDS] = kmeans(meas(:, [1, 2]), centroids, 'dist', 'sqEuclidean', 'display', 'off', 'start', 'uniform', 'onlinephase', 'off');
    end

    % GLOBAL PARAMETERS (the paper reports this values 0.72;1.49;1.49)
    w  = 0.72; % INERTIA
    c1 = 1.49; % COGNITIVE
    c2 = 1.49; % SOCIAL

    % SETTING UP PSO DATA STRUCTURES
    swarm_vel = rand(centroids, dimensions, particles) * 0.1;
    swarm_pos = rand(centroids, dimensions, particles);
    ranges = max(meas(:, [1, 2])) - min(meas(:, [1, 2])); % scale
    swarm_pos = swarm_pos .* repmat(ranges, centroids, 1, particles) + repmat(min(meas(:, [1, 2])), centroids, 1, particles);
    swarm_fitness = Inf(particles, 1);
    if hybrid_pso
        swarm_pos(:,:,1) = KMEANS_CENTROIDS;
    end

    % Initialize Firefly algorithm parameters
    alpha = 0.2;        % Attraction coefficient base value
    beta0 = 2;          % Attraction coefficient base value
    gamma = 1;          % Light absorption coefficient
    m = 2;              % Exponent in the attraction formula

    global_fitness_values = zeros(iterations, 1);

    for iteration = 1:iterations
        % CALCULATE EUCLIDEAN DISTANCES TO ALL CENTROIDS
        distances = zeros(dataset_size(1), centroids, particles);
        for particle = 1:particles
            for centroid = 1:centroids
                for data_vector = 1:dataset_size(1)
                    distances(data_vector, centroid, particle) = norm(swarm_pos(centroid, :, particle) - meas(data_vector, [1, 2]));
                end
            end
        end

        % ASSIGN MEASURES with CLUSTERS    
        c = zeros(dataset_size(1), particles);
        for particle = 1:particles
            [~, index] = min(distances(:, :, particle), [], 2);
            c(:, particle) = index;
        end

        % CALCULATE GLOBAL FITNESS and LOCAL FITNESS
        average_fitness = zeros(particles, 1);
        for particle = 1:particles
            for centroid = 1:centroids
                if any(c(:, particle) == centroid)
                    local_fitness = fitness_function_for_CH(SN, distances, floor(meas(:, 3)), c, particle, centroid);
                    average_fitness(particle, 1) = average_fitness(particle, 1) + local_fitness;
                end
            end
            average_fitness(particle, 1) = average_fitness(particle, 1) / centroids;
            if average_fitness(particle, 1) < swarm_fitness(particle)
                swarm_fitness(particle) = average_fitness(particle, 1);
                swarm_best(:, :, particle) = swarm_pos(:, :, particle); % LOCAL BEST FITNESS
            end
        end    
        [global_fitness, index] = min(swarm_fitness); % GLOBAL BEST FITNESS
        global_fitness_values(iteration) = global_fitness; % Store global fitness value
        swarm_overall_pose = swarm_pos(:, :, index);  % GLOBAL BEST POSITION

        % SAMPLE r1 AND r2 FROM UNIFORM DISTRIBUTION [0..1]
        r1 = rand;
        r2 = rand;

        % UPDATE CLUSTER CENTROIDS USING HFPSO
        for particle = 1:particles
            % Apply PSO updates
            inertia = w * swarm_vel(:, :, particle);
            cognitive = c1 * r1 * (swarm_best(:, :, particle) - swarm_pos(:, :, particle));
            social = c2 * r2 * (swarm_overall_pose - swarm_pos(:, :, particle));
            vel = inertia + cognitive + social;

            % Apply Firefly updates if necessary
            if iteration > 2 && average_fitness(particle) <= global_fitness
                rij = norm(swarm_pos(:, :, particle) - swarm_overall_pose) ./ max(ranges);
                beta = beta0 * exp(-gamma * rij.^m);  
                e = rand(1, dimensions) - 0.5;
                prev_pos = swarm_pos(:, :, particle);
                swarm_pos(:, :, particle) = swarm_pos(:, :, particle) + beta .* (swarm_pos(:, :, particle) - swarm_overall_pose) + alpha .* e;

                % Ensure boundaries are respected
                for k = 1:dimensions
                    if swarm_pos(:, k, particle) > max(meas(:, k))
                        swarm_pos(:, k, particle) = max(meas(:, k));
                    elseif swarm_pos(:, k, particle) < min(meas(:, k))
                        swarm_pos(:, k, particle) = min(meas(:, k));
                    end
                end

                % Update velocity
                swarm_vel(:, :, particle) = swarm_pos(:, :, particle) - prev_pos;
            else
                % Ensure boundaries are respected
                for k = 1:dimensions
                    if swarm_pos(:, k, particle) > max(meas(:, k))
                        swarm_pos(:, k, particle) = max(meas(:, k));
                    elseif swarm_pos(:, k, particle) < min(meas(:, k))
                        swarm_pos(:, k, particle) = min(meas(:, k));
                    end
                end

                % Update position and velocity
                swarm_pos(:, :, particle) = swarm_pos(:, :, particle) + vel;
                swarm_vel(:, :, particle) = vel;
            end
        end
    end

    % PLOT THE ASSOCIATIONS WITH RESPECT TO THE CLUSTER
    particle = index; 
    [~, ind] = min(distances(:, :, particle));
    CH_array = floor(meas(ind, 3));
     % Plot the fitness values
    % figure;
    % plot(1:iteration, global_fitness_values, 'o', 'MarkerSize', 2, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');
end




% function CH_array = HFPSO_algo(SN, centroids)
%     meas2 = [vertcat(SN(:).x), vertcat(SN(:).y)];
%     dataset_size = size(meas2);
%     meas = [];
%     for i = 1:dataset_size(1)
%         if SN(i).E > 0
%             meas = [meas; meas2(i, :), SN(i).id];
%         end
%     end
%     dataset_size = size(meas);
% 
%     % INIT PARTICLE SWARM
%     centroids = max(1, 5);    % number of centroids
%     dimensions = 2;                   % number of dimensions
%     particles = 30;                   % number of particles in swarm
%     iterations = 300;                 % number of iterations
% 
%     % EXECUTE K-MEANS
%     hybrid_pso = true;                % activate hybrid PSO
%     if hybrid_pso
%         [idx, KMEANS_CENTROIDS] = kmeans(meas(:, [1, 2]), centroids, 'dist', 'sqEuclidean', 'display', 'off', 'start', 'uniform', 'onlinephase', 'off');
%     end
% 
%     % GLOBAL PARAMETERS
%     w_max = 0.9;
%     w_min = 0.5;
%     c1 = 1.49445; % COGNITIVE
%     c2 = 1.49445; % SOCIAL
% 
%     % SETTING UP PSO DATA STRUCTURES
%     swarm_vel = rand(centroids, dimensions, particles) * 0.1;
%     swarm_pos = rand(centroids, dimensions, particles);
%     ranges = max(meas(:, [1, 2])) - min(meas(:, [1, 2])); % scale
%     swarm_pos = swarm_pos .* repmat(ranges, centroids, 1, particles) + repmat(min(meas(:, [1, 2])), centroids, 1, particles);
%     swarm_fitness = Inf(particles, 1);
%     if hybrid_pso
%         swarm_pos(:,:,1) = KMEANS_CENTROIDS;
%     end
% 
%     % Initialize Firefly algorithm parameters
%     alpha = 0.2;        % Attraction coefficient base value
%     beta0 = 2;          % Attraction coefficient base value
%     gamma = 1;          % Light absorption coefficient
%     m = 2;              % Exponent in the attraction formula
% 
%     global_fitness_values = zeros(iterations, 1);
% 
%     for iteration = 1:iterations
%         % Calculate inertia weight linearly decreasing
%         w = w_max - ((w_max - w_min) / iterations) * iteration;
% 
%         % CALCULATE EUCLIDEAN DISTANCES TO ALL CENTROIDS
%         distances = zeros(dataset_size(1), centroids, particles);
%         for particle = 1:particles
%             for centroid = 1:centroids
%                 for data_vector = 1:dataset_size(1)
%                     distances(data_vector, centroid, particle) = norm(swarm_pos(centroid, :, particle) - meas(data_vector, [1, 2]));
%                 end
%             end
%         end
% 
%         % ASSIGN MEASURES with CLUSTERS    
%         c = zeros(dataset_size(1), particles);
%         for particle = 1:particles
%             [~, index] = min(distances(:, :, particle), [], 2);
%             c(:, particle) = index;
%         end
% 
%         % CALCULATE GLOBAL FITNESS and LOCAL FITNESS
%         average_fitness = zeros(particles, 1);
%         for particle = 1:particles
%             for centroid = 1:centroids
%                 if any(c(:, particle) == centroid)
%                     local_fitness = fitness_function_for_CH(SN, distances, floor(meas(:, 3)), c, particle, centroid);
%                     average_fitness(particle, 1) = average_fitness(particle, 1) + local_fitness;
%                 end
%             end
%             average_fitness(particle, 1) = average_fitness(particle, 1) / centroids;
%             if average_fitness(particle, 1) < swarm_fitness(particle)
%                 swarm_fitness(particle) = average_fitness(particle, 1);
%                 swarm_best(:, :, particle) = swarm_pos(:, :, particle); % LOCAL BEST FITNESS
%             end
%         end    
%         [global_fitness, index] = min(swarm_fitness); % GLOBAL BEST FITNESS
%         global_fitness_values(iteration) = global_fitness; % Store global fitness value
%         swarm_overall_pose = swarm_pos(:, :, index);  % GLOBAL BEST POSITION
% 
%         % UPDATE CLUSTER CENTROIDS USING HFPSO
%         for particle = 1:particles
%             % Apply PSO updates
%             r1 = rand;
%             r2 = rand;
%             inertia = w * swarm_vel(:, :, particle);
%             cognitive = c1 * r1 * (swarm_best(:, :, particle) - swarm_pos(:, :, particle));
%             social = c2 * r2 * (swarm_overall_pose - swarm_pos(:, :, particle));
%             vel = inertia + cognitive + social;
% 
%             % Apply Firefly updates if necessary
%             if iteration > 2 && average_fitness(particle) <= global_fitness
%                 rij = norm(swarm_pos(:, :, particle) - swarm_overall_pose) ./ max(ranges);
%                 beta = beta0 * exp(-gamma * rij.^m);  
%                 e = rand(1, dimensions) - 0.5;
%                 prev_pos = swarm_pos(:, :, particle);
%                 swarm_pos(:, :, particle) = swarm_pos(:, :, particle) + beta .* (swarm_pos(:, :, particle) - swarm_overall_pose) + alpha .* e;
% 
%                 % Ensure boundaries are respected
%                 for k = 1:dimensions
%                     if swarm_pos(:, k, particle) > max(meas(:, k))
%                         swarm_pos(:, k, particle) = max(meas(:, k));
%                     elseif swarm_pos(:, k, particle) < min(meas(:, k))
%                         swarm_pos(:, k, particle) = min(meas(:, k));
%                     end
%                 end
% 
%                 % Update velocity
%                 swarm_vel(:, :, particle) = swarm_pos(:, :, particle) - prev_pos;
%             else
%                 % Ensure boundaries are respected
%                 for k = 1:dimensions
%                     if swarm_pos(:, k, particle) > max(meas(:, k))
%                         swarm_pos(:, k, particle) = max(meas(:, k));
%                     elseif swarm_pos(:, k, particle) < min(meas(:, k))
%                         swarm_pos(:, k, particle) = min(meas(:, k));
%                     end
%                 end
% 
%                 % Update position and velocity
%                 swarm_pos(:, :, particle) = swarm_pos(:, :, particle) + vel;
%                 swarm_vel(:, :, particle) = vel;
%             end
%         end
%     end
% 
%     % PLOT THE ASSOCIATIONS WITH RESPECT TO THE CLUSTER
%     particle = index; 
%     [~, ind] = min(distances(:, :, particle));
%     CH_array = floor(meas(ind, 3));
% 
%     % Plot the fitness values
%     figure;
%     plot(1:iterations, global_fitness_values, 'o', 'MarkerSize', 2, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b');
%     xlabel('Iteration');
%     ylabel('Global Fitness');
%     title('Global Fitness vs. Iteration for HFPSO');
% end


% function CH_array = HFPSO_algo(SN, centroids)
%     meas2 = [vertcat(SN(:).x), vertcat(SN(:).y)];
%     dataset_size = size(meas2);
%     meas = [];
%     for i = 1:dataset_size(1)
%         if SN(i).E > 0
%             meas = [meas; meas2(i, :), SN(i).id];
%         end
%     end
%     dataset_size = size(meas);
% 
%     % INIT PARTICLE SWARM
%     centroids = max(1, centroids);    % number of centroids
%     dimensions = 2;                   % number of dimensions
%     particles = 30;                   % number of particles in swarm
%     iterations = 300;                 % number of iterations
% 
%     % EXECUTE K-MEANS
%     hybrid_pso = true;                % activate hybrid PSO
%     if hybrid_pso
%         [idx, KMEANS_CENTROIDS] = kmeans(meas(:, [1, 2]), centroids, 'dist', 'sqEuclidean', 'display', 'off', 'start', 'uniform', 'onlinephase', 'off');
%     end
% 
%     % Global PSO parameters
%     w  = 0.72; % Inertia
%     c1 = 1.49; % Cognitive
%     c2 = 1.49; % Social
% 
%     % Setting up PSO data structures
%     ranges = max(meas(:, [1, 2])) - min(meas(:, [1, 2])); % scale
%     LB = repmat(min(meas(:, [1, 2])), 1, centroids);
%     UB = repmat(max(meas(:, [1, 2])), 1, centroids);
% 
%     % Define fitness function for the clustering problem
%     fitness_function = @(x) fitness_function_for_CH(SN, x, floor(meas(:, 3)), 5);
% 
%     % Execute HFPSO
%     [global_fitness, global_best_position] = hfpso_v4(iterations, particles, c1, c2, LB, UB, centroids * dimensions, 0.2, fitness_function);
% 
%     % Plot the global fitness values over iterations
%     figure;
%     plot(1:iterations, global_fitness, 'o', 'MarkerSize', 2, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b'); % 'o' to plot points, 2 is the size of the points
%     title('Global Fitness Values Over Iterations HFPSO Cluster head selection');
%     xlabel('Iteration');
%     ylabel('Fitness Value');
% 
%     % Assign cluster heads based on the global best position
%     [~, ind] = min(pdist2(global_best_position, meas(:, [1, 2])));
%     CH_array = floor(meas(ind, 3));
% end
% 
% 
% function [f_valc, xc] = hfpso_v4(iter, swarm_size, c1, c2, LB, UB, D, vmax_coef, fhd)
%     rand('state', sum(100*clock));
%     v_max = vmax_coef * (UB - LB);
%     v_min = -v_max;
% 
%     particles_x = zeros(swarm_size, D);
%     particles_v = zeros(swarm_size, D);
% 
%     for sw = 1:swarm_size
%         for ds = 1:D
%             particles_x(sw, ds) = LB(ds) + rand * (UB(ds) - LB(ds));
%             particles_v(sw, ds) = v_min(ds) + rand * (v_max(ds) - v_min(ds));
%         end
%     end
% 
%     f_val = zeros(swarm_size, 1);
%     for piiz = 1:swarm_size
%         f_val(piiz, 1) = feval(fhd, particles_x(piiz, :));
%     end
% 
%     p_best = particles_x;
%     p_best_val = f_val;
%     [~, index] = min(f_val(:, 1));
%     g_best = particles_x(index, :);
%     g_best_val = f_val(index, 1);
%     dmax = (UB - LB) * sqrt(D);
% 
%     g_best_t = zeros(iter, D);
%     g_best_val_t = zeros(iter, 1);
%     variabless = zeros(iter, D);
%     valuess = zeros(iter, 1);
% 
%     for i = 1:iter
%         w_linear = 0.9 - ((0.9 - 0.5) / iter) * i; % Linear Decreasing Inertia Weight
%         w = w_linear;
% 
%         for j = 1:swarm_size
%             if (i > 2) && (f_val(j, 1) <= g_best_val_t(i-2, :))
%                 rij = norm(particles_x(j, :) - g_best_t(i-2, :)) ./ dmax;
%                 alpha = 0.2;
%                 beta0 = 2; 
%                 m = 2;
%                 gamma = 1;
%                 beta = beta0 * exp(-gamma * rij.^m);  
%                 e = rand(1, D) - 0.5;
%                 prev_pos = particles_x(j, :);
%                 particles_x(j, :) = particles_x(j, :) + beta .* (particles_x(j, :) - g_best_t(i-2, :)) + alpha .* e;
% 
%                 for k = 1:D
%                     if particles_x(j, k) > UB(k)
%                         particles_x(j, k) = UB(k); 
%                     end
%                     if particles_x(j, k) < LB(k)
%                         particles_x(j, k) = LB(k); 
%                     end
%                 end
% 
%                 particles_v(j, :) = particles_x(j, :) - prev_pos;
% 
%                 for k = 1:D
%                     if particles_v(j, k) > v_max(k)
%                         particles_v(j, k) = v_max(k); 
%                     end
%                     if particles_v(j, k) < v_min(k)
%                         particles_v(j, k) = v_min(k);
%                     end 
%                 end
%             else
%                 for k = 1:D
%                     r1 = rand();
%                     r2 = rand();
%                     particles_v(j, k) = w * particles_v(j, k) + c1 * r1 * (p_best(j, k) - particles_x(j, k)) + c2 * r2 * (g_best(1, k) - particles_x(j, k)); 
%                 end
% 
%                 for k = 1:D
%                     if particles_v(j, k) > v_max(k)
%                         particles_v(j, k) = v_max(k); 
%                     end
%                     if particles_v(j, k) < v_min(k)
%                         particles_v(j, k) = v_min(k);
%                     end
%                 end
% 
%                 particles_x(j, :) = particles_x(j, :) + particles_v(j, :);
% 
%                 for k = 1:D
%                     if particles_x(j, k) > UB(k)
%                         particles_x(j, k) = UB(k); 
%                     end
%                     if particles_x(j, k) < LB(k)
%                         particles_x(j, k) = LB(k); 
%                     end
%                 end
%             end
%         end
% 
%         for piiz = 1:swarm_size
%             f_val(piiz, 1) = feval(fhd, particles_x(piiz, :));
%         end
% 
%         for j = 1:swarm_size
%             if f_val(j, 1) < p_best_val(j, 1)
%                 p_best(j, :) = particles_x(j, :);
%                 p_best_val(j, 1) = f_val(j, 1);
%             end
%             if p_best_val(j, 1) < g_best_val
%                 g_best = particles_x(j, :);
%                 g_best_val = p_best_val(j, 1);
%             end
%         end
% 
%         g_best_t(i, :) = g_best;
%         g_best_val_t(i, :) = g_best_val;
%         variabless(i, :) = g_best;
%         valuess(i, 1) = g_best_val;
%         disp(['Iteration ' num2str(i) ': Best Cost = ' num2str(valuess(i, 1))]);
% 
%     end
%     figure;
%     plot(1:iter, valuess, 'o', 'MarkerSize', 2, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b'); % 'o' để vẽ các điểm, 4 là kích thước của điểm
%     xlabel('Iteration');
%     ylabel('Best Cost');
%     title('Best Cost vs. Iteration HFPSO ');
%     xc = variabless(iter, :);
%     f_valc = valuess(iter, :);
% end