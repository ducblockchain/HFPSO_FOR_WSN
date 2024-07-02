function CH_array = PSO_algo(SN, centroids)
    meas2 = [vertcat(SN(:).x), vertcat(SN(:).y)];
    dataset_size = size(meas2);
    meas = [];
    for i = 1:dataset_size
        if SN(i).E > 0
            meas = [meas; meas2(i, :), SN(i).id];
        end
    end
    dataset_size = size(meas);

    % INIT PARTICLE SWARM
    centroids = max(1, centroids);    % số cụm (centroids)
    dimensions = 2;                   % số chiều của mỗi centroid
    particles = 20;                   % số hạt trong swarm
    iterations = 50;                  % số lần lặp của thuật toán

    % EXECUTE K-MEANS
    hybrid_pso = true;                % kích hoạt hybrid PSO
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
    ranges = max(meas(:, [1, 2])) - min(meas(:, [1, 2])); %% scale
    swarm_pos = swarm_pos .* repmat(ranges, centroids, 1, particles) + repmat(min(meas(:, [1, 2])), centroids, 1, particles);
    swarm_fitness = Inf(particles, 1);
    if hybrid_pso
        swarm_pos(:,:,1) = KMEANS_CENTROIDS;
    end

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
        swarm_overall_pose = swarm_pos(:, :, index);  % GLOBAL BEST POSITION

        % SAMPLE r1 AND r2 FROM UNIFORM DISTRIBUTION [0..1]
        r1 = rand;
        r2 = rand;

        % UPDATE CLUSTER CENTROIDS
        for particle = 1:particles        
            inertia = w * swarm_vel(:, :, particle);
            cognitive = c1 * r1 * (swarm_best(:, :, particle) - swarm_pos(:, :, particle));
            social = c2 * r2 * (swarm_overall_pose - swarm_pos(:, :, particle));
            vel = inertia + cognitive + social;
            swarm_pos(:, :, particle) = swarm_pos(:, :, particle) + vel; % UPDATE PARTICLE POSE
            swarm_vel(:, :, particle) = vel;                            % UPDATE PARTICLE VEL
        end
    end

    % PLOT THE ASSOCIATIONS WITH RESPECT TO THE CLUSTER
    particle = index; 
    [~, ind] = min(distances(:, :, particle));
    CH_array = floor(meas(ind, 3));
end


