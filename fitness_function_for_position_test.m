function score = fitness_function_for_position_test(positions, Trange, area)
    Nodes_pos = reshape(positions, [numel(positions)/2, 2]);
    pts = 10000;
    pointspos = area(1) .* rand(pts, 2);
    coveredpt = zeros(pts, 1);
    for ii = 1:pts
        for jj = 1:size(Nodes_pos, 1)
            dist = sqrt((pointspos(ii, 1) - Nodes_pos(jj, 1))^2 + (pointspos(ii, 2) - Nodes_pos(jj, 2))^2);
            if dist < Trange || dist == Trange
                coveredpt(ii) = 1;
                break;
            end
        end
    end
    coverage = numel(find(coveredpt == 1)) / pts;
    score = 1 / coverage;
end
