function fit=fitness_function_for_CH(SN, distances, indexArray, c, particle, centroid)

[val, ind]=min(distances(:,:,particle));
CH_array=indexArray(ind);
EnergySum=0;
BSDist=0;
for i=1:size(CH_array)
    EnergySum=EnergySum+SN(CH_array(i)).E;

end

fit=mean(distances(c(:,particle)==centroid,centroid,particle))/EnergySum;

% fit=mean(distances(c(:,particle)==centroid,centroid,particle));

% function fitness = fitness_function_for_CH(SN, x, meas_ids, centroids)
%     % Define the fitness function for the clustering problem here.
%     % For now, let's assume it calculates the sum of squared distances
%     % between the centroids and the nodes.
%     fitness = 0;
%     for i = 1:length(meas_ids)
%         cluster_id = meas_ids(i);
%         centroid = x((cluster_id-1)*2+1:cluster_id*2);
%         node = [SN(i).x, SN(i).y];
%         fitness = fitness + sum((node - centroid).^2);
%     end
% end