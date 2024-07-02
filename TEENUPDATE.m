close all;
clear;
clc;

%%%%%%%%%%%%%%%%%%%% Network Establishment Parameters %%%%%%%%%%%%%%%%%%%%
%%% Area of Operation %%%
% Field Dimensions in meters %
xm = 120;
ym = 120;
x = 0; % added for better display results of the plot
y = 0; % added for better display results of the plot
% Number of Nodes in the field %
n = 100;
% Number of Dead Nodes in the beginning %
dead_nodes = 0;
% Coordinates of the Sink (location is predetermined in this simulation) %
sinkx = 300;
sinky = 300;
%%% Energy Values %%%
% Initial Energy of a Node (in Joules) % 
Eo = 1; % units in Joules
% Energy required to run circuitry (both for transmitter and receiver) %
Eelec = 50 * 10^(-9); % units in Joules/bit
ETx = 50 * 10^(-9); % units in Joules/bit
ERx = 50 * 10^(-9); % units in Joules/bit
% Transmit Amplifier Types %
Eamp = 100 * 10^(-12); % units in Joules/bit/m^2 (amount of energy spent by the amplifier to transmit the bits)
% Data Aggregation Energy %
EDA = 5 * 10^(-9); % units in Joules/bit
% Size of data package %
k = 4000; % units in bits
% Suggested percentage of cluster head %
p = 0.05; % a 5 percent of the total amount of nodes used in the network is proposed to give good results
% Number of Clusters %
No = p * n;
% Round of Operation %
rnd = 0;
% Current Number of operating Nodes %
operating_nodes = n;
transmissions = 0;
temp_val = 0;
flag1stdead = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%% End of Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Creation of the Wireless Sensor Network %%%
% Plotting the WSN %
for i = 1:n
    SN(i).id = i; % sensor's ID number
    SN(i).x = rand(1, 1) * xm; % X-axis coordinates of sensor node
    SN(i).y = rand(1, 1) * ym; % Y-axis coordinates of sensor node
    SN(i).E = Eo; % nodes energy levels (initially set to be equal to "Eo"
    SN(i).role = 0; % node acts as normal if the value is '0', if elected as a cluster head it  gets the value '1' (initially all nodes are normal)
    SN(i).cluster = 0; % the cluster which a node belongs to
    SN(i).cond = 1; % States the current condition of the node. when the node is operational its value is =1 and when dead =0
    SN(i).rop = 0; % number of rounds node was operational
    SN(i).rleft = 0; % rounds left for node to become available for Cluster Head election
    SN(i).dtch = 0; % nodes distance from the cluster head of the cluster in which he belongs
    SN(i).dts = 0; % nodes distance from the sink
    SN(i).tel = 0; % states how many times the node was elected as a Cluster Head
    SN(i).rn = 0; % round node got elected as cluster head
    SN(i).chid = 0; % node ID of the cluster head which the "i" normal node belongs to
    
    % Initialize last_value and value for TEEN protocol
    SN(i).last_value = 0; % Initial last sensed value
    SN(i).value = rand() * 10; % Random initial sensed value
    
    hold on;
    figure(1)
    plot(x, y, xm, ym, SN(i).x, SN(i).y, 'ob', sinkx, sinky, '*r');
    title 'Wireless Sensor Network';
    xlabel '(m)';
    ylabel '(m)';
end

% Setting Thresholds for TEEN Protocol
HT = 5; % Hard Threshold
ST = 2; % Soft Threshold

                      %%%%%% Set-Up Phase %%%%%% 
first_dead = -1;
middle_dead = -1;
last_dead = -1;
remaining_energy=[];                      
total_energy = [];  
energy_per_round=[];
while operating_nodes > 0
    % Displays Current Round %     
    rnd     
    % Threshold Value %
    t = (p / (1 - p * (mod(rnd, 1 / p))));
    
    % Re-election Value %
    tleft = mod(rnd, 1 / p);
 
    % Resetting Previous Amount Of Cluster Heads In the Network %
    CLheads = 0;
    
    % Resetting Previous Amount Of Energy Consumed In the Network on the Previous Round %
    energy = 0;
 
    % Cluster Heads Election %
    for i = 1:n
        SN(i).cluster = 0; % resetting cluster in which the node belongs to
        SN(i).role = 0; % resetting node role
        SN(i).chid = 0; % resetting cluster head id
        if SN(i).rleft > 0
            SN(i).rleft = SN(i).rleft - 1;
        end
        if (SN(i).E > 0) && (SN(i).rleft == 0)
            generate = rand;    
            if generate < t
                SN(i).role = 1; % assigns the node role of a cluster head
                SN(i).rn = rnd; % Assigns the round that the cluster head was elected to the data table
                SN(i).tel = SN(i).tel + 1;   
                SN(i).rleft = 1 / p - tleft; % rounds for which the node will be unable to become a CH
                SN(i).dts = sqrt((sinkx - SN(i).x)^2 + (sinky - SN(i).y)^2); % calculates the distance between the sink and the cluster head
                CLheads = CLheads + 1; % sum of cluster heads that have been elected 
                SN(i).cluster = CLheads; % cluster of which the node got elected to be cluster head
                CL(CLheads).x = SN(i).x; % X-axis coordinates of elected cluster head
                CL(CLheads).y = SN(i).y; % Y-axis coordinates of elected cluster head
                CL(CLheads).id = i; % Assigns the node ID of the newly elected cluster head to an array
            end
        end
    end
        
    % Fixing the size of "CL" array %
    CL = CL(1:CLheads);
  
    % Grouping the Nodes into Clusters & calculating the distance between node and cluster head %
    for i = 1:n
        if (SN(i).role == 0) && (SN(i).E > 0) && (CLheads > 0) % if node is normal
            for m = 1:CLheads
                d(m) = sqrt((CL(m).x - SN(i).x)^2 + (CL(m).y - SN(i).y)^2);
            end
            d = d(1:CLheads); % fixing the size of "d" array
            [M, I] = min(d(:)); % finds the minimum distance of node to CH
            [Row, Col] = ind2sub(size(d), I); % displays the Cluster Number in which this node belongs too
            SN(i).cluster = Col; % assigns node to the cluster
            SN(i).dtch = d(Col); % assigns the distance of node to CH
            SN(i).chid = CL(Col).id;
        end
    end
        
                       %%%%%% Steady-State Phase %%%%%%
    % Energy Dissipation for normal nodes %
    for i = 1:n
        if (SN(i).cond == 1) && (SN(i).role == 0) && (CLheads > 0)
            if SN(i).E > 0
                % Check thresholds for TEEN
                if (abs(SN(i).last_value - SN(i).value) >= ST) && (SN(i).value > HT)
                    ETx = Eelec * k + Eamp * k * SN(i).dtch^2;
                    SN(i).E = SN(i).E - ETx;
                    energy = energy + ETx;
                    
                    % Dissipation for cluster head during reception
                    if SN(SN(i).chid).E > 0 && SN(SN(i).chid).cond == 1 && SN(SN(i).chid).role == 1
                        ERx = (Eelec + EDA) * k;
                        SN(SN(i).chid).E = SN(SN(i).chid).E - ERx;
                        energy = energy + ERx;
                        transmissions = transmissions + 1;
                        % Cluster heads energy depletion
                        if SN(SN(i).chid).E <= 0
                            SN(SN(i).chid).cond = 0;
                            SN(SN(i).chid).rop = rnd;
                            dead_nodes = dead_nodes + 1;
                            operating_nodes = operating_nodes - 1;
                        end
                    end
                    % Update last sensed value
                    SN(i).last_value = SN(i).value;
                end
            end
            if SN(i).E <= 0
                dead_nodes = dead_nodes + 1;
                operating_nodes = operating_nodes - 1;
                SN(i).cond = 0;
                SN(i).chid = 0;
                SN(i).rop = rnd;
            end
        end
    end

    % Energy Dissipation for cluster head nodes %
    for i = 1:n
        if (SN(i).cond == 1) && (SN(i).role == 1)
            if SN(i).E > 0
                ETx = (Eelec + EDA) * k + Eamp * k * SN(i).dts^2;
                SN(i).E = SN(i).E - ETx;
                energy = energy + ETx;
            end
            if SN(i).E <= 0
                dead_nodes = dead_nodes + 1;
                operating_nodes = operating_nodes - 1;
                SN(i).cond = 0;
                SN(i).rop = rnd;
            end
        end
    end
   
    if operating_nodes < n && temp_val == 0
        temp_val = 1;
        flag1stdead = rnd;
    end
    if rnd>=1
      remaining_energy(rnd) = sum([SN([SN.cond] == 1).E]);
    end
% Display Number of Cluster Heads of this round %
%CLheads;
    if rnd == 1
         total_energy(rnd) = energy;
    elseif rnd >= 2
         total_energy(rnd) = total_energy(rnd - 1) + energy; % accumulate energy
    end
    if rnd == 1
       energy_per_round(rnd) = energy;
    elseif rnd >= 2
       energy_per_round(rnd) = energy; % accumulate energy
    end
    % Lưu trữ vòng chết của nút đầu tiên, giữa và cuối cùng
    if first_dead == -1 && dead_nodes > 0
        first_dead = rnd;
    end
    if middle_dead == -1 && dead_nodes >= 50
        middle_dead = rnd;
    end
    if last_dead == -1 && dead_nodes == 100
        last_dead = rnd;
    end
    % Next Round %
    rnd = rnd + 1;

    tr(transmissions) = operating_nodes;
    op(rnd) = operating_nodes;

    if energy > 0
        nrg(transmissions) = energy;
    end
end

sum = 0;
for i = 1:flag1stdead
    if i <= length(nrg)
        sum = nrg(i) + sum;
    end
end
temp1 = sum / flag1stdead;
temp2 = temp1 / n;
for i = 1:flag1stdead
    avg_node(i) = temp2;
end
    
% Plotting Simulation Results "Operating Nodes per Round" %
figure(2)
plot(1:rnd, op(1:rnd), '-r', 'LineWidth', 2);

% Transpose the data to write in vertical format
T = table((1:rnd)', op(1:rnd)');
writetable(T, 'teen.txt')

title({'TEEN'; 'Operating Nodes per Round'});
xlabel('Rounds');
ylabel('Operational Nodes');
hold on;

% Plotting Simulation Results  %
figure(3)
plot(1:transmissions, tr(1:transmissions), '-r', 'Linewidth', 2);
title({'TEEN'; 'Operational Nodes per Transmission';})
xlabel 'Transmissions';
ylabel 'Operational Nodes';
hold on;

% Plotting Simulation Results  %
figure(4)
%plot(1:flag1stdead, nrg(1:flag1stdead), '-r', 'Linewidth', 2);
title({'TEEN'; 'Energy consumed per Transmission';})
xlabel 'Transmission';
ylabel 'Energy ( J )';
hold on;

% Plotting Simulation Results  %
figure(5)
plot(1:flag1stdead, avg_node(1:flag1stdead), '-r', 'Linewidth', 2);
title({'TEEN'; 'Average Energy consumed by a Node per Transmission';});
xlabel 'Transmissions';
ylabel 'Energy ( J )';
hold on;
figure(6)
plot(1:rnd-1, total_energy(1:rnd-1), '-r', 'LineWidth', 2);
A = table((1:rnd-1)', total_energy(1:rnd-1)');
writetable(A, 'Teen_total_energy.txt')
title({'TEEN'; 'Total Energy used by all Active Nodes per Round'});
xlabel('Rounds');
ylabel('Total Energy ( J )');
hold on;

% Plotting total remaining energy of all active nodes per round
figure(7) 
plot(1:rnd-1, remaining_energy(1:rnd-1), '-b', 'LineWidth', 2);
B = table((1:rnd-1)', remaining_energy(1:rnd-1)');
writetable(B, 'Teen_remain_energy.txt')
title({'TEEN'; 'Total Remaining Energy of All Active Nodes per Round'});
xlabel('Rounds');
ylabel('Total Remaining Energy ( J )');
hold on;
figure(8)
plot(1:rnd-1,energy_per_round(1:rnd-1),'-b','Linewidth',2);
C = table((1:rnd-1)', energy_per_round(1:rnd-1)');
writetable(C, 'Teen_energy_per_round.txt')
title ({'PEGASIS'; ' Energy per round All Nodes per Transmission';})
xlabel 'Transmissions';
ylabel 'Total Energy ( J )';
hold on;
disp(['Vòng mà nút đầu tiên chết: ', num2str(first_dead)]);
disp(['Vòng mà nút giữa chết: ', num2str(middle_dead)]);
disp(['Vòng mà nút cuối cùng chết: ', num2str(last_dead)]);
