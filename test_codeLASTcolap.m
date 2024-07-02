clear;


% Field Dimensions in meters %
xm=120;
ym=120;

x=0; % added for better display results of the plot
y=0; % added for better display results of the plot

% Number of Nodes in the field %
n=100;

% Number of Dead Nodes in the beggining %
dead_nodes=0;

% Coordinates of the Sink (location is predetermined in this simulation) %
sinkx=300;
sinky=300;

%%% Energy Values %%%
% Initial Energy of a Node (in Joules) % 
Eo=1; % units in Joules

% Energy required to run circuity (both for transmitter and receiver) %
Eelec=50*10^(-9); % units in Joules/bit
ETx=50*10^(-9); % units in Joules/bit
ERx=50*10^(-9); % units in Joules/bit

% Transmit Amplifier Types %
Eamp=100*10^(-12); % units in Joules/bit/m^2 (amount of energy spent by the amplifier to transmit the bits)

% Data Aggregation Energy %
EDA=5*10^(-9); % units in Joules/bit

% Size of data package %
k=4000;

% Suggested percentage of cluster head %
p=0.05; % a 5 percent of the total amount of nodes used in the network is proposed to give good results
% Number of Clusters %
No=p*n; 

% Round of Operation %
rnd=0;

% Current Number of operating Nodes %
operating_nodes=n;
transmissions=0;
temp_val=0;
flag1stdead=0;
totalenergy=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%% End of Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid = fopen('somefile.txt', 'r') ;              
fgetl(fid) ;                                  
                                   
buffer = fread(fid, Inf) ;                    
fclose(fid)
fid = fopen('somefile.txt', 'w')  ;   
fwrite(fid, buffer) ;                         
fclose(fid) ;

%% Creation of the optimized Sensor Network
A = load('somefile.txt');
%%% Creation of the Wireless Sensor Network %%%
% Plotting the WSN %
for i=1:n
    S(i).xd=A(i,1);
    SN(i).x=S(i).xd;
    S(i).yd=A(i,2);
    SN(i).y=S(i).yd;
    SN(i).id=i;	% sensor's ID number
    %SN(i).E=Eo;     % nodes energy levels (initially set to be equal to "Eo"
     
    SN(i).E=Eo; %generate random numbers between 2 and 5;
    %To generate random  numbers b/w a and b r = a + (b-a)*rand();

    SN(i).role=0;   % node acts as normal if the value is '0', if elected as a cluster head it  gets the value '1' (initially all nodes are normal)
    SN(i).cluster=0;	% the cluster which a node belongs to
    SN(i).cond=1;	% States the current condition of the node. when the node is operational its value is =1 and when dead =0
    SN(i).rop=0;	% number of rounds node was operational
    SN(i).dtch=0;	% nodes distance from the cluster head of the cluster in which he belongs
    SN(i).dts=0;    % nodes distance from the sink
    SN(i).tel=0;	% states how many times the node was elected as a Cluster Head
    SN(i).rn=0;     % round node got elected as cluster head
    SN(i).chid=0;   % node ID of the cluster head which the "i" normal node belongs to
    
    hold on;
    figure(9)
    plot(x,y,xm,ym,SN(i).x,SN(i).y,'ob',sinkx,sinky,'*r');
    title 'Wireless Sensor Network';
    xlabel '(m)';
    ylabel '(m)';
    
end

meas = [vertcat(SN(:).x), vertcat(SN(:).y)];
dataset_size = size (meas);
CH_array=zeros(1, No);
first_dead = -1;
middle_dead = -1;
last_dead = -1;
                      %%%%%% Set-Up Phase %%%%%% 
remaining_energy=[];                      
total_energy = []; 
energy_per_round=[];

while operating_nodes>0
    % Displays Current Round %     
    rnd     
	% Reseting Previous Amount Of Cluster Heads In the Network %
	CLheads=0;
    
    % Reseting Previous Amount Of Energy Consumed In the Network on the Previous Round %
    energy=0;

    No=floor(p*operating_nodes);
 
    CHdead=0;

% Cluster Heads Election %
    if(mod(rnd,10)==0 || CHdead==1)
        CH_array=HFPSO_algo(SN,No);
        CHdead=0;
    end

        for i=1:n
            SN(i).cluster=0;    % reseting cluster in which the node belongs to
            SN(i).role=0;       % reseting node role
            SN(i).chid=0;       % reseting cluster head id
            
            %RUN PSO algorithm%
            if(any(CH_array(:)==(SN(i).id) ))    
            
                SN(i).role=1;	% assigns the node role of acluster head
                SN(i).rn=rnd;	% Assigns the round that the cluster head was elected to the data table
                SN(i).tel=SN(i).tel + 1;
                SN(i).dts=sqrt((sinkx-SN(i).x)^2 + (sinky-SN(i).y)^2); % calculates the distance between the sink and the cluster head
                CLheads=CLheads+1;	% sum of cluster heads that have been elected 
                SN(i).cluster=CLheads; % cluster of which the node got elected to be cluster head
                CL(CLheads).x=SN(i).x; % X-axis coordinates of elected cluster head
                CL(CLheads).y=SN(i).y; % Y-axis coordinates of elected cluster head
                CL(CLheads).id=i; % Assigns the node ID of the newly elected cluster head to an array
     
            end
            

        end
    
	% Fixing the size of "CL" array %
	CL=CL(1:CLheads);
  
    
    
    
% Grouping the Nodes into Clusters & calculating the distance between node and cluster head %
     
       for i=1:n
        if  (SN(i).role==0) && (SN(i).E>0) && (CLheads>0) % if node is normal
            for m=1:CLheads
            d(m)=sqrt((CL(m).x-SN(i).x)^2 + (CL(m).y-SN(i).y)^2);
            % we calculate the distance 'd' between the sensor node that is
            % transmitting and the cluster head that is receiving with the following equation+ 
            % d=sqrt((x2-x1)^2 + (y2-y1)^2) where x2 and y2 the coordinates of
            % the cluster head and x1 and y1 the coordinates of the transmitting node
            end
        d=d(1:CLheads); % fixing the size of "d" array
        [M,I]=min(d(:)); % finds the minimum distance of node to CH
        [Row, Col] = ind2sub(size(d),I); % displays the Cluster Number in which this node belongs too
        SN(i).cluster=Col; % assigns node to the cluster
        SN(i).dtch= d(Col); % assigns the distance of node to CH
        SN(i).chid=CL(Col).id;
        end
       end
                           %%%%%% Steady-State Phase %%%%%%
    if(mod(rnd,10)==2 || CHdead==1)
        figure(10)
        cla
      for i=1:n
           hold on;
           figure(10)
           plot(sinkx,sinky,'*r');
          if  (SN(i).role==0) && (SN(i).E>0) && (CLheads>0)
           plot(x,y,xm,ym,SN(i).x,SN(i).y,'ob');
          end
          if  (SN(i).role==1) && (SN(i).E>0) 
           plot(x,y,xm,ym,SN(i).x,SN(i).y,'xr');
          end
          if (SN(i).role == 1) && (SN(i).E > 0)
           plot([SN(i).x, sinkx], [SN(i).y, sinky], '-r'); % Đường màu đỏ
          end
          %for i = 1:n
            if (SN(i).role == 0) && (SN(i).E > 0) && (CLheads > 0) % if node is normal and has energy and there are Cluster Heads
            % Find the Cluster Head index corresponding to this node
            ch_index = SN(i).cluster;
        % Check if the node has a valid Cluster Head index
               if (ch_index > 0) && (ch_index <= CLheads)
            % Plot a black line between the node and its Cluster Head
               plot([SN(i).x, CL(ch_index).x], [SN(i).y, CL(ch_index).y], '-k');
               end
            end
          %end

           title 'Wireless Sensor Network';
           xlabel '(m)';
           ylabel '(m)';
      end
    end           
% Energy Dissipation for normal nodes %
    
    for i=1:n
       if (SN(i).cond==1) && (SN(i).role==0) && (CLheads>0)
       	if SN(i).E>0
            ETx= Eelec*k + Eamp * k * SN(i).dtch^2;
            SN(i).E=SN(i).E - ETx;
            energy=energy+ETx;
            
        % Dissipation for cluster head during reception
        if SN(SN(i).chid).E>0 && SN(SN(i).chid).cond==1 && SN(SN(i).chid).role==1
            ERx=(Eelec+EDA)*k;
            energy=energy+ERx;
            SN(SN(i).chid).E=SN(SN(i).chid).E - ERx;
             if SN(SN(i).chid).E<=0  % if cluster heads energy depletes with reception
                SN(SN(i).chid).cond=0;
                SN(SN(i).chid).rop=rnd;
                dead_nodes=dead_nodes +1;
                operating_nodes= operating_nodes - 1
                CHdead=1;
             end
        end
        end
        
        
        if SN(i).E<=0       % if nodes energy depletes with transmission
        dead_nodes=dead_nodes +1;
        operating_nodes= operating_nodes - 1
        SN(i).cond=0;
        SN(i).chid=0;
        SN(i).rop=rnd;
        end
        
      end
    end            
    
    
    
% Energy Dissipation for cluster head nodes %
   
   for i=1:n
     if (SN(i).cond==1)  && (SN(i).role==1)
         if SN(i).E>0
            ETx= (Eelec+EDA)*k + Eamp * k * SN(i).dts^2;
            SN(i).E=SN(i).E - ETx;
            energy=energy+ETx;
         end
         if  SN(i).E<=0     % if cluster heads energy depletes with transmission
         dead_nodes=dead_nodes +1;
         operating_nodes= operating_nodes - 1
         SN(i).cond=0;
         SN(i).rop=rnd;
         CHdead=1;
         end
     end
   end
   
  
    if operating_nodes<n && temp_val==0
        temp_val=1;
        flag1stdead=rnd
    end

    % Display Number of Cluster Heads of this round %
    %CLheads;
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
    transmissions=transmissions+1;
    if CLheads==0
    transmissions=transmissions-1;
    end
    
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
    rnd= rnd +1;
    
    tr(transmissions)=operating_nodes;
    op(rnd)=operating_nodes;
    
    totalenergy=totalenergy+energy;
    nrg(transmissions)=totalenergy;
    

end
sum=0;

for i=1:flag1stdead
    sum=nrg(i) + sum;
end
temp1=sum/flag1stdead;
temp2=temp1/n;
for i=1:flag1stdead
avg_node(i)=temp2;
end
    
% Plotting Simulation Results "Operating Nodes per Round" %figure(2)
plot(1:rnd, op(1:rnd), '-r', 'LineWidth', 2);

% Transpose the data to write in vertical format
T = table((1:rnd)', op(1:rnd)');
writetable(T, 'hfpso.txt')

title({'TEEN'; 'Operating Nodes per Round'});
xlabel('Rounds');
ylabel('Operational Nodes');
hold on;
    % Plotting Simulation Results  %
    figure(13)
    plot(1:transmissions,tr(1:transmissions),'-r','Linewidth',2);
    title ({'PSO'; 'Operational Nodes per Transmission';})
    xlabel 'Transmissions';
    ylabel 'Operational Nodes';
    hold on;
    
    % Plotting Simulation Results  %
    figure(14)
    plot(1:flag1stdead,nrg(1:flag1stdead),'-r','Linewidth',2);
    title ({'PSO'; 'Energy consumed per Transmission';})
    xlabel 'Transmission';
    ylabel 'Energy ( J )';
    hold on;
    
    % Plotting Simulation Results  %
    figure(15)
    plot(1:flag1stdead,avg_node(1:flag1stdead),'-r','Linewidth',2);
    title ({'PSO'; 'Average Energy consumed by a Node per Transmission';})
    xlabel 'Transmissions';
    ylabel 'Energy ( J )';
    hold on;

figure(16)
plot(1:rnd-1, total_energy(1:rnd-1), '-r', 'LineWidth', 2);
A = table((1:rnd-1)', total_energy(1:rnd-1)');
writetable(A, 'HFPSO_total_energy.txt')
title({'LEACH'; 'Total Energy used by all Active Nodes per Round'});
xlabel('Rounds');
ylabel('Total Energy ( J )');
hold on;

% Plotting total remaining energy of all active nodes per round
figure(17)
plot(1:rnd-1, remaining_energy(1:rnd-1), '-b', 'LineWidth', 2);
B = table((1:rnd-1)', remaining_energy(1:rnd-1)');
writetable(B, 'HFPSO_remain_energy.txt')
title({'LEACH'; 'Total Remaining Energy of All Active Nodes per Round'});
xlabel('Rounds');
ylabel('Total Remaining Energy ( J )');
hold on;
figure(18)
plot(1:rnd-1,energy_per_round(1:rnd-1),'-b','Linewidth',2);
C = table((1:rnd-1)', energy_per_round(1:rnd-1)');
writetable(C, 'HFPSO_energy_per_round.txt')
title ({'PEGASIS'; ' Energy per round All Nodes per Transmission';})
xlabel 'Transmissions';
ylabel 'Total Energy ( J )';
hold on;
disp(['Vòng mà nút đầu tiên chết: ', num2str(first_dead)]);
disp(['Vòng mà nút giữa chết: ', num2str(middle_dead)]);
disp(['Vòng mà nút cuối cùng chết: ', num2str(last_dead)]);