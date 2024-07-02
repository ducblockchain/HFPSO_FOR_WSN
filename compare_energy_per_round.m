fid = fopen('HFPSO_energy_per_round.txt', 'r') ;              
fgetl(fid) ;                                               
buffer = fread(fid, Inf) ;                    
fclose(fid)
fid = fopen('HFPSO_energy_per_round.txt', 'w')  ;   
fwrite(fid, buffer) ;                         
fclose(fid) ;
hfpso = load('HFPSO_energy_per_round.txt');
fid = fopen('Teen_energy_per_round.txt', 'r') ;              
fgetl(fid) ;                                                        
buffer = fread(fid, Inf) ;                    
fclose(fid)
fid = fopen('Teen_energy_per_round.txt', 'w')  ;   
fwrite(fid, buffer) ;                         
fclose(fid) ;
teen = load('Teen_energy_per_round.txt');
fid = fopen('Leach_energy_per_round.txt', 'r') ;              
fgetl(fid) ;                                                            
buffer = fread(fid, Inf) ;                    
fclose(fid)
fid = fopen('Leach_energy_per_round.txt', 'w')  ;   
fwrite(fid, buffer) ;                         
fclose(fid) ;
leach = load('Leach_energy_per_round.txt');
fid = fopen('Pegasis_energy_per_round.txt', 'r') ;              
fgetl(fid) ;                                                     
buffer = fread(fid, Inf) ;                    
fclose(fid)
fid = fopen('Pegasis_energy_per_round.txt', 'w')  ;   
fwrite(fid, buffer) ;                         
fclose(fid) ;
pegasis = load('Pegasis_energy_per_round.txt');
plot(hfpso(:,1), hfpso(:,2), '-b', 'LineWidth', 2);
hold on;
plot(teen(:,1), teen(:,2), '-r', 'LineWidth', 2);
plot(leach(:,1), leach(:,2), '-g', 'LineWidth', 2);
plot(pegasis(:,1), pegasis(:,2), '-m', 'LineWidth', 2);
plot(hfpso(:,1), hfpso(:,2), '-b', 'LineWidth', 2);
title('So sánh năng lượng tiêu thụ mỗi vòng');
xlabel('Vòng');
ylabel('Năng lượng tiêu thụ mỗi vòng');
legend('HFPSO', 'TEEN', 'LEACH', 'PEGASIS', 'Location', 'Best');
grid on;
hold off;