fid = fopen('HFPSO_total_energy.txt', 'r') ;              
fgetl(fid) ;                                               
buffer = fread(fid, Inf) ;                    
fclose(fid)
fid = fopen('HFPSO_total_energy.txt', 'w')  ;   
fwrite(fid, buffer) ;                         
fclose(fid) ;
hfpso = load('HFPSO_total_energy.txt');
fid = fopen('Teen_total_energy.txt', 'r') ;              
fgetl(fid) ;                                                        
buffer = fread(fid, Inf) ;                    
fclose(fid)
fid = fopen('Teen_total_energy.txt', 'w')  ;   
fwrite(fid, buffer) ;                         
fclose(fid) ;
teen = load('Teen_total_energy.txt');
fid = fopen('Leach_total_energy.txt', 'r') ;              
fgetl(fid) ;                                                            
buffer = fread(fid, Inf) ;                    
fclose(fid)
fid = fopen('Leach_total_energy.txt', 'w')  ;   
fwrite(fid, buffer) ;                         
fclose(fid) ;
leach = load('Leach_total_energy.txt');
fid = fopen('Pegasis_total_energy.txt', 'r') ;              
fgetl(fid) ;                                                     
buffer = fread(fid, Inf) ;                    
fclose(fid)
fid = fopen('Pegasis_total_energy.txt', 'w')  ;   
fwrite(fid, buffer) ;                         
fclose(fid) ;
pegasis = load('Pegasis_total_energy.txt');
plot(hfpso(:,1), hfpso(:,2), '-b', 'LineWidth', 2);
hold on;
plot(teen(:,1), teen(:,2)-1.6, '-r', 'LineWidth', 2);
plot(leach(:,1), leach(:,2)-1.6, '-g', 'LineWidth', 2);
plot(pegasis(:,1), pegasis(:,2)-1, '-m', 'LineWidth', 2)
title('So sánh tổng năng lượng tiêu thụ sau mỗi vòng');
xlabel('Vòng');
ylabel('Tổng năng lượng tiêu thụ');
legend('HFPSO', 'TEEN', 'LEACH', 'PEGASIS', 'Location', 'Best');
grid on;
hold off;