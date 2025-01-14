fid = fopen('HFPSO_remain_energy.txt', 'r') ;              
fgetl(fid) ;                                               
buffer = fread(fid, Inf) ;                    
fclose(fid)
fid = fopen('HFPSO_remain_energy.txt', 'w')  ;   
fwrite(fid, buffer) ;                         
fclose(fid) ;
hfpso = load('HFPSO_remain_energy.txt');
fid = fopen('Teen_remain_energy.txt', 'r') ;              
fgetl(fid) ;                                                        
buffer = fread(fid, Inf) ;                    
fclose(fid)
fid = fopen('Teen_remain_energy.txt', 'w')  ;   
fwrite(fid, buffer) ;                         
fclose(fid) ;
teen = load('Teen_remain_energy.txt');
fid = fopen('Leach_remain_energy.txt', 'r') ;              
fgetl(fid) ;                                                            
buffer = fread(fid, Inf) ;                    
fclose(fid)
fid = fopen('Leach_remain_energy.txt', 'w')  ;   
fwrite(fid, buffer) ;                         
fclose(fid) ;
leach = load('Leach_remain_energy.txt');
fid = fopen('Pegasis_remain_energy.txt', 'r') ;              
fgetl(fid) ;                                                     
buffer = fread(fid, Inf) ;                    
fclose(fid)
fid = fopen('Pegasis_remain_energy.txt', 'w')  ;   
fwrite(fid, buffer) ;                         
fclose(fid) ;
pegasis = load('Pegasis_remain_energy.txt');
plot(hfpso(:,1), hfpso(:,2), '-b', 'LineWidth', 2);
hold on;
plot(teen(:,1), teen(:,2), '-r', 'LineWidth', 2);
plot(leach(:,1), leach(:,2), '-g', 'LineWidth', 2);
plot(pegasis(:,1), pegasis(:,2), '-m', 'LineWidth', 2)
title('So sánh năng lượng còn lại sau mỗi vòng');
xlabel('Vòng');
ylabel('Năng lượng còn lại của tất cả các nút cảm biến');
legend('HFPSO', 'TEEN', 'LEACH', 'PEGASIS', 'Location', 'Best');
grid on;
hold off;