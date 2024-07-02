fid = fopen('hfpso.txt', 'r') ;              
fgetl(fid) ;                                               
buffer = fread(fid, Inf) ;                    
fclose(fid)
fid = fopen('hfpso.txt', 'w')  ;   
fwrite(fid, buffer) ;                         
fclose(fid) ; 
hfpso = load('hfpso.txt');
fid = fopen('teen.txt', 'r') ;              
fgetl(fid) ;                                                        
buffer = fread(fid, Inf) ;                    
fclose(fid)
fid = fopen('teen.txt', 'w')  ;   
fwrite(fid, buffer) ;                         
fclose(fid) ;
teen = load('teen.txt');
fid = fopen('leach.txt', 'r') ;              
fgetl(fid) ;                                                            
buffer = fread(fid, Inf) ;                    
fclose(fid)
fid = fopen('leach.txt', 'w')  ;   
fwrite(fid, buffer) ;                         
fclose(fid) ;
leach = load('leach.txt');
fid = fopen('pegasis.txt', 'r') ;              
fgetl(fid) ;                                                     
buffer = fread(fid, Inf) ;                    
fclose(fid)
fid = fopen('pegasis.txt', 'w')  ;   
fwrite(fid, buffer) ;                         
fclose(fid) ;
pegasis = load('pegasis.txt');
plot(hfpso(:,1), hfpso(:,2), '-b', 'LineWidth', 2);
hold on;
plot(teen(:,1), teen(:,2), '-r', 'LineWidth', 2);
plot(leach(:,1), leach(:,2), '-g', 'LineWidth', 2);
plot(pegasis(:,1), pegasis(:,2), '-m', 'LineWidth', 2)
title('So sánh cảm biến còn hoạt động sau mỗi vòng');
xlabel('Vòng');
ylabel('Số nút cảm biến đang hoạt động');
legend('HFPSO', 'TEEN', 'LEACH', 'PEGASIS', 'Location', 'Best');
grid on;
hold off;