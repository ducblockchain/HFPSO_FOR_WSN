N=100;                       % number of nodes
area=[100,100];              % nodes deployment area in meter
Trange=7;                   % transmission range of sensor node in meter
nodes.pos=area(1).*rand(N,2);% nodes geographical locations
lambda=0.125;                % signal wavelength in meter
nodes.major = Trange;        % major axis for ellpitical range in meter
nodes.minor = lambda*Trange;  % minro axis for ellipitical range in meter
% redundantNo=9;               % number of healing nodes   
redundantNo=round(10*N/100);
data = readtable('random.txt');
nodes.pos = [data.Var1, data.Var2];
%% plot the nodes deployment
cnt=1;
for ii=1:N      
    for jj=1:N
        if ii~=jj
            nodes.distance(ii,jj)=pdist([nodes.pos(ii,:);nodes.pos(jj,:)]);
            if nodes.distance(ii,jj)<Trange || nodes.distance(ii,jj)==Trange
                nodes.inrange(ii,jj)=1;
            else
                nodes.inrange(ii,jj)=0;
            end
        end
    end
end
figure;
F5=plot(nodes.pos(:,1),nodes.pos(:,2),'.','color','r');
hold on
