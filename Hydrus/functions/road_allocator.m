function RoadCostNOK = road_allocator(Dis2Road_km,P_MW)
% estimates cost of temporary access road based on NVE Veileder2 (2.7.1)
% and Veileder3 (B.12.1) assuming <10MW is low quality road while >10MW is
% high quality

% D = distance in km
% PW = power capacity in MW
% RoadCostNOK are in NOK as rates are NOK2000/m and NOK 1000/m
for i=1:length(P_MW)
    if P_MW(i)>10
        RoadCostNOK(i) = Dis2Road_km*2000*1000;  %High quality difficult terrain
    else
        RoadCostNOK(i) = Dis2Road_km*1000*1000; %Low quality difficult terrain
    end
end
end