%% Calculate pipeplenght
% Based on columns rows and elevation diff

for  k = 1:numel(outlets)
    
    if isnan(ro(k))==1; continue; end;
    
    for l = 1:ni
        for m = 1:numel(Pinlet_win{k}{l})
            
            PL{k}{l}(m) = Lenght_finder(acc_inlet_win{k},Zout,Pinlet_win{k}{l}(m),latOut(k),Zoutlets(k),ZPinlet{k}{l}(m));
            
            [latP{k}{l}(m), lonP{k}{l}(m)] = setltln(acc, Rw, rPin{k}{l}(m), cPin{k}{l}(m)); %Coordinates inlets
            %fprintf('Total pipelength = %.2f km.\n',PL{k}{l}(m)/1000);
            
        end
    end
end