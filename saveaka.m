function saveaka(akaname, aka, pmin, pmax)

% akaname = 'temp_002_100_100.akc';
% pmax = strdata.params.pmax;
% pmin = 2;
% aka = strout.aka;

fichw = fopen(akaname,'w');

for i=0:pmax
    akam = aka(:,1+i);
    for j=pmin:pmax
        fprintf(fichw,'%f ', akam(j-pmin+1));
    end
    fprintf(fichw,'\n');
end

fclose(fichw);
