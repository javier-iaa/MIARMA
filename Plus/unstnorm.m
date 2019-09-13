function output = unstnorm(input,m,s)
% function output = unstnorm(input,m,s)
% Undo of stnorm. Input can be a filename of data processed with stnorm.m
% or a variable containing the time series.
% m and s are the original mean and standard deviation respectively.

if ischar(input),
    filename = input;
    strdata = importdata(filename);
    time = strdata.data(:,1);
    data = strdata.data(:,2);
    flag = strdata.data(:,3);
    
    output = data*s+m;
    
    ind = find(filename=='.');
    fout = cat(2,filename(1:ind-5),filename(ind:end));
    fich = fopen(fout,'w');
    fprintf(fich,'x y z\n');
    for i=1:length(time),
        fprintf(fich,'%16.12f %16.13f %d\n',time(i),output(i),flag(i));
    end
    fclose(fich);
    return;
else
    data = input;
end

output = data*s+m;
