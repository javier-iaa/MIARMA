function [output,med,sigma] = stnorm(input)
% function [output,med,sigma] = stnorm(input)
% Statistical normalization or gaussianization of a time series. Input can
% be a filename containing 3 cols: time, data and status, or a variable 
% containing only the data.
% med is the mean and sigma is the standard deviation of the input
% Modified 10-sep-2015

if ischar(input),
    filename = input;
    strdata = importdata(filename);
    [~,m] = size(strdata.data);
    if m>3,
        fprintf('Warning: wrong format. Removing columns\n');
        strdata.data(:,3:m-1) = [];
    elseif m<3,
        fprintf('Warning: wrong format. Aborting program.\n');
        return;
    end
    time = strdata.data(:,1);
    data = strdata.data(:,2);
    flag = strdata.data(:,3);
    
    med = mean(data);
    sigma = std(data);
    output = (data-med)./sigma;
    
    fout = cat(2,filename(1:end-4),'_std.dat');
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

med = mean(data);
sigma = std(data);
output = (data-med)./sigma;
