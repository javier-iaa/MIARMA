function saveout(strout)
% function saveout(strout)
% Save the output of MIARMA in ASCII file
%
% By Javier Pascual-Granado
% <a href="matlab:web http://www.iaa.es;">IAA-CSIC, Spain</a>
%
% Version: 0.1
%
% Date: 16/03/2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isfield(strout, 'filename')
    filename = strout.filename;
    fout = [ filename(1:end-4) '.agfs'];
else
    fout = 'output.agfs';
end

fich = fopen(fout,'w');

fprintf(fich,'# code: MIARMA\n');
fprintf(fich,'# version: %s\n', strout.numvers);

if isfield(strout, 'ord')
    fprintf(fich,'# model (p,q): %d %d\n',strout.ord(1),strout.ord(2));
end

fprintf(fich,'# length: %d\n', strout.L);
fprintf(fich,'# gaps_arma: %d\n', strout.lgaps0-strout.Llin);
fprintf(fich,'# gaps_linear: %d\n', strout.Llin);
fprintf(fich,'# facmin: %d\n', strout.params.facmin);
fprintf(fich,'# facmax: %d\n', strout.params.facmax);
fprintf(fich,'# facint: %d\n', strout.params.facint);
fprintf(fich,'# npi: %d\n', strout.params.npi);
fprintf(fich,'# npz: %d\n', strout.params.npz);
fprintf(fich,'# mseg: %d\n', strout.params.mseg);
fprintf(fich,'# ft_corr: %d\n\n', strout.ft_corr);
fprintf(fich,'x y z\n');

for i=1:strout.L
    fprintf(fich,'%16.12f %16.13f %d\n',...
        strout.timeout(i),strout.datout(i),strout.statout(i));
end

fclose(fich);

fprintf('\n  Interpolation finished successfully.  \n');