function datout = locdetrend(datout, datin, reco, seg1, seg2, interp)
% function datout = locdetrend(datout, datin, reco, seg1, seg2, interp)
% Fix local trends that might be not modeled properly introducing jumps

nd = 10;
L = length(datin);
datin = reshape(datin,L,1);
len = length(seg1);

if len<nd
    difs = [diff(seg1)' diff(seg2)'];
else
    difs = [diff(seg1((len-nd+1):len))' diff(seg2(1:nd))'];
end
lim_jump = 3*mean( abs(difs) );

jump = find(abs(datout(reco)-datin(reco))>lim_jump, 1);

if jump
    % Fit model
    if len>npint
        ssubi1 = 1:npint;
        sseg1 = seg1((len-npint+1):len);
        ssubi2ini = npint + np + 1;
        ssubi2fin = ssubi2ini + npint - 1;
        ssubi2 = ssubi2ini:ssubi2fin;
        sseg2 = seg2(1:npint);
    else
        ssubi1 = 1:len;
        sseg1 = seg1;
        ssubi2ini = len + np + 1;
        ssubi2fin = ssubi2ini + len - 1;
        ssubi2 = ssubi2ini:ssubi2fin;
        sseg2 = seg2;
    end
    data = [sseg1; datin(reco); sseg2];
    treco = npint + reco0;
    tdata = [ssubi1 treco ssubi2]';
    pp_data = polyfit(tdata, data, 3);
    t_tot = 1:ssubi2fin;
    pval_tot = polyval(pp_data, t_tot);
    subiint = (ssubi1(end)+1):(ssubi2(1)-1);
    pval_int = pval_tot( subiint) ;

    % detrending of interp
    pp_interp = polyfit( subiint, interp', 2);
    pval_interp = polyval(pp_interp, subiint);
    interp = pval_int' + interp - pval_interp';
    datout(ind1(1):ind1(2)) = interp;
end