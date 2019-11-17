function [datout,flagout] = armafill(datin,flagin,aka,ind1,params,varargin)
% [datout,flagout] = gapfill(datin,flagin,aka,varargin) interpolates in the
% gaps contained in data using ARMA models for the predictors
% Input:    datin - input data array
%           flagin - status array. The gaps must be correctly flagged.
%           aka - Akaike coefficient matrix
%           ind1 - gap indexes
%           params - parameter list, is an array of 3 elements:
%              facmin - min. ratio between segment length and number of 
%                parameters for the model
%              facmax - 
%              npi - inf. limit in gap length for the ARMA interpolation
%                (below this limit linear interpolation is used)
%              pmin - inf. limit por the AR order
% Output:   datout - ARMA interpolated data series
%           flagout - residual status array

% Calls:   armaint.m v1.3.4
% Version: 
% Changes: minor bugs fixed
%  Author: Javier Pascual-Granado
%  $Date: 17/11/2019$
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L = length(datin);
l0 = length(ind1);

% convert data into column vector
datout = reshape(datin,L,1);

%% Parameters %%

facmin = params(1);
facmax = params(2);
npi = params(3);
pmin = params(4);
% rstd = params(5);

% Internal parameter S is related to the efficiency. When the length of the
% gap is > S times the length of any of the segments the algorithm is no 
%longer efficient for filling such a long gap and the gap is left unfilled.
S = 4;

% Akaike matrix is used to select the optimum (p,q) order. The lowest p is
% prioritized
minaka = min(min(aka));

[cp,cq] = find(aka==minaka);
q = cq - 1;
p = cp + pmin - 1;
ord = [p q];

% Static variables
aka1 = aka;
minaka1 = minaka;
ord1 = ord;
flagout = flagin;
l1 = l0;
%% Gap-filling process %% 

ind1f = l0-1;
i = 1;

if nargin==5,
    h = waitbar(0,sprintf(...
        'Gap filling iteration 1\n     Please wait...'));
else
    h = waitbar(0,sprintf(...
        'Gap filling iteration %d\n     Please wait...',varargin{1}));
end

% Here begins the gap-filling process
while i<=ind1f,
    
    % Data segments selection
    if l1==2,   % only one gap
        if ind1(i)==1,
            seg2 = datout(ind1(i+1)+1:end);
            seg1 = NaN;
        else
            if ind1(i+1)==L,
                seg2 = NaN;
            else
                seg2 = datout(ind1(i+1)+1:end);
            end
            seg1 = datout(1:ind1(i)-1);
        end
    else
        if ind1(i)==1,
            seg2 = datout(ind1(i+1)+1:ind1(i+2)-1);
            seg1 = NaN;
        elseif ind1(i+1)==L,
            seg1 = datout((ind1(i-1)+1):(ind1(i)-1));
            seg2 = NaN;
        else
            % First gap
            if i==1,
                seg1 = datout(1:ind1(i)-1);
            else
                seg1 = datout((ind1(i-1)+1):(ind1(i)-1));              
            end
            % Last gap
            if i==(l0-1),
                seg2 = datout(ind1(i+1)+1:end);
            else
                seg2 = datout(ind1(i+1)+1:ind1(i+2)-1);
            end
        end      
    end
    
    % number of loss datapoints in the gap
    np = ind1(i+1)-ind1(i)+1;
    
    waitbar(ind1(i)/L,h)

    lseg1 = length(seg1);
    lseg2 = length(seg2);

    if isempty(find(isnan(seg1),1)) && isempty(find(isnan(seg2),1)),
        % Small gaps are linearly interpolated
        if (np <= npi),
            lseg = lseg1+np+lseg2;
            interp = interp1([1:lseg1 (lseg1+np+1):lseg],...
                    [seg1' seg2'],lseg1+(1:np),'linear');
            if i==1,                       
                datout(ind1(1):ind1(2)) = interp;
    %             flagout(ind1(1):ind1(2)) = 1;
                flagout(ind1(1):ind1(2)) = 0;
                ind1(1:2)=[];
            else                                    
                datout(ind1(i):ind1(i+1)) = interp;
    %             flagout(ind1(i):ind1(i+1)) = 1;
                flagout(ind1(i):ind1(i+1)) = 0;
                ind1(i:(i+1))=[];
            end
            l1 = length(ind1);
            ind1f = l1-1;
            continue;
        end
    end
    
    d = sum(ord);
    fac1 = lseg1 / d;
    fac2 = lseg2 / d;

    if isempty(find(isnan(seg1),1)), 
        if isempty(find(isnan(seg2),1)),
        % Model order fix for small segments
            while (fac1 <= facmin || fac2 <= facmin),              
                aka(cp,cq) = inf;
                minaka = min(min(aka));
                if isinf(minaka)==1,
                    break;    
                end
                [cp,cq] = find(aka==minaka);
                q = cq - 1;
                p = cp + pmin - 1;
                ord = [p q];
                d = sum(ord);
                if length(d)>1,                   
                    if d(1)==d(2),
                        [~,iq] = min(q);
                        q = q(iq);
                        p = p(iq);
                    else
                        [~,id] = min(d);
                        q = q(id);
                        p = p(id);
                    end
                    ord = [p q];
                    d = sum(ord);
                end           
                fac1 = lseg1 / d;
                fac2 = lseg2 / d;
            end
        else
            while (fac1 <= facmin),              
                aka(cp,cq) = inf;
                minaka = min(min(aka));
                if isinf(minaka)==1,
                    break;    
                end
                [cp,cq] = find(aka==minaka);
                q = cq - 1;
                p = cp + pmin - 1;
                ord = [p q];
                d = sum(ord);
                if length(d)>1,                   
                    if d(1)==d(2),
                        [~,iq] = min(q);
                        q = q(iq);
                        p = p(iq);
                    else
                        [~,id] = min(d);
                        q = q(id);
                        p = p(id);
                    end
                    ord = [p q];
                    d = sum(ord);
                end           
                fac1 = lseg1 / d;
            end
        end
    elseif isempty(find(isnan(seg2),1)),
        while (fac2 <= facmin),              
            aka(cp,cq) = inf;
            minaka = min(min(aka));
            if isinf(minaka)==1,
                break;    
            end
            [cp,cq] = find(aka==minaka);
            q = cq - 1;
            p = cp + pmin - 1;
            ord = [p q];
            d = sum(ord);
            if length(d)>1,                   
                if d(1)==d(2),
                    [~,iq] = min(q);
                    q = q(iq);
                    p = p(iq);
                else
                    [~,id] = min(d);
                    q = q(id);
                    p = p(id);
                end
                ord = [p q];
                d = sum(ord);
            end           
            fac2 = lseg2 / d;
        end
    end
        
    if isinf(minaka)==1,
        if i==1,
%             i = 3;
            flagout(ind1(1):ind1(2)) = 1;
            ind1(1:2)=[];
            l1 = length(ind1);
            ind1f = l1-2;
        else
%             i = i+2;
            flagout(ind1(i):ind1(i+1)) = 1;
            ind1(i:(i+1))=[];
            l1 = length(ind1);
            ind1f = l1-2;
        end
        ord = ord1;
        cq = ord(2)+1;
        cp = ord(1)+1-pmin;
        aka = aka1;
        minaka = minaka1;
        continue;
    end

    % Too long segments are reduced by facmax for efficiency
    if isempty(find(isnan(seg1),1)),
        if (fac1 > facmax) && (floor(facmax*d)>np/S),
            newi = 1+floor(fac1-facmax)*d;
            seg1 = seg1(newi:end);
        end
        if isempty(find(isnan(seg2),1)),          
            if (fac2 > facmax) && (floor(facmax*d)>np/S),             
                newi = 1+facmax*d;
                seg2 = seg2(1:newi);
            end
        end
    elseif isempty(find(isnan(seg2),1)),
        if (fac2 > facmax) && (floor(facmax*d)>np/S),             
            newi = 1+facmax*d;
            seg2 = seg2(1:newi);
        end
    end 
        
    % Interpolation algorithm
    [int,go] = armaint(seg1,seg2,ord,np);
    if go,
        interp = int;
    else
        if i==1,
            flagout(ind1(1):ind1(2)) = 1;
            ind1(1:2)=[];
            l1 = length(ind1);
            ind1f = l1-2;
%             i = 3;
        else
            flagout(ind1(i):ind1(i+1)) = 1;
            ind1(i:(i+1))=[];
            l1 = length(ind1);
            ind1f = l1-2;
%             i = i+2;
        end
        % reinicialization
        ord = ord1;
        cq = ord(2)+1;
        cp = ord(1)+1-pmin;
        aka = aka1;
        minaka = minaka1;
        continue;
    end

    if i==1,
        datout(ind1(1):ind1(2)) = interp;
        flagout(ind1(1):ind1(2)) = 0;
        ind1(1:2)=[];
        l1 = length(ind1);
        ind1f = l1-2;
    else
        datout(ind1(i):ind1(i+1)) = interp;
        flagout(ind1(i):ind1(i+1)) = 0;
        ind1(i:(i+1))=[];
        l1 = length(ind1);
        ind1f = l1-2;
    end
    
    % reinicialization
    aka = aka1;
    minaka = minaka1;
    ord = ord1;
    cq = ord(2)+1;
    cp = ord(1)+1-pmin;
end

close(h);
