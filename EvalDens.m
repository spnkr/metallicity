function [fval]= EvalDens(x,y,f,xgrid,ygrid,T,xn,yn,gap)
n = length(x);
xmin = min(xgrid); xmax = max(xgrid);
ymin = min(ygrid); ymax = max(ygrid);

fval= zeros(n,T);

for j = 1:n,
    kx = 30; ky = 17;
    for ix = 1:xn,
        if(x(j) < xgrid(ix) + gap) kx = ix; break; end
    end    
    for iy = 1:yn,
        if(y(j) < ygrid(iy) + gap) ky = iy; break; end
    end
    fval(j,:) = f(:,kx,ky);
end
    

    
end