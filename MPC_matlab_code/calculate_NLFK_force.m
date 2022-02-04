%% Calculating tyhe NLFK force:
%
%  F_exc = sum( pressure * n * dS ); on the wetted surface of the cylinder.
%

function F_fk = calculate_NLFK_force(var, time, dX)

global wave

NDIM = length(dX);
F_fk = zeros(NDIM, 1);

cell_volume = prod(dX);

for DIM = 3 : NDIM
    
    x1 = min(var(:,DIM));
    x2 = max(var(:,DIM));
        
    x_sweep = x1 : dX(DIM) : x2;
    
    dS = cell_volume / dX(DIM);
    
    for i = 1 : length(x_sweep)-1

        idx1 = find(var(:,DIM) == x_sweep(i));
        idx2 = find(var(:,DIM) == x_sweep(i+1));

        for k = 1 : length(idx1)

            if ( var(idx1(k),4) >= 0 && var(idx2(k),4) < 0 && var(idx1(k),5) <= 0.0 )
                
                n = sign(var(idx2(k),4) - var(idx1(k),4));
                pressure_1 = calculate_pressure(var(idx1(k), 1), time, var(idx1(k),5));
                pressure_2 = calculate_pressure(var(idx2(k), 1), time, var(idx2(k),5));
                
                a = (pressure_1 - pressure_2)/(var(idx1(k),4) - var(idx2(k),4));
                b = pressure_1 - a * var(idx1(k),4);
                
                pressure_face = b;
                
                F_fk(DIM, 1) = F_fk(DIM, 1) + n * pressure_face * dS;
                
            elseif ( var(idx1(k),4) < 0 && var(idx2(k),4) >= 0 && var(idx1(k),5) <= 0.0 )
                
                n = sign(var(idx2(k),4) - var(idx1(k),4)); 
                pressure_1 = calculate_pressure(var(idx1(k), 1), time, var(idx1(k),5));
                pressure_2 = calculate_pressure(var(idx2(k), 1), time, var(idx2(k),5));
                
                a = (pressure_1 - pressure_2)/(var(idx1(k),4) - var(idx2(k),4));
                b = pressure_1 - a * var(idx1(k),4);
                
                pressure_face = b;
                
                F_fk(DIM, 1) = F_fk(DIM, 1) + n * pressure_face * dS;
                
            end
        end
    end
end