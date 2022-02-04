function y_predicted = AR_prediction(y, t, Np, order, dt_c)

%
%   Auto-regressive model to predict values in the future
%
%      y_predicted = Cp . y
%
%   |          |           |                 |
%   |   y(k)   |           |  y(k - order+1) |
%   |  y(k+1)  |           |  y(k - order+2) |
%   |     .    |  = [Cp] . |        .        |    
%   |     .    |           |        .        |
%   |     .    |           |      y(k-1)     |
%   |  y(k+Np) |           |       y(k)      |
%   |          |           |                 |
%
%   where,
%     y     --> available past data
%     t     --> time 
%     Np    --> values in future to be predicted (prediction horizon)
%     order --> AR model order (number of values from past spaced at intervals of dt_c)
%     dt_c  --> controller time
%     Cp    --> Prediction matrix of size( (Np+1)x(order) )
%

%% Identify the model:

t1 = t(end) : -dt_c : t(1);
if(t1(end) ~= t(1))
   t1(end) = t(1); 
end
time = flip(t1);

F = spline(t, y, time);       % Interpolate at equal intervals of dt_c.

past = 500;                   % number of data points to use, should be < length(F).
if( past >= length(F) )
    past = length(F) - 1;
end

dat = iddata(F(end-past:past)');
mod = ar(dat, order);

A = mod.A;                    % getting coeffs. from the AR-model.

%% Construct the prediction matrix:

B = [-A(2:order+1); [eye(order-1) zeros(order-1,1)]];
C = [1 zeros(1,order-1)];

Cptemp = [ 1  zeros(1,order-1) ];   % First row is for y(k)

for i = 1 : Np
    Cptemp = [ Cptemp ; C * B^i ];  % For y(k+1) to y(k+Np)
end

%% Re-order the matrix:

[l ,w] = size(Cptemp);        
Cp = [];                              
                               
for i = 0 : (w - 1)            
    Cp = [ Cp Cptemp(:,w-i) ];  
end                             
    
%% Predict the values:

[a ,b] = size(Cp);              % Here, a = (Np+1) and b = order.

y_past = F(end - b+1 : end);    % Taking past order number of data points.

y_predicted = Cp * y_past';     % predicting Np values in future.

end