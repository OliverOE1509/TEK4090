function [d_bar] = getDisturbanceVector(ind_start,pred_h,AmbientCond,BuildingDemand)
%GETDISTURBANCEVECTOR inputs the predicted disturbance over the prediction
%horizon pred_hor

d_bar = zeros(11*pred_h,1);
jj = ind_start;

d_size = length([AmbientCond{jj,2:end}, BuildingDemand{jj,2}]);

% first entry is zero
for ii = 2:pred_h-1
    d_bar((ii-1)*d_size +1:ii*d_size) = ...
        [AmbientCond{jj,2:end}, BuildingDemand{jj,2}]';
    jj = jj+1;
end



end

