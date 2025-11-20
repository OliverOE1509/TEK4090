function [t_lb, t_ub] = get_temp_const(start_time, pred_h, params, AmbientCond)
%GET_TEMP_CONST: input the start and end time of a horizon and outputs the
%temperature constraints at that time

% index of start time
jj = find(AmbientCond.timestamp==start_time);



if params.ctrl_umar == true && params.ctrl_dfab == true
    nx_u = 3 + 7;
    nx_t = nx_u + 1;
    t_lb = zeros(pred_h*nx_t,1);
    t_ub = t_lb;  
    
elseif params.ctrl_umar == false && params.ctrl_dfab == true
    nx_u = 7;
    nx_t = nx_u + 1;
    t_lb = zeros(pred_h*nx_t);
    t_ub = t_lb;  
elseif params.ctrl_umar == true && params.ctrl_dfab == false
    nx_u = 3;
    nx_t = nx_u + 1;
    t_lb = zeros(pred_h*nx_t);
    t_ub = t_lb;  
else
    error('You arent controlling any units!')
end


for ii = 1:pred_h
    date = AmbientCond{jj,1};
    hr = hour(date);
    if hr >=8 && hr <=20
        t_lb((ii-1)*nx_t +1:ii*nx_t-1) = params.T_min_day*ones(nx_u,1);
        t_ub((ii-1)*nx_t +1:ii*nx_t-1) = params.T_max_day*ones(nx_u,1);
    else
        t_lb((ii-1)*nx_t +1:ii*nx_t-1) = params.T_min_ngt*ones(nx_u,1);
        t_ub((ii-1)*nx_t +1:ii*nx_t-1) = params.T_min_ngt*ones(nx_u,1);
    end
    t_lb(ii*nx_t) = params.tank_T_min;
    t_ub(ii*nx_t) = params.tank_T_max;
    jj = jj+1;
end

end

