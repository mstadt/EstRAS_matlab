function est = get_estrogen(t)
    a_e = 50; % age of estrogen decline
    tau_e = 2.6; % slope of decline
    AGE = (t/24)/365; % t is in hours, convert to years
    if AGE < (a_e)
        est = 1.0;
    elseif AGE >= a_e % estrogen decline
        est = 1 / (1.0 + (AGE-a_e)/tau_e);
    else
        fprintf('something went wrong... AGE = %f \n', AGE)
    end
end