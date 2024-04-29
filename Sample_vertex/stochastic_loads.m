function [Pd] = stochastic_loads(mpc)
        % generate data for each loads(WT, PV), PD = 3
        Nbus = size(mpc.bus,1);
        Pd = mpc.bus(:,3);
        % mean value = forecast, standard deviation = 0.1*mean value
        for i = 1:Nbus
            Pd(i) = normrnd(Pd(i), 0.25 * abs(Pd(i)));
        end
end

