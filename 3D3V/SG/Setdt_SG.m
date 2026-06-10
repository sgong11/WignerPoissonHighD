function [dt, t_new] = Setdt_SG(CFL,vmax,h_mesh,t_now,T)

    dt = CFL*min([h_mesh(1)/vmax, h_mesh(3)/vmax, h_mesh(5)/vmax]);

    % flag = false;
    if t_now + dt > T
        dt = T - t_now;
        flag = true;
    end

    t_new = t_now + dt;

end