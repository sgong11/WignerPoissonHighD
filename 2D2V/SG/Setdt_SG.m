function [dt, t_new] = Setdt_SG(CFL,vmax,h_mesh,t_now,T)

    % dt = CFL/(vmax/h_mesh(1)+vmax/h_mesh(3));

    % E1_aver = norm(E1_ht)/sqrt(N_mesh(1)*N_mesh(3));
    % E2_aver = norm(E2_ht)/sqrt(N_mesh(1)*N_mesh(3));
    % dt = CFL/(vmax/h_mesh(1)+vmax/h_mesh(3)+E1_aver/h_mesh(2)+ E2_aver/h_mesh(4));
    % dt = CFL*min([h_mesh(1)/vmax, h_mesh(3)/vmax, h_mesh(2)/E1_aver, h_mesh(4)/E2_aver]);
    dt = CFL*min([h_mesh(1)/vmax, h_mesh(3)/vmax]);

    % flag = false;
    if t_now + dt > T
        dt = T - t_now;
        flag = true;
    end

    t_new = t_now + dt;

end