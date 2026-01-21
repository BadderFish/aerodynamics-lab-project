function [induced_U] = induced_velocity_fun(ControlPoint, Extreme_1, Extreme_2,Gamma)


toll = 1e-10;


r0 = Extreme_2 - Extreme_1;
r1 = ControlPoint - Extreme_1;
r2 = ControlPoint - Extreme_2;

CP = cross(r1, r2);

r1_norm = r1 / norm(r1);
r2_norm = r2 / norm(r2);


CP_norm = dot(CP, CP);

if(CP_norm < toll)
    CP_norm = toll;
end

dot_result = dot(r0, r1_norm - r2_norm);

induced_U = (Gamma/(4*pi)) * dot(r0, r1_norm - r2_norm) .* CP ./ (CP_norm);

