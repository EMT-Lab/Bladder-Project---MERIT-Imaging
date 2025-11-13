C_in = Waist_mean;
major_diameter_in = SAD_mean;
a = major_diameter_in/2;
% minor_rad_in = 2;
% b = minor_rad_in;

% C_calc = pi*(3*(a+b)-sqrt((3*a+b)*(a+3*b)));
% C_in = C_calc;

b_calc_p = (-1*a + sqrt(a^2 - 6*(1.5*a^2-C_in^2/pi^2)))/3;

b_calc_m = (-1*a - sqrt(a^2 - 6*(1.5*a^2-C_in^2/pi^2)))/3;

C_calc = pi*(3*(a+b_calc_p)-sqrt((3*a+b_calc_p)*(a+3*b_calc_p)));