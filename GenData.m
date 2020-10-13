rng('shuffle');
pressure_drop_1 = input('Pressure 1 (psi) [10]:');
if isempty(pressure_drop_1)
    pressure_drop_1 = 10;
end
pressure_drop_2 = input('Pressure 2 (psi) [15]:');
if isempty(pressure_drop_2)
    pressure_drop_2 = 15;
end
pressure_drop_3 = input('Pressure 3 (psi) [20]:');
if isempty(pressure_drop_3)
    pressure_drop_3 = 20;
end

pressure_drop_array = [pressure_drop_1, pressure_drop_2, pressure_drop_3];
pressure_drop_array_psi = pressure_drop_array;
pressure_drop_array = pressure_drop_array * 6894.76; % Convert PSI to Pa

specific_cake_resistance = input('Target Cake Resistance (kg/m) [3e10]:');
if isempty(specific_cake_resistance)
    specific_cake_resistance = 3e10;
end
% Limits on the swing of input variables (%):
a = -5;
b = 5;
specific_cake_resistance = specific_cake_resistance+specific_cake_resistance*((b-a).*rand() + a)/100;
fprintf('Actual Cake Resistance Used: %.2e kg/m\n', specific_cake_resistance);

medium_resistance = input('Target Medium Resistance (1/m) [8e10]:');
if isempty(medium_resistance)
    medium_resistance = 8e10;
end
medium_resistance = medium_resistance+medium_resistance*((b-a).*rand() + a)/100;
fprintf('Actual Medium Resistance Used: %.2e 1/m\n', medium_resistance);

slurry_concentration = input('Target Slurry Concentration (kg/m3) [15]:');
if isempty(slurry_concentration)
    slurry_concentration = 15;
end
slurry_concentration_array = slurry_concentration+slurry_concentration*((b-a).*rand(1,3) + a)/100;
fprintf('Actual Slurry Concentration Used: %.2f kg/m^3\n', slurry_concentration_array);

% These are constant variables that are not going to change.
filtration_area = 0.0628; 
mu_water = 8.9e-4; 
rho_water = 997; 
rho_chalk = 2710;
volume_small = 0.00151;
large_diameter = 30/100;

value_holder = zeros(30,3);

for i = 1:3
a = 0;
b = 5;
r = (b-a).*rand() + a;

%dt_dv = (mu_water/(filtration_area*pressure_drop))*((specific_cake_resistance*slurry_concentration*volume_small*15/filtration_area)+medium_resistance);

    pressure_drop = pressure_drop_array(i);
    slurry_concentration = slurry_concentration_array(i);
time_mat = zeros(25,1);
final_pressure_drop = pressure_drop;
for time_step = 1:25
    if time_step < ceil(r)
        pressure_drop = final_pressure_drop+final_pressure_drop*((-20-0).*rand() + 0)/100; % For the first x points allow a big swing down in pressure
    else
        pressure_drop = final_pressure_drop+final_pressure_drop*((4).*rand() -2 )/100; % For subsequent points allow a only allow +/- 2% swing
    end
    k_c = (mu_water*slurry_concentration*specific_cake_resistance)/(pressure_drop*filtration_area^2);
    one_q0 = (mu_water*medium_resistance)/(filtration_area*pressure_drop);
    time = k_c/2 * (time_step * volume_small)^2 + one_q0 * (time_step * volume_small) - (k_c/2 * ((time_step-1) * volume_small)^2 + one_q0 * ((time_step-1) * volume_small));
    %time = time + 0.025*time*((b-a).*rand() + a);
    time_mat(time_step) = round(time,2);
end
%disp(time_mat);
% Plot the data so we can examine it.
scatter(cumsum(volume_small*ones(size(time_mat))), time_mat./volume_small)
hold on
% We need to back-work the concentrations: 
% Total volumetric flow will be 25 times total volume filled in small tank
% + margin.
volume_big = ((0.005).*rand() +1 )*25*(32/30)*volume_small;
height_big = round(volume_big/(pi()*(large_diameter^2/4))*100,1);

mass_solids = volume_big*slurry_concentration;
porosity = 0.1*rand()+0.60;

volume_cake = (mass_solids/(rho_chalk*(1-porosity)));
rho_cake = porosity*rho_water+(1-porosity)*rho_chalk;
mass_cake = rho_cake*volume_cake;

% We measured the volume in a 1L measuring cylinder!
volume_cake = volume_cake*1e6; % Convert to ML
max_water = 1000 - volume_cake;
volume_water = round(0.4*rand()*max_water+0.6*max_water, -1);
volume_total = round(volume_water+volume_cake, -1);

mass_tray = 0.078;

mass_total = round(mass_tray+mass_cake, 3);

% DUMP in a matrix:
value_holder(1,i) = pressure_drop_array_psi(i);
value_holder(2:26,i) = time_mat;
value_holder(27,i) = height_big;
value_holder(28,i) = mass_total;
value_holder(29,i) = volume_water;
value_holder(30,i) = volume_total;

end

% Export Data to Excel Sheet: 

% Pressure
xlswrite('data.xlsx',value_holder(1,1),1,'C4');
xlswrite('data.xlsx',value_holder(1,2),1,'F4');
xlswrite('data.xlsx',value_holder(1,3),1,'I4');

% Times
xlswrite('data.xlsx',value_holder(2:26,1),1,'C8');
xlswrite('data.xlsx',value_holder(2:26,2),1,'F8');
xlswrite('data.xlsx',value_holder(2:26,3),1,'I8');


% Filtrate Height
xlswrite('data.xlsx',value_holder(27,1),1,'C35');
xlswrite('data.xlsx',value_holder(27,2),1,'F35');
xlswrite('data.xlsx',value_holder(27,3),1,'I35');

% Cake + Tray Mass
xlswrite('data.xlsx',value_holder(28,1),1,'C38');
xlswrite('data.xlsx',value_holder(28,2),1,'F38');
xlswrite('data.xlsx',value_holder(28,3),1,'I38');

% Volumes
xlswrite('data.xlsx',value_holder(29:30,1),1,'C42');
xlswrite('data.xlsx',value_holder(29:30,2),1,'F42');
xlswrite('data.xlsx',value_holder(29:30,3),1,'I42');

xlswrite('data.xlsx',round(4*rand()+58,2),1,'L11');
