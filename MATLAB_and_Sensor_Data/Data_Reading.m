%% 
% =========================================================================
% 9 DOF AHRS & FULL KALMAN FILTER SIMULATION
% =========================================================================
clear; clc; close all;

%% 1. READ DATA
disp("1. Reading CSV File...");
file_name = 'sensor_data.csv';
data_matrix = readmatrix(file_name);

% STM32 Format: ax, ay, az, gx, gy, gz, mx, my, mz, pressure
ax = data_matrix(:, 1);  ay = data_matrix(:, 2);  az = data_matrix(:, 3);
gx = data_matrix(:, 4);  gy = data_matrix(:, 5);  gz = data_matrix(:, 6);
mx = data_matrix(:, 7);  my = data_matrix(:, 8);  mz = data_matrix(:, 9);
pressure = data_matrix(:, 10);

num_samples = length(ax);
dt = 0.62; % Sampling time 

% MPU6050 axis correction (Breadboard orientation)
temp_ax = ax; ax = ay; ay = temp_ax;
temp_gx = gx; gx = gy; gy = temp_gx;

% Fix Right-Hand Rule mismatch for Gyro
gx = -gx;
gy = -gy;

%% --- HARDWARE AXIS ALIGNMENT ---
mx_aligned = my;   
my_aligned = -mx;  
mz_aligned = mz;   

mx = mx_aligned;
my = my_aligned;
mz = mz_aligned;

%% 2. AUTO CALIBRATION
disp("2. Calibrating Sensors...");
% Dinamik Kalibrasyon Sınırı (Kısa kayıtlarda kodun çökmesini engeller)
offset_limit = min(50, num_samples); 

ax_offset = mean(ax(1:offset_limit)); ay_offset = mean(ay(1:offset_limit)); az_offset = mean(az(1:offset_limit)) - 1.0;
gx_offset = mean(gx(1:offset_limit)); gy_offset = mean(gy(1:offset_limit)); gz_offset = mean(gz(1:offset_limit));

ax_cal = ax - ax_offset; ay_cal = ay - ay_offset; az_cal = az - az_offset;
gx_cal = gx - gx_offset; gy_cal = gy - gy_offset; gz_cal = gz - gz_offset;

%% 3. CALCULATE RAW ANGLES
disp("3. Calculating Angles...");
% Accelerometer angles
roll_acc = atan2d(ay_cal, az_cal);
pitch_acc = atan2d(-ax_cal, sqrt(ay_cal.^2 + az_cal.^2));

% Gyroscope angles
roll_gyro = cumsum(gx_cal) * dt;
pitch_gyro = cumsum(gy_cal) * dt;

%% 4. KALMAN FILTER (ROLL & PITCH)
disp("4. Applying Kalman Filter for X and Y Axes...");
Q_angle = 0.005; 
Q_bias = 0.003; 
R_measure = 0.1;

roll_kalman = zeros(1, num_samples);
pitch_kalman = zeros(1, num_samples);

roll_kalman(1) = roll_acc(1);
pitch_kalman(1) = pitch_acc(1);

P_roll = eye(2);
P_pitch = eye(2);

x_roll = [roll_acc(1); 0];
x_pitch = [pitch_acc(1); 0];

for k = 2:num_samples
    % -- ROLL FILTER --
    % Prediction
    rate_roll = gx_cal(k) - x_roll(2);
    x_roll(1) = x_roll(1) + dt * rate_roll;
    P_roll(1,1) = P_roll(1,1) + dt * (dt*P_roll(2,2) - P_roll(1,2) - P_roll(2,1) + Q_angle);
    P_roll(1,2) = P_roll(1,2) - dt * P_roll(2,2);
    P_roll(2,1) = P_roll(2,1) - dt * P_roll(2,2);
    P_roll(2,2) = P_roll(2,2) + Q_bias * dt;
    
    % Update
    y_roll = roll_acc(k) - x_roll(1);
    S_roll = P_roll(1,1) + R_measure;
    K_roll = [P_roll(1,1) / S_roll; P_roll(2,1) / S_roll];
    x_roll = x_roll + K_roll * y_roll;
    
    P00 = P_roll(1,1); P01 = P_roll(1,2); P10 = P_roll(2,1); P11 = P_roll(2,2);
    P_roll(1,1) = P00 - K_roll(1)*P00;
    P_roll(1,2) = P01 - K_roll(1)*P01;
    P_roll(2,1) = P10 - K_roll(2)*P00;
    P_roll(2,2) = P11 - K_roll(2)*P01;
    
    roll_kalman(k) = x_roll(1);
    
    % -- PITCH FILTER --
    % Prediction
    rate_pitch = gy_cal(k) - x_pitch(2);
    x_pitch(1) = x_pitch(1) + dt * rate_pitch;
    P_pitch(1,1) = P_pitch(1,1) + dt * (dt*P_pitch(2,2) - P_pitch(1,2) - P_pitch(2,1) + Q_angle);
    P_pitch(1,2) = P_pitch(1,2) - dt * P_pitch(2,2);
    P_pitch(2,1) = P_pitch(2,1) - dt * P_pitch(2,2);
    P_pitch(2,2) = P_pitch(2,2) + Q_bias * dt;
    
    % Update
    y_pitch = pitch_acc(k) - x_pitch(1);
    S_pitch = P_pitch(1,1) + R_measure;
    K_pitch = [P_pitch(1,1) / S_pitch; P_pitch(2,1) / S_pitch];
    x_pitch = x_pitch + K_pitch * y_pitch;
    
    P00 = P_pitch(1,1); P01 = P_pitch(1,2); P10 = P_pitch(2,1); P11 = P_pitch(2,2);
    P_pitch(1,1) = P00 - K_pitch(1)*P00;
    P_pitch(1,2) = P01 - K_pitch(1)*P01;
    P_pitch(2,1) = P10 - K_pitch(2)*P00;
    P_pitch(2,2) = P11 - K_pitch(2)*P01;
    
    pitch_kalman(k) = x_pitch(1);
end

%% 4.5 MAGNETOMETER TILT COMPENSATION
disp("Calculating Magnetometer Tilt Compensation...");
% Hard-Iron calibration 
offset_mx = -440.5000;
offset_my = -1197.0;
offset_mz = 999.5000;

mx_cal = mx - offset_mx;
my_cal = my - offset_my;
mz_cal = mz - offset_mz;

mx_norm = sqrt(mx_cal.^2 + my_cal.^2 + mz_cal.^2);
mx_cal = mx_cal ./ mx_norm;
my_cal = my_cal ./ mx_norm;
mz_cal = mz_cal ./ mx_norm;

yaw_mag = zeros(1, num_samples);

for k = 1:num_samples
    phi = deg2rad(roll_kalman(k));
    theta = deg2rad(pitch_kalman(k));
    
    X_h = mx_cal(k) * cos(theta) + my_cal(k) * sin(phi) * sin(theta) + mz_cal(k) * cos(phi) * sin(theta);
    Y_h = my_cal(k) * cos(phi) - mz_cal(k) * sin(phi);
    
    yaw_mag(k) = atan2d(Y_h, X_h);
    
    if yaw_mag(k) > 180
        yaw_mag(k) = yaw_mag(k) - 360;
    elseif yaw_mag(k) < -180
        yaw_mag(k) = yaw_mag(k) + 360;
    end
end

%% YAW KALMAN (GYRO + MAGNETOMETER)
R_measure_yaw = 50.0;
disp("Applying Yaw Kalman Filter...");

yaw_kalman = zeros(1, num_samples);
P_yaw = eye(2);
x_yaw = [yaw_mag(1); 0];
yaw_kalman(1) = yaw_mag(1);

for k = 2:num_samples
    % Prediction
    rate = gz_cal(k) - x_yaw(2);
    x_yaw(1) = x_yaw(1) + dt * rate;
    P_yaw(1,1) = P_yaw(1,1) + dt*(dt*P_yaw(2,2) - P_yaw(1,2) - P_yaw(2,1) + Q_angle);
    P_yaw(1,2) = P_yaw(1,2) - dt*P_yaw(2,2);
    P_yaw(2,1) = P_yaw(2,1) - dt*P_yaw(2,2);
    P_yaw(2,2) = P_yaw(2,2) + Q_bias*dt;
    
    % Measurement Update
    innovation = yaw_mag(k) - x_yaw(1);
    
    % Angle wrap correction
    if innovation > 180
        innovation = innovation - 360;
    elseif innovation < -180
        innovation = innovation + 360;
    end
    
    S = P_yaw(1,1) + R_measure_yaw;
    K = [P_yaw(1,1)/S; P_yaw(2,1)/S];
    x_yaw = x_yaw + K * innovation;
    
    P00 = P_yaw(1,1); P01 = P_yaw(1,2); P10 = P_yaw(2,1); P11 = P_yaw(2,2);
    P_yaw(1,1) = P00 - K(1)*P00;
    P_yaw(1,2) = P01 - K(1)*P01;
    P_yaw(2,1) = P10 - K(2)*P00;
    P_yaw(2,2) = P11 - K(2)*P01;
    
    yaw_kalman(k) = x_yaw(1);
    
    if yaw_kalman(k) > 180
        yaw_kalman(k) = yaw_kalman(k) - 360;
    elseif yaw_kalman(k) < -180
        yaw_kalman(k) = yaw_kalman(k) + 360;
    end
end

yaw_gyro = yaw_mag(1) + cumsum(gz_cal) * dt;
yaw_gyro = mod(yaw_gyro + 180, 360) - 180;

%% 5. ALTITUDE CALCULATION AND Z-AXIS KALMAN FILTER
disp("5. Calculating Barometric Altitude and Filtering Z Axis...");
P0 = mean(pressure(1:offset_limit)); % Ground pressure reference (Dinamik sınır eklendi)
alt_raw = 44330 * (1 - (pressure / P0).^(1/5.255));

alt_kalman = zeros(1, num_samples);
P_z = 1.0; 
Q_z = 0.0001; 
R_z = 100.0;  

x_z = 0; 

for k = 2:num_samples
    % Predict
    P_z = P_z + Q_z;
    
    % Update
    K_z = P_z / (P_z + R_z);
    x_z = x_z + K_z * (alt_raw(k) - x_z);
    P_z = (1 - K_z) * P_z;
    
    alt_kalman(k) = x_z;
end

%% 6. PLOTS (FILTER PERFORMANCE)
figure('Name', 'Drone Sensor Fusion Performance', 'Position', [100, 100, 1000, 800]);

subplot(3,1,1);
plot(roll_acc, 'Color', [0.7 0.7 1], 'DisplayName', 'Accelerometer'); hold on;
plot(roll_gyro, 'r', 'DisplayName', 'Gyroscope');
plot(roll_kalman, 'g', 'LineWidth', 2.5, 'DisplayName', 'Kalman');
title('Drone Roll Angle'); ylabel('Angle (Degrees)');
legend; grid on;

subplot(3,1,2);
plot(pitch_acc, 'Color', [0.7 0.7 1], 'DisplayName', 'Accelerometer'); hold on;
plot(pitch_gyro, 'r', 'DisplayName', 'Gyroscope');
plot(pitch_kalman, 'g', 'LineWidth', 2.5, 'DisplayName', 'Kalman');
title('Drone Pitch Angle'); ylabel('Angle (Degrees)');
legend; grid on;

subplot(3,1,3);
plot(alt_raw * 100, 'Color', [0.7 0.7 1], 'DisplayName', 'Barometer'); hold on;
plot(alt_kalman * 100, 'g', 'LineWidth', 2.5, 'DisplayName', 'Kalman Z-Axis');
title('Drone Altitude (Z-Axis)'); xlabel('Sample'); ylabel('Centimeters (cm)');
legend; grid on;

figure;
plot(yaw_mag, 'Color', [0.7 0.7 1], 'DisplayName', 'Magnetometer'); hold on;
plot(yaw_gyro, 'r', 'DisplayName', 'Gyroscope');
plot(yaw_kalman, 'g', 'LineWidth', 2, 'DisplayName', 'Kalman');
legend; title("Yaw Kalman"); ylabel("Degrees"); xlabel("Sample");
grid on; ylim([-180 180]);

%% 7. 3D DRONE DIGITAL TWIN (FLIGHT SIMULATION)
disp("6. Starting 3D Drone Flight Simulation...");
figure('Name', 'Drone 3D Flight Simulation', 'Position', [150, 150, 800, 600]);
view(3);

axis([-0.3 0.3 -0.3 0.3 -0.2 0.2]);
grid on;
xlabel('X Axis (m)'); ylabel('Y Axis (m)'); zlabel('Altitude (m)');
title('Drone 3D Orientation and Flight Simulation');

% BREADBOARD FİZİKSEL ÖLÇÜLERİ (Metre Cinsinden)
L = 0.16;   % 16 cm Uzunluk
W = 0.055;  % 5.5 cm Genişlik
H = 0.005;  % 0.5 cm Yükseklik/Kalınlık

vertices_base = [
    -L/2, -W/2, -H/2;  L/2, -W/2, -H/2;  L/2,  W/2, -H/2; -L/2,  W/2, -H/2;
    -L/2, -W/2,  H/2;  L/2, -W/2,  H/2;  L/2,  W/2,  H/2; -L/2,  W/2,  H/2
];

faces = [1 2 3 4; 5 6 7 8; 1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8];

drone_patch = patch('Vertices', vertices_base, 'Faces', faces, ...
    'FaceColor', [0.8 0.3 0.2], 'EdgeColor', 'k', 'FaceAlpha', 0.9);

% Animation Loop
for k = 1:num_samples
    % Apply negative sign to fix mirror image orientations
    phi = deg2rad(-roll_kalman(k));
    theta = deg2rad(-pitch_kalman(k));
    psi = deg2rad(yaw_kalman(k) - yaw_kalman(1));
    
    z_trans = alt_kalman(k);
    
    Rx = [1, 0, 0;
        0, cos(phi), -sin(phi);
        0, sin(phi), cos(phi)];
        
    Ry = [cos(theta), 0, sin(theta);
        0, 1, 0;
        -sin(theta), 0, cos(theta)];
        
    Rz = [cos(psi), -sin(psi), 0;
        sin(psi), cos(psi),  0;
        0,        0,         1];
        
    R = Rz * Ry * Rx;
    rotated_vertices = (R * vertices_base')';
    translated_vertices = rotated_vertices + [0, 0, z_trans];
    
    set(drone_patch, 'Vertices', translated_vertices);
    drawnow limitrate;
    pause(0.5);
end

disp("Drone Flight Simulation Completed!");