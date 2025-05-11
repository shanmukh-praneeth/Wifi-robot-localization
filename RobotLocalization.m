function room = createWiFiEnvironment()
    % Room Dimensions
    room.length = 10;    % Room length in meters
    room.width = 8;      % Room width in meters
    room.height = 3;     % Room height in meters

    num_rows = 2;
    num_cols = 2;
    lambda = 0.125;
    d = lambda/2;
    
    room.routers = struct('id',{},'position',{},'transmitPower',{},'antennaPositions',{});
    router_positions = [2.5,2,2.97;-2.5,2,2.97;-2.5,-2,2.97;2.5,-2,2.97];

    for i = 1:4
        antenna_positions = zeros(num_rows * num_cols, 3);
        index = 1;
        for x = 0:num_rows-1
            for y = 0:num_cols-1
                x_offset = (y - (num_cols - 1)/2)*d;
                y_offset = (x - (num_rows - 1)/2)*d;
                antenna_positions(index, :) = router_positions(i,:) + [x_offset, y_offset, 0];
                index = index + 1;
            end
        end
        room.routers(i) = struct('id', i, 'position', router_positions(i,:), 'transmitPower', 30, 'antennaPositions', antenna_positions);
    end
    room.robot.position = [0, 0, 0];
end

function visualizeEnvironment(room, estimated_pos)
    figure;
    hold on;
    grid on;
    
    % Plot routers and antennas
    for i = 1:length(room.routers)
        router = room.routers(i);
        plot3(router.position(1), router.position(2), router.position(3), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
        text(router.position(1), router.position(2), router.position(3), sprintf('Router %d', router.id));
        
        % Plot antennas
        for j = 1:size(router.antennaPositions, 1)
            plot3(router.antennaPositions(j, 1), router.antennaPositions(j, 2), router.antennaPositions(j, 3), 'kx', 'MarkerSize', 8);
        end
    end

    plot3(estimated_pos(1), estimated_pos(2), estimated_pos(3), 'gd', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
    text(estimated_pos(1), estimated_pos(2), estimated_pos(3), 'Estimated Position', 'VerticalAlignment','bottom');
    
    % Plot robot
    plot3(room.robot.position(1), room.robot.position(2), room.robot.position(3), 'bs', 'MarkerSize', 10, 'MarkerFaceColor', 'b');
    text(room.robot.position(1), room.robot.position(2), room.robot.position(3), 'Robot');
    xlabel('X (m)');
    ylabel('Y (m)');
    zlabel('Z (m)');
    title('Wi-Fi Localization Environment with URA Antennas');
    axis equal;
end

function rssi = computeRSSI(room)
    d0 = 1;  % Reference Distance
    PL0 = 40; % Path loss at the reference distance
    n = 3;  % Path loss exponent for indoor environment
    sigma = 1; 
    materialLoss = 0;
    robot_pos = room.robot.position;
    rssi = zeros(1, length(room.routers));
    for i = 1:length(room.routers)
        router = room.routers(i);
        d = norm(robot_pos - router.position);
        if d < d0
            d = d0;
        end
        noise = sqrt(sigma/2) * randn();  
        rssi(i) = router.transmitPower - PL0 - 10 * n * log10(d / d0) - materialLoss + noise;
        fprintf('RSSI from Router %d: %.2f dBm (Distance: %.2f m)\n', router.id, rssi(i), d);
    end
end

function aoa_estimates = estimateAoA_2DMUSIC(room)
    %Parameters
    fc = 2.4e9;             
    c = 3e8;                
    lambda = c/fc;     
    tx_power = 30;
    linear_power = 10^(tx_power/10);
    cfgWLAN = wlanNonHTConfig("ChannelBandwidth",'CBW20',...
        'MCS',4,...
        'PSDULength',1000);
    rxPos = room.robot.position';
    element = phased.CrossedDipoleAntennaElement("FrequencyRange",[2.39e9,2.41e9]);
    fprintf('AoA estimates from all routers: \n');
    for i=1:length(room.routers)
        txPos = room.routers(i).position';
        direction = rxPos-txPos;
        az = atan2d(direction(2),direction(1));
        el = asind(direction(3)/norm(direction));
        aoa_estimates = zeros(length(room.routers),2);

        normals = repmat([az; el], 1, size(room.routers(i).antennaPositions,1));
        array = phased.ConformalArray(...
            'ElementPosition',room.routers(i).antennaPositions',...
            'Element',element,...
            'ElementNormal',normals);

        fprintf('Router %d -> Expected: Azimuth = %.2f°, Elevation = %.2f°\n', ...
            i,az,el);

        wifiWaveform = sqrt(linear_power) * wlanWaveformGenerator(randi([0 1], 8*cfgWLAN.PSDULength, 1), cfgWLAN);
        collector = phased.WidebandCollector('Sensor',array,...
            'SampleRate',20e6,...
            'ModulatedInput',true,...
            'CarrierFrequency',fc,...
            'PropagationSpeed',c,...
            'NumSubbands',64);

        rx_collected = collector(wifiWaveform,[az;el]);

        noise = sqrt(0.5)*(randn(size(rx_collected)) + 1j*randn(size(rx_collected)));
        noisy_rx_collected = rx_collected + noise;

        snapshot_count = size(noisy_rx_collected, 1);
        num_antennas = size(noisy_rx_collected, 2);
        rx_sig = noisy_rx_collected.';

        R = (rx_sig * rx_sig')/snapshot_count;

        [V, D] = eig(R);
        [~, idx] = sort(diag(D), 'descend');
        V = V(:, idx);
        En = V(:, 2:end);

        el_grid = -90:0.1:0;
        az_grid = -180:0.1:180;

        elementPos = room.routers(i).antennaPositions';
        elementPos = elementPos.';

        Pmusic = zeros(length(el_grid), length(az_grid));
        k = 2 * pi/lambda;

        for aIdx = 1:length(az_grid)
            az_ang = az_grid(aIdx);
            for eIdx = 1:length(el_grid)
                el_ang = el_grid(eIdx);
                theta = deg2rad(90 - el_ang);
                phi = deg2rad(az_ang);
                u = [sin(theta)*cos(phi); sin(theta)*sin(phi); cos(theta)];
                phaseShifts = exp(1j * k * (elementPos * u));
                a = phaseShifts / norm(phaseShifts);
                Pmusic(eIdx, aIdx) = 1 / (a' * (En * En') * a);
            end
        end

        % Normalize and find peak
        Pmusic_dB = 10*log10(abs(Pmusic) / max(abs(Pmusic(:))));
        [~, maxIdx] = max(Pmusic_dB(:));
        [elEstIdx, azEstIdx] = ind2sub(size(Pmusic_dB), maxIdx);
        azEst = az_grid(azEstIdx);
        elEst = el_grid(elEstIdx);

        azEst = mod(azEst + 360, 360); % Normalize to [0, 360°)
        if azEst > 180
            azEst = azEst - 360; % Convert to [-180, 180°]
        end

        if (room.routers(i).position(1) == room.robot.position(1) & room.routers(i).position(2) == room.robot.position(2))
            azEst = 0;
        end

        aoa_estimates(i,:) = [azEst, elEst];

        fprintf('Router %d -> Result: Azimuth = %.2f°, Elevation = %.2f°\n', ...
        i, azEst, elEst);
    end
end

function est_robot_pos = robot_localization(room, rssi_cvalues, aoa_estimates)
    num_routers = length(room.routers);
    d0 = 1;   
    PL0 = 40; 
    n = 3;    

    d = zeros(num_routers, 1);
    dir_vectors = zeros(num_routers, 3);
    router_positions = zeros(num_routers, 3);

    for i = 1:num_routers
        router_positions(i, :) = room.routers(i).position;
        d(i) = d0 * 10^((room.routers(i).transmitPower - rssi_cvalues(i) - PL0) / (10 * n));
        az = deg2rad(aoa_estimates(i, 1));
        el = deg2rad(aoa_estimates(i, 2));
        dir_vectors(i, :) = [cos(el) * cos(az), cos(el) * sin(az), sin(el)];
    end

    function residuals = objective_function(robot_pos)
        residuals = zeros(2 * num_routers, 1);
        for j = 1:num_routers
            % Trilateration residual
            estimated_distance = norm(robot_pos - router_positions(j, :));
            residuals(j) = estimated_distance - d(j);

            % Triangulation residual
            vector_to_robot = robot_pos - router_positions(j, :);
            if norm(vector_to_robot) > 0
                cos_angle = dot(vector_to_robot / norm(vector_to_robot), dir_vectors(j, :));
                residuals(num_routers + j) = 0.8*(1 - cos_angle);
            else
                residuals(num_routers + j) = 1;
            end
        end
    end

    initial_guess = mean(router_positions, 1);
    options = optimoptions('lsqnonlin', 'Display', 'off', 'MaxIterations', 1000);
    est_robot_pos = lsqnonlin(@objective_function, initial_guess, [], [], options);

    error_pos = sqrt(sum(est_robot_pos - room.robot.position).^2);
    fprintf('Estimated Robot Position: [%.2f, %.2f, %.2f]\n', est_robot_pos(1), est_robot_pos(2), est_robot_pos(3));
    fprintf('Position Estimation Error: %.2f m\n', error_pos);
end

clc;
clear all;

room = createWiFiEnvironment();
rssi_values = computeRSSI(room);
aoa_estimates = estimateAoA_2DMUSIC(room);
robot_pos = robot_localization(room, rssi_values, aoa_estimates);
visualizeEnvironment(room, robot_pos);