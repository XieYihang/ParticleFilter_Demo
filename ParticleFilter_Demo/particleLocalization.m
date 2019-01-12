% Robotics: Estimation and Learning 
% WEEK 4
% 
% Complete this function following the instruction. 
function myPose = particleLocalization(ranges, scanAngles, map, param)

% Number of poses to calculate
N = size(ranges, 2);
% Output format is [x1 x2, ...; y1, y2, ...; z1, z2, ...]
myPose = zeros(3, N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Map Parameters 
% 
% % the number of grids for 1 meter.
myResolution = param.resol;
% % the origin of the map in pixels
myOrigin = param.origin; 

% The initial pose is given
myPose(:,1) = param.init_pose;
% You should put the given initial pose into myPose for j=1, ignoring the j=1 ranges. 
% The pose(:,1) should be the pose when ranges(:,j) were measured.

% Decide the number of particles, M.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M = 100;                          % Please decide a reasonable number of M, 
                               % based on your experiment using the practice data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create M number of particles
P = repmat(myPose(:,1), [1, M]);
W = ones(1, M) / M;
for j = 2:N % You will start estimating myPose from j=2 using ranges(:,2).
   
    % 1) Propagate the particles 
    Sigma = [0.1;0.1;0.04];
    correlation=[];
    for k = 1:M
%         P(:, k) = P(:, k) + Sigma .* randn(3, 1);
%         randomCita=currentCita+0.25*randn(1,1);
        %dist=0.0125+0.001*randn(1,1);
%         dist=0.0125;
%         P(:,k)=P(:,k)+dist*[cos(randomCita) -sin(randomCita) 0]'; 

        dist = 0.025 + randn(1,1) * 0.1;
        heading = P(3, k) + randn(1, 1) * 0.15; %heading is the sum of the previous one plus a random value
        P(1, k) = P(1, k) + dist * cos(heading);
        P(2, k) = P(2, k) - dist * sin(heading);
        P(3, k) = heading;
        
    end
      
    % 2) Measurement Update 
    %   2-1) Find grid cells hit by the rays (in the grid map coordinate frame)    
    for k = 1:M
        
        i_pos_x = ceil(myResolution * P(1, k)) + myOrigin(1);
        i_pos_y = ceil(myResolution * P(2, k)) + myOrigin(2);
        
        ids = randsample(size(ranges, 1), 100);
        phi = P(3, k) + scanAngles(ids);
        
        x_occ = ranges(ids, j) .* cos(phi) + P(1, k);
        y_occ = -ranges(ids, j) .* sin(phi) + P(2, k);
        
        i_x_occ = ceil(myResolution * x_occ) + myOrigin(1);
        i_y_occ = ceil(myResolution * y_occ) + myOrigin(2);
         N_occ = size(i_x_occ, 1);
        corr = 0;
        for n = 1:N_occ
            if (i_x_occ(n) > size(map, 2)) || (i_x_occ(n) < 1)
                continue;
            end
            if (i_y_occ(n) > size(map, 1)) || (i_y_occ(n) < 1)
                continue;
            end
            
            if map(i_y_occ(n), i_x_occ(n)) <= 0
                %corr = corr - 5 * map(i_y_occ(k), i_x_occ(k));
                corr = corr + map(i_y_occ(n), i_x_occ(n)) + 0.5;
            end
            if map(i_y_occ(n), i_x_occ(n)) >= 1
                %corr = corr + 10 * map(i_y_occ(k), i_x_occ(k));
                corr = corr + map(i_y_occ(n), i_x_occ(n)) + 0.5;
            end
         
        end
        W(k) = W(k) * corr; 
    end
    W_norm = W / sum(W);
    N_eff = floor( sum(W_norm) / sum(W_norm .^ 2) );
    %   2-2) For each particle, calculate the correlation scores of the particles
    [W_sorted, sort_ids] = sort(W_norm, 'descend');
    W_sample = W_sorted(1:N_eff); 
    P_sorted = P(:, sort_ids);
    P_sample = P_sorted(:, 1:N_eff); % 
    %   2-3) Update the particle weights         
    p_ids = randsample(N_eff, M, 'true', W_sample); % 
    W = W_sample(p_ids');
    P = P_sample(:, p_ids');
    %   2-4) Choose the best particle to update the pose
     myPose(:, j) = mean(P(:, 1:10), 2);
    % 3) Resample if the effective number of particles is smaller than a threshold

    % 4) Visualize the pose on the map as needed
   

end

end

