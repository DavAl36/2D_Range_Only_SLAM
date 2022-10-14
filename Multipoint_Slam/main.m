clc 
clear
warning off;

addpath '../'
addpath '../tools/g2o_wrapper'
source "../tools/utilities/geometry_helpers_2d.m"
addpath "./Multipoint_Slam"


disp("\n\nExtracting data from ground truth and initial guess\n")

#----------------GROUND TRUTH-----------------------
[lan_truth, poses, transitions_truth, observations] = loadG2o('datasets/slam2d_range_only_ground_truth.g2o');
poses_length = size(poses,2);
landmark_length = size(lan_truth,2);
transitions_length = size(transitions_truth,2);
observations_length = size(observations,2);
complete_poses_truth = zeros(3,1,poses_length);
complete_landmark_truth = zeros(2,1,landmark_length);
complete_observations_truth = zeros(1,1,observations_length);
for (n=1:poses_length)
    complete_poses_truth(1,1,n) = poses(n).x;
    complete_poses_truth(2,1,n) = poses(n).y;
    complete_poses_truth(3,1,n) = poses(n).theta;
endfor;
for (n=1:landmark_length)
    complete_landmark_truth(1,1,n) = lan_truth(n).x_pose;
    complete_landmark_truth(2,1,n) = lan_truth(n).y_pose;
endfor;
for (n=1:observations_length)
    complete_observations_truth(1,1,n) =observations(n).observation.bearing;
endfor;
#----------------INITIAL GUESS-----------------------
[_, poses, transitions, observations] = loadG2o("datasets/slam2d_range_only_initial_guess.g2o");
poses_length = size(poses,2);
observations_length = size(observations,2);
transitions_length = size(transitions,2);
land_list = [ ];
for (n=1:observations_length)
    meas_num = size(observations(n).observation,2);
    for(c=1:meas_num)
        value = observations(n).observation(c).id; 
        if (any (land_list == value) == 0) 
           land_list = [land_list ; value];
        endif
    endfor;
endfor;
landmark_length = size(land_list,1);
land_list = sort(land_list);
complete_poses = zeros(3,1,poses_length);
complete_landmark = zeros(2,1,landmark_length);
complete_observations = zeros(1,1,observations_length);
for (n=1:poses_length)
    complete_poses(1,1,n) = poses(n).x;
    complete_poses(2,1,n) = poses(n).y;
    complete_poses(3,1,n) = poses(n).theta;
endfor;
for (n=1:observations_length)
    complete_observations(1,1,n) =observations(n).observation.bearing;
endfor;

#----------------FUNCTIONS-----------------------

#take in input landmark ID and return a list of robot poses observers
function [XY,id_poses_of_landmark]=whoSeeLandmark(land_id,poses,observations,observations_length ,poses_length)
  XY = [];
  id_poses_of_landmark = [];
  for (n=1:observations_length) 
      pose_id =  observations(n).pose_id ; 
      meas_num = size(observations(n).observation,2); 
      for(c=1:meas_num) 
          value = observations(n).observation(c).id; 
          range = observations(n).observation(c).bearing; 
          if (value == land_id) 
             id_poses_of_landmark = [id_poses_of_landmark ; [pose_id range value] ];
             for (z=1:poses_length)  
                 if(poses(z).id == pose_id)
                     XY = [XY ; [ poses(z).x poses(z).y  poses(z).theta]];
                 endif;
             endfor; 
          endif;
      endfor;
  endfor;
endfunction;

function res = landmark_positions(poses_robot,range_measurements,id_land,n_p)
	first_x = poses_robot(1,1);
	first_y = poses_robot(2,1);
	first_r = range_measurements(1);
  identity=eye(2)*0.001;
	A = zeros(n_p-1,2);
	b = zeros(n_p-1,1);
	for(i=2:n_p)
		x_n = poses_robot(1,i);
		y_n = poses_robot(2,i);
		r_n = range_measurements(i);
		A(i,1) = first_x - x_n;
		A(i,2) = first_y - y_n;
		b(i) = 0.5 * (-first_r^2 + r_n^2 + first_x^2 - x_n^2 + first_y^2 - y_n^2 );
	end
  pseudoinverse = inv(A'*A+identity)*A'; %slide "Calibration" Rob_2 full column rank
	res = pseudoinverse*b;%b/A
end





%----------------SETTINGS-----------------
num_iterations = 65;
damping = 0.01;
kernel_threshold = 1;
num_landmarks = size(complete_landmark,3);
num_poses = size(poses,2);
associations_range = [];
Z_range = [];
XR_guess=zeros(3,3,num_poses);
XL_guess = zeros(2,1,num_landmarks);
id_threshold = min([poses.id])-1; %pose robot 1100 become 1 and so on... 



Z_odometry = zeros(3,3,transitions_length);
associations_odom = zeros(2,transitions_length);
for (n=1:transitions_length)
	associations_odom(1,n) = transitions(n).id_from - id_threshold;
	associations_odom(2,n) = transitions(n).id_to - id_threshold;
	displacement = transitions(n).v;
	odometry_transformation = v2t(displacement);
	Z_odometry(:,:,n) = odometry_transformation;
endfor





%--------------COMPUTE INITIAL GUESS LANDMARKS--------------
disp("\nCompute landmarks positions\n")
for (n=1:landmark_length) 
    [A,B] = whoSeeLandmark(land_list(n),poses,observations,observations_length ,poses_length);
    B(:,1) = B(:,1)-id_threshold;
    setPoints = [B A];
    associations_range = [ associations_range , [B(:,1)' ; B(:,3)']];
    Z_range = [Z_range , B(:,2)' ];
    howManyPosesSeenLandmark = size(B,1);
    if(howManyPosesSeenLandmark == 1 )%only one pose see the landmark
        X = setPoints(1,4); 
        Y = setPoints(1,5); 
    elseif(howManyPosesSeenLandmark >= 2 && howManyPosesSeenLandmark <= 3)%take nearest pose
         [minval, idx] = min(setPoints(:,2)', [], 2);%take Id of the smallest range
         X = setPoints(idx,4);
         Y =setPoints(idx,5); 
    else
        sol = landmark_positions(A(:,1:2)',B(:,2)',B(1,3),howManyPosesSeenLandmark);
        X = sol(1);
        Y = sol(2);
    endif
    XL_guess(:,:,n) =  [X Y]';
    %printf("Land: %d Coordinates: %f , %f N_Poses:%d \n",land_list(n), X,Y,howManyPosesSeenLandmark)
endfor
%transform input pose x,y,theta in a 3x3 matrix ([R|t])
for(i = 1:num_poses)
   ang = poses(i).theta;
   R = [cos(ang) -sin(ang);
        sin(ang)  cos(ang)];

   t = [poses(i).x poses(i).y];
   T = [ R(1,1) R(1,2) t(1);
         R(2,1) R(2,2) t(2);
            0     0     1  ];
   XR_guess(:,:,i) = T;
endfor
%Create an array useful in the association part
landmark_id_to_array_index = zeros(1, max(land_list));
for(i = 1:num_landmarks)
  for(ii = 1:max(land_list))
      if(land_list(i) == ii)
         landmark_id_to_array_index(ii) = i;
      endif
  endfor
endfor



pert_deviation=0.0001; 
pert_scale=eye(3)*pert_deviation;
for (pose_num=2:num_poses) 
    xr=rand(3,1)-0.5;
    dXr=v2t(pert_scale*xr);
    XR_guess(:,:,pose_num)=dXr*XR_guess(:,:,pose_num);
endfor;
dXl=(rand(2,1,num_landmarks)-0.5)*pert_deviation;
XL_guess+=dXl;



%plot before ICP
figure(1, 'position',[0,0,800,800]);
subplot(2,2,1);
plot(complete_poses_truth(1,:,:),complete_poses_truth(2,:,:),'ro',"linewidth",2)
hold on;
plot(complete_landmark_truth(1,:,:),complete_landmark_truth(2,:,:),'bo',"linewidth",2)
title ("[TRUTH] Poses+Landmarks (x,y)");
hold on;
subplot(2,2,2);
plot( XR_guess(1,3,:), XR_guess(2,3,:),'ro',"linewidth",2)
hold on; 
subplot(2,2,2);
plot( XL_guess(1,:,:), XL_guess(2,:,:),'bo',"linewidth",2)
title ("[GUESS] Poses+Landmarks (x,y)");
hold on;

subplot(2,2,3);
plot(complete_poses_truth(1,:,:),complete_poses_truth(2,:,:),'ro',"linewidth",2)
hold on;
plot( XR_guess(1,3,:), XR_guess(2,3,:),'bo',"linewidth",2)
title ("[TRUTH/GUESS] Poses (x,y)");
hold on; 
subplot(2,2,4);
plot(complete_landmark_truth(1,:,:),complete_landmark_truth(2,:,:),'ro',"linewidth",2)
hold on;
plot( XL_guess(1,:,:), XL_guess(2,:,:),'bo',"linewidth",2)
title ("[TRUTH/GUESS] Landmark (x,y)");
hold on;

%--------------------------ICP-----------------------------
source "./slam.m"

disp("\nCompute landmarks positions with ICP\n")

[XR, XL, chi_stats1,chi_stats2, num_inliers]=doMultiICP(XR_guess,XL_guess , 
              Z_range,  
					  	associations_range,
              Z_odometry, 
							associations_odom,
              num_poses, 
							num_landmarks, 
							num_iterations, 
							damping, 
							kernel_threshold,
              landmark_id_to_array_index);
  
%printf("Chi stats1: %f\n",chi_stats1)   
%printf("Chi stats2: %f\n",chi_stats2)   

%plot after ICP

figure(3, 'position',[0,0,800,800]);

subplot(2,2,1);
plot(complete_poses_truth(1,:,:),complete_poses_truth(2,:,:),'g*',"linewidth",2)
hold on;
plot(XR(1,3,:),XR(2,3,:),'b*',"linewidth",2)
title ("Poses (x,y)");
hold on; 


subplot(2,2,2);
plot(complete_landmark_truth(1,:,:),complete_landmark_truth(2,:,:),'g*',"linewidth",2)
hold on;
plot(XL(1,:,:),XL(2,:,:),'b*',"linewidth",2)
title ("Landmarks (x,y)");
legend("ground truth","estimate")

hold on; 

subplot(2,2,3);
plot(chi_stats1,"linewidth",4)
title ("chi stats range");
hold on;

subplot(2,2,4);
plot(chi_stats2,"linewidth",4)
title ("chi stats odometry");
hold on;
pause (10);
