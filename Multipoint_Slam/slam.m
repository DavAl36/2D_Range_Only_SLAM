

source "../tools/utilities/geometry_helpers_2d.m"


function v_flatten=flatten_transform(X)
	
	R=X(1:2,1:2);
	r1=R(1:2,1);
	r2=R(1:2,2);
    t=X(1:2,3);
	v_flatten=[r1; r2; t];
	
endfunction 

function v_flat=flatten_rot(R)
	
	r1=R(1:2,1);
	r2=R(1:2,2);
	v_flat=[r1;r2];

endfunction

function v_idx=poseMatrixIndex(pose_index, num_poses, num_landmarks)
  pose_dim=3;
  landmark_dim=2;

  if (pose_index>num_poses)
    v_idx=-1;
    return;
  endif;
  v_idx=1+(pose_index-1)*pose_dim;
endfunction;


function v_idx=landmarkMatrixIndex(landmark_index, num_poses, num_landmarks)
  pose_dim=3;
  landmark_dim=2;
  if (landmark_index>num_landmarks)
    v_idx=-1;
    return;
  endif;
  v_idx=1 + (num_poses)*pose_dim + (landmark_index-1) * landmark_dim;
endfunction;

function [e,Jr,Jl]=errorAndJacobianRange(Xr,Xl,z)
   
   xr = Xr(1,3); 
   yr = Xr(2,3); 
   cost = Xr(1,1); 
   sint = Xr(2,1);
   xl = Xl(1); 
   yl = Xl(2); 
   
   %without rotation
   
   %z_hat = sqrt( (xr-xl)^2 + (yr-yl)^2 ) ;
   %e=z_hat-z;
   %Jr = [ (xr-xl)/z_hat    (yr-yl)/z_hat    0];
   %Jl = [-(xr-xl)/z_hat   -(yr-yl)/z_hat     ];
   
   %with rotation   
   
   dR = [-sint  cost  
         -cost -sint];
   R=Xr(1:2,1:2);
   t=Xr(1:2,3);
   
   delta_t = Xl-t;
   r_Xl = R'*(delta_t);
   
   z_hat=norm(r_Xl);
   e=z_hat-z;
   
   
   Jr = (r_Xl'/norm(r_Xl))*[-R' dR*delta_t];
   Jl = (r_Xl'/norm(r_Xl))*[R'];

endfunction;


function [e,Jr_j,Jr_i]=errorAndJacobianOdometry(Xr_i,Xr_j,z)
	
	Rj=Xr_j(1:2,1:2);
	tj=Xr_j(1:2,3);
	Ri=Xr_i(1:2,1:2);
	ti=Xr_i(1:2,3);
	R_dot_zero=[0 -1;1 0];
	z_hat = flatten_transform(inv(Xr_i)*Xr_j);%chordal+flatten
	e = z_hat - flatten_transform(z);
  Jr_j=zeros(6,3);  
  Jr_j(1:4,1:2)=zeros(4,2);
  Jr_j(5:6,1:2)=Ri';
  Jr_j(1:4,3)=flatten_rot(Ri'*R_dot_zero*Rj);
  Jr_j(5:6,3)=Ri'*R_dot_zero*tj;  
	Jr_i=zeros(6,3);
  Jr_i = -Jr_j;
    

endfunction




function [XR, XL]=boxPlus(XR, XL, num_poses, num_landmarks, dx)
  pose_dim=3;
  landmark_dim=2;
  for(pose_index=1:num_poses)
    pose_matrix_index=poseMatrixIndex(pose_index, num_poses, num_landmarks);
    dxr=dx(pose_matrix_index:pose_matrix_index+pose_dim-1);
    XR(:,:,pose_index)=v2t(dxr)*XR(:,:,pose_index);
  endfor;
  for(landmark_index=1:num_landmarks)
    landmark_matrix_index=landmarkMatrixIndex(landmark_index, num_poses, num_landmarks);
    dxl=dx(landmark_matrix_index:landmark_matrix_index+landmark_dim-1,:);
    XL(:,landmark_index)+=dxl;
  endfor;
endfunction;


function [XR, XL, chi_stats, chi_stats2, num_inliers]=doMultiICP(XR, XL, Z_range, 
							associations_range,
              Z_odometry, 
							associations_odom, 
							num_poses, 
							num_landmarks, 
							num_iterations, 
							damping, 
							kernel_threshold,
              landmark_id_to_array_index)
   
  pose_dim=3;
  landmark_dim=2;  
  omega= 0.1*diag(6);

  
  chi_stats=zeros(1,num_iterations);
  chi_stats2=zeros(1,num_iterations);
  
  num_inliers=zeros(1,num_iterations);
  system_size=pose_dim*num_poses+landmark_dim*num_landmarks;
  
  for (iteration=1:num_iterations)
    H=zeros(system_size, system_size);
    b=zeros(system_size,1);
    chi_stats(iteration)=0;
    %RANGE PART
    
   for (measurement_num=1:size(Z_range,2))
      pose_index=associations_range(1,measurement_num);
      landmark_index=associations_range(2,measurement_num);
      landmark_index=landmark_id_to_array_index(landmark_index);
      z=Z_range(:,measurement_num);
      Xr=XR(:,:,pose_index);
      Xl=XL(:,landmark_index);
      [e,Jr,Jl] = errorAndJacobianRange(Xr, Xl, z);
      chi=e'*e;
      if (chi>kernel_threshold)
      	e*=sqrt(kernel_threshold/chi);
      	chi=kernel_threshold;
      else
      	num_inliers(iteration)++;
      endif;
      chi_stats(iteration)+=chi;
      
      Hrr=Jr'*Jr;
      Hrl=Jr'*Jl;
      Hll=Jl'*Jl;
      br=Jr'*e;
      bl=Jl'*e;

      %pose_matrix_index = pmi
      pmi=poseMatrixIndex(pose_index, num_poses, num_landmarks);
      %landmark_matrix_index = lmi
      lmi=landmarkMatrixIndex(landmark_index, num_poses, num_landmarks);

      H(pmi:pmi+pose_dim-1,
	      pmi:pmi+pose_dim-1)+=Hrr;

      H(pmi:pmi+pose_dim-1,
	      lmi:lmi+landmark_dim-1)+=Hrl;

      H(lmi:lmi+landmark_dim-1,
	      lmi:lmi+landmark_dim-1)+=Hll;

      H(lmi:lmi+landmark_dim-1,
	      pmi:pmi+pose_dim-1)+=Hrl';

      b(pmi:pmi+pose_dim-1)+=br;
      b(lmi:lmi+landmark_dim-1)+=bl;

  endfor
  %ODOMETRY PART
 
 for (odometry_measurement_num=1:size(Z_odometry,2))
		
		pose_i_index=associations_odom(1,odometry_measurement_num);
		pose_j_index=associations_odom(2,odometry_measurement_num);
		odometry = Z_odometry(:,:,odometry_measurement_num);
		Xr_i=XR(:,:,pose_i_index);
		Xr_j=XR(:,:,pose_j_index);
		
		[e_odom,Jr_j,Jr_i] = errorAndJacobianOdometry(Xr_i, Xr_j, odometry);
        chi2=e_odom'*e_odom;
        chi_stats2(iteration)+=chi2;
        
        Hri=Jr_i'*omega*Jr_i;
		Hrij=Jr_i'*omega*Jr_j;
		Hrj=Jr_j'*omega*Jr_j;
		bri=Jr_i'*omega*e_odom;
		brj=Jr_j'*omega*e_odom;
      
    %pose_matrix_index_i = pmi
    pmi=poseMatrixIndex(pose_i_index, num_poses, num_landmarks);
    %pose_matrix_index_j = lmi
    lmi=poseMatrixIndex(pose_j_index, num_poses, num_landmarks);
 
    H(pmi:pmi+pose_dim-1,
	    pmi:pmi+pose_dim-1)+=Hri;
    H(pmi:pmi+pose_dim-1,
      lmi:lmi+pose_dim-1)+=Hrij;
    H(lmi:lmi+pose_dim-1,
      lmi:lmi+pose_dim-1)+=Hrj;
    H(lmi:lmi+pose_dim-1,
      pmi:pmi+pose_dim-1)+=Hrij';
    b(pmi:pmi+pose_dim-1)+=bri;
    b(lmi:lmi+pose_dim-1)+=brj;

	endfor

  
  
  
    H+=eye(system_size)*damping;
    dx=zeros(system_size,1);

    dx(pose_dim+1:end)=-(H(pose_dim+1:end,pose_dim+1:end)\b(pose_dim+1:end,1));
    [XR, XL]=boxPlus(XR,XL,num_poses, num_landmarks, dx);

  endfor
endfunction




 