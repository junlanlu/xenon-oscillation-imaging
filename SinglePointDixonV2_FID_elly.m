function [shifted_FID, Delta_angle_deg] = SinglePointDixonV2_FID_elly(DisFID,RbcBarrierRatio,GasFID)

% Calculate Parameters from User Defined Parameters
Desired_angle_rad=atan2(1/RbcBarrierRatio,1);
Desired_angle_deg=rad2deg(Desired_angle_rad);

%% Calculate Initial Phase
Initial_angle_rad=angle(sum(DisFID(1,:)));
Initial_angle_deg=rad2deg(Initial_angle_rad);

%% Calculate Delta Angle
Delta_angle_rad=Desired_angle_rad-Initial_angle_rad;
Delta_angle_deg=-(Desired_angle_deg-Initial_angle_deg);

%% Global Phase Shift to Align 
shifted_FID = DisFID.*exp(1i*Delta_angle_rad);
