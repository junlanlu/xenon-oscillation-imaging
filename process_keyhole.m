
    filelist_keyhole;

for k = 3%:22%38:length(all_keyhole)
    close all
    clear dissolved_data, clear gas_data
%     [dissolved_data, gas_data,  dwell_time] = prep_data_for_keyhole(all_keyhole_GE{1,k});
    [dissolved_data, gas_data,  TR] = prep_data_for_keyhole(all_keyhole{1,k});
    load(all_keyhole{4,k})
    Subject = all_keyhole{2,k};
   
    ImSize = 128;
    % Use only the first 1/2 of FID (higher SNR?)
%     dissolved_data(65:128,:) = [];
%     gas_data(65:128,:) = [];
%     traj(:,65:128,:) = [];
%     subjName = [all_keyhole{2,k},''];

    %% Load Proton Mask and RBC to Barrier
    if contains(Subject,'_Pre')
        proton_mask_loc = ['Y:\01_ClinicalOutput\',Subject(1:7),'\Gas_Exchange\',erase(Subject(1:7),'-'),'_pre.mat'];
    elseif contains(Subject,'_Post')
        proton_mask_loc = ['Y:\01_ClinicalOutput\',Subject(1:7),'\Gas_Exchange\',erase(Subject(1:7),'-'),'_post.mat'];
    elseif contains(Subject,'_3hr')
        proton_mask_loc = ['Y:\01_ClinicalOutput\',Subject(1:7),'\Gas_Exchange\',erase(Subject(1:7),'-'),'_final.mat'];
    elseif contains(Subject,'_')
        tempSubject = extractBefore(Subject,'_');
        proton_mask_loc = ['Y:\01_ClinicalOutput\',tempSubject,'\Gas_Exchange\',erase(tempSubject,'-'),'_highBW.mat'];
    else 
        tempSubject = Subject;
        proton_mask_loc = ['Y:\01_ClinicalOutput\',tempSubject,'\Gas_Exchange\',erase(tempSubject,'-'),'_highBW.mat'];
    end 

    if ~exist(proton_mask_loc,'file') && contains(proton_mask_loc,'_highBW')
        proton_mask_loc = erase(proton_mask_loc,'_highBW');
    end 

    load(proton_mask_loc,'mask_reg')
    protonMask = permute((mask_reg),[1, 3, 2]);

    if ImSize == 64
        protonMask = imresize3(double(protonMask),[ImSize ImSize ImSize],'nearest');
        protonMask(protonMask <= 0.75) = 0; protonMask(protonMask > 0.75) = 1; protonMask = logical(protonMask);
    end 
    % protonMask = High_Res_Gas_Mask;
    
   %%
    
    dissolved_phase_detrend_keyhole_analysis_elly(traj,gas_data,...
        dissolved_data,TR,all_keyhole{3,k},'Siemens',Subject,...
        'D:\OneDrive\Documents\Duke\CIVM Research\keyhole\',protonMask)
    
%         dissolved_phase_detrend_keyhole_analysis_elly(traj,gas_data,...
%         dissolved_data,0.015,all_keyhole{3,k},'GE',subjName,...
%         'D:\OneDrive\Documents\Duke\CIVM Research\keyhole\')
end 