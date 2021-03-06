% (1) All 'coi'.
coi = {'sequencing_end_time (seconds)', 'sequencing_lifetime (seconds)', 'lifetime_after_tag_flow (seconds)', 'len_starting_single_pore_run (#)', ...
       'level_call_transition_rate (1/seconds)', 'ttt_median (milliseconds)', 'ttt_rate (1/milliseconds)', ...
       'align_copies (#)', ...
       'align_align_length (#)', 'align_num_delete (#)', 'align_num_ident (#)', 'align_num_insert (#)', 'align_num_mismatch (#)', 'align_procession_length (#)', 'align_read_length (#)', ...
       'align_homo_align_length (#)', 'align_homo_num_delete (#)', 'align_homo_num_ident (#)', 'align_homo_num_insert (#)', 'align_homo_num_mismatch (#)', 'align_homo_procession_length (#)', 'align_homo_read_length (#)', ...
       'level_call_A_counts (#)', 'level_call_A_k_cat_rate (1/seconds)', 'level_call_A_k_off_rate (1/seconds)', 'level_call_A_lower_bound (OC fraction)', 'level_call_A_mean_dwell_time (seconds)', 'level_call_A_ttt_median (milliseconds)', 'level_call_A_ttt_rate (1/milliseconds)', 'level_call_A_upper_bound (OC fraction)', ...
       'level_call_C_counts (#)', 'level_call_C_k_cat_rate (1/seconds)', 'level_call_C_k_off_rate (1/seconds)', 'level_call_C_lower_bound (OC fraction)', 'level_call_C_mean_dwell_time (seconds)', 'level_call_C_ttt_median (milliseconds)', 'level_call_C_ttt_rate (1/milliseconds)', 'level_call_C_upper_bound (OC fraction)', ...
       'level_call_G_counts (#)', 'level_call_G_k_cat_rate (1/seconds)', 'level_call_G_k_off_rate (1/seconds)', 'level_call_G_lower_bound (OC fraction)', 'level_call_G_mean_dwell_time (seconds)', 'level_call_G_ttt_median (milliseconds)', 'level_call_G_ttt_rate (1/milliseconds)', 'level_call_G_upper_bound (OC fraction)', ...
       'level_call_T_counts (#)', 'level_call_T_k_cat_rate (1/seconds)', 'level_call_T_k_off_rate (1/seconds)', 'level_call_T_lower_bound (OC fraction)', 'level_call_T_mean_dwell_time (seconds)', 'level_call_T_ttt_median (milliseconds)', 'level_call_T_ttt_rate (1/milliseconds)', 'level_call_T_upper_bound (OC fraction)', ...
       'level_call_super_dwell_waiting_time_median'};

% (2) No time 'coi'.
coi = {'level_call_transition_rate (1/seconds)', 'ttt_median (milliseconds)', 'ttt_rate (1/milliseconds)', ...
       'align_copies (#)', ...
       'align_align_length (#)', 'align_num_delete (#)', 'align_num_ident (#)', 'align_num_insert (#)', 'align_num_mismatch (#)', 'align_procession_length (#)', 'align_read_length (#)', ...
       'align_homo_align_length (#)', 'align_homo_num_delete (#)', 'align_homo_num_ident (#)', 'align_homo_num_insert (#)', 'align_homo_num_mismatch (#)', 'align_homo_procession_length (#)', 'align_homo_read_length (#)', ...
       'level_call_A_counts (#)', 'level_call_A_k_cat_rate (1/seconds)', 'level_call_A_k_off_rate (1/seconds)', 'level_call_A_lower_bound (OC fraction)', 'level_call_A_mean_dwell_time (seconds)', 'level_call_A_ttt_median (milliseconds)', 'level_call_A_ttt_rate (1/milliseconds)', 'level_call_A_upper_bound (OC fraction)', ...
       'level_call_C_counts (#)', 'level_call_C_k_cat_rate (1/seconds)', 'level_call_C_k_off_rate (1/seconds)', 'level_call_C_lower_bound (OC fraction)', 'level_call_C_mean_dwell_time (seconds)', 'level_call_C_ttt_median (milliseconds)', 'level_call_C_ttt_rate (1/milliseconds)', 'level_call_C_upper_bound (OC fraction)', ...
       'level_call_G_counts (#)', 'level_call_G_k_cat_rate (1/seconds)', 'level_call_G_k_off_rate (1/seconds)', 'level_call_G_lower_bound (OC fraction)', 'level_call_G_mean_dwell_time (seconds)', 'level_call_G_ttt_median (milliseconds)', 'level_call_G_ttt_rate (1/milliseconds)', 'level_call_G_upper_bound (OC fraction)', ...
       'level_call_T_counts (#)', 'level_call_T_k_cat_rate (1/seconds)', 'level_call_T_k_off_rate (1/seconds)', 'level_call_T_lower_bound (OC fraction)', 'level_call_T_mean_dwell_time (seconds)', 'level_call_T_ttt_median (milliseconds)', 'level_call_T_ttt_rate (1/milliseconds)', 'level_call_T_upper_bound (OC fraction)', ...
       'level_call_super_dwell_waiting_time_median'};

% (3) No base counts 'coi'.   
coi = {'level_call_transition_rate (1/seconds)', 'ttt_median (milliseconds)', 'ttt_rate (1/milliseconds)', ...
       'align_copies (#)', ...
       'align_align_length (#)', 'align_num_delete (#)', 'align_num_ident (#)', 'align_num_insert (#)', 'align_num_mismatch (#)', 'align_procession_length (#)', 'align_read_length (#)', ...
       'align_homo_align_length (#)', 'align_homo_num_delete (#)', 'align_homo_num_ident (#)', 'align_homo_num_insert (#)', 'align_homo_num_mismatch (#)', 'align_homo_procession_length (#)', 'align_homo_read_length (#)', ...
       'level_call_A_k_cat_rate (1/seconds)', 'level_call_A_k_off_rate (1/seconds)', 'level_call_A_lower_bound (OC fraction)', 'level_call_A_mean_dwell_time (seconds)', 'level_call_A_ttt_median (milliseconds)', 'level_call_A_ttt_rate (1/milliseconds)', 'level_call_A_upper_bound (OC fraction)', ...
       'level_call_C_k_cat_rate (1/seconds)', 'level_call_C_k_off_rate (1/seconds)', 'level_call_C_lower_bound (OC fraction)', 'level_call_C_mean_dwell_time (seconds)', 'level_call_C_ttt_median (milliseconds)', 'level_call_C_ttt_rate (1/milliseconds)', 'level_call_C_upper_bound (OC fraction)', ...
       'level_call_G_k_cat_rate (1/seconds)', 'level_call_G_k_off_rate (1/seconds)', 'level_call_G_lower_bound (OC fraction)', 'level_call_G_mean_dwell_time (seconds)', 'level_call_G_ttt_median (milliseconds)', 'level_call_G_ttt_rate (1/milliseconds)', 'level_call_G_upper_bound (OC fraction)', ...
       'level_call_T_k_cat_rate (1/seconds)', 'level_call_T_k_off_rate (1/seconds)', 'level_call_T_lower_bound (OC fraction)', 'level_call_T_mean_dwell_time (seconds)', 'level_call_T_ttt_median (milliseconds)', 'level_call_T_ttt_rate (1/milliseconds)', 'level_call_T_upper_bound (OC fraction)', ...
       'level_call_super_dwell_waiting_time_median'};   
   
% (4) No alignment-related 'coi'.   
coi = {'level_call_transition_rate (1/seconds)', 'ttt_median (milliseconds)', 'ttt_rate (1/milliseconds)', ...
       'level_call_A_k_cat_rate (1/seconds)', 'level_call_A_k_off_rate (1/seconds)', 'level_call_A_lower_bound (OC fraction)', 'level_call_A_mean_dwell_time (seconds)', 'level_call_A_ttt_median (milliseconds)', 'level_call_A_ttt_rate (1/milliseconds)', 'level_call_A_upper_bound (OC fraction)', ...
       'level_call_C_k_cat_rate (1/seconds)', 'level_call_C_k_off_rate (1/seconds)', 'level_call_C_lower_bound (OC fraction)', 'level_call_C_mean_dwell_time (seconds)', 'level_call_C_ttt_median (milliseconds)', 'level_call_C_ttt_rate (1/milliseconds)', 'level_call_C_upper_bound (OC fraction)', ...
       'level_call_G_k_cat_rate (1/seconds)', 'level_call_G_k_off_rate (1/seconds)', 'level_call_G_lower_bound (OC fraction)', 'level_call_G_mean_dwell_time (seconds)', 'level_call_G_ttt_median (milliseconds)', 'level_call_G_ttt_rate (1/milliseconds)', 'level_call_G_upper_bound (OC fraction)', ...
       'level_call_T_k_cat_rate (1/seconds)', 'level_call_T_k_off_rate (1/seconds)', 'level_call_T_lower_bound (OC fraction)', 'level_call_T_mean_dwell_time (seconds)', 'level_call_T_ttt_median (milliseconds)', 'level_call_T_ttt_rate (1/milliseconds)', 'level_call_T_upper_bound (OC fraction)', ...
       'level_call_super_dwell_waiting_time_median'};   
   
% (5) No cumulative 4 base 'coi'.   
coi = {'level_call_A_k_cat_rate (1/seconds)', 'level_call_A_k_off_rate (1/seconds)', 'level_call_A_lower_bound (OC fraction)', 'level_call_A_mean_dwell_time (seconds)', 'level_call_A_ttt_median (milliseconds)', 'level_call_A_ttt_rate (1/milliseconds)', 'level_call_A_upper_bound (OC fraction)', ...
       'level_call_C_k_cat_rate (1/seconds)', 'level_call_C_k_off_rate (1/seconds)', 'level_call_C_lower_bound (OC fraction)', 'level_call_C_mean_dwell_time (seconds)', 'level_call_C_ttt_median (milliseconds)', 'level_call_C_ttt_rate (1/milliseconds)', 'level_call_C_upper_bound (OC fraction)', ...
       'level_call_G_k_cat_rate (1/seconds)', 'level_call_G_k_off_rate (1/seconds)', 'level_call_G_lower_bound (OC fraction)', 'level_call_G_mean_dwell_time (seconds)', 'level_call_G_ttt_median (milliseconds)', 'level_call_G_ttt_rate (1/milliseconds)', 'level_call_G_upper_bound (OC fraction)', ...
       'level_call_T_k_cat_rate (1/seconds)', 'level_call_T_k_off_rate (1/seconds)', 'level_call_T_lower_bound (OC fraction)', 'level_call_T_mean_dwell_time (seconds)', 'level_call_T_ttt_median (milliseconds)', 'level_call_T_ttt_rate (1/milliseconds)', 'level_call_T_upper_bound (OC fraction)'};    
   
% (6) No upper/lower bound 'coi'.   
coi = {'level_call_A_k_cat_rate (1/seconds)', 'level_call_A_k_off_rate (1/seconds)', 'level_call_A_mean_dwell_time (seconds)', 'level_call_A_ttt_median (milliseconds)', 'level_call_A_ttt_rate (1/milliseconds)', ...
       'level_call_C_k_cat_rate (1/seconds)', 'level_call_C_k_off_rate (1/seconds)', 'level_call_C_mean_dwell_time (seconds)', 'level_call_C_ttt_median (milliseconds)', 'level_call_C_ttt_rate (1/milliseconds)', ...
       'level_call_G_k_cat_rate (1/seconds)', 'level_call_G_k_off_rate (1/seconds)', 'level_call_G_mean_dwell_time (seconds)', 'level_call_G_ttt_median (milliseconds)', 'level_call_G_ttt_rate (1/milliseconds)', ...
       'level_call_T_k_cat_rate (1/seconds)', 'level_call_T_k_off_rate (1/seconds)', 'level_call_T_mean_dwell_time (seconds)', 'level_call_T_ttt_median (milliseconds)', 'level_call_T_ttt_rate (1/milliseconds)'}; 
 
% (7) No ttt 'coi'.   
coi = {'level_call_A_k_cat_rate (1/seconds)', 'level_call_A_k_off_rate (1/seconds)', 'level_call_A_mean_dwell_time (seconds)', ...
       'level_call_C_k_cat_rate (1/seconds)', 'level_call_C_k_off_rate (1/seconds)', 'level_call_C_mean_dwell_time (seconds)', ...
       'level_call_G_k_cat_rate (1/seconds)', 'level_call_G_k_off_rate (1/seconds)', 'level_call_G_mean_dwell_time (seconds)', ...
       'level_call_T_k_cat_rate (1/seconds)', 'level_call_T_k_off_rate (1/seconds)', 'level_call_T_mean_dwell_time (seconds)'}; 

% (8) No k_off 'coi'.   
coi = {'level_call_A_k_cat_rate (1/seconds)', 'level_call_A_mean_dwell_time (seconds)', ...
       'level_call_C_k_cat_rate (1/seconds)', 'level_call_C_mean_dwell_time (seconds)', ...
       'level_call_G_k_cat_rate (1/seconds)', 'level_call_G_mean_dwell_time (seconds)', ...
       'level_call_T_k_cat_rate (1/seconds)', 'level_call_T_mean_dwell_time (seconds)'}; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% (9) K_cat only.   
coi = {'level_call_A_k_cat_rate (1/seconds)', ...
       'level_call_C_k_cat_rate (1/seconds)', ...
       'level_call_G_k_cat_rate (1/seconds)', ...
       'level_call_T_k_cat_rate (1/seconds)'}; 

% (10) K-off only. 
coi = {'level_call_A_k_off_rate (1/seconds)', ...
       'level_call_C_k_off_rate (1/seconds)', ...
       'level_call_G_k_off_rate (1/seconds)', ...
       'level_call_T_k_off_rate (1/seconds)'}; 
   
% (11) Dwell time only.   
coi = {'level_call_A_mean_dwell_time (seconds)', ...
       'level_call_C_mean_dwell_time (seconds)', ...
       'level_call_G_mean_dwell_time (seconds)', ...
       'level_call_T_mean_dwell_time (seconds)'}; 
   
% (12) Alignment-related 'coi' only.   
coi = {'align_copies (#)', ...
       'align_align_length (#)', 'align_num_delete (#)', 'align_num_ident (#)', 'align_num_insert (#)', 'align_num_mismatch (#)', 'align_procession_length (#)', 'align_read_length (#)', ...
       'align_homo_align_length (#)', 'align_homo_num_delete (#)', 'align_homo_num_ident (#)', 'align_homo_num_insert (#)', 'align_homo_num_mismatch (#)', 'align_homo_procession_length (#)', 'align_homo_read_length (#)'}; 
   
% (13) Cumulative 4 base 'coi' only.
coi = {'level_call_transition_rate (1/seconds)', 'ttt_median (milliseconds)', 'ttt_rate (1/milliseconds)', ...
       'level_call_super_dwell_waiting_time_median'}; 

% (14) K_cat and dwell time 'coi' only.
coi = {'level_call_A_k_cat_rate (1/seconds)', ...
       'level_call_C_k_cat_rate (1/seconds)', ...
       'level_call_G_k_cat_rate (1/seconds)', ...
       'level_call_T_k_cat_rate (1/seconds)', ...
       'level_call_A_mean_dwell_time (seconds)', ...
       'level_call_C_mean_dwell_time (seconds)', ...
       'level_call_G_mean_dwell_time (seconds)', ...
       'level_call_T_mean_dwell_time (seconds)'}; 
