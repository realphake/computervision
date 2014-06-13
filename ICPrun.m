% this file is used to gather data for the first assignment

%[out_est_uni100, R_est_uni100, T_est_uni100, iter_est_uni100] = merge_pointclouds([0:99], 100, 'uniform', 'estimate', 400);

%[out_est_ran100, R_est_ran100, T_est_ran100, iter_est_ran100] = merge_pointclouds([0:99], 100, 'random', 'estimate', 400);

%[out_est_nor100, R_est_nor100, T_est_nor100, iter_est_nor100] = merge_pointclouds([0:99], 100/20, 'normals', 'estimate', 400);

[out_est_uni500, R_est_uni500, T_est_uni500, iter_est_uni500] = merge_pointclouds([0:99], 500, 'uniform', 'estimate', 400);

[out_est_ran500, R_est_ran500, T_est_ran500, iter_est_ran500] = merge_pointclouds([0:99], 500, 'random', 'estimate', 400);

[out_est_nor500, R_est_nor500, T_est_nor500, iter_est_nor500] = merge_pointclouds([0:99], 500/20, 'normals', 'estimate', 400);

[out_est_uni2000, R_est_uni2000, T_est_uni2000, iter_est_uni2000] = merge_pointclouds([0:99], 2000, 'uniform', 'estimate', 400);

[out_est_ran2000, R_est_ran2000, T_est_ran2000, iter_est_ran2000] = merge_pointclouds([0:99], 2000, 'random', 'estimate', 400);

[out_est_nor2000, R_est_nor2000, T_est_nor2000, iter_est_nor2000] = merge_pointclouds([0:99], 2000/20, 'normals', 'estimate', 400);

[out_est_uni5000, R_est_uni5000, T_est_uni5000, iter_est_uni5000] = merge_pointclouds([0:99], 5000, 'uniform', 'estimate', 400);

[out_est_ran5000, R_est_ran5000, T_est_ran5000, iter_est_ran5000] = merge_pointclouds([0:99], 5000, 'random', 'estimate', 400);

[out_est_nor5000, R_est_nor5000, T_est_nor5000, iter_est_nor5000] = merge_pointclouds([0:99], 5000/20, 'normals', 'estimate', 400);
