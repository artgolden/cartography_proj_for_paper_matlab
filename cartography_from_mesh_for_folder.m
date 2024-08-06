% title  : GeodesicTissueAnalysisPipeline
% date   : July 2024
% author : Zoë Lange (zlange@fias.uni-frankfurt.de) modified by Artemiy
% Golden (artemiy.golden@physikalischebiologie.de)
% Geodesic projections, segmentation, Voronoi diagramms and 
% inference of forces. Morphological analysis of gastrulation
% in the extra-embryonic serosa tissue in Tribolium castaneum.
% Geodesic projections after Heemskerk & Streichan 2015.
% Segmentation using StarDist.
% Baysian inference of forces after Ishihara & Sugimura 2012.

clear; close all; clc;

%% Start Imsane to create geodesic projections and add all relevant paths
restoredefaultpath;
run("imsaneV1.2/setup.m");

%% SerosaDynamics_DS0008__514__20_02_11 single file test


parameters.timepoint       = 1;
parameters.stackSize       = [481 985 505 1];
parameters.createMesh      = "no"; % yes or no
parameters.performClosure  = "no"; % yes or no (yes for after gastrulation)
parameters.list_of_charts  = ["cylinder1"]; %, "cylinder2" , "equidistant_1", "equidistant_2"
parameters.semseg.beadradius            =  100; % size in pixels
parameters.semseg.beadthresh            =  0.3; % relative intensity between 0 and 1 (~ 0.7)
parameters.semseg.BWthreshMan           =  0.5; % additive (can be zero)
parameters.semseg.alpha                 =   70; % inf is convex hull, too small makes holes (~ 245)
parameters.semseg.alphaRegionThreshold  = 1000; % max volume of regions to suppress (~ bead volume)
parameters.layers.number                =   61; %41
parameters.layers.distance              =    1; %in pixels
parameters.layers.ref_layer_cylinder    =   13; %23
parameters.layers.ref_layer_equidistant =   11; %21

%% Kraemer2024A_DS0008TP0005DR_D1_CH0001PL_NS control volume

% tifFileName = "fused_tp_233_ch_0.tif";
% tifDir         = custom_fullfile("F:\Mariia\SerosaDynamics_DS0008__514__20_02_11\fused\volumes");
% tifFilePath        = custom_fullfile(tifDir, tifFileName);
% 
% datasetScratchDir = "_flippedMesh";
% 
% scratchDir = custom_fullfile("F:\Mariia\SerosaDynamics_DS0008__514__20_02_11\fused\_scratch", datasetScratchDir);
% projectDir = scratchDir;
% chartDir   = custom_fullfile(projectDir, "geodesicProjections", "fields", "data");
% meshDir = custom_fullfile("F:\Mariia\SerosaDynamics_DS0008__514__20_02_11\fused\smooth_meshes");
% meshFilePath = custom_fullfile(meshDir, "meshlab_smoothed_fused_tp_233.obj");
% 
% outputDir = custom_fullfile("F:\Mariia\SerosaDynamics_DS0008__514__20_02_11\fused\cartography_projections");
% 
% mkdir(projectDir)
% mkdir(chartDir)
% 
% parameters.timepoint       = 1;
% parameters.stackSize       = [481 985 505 1];
% parameters.createMesh      = "no"; % yes or no
% parameters.performClosure  = "no"; % yes or no (yes for after gastrulation)
% parameters.list_of_charts  = ["cylinder1", "cylinder2", "equidistant_1", "equidistant_2"];
% parameters.semseg.beadradius            =  100; % size in pixels
% parameters.semseg.beadthresh            =  0.3; % relative intensity between 0 and 1 (~ 0.7)
% parameters.semseg.BWthreshMan           =  0.5; % additive (can be zero)
% parameters.semseg.alpha                 =   70; % inf is convex hull, too small makes holes (~ 245)
% parameters.semseg.alphaRegionThreshold  = 1000; % max volume of regions to suppress (~ bead volume)
% parameters.layers.number                =   11; %41
% parameters.layers.distance              =    1; %in pixels
% parameters.layers.ref_layer_cylinder    =   13; %23
% parameters.layers.ref_layer_equidistant =   11; %21

    
%% Geodesic projection and creation of the charts
% NB! if there is some issue here make sure that imsane/setup.m has run to initialize
% NB! for actin data: copy xp-file and mesh-obj-file from nuclei segmentation folder and change name
% timepoint = 233;
for timepoint = 234:332
    tifFileName = "fused_tp_" + timepoint + "_ch_0.tif";
    tifDir         = custom_fullfile("F:\Mariia\SerosaDynamics_DS0008__514__20_02_11\fused\rotated_volumes");
    tifFilePath        = custom_fullfile(tifDir, tifFileName);
    
    datasetScratchDir = "_rot_tp_"+timepoint;
    
    scratchDir = custom_fullfile("F:\Mariia\SerosaDynamics_DS0008__514__20_02_11\fused\_scratch", datasetScratchDir);
    projectDir = scratchDir;
    chartDir   = custom_fullfile(projectDir, "geodesicProjections", "fields", "data");
    meshDir = custom_fullfile("F:\Mariia\SerosaDynamics_DS0008__514__20_02_11\fused\embryo_mesh");
    meshFilePath = custom_fullfile(meshDir, "rot_flipxy_meshlab_smoothed_fused_tp_233.obj");
    
    outputDir = custom_fullfile("F:\Mariia\SerosaDynamics_DS0008__514__20_02_11\fused\cartography_projections_2_layers");
    
    mkdir(projectDir)
    mkdir(chartDir)

    GeodesicProjections(tifDir, projectDir, ...
            tifFileName, meshFilePath, scratchDir, outputDir, ...
            parameters.stackSize, parameters.timepoint, ...
            parameters.layers.number, parameters.layers.distance);

    % cylinderMIPFilePath = custom_fullfile(scratchDir, "geodesicProjections", "fields", "data_MIP", "cylinder1_index", "cylinder1", "cmp_1_1_T0001.tif");
    % cylinderMIPOutputPath = custom_fullfile(outputDir, "cylinder_1_MIP_TP_"+timepoint+".tif");
    % movefile(cylinderMIPFilePath, cylinderMIPOutputPath);
    % cylinderMIPFilePath = custom_fullfile(scratchDir, "geodesicProjections", "fields", "data_MIP", "cylinder2_index", "cylinder2", "cmp_1_1_T0001.tif");
    % cylinderMIPOutputPath = custom_fullfile(outputDir, "cylinder_2_MIP_TP_"+timepoint+".tif");
    % movefile(cylinderMIPFilePath, cylinderMIPOutputPath);
    serosaLayersRange = 28:37;
    embryoLayersRange = 56:61;


    nr_of_layers = (parameters.layers.number - 1) /2;

    layer_paths_combined = {};
    % Load outer layers
    for i = nr_of_layers:-1:1
        layer_paths_combined{end+1} = custom_fullfile(scratchDir, "geodesicProjections", "fields","data_layer_m" + i); 
    end

    % Load middle layer
    layer_paths_combined{end+1} = custom_fullfile(scratchDir, "geodesicProjections", "fields","data");


    % Load inner layers
    for i = 1:nr_of_layers
        layer_paths_combined{end+1} = custom_fullfile(scratchDir, "geodesicProjections", "fields","data_layer_p" + i); 
    end

    layersStackDir = custom_fullfile(outputDir, "layersStack", "layersStack_cyl_"+timepoint+"_tp");
    mkdir(layersStackDir);
    for i = 1:parameters.layers.number
        layersStackPath = custom_fullfile(layersStackDir, "layer_" + i + ".tif");
        copyfile(custom_fullfile(layer_paths_combined{i},  "cylinder1_index","cylinder1", "cmp_1_1_T0001.tif"), layersStackPath);
    end

    % Load and do maximum projection of serosa layers
    serosa_mip_dir = custom_fullfile(outputDir, "extraembryonic_membranes");
    mkdir(serosa_mip_dir);

    cyl_1_init_path = custom_fullfile(layer_paths_combined{serosaLayersRange(1)}, "cylinder1_index","cylinder1", "cmp_1_1_T0001.tif");
    cyl_2_init_path = custom_fullfile(layer_paths_combined{serosaLayersRange(1)}, "cylinder2_index","cylinder2", "cmp_1_1_T0001.tif");
    MIP_cyl_1 = tiffreadVolume(cyl_1_init_path);    
    MIP_cyl_2 = tiffreadVolume(cyl_2_init_path);    
    for i = serosaLayersRange
        cyl_1_path = custom_fullfile(layer_paths_combined{i}, "cylinder1_index","cylinder1", "cmp_1_1_T0001.tif");
        cyl_2_path = custom_fullfile(layer_paths_combined{i}, "cylinder2_index","cylinder2", "cmp_1_1_T0001.tif");
        layer_cyl_1 = tiffreadVolume(cyl_1_path);
        layer_cyl_2 = tiffreadVolume(cyl_2_path);
        MIP_cyl_1 = max(MIP_cyl_1, layer_cyl_1);
        MIP_cyl_2 = max(MIP_cyl_2, layer_cyl_2);
    end
    serosa_cyl_1_mip_path = custom_fullfile(serosa_mip_dir, "extr_memb_cyl_1_MIP_tp_" + timepoint + ".tif");
    serosa_cyl_2_mip_path = custom_fullfile(serosa_mip_dir, "extr_memb_cyl_2_MIP_tp_" + timepoint + ".tif");
    imwrite(MIP_cyl_1, serosa_cyl_1_mip_path, 'tif');
    imwrite(MIP_cyl_2, serosa_cyl_2_mip_path, 'tif');

    % Load and do maximum projection of embryo layers
    embryo_mip_dir = custom_fullfile(outputDir, "embryo");
    mkdir(embryo_mip_dir);

    cyl_1_init_path = custom_fullfile(layer_paths_combined{embryoLayersRange(1)}, "cylinder1_index","cylinder1", "cmp_1_1_T0001.tif");
    cyl_2_init_path = custom_fullfile(layer_paths_combined{embryoLayersRange(1)}, "cylinder2_index","cylinder2", "cmp_1_1_T0001.tif");
    MIP_cyl_1 = tiffreadVolume(cyl_1_init_path);    
    MIP_cyl_2 = tiffreadVolume(cyl_2_init_path);    
    for i = embryoLayersRange
        cyl_1_path = custom_fullfile(layer_paths_combined{i}, "cylinder1_index","cylinder1", "cmp_1_1_T0001.tif");
        cyl_2_path = custom_fullfile(layer_paths_combined{i}, "cylinder2_index","cylinder2", "cmp_1_1_T0001.tif");
        layer_cyl_1 = tiffreadVolume(cyl_1_path);
        layer_cyl_2 = tiffreadVolume(cyl_2_path);
        MIP_cyl_1 = max(MIP_cyl_1, layer_cyl_1);
        MIP_cyl_2 = max(MIP_cyl_2, layer_cyl_2);
    end
    embryo_cyl_1_mip_path = custom_fullfile(embryo_mip_dir, "embryo_cyl_1_MIP_tp_" + timepoint + ".tif");
    embryo_cyl_2_mip_path = custom_fullfile(embryo_mip_dir, "embryo_cyl_2_MIP_tp_" + timepoint + ".tif");
    imwrite(MIP_cyl_1, embryo_cyl_1_mip_path, 'tif');
    imwrite(MIP_cyl_2, embryo_cyl_2_mip_path, 'tif');
end

%% Combine cartography layers to stack to check
    serosaLayersRange = 28:37;
    embryoLayersRange = 56:61;


    nr_of_layers = (parameters.layers.number - 1) /2;

    layer_paths_combined = {};
    % Load outer layers
    for i = nr_of_layers:-1:1
        layer_paths_combined{end+1} = custom_fullfile(scratchDir, "geodesicProjections", "fields","data_layer_m" + i); 
    end

    % Load middle layer
    layer_paths_combined{end+1} = custom_fullfile(scratchDir, "geodesicProjections", "fields","data");


    % Load inner layers
    for i = 1:nr_of_layers
        layer_paths_combined{end+1} = custom_fullfile(scratchDir, "geodesicProjections", "fields","data_layer_p" + i); 
    end

    layersStackDir = custom_fullfile(outputDir, "layersStack", "layersStack_cyl_1_tp");
    mkdir(layersStackDir);
    for i = 1:parameters.layers.number
        layersStackPath = custom_fullfile(layersStackDir, "layer_" + i + ".tif");
        copyfile(custom_fullfile(layer_paths_combined{i},  "cylinder1_index","cylinder1", "cmp_1_1_T0001.tif"), layersStackPath);
    end

%% MIPs for serosa and embryo
    % Load and do maximum projection of serosa layers
    serosa_mip_dir = custom_fullfile(outputDir, "extraembryonic_membranes");
    mkdir(serosa_mip_dir);

    cyl_1_init_path = custom_fullfile(layer_paths_combined{serosaLayersRange(1)}, "cylinder1_index","cylinder1", "cmp_1_1_T0001.tif");
    cyl_2_init_path = custom_fullfile(layer_paths_combined{serosaLayersRange(1)}, "cylinder2_index","cylinder2", "cmp_1_1_T0001.tif");
    MIP_cyl_1 = tiffreadVolume(cyl_1_init_path);    
    MIP_cyl_2 = tiffreadVolume(cyl_2_init_path);    
    for i = serosaLayersRange
        cyl_1_path = custom_fullfile(layer_paths_combined{i}, "cylinder1_index","cylinder1", "cmp_1_1_T0001.tif");
        cyl_2_path = custom_fullfile(layer_paths_combined{i}, "cylinder2_index","cylinder2", "cmp_1_1_T0001.tif");
        layer_cyl_1 = tiffreadVolume(cyl_1_path);
        layer_cyl_2 = tiffreadVolume(cyl_2_path);
        MIP_cyl_1 = max(MIP_cyl_1, layer_cyl_1);
        MIP_cyl_2 = max(MIP_cyl_2, layer_cyl_2);
    end
    serosa_cyl_1_mip_path = custom_fullfile(serosa_mip_dir, "extr_memb_cyl_1_MIP_tp_" + timepoint + ".tif");
    serosa_cyl_2_mip_path = custom_fullfile(serosa_mip_dir, "extr_memb_cyl_2_MIP_tp_" + timepoint + ".tif");
    imwrite(MIP_cyl_1, serosa_cyl_1_mip_path, 'tif');
    imwrite(MIP_cyl_2, serosa_cyl_2_mip_path, 'tif');

    % Load and do maximum projection of embryo layers
    embryo_mip_dir = custom_fullfile(outputDir, "embryo");
    mkdir(embryo_mip_dir);

    cyl_1_init_path = custom_fullfile(layer_paths_combined{embryoLayersRange(1)}, "cylinder1_index","cylinder1", "cmp_1_1_T0001.tif");
    cyl_2_init_path = custom_fullfile(layer_paths_combined{embryoLayersRange(1)}, "cylinder2_index","cylinder2", "cmp_1_1_T0001.tif");
    MIP_cyl_1 = tiffreadVolume(cyl_1_init_path);    
    MIP_cyl_2 = tiffreadVolume(cyl_2_init_path);    
    for i = embryoLayersRange
        cyl_1_path = custom_fullfile(layer_paths_combined{i}, "cylinder1_index","cylinder1", "cmp_1_1_T0001.tif");
        cyl_2_path = custom_fullfile(layer_paths_combined{i}, "cylinder2_index","cylinder2", "cmp_1_1_T0001.tif");
        layer_cyl_1 = tiffreadVolume(cyl_1_path);
        layer_cyl_2 = tiffreadVolume(cyl_2_path);
        MIP_cyl_1 = max(MIP_cyl_1, layer_cyl_1);
        MIP_cyl_2 = max(MIP_cyl_2, layer_cyl_2);
    end
    embryo_cyl_1_mip_path = custom_fullfile(embryo_mip_dir, "embryo_cyl_1_MIP_tp_" + timepoint + ".tif");
    embryo_cyl_2_mip_path = custom_fullfile(embryo_mip_dir, "embryo_cyl_2_MIP_tp_" + timepoint + ".tif");
    imwrite(MIP_cyl_1, embryo_cyl_1_mip_path, 'tif');
    imwrite(MIP_cyl_2, embryo_cyl_2_mip_path, 'tif');