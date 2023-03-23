Matlab scripts enclosed in this repository are those used to generate the data and the figures for the manuscript ‘False friends: How mutualism impedes community assembly and creates empty niches for invasion’, by HO Minoarivelo, U. Dieckmann, and C. Hui
-------------------------------------------------------------

The codes were run using Matlab 2020b on a device with the following specifications:
OS:  Windows 10 64-bit, Processor: Intel Core i5-10210U CPU @ 1.60GHz 2.10GHz, RAM: 8Gb
-------------------------------------------------------------

To generate the simulated communities, the scripts should be run in the following order:

- Run 'main_model_intra_guild_resources.m'.
This will generate the communities until an ESS is reached for the case where the two ancestral traits evolve to optimize their intake of intra-guild resources.
This code runs over all ranges of parameters. A folder named 'Results_intra_guild' is created and data on population abundances and trait values are stored in there.
This code took approximately 10 hours to run.

- Run 'main_model_cross_guild_resources.m'. 
This will generate the communities until an ESS is reached for the case where the two ancestral traits evolve to optimize their intake of cross-guild resources.
This code runs over all ranges of parameters. A folder named 'Results_cross_guild' is created and data on population abundances and trait values are stored in there.
This code took approximately 13 hours to run.

- Run 'invasion_EN_after_ESS.m'. 
This will take the simulations after an ESS and continue with the introduction of new species and traits evolution until ESS for those communities presenting empy niches.
This code runs over all ranges of parameters for both the two cases of resource intake optimization (intra or cross guilds). Folders named 'Results_invasion_intra_guild' and 'Results_invasion_cross guild'
are created and data on trait values and population densities are stored in there.
This code took approximately 9 hours to run.

-Run 'main_resource_competition_model.m' located within the folder 'Resource_competition_model' to simulate under a model where there is only competition and no mutualism.
This will create a folder named 'Results_RC' where trait values and population densities are stored. It runs only over one set of parameter, the one needed to plot figure2.
This code takes approximately 40 seconds to run. 
 
- To plot the figures, run 'Plot_Fig(figure number).m'. For example, run 'Plot_Fig1.m' to plot figure 1, or 'Plot_Fig2.m' to plot figure 2. All plots are stored in the folder 'Plots' which will be generated.

- Plotting figure 1 requires that the needed data are present in the folders 'Results_intra_guild' and 'Results_cross_guild'.
'Plot_Fig1.m' takes approximately 1mn to run.
- Plotting figure 2 requires that the needed data are present in the folders 'Results_intra_guild', 'Results_cross_guild' and 'Results_RC' (located in the forlder 'Resource_competition_model'). 
'Plot_Fig2.m' takes approximately 6mn to run.
- Plotting figure 3 requires that all data are present in the folders 'Results_intra_guild', 'Results_cross_guild', 'Results_invasion_intra_guild' and 'Results_invasion_cross_guild'.
'Plot_Fig3.m' takes approximately 13mn to run.
- Plotting figure 4 requires that the needed data are present in the folders 'Results_intra_guild', 'Results_cross_guild', 'Results_invasion_intra_guild' and 'Results_invasion_cross_guild'.
'Plot_Fig4.m' takes approximately 5mn to run.
- Plotting figure 5 requires that all data are present in the folders 'Results_intra_guild', 'Results_cross_guild', 'Results_invasion_intra_guild' and 'Results_invasion_cross_guild'.
'Plot_Fig5.m' takes approximately 36mn to run.
- Plotting figure 6 requires that all data are present in the folders 'Results_intra_guild', 'Results_cross_guild', 'Results_invasion_intra_guild' and 'Results_invasion_cross_guild'.
'Plot_Fig6.m' creates two folders: 'Selection_pressure_intra_guild' and 'Selection_pressure_cross_guild' which store the values of the selection gradients over time, before outputing figure 6.
'Plot_Fig6.m' takes approximately 1h20mn to run.

- To generate the figures in the supplementary material B (Figures B1 to B11), i.e. before running 'Plot_FigB...', it is required that data are available in the folders 'Results_intra_guild', 'Results_cross_guild', 'Results_invasion_intra_guild' and 'Results_invasion_cross_guild'.  
_ To generate the figures in the supplementary materia C (Figures C1 to C5), i.e. before running 'Plot_FigC...', the code 'robustness_test.m' shoudl be run first. This code will generate the additional data needed for the robustness analysis. 


