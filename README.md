## Matlab scripts for simulation and visualisation in the manuscript "False friends: How mutualism impedes community assembly and creates empty niches for invasion" 
### by H.O. Minoarivelo, U. Dieckmann, and C. Hui
(submitted for review in The American Naturalist)
-------------------------------------------------------------
The codes were run using Matlab 2020b on a device with the following specifications:
OS:  Windows 10 64-bit, Processor: Intel Core i5-10210U CPU @ 1.60GHz 2.10GHz, RAM: 8Gb.

To simulate the assembly of ecological communities wth trait-mediated competitive and mutualistic interactions facing biological invasions, one can run the scripts in the following order:

- Run 'main_model_intra_guild_resources.m'
This will generate pre-invasion communities that the initial pair of resident species has reached the evolutionarily stable strategy (ESS) when the adaptive traits of the two ancestral species evolved to optimize their intakes of intra-guild resources. This code runs over all ranges of parameters. A folder named 'Results_intra_guild' is created, and data on population densities and trait values are stored therein. This code takes approximately 10 hours to run in the specified computational system.

- Run 'main_model_cross_guild_resources.m' 
This will generate pre-invasion communities that the initial pair of resident specis has reached the ESS when the adaptive traits of the two ancestral species evolved to optimize their intakes of cross-guild resources. This code runs over all ranges of parameters. A folder named 'Results_cross_guild' is created, and data on population densities and trait values are stored therein. This code takes approximately 13 hours to run in the specified computational system.

- Run 'invasion_EN_after_ESS.m'
This will simulate biological invasions into communities with empty niches. In one round it introduces one species and simulates trait evolution in the invaded community till a new ESS has been reached. Multiple rounds are run till the final ESS where no empty niche in the trait space is present. This code runs over all ranges of parameters for both intra- and cross-guild cases of resource intake optimization. Folders named 'Results_invasion_intra_guild' and 'Results_invasion_cross guild' are created, and data on trait values and population densities are stored therein. This code takes approximately 9 hours to run in the specifid computational system.

- Run 'main_resource_competition_model.m' 
This is located within the folder 'Resource_competition_model' and is used to simulate a simplifed model where the community contains only competitive interactions and lacks mutualistic interactions. This will create a folder named 'Results_RC' where trait values and population densities are stored. It runs only over one set of parameters and is needed for plotting figure2. This code takes approximately 40 seconds to run. 
 
- Run 'Plot_Fig(figure number).m' to produce figures
For example, run 'Plot_Fig1.m' to plot figure 1, or 'Plot_Fig2.m' to plot figure 2. All plots are stored in a folder created and named 'Plots'. Note the following specifics for producing specific figures:

  * Plotting figure 1 requires that the needed data are present in the folders 'Results_intra_guild' and 'Results_cross_guild'. 'Plot_Fig1.m' takes approximately 1 minute to run.
  * Plotting figure 2 requires that the needed data are present in the folders 'Results_intra_guild', 'Results_cross_guild' and 'Results_RC', located within the forlder 'Resource_competition_model'. 'Plot_Fig2.m' takes approximately 6 minutes to run.
  * Plotting figure 3 requires that all data are present in the folders 'Results_intra_guild', 'Results_cross_guild', 'Results_invasion_intra_guild' and 'Results_invasion_cross_guild'. 'Plot_Fig3.m' takes approximately 13 minutes to run.
  * Plotting figure 4 requires that the needed data are present in the folders 'Results_intra_guild', 'Results_cross_guild', 'Results_invasion_intra_guild' and 'Results_invasion_cross_guild'. 'Plot_Fig4.m' takes approximately 5 minutes to run.
  * Plotting figure 5 requires that all data are present in the folders 'Results_intra_guild', 'Results_cross_guild', 'Results_invasion_intra_guild' and 'Results_invasion_cross_guild'. 'Plot_Fig5.m' takes approximately 36 minutes to run.
  * Plotting figure 6 requires that all data are present in the folders 'Results_intra_guild', 'Results_cross_guild', 'Results_invasion_intra_guild' and 'Results_invasion_cross_guild'. 'Plot_Fig6.m' creates two folders: 'Selection_pressure_intra_guild' and 'Selection_pressure_cross_guild' which store the values of the selection gradients over time, before producing figure 6. 'Plot_Fig6.m' takes approximately 1h20mn to run.

- Run 'Plot_FigB(figure number).m' to produce figures in the supplementary material B (Figures B1 to B10). It requires that data are available in the folders 'Results_intra_guild', 'Results_cross_guild', 'Results_invasion_intra_guild' and 'Results_invasion_cross_guild'.

- Run 'Plot_FigC(figure number).m' to produce figures in the supplementary material C (Figures C1 to C5). Note, the code 'robustness_test.m' should be run first before producing these figures. This code will generate additional data needed for the robustness analysis.
