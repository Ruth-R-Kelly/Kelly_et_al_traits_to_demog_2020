#### Code and data files for "Climatic and evolutionary contexts are required to infer plant life history strategies from functional traits at a global scale"

##### Supplementary material for Kelly et al. Ecology Letters 2021
<br/>

This folder contains the code and data files required for the analysis presented in the paper - "Climatic and evolutionary contexts are required to infer plant life history strategies from functional traits at a global scale" 
For more details on the hypothesis, methods and results please refer to the main paper. 

**Code and analysis**
There are two main stages in the analysis, 

   -  The calculation of the Life history metrics from the Compadre database,
       and it's linkage to the functional traits databases.
   -  Attaching the environmental and phylogenetic data to this data. 
    
These are stored in the folders - "1_Calculating_Life_History_Metrics" and   "2_Joining_Env_Phylo_and_GLMMing" respectively.  

Inside each of these folders there is a file called overview code, this is an rmarkdown file that describes and runs all of the analyses in that folder. 

I have also added a third folder 3 - Working_with_glmm_outputs, this folder contains a few bits of code that were used after the fitting of the main multiple response model presented in the paper. Relating to model variance calculation and some plotting, which may be of interest to some readers.

**Data**

If you'd prefer to just have the data for the derived demography metrics with matching climate and trait data so you can do your own thing, you will find it in the file 'data_with_needed_env_vars_for_GLMM_analysis' in the home folder. Similarly, 'reduced_tree_for_traits_03_2020.nex' contains the matching phylogeny. 

Please keep in mind that the databases we took our original data from are constantly growing and being updated. Therefore, if you want to do further more detailed research in this area we refer you to the following datasources which were invaluable to us.
<br/>

##### Life history
Compadre - https://compadre-db.org/

Salguero-Gomez, R., Jones, O.R., Archer, C.R., Buckley, Y.M., Salguero-g, R., Che-castaldo, J., et al. (2015). The COMPADRE Plant Matrix Database : an open online repository for plant demography. J. Ecol., 202–218.

##### Traits
TRY - https://www.try-db.org/
Kattge, J., Bönisch, G., Díaz, S., Lavorel, S., Prentice, I.C., Leadley, P., et al. (2020). TRY plant trait database – enhanced coverage and open access. Glob. Chang. Biol., 26, 119–188.

Bien - https://bien.nceas.ucsb.edu/bien/
Enquist, B., Condit, R., Peet, R., Schildhauer, M. & Thiers, B. (2016). Cyberinfrastructure for an integrated botanical information network to investigate the ecological impacts of global climate change on plant biodiversity. PeerJ Prepr., e2615v2.

##### Climate
http://www.csi.cgiar.org
Trabucco, A. & Zomer, R.J. (2009). Global Aridity Index (Global-Aridity) and Global Potential Evapo-Transpiration (Global-PET) Geospatial Database. CGIAR Consort. Spat. Information.

http://www.worldclim.org
Fick, S.E. & Hijmans, R.J. (2017). WorldClim 2: new 1‐km spatial resolution climate surfaces for global land areas. Int. J. Climatol., 37, 4302–4315.

##### General
Dryad - https://datadryad.org

##### The phylogeny used here is an edited subset of: 
Zanne, Amy E. et al. (2014), Data from: Three keys to the radiation of angiosperms into freezing environments, Dryad, Dataset, https://doi.org/10.5061/dryad.63q27


