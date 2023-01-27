
This project implements the experiments section of the paper K. Musayeva and M. Binois, "Revisiting multi-label propagation: the case of small data", 2023.

The repository is as follows:
    
    main.R
    |
    src: contains the R code of the project
    |
    data: contains the data used in the paper
    |
    results: the results of the experiments are saved into this folder


Running instructions:

    main.R : this is where the project should be run. This file contains the examples for running each method considered in the paper.

    The following variables should be specified: the name of the method, the decision function used and the hyperparameter list. 
            
    The evaluation measures used and the parameters of the cross-validation should be modified in the global_var.R file.

    To run the harmonic function as a stacking method, you need to first run the method, the predictions of which you want to 

    use in the stacking (the results are saved automatically into results/name_data_*.obj file). Then the compromise solution, 

    more precisely, the hyperparameter configuration corresponding to the compromise solution should be computed from this .obj file. 

    This hyperparameter configuration is then used in the stacking. 

    The reported results are obtained by applying the ks_solution to the results/name_data_*.obj file.


Some reproducibiltiy notes:

    In the consistency method ("Learning with local and global consistency," D. Zhou et al., 2004) and  in the locally linear neighbourhood method

    (“Label propagation through linear neighborhoods,” Wang and C. Zhang, 2007) the regularization parameter alpha is defaulted to 0.99. In our work,

    we experimented with different values for this parameter and for the scene dataset, we found that very small values of alpha, for instance, {0.05, 0.1} 

    (this implies smaller input from the propagation matrix) coupled with the class-mass normalization thresholding strategy produced superior results,

    however, if we couple it with the label cardinality optimizer the results are meaningless. 
    

        






