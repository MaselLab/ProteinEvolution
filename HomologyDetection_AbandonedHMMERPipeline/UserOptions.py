def UserOptions():


    # Max Index = 23
    ###########################################################
    ############### Input and Temporary Files  ################
    ###########################################################

    # Databases
    FullGenesDatabase = './Databases/ProteinCodingMouseGenes.fasta' #[0]
    FullGenesDictionary = './Databases/Genes_Dictionary.txt'#[1]

    FragmentedGenesDatabase = './Databases/PartitionedGenesDatabase.fasta' #[2]
    FragmentedGenesDictionary = './Databases/PartitionedGenesDatabase.dictionary' #[3]

    # Temporary Files
    inputGeneFile = './TemporaryFiles/InputSequence.fasta' #[4]
    pHMMERFile = './TemporaryFiles/pHMMER.out' #[5]
    HMMFile = './TemporaryFiles/HMMFile.hmm' #[6]
    HitsFastaFile = './TemporaryFiles/Hits.fasta' #[7]
    MSAFile = './TemporaryFiles/Hits.msa' #[8]
    HMMSearchFile = './TemporaryFiles/HMMSearchOutput.txt' #[9]
    TemporaryFastaFile = './TemporaryFiles/Temp.fasta' #[21]
    TemporaryMSA = './TemporaryFiles/Temp.msa' #[22]

    # Input UIDs (arranged from shortest to longest proteins)
    InputUIDs = './Databases/UIDInputs_AscendingLengths.txt' #[23]

    ###########################################################
    ###########################################################
    ###########################################################



    ###########################################################
    ################## HMMER Parameters  ######################
    ###########################################################
    maximumIterations = 10 #[10]

    EThreshold = float(1e-6) #[11]

    domEThreshold = float(1e-6) #[12]

    incEThreshold = float(1e-20) #[13]

    incdomEThreshold = float(1e-20) #[14]

    bitBiasRatio = .5 #[15]

    ###########################################################
    ###########################################################
    ###########################################################

    
    ###########################################################
    ################## Other Parameters  ######################
    ###########################################################

    changePointPenalty = 10 #[16]

    maskIfInitialBiasExceeded = False #[17]

    minDistanceBetweenChangePoints = 20 #[18]

    ViewChangePoints = False #[20]

    ###########################################################
    ###########################################################
    ###########################################################



    ###########################################################
    ####################### Output Files ######################
    ###########################################################

    outputGenesDictionary = '../Output/Results.dictionary' #[19]


    ###########################################################
    ###########################################################
    ###########################################################

    
    return FullGenesDatabase, FullGenesDictionary, FragmentedGenesDatabase, FragmentedGenesDictionary, inputGeneFile, pHMMERFile, HMMFile, HitsFastaFile, MSAFile, HMMSearchFile, maximumIterations, EThreshold, domEThreshold, incEThreshold, incdomEThreshold, bitBiasRatio, changePointPenalty, maskIfInitialBiasExceeded, minDistanceBetweenChangePoints, outputGenesDictionary, ViewChangePoints, TemporaryFastaFile, TemporaryMSA, InputUIDs
