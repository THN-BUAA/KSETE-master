# KSETE-master
Source code and all the datasets used in the paper (Kernel Spectral Embedding Transfer Ensemble for Heterogeneous Defect Prediction).

## There are the following three main subfolders in KSETE-master
- Datasets -- This folder includes three subfolders (AEEEM, JIRA, and PROMISE), each subfolder contains the corresponding datasets used in the paper. Each dataste is a ARFF file.
- Datasets-inOne -- All the datasets are put in one folder.
- KSETE -- This folder includes the source code of KSETE and an example to use it(i.e., Experiment.m).

## Usage of KSETE
- Step1: Download KSETE-master;
- Step2: Open MATLAB, add the path of KSETE-master and its subfolders;
- Step3: Add the absolute path of weka.jar and SMOTE.jar in the folder KSETE into the classpath.txt of MATLAB.
- Step4: Enter folder KSETE and open the script Experiment.m in MATLAB;
- Step5: Click the run button in the manu of MATLAB Editor to run the example; 

Note that owing to the randomness in the KSETE, the results of RQ1 and RQ2 may have a slight fluctuation compared with that in the paper.



## Contact
If you have any problem about our code, feel free to contact
- tonghaonan@buaa.edu.cn
or describe your problem in Issues.
