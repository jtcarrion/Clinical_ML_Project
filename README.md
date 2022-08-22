# Clinical_ML_Project
Developed and analyzed various ML models, Decision Tree, to classify breast cancer tumor cells as malignant or benign. Code provided in this repository is responsible for testing an array of Decision Tree models to find the most accurate and useful model for clinical use. After 10-fold cross validation, the tree with a minumun split of 25 and a max depth of 7 resulted in a model that can classify maligant and benign tumor cells with 92% accuracy 

DATA:
The data provided consists of 700 samples each with 13 measurments. After cleaning up the data, PCA was used to reduce the complexity of dataset and allowed us to see a clear distinction between two groups using the top 2-3 PCs

Model:
Preliminary data using Weka suggested a Decsion tree will be the most effective appraoch for this dataset so we developed and tested various deciosn tree models in order to test the effects of variables like Min_Split and Max_Depth

Results:
After 10-fold cross validation, the tree with a minumun split of 25 and a max depth of 7 resulted in a model that can classify maligant and benign tumor cells with 92% accuracy
