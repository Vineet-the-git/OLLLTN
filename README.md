# OLLLTN
In several real world problems, we need to model data with a variety of shapes such as data having two peaks, or bimodal data. We often have to resort to mixture distributions in order to model them. The paper [1] proposes a single distribution which can account for the bimodality of our data. In this report, we have established some theoretical properties and practical applications of the proposed OLLLTN distribution. We began by defining its PDF, CDF, moments, skewness, and kurtosis, and then illustrated some of them using plots. We went ahead and implemented functions to sample from the distribution and return values of the above properties for the distribution. We then described OLLLTN regression and tested the same on the dataset we generated using our sampler. Finally we modelled the Human Development Index (HDI) data of two Brazilian cities using OLLLTN distribution and compared the quality of fit with other distributions to establish its practicality.
We were able to reproduce most of the results and findings that are proposed in the paper and implement the same in R. The pdf of our presentation slides and the Brazilian HDI dataset is uploaded here.

## To Run the code 
`Run pitchers.r`  

## You can see our results in-
1. [Project Report](https://github.com/Vineet-the-git/OLLLTN/blob/main/pitchers.pdf)  
2. [Presentation](https://github.com/Vineet-the-git/OLLLTN/blob/main/Presentation.pdf)
