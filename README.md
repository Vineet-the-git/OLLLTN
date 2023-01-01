# OLLLTN
In several real world problems, we need to model data with a variety of shapes such as data having two peaks, or bimodal data. We often have to resort to mixture distributions in order to model them. The paper [1] proposes a single distribution which can account for the bimodality of our data. In this report, we have established some theoretical properties and practical applications of the proposed OLLLTN distribution. We began by defining its PDF, CDF, moments, skewness, and kurtosis, and then illustrated some of them using plots. We went ahead and implemented functions to sample from the distribution and return values of the above properties for the distribution. We then described OLLLTN regression and tested the same on the dataset we generated using our sampler. Finally we modelled the Human Development Index (HDI) data of two Brazilian cities using OLLLTN distribution and compared the quality of fit with other distributions to establish its practicality.
We were able to reproduce most of the results and findings that are proposed in the paper and implement the same in R. The pdf of our presentation slides and the Brazilian HDI dataset is uploaded here.  

-> The code is in the pitchers.R file.   
-> It requires two libraries, "gamlss" and "pracma"  
Install them if requiered by the install.packages("gamlss") and install.packages("pracma")  

-> The code is divided section wise according to the report  
-> In some places, a path to the output location for plots needs to be specified. The plots will be saved at those locations   
-> To implement code in section 7, the brazilian cities dataset is required (the exact location is commented in the code) (line 669)  
-> The dataset is added to the zip folder and is also uploaded to this drive link: https://drive.google.com/drive/folders/1jZX9k-o4lNetEXvQr-wnQeAEX3RheZsW  
-> Our presentation is in the "Presentation.pdf" file  
-> Our report is in the "pitchers.pdf" file  
-> The latex folder contains the raw latex files for generating the report  
-> The plots have been added to a pdf file separately as well in "plots.pdf"  

## To Run the code 
`Run pitchers.r`  

## You can see our results in-
1. [Project Report](https://github.com/Vineet-the-git/OLLLTN/blob/main/pitchers.pdf)  
2. [Presentation](https://github.com/Vineet-the-git/OLLLTN/blob/main/Presentation.pdf)
