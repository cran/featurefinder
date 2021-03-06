% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/featurefinder.r
\name{findFeatures}
\alias{findFeatures}
\title{findFeatures}
\usage{
findFeatures(outputPath = "NoPath", fcsv, exclusionVars, factorToNumericList,
  treeGenerationMinBucket = 20, treeSummaryMinBucket = 50,
  treeSummaryResidualThreshold = 0,
  treeSummaryResidualMagnitudeThreshold = 0, doAllFactors = TRUE,
  maxFactorLevels = 20, useSubDir = TRUE, tempDirFolderName = "")
}
\arguments{
\item{outputPath}{A string containing the location of the input csv file. Results are also stored in this location. Set to "NoPath" to use tempdir() or leave blank}

\item{fcsv}{A string containing the name of a csv file}

\item{exclusionVars}{A string consisting of a list of variable names with double quotes around each variable}

\item{factorToNumericList}{A list of variable names as strings}

\item{treeGenerationMinBucket}{Desired minimum number of data points per leaf (default 20)}

\item{treeSummaryMinBucket}{Minimum number of data points in each leaf for the summary (default 50)}

\item{treeSummaryResidualThreshold}{Minimum residual in the summary (default 0 for positive residuals)}

\item{treeSummaryResidualMagnitudeThreshold}{Minimum residual magnitude in the summary (default 0 i.e. no restriction)}

\item{doAllFactors}{Flag to indicate whether to analyse the levels of all factor variables (default TRUE)}

\item{maxFactorLevels}{Maximum number of levels per factor before it is converted to numeric (default 20)}

\item{useSubDir}{Flag to specify whether the partition trees should be saved in the current directory or a subdirectory}

\item{tempDirFolderName}{specify a subfolder name if writing multiple scans to the temporary directory}
}
\value{
outputPath returns the location of the output for reference in addFeatures and for any other purpose. Saves residual CART trees and associated highlighted residuals for each to the path provided.
}
\description{
Perform analysis of residuals grouped by factor to identify features which explain the target variable
}
\examples{

require(featurefinder)
data(futuresdata)
data=futuresdata
data$SMIfactor=paste("smi",as.matrix(data$SMIfactor),sep="")
n=length(data$DAX)
nn=floor(length(data$DAX)/2)

# Can we predict the relative movement of DAX and SMI?
data$y=data$DAX*0 # initialise the target to 0
data$y[1:(n-1)]=((data$DAX[2:n])-(data$DAX[1:(n-1)]))/
  (data$DAX[1:(n-1)])-(data$SMI[2:n]-(data$SMI[1:(n-1)]))/(data$SMI[1:(n-1)])

# Fit a simple model
thismodel=lm(formula=y ~ .,data=data)
expected=predict(thismodel,data)
actual=data$y
residual=actual-expected
data=cbind(data,expected, actual, residual)

CSVPath=tempdir()
fcsv=paste(CSVPath,"/futuresdata.csv",sep="")
write.csv(data[(nn+1):(length(data$y)),],file=fcsv,row.names=FALSE)
exclusionVars="\\"residual\\",\\"expected\\", \\"actual\\",\\"y\\""
factorToNumericList=c()

# Now the dataset is prepared, try to find new features
findFeatures(outputPath="NoPath", fcsv, exclusionVars,factorToNumericList,                     
         treeGenerationMinBucket=50,
         treeSummaryMinBucket=20,
         useSubDir=FALSE)  
         
newfeat1=((data$SMIfactor==0) & (data$CAC < 2253) & (data$CAC< 1998) & (data$CAC>=1882)) * 1.0
newfeat2=((data$SMIfactor==1) & (data$SMI < 7837) & (data$SMI >= 7499)) * 1.0
newfeatures=cbind(newfeat1, newfeat2) # create columns for the newly found features
datanew=cbind(data,newfeatures)
thismodel=lm(formula=y ~ .,data=datanew)
expectednew=predict(thismodel,datanew)

requireNamespace("Metrics")
OriginalRMSE = Metrics::rmse(data$y,expected)
NewRMSE = Metrics::rmse(data$y,expectednew)

print(paste("OriginalRMSE = ",OriginalRMSE))
print(paste("NewRMSE = ",NewRMSE))
}
\keyword{findFeatures}
