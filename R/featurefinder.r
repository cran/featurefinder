
#' @importFrom grDevices dev.off
#' @importFrom stats as.formula na.omit
#' @importFrom utils read.csv

#' @title saveTree
#' @description Generate a residual tree on a subset of the data specified by the factor level mainfaclev (main factor level)
#' @param data A dataframe containing the residual and some predictors
#' @param vars A list of candidate predictors
#' @param expr A expression to be modelled by the RPART tree
#' @param i An integer corresponding to the factor level
#' @param varname A string corresponding to the name of the factor variable being analysed
#' @param mainfaclev A level of the mainfac factor
#' @param treeGenerationMinBucket Minimum size for tree generation
#' @param ... and parameters to be passed through
#' @return A tree object
#' @keywords saveTree
#' @export 
#' @name saveTree

saveTree<-function(data,vars,expr,i,varname,mainfaclev,treeGenerationMinBucket,...) {
    
  dat=data[data$mainfac==mainfaclev,]
  print(treeGenerationMinBucket)
  
  # split out train and test ;
  ttd_flag=dat$ttd_flag
  nares=dat$residual[is.na(dat$residual)]
  if (length(nares)>0) {dat$residual[is.na(dat$residual)]=0}
  residual=as.numeric(as.vector(dat$residual))
  residual[which(is.na(residual))]=0
  dat=cbind(dat[,vars],ttd_flag)
  dat[is.na(dat)]=0
  dat[dat=="."]=0
  
  # convert to numeric and factor
  nn=names(dat)
  #class(dat[,nn[1]])
  getclasses<-function(nname) {
    return(class(dat[1,nname]))
  }
  dclasses=lapply(nn,FUN=getclasses)
  factorlist=which(dclasses=="factor")
  numerics=which(dclasses=="numeric")
  intlist=which(dclasses=="integer")
  
  dat[,intlist]=apply(FUN=as.numeric,as.data.frame(dat[,intlist]),MARGIN=2)
  dat[is.na(dat[,intlist]),intlist]=0
  dat[is.na(dat[,numerics]),numerics]=0
  dat=as.data.frame(dat,stringsAsFactors=FALSE)  
  
  any(is.na(dat[,intlist]))  
  any(is.na(dat[,numerics]))    
  any(is.na(dat[,factorlist]))  
  
  dat=cbind(residual,dat)
  dat_train=dat[dat$ttd_flag==0,]
  dat_test=dat[dat$ttd_flag==1,]

  fit <- rpart::rpart(expr, data = dat_test, control=rpart::rpart.control(minbucket=treeGenerationMinBucket))

  improve<-NULL # this fixes the NOTE in R CMD check, without throwing other errors
  try({
    datfit=fit$splits
    name=rownames(fit$splits)
    datfit=as.data.frame(cbind(name,datfit))
    fitarrange=datfit
    if (dim(datfit)[[1]]>0) {
      fitarrange=plyr::arrange(datfit,plyr::desc(improve)) # ~improve fixes the note but throws errors in the example
    } 
  })
  grDevices::png(paste(varname,'_',mainfaclev,'.png',sep=""), width=750, height = 500)
  rpart.plot::rpart.plot(fit, varlen=0, faclen=0, extra=1, digits=5)
  dev.off()
  grDevices::png(paste(varname,'_',mainfaclev,'_medres.png',sep=""), width=1800, height = 1200)
  rpart.plot::rpart.plot(fit, varlen=0, faclen=0, extra=1, digits=5)
  dev.off()
  grDevices::png(paste(varname,'_',mainfaclev,'_highres.png',sep=""), width=3000, height = 2000)
  rpart.plot::rpart.plot(fit, varlen=0, faclen=0, extra=1, digits=5)
  dev.off()
  return(list(fitarrange,fit))
}

#' @title generateTrees
#' @description  Generate a residual tree for each level of factor mainfac
#' @param data A dataframe
#' @param vars A list of candidate predictors
#' @param expr A expression to be modelled by the RPART tree
#' @param runname A string corresponding to the name of the variable being modelled
#' @param ... and parameters to be passed through
#' @keywords generateTrees
#' @return A list of residual trees for each level of the mainfac factor provided
#' @export 
 
generateTrees<-function(data,vars,expr,runname,...) {
  
  mainfaclevels=levels(data$mainfac)
  fits=list()
  for (i in 1:length(mainfaclevels)) {
    print(paste("Doing level ",i,": ",mainfaclevels[i],sep=""))
    try({fits[[i]]=saveTree(data,vars,expr,i,runname,mainfaclevels[i],...)})
  }
  save(file=paste("mainfactors ",runname,".Rdata",sep=""),fits)
  return(fits)
}

#' @title parseSplits
#' @description Extract information relating to the paths and volume of data in the leaves of the tree 
#' @param thistree A tree
#' @return A list of parsed splits.
#' @keywords saveTree
#' @export 
 
parseSplits<-function(thistree) {
  nodes <- as.numeric(row.names(thistree$frame))
  leafnodes=thistree$frame$var=="<leaf>"
  allpaths=rpart::path.rpart(thistree, node = nodes[leafnodes])
  allresiduals=as.matrix(thistree[[1]]$yval)[leafnodes]
  totalvol=as.matrix(thistree[[1]]$n)[1]
  allvolumes=as.matrix(thistree[[1]]$n)[leafnodes]/totalvol
  allnumbers=as.matrix(thistree[[1]]$n)[leafnodes]
  all=cbind(allpaths,allresiduals,allvolumes,allnumbers,totalvol)
  return(all)
}
#all=parsesplits(tree)
#mytree=(treesThisvar[[2]][[2]])
#pp=parsesplits(mytree)

#' @title getVarAv
#' @description This function generates a residual tree on a subset of the data
#' @param dd A dataframe
#' @param varAv A string corresponding to the numeric field to be averaged within each leaf node
#' @param varString A string
#' @return An average of the numeric variable varString in the segment
#' @keywords saveTree
#' @export 

getVarAv<-function(dd,varAv,varString)
{
  if (is.na(varString[1])) {return(0)}
  thisRule=varString[1]
  thisRule=gsub("=", "==", thisRule)
  thisRule=gsub(">==", ">=", thisRule)
  thisRule=gsub("<==", "<=", thisRule)
  thisVar=gsub( "[><=][ =].*$", "", thisRule )
  thisVals=gsub( ".*[><=][ =]", "", thisRule )
  multiflag=FALSE;
  if ( grepl(",",thisRule)) {
    multiflag=TRUE;
  }
  if (length(varString)>1) {
    
    if (multiflag==FALSE) {
      evalstr=paste("dd[dd$",thisRule,",]",sep="")
    } else {
      # handle multiple value cases
      if (inherits(dd[, thisVar], "factor")) {
        evalstr=paste("dd[which(as.character(dd$",gsub("==",") %in% as.character(c(",thisRule),"))),]",sep="",collapse="")
      } else {
        evalstr=paste("dd[which(dd$",gsub("==","==c(",thisRule),")),]",sep="",collapse="")
      }
    }
    evalstr=gsub("==,", "==\'\',", evalstr) # missing strings need values
    if (length(grep("[A-Za-z]", thisVals))>0 ) {
      thisValCommas=gsub(",","\',\'",thisVals) # insert quotes around values
      thisValCommas=paste("\'",thisValCommas,"\'",sep="",collapse="")
      evalstr=gsub(thisVals,thisValCommas,evalstr)
    }
    tryCatch({dd2=eval(parse(text=evalstr))},warning=function(){},error=function(){print(paste("Error in expectation: ",evalstr)); browser()})
    res=getVarAv(dd2,varAv,varString[2:length(varString)])
  } else { # evaluate varAv here
    if (multiflag==FALSE) {
      evalstr=paste("dd[dd$",thisRule,",]$",varAv,collapse="",sep="")  
    } else {
      if (inherits(dd[, thisVar], "factor")) { # check the class is factor
        evalstr=paste("dd[which(as.character(dd$",gsub("==",") %in% as.character(c(",thisRule),"))),]$",varAv,sep="",collapse="")
      } else {
        evalstr=paste("dd[which(dd$",gsub("==","==c(",thisRule),")),]$",varAv,sep="",collapse="")
      }
    }
    evalstr=gsub("==,", "==\'\',", evalstr) # missing strings need values
    if (length(grep("[A-Za-z]", thisVals))>0 ) {
      thisValCommas=gsub(",","\',\'",thisVals) # insert quotes around values
      thisValCommas=paste("\'",thisValCommas,"\'",sep="",collapse="")
      evalstr=gsub(thisVals,thisValCommas,evalstr)
    }
    tryCatch({res=eval(parse(text=evalstr))},warning=function(){},error=function(){print(paste("Error in expectation: ",evalstr)) })
    res=res[res!="."]
    if (inherits(res, "factor")) {res=as.numeric(levels(res))[res]}
    res=na.omit(res)
    res=mean(res)
    if (is.na(res)) {
      print(paste("Mean residual is NA - check evalstr: ", evalstr))}
  }
  return(res)
}

#' @title printResiduals
#' @description This function generates a residual tree on a subset of the data
#' @param fileConn A file connection
#' @param all A dataframe
#' @param dat The dataset
#' @param runname A string corresponding to the name of the factor being analysed
#' @param levelname A string corresponding to the factor level being analysed
#' @param treeSummaryResidualThreshold The minimum residual threshold
#' @param treeSummaryMinBucket The minumum volume per leaf
#' @param treeSummaryResidualMagnitudeThreshold Minimun residual magnitude
#' @param ... and parameters to be passed through
#' @return Residuals are printed and also saved in a simplified format.
#' @keywords saveTree
#' @export 

printResiduals<-function(fileConn,all, dat, runname, levelname,treeSummaryResidualThreshold,treeSummaryMinBucket,treeSummaryResidualMagnitudeThreshold,...) {
  for (i in 1:dim(all)[[1]]) {
    #print("\n")
    if (all[[i,2]]>treeSummaryResidualThreshold & all[[i,4]]>treeSummaryMinBucket & abs(all[[i,2]])>treeSummaryResidualMagnitudeThreshold) { 
      s1=as.character(signif(all[[i,2]], digits = 3))
      s3=as.character(signif(all[[i,3]]*100, digits = 3))
      s4=as.character(signif(all[[i,4]], digits = 3))
      s5=as.character(signif(all[[i,5]], digits = 3))
      pathterms=as.vector(all[[i,1]])[2:length(all[[i,1]])]
      s2=paste(pathterms,collapse=" and ")   
      s6=as.character(signif(getVarAv(dat,"expected",pathterms), digits = 3))
      s7=as.character(signif(getVarAv(dat,"actual",pathterms), digits = 3))
      s8=as.character(signif(getVarAv(dat,"residual",pathterms), digits = 3))
      print(i)
      print(paste("RESIDUAL:: ", runname, ":", levelname, ' :: ',s1,"(",s3,"%: ", s4," of ",s5," in tree, E=",s6,", A=",s7,", residual=",s8,") :: ",s2,collapse="",sep=""))
      writeLines(paste("RESIDUAL: ",runname, ",", levelname,",", s1,",",s3,",", s4,",",s5,",",s6,",",s7,",",s8,",",s2,collapse="",sep=""), con=fileConn)
      
      # Output may be reconfigured as desired, here are some examples:
      # writeLines(s2, con=fileConn) 
      # print(paste("RESIDUAL: ",s4,",",s6,",",s7,",",s5,",",s1,",",s2,",",s3,",",s2,collapse="",sep=""))      
      # print(paste("RESIDUAL: ",s1,",",s3,",", s4,",",s5,",",s6,",",s7,",",s8,",",s2,collapse="",sep=""))      
      # writeLines(paste("RESIDUAL: ",s1,"(",s3,"%: ", s4," of ",s5," in tree, E=",s6,", A=",s7,", residual=",s8,") :: ",s2,collapse="",sep=""), con=fileConn)
      # print(paste("RESIDUAL: ",s1,"( ",s3,"%: ", s4," of ",s5," in tree, E=",s6,", A=",s7,") :: ",s2,collapse=""))
      # writeLines(s2, con=fileConn) 
      # writeLines(paste("RESIDUAL: ",s1,"(",s3,"%: ", s4," of ",s5," in tree, E=",s6,", A=",s7," residual=",s8,") :: ",s2,collapse=""), con=fileConn)      
    }
  }
}
#printResiduals(filename="output2.txt",all)

#' @title generateResidualCutoffCode
#' @description For each tree print a summary of the significant residuals as specified by the user
#' @param data A dataframe
#' @param filename A string
#' @param trees A list of trees generated by saveTree
#' @param names A list of level names
#' @param runname A string corresponding to the name of the factor variable being analysed
#' @param ... and parameters to be passed through
#' @return A list of residuals for each tree provided.
#' @keywords saveTree
#' @export 

generateResidualCutoffCode<-function(data,filename,trees,names, runname,...){
  fileConn<-file(filename,"w")
  splitlist=list()
  for (t in 1:length(trees)){
    #writeLines("", con=fileConn) 
    writeLines(paste("Tree ",t,": ",names[t],collapse=""), con=fileConn) 
    try({tree=trees[[t]][[2]]
         splitlist[t]=list(parseSplits(tree))
         dat=data[data$mainfac==names[t],]
         # Add any other filters being applied to the data in the function saveTree here    
         printResiduals(fileConn,splitlist[t][[1]],dat, runname, names[t],...)
    })
  }
  close(fileConn)
  return(splitlist)
}
#generateResidualCutoffCode("treesTest.txt",treesProducts,mainfaclevels)

#' @title findFeatures
#' @description Perform analysis of residuals grouped by factor to identify features which explain the target variable
#' @param OutputPath A string containing the location of the input csv file. Results are also stored in this location.
#' @param fcsv A string containing the name of a csv file
#' @param ExclusionVars A string consisting of a list of variable names with double quotes around each variable
#' @param FactorToNumericList A list of variable names as strings
#' @param treeGenerationMinBucket Desired minimum number of data points per leaf (default 50)
#' @param treeSummaryMinBucket Minimum number of data points in each leaf for the summary (default 20)
#' @param treeSummaryResidualThreshold Minimum residual in the summary (default 0 for positive residuals)
#' @param treeSummaryResidualMagnitudeThreshold Minimum residual magnitude in the summary (default 0 i.e. no restriction)
#' @param doAllFactors Flag to indicate whether to analyse the levels of all factor variables (default TRUE)
#' @param maxFactorLevels (maximum number of levels per factor before it is converted to numeric (default 20)
#' @return Saves residual CART trees and associated highlighted residuals for each to the path provided.
#' @keywords findFeatures
#' @export 
#' @examples
#' 
#' require(featurefinder)
#' data(mycsv)
#' data$SMIfactor=paste("smi",as.matrix(data$SMIfactor),sep="")
#' nn=floor(length(data$DAX)/2)
#' 
#' # Can we predict the relative movement of DAX and SMI?
#' data$y=data$DAX*0
#' data$y[1:(nn-1)]=((data$DAX[2:nn])-(data$DAX[1:(nn-1)]))/
#'                   (data$DAX[1:(nn-1)])-(data$SMI[2:nn]-(data$SMI[1:(nn-1)]))/(data$SMI[1:(nn-1)])
#' 
#' thismodel=lm(formula=y ~ .,data=data)
#' expected=predict(thismodel,data)
#' actual=data$y
#' residual=actual-expected
#' data=cbind(data,expected, actual, residual)
#' 
#' # setwd('.\\test')
#' write.csv(data[(nn+1):(length(data$y)),],file='mycsv.csv',row.names=FALSE)
#' 
#' OutputPath="."
#' fcsv="mycsv.csv"
#' ExclusionVars="\"residual\",\"expected\", \"actual\",\"y\""
#' FactorToNumericList=c()
#' findFeatures(OutputPath, fcsv, ExclusionVars,FactorToNumericList,                     
#'          treeGenerationMinBucket=50,
#'          treeSummaryMinBucket=20)  
          
findFeatures <- function(OutputPath, fcsv, ExclusionVars,FactorToNumericList, 
                     treeGenerationMinBucket=50,
                     treeSummaryMinBucket=20,
                     treeSummaryResidualThreshold=0,
                     treeSummaryResidualMagnitudeThreshold=0,
                     doAllFactors=TRUE,
                     maxFactorLevels=20) {

  print(treeGenerationMinBucket)  
  print(treeSummaryMinBucket)
  print(treeSummaryResidualThreshold)
  print(treeSummaryResidualMagnitudeThreshold)
  print(doAllFactors)
  print(maxFactorLevels)
  
  currdir=getwd()
  setwd(OutputPath)
  
  dat=read.csv(paste(fcsv,sep="",collapse=""))
  
  vars=names(dat)
  eval(parse(text=paste("vars=vars[!(vars %in% c(",ExclusionVars,"))]")))
  vars=vars[]
  varsexpr=paste(vars,collapse=" + ")
  expr=as.formula(paste("residual ~", eval(varsexpr)))
  dat0=dat
  dat0$ttd_flag=1
  
  # Handle factor to numeric conversion list
  if (!is.null(FactorToNumericList)) {
    for (v in 1:length(FactorToNumericList)) {
      dat0[,FactorToNumericList[[v]]]=as.numeric(dat0[,FactorToNumericList[[v]]])
    }
  }
  
  # any factor data with more than 10 levels, convert to numeric
  nn=names(dat0)
  getclasses<-function(nname) {
    return(class(dat0[1,nname]))
  }
  dclasses=lapply(nn,FUN=getclasses)
  factorlist=which(dclasses=="factor")  
  if (length(factorlist)>0) {
    for (f in 1:length(factorlist)) {
      if (length(levels(dat0[,factorlist[f]]))>maxFactorLevels) {
        print(paste("Factor ",names(dat0[,factorlist])[f]," converted to numeric as it has ",length(levels(dat0[,factorlist[f]])), " levels.",sep="",collapse=""))
        vals=as.numeric(dat0[,factorlist[f]]) 
        vals[is.na(vals)]=0
        dat0[,factorlist[f]]=vals
      } else {
        # any factor data with fewer than 10 levels which are all numeric, retain format
        # print(class(dat[,factorlist[f]]))
      }
    }
    }
  #dat0=as.data.frame(dat0,stringsAsFactors=FALSE)          
  dat0[] <- lapply(dat0, function(x) if (is.character(x)) as.factor(x) else x)

  # An example of how to generate residual trees for all levels of a categorical variable, 
  # in this example just one level ("ALL"). The variable must be formatted as factor.
  dat0$mainfac=as.factor("ALL") # can replace with dat0$mainfac=dat0$segmentlabel, for example
  runname="ALL" # can replace with "SEGMENTLABEL"
  print('done1')

  treesAll=generateTrees(data=dat0,vars=vars,expr=expr,runname=runname, 
                         treeGenerationMinBucket=treeGenerationMinBucket,
                         treeSummaryMinBucket=treeSummaryMinBucket,
                         treeSummaryResidualThreshold=treeSummaryResidualThreshold,
                         treeSummaryResidualMagnitudeThreshold=treeSummaryResidualMagnitudeThreshold,
                         doAllFactors=doAllFactors,
                         maxFactorLevels=maxFactorLevels)

  #bestCustScoreAll=optimiseCustScore("ALL")
  #plot(x=cutofflist,bestCustScoreAll[[1]],xlab='Cutoff',ylab='Mean positive residuals over entire dataset')
  mainfaclevels=levels(dat0$mainfac)
  generateResidualCutoffCode(data=dat0,"treesAll.txt",treesAll,mainfaclevels, runname,
                             treeGenerationMinBucket=treeGenerationMinBucket,
                             treeSummaryMinBucket=treeSummaryMinBucket,
                             treeSummaryResidualThreshold=treeSummaryResidualThreshold,
                             treeSummaryResidualMagnitudeThreshold=treeSummaryResidualMagnitudeThreshold,
                             doAllFactors=doAllFactors,
                             maxFactorLevels=maxFactorLevels)
  
  # Run this section if you want to generate residuals for every level of every categorical variable  
  if (doAllFactors) {
    # Loop over all categorical variables. These can include user-defined variables based on numerical variables 
    dir.create(file.path(getwd(), "allfactors"), showWarnings = FALSE)
    setwd("./allfactors")
    catlist = lapply(dat0[,vars],FUN=class)=="factor"
    lens=lapply(dat0[,vars],FUN=function(x){length(levels(x))})
    
    lengthlist = lens>=2 & lens<=maxFactorLevels # only plot factors with at most
    faclist = vars[catlist & lengthlist]
    lenlist = lens[catlist & lengthlist]
    length(lenlist)
    length(faclist)
    if (length(faclist)>0) {
      for (v in 1:length(faclist)) {
        thisvar=faclist[v]
        dat0$mainfac=dat0[,thisvar]
        mainfaclevels=levels(dat0$mainfac)
        print(paste("Variable ", v, ": ", thisvar, " has levels:", sep=""))
        print(mainfaclevels)
        runname=thisvar
        treesThisvar=generateTrees(data=dat0,vars,expr,runname,
                                   treeGenerationMinBucket=treeGenerationMinBucket,
                                   treeSummaryMinBucket=treeSummaryMinBucket,
                                   treeSummaryResidualThreshold=treeSummaryResidualThreshold,
                                   treeSummaryResidualMagnitudeThreshold=treeSummaryResidualMagnitudeThreshold,
                                   doAllFactors=doAllFactors,
                                   maxFactorLevels=maxFactorLevels)
        generateResidualCutoffCode(data=dat0,paste("trees",thisvar,".txt",sep=""),treesThisvar,mainfaclevels, runname, 
                                   treeGenerationMinBucket=treeGenerationMinBucket,
                                   treeSummaryMinBucket=treeSummaryMinBucket,
                                   treeSummaryResidualThreshold=treeSummaryResidualThreshold,
                                   treeSummaryResidualMagnitudeThreshold=treeSummaryResidualMagnitudeThreshold,
                                   doAllFactors=doAllFactors,
                                   maxFactorLevels=maxFactorLevels)
      }  
    }
    setwd("..")
  }
  setwd(currdir)
}
#findFeatures(OutputPath, fcsv, ExclusionVars)

#' @title data
#' @description Sample data based on dataset EuStockMarkets in the datasets package. 
#' @docType data
#' @keywords sample data
#' @name data
#' @author Richard Davis \email{richard.davis@cba.com.au}
#' @format A data frame with 1860 rows and 4 variables
#' @source \url{https://stat.ethz.ch/R-manual/R-devel/library/datasets/html/00Index.html}
#' @examples
#' data(mycsv)
#' thismodel=lm(formula=DAX ~ .,data=data)
#' expectedprob=predict(thismodel,data)
#' actualprob=data$DAX
#' residual=actualprob-expectedprob
#' data=cbind(data,expectedprob, actualprob, residual)
NULL
