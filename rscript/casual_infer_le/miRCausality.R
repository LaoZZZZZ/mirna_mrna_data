#Inferring miRNA-mRNA causal regulatory relationships from expression data.
#Thuc Duy Le, Lin Liu, Anna Tsykin, Gregory J Goodall, Bing Liu and Jiuyong Li.
#
#Created by Thuc Duy Le
#March 2012
#1.Save miRCausality.R into your local machine
#2.Follow the steps in the RunningCausality file to run the program.

#Main function
#This function is for testing purpose.
#param: dataset (.csv), tuning parameter for PC algorithm.
	run<-function(dataset,alp){	
		data<-Read(dataset)	
		stdData<-Standardise(data)
		caef<-Causaleffects(stdData,alp, nmiRs)
		caef
	}
#Bootstrapping
#This function is to bootstrap the data nbootstr times, and estimate the causal effects between the miRNAs and mRNAs accordingly.
#Param: stdData: the standardised dataset. nbootstr: the number of bootstraps. alp: the tuning parameter for the PC algorithm.
#Param: nmiRs: the number of miRNAs in the dataset. 
#return: bstrResult. The result is the array that contains nbootstr table. Each table represents the miRNAs-mRNAs causal effects.
bootstr<-function(stdData, nbootstr, alp, nmiRs){
	bstrResult<-lapply(1:nbootstr, function(i) Causaleffects(stdData[sample(nrow(stdData), nrow(stdData), replace=T), ], alp, nmiRs))
	return(bstrResult)
	
}

#Read  the dataset with .csv format, the first column should be the #rownames (sample names), the first row is the column names (gene #names).
#param: filename (.csv)
#return: data table
	Read<-function(dataset){
		data<-read.csv(dataset, header=TRUE, sep=",")
		return(data)
	}
	
	
#Standardise the dataset to have unit standard deviation.
#param: data table
#return: standardised data set
	Standardise<-function(data){
		ncol<-ncol(data)
		nrow<-nrow(data)
		stdData<-matrix(nrow=nrow,ncol=ncol-1)
		rownames(stdData)<-data[,1]
		for (i in 2:ncol){
			stdData[,i-1]<-scale(data[i],center=TRUE, scale=TRUE)
			
			}
		return(stdData)
	}

#Learning PCDAG from data then inferring causal effects
#required: pcalg package
#param: standardised dataset, turning parameter for PC algorithm, number of nodes need to
#intervene (number of miRNAs)
#return: The causal effect values between miRNAs and mRNAs.
	Causaleffects<-function(stdData,alp, nmiRs){
		library(pcalg)
		#Learning PCDAG
		suffStat<-list(C=cor(stdData), n=nrow(stdData))  
		pc.fit<-pc(suffStat, indepTest=gaussCItest, p=ncol(stdData), alpha=alp)	
		nmRs<-ncol(stdData)-nmiRs
		result<-matrix(nrow=nmRs, ncol=nmiRs)
		for (l in 1:nmiRs){
				
			#Inferring causal effects
			caef<-idaFast(l,(nmiRs+1):ncol(stdData),cov(stdData), pc.fit@graph )
		
			#min of absolute values.
			caef1<-matrix(nrow=nmRs,ncol=1)
			for (k in 1:nmRs){
				caefabs<-abs(caef)
				index<-which(caefabs==min(caefabs[k,]), arr.ind=TRUE)
				pos<-index[1,2]
				caef1[k,]<-caef[k,pos]
			}
			result[,l]<-caef1
		}
		temp<-"Finish!"
		print(temp)
		return(result)
	}
	