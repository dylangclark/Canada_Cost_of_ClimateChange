
#####Lasso climate model analysis Canada
#This script is used to: 1) pull all 24 statistically downscaled GCMs into R (netCDF format); 2) cut the GCM into a region (select region option);
# 3) calculate region mean temp and mean prcp for each period and RCP; 4) plot the mean temp vs. mean prcp for each period and RCP; 5) create a convex hull to select the outermost points on each plot


#Statistically downscaled GCMs have been downloaded from Environment and Climate Change Canada http://climate-scenarios.canada.ca/?page=bccaqv2-data

##Load libraries
library(RNetCDF)
library(raster)
library(ncdf4)
library(chron)
library(lattice)
library(RColorBrewer)
library(rgdal)
library(dplyr)
library(rgeos)
library(purrr)
library(ggplot2)


######## Pull all 24 netCDFs into R




# read netCDFs [needs to be a loop]
netCDFdirectory<-"/Users/cccgi/Documents/GIS/GCMs/"

# read shapefile polygon for reagional boundaries
AdminBound<-readOGR(dsn=boundaryFolder, layer = boundaryName)
boundaryFolder<-"/Users/cccgi/Documents/GIS/CanMap/GCM AnalysisLayer/"
boundaryName<-"SimpleAdminArea"

####set netcdf variable
varName<-"tg_mean"
#OR
varName<-"prcptot"
##


#######Function. Install this function and the one below.

## Remember to change the varName file between prcp and temp.

ClimatePlots<-function(varPath,varName,boundaryFolder,boundaryName,netCDFdirectory){
  
  ###create table of all statistically downscaled GCMs of the same variable type and note the classes
  GCMname<-c("ACCESS1-0","bcc-csm1-1","bcc-csm1-1-m","NorESM1-M","NorESM1-ME","MRI-CGCM3","MPI-ESM-MR","MPI-ESM-LR","MIROC5","MIROC-ESM","MIROC-ESM-CHEM","IPSL-CM5A-MR",
             "IPSL-CM5A-LR","HadGEM2-ES","HadGEM2-AO","GFDL-ESM2M","GFDL-ESM2G","GFDL-CM3","FGOALS-g2","CanESM2","CSIRO-Mk3-6-0","CNRM-CM5","CESM1-CAM5","CCSM4","BNU-ESM")
  Eras<-c("1971-2000", "2041-2070","2071-2100")
  RCP<-c("rcp45","rcp85")
  PT<-c("BC","AB","YT","SK","QC","PE","ON","NU","NS","NT","NL","NB","MB")
  counter<-0
  
  gcmTable<-GCMTableFun(gcmname=GCMname,era=Eras,rcp=RCP,var=varName)
  runs<-13*(nrow(gcmTable))
  print("wrote GCMTable!")
  
  for (i in PT){
    pt<-i
    
    print(paste0("working on ",pt))
    ###Read shapefile that will be used to cut GCMs for regional analysis
    selectPT<-subset(AdminBound, PROV==pt)
    selectPT<-gSimplify(selectPT,tol=0.001)
    
    ###Iterate through GCMs
    n<-1:nrow(gcmTable)
    for (i in n){
      counter<- counter + 1
      y<-gcmTable[i,]
      gcm<-y[,2]
      Rcp<-y[,4]
      era<-y[,3]
      
      if(era=="1971-2000"){
        print("baseline era found, skipping row")
        print(paste0("Run number ",counter, " of ",runs))
        
      }else if (era=="2041-2070"|era=="2071-2100"){
        GCMpath<-as.character(paste0(netCDFdirectory,(y[,5])))
        rastGCM<-raster(GCMpath,var=varName, snap='in')
        Rast<-crop(rastGCM, selectPT)
        check<- data.frame(x.mean=cellStats(Rast, "mean"))
        if (check<0 | check>0){}else if(check==0){errorCondition(paste0("Error, Raster has mean value of 0 ", gcm," ",Rcp," ",era))}
        
          ####Find the correct baseline GCM
        df<-gcmTable %>% 
          dplyr::select(var, model, era, rcp, url) %>%
          filter(var==varName & model==gcm & rcp==Rcp & era=="1971-2000")
        
        urlBase<-df[,5]
        urlBase<-as.character(paste0(netCDFdirectory,urlBase))
        rastbase<-raster(urlBase,var=varName, snap='in')
        RastBase<-crop(rastbase,selectPT)
        check<- data.frame(x.mean=cellStats(RastBase, "mean"))
        if (check<0 | check>0){}else if(check==0){errorCondition(paste0("Error, Raster base has mean value of 0 ", gcm," ",Rcp," ",era))}
        
        ##Rast math (Rast - RastBase)
        if (varName=="prcptot"){Delta<-overlay(Rast,RastBase,fun=function(r1,r2){return(r1/r2)})}
        else if (varName=="tg_mean"){Delta<-overlay(Rast,RastBase,fun=function(r1,r2){return(r1-r2)})}
        
        m<- data.frame(x.mean=cellStats(Delta, "mean"))
        sd<-data.frame(x.sd=cellStats(Delta, "sd"))
        data<-cbind(pt,varName,gcm,Rcp,era,m,sd)
        
        ##write data
        if(exists("table_iterate")){table_iterate<-rbind(table_iterate,data)}else{table_iterate<-data}
        
        
        
        ##report status of run
        print(paste0("Run number ",counter, " of ",runs))
        
        
        ##End else if loop
      }
      
      ##END spatial analysis for loop 
    }
    
    ##End PT for loop
  }
  
  Out<-return(table_iterate)
  ##END function 
}





#######################################################################################################################################GCM Table function
#######################Function for file lookup

GCMTableFun<-function(gcmname,era,rcp,var){
  for(i in gcmname){
    gcm<-i
    for (i in era){
      time<-i
      tbl<-lapply(rcp, function(x){
        URL<-as.character((paste0("BCCAQv2+ANUSPLIN300_",gcm,"_",x,"_",var,"_",time,"_mean.nc")))
        paste0(c(var,gcm,time,x,URL))
        
      })
      
      if(exists("GCMtab")){GCMtab<-rbind(tbl,GCMtab)}else{GCMtab<-tbl}
      
    }
  }
  GCMtab<-do.call(rbind.data.frame, GCMtab)
  colnames(GCMtab) <- c("var", "model","era","rcp","url")
  return(GCMtab)
}




#################################################################################################   Run model with this:
Out<-ClimatePlots(varPath = varPath,varName=varName,boundaryFolder = boundaryFolder,boundaryName = boundaryName,netCDFdirectory=netCDFdirectory)

write.csv(Out,file="/Users/cccgi/Documents/GIS/temp.csv")
#Or
write.csv(Out,file="/Users/cccgi/Documents/GIS/prcp.csv")




#####Analysis of the data

Data<-read.csv(file="/Users/cccgi/Documents/GIS/DataRaw.csv",header=T)

PT<-levels(Data$pt)
RCP<-levels(Data$Rcp)
ERA<-levels(Data$era)

for (i in PT){
  pT<-i
  for (i in RCP){
    rcP<-i
      for (i in ERA){
        erA<-i
        
        pull<-Data %>% 
          dplyr::select(pt, gcm, Rcp, era, Temp.mean, Temp.sd, Prcp.mean, Prcp.sd) %>%
          filter(pt==pT & Rcp==rcP & era==erA)
          gcm<-pull$gcm
          pull$x<-pull$Temp.mean
          pull$y<-((pull$Prcp.mean)-1)*100
          index<-chull(pull$x,pull$y)
          boundaryGCM<-(paste(as.character(gcm[index]), sep=", ", collapse=", "))
          Title<-paste0(rcP, " ", pT, " ", erA)
          ###Save output
          save<-cbind(pT,rcP,erA,boundaryGCM)
          if(exists("chullOut")){chullOut<-rbind(save,chullOut)}else{chullOut<-save}
          
          
          ###Draw diagram
          rownames(pull)<-pull$gcm
          xAxis<-"temp change (deg C)"
          yAxis<-"% precp change"
          rownames(pull) <- pull[,2]
          hpts <- chull(pull$x, pull$y)
          hpts <- c(hpts, hpts[1])
          tempPlot<-ggplot(d=pull, aes(x,y))+ geom_point(shape=1,color="blue") +
            geom_text(aes(x,y,label=rownames(pull)),
                      position = position_dodge(width=0.5), size=2.5)+
            geom_path(data=pull[hpts,], color="blue")+
            labs(title = Title, x=xAxis,y=yAxis)
        
          ggsave(tempPlot, file=paste0("/Users/cccgi/Documents/GIS/plots/", Title,".png"),
                 width = 14, height = 10, units = "cm")
          
    }
  }
}   


write.csv(chullOut,file="/Users/cccgi/Documents/GIS/boundaries.csv")

remove(chullOut)







