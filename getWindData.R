getWindData <- function(indexs='S0031529,S0070127',sTime='2015-07-13',eTime='2016-07-12',filename=NULL){

        # wind 实时数据获取, 以天为基础
        library(WindR)
        w.start()
        w_edb_data<-w.edb(indexs,sTime,eTime,'')
        newData <- w_edb_data$Data
        if(!is.null(filename)) write.csv(newData,file=paste(filename,day,".csv",sep=""),row.names=FALSE,quote=FALSE)
        w.stop()
        
        newData
}