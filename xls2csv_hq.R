xls2csv_hq <- function(filenames,wpath=""){
        library(xlsx)
        options(warn=-1)
        
        for(i in 1:length(filenames)){
                print(i)
                aa <- read.xlsx2(filenames[i],1,as.data.frame = TRUE, header=FALSE, colClasses="character")
                tmp <- as.Date(as.numeric(aa[,1]),origin="1899-12-30")
                tmp1 <- grepl("-",aa[,1])
                
                if( !all(is.na(tmp)) | any(tmp1) ){
                        subhead <- which(aa[,1]=="ָ������")
                        cols <- aa[subhead, ]
                        subref <-  which(grepl("������Դ",aa[,1]))
                        aa <- aa[-subref, ]
                        
                        if( !any(grepl("-",aa[,1])) ){
                                tmp <- as.Date(as.numeric(aa[,1]),origin="1899-12-30")
                                aa <- aa[!is.na(tmp), ]
                                aa[,1] <- tmp[!is.na(tmp)]
                        }else{
                                aa <- aa[grepl("-",aa[,1]), ]
                        }
                        colnames(aa) <- cols
                        
                        if(wpath==""){
                                write.csv(aa,file=gsub(".xls",".csv",filenames[i]),quote=FALSE,row.names = FALSE)     
                        }else{
                                write.csv(aa,file=paste(wpath,"/", gsub(".xls",".csv",basename(filenames[i])),sep=""),quote=FALSE,row.names = FALSE)      
                        }
                }else{
                        print(paste(filenames[i],": No validate records.",sep=""))
                }
        }
        options(warn=0)

}

test <- function(){
        source("xls2csv_hq.R")
        path <- "D:/data/����/���ݹ���/PX��������/��һ�׶�"
        filenames <- list.files(path,pattern=".xls",recursive = TRUE,full.names = TRUE)
        filenames <- filenames[1:56]
        xls2csv_hq(filenames,wpath="D:/data/����/���ݹ���/PX��������/��һ�׶�/tmp")
}
