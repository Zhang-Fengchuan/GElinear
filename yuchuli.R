
name_x <- c('IID','near','out');name_x###############输入规范
name_y <- c('WEINO','AGE','AGE_1','AGE_2','GENDER','RA','RA_1','RA_2')#########输入规范
number_repeat###########堆叠次数
response_variable <-'RA' 
name_year <- 'AGE'
rangesnp <- c(9:33)###snp表中范围
c('全一','gender','year','near','out')
range_zscore <- c(3,4,5)

yuchuli <- function(position_snp,position_phenotype,name_x,name_y,number_repeat,response_variable,name_year,rangesnp,range_zscore,grouping){
  library(xlsx)#加载xlsx包
  library(readxl)
  data_biaotou_x <- read_excel(position_snp,col_names = F,range = cell_rows(1));data_biaotou_x###########默认第一行为表头名称
  data_biaotou_x <- as.matrix(data_biaotou_x);data_biaotou_x
  lieshu_x <- matrix(0,length(name_x),1);lieshu_x#存储选取的自变量所在的位置索引
  for(i in c(1:length(lieshu_x))){
    lieshu_x[i,1] <- which(data_biaotou_x==name_x[i], arr.ind = TRUE)[2]
  }
  
  data_biaotou_y <- read_excel(position_phenotype,col_names = F,range = cell_rows(1));data_biaotou_y###########默认第一行为表头名称
  data_biaotou_y <- as.matrix(data_biaotou_y);data_biaotou_y
  data_biaotou_y==name_y[1]
  lieshu_y <- matrix(0,length(name_y),1);lieshu_y#存储选取的自变量所在的位置索引
  for(i in c(1:length(lieshu_y))){
    lieshu_y[i,1] <- which(data_biaotou_y==name_y[i], arr.ind = TRUE)[2]
  }
  lieshu_x <- as.integer(rbind(lieshu_x,as.matrix(rangesnp)));lieshu_x;#将snp与协变量索引进行结合
  lieshu_y <- as.integer(lieshu_y);lieshu_y#将协变量与响应变量索引结合
  
  
  data_snp<-as.matrix(read_excel(position_snp,sheet = "Sheet1"));data_snp#读取全部snp数据
  data_phenotype <- as.matrix(read_excel(position_phenotype,sheet = "Sheet1"));data_phenotype#读取全部表型加部分协变量数据
  data_snp_only <- data_snp[,lieshu_x];data_snp_only#仅包括本次实验的数据
  data_phenotype_only <- data_phenotype[,lieshu_y];data_phenotype_only#仅包括本次实验的数据
  
  ############################################################以上为数据导入，以下为数据预处理
  
  
  index <- matrix(0,length(data_snp_only[,1]),1);index#第一种预处理方案与剔除snp缺失值后的数据对齐
  for(i in c(1:length(index))){
    index[i,1] <- which(data_phenotype_only[,1]==data_snp_only[i,1],arr.ind = TRUE);
  }
  data_phenotype_only_duiqi <- apply(data_phenotype_only[index,],2,as.numeric);data_phenotype_only_duiqi#转为数值型矩阵
  data_snp_only_duiqi <- apply(data_snp_only,2,as.numeric);data_snp_only_duiqi#转为数值型矩阵
  
  
  number_na_1 <- as.matrix(which(is.na(data_phenotype_only_duiqi),arr.ind = TRUE));number_na_1#检测NA的个数(缺失值位置索引矩阵)
  number_na_2 <- as.matrix(which(is.na(data_snp_only_duiqi),arr.ind = TRUE));number_na_2#检测NA的个数(缺失值位置索引矩阵)
  na_matrix_1 <- matrix(0,1,(length(name_x)-1));na_matrix_1#每个变量的缺失值情况（缺失统计情况矩阵）
  na_matrix_2 <- matrix(0,1,(length(name_y)-1));na_matrix_2#每个变量的缺失值情况（缺失统计情况矩阵）
  for(i in c(1:length(na_matrix_1))){
    na_matrix_1[1,i] <- (length(data_snp_only_duiqi[,(i+1)])-length(na.omit(data_snp_only_duiqi[,(i+1)])))/length(data_snp_only_duiqi[,(i+1)]);#x部分协变量缺失值统计
  }
  for(i in c(1:length(na_matrix_2))){
    na_matrix_2[1,i] <- (length(data_phenotype_only_duiqi[,(i+1)])-length(na.omit(data_phenotype_only_duiqi[,(i+1)])))/length(data_phenotype_only_duiqi[,(i+1)]);#x部分协变量缺失
  }
  na_matrix <- data.frame(rbind(cbind(t(as.matrix(name_x)[-1]),t(as.matrix(name_y)[-1])),cbind(na_matrix_1,na_matrix_2)));na_matrix#######缺失值统计情况
  
  
  
  
  for(i in c(1:nrow(number_na_1))){
    data_phenotype_only_duiqi[number_na_1[i,1],number_na_1[i,2]] <- mean(na.omit(data_phenotype_only_duiqi[,(number_na_1[i,2])]))#去除缺失值后的均值
  }
  if(nrow(number_na_2)>0){
    for(i in c(1:nrow(number_na_2))){
      data_snp_only_duiqi[number_na_2[i,1],number_na_2[i,2]] <- mean(na.omit(data_snp_only_duiqi[,(number_na_2[i,2])]))#去除缺失值后的均值
    }
  }
  ########################################################以上按照需求选取对应的列，统计缺失值，并用均值代替
  #na_matrix 
  #data_phenotype_only_duiqi
  #data_snp_only_duiqi
  
  ##############################以下根据因变量在表几 构建位置索引并获得数据矩阵
  panduan_source <- 0;
  for (i in c(1:length(name_y))){
    if(name_y[i]==response_variable){
      panduan_source <- panduan_source+1;
    }
  }
  if(panduan_source>=1){
    index_y <- seq(which(name_y==response_variable,arr.ind = TRUE),by=1,length.out=number_repeat);index_y#
    index_year <- seq(which(name_y==name_year,arr.ind = TRUE),by=1,length.out=number_repeat);index_year#
    index_x_1 <- seq(2,by=1,length.out= (length(name_x)-1));index_x_1
    index_x_2 <- seq(2,by=1,length.out=(length(name_y)-1));index_x_2
    for(i in c(1:length(index_y))){
      index_x_2 <- index_x_2[-which(index_x_2==index_y[i],arr.ind = TRUE)];index_x_2
    }
    for(i in c(1:length(index_year))){
      index_x_2 <- index_x_2[-which(index_x_2==index_year[i],arr.ind = TRUE)];index_x_2
    }
    index_x_2
    index_w <- seq((length(name_x)+1),by=1,length.out=length(rangesnp));index_w#
    c <- number_repeat*nrow(data_snp_only_duiqi);c#多年数据堆叠维度
    x <- matrix(1,c,1);#截距项x值全一
    for (i in c(1:length(index_x_2))) {
      x <- cbind(x,as.matrix(rep(as.matrix(data_phenotype_only_duiqi[,index_x_2[i]]),number_repeat)))
    }
    x <- cbind(x,as.matrix(as.vector(data_phenotype_only_duiqi[,index_year])))
    for (i in c(1:length(index_x_1))) {
      x <- cbind(x,as.matrix(rep(as.matrix(data_snp_only_duiqi[,index_x_1[i]]),number_repeat)))
    }
    x[(x[,2]==1),2]=0#性别0-1二分类哑变量处理
    x[(x[,2]==2),2]=1#性别0-1二分类哑变量处理
    x <- apply(x, 2, as.numeric);x #处理成功之后的x
    if (grouping=='yes'){
      for (i in c(4:5)) {
        x[(x[,4]>median(x[,4])),4]=1;
        x[(x[,4]<=median(x[,4])),4]=0;
        x[(x[,5]>median(x[,5])),5]=1;
        x[(x[,5]<=median(x[,5])),5]=0;
      }
    }
    y <- apply(as.matrix(as.vector(data_phenotype_only_duiqi[,index_y])),2,as.numeric);y#处理成功之后的y
    w_org <- matrix(,c,1);#创建w空矩阵
    for (i in c(1:length(index_w))) {
      w_org <- cbind(w_org,rep(as.matrix(data_snp_only_duiqi[,index_w[i]]),number_repeat));w_org
    }
    w_org <-as.matrix(w_org[,-1]);w_org#处理好的原始的存储snp的矩阵
    w <- apply(matrix(0,c,(length(rangesnp)*ncol(x))),2,as.numeric);w
    for (i in c(1:length(rangesnp))) {
      for (j in c(1:ncol(x))) {
        w[,(((i-1)*ncol(x))+j)] <- as.matrix(w_org[,i]*x[,j])
      }
    }
    w #处理成功之后的w 
  }else{
    index_y <- seq(which(name_x==response_variable,arr.ind = TRUE),by=1,length.out=number_repeat);index_y
    index_x_1 <- seq(2,by=1,length.out=(length(name_x)-1));index_x_1
    for (i in c(1:length(index_y))) {
      index_x_1 <- index_x_1[-which(index_x_1==index_y[i],arr.ind = TRUE)];index_x_1
    }
    index_x_1
    index_x_2 <- seq(2,by=1,length.out=(length(name_y)-1));index_x_2
    index_w <- seq((length(name_x)+1),by=1,length.out(length(rangesnp)));index_w
    c <- number_repeat*nrow(data_snp_only_duiqi);c#多年数据堆叠维度
    x <- matrix(1,c,1);#截距项x值全一
    for (i in c(1:length(index_x_2))) {
      x <- cbind(x,as.matrix(rep(as.matrix(data_phenotype_only_duiqi[,index_x_2[i]]),number_repeat)))
    }
    x <- cbind(x,as.matrix(as.vector(data_phenotype_only_duiqi[,index_year])))
    for (i in c(1:length(index_x_1))) {
      x <- cbind(x,as.matrix(rep(as.matrix(data_snp_only_duiqi[,index_x_1[i]]),number_repeat)))
    }
    x[(x[,2]==1),2]=0#性别0-1二分类哑变量处理
    x[(x[,2]==2),2]=1#性别0-1二分类哑变量处理
    x <- apply(x, 2, as.numeric);x #处理成功之后的x
    y <- apply(as.matrix(as.vector(data_snp_only_duiqi[,index_y])),2,as.numeric);y#处理成功之后的y
    w_org <- matrix(,c,1);#创建w空矩阵
    for (i in c(1:length(index_w))) {
      w_org <- cbind(w_org,rep(as.matrix(data_snp_only_duiqi[,index_w[i]]),number_repeat));w_org
    }
    w_org <- w_org[,-1];w_org#处理好的原始的存储snp的矩阵
    w <- apply(matrix(0,c,(length(rangesnp)*ncol(x))),2,as.numeric);w
    for (i in c(1:length(rangesnp))) {
      for (j in c(1:ncol(x))) {
        w[,(((i-1)*ncol(x))+j)] <- as.matrix(w_org[,i]*x[,j])
      }
    }
    w #处理成功之后的w
    
  }
  
  ############################以上数据构建完成，以下归一化r (i in c(1:length(rangesnp))) {
  if(grouping=='yes'){
   range_zscore=range_zscore[1] 
  }
  for (j in range_zscore) {
    w[,(((i-1)*ncol(x))+j)] <- (w[,(((i-1)*ncol(x))+j)]-mean(w[,(((i-1)*ncol(x))+j)]))/sd(w[,(((i-1)*ncol(x))+j)])
  }

  

list_data <- list(x,w,y,na_matrix);list_data
names(list_data) <- c("x","w","y","na_matrix")
return(list_data)
}
