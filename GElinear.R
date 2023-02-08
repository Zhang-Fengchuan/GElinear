


#----------------------------------------------------------------------------------------------------------------------------------------------
#-------------------Gene and environment interaction linear regression ------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------
GElinear <-  function(x,w,y,lambda,geshu,kesai,start,
                      principle,catagory,response_variable,
                      site,leftpart,left,right,yinzi_left,
                      yinzi_right,eps1,eps2,eps3,max_iter1,max_iter2){
  
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## The name of the function: GElinear
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Description: 
  ##            Performing joint analysis based on linear regression with the sparse group penalty.
  ##            The regression coefficients under the given tuning parameters can be calculated, 
  ##            and the optimal regression coefficient is selected under the given criterion such as BIC or AIC.
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Required preceding functions or packages:
  ##            R functions: model()   qujian()  
  ##            Rcpp functions: gmcp_cpp()
  ##            R packages: Rcpp
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Input:
  ## @ y: n * 1 vector, the response variable
  ## @ x: n * (q+1) matrix, the design matrix corresponding environmental variables.
  ## @ w: n * ((q+1)*p) matrix, the design matrix corresponding gene variables and its interactive effects.
  ## @ kesai: a float value, regularization parameter in MCP, the default setting is 3.
  ## @ max_iter1: int, Maximum number of cycles of the algorithm, the default setting is 200.
  ## @ max_iter2: int, Maximum number of cycles of the algorithm, the default setting is 50.
  ## @ start: the type of starting value, which can be selected from "lassostart" and "olsstart" and "warmstart".
  ## @ geshu: a value (s) the number of the tuning parameters
  ## @ lambda: s * 2 matrix, the tuning parameters,!!!!!!!!(if you choose "lassostart" or "olsstart", you don't need to provide the lambda) 
  ##         the s rows represent s choices of the tuning parameters,
  ##         the 2 columns represent lambda1, lambda2,
  ##         (in this paper, to simplify, lambda1 = lambda2 is set, but this function is written in such a way that the two can be not equal.)
  ## @ eps1: a float value, algorithm termination threshold, the default setting is 1e-8
  ## @ eps2: a float value, algorithm termination threshold, the default setting is 1e-3
  ## @ eps3: a float value, algorithm termination threshold, the default setting is 1e-5
  ## @ principle: the parameter selection criterion, which can be selected from "AIC" and "BIC".
  ## @ catagory: the name matrix of covariates, such as [catagory <- as.matrix(c('','gender','year','near','out'))]
  ## @ response_variable: the name of the research response variable, such as ('AL') 
  ## @ site: the address of the saved consequence such as ('C:/Users/Desktop/consequence/')
  ## @ leftpart:
  ## @ left:
  ## @ right:
  ## @ yinzi_left:
  ## @ yinzi_right:
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Output:
  ## A list including:
  ## @ consequence: ( (q+1)*(p+1) ) * s vector, all regression coefficients corresponding s choices of given tuning parameters.
  ## @ BIC: the number of BIC values corresponding s choice of given tuning parameters.
  ## @ jieguo_best: the estimated regression coefficients under the principle of BIC or AIC 
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  res=list()
  #####################################################################################参数设置
  yanzhou <- y;yanzhou
  A<-cbind(x,w);A########得自变量的数据矩阵
  n<-length(yanzhou);n
  d<-ncol(x);d
  rangesnp <- c(1:(ncol(w)/d));rangesnp
  range_snp<-c(1:length(rangesnp));range_snp#25个snp,循环器
  
  #left和right可取同一个值
  #leftpart控制总系数的比例
  #yinzi是惩罚因子变动系数(0-1)
  #geshu是设置的惩罚因子个数
  #a1、b1为初值(最小二乘解);
  #kesai为常数; #n为样本量;
  #range_snp为snp范围的数组(例如c(1:25));
  #d为每个bj的维度,即一个snp与其交互项的维度和
  #x,w分别为协变量和经过处理的snp交互项数据矩阵
  #yanzhou为响应变量的数据向量
  #eps1为求解每个bj的收敛准则
  #eps2为求解所有的bj的收敛准则(每坐标下降所有snp一遍后,前后两步的系数差范数)
  #eps3为交替求解a,b时前后两步的系数差范数
  #max_iter1为针对一个当前的a,求解b时坐标下降的最大迭代步数
  #max_iter2为a,b交替求解时，最大交替迭代步数
  
  if (start=='lassostart'){
    ####################################################################################lasso初值
    library(caret);library(glmnet);library(corrplot)
    library(ggplot2)
    lambdas <- seq(0,1, length.out = 1000)####################################初值估计用lasso回归，设定惩罚因子的变化区间
    set.seed(1245)#设定种子
    lasso_model <- cv.glmnet(A,yanzhou,alpha = 1,lambda = lambdas,nfolds =5)#用交叉验证进行lasso回归选取最优因子lambda
    plot(lasso_model)#随最优因子的变化mse的变化情况
    plot(lasso_model$glmnet.fit, "lambda", label = T)#画出lasso的系数变化图
    lasso_min <- lasso_model$lambda.min;lasso_min #最优的惩罚因子
    lasso_best <- glmnet(A,yanzhou,alpha = 1,lambda = lasso_min)#利用最优因子进行lasso回归
    coef(lasso_best)#得到lasso估计的系数作为初值(提高速度)
    canshu_chushi<-coef(lasso_best);canshu_chushi#定义参数初始值变量
    b1<-as.matrix(canshu_chushi[(d+1+1):length(canshu_chushi),1]);b1#
    a1<-as.matrix(canshu_chushi[2:(d+1),1]);a1#从截距项后面一项开始取d个
    a1[1]<-canshu_chushi[1,1];a1#将a1的第一项改为截距项
  }else if (start=='olsstart'){
    #####################################################################################最小二乘估计初值
    canshu_chushi<-solve(t(A)%*%A)%*%t(A)%*%yanzhou;canshu_chushi
    #b1<-as.matrix(canshu_chushi[(d+1):length(canshu_chushi),1]);b1
    #canshu_chushi[(d+1):length(canshu_chushi),1] <- 0;canshu_chushi##############
    b1<-as.matrix(canshu_chushi[(d+1):length(canshu_chushi),1]);b1###############
    a1<-as.matrix(canshu_chushi[1:d,1]);a1
    a1 <- solve(t(x)%*%x)%*%t(x)%*%(yanzhou-w%*%b1);a1
  }else if (start=='warmstart'){
    #####################################################################################最小二乘估计初值
    canshu_chushi<-solve(t(A)%*%A)%*%t(A)%*%yanzhou;canshu_chushi
    #b1<-as.matrix(canshu_chushi[(d+1):length(canshu_chushi),1]);b1
    #canshu_chushi[(d+1):length(canshu_chushi),1] <- 0;canshu_chushi##############
    b1<-as.matrix(canshu_chushi[(d+1):length(canshu_chushi),1]);b1###############
    a1<-as.matrix(canshu_chushi[1:d,1]);a1
    a1 <- solve(t(x)%*%x)%*%t(x)%*%(yanzhou-w%*%b1);a1
  }
  
  
  if (start=='warmstart'){
    lamda1 <- lambda[,1];lamda1
    lamda2 <- lambda[,2];lamda2
    lamda_size<-length(lamda2);lamda_size #惩罚因子的个数
    BIC<-matrix(Inf,lamda_size,1);BIC #贝叶斯信息准则的取值，用来选取最优惩罚因子lambda2以及lambda1
    #z<-array(Inf,dim=c(3+length(a1)+length(b1)+3,lamda_size));z#存储每次估计的系数矩阵,1-3三个参数值(lamda1,lamda2,kesai),4a,5b,6time,7panding,8iter
    z <- list();z#设置储存空列表
  }else{
    ########################################################################################搜寻惩罚因子区间
    s <- Sys.time()#计时开始
    data_qujian <- qujian(left,right,leftpart,yinzi_left,yinzi_right,geshu,a1,b1,kesai,n
                          ,range_snp,d,x,w,yanzhou,eps2,eps1,max_iter1,max_iter2,eps3);data_qujian
    e <- Sys.time()
    print(e-s)#计时结束
    #################################################################################################存储矩阵构建
    lamda1 <- seq(data_qujian$left_end,data_qujian$right_end,length.out = data_qujian$geshu);lamda1 #设定惩罚因子lambda1的变化区间
    lamda2 <- seq(data_qujian$left_end,data_qujian$right_end,length.out = data_qujian$geshu);lamda2 #设定惩罚因子lambda2的变化区间
    lamda_size<-length(lamda2);lamda_size #惩罚因子的个数
    BIC<-matrix(Inf,lamda_size,1);BIC #贝叶斯信息准则的取值，用来选取最优惩罚因子lambda2以及lambda1
    #z<-array(Inf,dim=c(3+length(a1)+length(b1)+3,lamda_size));z#存储每次估计的系数矩阵,1-3三个参数值(lamda1,lamda2,kesai),4a,5b,6time,7panding,8iter
    z <- list();z#设置储存空列表
  }
  for (i in c(1:lamda_size)){
    print(i+10000)
    if (start=='warmstart'){
      ans<-model(a1,b1,lamda1[i],lamda2[i],kesai,n,range_snp,d,x,w,yanzhou
                 ,eps2,eps1,max_iter1,max_iter2,eps3);ans
      a1 <- ans$a1;a1
      b1 <- ans$b1;b1
    }else{
      ans<-model(a1,b1,lamda1[i],lamda2[i],kesai,n,range_snp,d,x,w,yanzhou
                 ,eps2,eps1,max_iter1,max_iter2,eps3);ans
    }
    nozero<-sum(rbind(as.matrix(ans$a1),as.matrix(ans$b1))!=0);nozero
    
    #sigma_L <- sum((residual-mean(residual))^2)/n;sigma_L
    #sse<-sum((yanzhou-x%*%as.matrix(a1)-w%*%as.matrix(b1))^2);sse
    #sse_<-sum((yanzhou-x%*%as.matrix(ans$a1)-w%*%as.matrix(ans$b1))^2);sse_
    residual <- (yanzhou-x%*%as.matrix(ans$a1)-w%*%as.matrix(ans$b1));residual
    sigma_L <- sum((residual)^2)/(n-1);sigma_L
    residual <- (residual)/sqrt(sigma_L);residual
    L <- n*log(2*pi*sqrt(sigma_L))+(1/sigma_L)*sum(residual^2);L
    ##BIC[i,1]<-log(n)*(nozero+2)+n*log(sse/n)#BIC准则
    if (principle=='BIC'){
      BIC[i,1]<-log(n)*(nozero)+L#BIC准则
      BIC[i,1] <- log(n)*(nozero)+sum((residual^2))+2*n*log(sqrt(sigma_L));
    }else if (principle=='AIC'){
      BIC[i,1] <- 2*(nozero)+sum((residual^2))+2*n*log(sqrt(sigma_L));
    }
    #jiehe<-matrix()
    #for (u in c(1:6)){
    #  jiehe<-rbind(jiehe,as.matrix(ans[[u]]))
    #}
    #jiehe<-jiehe[-1,];jiehe
    #z[,i]<-jiehe;
    z <- c(z,list(ans));z
  }
  
  #################################################################################################结果可视化
  BIC <- data.frame(BIC);BIC#转化为数据框
  lamda<-data.frame(lamda1);lamda1#转化为数据框
  canshu <- matrix(rbind(as.matrix(z[[1]]$a1),as.matrix(z[[1]]$b1)));canshu
  for (i in c(2:lamda_size)){
    canshu <- cbind(canshu,rbind(as.matrix(z[[i]]$a1),as.matrix(z[[i]]$b1)));
  }
  canshu <- t(canshu);canshu
  #canshu<-z[4:(3+length(a1)+length(b1)),];canshu#d+number(snp)*d行；length(lambda1)列
  #canshu <- t(canshu);canshu#length(lambda1)行；d+number(snp)*d列
  consequence <- cbind.data.frame(lamda,BIC,as.data.frame(canshu));consequence#数据框:第一列lamda,第二列BIC值，后面是每个snp
  ###############################################################################################画图前执行
  #catagory <- as.matrix(c('','gender','year','near','out'));catagory#协变量的种类
  c <- cbind(c('lamda'),c('BIC'));c#为conseque前两列署名
  for (i in c(1:d)) {
    if(i==1){
      c <- cbind(c,'intercept')
    }else{
      c <- cbind(c,catagory[i])
    }
  }
  for(i in c(1:length(range_snp))){
    for (j in c(1:d)) {
      
      gengxin<-paste('snp',i,sep = '')
      c <- cbind(c,paste(gengxin,catagory[j],sep = ''))
    }
  }
  c <- as.character(c);c##########################################################################
  colnames(consequence)<-c#为数据框consequence署名
  #View(consequence)
  
  weizhi_index<-which(BIC== min(BIC), arr.ind = TRUE);weizhi_index
  lamda_best<-lamda1[weizhi_index[1]];lamda_best#最优惩罚因子大小
  #jieguo <- consequence[weizhi_index[1],];jieguo#将最优模型从数据框中抽取出来
  #View(t(jieguo))#可视化最优模型系数
  #jieguo<-z[,weizhi_index[1]];jieguo
  #jieguo<-jieguo[4:(3+length(a1)+length(b1))];jieguo
  
  #a_best<-t(jieguo[2+1:d+1]);a_best
  #b_best<-t(jieguo[(d+3):(length(range_snp)*d+d+2)]);b_best
  #a_best#############################################最优a
  #b_best#############################################最优b
  
  
  canshu_chushi<-solve(t(A)%*%A)%*%t(A)%*%yanzhou;canshu_chushi
  b1<-as.matrix(canshu_chushi[(d+1):length(canshu_chushi),1]);b1
  #canshu_chushi[(d+1):length(canshu_chushi),1] <- 0;canshu_chushi##############
  b1<-as.matrix(canshu_chushi[(d+1):length(canshu_chushi),1]);b1###############
  a1<-as.matrix(canshu_chushi[1:d,1]);a1
  residual <- (yanzhou-x%*%as.matrix(a1)-w%*%as.matrix(b1));residual
  r_2_ols <- (sum((yanzhou-mean(yanzhou))^2)-sum(residual^2))/sum((yanzhou-mean(yanzhou))^2);r_2_ols################
  canshu_chushi<-solve(t(A)%*%A)%*%t(A)%*%yanzhou;canshu_chushi
  b1<-as.matrix(canshu_chushi[(d+1):length(canshu_chushi),1]);b1
  canshu_chushi[(d+1):length(canshu_chushi),1] <- 0;canshu_chushi##############
  b1<-as.matrix(canshu_chushi[(d+1):length(canshu_chushi),1]);b1###############
  a1<-as.matrix(canshu_chushi[1:d,1]);a1
  a1 <- solve(t(x)%*%x)%*%t(x)%*%(yanzhou-w%*%b1);a1
  residual <- (yanzhou-x%*%as.matrix(a1)-w%*%as.matrix(b1));residual
  r_2_ols_empty <- (sum((yanzhou-mean(yanzhou))^2)-sum(residual^2))/sum((yanzhou-mean(yanzhou))^2);r_2_ols_empty################
  jieguo_best <- t(consequence[weizhi_index[1],]);jieguo_best###################################最优信息汇总
  a1 <- jieguo_best[c(3:(3+d-1)),1];a1
  b1 <- jieguo_best[c((3+d):length(jieguo_best)),1];b1
  residual <- (yanzhou-x%*%as.matrix(a1)-w%*%as.matrix(b1));residual
  r_2_opt <- (sum((yanzhou-mean(yanzhou))^2)-sum(residual^2))/sum((yanzhou-mean(yanzhou))^2);r_2_opt#################
  
  
  
  ###############################################################################################结果呈现
  jieguo_best <- t(consequence[weizhi_index[1],]);jieguo_best###################################最优信息汇总
  jieguo_best <- rbind(jieguo_best,r_2_ols,r_2_ols_empty,r_2_opt,miu,sigma);jieguo_best
  consequence####################################数据框
  
  
  
  ###############################################################################################结果保存
  consequence
  jieguo_best
  BIC
  #site <- 'C:/Users/张丰川/Desktop/联合模型BIC程序/运行结果文件/'
  dd <- paste(site,response_variable,sep = '');dd
  dizhi <- paste(dd,".rda",sep = '');dizhi
  save(consequence,jieguo_best,BIC,file = dizhi)
  
  ###############################################################################################图像展示
  for (i in c(1:(length(range_snp)*d))){
    if(i==1){
      biaohao <- paste('consequence$',c[7+i],sep='');biaohao
      eval(biaohao)
      x_huatu <- as.list.data.frame(consequence$lamda);x_huatu
      y_huatu <- as.list.data.frame(eval(parse(text=biaohao)));y_huatu
      plot(x_huatu,y_huatu,pch=20,cex=0.5,col=i,ylim=c(-1,1))
      lines(x_huatu,y_huatu,lwd=0.6,lty=1.2,col=i,ylim=c(-1,1));
    }else{
      biaohao <- paste('consequence$',c[7+i],sep='');biaohao
      eval(biaohao)
      x_huatu <- as.list.data.frame(consequence$lamda);x_huatu
      y_huatu <- as.list.data.frame(eval(parse(text=biaohao)));y_huatu
      lines(x_huatu,y_huatu,lwd=0.6,lty=1.2,col=i,ylim=c(-1,1));
    }
  }
  #legend("top",legend=c[7+1:(length(range_snp)*d)], ncol=10, cex=0.5, bty="n", col=c(1:(length(range_snp)*d)), lty=1,lwd=1)
  
  res <- list(consequence,jieguo_best,BIC)
  return(res)
}

#rangesnp <- c(9:33)
#range_snp<-c(1:length(rangesnp));range_snp#25个snp,循环器
#d <- 5;



################################################################################################## R方计算];b1
#residual <- (yanzhou-x%*%as.matrix(a1)-w%*%as.matrix(b1));residual
#r_2 <- (sum((yanzhou-mean(yanzhou))^2)-sum(residual^2))/sum((yanzhou-mean(yanzhou))^2);r_2
