model<-function(a1,b1,lamda1,lamda2,kesai,n,range_snp,d,x,w,yanzhou
                ,eps2,eps1,max_iter1,max_iter2,eps3){
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
  
  #cholesky分解
  Q <- matrix(,n,ncol(w));Q #预先分配内存，定义空矩阵
  R_inv <- matrix(,d,ncol(w));R_inv #预先分配内存，定义空矩阵
  for (i in range_snp) {
    XX <- (t(w[,((i-1)*d+1):(i*d)])%*%w[,((i-1)*d+1):(i*d)])/n;XX
    RR <- chol(XX);RR
    c <- w[,((i-1)*d+1):(i*d)]%*%solve(RR);c
    Q[,((i-1)*d+1):(i*d)] <- c;
    R_inv[,((i-1)*d+1):(i*d)] <- solve(RR)
  }
  
  #计算一个惩罚因子对应的最优解
  s <- Sys.time()
  ans <- FFcpp(range_snp,b1,b1,w,yanzhou,x,a1,n,lamda1,lamda2
               ,kesai,d,Q,R_inv,eps2,eps1,max_iter1,max_iter2,eps3);ans
  e <- Sys.time()
  time <- e-s;time
  
  a1 <- ans$a1;a1#提取协变量系数
  b1 <- ans$b1;b1#提取snp及交互项系数
  panding <- ans$panding;panding#最外层a与b的前后两步系数差二范数
  iter <- ans$iter;iter#最外层a与b的迭代步数
  list_data <- list(lamda1,lamda2,kesai,a1,b1,time,panding,iter);list_data
  names(list_data) <- c("lamda1","lamda2","kesai","a1","b1","time","panding","iter")
  return(list_data)
  
}
