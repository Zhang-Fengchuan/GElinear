qujian <- function(left,right,leftpart,yinzi_left,yinzi_right,geshu,a1,b1,kesai,n,range_snp,d,x,w,yanzhou
                   ,eps2,eps1,max_iter1,max_iter2,eps3){
  #left和right可取同一个值，leftpart控制总系数的比例，因子是惩罚因子变动系数（0-1），geshu是设置的惩罚因子个数
  number_coefficent <- max(range_snp)*d+d;number_coefficent
  left_end <- left;
  right_end <- right;
  
  
  ans <- model(a1,b1,left,left,kesai,n,range_snp,d,x,w,yanzhou
               ,eps2,eps1,max_iter1,max_iter2,eps3);ans
  nozero<-sum(rbind(as.matrix(ans$a1),as.matrix(ans$b1))!=0);nozero
  while(nozero<leftpart*number_coefficent){
    left_end <- left_end*yinzi_left;left_end
    ans <- model(a1,b1,left_end,left_end,kesai,n,range_snp,d,x,w,yanzhou
                 ,eps2,eps1,max_iter1,max_iter2,eps3);ans
    
    nozero<-sum(rbind(as.matrix(ans$a1),as.matrix(ans$b1))!=0);nozero
  }
  while(nozero>=leftpart*number_coefficent){
    left_end <- left_end/yinzi_left;left_end
    ans <- model(a1,b1,left_end,left_end,kesai,n,range_snp,d,x,w,yanzhou
                 ,eps2,eps1,max_iter1,max_iter2,eps3);ans
    
    nozero<-sum(rbind(as.matrix(ans$a1),as.matrix(ans$b1))!=0);nozero
  }
  while(nozero>d){
    right_end <- right_end*(1/yinzi_right);right_end
    ans <- model(a1,b1,right_end,right_end,kesai,n,range_snp,d,x,w,yanzhou
                 ,eps2,eps1,max_iter1,max_iter2,eps3);ans
    nozero<-sum(rbind(as.matrix(ans[[4]]),as.matrix(ans[[5]]))!=0);nozero
  }
  while(nozero<=d){
    right_end <- right_end*yinzi_right;right_end
    ans <- model(a1,b1,right_end,right_end,kesai,n,range_snp,d,x,w,yanzhou
                 ,eps2,eps1,max_iter1,max_iter2,eps3);ans
    nozero<-sum(rbind(as.matrix(ans$a1),as.matrix(ans$b1))!=0);nozero
  }
  left_end <- left_end*yinzi_left;left_end
  right_end <- right_end/yinzi_right;right_end
  list_data <- list(left_end,right_end,geshu);list_data
  names(list_data) <- c("left_end","right_end","geshu")
  return(list_data)
}




