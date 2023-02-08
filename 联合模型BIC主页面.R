#----------------------------------------------------------------------------------------------------------------------------------------------
#---------------------------数据准备-----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------
name_x <- c('IID','near','out');name_x###############输入规范
name_y <- c('WEINO','AGE','AGE_1','AGE_2','GENDER','AL','AL_1','AL_2')###############输入规范
number_repeat <- 3###########堆叠次数
response_variable <-'AL' 
name_year <- 'AGE'
rangesnp <- c(9:33)###snp表中范围 ##############################################################边际联合的w开关(9:33)
c('全一','gender','year','near','out')
range_zscore <- c(3,4,5)
#####################################################################################数据加载
position_snp <- 'C:/Users/张丰川/Desktop/联合模型BIC程序原始/data.xlsx';position_snp
position_phenotype <- 'C:/Users/张丰川/Desktop/联合模型BIC程序原始/middle_school.xlsx';position_phenotype
grouping='yes'######yes or no是否将协变量的近距离和远距离工作时间按照中位数二分类
data <- yuchuli(position_snp,position_phenotype,name_x,name_y,number_repeat,response_variable,name_year,rangesnp,range_zscore,grouping)
x <- as.matrix(data$x);x   
w <- as.matrix(data$w);w    
yanzhou <- as.matrix(data$y);yanzhou    
na_matrix <-as.matrix(data$na_matrix);na_matrix
miu <- mean(yanzhou);miu
sigma <- sd(yanzhou);sigma
###############################################################################################因变量y标准化开关
yanzhou <- (yanzhou-mean(yanzhou))/sd(yanzhou);yanzhou
y <- yanzhou;y
kesai <- 3;kesai
leftpart <- 1;leftpart
left <- 0.001;left
right <- 0.25;right
yinzi_left <- 0.5;yinzi_left
yinzi_right <- 0.8;yinzi_right
geshu <- 200;geshu#惩罚因子设置100个点 ##############################################################  ###################################
eps1 <- 1e-8;eps1
eps2 <- 1e-3;eps2
eps3 <- 1e-5;eps3
max_iter1 <- 200;max_iter1
max_iter2 <- 50;max_iter2
start <- 'lassostart'
response_variable <- 'AL'
site <- 'C:/Users/张丰川/Desktop/结果/'
lambda <- cbind(seq(0,0.4,length.out=100),seq(0,0.4,length.out=100));lambda
catagory <- as.matrix(c('','gender','year','near','out'));catagory

#----------------------------------------------------------------------------------------------------------------------------------------------
#----------------------执行命令行--------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------
s <- Sys.time()#计时开始
output <- GElinear(x,w,y,lambda,200,3,start,'BIC',catagory,response_variable,site,leftpart,
                   left,right,yinzi_left,yinzi_right,eps1,eps2,eps3,max_iter1,max_iter2);output
e <- Sys.time()
print(e-s)#计时结束







