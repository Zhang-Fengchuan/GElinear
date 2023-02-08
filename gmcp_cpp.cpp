#include <RcppArmadillo.h>
#include <iostream>
#include <vector>
#include <string>
using namespace Rcpp;
using namespace std;
// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
double Sfun1(double u,double lamda) {
  if(u>lamda){
    return u-lamda;
  }else if(u<(-lamda)){
    return u+lamda;
  }else{
    return 0;
  }
}




// [[Rcpp::export]]
arma::vec Sfun2(arma::vec v,double lamda){
  double a = 1-lamda/norm(v,2);
  arma::vec b = zeros(size(v));
  if (a>0){
    return a*v;
  }else{
    return b;
  }
}




// [[Rcpp::export]]
arma::vec fcpp(arma::mat w, arma::vec b0, arma::uvec index,arma::vec r_, double n, double lamda1, double lamda2, double kesai, double d, arma::mat Q, arma::mat R_inv, double eps1){
  arma::mat A = Q.cols(index-1);
  arma::vec b1 = inv(R_inv.cols(index-1))*b0.rows(index-1); /*对前一步更新后的b取相应位置*/
  arma::vec b1_new;
  w = A; /*w、A都是基变换后的Q对应部分 */
  arma::mat wt = A.t();
  double pandingzhi = 0.1; /*第一步的判定值，保证进入循环进行运算*/
  arma::vec v; v.zeros(d);
  double g;
  while(pandingzhi>eps1){
    if(norm(b1,2)<=(kesai*sqrt(d)*lamda1)){
      if(norm(b1,2)==0){
        g = 1+(1/(norm(b1,2)))*sqrt(d)*lamda1-(1/kesai);
      }else{
        g = 1+(1/(norm(b1,2)))*(sqrt(d)*lamda1-norm(b1,2)/kesai);
      }
    }else{
      g =1;
    }
    arma::vec u = wt*r_/n; /*计算u*/
    for(int i = 0; i<d; i++){
      if(i==0){
        v(i) = u(i);
      }else{
        if (abs(u(i))<=(kesai*lamda2*g)){
          v(i) = Sfun1(u(i),lamda2)/(1-1/(kesai*g));
        }else{
          v(i) = u(i);
        }
      }
    }
    if (norm(v,2)<=(kesai*sqrt(d)*lamda1)){
      b1_new = (kesai/(kesai-1))*Sfun2(v,(sqrt(d)*lamda1));
    }else{
      b1_new = v;
    }
    pandingzhi = norm(abs(b1_new-b1),2);
    b1 = b1_new;
  }
  b1 = R_inv.cols(index-1)*b1;
  return b1;
}





//[[Rcpp::export]] 
/*生成int类型的从a到b的向量*/
arma::uvec create(double a, double b){
  int x = abs(b-a)+1;
  NumericVector ans(x);
  arma::vec consequense;
  for(int i = 0;i <= abs(b-a);i++){
    ans(i) = a+i;
  }
  consequense = ans;
  arma::uvec consequense_ = arma::conv_to<arma::uvec>::from(consequense);
  return consequense_;
}





//[[Rcpp::export]]
/*生成剔除index的索引向量*/
arma::uvec delete_index(arma::uvec index, double d, arma::vec range_snp){
  int number = range_snp.size();
  int total = number*d;
  NumericVector ans(total);
  for(int i = 0;i <=total-1;i++){
    ans(i) = 1+i;
  }
  int number_ans = index.size();
  NumericVector ans_(number_ans);
  for(int i = 0;i <=number_ans-1;i++){
    ans_(i) = index(i);
  }
  ans.erase(ans_(0)-1,ans_(d-1));
  arma::vec consequence = ans;
  arma::uvec consequence_ = arma::conv_to<arma::uvec>::from(consequence);
  return consequence_;
}




// [[Rcpp::export]]
arma::vec Fcpp(arma::vec range_snp, arma::vec b0_hou, arma::vec b0_qian, arma::mat wquan, arma::mat yanzhouquan,
               arma::mat xquan, arma::vec a0_qian, double n, double lamda1, double lamda2, double kesai, double d,
               arma::mat Q, arma::mat R_inv, double eps2, double eps1, double max_iter1){
  int length = range_snp.size();
  double b_panding = 100;/*保证第一步进入while循环*/
int iter = 0;
if(length==1){
  b0_hou = b0_qian;
  arma::uvec index = create(1,d);
  arma::mat w = wquan.cols(index-1);
  arma::vec r_ = yanzhouquan-xquan*a0_qian;
  arma::vec ans = fcpp(w,b0_qian,index,r_,n,lamda1,lamda2,kesai,d,Q,R_inv,eps1);
  b0_qian(index-1) = ans;   
}else{
  arma::uvec index; /*= create(((1-1)*d+1),(1*d));   提前在循环外声明变量*/
arma::mat w;   /*= wquan.cols(index-1);     提前在循环外声明变量*/
arma::mat w_;  /*= wquan.cols(delete_index(index,d,range_snp)-1);      提前在循环外声明变量*/
arma::vec b1_; /*= b0_qian(delete_index(index,d,range_snp)-1);     提前在循环外声明变量*/
arma::vec r_;  /*= yanzhouquan-xquan*a0_qian-w_*b1_;    提前在循环外声明变量*/
arma::vec ans; /*= fcpp(w,b0_qian,index,r_,n,lamda1,lamda2,kesai,d,Q,R_inv,eps1);     提前在循环外声明变量*/
while(b_panding > eps2 && iter < max_iter1){/*while循环误差大于eps2且迭代步数小于最大迭代步数时继续计算*/
b0_hou = b0_qian;/*存储前一步的b*/
for(int i = 1;i<=length;i++){
  index = create(((i-1)*d+1),(i*d));
  w = wquan.cols(index-1);
  w_ = wquan.cols(delete_index(index,d,range_snp)-1);
  b1_ = b0_qian(delete_index(index,d,range_snp)-1);
  r_ = yanzhouquan-xquan*a0_qian-w_*b1_;
  ans = fcpp(w,b0_qian,index,r_,n,lamda1,lamda2,kesai,d,Q,R_inv,eps1);
  b0_qian(index-1) = ans;
}
b_panding = norm(b0_qian-b0_hou,2);
/*b_panding = 0.5*eps2;*/
iter = iter+1;
}
}
/*arma::uvec index = create(1,d);
 arma::mat w = wquan.cols(index-1);
 arma::vec r_ = yanzhouquan-xquan*a0_qian;
 arma::vec ans = fcpp(w,b0_qian,index,r_,n,lamda1,lamda2,kesai,d,Q,R_inv,eps1);
 arma::mat w_ = wquan.cols(delete_index(index,d,range_snp)-1);
 arma::vec b1_ = b0_qian(delete_index(index,d,range_snp)-1);
 arma::vec r_2 = yanzhouquan-xquan*a0_qian-w_*b1_;
 arma::vec ans_ = fcpp(w,b0_qian,index,r_2,n,lamda1,lamda2,kesai,d,Q,R_inv,eps1);
 b0_qian(index-1) = ans;
 double b_panding = norm(b0_qian-b0_hou,2);*/
/*List res = List::create(Named("b0_qian")=b0_qian, _["b_panding"] = b_panding, _["iter_1"] = iter);*/
return b0_qian;
}





// [[Rcpp::export]]
List Fcpp2(arma::vec range_snp, arma::vec b0_hou, arma::vec b0_qian, arma::mat wquan, arma::mat yanzhouquan,
           arma::mat xquan, arma::vec a0_qian, double n, double lamda1, double lamda2, double kesai, double d,
           arma::mat Q, arma::mat R_inv, double eps2, double eps1, double max_iter1){
  int length = range_snp.size();
  double b_panding = 100;/*保证第一步进入while循环*/
int iter = 0;
if(length==1){
  b0_hou = b0_qian;
  arma::uvec index = create(1,d);
  arma::mat w = wquan.cols(index-1);
  arma::vec r_ = yanzhouquan-xquan*a0_qian;
  arma::vec ans = fcpp(w,b0_qian,index,r_,n,lamda1,lamda2,kesai,d,Q,R_inv,eps1);
  b0_qian(index-1) = ans;   
}else{
  arma::uvec index = create(((1-1)*d+1),(1*d));/*提前在循环外声明变量*/
arma::mat w = wquan.cols(index-1);/*提前在循环外声明变量*/
arma::mat w_ = wquan.cols(delete_index(index,d,range_snp)-1);/*提前在循环外声明变量*/
arma::vec b1_ = b0_qian(delete_index(index,d,range_snp)-1);/*提前在循环外声明变量*/
arma::vec r_ = yanzhouquan-xquan*a0_qian-w_*b1_;/*提前在循环外声明变量*/
arma::vec ans = fcpp(w,b0_qian,index,r_,n,lamda1,lamda2,kesai,d,Q,R_inv,eps1);/*提前在循环外声明变量*/
while(b_panding > eps2 && iter < max_iter1){/*while循环误差大于eps2且迭代步数小于最大迭代步数时继续计算*/
b0_hou = b0_qian;/*存储前一步的b*/
for(int i = 1;i<=length;i++){
  index = create(((i-1)*d+1),(i*d));
  w = wquan.cols(index-1);
  w_ = wquan.cols(delete_index(index,d,range_snp)-1);
  b1_ = b0_qian(delete_index(index,d,range_snp)-1);
  r_ = yanzhouquan-xquan*a0_qian-w_*b1_;
  ans = fcpp(w,b0_qian,index,r_,n,lamda1,lamda2,kesai,d,Q,R_inv,eps1);
  b0_qian(index-1) = ans;
}
b_panding = norm(b0_qian-b0_hou,2);
iter = iter+1;
}
}
List res = List::create(Named("b0_qian")=b0_qian, _["b_panding"] = b_panding, _["iter_1"] = iter);
return res;
}






//[[Rcpp::export]]
List FFcpp(arma::vec range_snp, arma::vec b0_hou, arma::vec b0_qian, arma::mat wquan, arma::mat yanzhouquan,
           arma::mat xquan, arma::vec a0_qian, double n, double lamda1, double lamda2, double kesai, double d,
           arma::mat Q, arma::mat R_inv, double eps2, double eps1, double max_iter1, double max_iter2, double eps3){
  double z_panding = 100;
  int iter = 0;
  arma::mat xquant = xquan.t();
  arma::vec b1 = b0_qian;
  arma::vec a1 = a0_qian;
  while(z_panding > eps3 && iter < max_iter2){
    b1 = b0_qian;
    b0_qian = Fcpp(range_snp,b0_qian,b0_qian,wquan,yanzhouquan,xquan,a0_qian,n,lamda1,lamda2,kesai,d,Q,R_inv,eps2,eps1,max_iter1);
    a1 = a0_qian;
    
    a0_qian = inv(xquant*xquan)*xquant*(yanzhouquan-wquan*b0_qian);
    z_panding = sqrt(norm(a0_qian-a1,2)*norm(a0_qian-a1,2)+norm(b0_qian-b1,2)*norm(b0_qian-b1,2));
    b1 = b0_qian;
    iter = iter+1;
  }
  List res = List::create(Named("a1")=a1, _["b1"] = b1, _["panding"] = z_panding, _["iter"] = iter);
  return res;
}











// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//




