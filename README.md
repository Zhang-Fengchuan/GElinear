# GElinear
描述：本项目为基于稀疏惩罚的线性回归联合分析模型，适用于基因环境交互风险因素的识别，本项目还依赖Rcpp、rJava等R包，使用时请提前安装并配置环境。
作用：可以计算给定调优参数下的回归系数，并在给定的准则(如BIC或AIC)下选择最优回归系数。
gmcp_cpp.cpp是c++文件，用于计算算法最内层的迭代方程。
model.R用于计算在一个调优参数下的最优回归系数。
qujian.R用于计算调优参数取值范围的函数。
yuchuli.R用于处理excel格式的基因环境数据。
GE_linear是整体的联合分析模型，给定参数条件下可选择最优回归系数。
联合模型BIC主页面为调用GE_linear进行具体分析的页面。
data_snp.xlsx和middle_school.xlsx是基因环境数据。
