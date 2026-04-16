## version 0.2
# source("https://bioconductor.org/biocLite.R") 通过bioconductor官网下载相应的R包
# biocLite("edgeR")
args = commandArgs(TRUE)
a=length(args)
if(a<4){
	cat("Usage: Rscript edgeR.r <count.xls> <col_num1,col_num2> <outpfx> <dispersion>\n")
	q()
}
## input prepare
infile = args[1]
cols=as.numeric((unlist(strsplit(args[2],','))))
group_colnames=(unlist(strsplit(args[3],',')))
#cols
allcols=sum(cols)
outpfx = args[4]
dispersion = args[5]  # 没有重复时使用，人为设定生物学变异系数，edgeR建议：人类实验0.4，遗传背景相似的模式生物0.1，技术重复0.01 
pq = as.numeric(args[6])  ## 1=> p , 2=>q
cutp = as.numeric(args[7])
cutfc = as.numeric(args[8])
cutfc = log2(cutfc)
norm = as.character(args[9])
fpkm = as.character(args[10])
pair_test = args[11] ## [no] / yes 

if(allcols == 2) pair_test = "no" ## 1-vs-1 always no
#cutfc = 1 ## log2(fc) filter cutoff
#cutp = 0.05 ## pvalue filter cutoff

pq1=c("PValue","FDR")
pq=pq1[pq]
#pq
##1. 导入文件，并计算校正因子
library(edgeR,quietly=TRUE)   
data = read.table(args[1], header=T, row.names=1, sep='\t',check.names = F,quote = "",comment.char = "",na.strings="")
x = data[, 1:ncol(data)] #读入所有列；
all_list = DGEList(counts=x, group=colnames(x)) # 生成每个样本差异表达矩阵
all_list <- calcNormFactors(all_list,method="TMM")    #从所有样本预估校正因子
#all_list  # 查看下y数组的内容，注意校正因子的数值； 
#plotMDS(all_list)  # 使用mutlidimensional scaling的方法进行聚类，MDS在算法上与PCA分析等效。

##2. 准备数据，并进行差异分析
C_number=cols[1];   # 对照组样本数
T_number=cols[2];   # 处理组样本数
col_ordering = 1:allcols # 设定 用于差异分析的样本的编号，这里选择A组和C组
#head(data)
DiffMatrix = data[,col_ordering] # 从总表中筛选用于差异的样本数据
## filter row and num
filter = 0.001
f1=(rowSums(DiffMatrix)>filter)
f2=(rowSums(DiffMatrix)<=filter)
f1n=nrow(DiffMatrix[rowSums(DiffMatrix)>filter,])
f2n=nrow(DiffMatrix[rowSums(DiffMatrix)<=filter,])
#f2n
f2name = rownames(DiffMatrix[rowSums(DiffMatrix)<=filter,])

DiffMatrix = DiffMatrix[rowSums(DiffMatrix)>filter,] # 只挑选readcount均值大于1的基因，即去除不表达的基因;
conditions = factor(c(rep("control", C_number), rep("treat", T_number))) # 设定分组信息名称，注意两个变量的数值是提前设定好的
Diff_list = DGEList(counts=DiffMatrix, group=conditions) # 生成两个比较组的差异表达矩阵
#Diff_list$samples[3]=rep(1,ncol(DiffMatrix))  # 校正因子都使用默认的1
# # 以下2行的循环命令，则是将校正因子替换为edgeR计算出的校正因子； 
Diff_list = calcNormFactors(Diff_list,method="TMM")
if(norm == "no"){
	Diff_list$samples[3]=rep(1,ncol(DiffMatrix))
}
Diff_list$samples[3]
# name= rownames(Diff_list$samples)
# for (i in 1:length(name)){Diff_list$samples[name[i],3]= all_list$samples[name[i],3]}

if(pair_test == "yes"){
	Patient <- factor(c(1:C_number,1:T_number))                    # 实验因素1：样本
	Tissue <- factor(rep(c("N","T"),each=C_number))             # 实验因素2：癌和癌旁
	design <- model.matrix(~Patient+Tissue)
	rownames(design) <- colnames(Diff_list)
	# Estimating the dispersion
	Diff_list <- estimateDisp(Diff_list, design, robust=TRUE,trend.method="locfit")
	fit <- glmFit(Diff_list, design)    # 构建广义线性回归模型
	etest <- glmLRT(fit,coef=ncol(design))
	print(design)
}else if(C_number > 1 && T_number > 1 && dispersion == "auto"){

#if (dispersion == 'no' & (C_number > 1 || T_number > 1)){             #或的关系，只要求一组有重复；
#if(C_number > 1 && T_number > 1){
    Diff_list = estimateCommonDisp(Diff_list)  # 计算组的离散度
    Diff_list = estimateTagwiseDisp(Diff_list) # 计算每个基因的离散度
#    plotBCV(Diff_list)
    etest = exactTest(Diff_list)                  # 精确检验
	cat("group\n")
}else{
	if(dispersion == "auto") dispersion = 0.01
    etest = exactTest(Diff_list, pair = c("control", "treat"), dispersion=as.numeric(dispersion)^2)  # 选择组名称，认为设定bcv
#	cat("dispersion\n")
}
tTags = topTags(etest,n=NULL) # 排序矩阵,并计算FDR


#head(data)

## calculate cpm
data2=data

if(fpkm == "none"){
	for(i in 1:allcols)
	{
		data2[,i] = round(data[,i]*1000000/sum(as.numeric(data[,i])),2)
		colnames(data2)[i] = paste(colnames(data)[i],"_CPM",sep="")
	}
}else{
	data2 = read.table(fpkm, header=T, row.names=1, sep='\t',check.names = F,quote = "",comment.char = "")
}
#head(data2)

## mean 
data3=data2[,1:2]
data3[,1] = apply(data.frame(data2[,1:C_number]), 1, mean )
data3[,2] = apply(data.frame(data2[,(C_number+1):allcols]), 1, mean)
data3[which(data3[,1]==0),1]=0.001
data3[which(data3[,2]==0),2]=0.001
log2fc=log2(data3[,2]/data3[,1])
colnames(data3)=paste(group_colnames,"_mean",sep="")
#head(data3)

## ouput file, id,count,log2(fc),pvalue,qvalue
data4=data.frame(rownames(data),data,data2,data3,check.names=FALSE)
colnames(data4)[1] = "id"
data4$'log2(fc)' = log2fc
#head(tTags@.Data[[1]])
# pair.test => logFC    logCPM       LR       PValue          FDR 
# normal => logFC    logCPM PValue          FDR
pq_col_list = c(3:4) + ifelse(pair_test == "yes" , 1, 0)
datat = tTags@.Data[[1]][,pq_col_list]
datat$id = rownames(datat)
if(f2n != 0){
	n = nrow(datat)
	datat[(n+1):(n+f2n),1:3]=rep(1,f2n)
	datat[(n+1):(n+f2n),2]=rep(1,f2n)
	datat[(n+1):(n+f2n),3]=f2name
}

data5 = merge(data4,datat,by="id")
#head(data5)

#write.table(data2, file=paste(outpfx,".edgeR.cpm.xls",sep=""), sep='\t', quote=F, row.names=T)  # 鍐欏嚭鏁版嵁
#write.table(data3, file=paste(outpfx,".edgeR.mean.xls",sep=""), sep='\t', quote=F, row.names=T)  # 鍐欏嚭鏁版嵁
#write.table(tTags, file=paste(outpfx,".edgeR.DE_results.xls",sep=""), sep='\t', quote=F, row.names=T)  # 鍐欏嚭鏁版嵁
write.table(data5, file=paste(outpfx,".all.xls",sep=""), sep='\t', quote=F, row.names=F)  # 写出数据
data6 = data5[which(abs(data5[,'log2(fc)'])>cutfc & data5[,pq]<cutp),]
#head(data6)
write.table(data6, file=paste(outpfx,".filter.xls",sep=""), sep='\t', quote=F, row.names=F)

1
## 3. 绘图 
#summary(de <- decideTestsDGE(etest,p.value=0.05,lfc=1)) # 设定差异基因 FDR 和 log2FC的阈值，lfc的设置只对R 3.2.3的版本有效；
#detags <- rownames(Diff_list)[as.logical(de)]  # 设定一个是否显著差异的逻辑判断向量
#plotSmear(etest, de.tags=detags) # 绘制火山图
#abline(h=c(-1, 1), col="blue")   # 辅助线
