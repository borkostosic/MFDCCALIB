rm(list = ls(all.names = TRUE))
dirname <- "D:\\laptop_07_02_2010\\D\\papers\\MDFA\\MFDCCALIB\\GITHUB\\R"
setwd(dirname)
fname <- "..\\DATA\\sug_eta_returns.dat"

seq <- read.table(fname, header=TRUE)				# load data

MAX_BOX<-200		# maximum 200 points on logarithmic scale...
MAXQ<-201		# max q resolution -10,...,10 dq=0.1

data1<-as.numeric(unlist(seq[1]))
data2<-as.numeric(unlist(seq[2]))
#p<-0.3

total <- length(data1)	# data length
qmin <- -10.000001	# small offset to avoid 0
qmax <- 10.0
dq <- 0.1
qseq <- seq(qmin, qmax, by = dq)	# sequence of q's
nq<-length(qseq)			# length of q sequence
H <- numeric(length(qseq))		# generalized Hurst exponet
tau <- numeric(length(qseq))		# Renyi exponet
f <- numeric(length(qseq))		# f(alpha)
alpha <- numeric(length(qseq))		#alpha
dmse<-as.numeric(matrix(0, nrow=MAXQ, ncol=MAX_BOX))	# fluctuations logarithm (base 10)
nfit<-1		# polynomial degree 1 (linear), 2 or 3
sw<-0		# sw=1 for sliding segments, takes longer, use with care... 
perc<-0

# version of MFDCCA algorithm:
# 1 - Original MF-DXA W.-X. Zhou, Phys. Rev. E 77, 066211 (2008).
# 2 - ABS
# 3 - MFCCA Phys. Rev. E 89, 023305 (2014)
# 4 - Plus sum
# 5 - Minus sum
# 6 - Plus Box
# 7 - Minus box
# 8 - PP
# 9 - PM
# 10 - MP
# 11 - MM
ver0<- 1	# used for MFDFA of indivudual series
ver <- 4	# current choice of MFDCCA flavor

rs<-integer(MAX_BOX)			# segment size list
# prepare multiplicative scale (equidistant on log-log plot)
   nrs<-0
   minseg<-4		# minimum segment size
   maxseg<-total/4	# maximum segment size
   boxratio<-2^(1/8)	# segment ratio (multiplicative scale)
#...OR prepare your own segment scale
   #rs<-10:100 	
   #nrs<-length(rs)

dllname <- paste(dirname,"\\mfdcca64.dll",sep="")	# dll name
dyn.load(dllname)					# load the dynamic library

ans1 <- .C("mfdcca", as.integer(ver0),
	as.numeric(data1), as.numeric(data1), as.integer(total), as.numeric(H), as.numeric(tau), as.numeric(f), as.numeric(alpha),
	as.numeric(dmse), as.integer(rs), as.integer(nrs),as.numeric(perc), as.numeric(qmin), as.numeric(qmax), as.numeric(dq),
	as.integer(minseg), as.integer(maxseg), as.numeric(boxratio), as.integer(nfit), as.integer(sw)
	)
H1<-ans1[[5]]
tau1<-ans1[[6]]
f1<-ans1[[7]]
alpha1<-ans1[[8]]

ans2 <- .C("mfdcca", as.integer(ver0),
	as.numeric(data2), as.numeric(data2), as.integer(total), as.numeric(H), as.numeric(tau), as.numeric(f), as.numeric(alpha),
	as.numeric(dmse), as.integer(rs), as.integer(nrs),as.numeric(perc), as.numeric(qmin), as.numeric(qmax), as.numeric(dq),
	as.integer(minseg), as.integer(maxseg), as.numeric(boxratio), as.integer(nfit), as.integer(sw)
	)
H2<-ans2[[5]]
tau2<-ans2[[6]]
f2<-ans2[[7]]
alpha2<-ans2[[8]]

ans <- .C("mfdcca", as.integer(ver),
	as.numeric(data1), as.numeric(data2), as.integer(total), as.numeric(H), as.numeric(tau), as.numeric(f), as.numeric(alpha),
	as.numeric(dmse), as.integer(rs), as.integer(nrs),as.numeric(perc), as.numeric(qmin), as.numeric(qmax), as.numeric(dq),
	as.integer(minseg), as.integer(maxseg), as.numeric(boxratio), as.integer(nfit), as.integer(sw)
	)

dyn.unload(dllname)					# unload the dll

H<-ans[[5]]
tau<-ans[[6]]
f<-ans[[7]]
alpha<-ans[[8]]
dmse<-ans[[9]]
rs<-ans[[10]]
nr<-ans[[11]] 
perc<-ans[[12]] 

# 4 figures arranged in 2 rows and 2 columns
par(mfrow=c(2,2))
plot(qseq,tau, main="Rényi exponent",cex=1, cex.lab=1.6, cex.axis=1.6,xlab = "q", ylab=expression(tau))
lines(qseq,(tau1+tau2)/2,type="l", col=3, lwd=2)

plot(log10(rs[1:nr]),log10(dmse[1:nr]), main="Fluctuations",
	cex=1, cex.lab=1.6, cex.axis=1.6, xlab = "log(s)", ylab = "logF(s)",
	ylim=c(min(log10(dmse[dmse>0])),max(log10(dmse[is.finite(dmse)]))))
for (i in 2:nq){ 
  points(log10(rs[1:nr]), log10(dmse[(i*MAX_BOX+1):(i*MAX_BOX+nr)]), col=i,pch=1) 
}

plot(alpha,f, main="Multifractal spectrum",cex=1, cex.lab=1.6, cex.axis=1.6, col=4, xlab=expression(alpha))
lines((alpha1+alpha2)/2,(f1+f2)/2,type="l", col=3, lwd=2)

#plot(qseq,H, main="Generalized Hurst exponent",cex=1, cex.lab=1.6, cex.axis=1.6, xlab = "q", ylim=c(0.5,1.5))
plot(qseq,H, main="Generalized Hurst exponent",cex=1, cex.lab=1.6, cex.axis=1.6, xlab = "q")
lines(qseq,(H1+H2)/2,type="l", col=3, lwd=2)

print(sprintf("used %.2f%% data pairs",perc))
