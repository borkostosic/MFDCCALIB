rm(list = ls(all.names = TRUE))
dirname <- "D:\\laptop_07_02_2010\\D\\papers\\MDFA\\MFDCCALIB\\GITHUB\\R"	# adjust for your own directory path
setwd(dirname)
fname <- "..\\DATA\\BIOMIAL_MULTIFRACTAL_1048576_0.3_0.4.dat"

seq <- read.table(fname, header=TRUE)				# load data

MAX_BOX<-200		# maximum 200 points on logarithmic scale...
MAXQ<-201		# max q resolution -10,...,10 dq=0.1

data1<-as.numeric(unlist(seq[1]))
data2<-as.numeric(unlist(seq[1]))
p<-0.3

total <- length(data1)	# data length
qmin <- -10.000001	# small offset to avoid 0
qmax <- 10.0
dq <- 1.0
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
flavor <- 1

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
start.time <- Sys.time()	# call mfdfa_calc...
ans <- .C("mfdcca", as.integer(flavor),
	as.numeric(data1), as.numeric(data2), 
	as.integer(total), 
	as.numeric(H), as.numeric(tau), 	
	as.numeric(f), as.numeric(alpha),
	as.numeric(dmse), as.integer(rs), as.integer(nrs),as.numeric(perc),  
	as.numeric(qmin), as.numeric(qmax), as.numeric(dq),
	as.integer(minseg), as.integer(maxseg), as.numeric(boxratio),
	as.integer(nfit), as.integer(sw)
	)
end.time <- Sys.time()
time.taken <- end.time - start.time


dyn.unload(dllname)					# load the dll

H<-ans[[5]]
tau<-ans[[6]]
f<-ans[[7]]
alpha<-ans[[8]]
dmse<-ans[[9]]
rs<-ans[[10]]
nr<-ans[[11]] 
perc<-ans[[12]] 

#calculate theoretical curves
tot<-100
H_t <- numeric(tot)
tau_t <- numeric(tot)
f_t <- numeric(tot)
alpha_t <- numeric(tot)
q_t <- numeric(tot)
for (i in 1:tot){
q<-qmin+i*(qmax-qmin)/tot
q_t[i]<-q
H_t[i]<-0.1e1 / q - log((p)^q + (0.1e1 - p)^q) / q / log(0.2e1);
tau_t[i]<-q * (0.1e1 / q - log((p)^q + (0.1e1 - p)^q) / q / log(0.2e1)) - 0.1e1;
f_t[i]<- -q * ((p)^q * log(p) + 
	(0.1e1 - p)^q * log(0.1e1 - p)) / ((p)^q + (0.1e1 - p)^q) / log(0.2e1) - 
	q * (0.1e1 / q - log((p)^q + (0.1e1 - p)^q) / q / log(0.2e1)) + 0.1e1
alpha_t[i]<- -((p)^q * log(p) + (0.1e1 - p)^q * log(0.1e1 - p)) / ((p)^q + 
	(0.1e1 - p)^q) / log(0.2e1)
}

# 4 figures arranged in 2 rows and 2 columns
par(mfrow=c(2,2))
plot(qseq,tau, main="Rényi exponent",cex=1, cex.lab=1.6, cex.axis=1.6,xlab = "q", ylab=expression(tau))
lines(q_t,tau_t,type="l", col=2, lwd=2)

plot(log10(rs[1:nr]),log10(dmse[1:nr]), main="Fluctuations",
	cex=1, cex.lab=1.6, cex.axis=1.6,xlab = "log(s)", ylab = "logF(s)",
	ylim=c(min(log10(dmse[dmse>0])),max(log10(dmse[is.finite(dmse)]))))
for (i in 2:nq){ 
  points(log10(rs[1:nr]), log10(dmse[(i*MAX_BOX+1):(i*MAX_BOX+nr)]), col=i,pch=1) 
}
plot(alpha,f, main="Multifractal spectrum",cex=1, cex.lab=1.6, cex.axis=1.6, col=4, xlab=expression(alpha))
lines(alpha_t,f_t,type="l", col=2, lwd=2)

plot(qseq,H, main="Generalized Hurst exponent",cex=1, cex.lab=1.6, cex.axis=1.6, xlab = "q")
lines(q_t,H_t,type="l", col=2, lwd=2)

print(time.taken)
print(sprintf("used %.2f%% data pairs",perc))
