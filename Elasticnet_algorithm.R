#function for creating random cities
random_cities<-function(cities){
  x_coor <- runif(cities, 0.2, 1.1)
  y_coor <- runif(cities, 0.2, 1.1)
  coors <- cbind(x_coor,y_coor) 
  return(coors)
}

#function for creating a plot of the randomised tours
plotdata<-function(coors,neurons){
  plot(coors[,1],coors[,2], ylim=c(0,1.2), xlim=c(0,1.2),
        xlab="", ylab="")
  points(neurons[,1],neurons[,2], pch=".", col="red", cex=3)
  neurons<-rbind(neurons,neurons[1,])
  for (i in 1:(nrow(neurons)-1)){
    segments(neurons[i,1],neurons[i,2],neurons[i+1,1],
             neurons[i+1,2], col="red")}
  grid(12,12,lwd=1.1)
}

#function for creating neurons of the elastic net
neuron<-function(c,coors,r=0.1){
  c=c+10
  if (c == 2){
    t<-seq(0, (2*pi), by = pi)
  } else {
    t<-seq(0, (1.9*pi), length = c)
  }
  x_coor<-c()
  y_coor<-c()
  
  #position elastic net on the centroid of data
  for (i in 1:c){
    x_coor[i]<-r*cos(t[i])+mean(coors[,1])#mean_coorsx
    y_coor[i]<-r*sin(t[i])+mean(coors[,2])#mean_coorsy
  }
  neuron <- cbind(x_coor,y_coor) 
  return(neuron)
}

#function that calculates distance between points
distance = function(coor,neuron) {
  line_len<-round(sqrt((coor[1]-neuron[1])^2+
                         (coor[2]-neuron[2])^2),8)
  return(as.numeric(line_len))
}

#create a function to calculate energy
energy<-function(k,y){
  a=alpha; b=beta;c=cities+10
  d=cities; x=coors
  #y=neurons; k=0.2
  all_exp<-list()
  for (i in 1:c){
    exps<-c()
    for (n in 1:d){
      exps[n]<-exp(-1*(distance(x[n,],y[i,])^2)/(2*k^2))
    }
    all_exp[[i]]<-exps
  }
  A1<-c()
  for(i in 1:length(all_exp)){
    A1[i]<-log(sum(all_exp[[i]]))
  }
  A<--1*a*k*sum(A1)
  y_diff_term<-c()
  for (i in 1:c){
    if(i==c){
      term1=y[1,]
      term2=y[i,]
    }else{
      term1=y[(i+1),]
      term2=y[i,]
    }
    y_diff_term[i]<-distance(term1,term2)^2
  }
  B<-b*sum(y_diff_term)
  E<-A+B
  return(E)
}

#create weight update function
delta_y<-function(k,y){
  c=cities+10; a=alpha;x=coors;b=beta;
  #y=neurons; k=0.2;
  d=cities
  all_exp<-list()
  xi_yj<-list()
  for (i in 1:d){
    exps<-c()
    diff_matrix<-matrix(0,c,2)
    for (n in 1:c){
      exps[n]<-exp(-1*(distance(x[i,],y[n,])^2)/(2*k^2))
      diff_matrix[n,]<-(x[i,]-y[n,])
    }
    all_exp[[i]]<-exps
    xi_yj[[i]]<-diff_matrix
  }
  
  wij<-list()
  for(i in 1:length(all_exp)){
    wij[[i]]<-t(all_exp[[i]])/sum(all_exp[[i]])
  }
  
  xi_yj_as_matrix<-do.call(cbind, xi_yj)
  wij_as_matrix<-matrix(0,c,c)
  allexp_as_matrix<-matrix(0,c,c)
  for (i in 1:length(wij)){
    wij_as_matrix[i,]<-as.vector(t(wij[[i]]))
    allexp_as_matrix[i,]<-as.vector(t(all_exp[[i]]))
  }
  
  start <- seq(1, by = 2, length = ncol(xi_yj_as_matrix) / 2)
  prod_matrix_aslist<-list()
  count=0
  for (i in start){
    count=count+1
    prod_matrix_aslist[[count]]<-wij_as_matrix[count,]*xi_yj_as_matrix[,i:(i+1)]
  }
  
  summation<-Reduce('+', prod_matrix_aslist) #summation over neurons
  
  y_diff_term<-matrix(0,c,2)
  for (i in 1:c){
    if (i == 1){
      term1=y[(i+1),]
      term2=2*y[i,]
      term3=y[c,]
    } else if(i==c){
      term1=y[1,]
      term2=2*y[i,]
      term3=y[(i-1),]
    }else{
      term1=y[(i+1),]
      term2=2*y[i,]
      term3=y[(i-1),]
    }
    y_diff_term[i,]<-term1-term2+term3
  }
  
  dy=(a*summation)+(b*k*y_diff_term)
  return(dy)
}

cities <-20 #number of cities
coors<-random_cities(cities)
neurons<-neuron(cities,coors)
dev.off();plotdata(coors,neurons)
K=0.25;alpha=0.3;beta=2;iter=10000; n=25
M=0.25*cities
E<-c(); temp<-c()

for (q in 1:iter){
  temp[q]<-K
  E[q]=energy(K,neurons)
  delta_neurons<-delta_y(K,neurons)
  neurons<-neurons+delta_neurons
  if (q%%n==0){
    K=0.99*K
  }
  if (q%%1000==0){
    dev.off()
    plotdata(coors,neurons)
    Sys.sleep(1)
  }
}

plot(seq(1,iter,1),E, type="l", col="indianred", 
     main= "Energy Minimisation",
     xlab= "Iteration",
     ylab= "Energy",
     cex.main=1.7, cex.axis=1,cex.lab=1.4)
grid(13,13,lwd=1.3)
plot(seq(1,iter,1),temp, type="l", col="indianred", 
     main= "Length Parameter (K) Variation", lwd = 3,
     xlab= "Iteration",
     ylab= "K",
     cex.main=1.7, cex.axis=1,cex.lab=1.4)
grid(13,13,lwd=1.3)
