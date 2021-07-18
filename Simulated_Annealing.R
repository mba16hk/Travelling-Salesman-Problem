#Simulated Annealing
#setwd("C:/Users/Heba/Desktop/MPhil/SP/SP3/Sim_annealing") #directory to save plots
cities <-12 #number of cities

#Different Tours
random_cities<-function(cities){
  x_coor <- runif(cities, 0.2, 1.1)
  y_coor <- runif(cities, 0.2, 1.1)
  positions<-c(1:length(x_coor))
  coors <- cbind(x_coor,y_coor,positions) 
  return(coors)
}

circle_cities<-function(cities,h=0.5,r=0.5){
  theta = runif(cities) * 2 * pi
  x_coor<-h+cos(theta) * r
  y_coor<-h+sin(theta) * r
  positions<-c(1:length(x_coor))
  coors <- cbind(x_coor,y_coor,positions) 
  return(coors)
}

square_cities<-function(cities){
  if (cities == 4){
    x_coor<-c(0.2,0.2,1,1)
    y_coor<-c(0.2,1,0.2,1)
  }else{
    x_coor<-c(0.2,0.2,1,1)
    y_coor<-c(0.2,1,0.2,1)
    switcher<-5
    for (i in 1:(cities-4)){
      if ((switcher%%5)==0){
        x_coor<-append(x_coor,runif(1,0.2,1))
        y_coor<-append(y_coor,0.2)
        switcher=switcher-4
      } else if(switcher==1){
        x_coor<-append(x_coor,runif(1,0.2,1))
        y_coor<-append(y_coor,1)
        switcher= switcher+1
      } else if ((switcher%%2)==0){
        y_coor<-append(y_coor,runif(1,0.2,1))
        x_coor<-append(x_coor,1)
        switcher=switcher+1
      }else if ((switcher%%3)==0){
        y_coor<-append(y_coor,runif(1,0.2,1))
        x_coor<-append(x_coor,0.2)
        switcher=switcher+2
      }
    }
  }
  positions<-c(1:length(x_coor))
  coors <- cbind(x_coor,y_coor,positions) 
  return(coors)
}

#create a distance function that measures journey lengths between each city
distance = function(coors2) {
  coors2<-rbind(coors2,coors2[1,])
  for (i in 1:(nrow(coors2)-1)){
    line_len[i]<-round(sqrt((coors2[i,1]-coors2[i+1,1])^2+
                              (coors2[i,2]-coors2[i+1,2])^2),7)
  }
  return(line_len)
}

#function for creating a plot of the randomised tours
plotdata<-function(coors2){
  plot(coors2[,1],coors2[,2],ylim=c(min(coors2[,2])-0.05,max(coors2[,2])+0.05), 
       xlim=c(min(coors2[,1])-0.05,max(coors2[,1])+0.05),xlab="", ylab="")
  coors2<-rbind(coors2,coors2[1,])
  grid(12,12,lwd=1.1)
  for (i in 1:(nrow(coors2)-1)){
    segments(coors2[i,1],coors2[i,2],coors2[i+1,1],
             coors2[i+1,2], col="red")}
}


#coors2<-square_cities(cities)
coors2<-read.table("pr439.tsp", sep="") #could import coordinates from directory
coors2<-(coors2[,2:3])
positions<-c(1:nrow(coors2))
coors2 <- as.matrix(cbind(coors2,positions))
cities<-nrow(coors2)
dev.off()
plot(coors2[,1],coors2[,2], ylim=c(min(coors2[,2])-0.05,max(coors2[,2])+0.05),
     xlim=c(min(coors2[,1])-0.05,max(coors2[,1])+0.05),
     xlab="", ylab="")
grid(12,12,lwd=1.1)

#draw lines between coordinates
locations<-matrix(0,nrow=length(coors2[,1]),ncol=3)
start<-as.numeric(coors2[sample(nrow(coors2),1),])
x2 <- coors2[-which(coors2[,3]==start[3]),]

#create matrix of paths to follow through sampled city point
locations[1,]=coors2[which(coors2[,3]==start[3]),]

#all_cities<-seq(2,nrow(coors),1) #remaining unsampled cities
all_cities<-seq(2,nrow(coors2),1)
for (i in all_cities){
    if (i==length(all_cities)+1){
      locations[i,]=as.numeric(x2)
      break
    }
    follow<-as.numeric(x2[sample(nrow(x2),1),])
    locations[i,]=x2[which(x2[,3]==follow[3]),]
    x2 <- x2[-which(x2[,3]==follow[3]),]
}

#create a random tour through city coordinates
tour<-locations[,3]
rand_tour<-rbind(locations,as.numeric(locations[1,]))

#measure the euclidean distance between each random walk
line_len<-c() #journey lengths between each city stored here
line_len<-distance(coors2)

for (i in 1:(nrow(rand_tour)-1)){
  segments(rand_tour[i,1],rand_tour[i,2],rand_tour[i+1,1],rand_tour[i+1,2],
           col="red")
  line_len[i]<-round(sqrt((rand_tour[i,1]-rand_tour[i+1,1])^2+
                            (rand_tour[i,2]-rand_tour[i+1,2])^2),7)
}

#Initialise params
max_iter=35000;loops = 1000000; orig_temp = 2000

#other initialised matrices and vectors
ratios<-c() #to observed decaying probability of accepting new tours
temperatures<-c() # to observe temperature changes
tour_len<-sum(line_len) #a vector of tour lengths of all tours
candidates<-list(rand_tour) #list of best tour length and corresponding tour
names(candidates)<-c(tour_len) #name list by the best tour length
tour_distance<-tour_len[1] #initialise tour distance as the firs tour length
best_distance<-tour_distance #the best distance is the shortest tour length
best_tour<-tour #best tour is tour with shortest tour length

#iterate until GLOBAL minimum is found
for (m in 1:loops) {
   iteration=m
   temp=orig_temp*(1/(1+exp((iteration-max_iter)/max_iter)))#annealing schedule
   temperatures[m]<-temp
   candidate_tour=tour
   swap = sample(candidate_tour, 2)
   candidate_tour[swap[1]:swap[2]] = rev(candidate_tour[swap[1]:swap[2]])
   n<-0; locations2<-matrix(0,nrow=length(coors2[,1]),ncol=3)
   for(i in candidate_tour){
     n=n+1
     locations2[n,]=locations[which(locations[,3]==i),]
   }
   rand_tour2<-rbind(locations2,as.numeric(locations2[1,]))
   line_len<-c() #empty the vector
   line_len<-distance(rand_tour2)
   
   #change plot every given number of iterations
   if (iteration %% 50000 == 0){
     plotdata(candidates[[1]])
     print(best_distance)
     Sys.sleep(1)
   }
   this_tour<-sum(line_len)
   tour_len[length(tour_len)+1]<-sum(line_len)
   
   if (temp > 0) {
     delta<-(this_tour-tour_distance)
     ratio = exp(-(delta)/ temp)
     ratios[length(ratios)+1]<-ratio
   } else {
     ratio = 0
     ratios[length(ratios)+1]<-ratio
   }
   
   if (this_tour<best_distance){
     tour<-candidate_tour
     tour_distance<-this_tour
     best_tour<-tour
     best_distance<-tour_distance
     candidates<-c(list(rand_tour2))
     names(candidates)<-c(best_distance)
   } else{
     
     if (runif(1)<ratio) {
           tour<-candidate_tour
           tour_distance<-this_tour
      }
   }
   rand_tour2<-c() #empty the vector
 }


#plot Acceptance probability and annealing temperature
x=seq(1,length(tour_len),100)
plot(x,ratios[x],ylim=c(0,1), pch=".", col="indianred", 
     main= "Acceptance Probability",
     xlab= "Iteration",
     ylab= "Probability",
     cex.main=1.7, cex.axis=1,cex.lab=1.4)
grid(10,10,lwd=1.3)
plot(x,temperatures[x], ylim=c(0,1500), type="l", col="blue", 
     main= "Temperature Annealing",
     xlab= "Iteration",
     ylab= "Temperature", cex.main=1.7, cex.axis=1,cex.lab=1.4,
     lwd=2)
grid(10,10,lwd=1.3)
plot(x,tour_len[x], pch=".", col="seagreen",
     main= "Variation of Tour Length",
     xlab= "Iteration",
     ylab= "Tour Distance", cex.main=1.7, cex.axis=1,cex.lab=1.4)
grid(10,10,lwd=1.3)
