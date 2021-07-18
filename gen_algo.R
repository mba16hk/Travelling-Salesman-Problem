#function to plot data
plotdata<-function(coors){
  plot(coors[,1],coors[,2],ylim=c(min(coors[,2])-0.05,max(coors[,2])+0.05), 
       xlim=c(min(coors[,1])-0.05,max(coors[,1])+0.05),xlab="", ylab="",
       xaxt="n", yaxt="n")
  coors<-rbind(coors,coors[1,])
  
  grid(12,12,lwd=1.1)
  for (i in 1:(nrow(coors)-1)){
    segments(coors[i,1],coors[i,2],coors[i+1,1],
             coors[i+1,2], col="red")}
}

#function to generate a population of tours of a certain size
generate_population<-function(population_size, city_positions){
  population<-matrix(0,nrow=population_size,ncol=length(city_positions))
  for (i in 1:nrow(population)){
    population[i,]<-sample(city_positions,length(city_positions))
  }
  return(population)
}

#function to calculate tour distance
tour_distance<-function(population,coors){
  distances<-c()
  for (chromosome in 1:nrow(population)){
    tour<- coors[population[chromosome,],,drop=FALSE]
    tour<-rbind(tour,tour[1,])
    #print(tour)
    line_len<-0
    for (i in 1:(nrow(tour)-1)){
      line_len[i]<-round(sqrt((tour[i,1]-tour[i+1,1])^2+
                                (tour[i,2]-tour[i+1,2])^2),7)
    }
    distances[chromosome]<-as.numeric(sum(line_len))
  }
  return(distances)
}

#function to introduce mutations by swapping city positions
mutate<-function(tour){
  swap = sample(1:ncol(tour), 2)
  tour[,swap[1]:swap[2]] = tour[,swap[2]:swap[1]]
  return(tour)
}

#function to make crossover children of two tour parents
crossover<-function(population){
  #initialise children matrix
  children<-matrix(0,nrow=nrow(population),ncol=ncol(population))
  start_pos <- seq(1, by = 2, length = nrow(children) / 2)
  for (k in start_pos){
    #specify populations
    tour1<-population[k,]
    tour2<-population[(k+1),]
    #select a section from first tour to keep in child
    keep_from_p1<-sample(tour1[1:(length(tour1)-3)],1)
    start<-which(tour1==keep_from_p1)
    keep_from_p1<-tour1[start:(start+2)]
    #initilise matrix of first child
    child1<-rep(0,length(tour1))
    from_p1_to_child1<-match(keep_from_p1,tour1)
    child1[from_p1_to_child1]<-keep_from_p1
    #fill in the points where the child matrix is 0 from the second tour
    from_p2<-which(child1==0)
    for (i in from_p2){
      for (j in 1:length(tour2)){
        if (!(tour2[j] %in% child1)){
          child1[i]<-tour2[j]
          break
        }
      }
    }
    #do the same for the second child
    keep_from_p2<-sample(tour2[1:(length(tour2)-3)],1)
    start<-which(tour2==keep_from_p2)
    keep_from_p2<-tour2[start:(start+2)]
    child2<-rep(0,length(tour2))
    from_p2_to_child2<-match(keep_from_p2,tour2)
    child2[from_p2_to_child2]<-keep_from_p2
    from_p1<-which(child2==0)
    for (i in from_p1){
      for (j in 1:length(tour1)){
        if (!(tour1[j] %in% child2)){
          child2[i]<-tour1[j]
          break
        }
      }
    }
    #fill in the crossed over children into the children matrix
    children[k:(k+1),]<-rbind(child1,child2)
  }

  return(children)
}

genes<-50
first_chromosome<-random_cities(genes)
first_chromosome<-cbind(first_chromosome,c(1:nrow(first_chromosome)))


#import tsp questions
#coors<-read.table("ulysses16.tsp", sep="")
#first_chromosome<-as.matrix(coors[,2:3])
#first_chromosome<-cbind(first_chromosome,c(1:nrow(first_chromosome)))
#genes<-nrow(first_chromosome)

#determine population size & iteration threshold based on number of cities/genes
if (genes<=15){
  population_size=50
  threshold_iteration<-500
} else if (genes >15 && genes <=23){
  population_size=70
  threshold_iteration<-500
} else if (genes >23 && genes <=60){
  population_size=210
  threshold_iteration<-1600
}else if (genes>60 && genes <=110){
  population_size=250
  threshold_iteration<-7500
}else if (genes>110){
  population_size=30
  threshold_iteration<-25000
}

#create initial population from the generate population function
initial_population<-generate_population(population_size,first_chromosome[,3])
plotdata(first_chromosome[initial_population[1,],,drop=FALSE])

#initialise variables
iteration<-0
best_distances<-c()

#eneter while loop that implements the genetic algorithm
while (iteration<=threshold_iteration){
  
  #create empty vectors for each loop to be filled
  all_distances<-c()
  populations<-c()
  all_pop_distances<-c()
  
  #compute the distances for each tour in the initial population
  population_distances<-round(tour_distance(initial_population, 
                                            first_chromosome),3)
  
  #mutation: mutate 25% of population and measure their tour distances
  mutation_population<-initial_population[1:round(0.25*population_size),]
  mutated_population<- mutate(mutation_population)
  mutated_distances<-round(tour_distance(mutated_population,first_chromosome),3)
  
  #crossover: cross 65% of population and measure their touor distances
  begin<-(nrow(mutated_population)+1)
  end<-(round(0.85*population_size)-1)
  if (((begin %% 2)==0 && (end%%2)==0)|((begin %% 2)==1 && (end%%2)==1)){
    end<-end-1
  }
  crossover_population<-initial_population[begin:end,]
  crossed_population<-crossover(crossover_population)
  cross_distances<-round(tour_distance(crossed_population,first_chromosome),3)
  
  #new chromosomes: introduce comepletely new tours to the population
  new_chromosomes<-generate_population((population_size-end),
                                       first_chromosome[,3])
  newchr_distances<-round(tour_distance(new_chromosomes,first_chromosome),3)
  
  #concatenate all new data: population size has doubled
  all_distances<-c(population_distances,mutated_distances,cross_distances,
                       newchr_distances)
  populations<-rbind(initial_population,mutated_population,crossed_population,
                     new_chromosomes)
  
  #remove the worst tours (with longest distances), and keep the best tours
  all_pop_distances<-cbind(all_distances,populations) #all populations & their distances
  orderd_popdists<-all_pop_distances[order(all_pop_distances[,1]),] #order populations based on distance
  initial_population<-orderd_popdists[1:population_size,2:ncol(orderd_popdists)]
  best_distances[iteration]<-as.numeric(all_pop_distances[1,1])
  #number of tours in the population is reduced initial population size

  if ((iteration %% 1000)==0){
    #pdf(paste0("439_",iteration,".pdf"))
    plotdata(first_chromosome[initial_population[1,],,drop=FALSE])
    print(paste0("here",best_distances[iteration]))
    dev.off()
    Sys.sleep(1)
  }
  
  if ((iteration %% 100)==0){
    print(best_distances[iteration])
  }
  
  iteration=iteration+1
  
}

plot(seq(1,length(best_distances),1),best_distances,
     type="l", col="seagreen", 
     main= "Variation of Tour Distance", lwd = 3,
     xlab= "Iteration",
     ylab= "Tour Distance",
     xaxt = "n",
     yaxt = "n",
     xlim = c(0,700),
     cex.main=1.7, cex.axis=1,cex.lab=1.4)
grid(13,13,lwd=1.3)
axis(1,                         # Define x-axis manually
     at = seq(0,700,116),
     labels = c("0","2000","4000", "6000", "8000", "10000","12000"))
axis(2,                         # Define x-axis manually
     at = seq(5.2,max(best_distances),1.8),
     labels = c("3e5","4e5","5e5", "6e5", "7e5", "8e5", "9e5", "1e6"))
