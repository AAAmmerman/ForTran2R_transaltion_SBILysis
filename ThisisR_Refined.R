library(tidyverse) # allows pipeline and neat wrangling
library(zeallot) # allows compact way to assign variables

#   N   O   T   E   S  /   W  rk  In   P  rg

# turns out source of error was from the randomchoice function
# an addition of an intermediate variable seems to have done the trick





# do _ = __,__ all turned into for loops. 
# Gillspie Algo is repeat until t > tfinal
# see ThisisR file for additional changes on first translation


#             F    U    N   C   T    I   O   N   S


# minvect : intention is to find the location of the minimum
#           entry of a vector excluding 0
#           ..............................
#           in the OG code, this pics the first minimum
#           since its basis is whether the next value is
#           .lt. (lesser than); whereas min[i] < min[i-1]
#           is false.   We can set it up to randomly pick
#           from all value == min(vector), but for now we
#           can superscede this by using the following r
#           commands:        
#                           
#                           thismin <- min(vector[which(vector>0)])
#                           which(vector==thismin)[1]
#                           or more nicely
#
#                           thismin <- which(tester>0 & 
#                                       is.na(tester)==F) %>% 
#                                       vector[.] %>%  
#                                       min()  
#                                         
#                           loca <- which(vector==thismin)[1]
#           this will become the new minvect function

minvect <- function(g){
  g2 <- c(g)
  
  thismin <- which(  g2 > 0   # &  is.na(g2)==F  
                              # seems which filters out alrdy
                     ) %>% 
              g2[.] %>% 
              min()
  loca <- which(g2==thismin)[1]
  return(loca)
}

# findint: intention is to find the number of times a value
#           shows up in a given vector.
#           This is supersceded by the following R code:
#
#                             sum(some_vector==the_value, 
#                                           na.rm = T)
#           hence:

findint <- function(g, the_value){
  some_vector <- c(g)
  
  this_many_times <- sum(some_vector==the_value, 
                          na.rm = T)
  return(this_many_times)
}


# findfirst: simply to find the first instance of a 
#             value in a vector. this can be done via
#             the followig r code:
#                                   
#                 loca <- which(vector==this_value)[1]    
#                                   
#             hence:

findfirst <- function(g,numspec){
  this_vector <- c(g)
  this_value <- numspec
  loca <- which(this_vector == this_value)[1]
  return(loca)
}


# randchoice: pick a random location of the specified value
#             likely can be supersceded by the following
#
#             locas <- which(vector==this_value) # all with this
#             randomed_to_this <- sample(locas,1)
#             
#             e.g. vector = c( 1 , 2 , 3 , 2, 5, 2)
#             which(blah = 3) saves a list of c( 2 , 4 ,6)
#             sample(which(blah = 3),1), picks one of them
#             each has an equal chance of being sampled...
#             
#             see:sample(which(tester==3),1000,replace = T) %>% table()/1000
# .
#            loca?     |   2    |  4    | 6 
#                       ..................
#         picked/1000  | 0.318 | 0.336 | 0.346  # made pretty
#             
#             hence

randchoice <- function(g,numspec){
  this_vector <- g
  this_value <- numspec
  
  locas <- which(this_vector == this_value) # all with this
  randomed_to_this <- sample(1:length(locas),size = 1, replace = T) #%>% locas[.]
  return_this <- locas[randomed_to_this]
  
  return(return_this)
}

# movetpa: this one i will not attempt to condense beyond the basics
#         From what I understand, it does the following:
#         Only output is p_rebind
#           
#                 find the number of values in this matrix equal to 5
#                 find the number of values in FIRST COLUMN equal to 1
#                 find the number of values in NOT IN FIRST COLUMN equal to 1
#
#                 exposed = number based on those outputs (and probabilities)
#                 Other hardcoded values? 602.2 oH micro M o l a r
#                 some calculation that results in p_rebind
#         
#                 thankfully we can superscede most of the forloops with this:
#
#                 sum(this_matrix==this_value,na.rm = T)
#                 example:
#
#                 ok <- matrix(c(1,2,3), nrow = 100,ncol = 3) #300 entires
#
#                 > sum(ok==3,na.rm = T)
#                 [1] 100

movetpa <- function(state){
  this_matrix <- state
  #necessary to direct to global env. to retrive set value
  prob_N02n <- .GlobalEnv$prob_N02n
  prob_N00n <- .GlobalEnv$prob_N00n
  prob_N02i <- .GlobalEnv$prob_N02i
  prob_N00i <- .GlobalEnv$prob_N00i
  param9 <- .GlobalEnv$param9
  ktPAon <- .GlobalEnv$ktPAon
  vol3 <- .GlobalEnv$vol3
  D <- .GlobalEnv$D
  Ntot<- .GlobalEnv$Ntot
  nodes<- .GlobalEnv$nodes
  countt <- .GlobalEnv$countt
  
  countstate6 <- sum(this_matrix == 5 ,na.rm = T)
  
  countstate7 <- sum(this_matrix[ ,1] == 1 ,na.rm = T) #first col only
  countstate8 <- sum(this_matrix[ ,-1] == 1 ,na.rm = T) #EXCLUDE col only
  
  exposed = countstate6 + 
    (prob_N02i + 2*prob_N00i)*countstate7 + 
    (prob_N02n + 2*prob_N00n)*countstate8
  
  bs0 = exposed/vol3/602.2
  
  r0 = 0.5    
  r1 = 22.5   
  k0 = ktPAon*bs0  #[F01]
  
  p_rebind = 1 - 
    r1/r0*
    (1/(1 + (r1 - r0)*sqrt(k0/D)*
          cosh(sqrt(k0/D)*r0)/sinh(sqrt(k0/D)*r0)))
              
  rr = runif(1)
  
  if(rr<=p_rebind){
    stop("Problem: random number smallerer than p_rebind. 
            must add rebinding to code. issue orig from funct: movetpa")
  }
  
  return(p_rebind)
}

# movepli: same issue as movetpa, wont condense unless im absolutely sure
#           code flow is:
#
#           for each connected ind, find the number of times
#           of the corresponding ROW in STATE == some value
#           
#           IF, the entry in the 1st column of the matrix STATE AND THE 
#           corresponding
#           ROW equal to THAT connected ind.== 1,  this value, else that value
#           
#           add up tot_sitesPLi based on those outputs and probabilities for 
#           each ind.
#           
#           then, pick a random site from tot_sitesPLi
#           this for-loop seems unnecessary but just in case I'll leave it be

movepli <- function(s, this, st){
  loc <- this
  state <- s
  statetemp <- st
  
  PLilocation<- .GlobalEnv$PLilocation
  PLilocationt <- .GlobalEnv$PLilocationt
  connect_ind <- .GlobalEnv$connect_ind
  prob_N02n <- .GlobalEnv$prob_N02n
  prob_N00n <- .GlobalEnv$prob_N00n
  prob_N02i <- .GlobalEnv$prob_N02i
  prob_N00i <- .GlobalEnv$prob_N00i
  
  alpha <- .GlobalEnv$alpha
  Ntot <- .GlobalEnv$Ntot
  nodes <- .GlobalEnv$nodes
  sizenonzero_connectind <- .GlobalEnv$sizenonzero_connectind # may no longer 
          
  newsite = 0
  newdoublet = 0
  newsite2 = 0
  newdoublet2 = 0
  cumprob = 0
  
  tot_sitesPLi <- c()#rep(dbl,6)
  prob_sitePLi <- c()#rep(dbl,6)
  
  sites_N03 <- sites_N10 <- sites_Nplg <- 
    sites_Nplgn <- sites_Nplgi <- sites_degrade <- 
    c()#rep(0,6)
  
  PLilocationt = PLilocation #not sure why its here 

  
  for (j in 1:length(connect_ind)) {
    # oktest <- connect_ind[j]
    sites_N10[j] <-     sum(state[connect_ind[j],  ]     ==  3, na.rm = T)
    sites_N03[j] <-     sum(state[connect_ind[j],  ]     ==  5, na.rm = T)
    sites_degrade[j] <- sum(state[connect_ind[j],  ]     == -1, na.rm = T)
    sites_Nplgn[j] <-   sum(state[connect_ind[j],-1]     ==  1, na.rm = T)
    
    
    if (state[connect_ind[j], 1] == 1) {
      sites_Nplgi[j]<- 1
    } else {
      sites_Nplgi[j]<- 0
    }
    
    tot_sitesPLi[j] <- sites_N10[j] + sites_N03[j] + 
      (prob_N02n + 2*prob_N00n)*
      sites_Nplgn[j] +
      (prob_N02i + 2*prob_N00i)*
      sites_Nplgi[j] + 
      alpha*sites_degrade[j]  #  total number of empty binding sites 4 PLi 
    #  at location j.
    
  }
  
  #if(sum())
  prob_sitePLi <- tot_sitesPLi/sum(tot_sitesPLi)
  newsite2 <- sample(1:length(prob_sitePLi),1,prob = prob_sitePLi, replace = T)
  
  P <- c()
  Pnsmall <- c()
  Pismall <- c()
  Pismall_temp <- c()
  Pnsmall_temp <- c()
  
  for (j in 1:length(connect_ind)) {
    if (newsite2 == j) {
      P[1] <- sites_N10[j]/tot_sitesPLi[j]
      P[2] <- (prob_N02n + 2*prob_N00n)*
        sites_Nplgn[j]/tot_sitesPLi[j]
      P[3] <- (prob_N02i + 2*prob_N00i)*
        sites_Nplgi[j]/tot_sitesPLi[j]
      P[4] <- sites_N03[j]/tot_sitesPLi[j]
      P[5] <- alpha*sites_degrade[j]/tot_sitesPLi[j]
    }else if (newsite2 > sizenonzero_connectind) {
      stop("Problem: newsite2 larger than possible")
      break
    }
  }
  
  P_temp <- P/sum(P)
  newdoublet2 <- sample(1:length(P_temp),1,prob = P_temp, replace = T)
  
  if (newdoublet2 == 1) {
    
    N10toN13 <- findfirst(
      state[connect_ind[newsite2], ] # of this
      , 3                            # find 1st this
    )
    state[connect_ind[newsite2], N10toN13] <- 4
    
    PLilocationt[loc,1] <- connect_ind[newsite2]
    PLilocationt[loc,2] <- N10toN13
    
  } else if (newdoublet2 == 2) {
    
    NplgtoN23 <- findfirst(
      state[connect_ind[newsite2],-1]
      ,1
    )
    NplgtoN23 = 1 + NplgtoN23  
    
    if (NplgtoN23 == 1){
      stop("Problem: NplgtoN23 == 1.")
    }
    
    Pnsmall[1] <- prob_N00n/(prob_N02n + prob_N00n)
    Pnsmall[2] <- prob_N02n/(prob_N02n + prob_N00n)
    
    Pnsmall_temp <- Pnsmall/sum(Pnsmall) #just in case
    temp <- sample(1:length(Pnsmall_temp),1,prob = Pnsmall_temp, replace = T)
    
    if (temp == 1) {
      state[connect_ind[newsite2],NplgtoN23] <- 5
    }else if (temp == 2) {
      state[connect_ind[newsite2],NplgtoN23] <- 6
    }
    
    PLilocationt[loc,1] <- connect_ind[newsite2]
    PLilocationt[loc,2] <- NplgtoN23
    
  } else if (newdoublet2 == 3) {
    
    NplgtoN23 = 1 # only place we have intact fibrin is in column 1
    Pismall[1] = prob_N00i/(prob_N02i + prob_N00i)
    Pismall[2] = prob_N02i/(prob_N02i + prob_N00i)
    
    Pismall_temp <- Pismall/sum(Pismall) #just in case
    temp <- sample(1:length(Pismall_temp),1,prob = Pismall_temp, replace = T)
    
    if (temp == 1) {
      state[connect_ind[newsite2],NplgtoN23] <- 5
    }else if(temp == 2) {
      state[connect_ind[newsite2],NplgtoN23] <- 6
    }
    
    PLilocationt[loc,1] <- connect_ind[newsite2]
    PLilocationt[loc,2] <- NplgtoN23
    
  } else if(newdoublet2 == 4) {
    
    N03toN33 <- findfirst(state[connect_ind[newsite2], ],5 )
    state[connect_ind[newsite2], N03toN33] <- 7
    
    PLilocationt[loc,1] <- connect_ind[newsite2]
    PLilocationt[loc,2] <- N03toN33
    
  }else if(newdoublet2 == 5) {
    
    degradetoPLi <- findfirst(state[connect_ind[newsite2], ], -1 )
    state[connect_ind[newsite2], degradetoPLi] <- 9
    
    PLilocationt[loc,1] <- connect_ind[newsite2]
    PLilocationt[loc,2] <- degradetoPLi
    
  } 
  
  statetemp <<- state
  PLilocationt <<- PLilocationt
  
}


move2pli <- function(s, this, this2, st){
  loc <- this
  loc2 <- this2
  state <- s
  statetemp <- st
  
  PLilocation<- .GlobalEnv$PLilocation
  PLilocationt <- .GlobalEnv$PLilocationt
  connect_ind <- .GlobalEnv$connect_ind
  prob_N02n <- .GlobalEnv$prob_N02n
  prob_N00n <- .GlobalEnv$prob_N00n
  prob_N02i <- .GlobalEnv$prob_N02i
  prob_N00i <- .GlobalEnv$prob_N00i
  
  alpha <- .GlobalEnv$alpha
  Ntot <- .GlobalEnv$Ntot
  nodes <- .GlobalEnv$nodes
  sizenonzero_connectind <- .GlobalEnv$sizenonzero_connectind # may no longer 
  
  newsite = 0
  newdoublet = 0
  newsite2 = 0
  newdoublet2 = 0
  cumprob = 0
  
  tot_sitesPLi <- c()#rep(dbl,6)
  prob_sitePLi <- c()#rep(dbl,6)
  
  sites_N03 <- sites_N10 <- sites_Nplg <- 
    sites_Nplgn <- sites_Nplgi <- sites_degrade <- 
    c()#rep(0,6)
  
  PLilocationt = PLilocation #not sure why its here 
  # be used
  
  for (j in 1:length(connect_ind)) {
    # oktest <- connect_ind[j]
    sites_N10[j] <-     sum(state[connect_ind[j],  ]     ==  3, na.rm = T)
    sites_N03[j] <-     sum(state[connect_ind[j],  ]     ==  5, na.rm = T)
    sites_degrade[j] <- sum(state[connect_ind[j],  ]     == -1, na.rm = T)
    sites_Nplgn[j] <-   sum(state[connect_ind[j],-1]     ==  1, na.rm = T)
    
    
    if (state[connect_ind[j], 1] == 1) {
      sites_Nplgi[j]<- 1
    } else {
      sites_Nplgi[j]<- 0
    }
    
    tot_sitesPLi[j] <- sites_N10[j] + sites_N03[j] + 
      (prob_N02n + 2*prob_N00n)*
      sites_Nplgn[j] +
      (prob_N02i + 2*prob_N00i)*
      sites_Nplgi[j] + 
      alpha*sites_degrade[j]  #  total number of empty binding sites 4 PLi 
    #  at location j.
    
  }
  
  prob_sitePLi <- tot_sitesPLi/sum(tot_sitesPLi)
  newsite <- sample(1:length(prob_sitePLi),1,prob = prob_sitePLi, replace = T)
  
  P <- c()
  Pnsmall <- c()
  Pismall <- c()
  Pismall_temp <- c()
  Pnsmall_temp <- c()
  
  for (j in 1:length(connect_ind)) {
    if (newsite == j) {
      P[1] <- sites_N10[j]/tot_sitesPLi[j]
      P[2] <- (prob_N02n + 2*prob_N00n)*
        sites_Nplgn[j]/tot_sitesPLi[j]
      P[3] <- (prob_N02i + 2*prob_N00i)*
        sites_Nplgi[j]/tot_sitesPLi[j]
      P[4] <- sites_N03[j]/tot_sitesPLi[j]
      P[5] <- alpha*sites_degrade[j]/tot_sitesPLi[j]
    }else if (newsite > sizenonzero_connectind) {
      stop("Problem: newsite larger than possible")
      break
    }
  }
  
  P_temp <- P/sum(P)
  newdoublet <- sample(1:length(P_temp),1,prob = P_temp, replace = T)
  
  if (newdoublet == 1) {
    
    N10toN13 <- findfirst(
      state[connect_ind[newsite], ] # of this
      , 3                            # find 1st this
    )
    state[connect_ind[newsite], N10toN13] <- 4
    
    PLilocationt[loc,1] <- connect_ind[newsite]
    PLilocationt[loc,2] <- N10toN13
    
  } else if (newdoublet == 2) {
    
    NplgtoN23 <- findfirst(
      state[connect_ind[newsite],-1]
      ,1
    )
    NplgtoN23 = 1 + NplgtoN23  
    
    if (NplgtoN23 == 1){
      stop("Problem: NplgtoN23 == 1.")
    }
    
    Pnsmall[1] <- prob_N00n/(prob_N02n + prob_N00n)
    Pnsmall[2] <- prob_N02n/(prob_N02n + prob_N00n)
    
    Pnsmall_temp <- Pnsmall/sum(Pnsmall) #just in case
    temp <- sample(1:length(Pnsmall_temp),1,prob = Pnsmall_temp, replace = T)
    
    if (temp == 1) {
      state[connect_ind[newsite],NplgtoN23] <- 5
    }else if (temp == 2) {
      state[connect_ind[newsite],NplgtoN23] <- 6
    }
    
    PLilocationt[loc,1] <- connect_ind[newsite]
    PLilocationt[loc,2] <- NplgtoN23
    
  } else if (newdoublet == 3) {
    
    NplgtoN23 = 1 # only place we have intact fibrin is in column 1
    Pismall[1] = prob_N00i/(prob_N02i + prob_N00i)
    Pismall[2] = prob_N02i/(prob_N02i + prob_N00i)
    
    Pismall_temp <- Pismall/sum(Pismall) #just in case
    temp <- sample(1:length(Pismall_temp),1,prob = Pismall_temp, replace = T)
    
    if (temp == 1) {
      state[connect_ind[newsite],NplgtoN23] <- 5
    }else if(temp == 2) {
      state[connect_ind[newsite],NplgtoN23] <- 6
    }
    
    PLilocationt[loc,1] <- connect_ind[newsite]
    PLilocationt[loc,2] <- NplgtoN23
    
  } else if(newdoublet == 4) {
    
    N03toN33 <- findfirst(state[connect_ind[newsite], ],5 )
    state[connect_ind[newsite], N03toN33] <- 7
    
    PLilocationt[loc,1] <- connect_ind[newsite]
    PLilocationt[loc,2] <- N03toN33
    
  }else if(newdoublet == 5) {
    
    degradetoPLi <- findfirst(state[connect_ind[newsite], ], -1 )
    state[connect_ind[newsite], degradetoPLi] <- 9
    
    PLilocationt[loc,1] <- connect_ind[newsite]
    PLilocationt[loc,2] <- degradetoPLi
    
  } 
  
  # ROUND 2 #########################################
  
  for (j in 1:length(connect_ind)) {
    # oktest <- connect_ind[j]
    sites_N10[j] <-     sum(state[connect_ind[j],  ]     ==  3, na.rm = T)
    sites_N03[j] <-     sum(state[connect_ind[j],  ]     ==  5, na.rm = T)
    sites_degrade[j] <- sum(state[connect_ind[j],  ]     == -1, na.rm = T)
    sites_Nplgn[j] <-   sum(state[connect_ind[j],-1]     ==  1, na.rm = T)
    
    
    if (state[connect_ind[j], 1] == 1) {
      sites_Nplgi[j]<- 1
    } else {
      sites_Nplgi[j]<- 0
    }
    
    tot_sitesPLi[j] <- sites_N10[j] + sites_N03[j] + 
      (prob_N02n + 2*prob_N00n)*
      sites_Nplgn[j] +
      (prob_N02i + 2*prob_N00i)*
      sites_Nplgi[j] + 
      alpha*sites_degrade[j]  #  total number of empty binding sites 4 PLi 
    #  at location j.
    
  }
  
  prob_sitePLi <- tot_sitesPLi/sum(tot_sitesPLi)
  newsite2 <- sample(1:length(prob_sitePLi),1,prob = prob_sitePLi, replace = T)
  
  P <- c()
  Pnsmall <- c()
  Pismall <- c()
  Pismall_temp <- c()
  Pnsmall_temp <- c()
  
  for (j in 1:length(connect_ind)) {
    if (newsite2 == j) {
      P[1] <- sites_N10[j]/tot_sitesPLi[j]
      P[2] <- (prob_N02n + 2*prob_N00n)*
        sites_Nplgn[j]/tot_sitesPLi[j]
      P[3] <- (prob_N02i + 2*prob_N00i)*
        sites_Nplgi[j]/tot_sitesPLi[j]
      P[4] <- sites_N03[j]/tot_sitesPLi[j]
      P[5] <- alpha*sites_degrade[j]/tot_sitesPLi[j]
    }else if (newsite2 > sizenonzero_connectind) {
      stop("Problem: newsite2 larger than possible")
      break
    }
  }
  
  P_temp <- P/sum(P)
  newdoublet2 <- sample(1:length(P_temp),1,prob = P_temp, replace = T)
  
  if (newdoublet2 == 1) {
    
    N10toN13 <- findfirst(
      state[connect_ind[newsite2], ] # of this
      , 3                            # find 1st this
    )
    state[connect_ind[newsite2], N10toN13] <- 4
    
    PLilocationt[loc2,1] <- connect_ind[newsite2]
    PLilocationt[loc2,2] <- N10toN13
    
  } else if (newdoublet2 == 2) {
    
    NplgtoN23 <- findfirst(
      state[connect_ind[newsite2],-1]
      ,1
    )
    NplgtoN23 = 1 + NplgtoN23  
    
    if (NplgtoN23 == 1){
      stop("Problem: NplgtoN23 == 1.")
    }
    
    Pnsmall[1] <- prob_N00n/(prob_N02n + prob_N00n)
    Pnsmall[2] <- prob_N02n/(prob_N02n + prob_N00n)
    
    Pnsmall_temp <- Pnsmall/sum(Pnsmall) #just in case
    temp <- sample(1:length(Pnsmall_temp),1,prob = Pnsmall_temp, replace = T)
    
    if (temp == 1) {
      state[connect_ind[newsite2],NplgtoN23] <- 5
    }else if (temp == 2) {
      state[connect_ind[newsite2],NplgtoN23] <- 6
    }
    
    PLilocationt[loc2,1] <- connect_ind[newsite2]
    PLilocationt[loc2,2] <- NplgtoN23
    
  } else if (newdoublet2 == 3) {
    
    NplgtoN23 = 1 # only place we have intact fibrin is in column 1
    Pismall[1] = prob_N00i/(prob_N02i + prob_N00i)
    Pismall[2] = prob_N02i/(prob_N02i + prob_N00i)
    
    Pismall_temp <- Pismall/sum(Pismall) #just in case
    temp <- sample(1:length(Pismall_temp),1,prob = Pismall_temp, replace = T)
    
    if (temp == 1) {
      state[connect_ind[newsite2],NplgtoN23] <- 5
    }else if(temp == 2) {
      state[connect_ind[newsite2],NplgtoN23] <- 6
    }
    
    PLilocationt[loc2,1] <- connect_ind[newsite2]
    PLilocationt[loc2,2] <- NplgtoN23
    
  } else if(newdoublet2 == 4) {
    
    N03toN33 <- findfirst(state[connect_ind[newsite2], ],5 )
    state[connect_ind[newsite2], N03toN33] <- 7
    
    PLilocationt[loc2,1] <- connect_ind[newsite2]
    PLilocationt[loc2,2] <- N03toN33
    
  }else if(newdoublet2 == 5) {
    
    degradetoPLi <- findfirst(state[connect_ind[newsite2], ], -1 )
    state[connect_ind[newsite2], degradetoPLi] <- 9
    
    PLilocationt[loc2,1] <- connect_ind[newsite2]
    PLilocationt[loc2,2] <- degradetoPLi
    
  } 
  
  statetemp <<- state
  PLilocationt <<- PLilocationt
  
  
}

#####################################################
#####################################################
######################################################


###############################

#############################


c(
  nodes , Nplginit , simulations , sV 
  , seed , lysunit , tPAunit , PLiunit , endunit 
  , lysunit2 , tPAunit2 , PLiunit2 , endunit2 , lysunit3 
  , tPAunit3 , PLiunit3 , endunit3 , lysunit4 , tPAunit4 
  , PLiunit4 , endunit4 , lysunit5 , tPAunit5 , PLiunit5 
  , endunit5 , lysunit6 , tPAunit6 , PLiunit6 , endunit6 
  , perdegunit , tunit , t2unit , perunit , ptunit 
  , countunit , tvectunit , sunit , plgunit , ctunit 
  , plgunbdunit , plgbdunit , plitimeunit , tPAunbdunit , tPAPLiunit 
  , s2unit , fpunit
) %<-% as.integer(c(
  7 , 1 , 1000 , 26  
  , 0 , 22, 23, 24, 25
  , 26, 27, 28, 29, 30
  , 31, 32, 33, 34, 35
  , 36, 37, 38, 39, 40
  , 41, 42, 43, 44, 45
  , 20, 21, 19, 18, 17
  , 16, 46, 47, 48, 49
  , 50, 51, 52, 53, 54
  , 55, 56)) # this is more compact... might try to find a way to keep
# it like this using an intermediate variable to change 
# respective parameters per "walk"

Ninit <-  as.integer(5*Nplginit)
Ntot <-  as.integer(Ninit + Nplginit )


c(radius , KdtPAyesplg , KdtPAnoplg 
  , KdPLGintact , KdPLGnicked , ktPAon 
  , kplgon , freeplg , kdeg 
  , kplioff , kapcat , kncat 
  , snap_proportion) %<-% c(72.7/2/1000 , 0.02 , 0.36 
                            , 38 , 2.2 , 0.1 
                            , 0.1 , 2 , 5 
                            , 57.6 , 0.1 , 5 ,
                            2.0/3.0 )


# est_mts <-      ## no longer needed, R can re-dim a vector per added entry
#   as.integer(1000000) 
#let this represent estimated max timesteps

dvl <- 0 %>% as.integer()

#c(V, 
#V_tmp,  # see note [2]
# q <-    # see note [3]
#)%<-% list(
#matrix(dvl,sV,9), 
#matrix(dvl,9*sV,1), # see note [2]
#matrix(dvl,9,9)
#)

# [2]   V_tmp was just a temporary matrix (list) s.t. it could be
#       be used to construct V, a matrix of dimensions 9, sV.
#       Let A = 
#         [ Nplg,N12,N10,N13,N03,
#         N23,N33,degradation,
#         degradation plus 1 PLi ] =
#
#         [v_1,...,v_2]
#
# Let each column, v_i, be a vector that denotes a state
# Let row represent a different reaction

#  [3]  q is currently already a 9x9 matrix of 0s (integrs)
#       the method of making a diagonal matrix has been superceeded
#       by the function "diag"
#
#  [3.b]  Diagonal from library(Matrix) may be more efficient
#         computation wise in the future

q <- diag(9) 

V_tmp <- as.integer(
  c(
    0, -1, 1, 0, 0, 0, 0, 0, 0
    ,0, 1, -1, 0, 0, 0, 0, 0, 0
    ,0, 0, 1, -1, 0, 0, 0, 0, 0
    ,0, 0, 0, -1, 1, 0, 0, 0, 0
    ,0, 0, 0, 0, -1, 1, 0, 0, 0
    ,0, 0, 0, 0, 1, -1, 0, 0, 0
    ,0, 0, 0, 0, 1, 0, -1, 0, 0
    ,1, -1, 0, 0, 0, 0, 0, 0, 0
    ,1, 0, -1, 0, 0, 0, 0, 0, 0
    ,1, 0, 0, 0, -1, 0, 0, 0, 0
    ,1, 0, 0, 0, 0, -1, 0, 0, 0
    ,0, 0, 0, 0, 0, 0, 0, 0, 0
    ,0, 0, 0, 0, 0, 0, 0, 1, -1
    ,0, 0, 0, -1, 0, 0, 0, 1, 0
    ,0, 0, 0, 0, -1, 0, 0, 1, 0
    ,0, 0, 0, 0, 0, -1, 0, 1, 0
    ,0, 0, 0, 0, 0, 0, -1, 1, 0
    ,0, 0, 0, 0, 0, 0, 0, 0, 0
    ,0, 0, 0, 0, 0, 0, 0, 0, 0
    ,0, 0, 0, 0, 0, 0, 0, 0, 0
    ,0, 0, 0, 0, 0, 0, 0, 0, 0
    ,0, 0, 0, 0, 0, 0, 0, 0, 0
    ,0, 0, -1, 0, 0, 0, 0, 1, 0
    ,0, -1, 0, 0, 0, 0, 0, 1, 0
    ,-1, 0, 0, 0, 0, 0, 0, 1, 0
    ,-1, 0, 0, 0, 0, 0, 0, 1, 0
  )
)

is.wholenumber <-
  function(x, tol = 
             .Machine$double.eps^0.5)  abs(x - round(x)) < tol
# checking if V_tmp is whole num when div.by 9
# (length(V_tmp)/9) %>% is.wholenumber() #TRUE
# so they want it to be each row is a reachtion and each column
# is a state

V <- matrix(V_tmp,sV,9, byrow = T)

#lets check how it looks
# print(format(V, digits = 3, justify = "right"))
#      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8]
# [1,] " 0" "-1" " 1" " 0" " 0" " 0" " 0" " 0"
# [2,] " 0" " 1" "-1" " 0" " 0" " 0" " 0" " 0"
# [3,] " 0" " 0" " 1" "-1" " 0" " 0" " 0" " 0"
# [4,] " 0" " 0" " 0" "-1" " 1" " 0" " 0" " 0" ...
#       ...  # yep there it is.

PLilocsmall <- 
  currentPLib <-                     # see [30]
  rowtPA <- coltPA<- rep(dvl,2)
#matrix(dvl,2,1) 

tPANplgn <- Nplgn <- rep(dvl,5)
#matrix(dvl,5,1) 




#connect_ind <- 
connect_ind2temp <- 
  # sites_N03 <- 
  # sites_N10 <- 
  # sites_Nplg <- 
  # sites_Nplgn <-  #see [35]
  # sites_Nplgi <- 
  # sites_degrade <- 
  # matrix(dvl,6,1) 
  rtest <- ctest <- testr <- testc <- 
  tot_sitestPA <- prob_sitetPA <- sitestPA_N03 <- sitestPA_N10 <- 
  sitestPA_Nplg <- sitestPA_Nplgn <- sitestPA_Nplgi <-
  #matrix(dvl,6,1) 
  rep(dvl,6)



newstate<- rep(dvl,9)
#matrix(dvl,9,1) 




# ones <-                 see [1]
# PLiwaste <-             #see [14]
#row <- col <-          # see [27]
rowPLi <- colPLi <- rep(dvl,300)
#matrix(dvl,300,1) 

# [1] ones is just a matrix of ones so let us 
#     recreate that matrix with the following
#     whereas ones is now an integer matrix
#     of 300 rows and 1 column

# ones <- rep(as.integer(1),300)#matrix(as.integer(1),300,1) #no longer used?



bigJ <- c()#rep(dvl,12000)
#matrix(dvl,12000,1) 




#reaction <- # see [9]
#tPA <- tPA_loc <-        #see [15] [16]
#PLi <- # see [11]
#per_lysis <- [18]
#   degrade_tot <- undegrade_tot <- undegrade_new <- 
# PLi_loc <- 
#   # Anumber <-            # see [10]
#                           matrix(dvl,est_mts,1) 
##  all rest [23]



# PLilocation <- #  see [13]
PLilocationt <- 
  matrix(dvl,300, 2) 



connect_ind2 <- 
  matrix(dvl,5, 6) 


#integer\,\s+dimension\(6\)\s{0,}:: 
#\n|\,
# <- \n

dbl <- 0 %>% as.double()


# PLitime <- see [7]
#                            matrix(dbl,100,1) 


tsave  <-  persave <- rep(dbl,96)
# matrix(dbl,96,1) 


Pnsmall <-  cumPnsmall <-  Pismall <-  cumPismall <- 
  rep(dbl,2)
#matrix(dbl,2,1) 


probnumP <- cumprobnumP <- rep(dbl,4)
#matrix(dbl,4,1) 


#P <-                           #see [29]
cumP <- rep(dbl,5)
#matrix(dbl,5,1) 



#tot_sitesPLi <- prob_sitePLi <-  cumprob <- 
# matrix(dbl,6,1)  see [34]



# PLideg <-  #see [12]
#tau2 <- [19]
#                        matrix(dbl,300,1) 



#tau <-                           # see [28]
#asum <- 
#matrix(dbl,12000,1) see[32]**? bug?



#tautot <- rep(dbl,12300) # no longer used i think
#matrix(dbl,12300,1) 



cumsumaa  <-  asumtemp <- rep(dbl,sV)
#matrix(dbl,sV,1) 


#timevect <- percent_degrade <- [20] [21]
#percent_degrade2 <-  # no longer used also rn
#PLG <- [17]
# tvals <-                 # see [8]
#matrix(dbl,1000000,1) 
#c()#rep(dbl,1000000)


qstate <-
  matrix(dbl,12000, 9)  

#aa <- 
# matrix(dbl,12000, sV)  
aa <-  rep(dbl,sV)#matrix(dbl,12000, sV) 

#Prob_rebind <- [22]
#                           matrix(dbl,2, 1000000)  

dvlNA <- NA %>% as.integer()
dblNA <- NA %>% as.double()

#state <-                   # see [4]
#statetemp <- 
#matrix(dvl,Ntot,nodes^2)

#init_state <- #            #see[5], [6]
#  cumsuminit <- 
#                            matrix(dblNA,nodes,1)


Lat <- 
  matrix(dvl,nodes^2, nodes^2)


row2 <- col2 <-  rep(dvl,nodes^2 * Ntot)
# matrix(dvl, nodes^2 * Ntot ,1)


lysis_time <- tPA_time <- firstPLi <- 
  #matrix(dblNA, simulations,1)
  #rep(dblNA,simulations)
  c()


tPAunbind <- tPAPLiunbd <- ltPA <- Plasmin <- max_Plg <- 
  lysiscomplete <- countvect <- PLGbd <- PLGunbd <- c()
#matrix(dvl, simulations,1)
#rep(dvl,simulations)





kplgoff = kplgon*KdPLGintact 
kplgoffnick = kplgon*KdPLGnicked 
kaoff12 = ktPAon*KdtPAyesplg    
kaoff10 = ktPAon*KdtPAnoplg    


#print(format(Lat, digits = 3, justify = "right"))
# Lat is [
# 0  0  0  0   ...
# 0  0  0 ...  ...
# 0  0... ...  ...
# 0...  ....   ...
# ... 
# ]

for (i in 1:nodes) {
  for (j in 1:nodes) {
    Lat[(i - 1)*nodes + j, (i - 1)*nodes + j] = 1 
    if (i != 1) Lat[(i - 1)*nodes + j, (i - 2)*nodes + j] = 1 
    if (j != 1) Lat[(i - 1)*nodes + j, (i - 1)*nodes + j - 1] = 1 
    if (i != nodes) Lat[(i - 1)*nodes + j, (i)*nodes + j] = 1 
    if (j != nodes) Lat[(i - 1)*nodes + j, (i - 1)*nodes + j + 1] = 1 
  }
}

#print(format(Lat, digits = 3, justify = "right"))
# Lat is [
# 1  1  0  0   ...
# 1  1  1  0.  ...
# 0  1  1  1 ...
# 0  0  1  1   ...
# 0  0  0  1
# ] essentialy each 1 on main diagonal connects to 1s on each side
#          1
#    1     1_main   1
#          1

seeder <- sample(1:50000,1)#350060996 %>% round(,0)
set.seed(seed)

#need to find our which form of IEE rounding
# if (seed == 0) seed = mscw()
# 
# ! seed = mscw() !randomly generate seed
# ! seed=-34041038
# write (*, *), ' seed=', seed
# 
# stater(1) = 129281
# stater(2) = 362436069
# stater(3) = 123456789
# stater(4) = seed
# call set_kiss32(stater)
# 
# call get_kiss32(stater) # no idea what this stuff is

Tdoublets <-  as.integer(Ntot*nodes**2 )

# pi exists as is in R as just "pi"
# ones have been made a matrix of 1s already  see [1]
# V has become a matrix from V_tmp already    see [2]
# q exists already                            see [3]

#lysiscomplete = 0     #exists already
#firstPLi = 0          #exists already
countfp = 0           #exists already

# For building purposes right now
# setting simulations to 1000 to avoid
# computational burden when building

# THIS begins the giant loop
# for conveniance, mark start and end point
# will be denoted as [SimStartLoop]

c2 = kplgoff
c3 = kplgon*freeplg
c5 = kplioff
c6 = kaoff12
c8 = kplgon*freeplg
c8a = kplgoffnick #!*freeplg
c9 = kplgoff
c11 = kplioff
c12 = kaoff12
c13 = kaoff10
c14 = kplioff
c15 = kplioff
c20 = kapcat
c25 = kplioff

#probability for fully intact fibrin
prob_N02i = 1/(1 + (0.5*KdPLGintact/freeplg) + (0.5*freeplg/KdPLGintact))
prob_N00i = 0.5*KdPLGintact/(freeplg*(1 + (0.5*KdPLGintact/freeplg) + (0.5*freeplg/KdPLGintact)))
prob_N22i = 0.5*freeplg/(KdPLGintact*(1 + (0.5*KdPLGintact/freeplg) + (0.5*freeplg/KdPLGintact)))

#probability for nicked fibrin
prob_N02n = 1/(1 + (0.5*KdPLGnicked/freeplg) + (0.5*freeplg/KdPLGnicked))
prob_N00n = 0.5*KdPLGnicked/(freeplg*(1 + (0.5*KdPLGnicked/freeplg) + (0.5*freeplg/KdPLGnicked)))
prob_N22n = 0.5*freeplg/(KdPLGnicked*(1 + (0.5*KdPLGnicked/freeplg) + (0.5*freeplg/KdPLGnicked)))


#test
# this_is_a_test <- function(){
#   erm <<-prob_N22n
#   ok <- erm
#   return(ok)
# } oh i see, if not spec. inside local env of function, searches global


SIM_params <- data.frame(
  fiber_radius = radius
  ,diss_const_tPA_wPLG = KdtPAyesplg
  ,diss_const_tPA_woPLG = KdtPAnoplg
  ,diss_const_PLG_intact = KdPLGintact
  ,diss_const_PLG_nicked = KdPLGnicked
  ,bind_rate_tPA = ktPAon
  ,bind_rate_PLG = kplgon
  ,conc_free_PLG = freeplg
  ,deg_rate_fibrin = kdeg
  ,unbind_rate_PLi = kplioff
  ,activation_rate_PLG = kapcat
  ,exposure_rate_binding_site = kncat
  ,nodes_in_micro_row = nodes
  ,snap_proportion = snap_proportion
  ,seed = seeder
)

#ORIGINAL SET MUST MATCH SUM STATS OF sp_f23.h5
# testing if model is signif difference from fortan

# sum stats for fiber degraded, plifirst time, pli_generated_num
#sim_final_time #tpa_final_num #tpa_leaving_time
#tpa_unbound_by_pli #tpa_unbound_kinetic #

#full_df %>% filter(., file_name == "sp_f23.h5") %>% summary()
 # fiber_degraded  pli_first_time    pli_generated_num sim_final_time     tpa_final_num    
 # Mode :logical   Min.   :  0.000   Min.   :0.0000    Min.   :  0.0012   Min.   :0.00000  
 # FALSE:43731     1st Qu.:  0.000   1st Qu.:0.0000    1st Qu.:  8.5097   1st Qu.:0.00000  
 # TRUE :6269      Median :  0.000   Median :0.0000    Median : 20.0300   Median :0.00000  
 #                 Mean   :  3.146   Mean   :0.1632    Mean   : 27.2868   Mean   :0.01798  
 #                 3rd Qu.:  0.000   3rd Qu.:0.0000    3rd Qu.: 37.8169   3rd Qu.:0.00000  
 #                 Max.   :235.517   Max.   :5.0000    Max.   :357.4574   Max.   :1.00000  
 
# tpa_leaving_time   tpa_unbound_by_pli tpa_unbound_kinetic  file_name        
 # Min.   :  0.0012   Mode :logical      Mode :logical       Length:50000      
 # 1st Qu.:  7.8655   FALSE:45731        FALSE:5168          Class :character  
 # Median : 18.4863   TRUE :4269         TRUE :44832         Mode  :character  
 # Mean   : 26.1304                                                            
 # 3rd Qu.: 36.2674                                                            
 # Max.   :357.4574  

# DO MIND THAT pli_first_time is not inherently filtered for fiber degraded == T
# here is the associated statistic for pli_firsttime WHEN Fiber degraded == T



# COMPARISON : fiber_degraded == T :: 
#
#                              FT: 0.12538
#                               R: 0.13

# COMPARISON : tpa_unbound_by_pli == T ::
#
#                              FT: 0.08538
#                               R: 0.083

# COMPARISON : tpa_unbound_kinetic == T ::
#
#                              FT: 0.89664
#                               R: 0.904

# COMPARISON : Pli_first_time WHEN fiber_degraded == T
#    min 	         q1 	     median   	mean 	    q3 	    max 	
#    0.0004892 	7.138617 	17.20760 	25.08966 	34.68722 	235.5169

#
#  summary(This_df)  ## OUTPUT SUM STATS OF THIS SIM ON ONLY 1000 SIMULATIONS
 # fiber_degraded pli_first_time     pli_generated_num sim_final_time      tpa_final_num  
 # 1   :130       Min.   :  0.2781   Min.   :0.000     Min.   :  0.04483   Min.   :0.000  
 # NA's:870       1st Qu.:  9.1780   1st Qu.:0.000     1st Qu.:  8.46621   1st Qu.:0.000  
 #                Median : 20.0403   Median :0.000     Median : 21.17805   Median :0.000  
 #                Mean   : 28.0465   Mean   :0.184     Mean   : 29.80391   Mean   :0.012  
 #                3rd Qu.: 37.1713   3rd Qu.:0.000     3rd Qu.: 45.03549   3rd Qu.:0.000  
 #                Max.   :137.1623   Max.   :4.000     Max.   :200.42667   Max.   :1.000  
 #                NA's   :870                                                             
 # tpa_leaving_time    tpa_unbound_by_pli tpa_unbound_kinetic
 # Min.   :  0.04483   0:917              0: 96              
 # 1st Qu.:  7.76653   1: 83              1:904              
 # Median : 18.46039                                         
 # Mean   : 26.99666                                         
 # 3rd Qu.: 37.64606                                         
 # Max.   :200.42667  



# Results seem similar from a subjective lense. I should run a statistical test
#         to confirm but I'm overall pretty confident these are the same


for (stats in 1:simulations) {        # [SimStartLoop]
  
  # Nsave = 0
  # Ninteger = 0
  # tsave = 0                     # to save 0 to 96%
  # idk if I mind this
  # so commenting all
  
  # [4] it looks like @ each loop, the matrix state
  #     must reset to the following each simulation.
  #     0 stands for false (not exposed) 1 stands for
  #     exposed. state is a matrix
  #     of the following dimensions: nodes^2 x Ntot
  
  state <- matrix(dvl,nodes^2,Ntot)
  statetemp <- matrix(dvl,nodes^2,Ntot)
  state[,1:Nplginit] <- 1
  # > state
  #       [,1] [,2] [,3] [,4] [,5] [,6]
  # [1,]    1    0    0    0    0    0
  # [2,]    1    0    0    0    0    0
  
  
  # resetting variables and vecotors here
  # potential bug previously?
  # tPA_loc was set to 0 after setting tPA_loc(1) to value
  # in the forloop called [find tpa start]
  noww <- c()
  countt = 1
  countPLi = 0
  PLGbind = 0
  PLGunbind = 0
  plotnumber = 1
  PLitime <- rep(dbl,100)             # [7]
  reaction2 = 0
  reaction3 = 0
  tvals <- c()#rep(dbl,1000000)           # [8]
  reaction <- c()#rep(dvl,est_mts)        # [9]
  Anumber<- c()#rep(dvl,est_mts)          # [10]
  PLi <- c()#rep(dvl,est_mts)             # [11]
  #matrix(dvl,est_mts,1)
  PLideg <- rep(0,300)                # [12]
  PLilocation <- matrix(dvl,300, 2)   #[13]
  PLiwaste <- rep(dvl,300)            # [14]
  # matrix(dvl,est_mts,1)
  tPA <- tPA_loc<- c()#(dvl,est_mts)   # [15] [16]
  # matrix(dbl,1000000,1) 
  PLG <- c()#rep(dbl,1000000)             #[17]
  per_lysis <- c()#rep(dvl,est_mts)       #[18]
  tau2 <- c()#rep(0,300)                  #[19] matrix(dbl,300,1)
  
  timevect <- percent_degrade <-c()#matrix(dbl,est_mts,1) 
  # [20] [21]
  
  Prob_rebind <- matrix(dbl,2, 1000000)  # [22]
  
  degrade_tot <- undegrade_tot <- undegrade_new <- 
    PLi_loc <- 
    c()#matrix(dvl,est_mts,1)
  #[23]
  undegrade_tot[1] = Nplginit*nodes^2
  undegrade_new[1] = 0
  Anumber[1] = 1
  
  temp <- 0
  PLi[1] = 0         
  tPA[1] = 1 
  
  countstate = 0
  countstate2 = 0
  
  rin = runif(1)    # random number to pick 1st tPA loc
  # note, the function cumsum now
  # has been integrated into base r
  # s.t. the following retrieves the
  # same effect desired as in fortran 
  
  # [5]   in state is a vector. its length is one row of lattice
  #       and is scalled s.t. sum(init_state) = 1
  # [6]   cumsuminit can now just diectly be made like so
  
  cumsuminit <- rep(1.0/nodes,nodes) %>% 
    cumsum(.) 
  # cumsuminit = 
  # [1] 0.1428571 0.2857143 0.4285714 0.5714286
  # [5] 0.7142857 0.8571429 1.0000000
  
  
  for (j in 1:nodes) {      # [find tpa start] see [6b]
    if (cumsuminit[j] >= rin) {
      break
    }
  }                         # end [find tpa start]
  
  # [6] init_entry = j # can just direct into tpa Loc
  # [6b]  first node in vector cumsuminit more than
  #       the random number rin then becomes tpa start
  init_entry <- j
  tPA_loc[1] <- init_entry
  
  r = runif(1)
  
  if (r<= prob_N02i) {      # [N02 oe N22] see [6b]
    state[init_entry, 1] <-  2
  } else {
    state[init_entry, 1] <-  3
  }
  
  t = 0
  tfinal = 24*60*60   
  
  countstate = 0
  countstate2 = 0
  
  for (i in 1:nodes^2) {
    for (j in 1:Ntot) {
      if (state[i,j] == 2 | state[i,j] == 6) {
        countstate = countstate + 1
      }
      
      if (state[i,j] == 1) {
        countstate2 = countstate2 + 1
      }
    }
  }
  
  PLG[1] = countstate + (prob_N02i + 2*prob_N22i)*countstate2 
  degrade_tot[1] = 0
  percent_degrade[1] = 0 
  
  D = 10**(7)             
  
  vol3 = 0.001*(radius)**2*pi   
  
  bs = Nplginit*nodes**2/vol3  
  
  bs = bs/602.2 
  
  repeat{
    
    countstate3 = 0
    row <- col <-  rep(dvl,300) #matrix(dvl,300,1)
    
    for (j in 1:Ntot) {
      for (i in 1:nodes^2) {
        if (state[i,j]>1) {
          countstate3 = countstate3 + 1   # if the state of the doublet is 
          # larger than 1 (i.e., any state other than Nplg),
          # countt it as an active site
          row[countstate3] = i
          col[countstate3] = j
        }
      }
    }
    
    A = countstate3 #number of activated states
    
    if (A == 0) {
      break
    }
    
    
    tau <- c()#rep(dbl,12000)   #matrix(dbl,12000,1)   # [28]     
    #reset tau to zero each time
    newsite = 0
    newsite2 = 0
    newdoublet = 0
    newdoublet2 = 0
    P <- rep(dbl,5)                   #matrix(dbl,5,1) [29]
    currentPLi = 0
    #<- rep(dvl,2)          #matrix(dvl,2,1) [30]
    countcurPLi = 0
    
    for (i in 1:A) {
      # qstate is  matrix(dbl,12000, 9) 
      # q is diag(9) 
      qstate[i,] <- q[state[row[i], col[i]], ] 
      #!qstate(i,) saves vector of row corresponding to the state of 
      # the active site, and 0's elsewhere
      # q[#,] would return whole row... since q is an identity matrix
      # doesnt seem to matter what is what
      qconvert = 0
      numN = 0
      numP2 = 0
      numP1 = 0
      numP1b = 0
      numPt = 0 #
      numP6 = 0
      n12 = 0
      n10 = 0
      n13 = 0
      n03 = 0
      n23 = 0
      n33 = 0
      nplgi = 0 #
      nplgn2 = 0 #
      
      numN <- findint(state[row[i],], 0)
      #####
      if (state[row[i],col[i]] >= 2 & state[row[i],col[i]]<=4) {
        numP2 <- findint(state[row[i],], 2)
        numP1 <- findint(state[row[i],c(2:Ntot)],1)
        numP1b <- findint(state[row[i],1],1)
        numP6 <- findint(state[row[i], ],6)
        
        numPt <-  (2*prob_N22n + prob_N02n)*numP1 + 
          numP2 + numP6 + 
          (2*prob_N22i + prob_N02i)*numP1b
        
        if (numPt>0) {
          qconvert <- 1
        }
      }
      
      if ( (state[row[i],col[i]] >= 4 & state[row[i],col[i]]<=7) |  
           state[row[i],col[i]] == 9) {
        n12 <- findint(state[row[i],],2)
        n10 <- findint(state[row[i],],3)
        n13 <- findint(state[row[i],],4)
        n03 <- findint(state[row[i],],5)
        n23 <- findint(state[row[i],],6)
        n33 <- findint(state[row[i],],7)
        nplgi <- findint(state[row[i],1],1)
        nplgn2 <- findint(state[row[i],c(2:Ntot)],1)
      }
      
      
      
      
      ## ?? we put into and then erase it right after? 
      ## we can thus save space by just making aa and asum a list??
      
      # if (col[i] == 1 ) {
      #   aa[i,1] = c2*qstate[2,i]  
      #   aa[i,6] = c2*qstate[6,i]
      # }else{
      #   if (col[i]>1) {
      #     aa[i,1] = c8a*qstate[2,i]  
      #     aa[i,6] = c8a*qstate[6,i]
      #   }
      # }
      
      # aa[i, 2] = c3*qstate[ 3,i]     
      # aa[i, 5] = c3*qstate[ 5,i]     
      # aa[i, 3] = c5*qstate[ 4,i]   
      # aa[i, 4] = c6*qstate[ 4,i]   
      # aa[i, 7] = c11*qstate[ 7,i]   
      # aa[i, 8] = c12*qstate[ 2,i]   
      # aa[i, 9] = c13*qstate[ 3,i]   
      # aa[i, 10] = c14*qstate[ 5,i]   
      # aa[i, 11] = c15*qstate[ 6,i]   
      # 
      # 
      # aa[i, 12] = c20*numPt*qconvert   
      # aa[i, 13] = c25*qstate[ 9,i]   
      # aa[i, 14] = kdeg*n13    
      # aa[i, 15] = kdeg*n03    
      # aa[i, 16] = kdeg*n23    
      # aa[i, 17] = kdeg*n33    
      # aa[i, 18] = kncat*numN*qstate[ 4,i] 
      # aa[i, 19] = kncat*numN*qstate[ 5,i] 
      # aa[i, 20] = kncat*numN*qstate[ 6,i] 
      # aa[i, 21] = kncat*numN*qstate[ 7,i] 
      # aa[i, 22] = kncat*numN*qstate[ 9,i] 
      # aa[i, 23] = kdeg*n10    
      # aa[i, 24] = kdeg*n12    
      # aa[i, 25] = kdeg*nplgi  
      # aa[i, 26] = kdeg*nplgn2 
      # 
      # #asum <- 
      #   #matrix(dbl,12000,1) 
      # asum <-  sum(aa[i, ])
      
      aa <-  rep(dbl,sV)
      asum <- 0
      
      if (col[i] == 1 ) {
        aa[1] = c2*qstate[i,2]  
        aa[6] = c2*qstate[i,6]
      } else if (col[i]>1) {
        aa[1] = c8a*qstate[i,2]  
        aa[6] = c8a*qstate[i,6]
      }
      
      
      aa[2] = c3*qstate[i, 3]     
      aa[5] = c3*qstate[i, 5]     
      aa[3] = c5*qstate[i, 4]   
      aa[4] = c6*qstate[i, 4]   
      aa[7] = c11*qstate[i, 7]   
      aa[8] = c12*qstate[i, 2]   
      aa[9] = c13*qstate[i, 3]   
      aa[10] = c14*qstate[i, 5]   
      aa[11] = c15*qstate[i, 6]   
      
      
      aa[12] = c20*numPt*qconvert   
      aa[13] = c25*qstate[i, 9]   
      aa[14] = kdeg*n13    
      aa[15] = kdeg*n03    
      aa[16] = kdeg*n23    
      aa[17] = kdeg*n33    
      aa[18] = kncat*numN*qstate[i, 4] 
      aa[19] = kncat*numN*qstate[i, 5] 
      aa[20] = kncat*numN*qstate[i, 6] 
      aa[21] = kncat*numN*qstate[i, 7] 
      aa[22] = kncat*numN*qstate[i, 9] 
      aa[23] = kdeg*n10    
      aa[24] = kdeg*n12    
      aa[25] = kdeg*nplgi  
      aa[26] = kdeg*nplgn2 
      
      #asum <- 
      #matrix(dbl,12000,1) 
      asum <-  sum(aa)
      asumtemp <- aa/asum
      cumsumaa <- 0
      cumsumaa <- cumsum(asumtemp) # vector of 26
      aa <-  rep(dbl,sV)
      
      r = runif(1)
      rin = runif(1)
      
      #tau <- rep(dbl,12000)
      tau[i] <- log(1/r)/asum # a single digit
      # saves this for all binding sites
      asum <- 0
      
      for (j in 1:sV) {
        if (cumsumaa[j] >= rin) {
          bigJ[i] <-  j # bigJ is vector of 26
          # saves location where cumsumaa >= rand num
          break
        }
        
      }
      
    }
    
    #minvect
    
    place <- minvect(tau) # a single integer, minimum LOCATION
    # binding location with minimum
    entry <- tau[place]
    
    #qstate <-matrix(dbl,12000, 9) refrenece comment, disregard
    # its first row has saved the identity matrix of q corresponding
    for (i in 1:A) {
      if (i != place) {
        qstate[i,] <- rep(0,9)
        bigJ[i] <- 0
      } #now bigJ has only the vector of 26 with 
    }
    #tot_sitesPLi <- prob_sitePLi <-  cumprob <- 
    # matrix(dbl,6,1)  see [34]
    
    # sites_N03 <- 
    # sites_N10 <- 
    # sites_Nplg <- 
    # sites_Nplgn <-  #see [35]
    # sites_Nplgi <- 
    # sites_degrade <- 
    # matrix(dvl,6,1) 
    
    tot_sitesPLi <- c()#rep(dbl,6)
    prob_sitePLi <- c()#rep(dbl,6)
    
    sites_N03 <- sites_N10 <- sites_Nplg <- 
      sites_Nplgn <- sites_Nplgi <- sites_degrade <- 
      c()#rep(dvl,6)
    
    alpha = 1   # integer
    
    # BIGPLACE IFS , goes until END BIGPLACE IFS
    
    if (bigJ[place] == 14) { ####################
      ############################################# 14
      
      findit <- findfirst(state[row[place], ],4)
      if (findit == 0) {
        stop("Problem, should have a N13, bigJ[place] == 14 
                and no 4 on state[row[place],]")
      }
      state[row[place], findit] = -1 #degraded
      
      for (iPLi in 1:countPLi) {
        if (PLilocation[iPLi,1]== row[place] & 
            PLilocation[iPLi,2]== findit ) {
          countcurPLi = countcurPLi + 1
          currentPLib[countcurPLi] = iPLi
        }
      }
      
      if (countcurPLi == 2) {
        stop("Problem: should only have one plasmin, 
                countcurPLi == 2.")
      } else if (countcurPLi == 1){
        loc = currentPLib[1]
      }else if (countcurPLi == 0) {
        stop("Problem: should be 1 plasmin. 
                countcurPLi == 0")
      }else{
        stop("SHOULDNT HAPPEN, review here [1]")
      }
      
      countt = countt + 1
      t = t + entry
      
      tvals[countt] = t
      Anumber[countt] = A
      reaction[countt] = bigJ[place]
      connect_ind <- c()#rep(dvl,6)
      countlat = 0
      
      for (j in 1:nodes^2) {
        if (Lat[PLilocation[loc,1],j] > 0){
          countlat = countlat + 1
          connect_ind[countlat] = j
        }
      }
      
      sizenonzero_connectind = countlat
      
      p_rebind <- movetpa(state)
      
      reaction3 = reaction3 + 1
      
      movepli(state,loc,statetemp)
      
      PLilocation <- PLilocationt
      state <-  statetemp #oh.... then cant we just save this to the new state?
      
      
    } else if (bigJ[place] == 15) {####################
      ############################################# 15
      
      findit <- findfirst(state[row[place], ],5)
      if (findit==0) {
        stop("Problem, should have a N03, bigJ[place] == 14 
                and no 5 on state[row[place],]")
      }
      
      state[row[place], findit] = -1 #degraded
      
      for (iPLi in 1:countPLi) {
        if (PLilocation[iPLi,1]== row[place] & 
            PLilocation[iPLi,2]== findit ) {
          countcurPLi = countcurPLi + 1
          currentPLib[countcurPLi] = iPLi
        }
      }
      
      if (countcurPLi == 2) {
        stop("Problem: should only have one plasmin, 
                countcurPLi == 2.")
      } else if (countcurPLi == 1){
        loc = currentPLib[1]
      }else if (countcurPLi == 0) {
        stop("Problem: should be 1 plasmin. 
                countcurPLi == 0")
      }else{
        stop("SHOULDNT HAPPEN, review here [2]")
      }
      
      countt = countt + 1
      t = t + entry
      
      tvals[countt] = t
      Anumber[countt] = A
      reaction[countt] = bigJ[place]
      connect_ind <- c()#rep(dvl,6)
      countlat = 0
      
      for (j in 1:nodes^2) {
        if (Lat[PLilocation[loc,1],j] > 0){
          countlat = countlat + 1
          connect_ind[countlat] = j
        }
      }
      
      sizenonzero_connectind = countlat
      
      movepli(state,loc,statetemp)
      
      PLilocation <- PLilocationt
      state <-  statetemp
      
    } else if (bigJ[place] == 16) {####################
      ############################################# 16
      
      findit <- findfirst(state[row[place], ],6)
      if (findit==0) {
        stop("Problem, should have a N23, bigJ[place] == 14 
                and no 6 on state[row[place],]")
      }
      
      state[row[place], findit] = -1 #degraded
      
      for (iPLi in 1:countPLi) {
        if (PLilocation[iPLi,1]== row[place] & 
            PLilocation[iPLi,2]== findit ) {
          countcurPLi = countcurPLi + 1
          currentPLib[countcurPLi] = iPLi
        }
      }
      
      if (countcurPLi == 2) {
        stop("Problem: should only have one plasmin, 
                countcurPLi == 2.")
      } else if (countcurPLi == 1){
        loc = currentPLib[1]
      }else if (countcurPLi == 0) {
        stop("Problem: should be 1 plasmin. 
                countcurPLi == 0")
      }else{
        stop("SHOULDNT HAPPEN, review here [3]")
      }
      
      countt = countt + 1
      t = t + entry
      
      tvals[countt] = t
      Anumber[countt] = A
      reaction[countt] = bigJ[place]
      connect_ind <- c()#rep(dvl,6)
      countlat = 0
      
      for (j in 1:nodes^2) {
        if (Lat[PLilocation[loc,1],j] > 0){
          countlat = countlat + 1
          connect_ind[countlat] = j
        }
      }
      
      sizenonzero_connectind = countlat
      
      movepli(state,loc,statetemp)
      
      PLilocation <- PLilocationt
      state <-  statetemp
      
      
    } else if (bigJ[place] == 17) {####################
      ############################################# 17
      
      findit <- findfirst(state[row[place], ],7)
      if (findit==0) {
        stop("Problem, should have a N33, bigJ[place] == 14 
                and no 7 on state[row[place],]")
      }
      
      state[row[place], findit] = -1 #degraded
      
      for (iPLi in 1:countPLi) {
        if (PLilocation[iPLi,1]== row[place] & 
            PLilocation[iPLi,2]== findit ) {
          countcurPLi = countcurPLi + 1
          currentPLib[countcurPLi] = iPLi
        }
      }
      
      if (countcurPLi == 2) {
        
        # warning("Problem: should only have one plasmin, 
        #         countcurPLi == 2.")
        loc <- currentPLib[1]
        loc2 <- currentPLib[2]
        
      } else if (countcurPLi == 1){
        
        # loc = currentPLib[1]
        stop("Problem: should only have two plasmin, 
                 countcurPLi == 1.")
        
      }else if (countcurPLi == 0) {
        
        stop("Problem: should be 2 plasmin. 
                countcurPLi == 0")
        
      }else{
        stop("SHOULDNT HAPPEN, review here [4]")
      }
      
      countt = countt + 1
      t = t + entry
      
      tvals[countt] = t
      Anumber[countt] = A
      reaction[countt] = bigJ[place]
      connect_ind <- c()#rep(dvl,6)
      countlat = 0
      
      for (j in 1:nodes^2) {
        if (Lat[PLilocation[loc,1],j] > 0){
          countlat = countlat + 1
          connect_ind[countlat] = j
        }
      }
      
      sizenonzero_connectind = countlat
      
      PLilocsmall[1] = PLilocation[loc,1]
      PLilocsmall[2] = PLilocation[loc,2]
      
      for (iPLi2 in 1:countPLi) {
        if(PLilocation[iPLi2,1] == PLilocsmall[1] &
           PLilocation[iPLi2,2] == PLilocsmall[2] & 
           iPLi2 != loc ){
          loc2b <- iPLi2
          
          if (loc2b != loc2) {
            stop(paste("Problem: loc2 and loc2b not the same. loc2 = ", loc2,
                       " ; and loc = ", loc))
          }
          
        }
      }
      
      move2pli(state, loc, loc2, statetemp)
      PLilocation <- PLilocationt
      state <- statetemp
      
    } else if (bigJ[place] == 23) {####################
      ############################################# 23
      
      findit <- findfirst(state[row[place], ],3)
      if (findit==0) {
        stop("Problem, should have a N10, bigJ[place] == 14 
                and no 3 on state[row[place],]")
      }
      state[row[place], findit] = -1
      countt = countt + 1
      t = t + entry
      
      tvals[countt] = t
      Anumber[countt] = A
      reaction[countt] = bigJ[place]
      
      p_rebind <- movetpa(state)
      reaction3 = reaction3 + 1
      
    } else if (bigJ[place] == 24) {####################
      ############################################# 24
      
      findit <- findfirst(state[row[place], ],2)
      if (findit==0) {
        stop("Problem, should have a N10, bigJ[place] == 14 
                and no 2 on state[row[place],]")
      }
      state[row[place], findit] = -1
      countt = countt + 1
      t = t + entry
      
      tvals[countt] = t
      Anumber[countt] = A
      reaction[countt] = bigJ[place]
      
      p_rebind <- movetpa(state)
      eaction3 = reaction3 + 1
      
      
    } else if (bigJ[place] == 25) {####################
      ############################################# 25
      
      if(state[row[place],1] == 1 ){
        state[row[place],1] <- -1
      } else{
        stop("Problem: should have an intact Nplg.")
      }
      
      countt = countt + 1
      t = t + entry
      
      tvals[countt] = t
      Anumber[countt] = A
      reaction[countt] = bigJ[place]
      
    } else if (bigJ[place] == 26) {####################
      ############################################# 26
      
      findit <- findfirst(state[row[place], -1],1)
      if (findit==0) {
        stop("Problem, should have a nicked Nplg, bigJ[place] == 26 
                and no 1 on state[row[place],]")
      } # stop / save error instance
      findit = findit + 1 
      state[row[place], findit] = -1
      countt = countt + 1
      t = t + entry
      
      tvals[countt] = t
      Anumber[countt] = A
      reaction[countt] = bigJ[place]
      
      
    } else if (
      bigJ[place] == 18 | bigJ[place] == 19 | bigJ[place] == 20 |
      bigJ[place] == 21 | bigJ[place] == 22
    ) { ####################
      ############################################# 18 19 20
      ############################################# 21 22
      
      newstate <- qstate[place, ] + V[bigJ[place], ]
      # newstate is a list of nine
      # qstate is  matrix(dbl,12000, 9) 
      # thus qstate[place, ] is vector of 9
      #
      qstate[place, ] <- 0 # oh this works
      
      for (j in 1:9) {
        if (newstate[j] > 0) {
          statetype = j
          break
        }
      }
      
      if (statetype == state[row[place], col[place]]){
        
        NtoNplg <- findfirst(state[row[place], ],0)
        state[row[place], NtoNplg] <- 1
        
        if (NtoNplg == 0) {
          stop("Problem: there must be Ns to expose if this rxn occurs")
        }
        
      }else{
        stop("Problem:this should be exposure;
                  so statetype should equal state(col(place),row(place)")
      }
      
      countt = countt + 1
      t = t + entry 
      
      tvals[countt] = t
      Anumber[countt] = A
      reaction[countt] = bigJ[place]
      
    } else if (bigJ[place] == 3 | bigJ[place] == 7 | bigJ[place] == 10 |
               bigJ[place] == 11 | bigJ[place] == 13) {####################
      ############################################# 3  7  10
      ############################################# 11 13
      
      newstate <- qstate[place, ] + V[bigJ[place], ]
      
      qstate[place, ] <- 0 # oh this works
      
      for (j in 1:9) {
        if (newstate[j] > 0) {
          statetype = j
          break
        }
      }
      #statetype <- findfirst(newstate,0)
      
      if (statetype == 8) {
        statetype <- -1
      }
      
      state[row[place], col[place]] <- statetype
      
      countt = countt + 1
      t = t + entry 
      
      tvals[countt] = t
      Anumber[countt] = A
      reaction[countt] = bigJ[place]
      
      for (iPLi in 1:countPLi) {
        if (PLilocation[iPLi,1]== row[place] & 
            PLilocation[iPLi,2]== col[place] ) {
          countcurPLi = countcurPLi + 1
          currentPLib[countcurPLi] = iPLi
        }
      }
      
      connect_ind <- c()#rep(dvl,6)
      countlat3 = 0
      
      for (j in 1:nodes^2) {
        if (Lat[row[place],j] > 0){
          countlat3 = countlat3 + 1
          connect_ind[countlat3] = j
        }
      }
      
      sizenonzero_connectind = countlat3
      
      tot_sitesPLi <- c()#rep(dbl,6)
      prob_sitePLi <- c()#rep(dbl,6)
      
      sites_N03 <- 
        sites_N10 <-
        sites_Nplg <- 
        sites_Nplgn <- 
        sites_Nplgi <- 
        sites_degrade <- 
        c()#rep(dvl,6)
      alpha = 1 
      
      if (countcurPLi == 2) {
        
        randPLi <- runif(1)
        if (randPLi <= 0.5) {
          currentPLi = currentPLib[1]
        } else {
          currentPLi = currentPLib[2]
        }
        
      } else if (countcurPLi == 1) {
        currentPLi = currentPLib[1]
        
        if(bigJ[place] == 7){
          stop("Problem: there should be 2 PLis if rxn 7 occurs")
        }
        
      }
      
      movepli(state,currentPLi,statetemp)
      
      PLilocation <-  PLilocationt
      state <-  statetemp
      
    } else if (bigJ[place] == 12) {####################
      ############################################# 12
      
      countPLi = countPLi + 1 
      PLilocation[countPLi, 1] = row[place] 
      
      numP2 <- findint(state[row[place],], 2)
      
      numP1 <-  findint(state[row[place], -1 ],1)
      numP1b <- findint(state[row[place],  1 ],1)
      
      numP6 <- findint(state[row[place], ],6)
      
      numPt <-  (2*prob_N22n + prob_N02n)*numP1 + 
        numP2 + numP6 + 
        (2*prob_N22i + prob_N02i)*numP1b
      
      probnumP[1] <- (2*prob_N22n + prob_N02n)*as.double(numP1)/numPt
      probnumP[2] <- as.double(numP2)/numPt 
      probnumP[3] <- as.double(numP6)/numPt 
      probnumP[4] <- (2*prob_N22i + prob_N02i)*as.double(numP1b)/numPt
      
      probnumP_temp <- probnumP/sum(probnumP)
      colPLGb <- sample(1:length(probnumP_temp),1,prob = probnumP_temp, replace = T)
      
      if (colPLGb == 3) {
        colPLGfind <- 6
        colPLG <- randchoice(state[row[place],],colPLGfind)
      } else if (colPLGb == 4) {
        colPLG = 1 
      } else if (colPLGb == 2) {
        colPLGfind = 2
        colPLG <- randchoice(state[row[place],],colPLGfind)
      } else if (colPLGb == 1) {
        colPLGfind = 1
        colPLG <- randchoice(state[row[place],-1],colPLGfind)
        colPLG <<-  colPLG + 1
        # somehow this is able to return zero as a value
        # I see, i made an error when making the location vector in
        # the function randchoice
        
      }
      
      PLilocation[countPLi, 2] = colPLG          
      PLitime[countPLi] = t
      
      if (countPLi == 1) {
        countfp = countfp + 1
        firstPLi[stats] = t + entry
      }
      
      qstate[place, ] = 0
      
      if (state[row[place], colPLG] == 1) {
        rP <- runif(1)
        if (colPLG == 1) {
          
          if (rP <= prob_N02i/(prob_N02i + prob_N22i)) {
            statetype = 5
          }else{
            statetype = 6
          }
          
        } else {
          
          if (rP <= prob_N02n/(prob_N02n + prob_N22n)) {
            statetype = 5
          } else{
            statetype = 6
          }
          
        }
      } else if (state[row[place], colPLG] == 2) {
        
        statetype = 4 
        
      } else if (state[row[place], colPLG] == 6) {
        
        statetype = 7 
        
      } else{
        
        stop("Problem: should be a doublet with PLG")
      
        }
      
      
      
      state[row[place], colPLG] = statetype
      
      countt = countt + 1
      t = t + entry 
      
      tvals[countt] = t
      Anumber[countt] = A
      reaction[countt] = bigJ[place]
      
      
    } else if (bigJ[place] == 4 | bigJ[place] == 8 | bigJ[place] == 9) {####################
      ############################################# 4  8  9
      
      reaction2 = reaction2 + 1
      newstate <- qstate[place, ] + V[bigJ[place], ]
      
      qstate[place, ] <- 0 # oh this works
      
      for (j in 1:9) {
        if (newstate[j] > 0) {
          statetype = j
          break
        }
      }
      
      if (statetype == 8) {
        statetype <- -1
        stop("Problem: degradation should not occur here")
      }
      
      state[row[place], col[place]] <- statetype
      
      countt = countt + 1
      t = t + entry 
      
      tvals[countt] = t
      Anumber[countt] = A
      reaction[countt] = bigJ[place]
      
      connect_ind <- c()#rep(dvl,6)
      countlat3 = 0
      
      for (j in 1:nodes^2) {
        if (Lat[row[place],j] > 0){
          countlat3 = countlat3 + 1
          connect_ind[countlat3] = j
        }
      }
      
      sizenonzero_connectind = countlat3
      
      tot_sitesPLi <- c()#rep(dbl,6)
      prob_sitePLi <- c()#rep(dbl,6)
      
      sites_N03 <- 
        sites_N10 <-
        sites_Nplg <- 
        sites_Nplgn <- 
        sites_Nplgi <- 
        sites_degrade <- 
        c()#rep(dvl,6)
      
      p_rebind <- movetpa(state)
      
    } else if (bigJ[place] == 2 | bigJ[place] == 5 | 
               bigJ[place] == 6 | bigJ[place] == 1) {####################
      #############################################  1  2  5  6
      
      newstate <- qstate[place, ] + V[bigJ[place], ]
      
      qstate[place, ] <- 0 # oh this works
      
      for (j in 1:9) {
        if (newstate[j] > 0) {
          statetype = j
          break
        }
      }
      
      if (statetype == 8) {
        statetype <- -1
        stop("Problem: degradation should not occur here")
      }
      
      state[row[place], col[place]] <- statetype
      
      countt = countt + 1
      t = t + entry 
      
      tvals[countt] = t
      Anumber[countt] = A
      reaction[countt] = bigJ[place]
      
    }
    
    # END BIGPLACE IFS
    noww[countt] <- statetype
    bigJ[place] <- 0
    
    tPA[countt] <- 0 
    #if(countt>1) {
    if (tPA[countt - 1] == 1) {
      countstate9 = 0
      for (i in 1:Ntot) {
        for (j in 1:nodes^2) {
          if (state[j,i] >= 2 & 
              state[j,i] <= 4) {
            countstate9 = countstate9 + 1
            rowtPA[countstate9] = j
          }
        }
      }
      
      tPA[countt] = countstate9
      if (tPA[countt]>1) {
        stop("should only have 1 tPA in system. something is wrong")
      }
      if (tPA[countt] == 1) {
        tPA_loc[countt] <-  rowtPA[countstate9]
    }
    }
    #}
    
    countstate13 = 0
    for (i in 1:Ntot) {
      for (j in 1:nodes^2) {
        if (state[j,i] == 2 | 
            state[j,i] == 6) {
          countstate13 = countstate13 + 1
        }
      }
    }
    
    countstate14 = 0
    for (j in 1:nodes^2) {
      if (state[j,1]==1) {
        countstate14 = countstate14 + 1
      }
    }
    
    countstate15 = 0
    for (i in 2:Ntot) {
      for (j in 1:nodes^2) {
        if (state[j,i] == 1 ) {
          countstate15 = countstate15 + 1
        }
      }
    }
    
    PLG[countt] = countstate13 + (prob_N02i + 2*prob_N22i)*countstate14 + (prob_N02n + 2*prob_N22n)*countstate15
    
    countstate10 = 0
    for (i in 1:Ntot) {
      for (j in 1:nodes^2) {
        if (state[j,i]  >= 4 &
            state[j,i]  <= 6) {
          countstate10 = countstate10 + 1
          rowPLi[countstate10] = j
          colPLi[countstate10] = i
        }
      }
    }
    
    countstate10a = 0
    for (i in 1:Ntot) {
      for (j in 1:nodes^2) {
        if (state[j,i]  == 9 ) {
          countstate10a = countstate10a + 1
          rowPLi[countstate10a] = j
          colPLi[countstate10a] = i
        }
      }
    }
    
    countstate11 = 0
    for (i in 1:Ntot) {
      for (j in 1:nodes^2) {
        if (state[j,i]  == 7 ) {
          countstate11 = countstate11 + 1
          rowPLi[countstate11] = j
          colPLi[countstate11] = i
        }
      }
    }
    
    PLi[countt] = countstate10 + 2*countstate11 + countstate10a
    
    if (countt>1) {
      
      if (PLi[countt] < PLi[countt-1]) {
        stop("Problem: should not lose any PLi from system. 
                something is wrong.")
        print(countt, PLi[countt], reaction[countt], PLi[countt - 1], reaction[countt-1])
      }
    }
    
    timevect[countt]= t
    
    countstate12 = 0
    for (i in 1:Ntot) {
      for (j in 1:nodes^2) {
        if (state[j,i]  > -1 &
            state[j,i]  <  9 ) {
          countstate12 = countstate12 + 1 
          rowPLi[countstate12] = j
          colPLi[countstate12] = i
        }
      }
    }
    
    Tdoublets2 = 6.0*nodes**2
    percent_degrade[countt] = 1 - countstate12/Tdoublets2
    
    if (countt > 50000) {
      
      stop("single simulation gone over 50k instances")
    }
    
    if (percent_degrade[countt] > snap_proportion) {
      lysiscomplete[stats] = 1
      break
    }
    
    
    if (t > tfinal) {
      break
    }
    
  } # end Gillespie loop (the repreat)
  
  
  lysis_time[stats] <- tvals[countt]
  if (tPA[countt] == 1) {
    tPA_time[stats] = tvals[countt]
  } else {
    loca <- findfirst(tPA,0)
    tPA_time[stats] = tvals[loca]
  }
  
  Plasmin[stats] <- PLi[countt]
  max_Plg[stats] <- max(PLG)
  countvect[stats] <- countt
  tPAunbind[stats] <- reaction2
  tPAPLiunbd[stats] <- reaction3
  ltPA[stats] <- tPA[countt]
  
}                                 # end [SimStartLoop]

if (length(lysiscomplete)<stats) {
  firstPLi[stats] <- NA
  lysiscomplete[stats] <- NA
}



This_df <- data.frame(
  fiber_degraded = as.factor(lysiscomplete)
  ,pli_first_time = firstPLi
  ,pli_generated_num = Plasmin
  ,sim_final_time = lysis_time
  ,tpa_final_num = ltPA
  ,tpa_leaving_time = tPA_time
  ,tpa_unbound_by_pli = as.factor(tPAPLiunbd)
  ,tpa_unbound_kinetic = as.factor(tPAunbind)
)

