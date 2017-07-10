#### model ####

stimulus <- function(t,s){
  if(0 <= t & t < 5){
    return(s)
  }
    return(0)
}

model_ode <-  function(p,y,t,s){
  dydt <- rep(0, times = 17)
  # dydt[1] =  p[4]*y[13]*2.0-((y[16]+y[17])*p[1]*y[1])/(p[6]+y[1])
  # dydt[2] =  p[2]*(y[2]*y[2])*-2.0+((y[16]+y[17])*p[1]*y[1])/(p[6]+y[1])
  # dydt[3] =  -p[3]*y[3]+p[2]*(y[2]*y[2])
  # dydt[4] =  p[3]*y[3]-p[4]*y[4]
  # dydt[5] =  p[4]*y[4]-p[4]*y[5]
  # dydt[6] =  p[4]*y[5]-p[4]*y[6]
  # dydt[7] =  p[4]*y[6]-p[4]*y[7]
  # dydt[8] =  p[4]*y[7]-p[4]*y[8]
  # dydt[9] =  p[4]*y[8]-p[4]*y[9]
  # dydt[10] =  p[4]*y[9]-p[4]*y[10]
  # dydt[11] =  p[4]*y[10]-p[4]*y[11]
  # dydt[12] =  p[4]*y[11]-p[4]*y[12]
  # dydt[13] =  p[4]*y[12]-p[4]*y[13]
  # dydt[14] =  p[3]*y[3]-p[4]*y[13]
  # dydt[15] =  p[8]*y[16]-stimulus(t,s)*p[5]*p[7]*y[15]
  # dydt[16] =  -p[8]*y[16]-p[9]*y[16]+stimulus(t,s)*p[5]*p[7]*y[15]
  # dydt[17] =  p[9]*y[16]-p[10]*y[17]
  
  dydt[ 1] =  p[4]*y[13]*2.0-((y[16]+y[17])*p[1]*y[1])/(p[6]+y[1]);
  dydt[ 2] =  p[2]*(y[2]*y[2])*-2.0+((y[16]+y[17])*p[1]*y[1])/(p[6]+y[1]);
  dydt[ 3] =  -p[3]*y[3]+p[2]*(y[2]*y[2]);
  dydt[ 4] =  p[3]*y[3]-p[4]*y[4];
  dydt[ 5] =  p[4]*y[4]-p[4]*y[5];
  dydt[ 6] =  p[4]*y[5]-p[4]*y[6];
  dydt[ 7] =  p[4]*y[6]-p[4]*y[7];
  dydt[ 8] =  p[4]*y[7]-p[4]*y[8];
  dydt[ 9] =  p[4]*y[8]-p[4]*y[9];
  dydt[ 10] =  p[4]*y[9]-p[4]*y[10];
  dydt[ 11] =  p[4]*y[10]-p[4]*y[11];
  dydt[ 12] =  p[4]*y[11]-p[4]*y[12];
  dydt[ 13] =  p[4]*y[12]-p[4]*y[13];
  dydt[ 14] =  p[3]*y[3]-p[4]*y[13];
  dydt[ 15] =  p[8]*y[16]-(stimulus(t,s)*p[5]*p[7]*(y[15]*y[15]))/(y[15]*y[15]+p[11]);
  dydt[ 16] =  -p[8]*y[16]-p[9]*y[16]+(stimulus(t,s)*p[5]*p[7]*(y[15]*y[15]))/(y[15]*y[15]+p[11]);
  dydt[ 17] =  p[9]*y[16]-p[10]*y[17];
  
  return(dydt)
}


model_trajectory <- function(data.trajectory, variables){
  data.trajectory.conditions <- data.trajectory %>% dplyr::distinct(priming, stimulation, time)
  data.derivatives <- data.frame(m = numeric(), time = numeric(), stimulation = numeric(), priming = numeric(), var = numeric())
  for(i in 1:nrow(data.trajectory.conditions)){
    data.derivatives.tmp <- data.frame(
      m = model_ode(p = parameters.model,
                   y = (data.trajectory %>% 
                          dplyr::filter(priming == as.numeric(data.trajectory.conditions[i,"priming"]),
                                        stimulation == as.numeric(data.trajectory.conditions[i,"stimulation"]),
                                        time == as.numeric(data.trajectory.conditions[i,"time"])))$m,
                   t = data.trajectory.conditions[i,]$time,
                   s = data.trajectory.conditions[i,]$stimulation),
      priming = data.trajectory.conditions[i,]$priming,
      stimulation = data.trajectory.conditions[i,]$stimulation,
      time = data.trajectory.conditions[i,]$time,
      var = 1:length(variables))
    data.derivatives <- rbind(data.derivatives, data.derivatives.tmp)  
  }
  return(data.derivatives)
}