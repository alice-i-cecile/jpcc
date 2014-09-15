# Libraries ####
library(foreach)


# Grabbing ages of trees ####
grab_age <- function(i, tra, sparse=TRUE)
{
  if (sparse)
  {
    age <- max(as.numeric(as.character(tra[tra$i==i,"a"])))
  } else{
    age_row <- min(which(!is.na(tra[i,,]), arr.ind=T)[[2]])
    age <- dimnames(tra)[[3]][age_row]
  }
  return(as.numeric(age))
}


# Grabbing birth years of trees ####
grab_birth_year <- function(i, tra, sparse=TRUE)
{
  if (sparse)
  {
    slice <- tra[tra$i==i,][1,]
    birth_year <- as.numeric(as.character(slice$t)) -  as.numeric(as.character(slice$a)) + 1
  } else{
    slice <- which(!is.na(tra[i,,]), arr.ind=T)[1,]
    birth_index <- slice[1]-slice[2]+1
    birth_year <- as.numeric(dimnames(tra)[[2]][birth_index])
  }
  
  return(birth_year)
}


get_birth_years <- function(tra, sparse=TRUE)
{
  
  if (sparse){
    birth_years <- sapply(levels(tra$i), grab_birth_year, tra=tra, sparse)
  } else
  {
    birth_years <- sapply(1:dim(tra)[1], grab_birth_year, tra=tra, sparse)
    names(birth_years) <- dimnames(tra)[[1]]
  }
  
  return (birth_years)
}


get_birth_index <- function(birth_year, tra, sparse)
{
  if (sparse){
    years <- sort(as.numeric(levels(tra$t)))
  } else {
    years <- sort(as.numeric(dimnames(tra)[[2]]))
  }
  
  year_index <- which(birth_year==years)
  age_index <- 1
  birth_index <- year_index - age_index
  
  return (birth_index)
}

# RWL to TRA and back again ####

# Get endpoints of a rwl series
getEndpoints <- function (series, side="start"){
  #side: "start" or "stop" of series
  rings <- subset(series, !is.na(series))
  if (side=="start"){
    return (head(rownames(rings),1))
  }
  if (side=="stop"){
    return (tail(rownames(rings),1))
  }      
}

# Convert a standard tree ring data frame (.rwl) into tree ring array form. Model is G(i, t, T) = Q*F*A
rwl.to.tra <- function (rwl, birth_years=NULL){
  
  if (is.null(birth_years)){
    # Determine birth for each tree
    birth_years <- foreach(i=colnames(rwl)) %do% {getEndpoints(rwl[i], side="start")}
    names (birth_years) <- names(rwl)
    
  }  else {
    # If rwl file is too short, add empty rows to it
    birth_limit <- min(sapply(birth_years, as.numeric))
    rwl_limit <- min(sapply(rownames(rwl), as.numeric))
    
    if (birth_limit < rwl_limit){
      num_missing <- rwl_limit-birth_limit
      
      empty_rows <- rbind(rwl[0,], matrix(NA, num_missing, ncol(rwl)))
      
      colnames(empty_rows) <- colnames(rwl)
      rownames(empty_rows) <- birth_limit:(rwl_limit-1)
      
      rwl <- rbind(empty_rows, rwl)
    }
  }
  
  # Find the indices for each year
  birth_index <- vector ()
  
  for (i in 1:length(birth_years)){
    birth_index[i] <- which(rownames(rwl)==toString(birth_years[i]))
  }
  
  names (birth_index) <- names (birth_years)
    
  # Compute dimension size for all elements
  i.size <- ncol (rwl)
  t.size <- nrow (rwl)
  T.size <- t.size
  
#   # Construct empty vectors of the appropriate length
#   Q.empty <- rep.int (NA, i.size)
#   F.empty <- rep.int (NA, t.size)
#   A.empty <- rep.int (NA, T.size)
#   
#   # Construct an empty tree ring array via tensor multiplication
#   tra <- Q.empty %o% F.empty %o% A.empty
  
  # Construct an array directly
  tra <- array(NA, c(i.size, t.size, T.size), list(colnames(rwl), rownames (rwl), 1:T.size))
  
  # Input data 
  for (tree in 1:ncol(rwl)){
    birth <- birth_index[tree]
    for (year in 1:nrow(rwl)){
      
      # Find age of ring data
      age <- year - birth + 1
      

      datum <- rwl [year, tree]
      
      # Only update filled values
      if (!is.na(datum)){        
        tra [tree,year,age] <- datum
      }
    }
  }
  
  # Truncate empty rows
  empty_tra <- is.na(tra)
  
  empty_Q <- apply(empty_tra, 1, all)
  empty_F <- apply(empty_tra, 2, all)
  empty_A <- apply(empty_tra, 3, all)
  
  tra <- tra[!empty_Q, !empty_F, !empty_A]
  
  return (tra)
}

# Convert a tree ring array back to a standard tree ring data frame
tra.to.rwl <- function (tra) {
  
  
  # Get dimension sizes
  i.size <- dim (tra)[1]
  t.size <- dim (tra)[2]
  
  # Make an empty data frame
  rwl <- as.data.frame(matrix (NA, t.size, i.size))
  
  colnames(rwl) <- dimnames(tra)[[1]]
  rownames(rwl) <- dimnames(tra)[[2]]
  
  # Fill data
  for (tree in 1:dim(tra)[1]){
    
    # Get lifespan to fill. Could be more general in style
    filled <- !is.na(tra[tree, ,])
    filledYears <- rownames(filled[rowSums(filled)>0,])
    names(filledYears) <- dimnames(tra)[[1]]
    
    # Extract data
    data <- tra[tree,,]
    data <- data[!is.na(data)]
    
    # Put data into appropriate position in data frame
    rwl[filledYears, tree] <- data
  }
  
  return (rwl)
}

# Handling sparse tree-ring arrays ####

# Compress a (full) tree-ring array to its sparse form
sparse_tra <- function(tra)
{
  full <- !is.na(tra)
  values <- tra[full]
  ind_pos <- which(full, arr.ind=T)
  
  i <- dimnames(full)[[1]][ind_pos[,1]]
  t <- dimnames(full)[[2]][ind_pos[,2]]
  a <- dimnames(full)[[3]][ind_pos[,3]]
  
  stra <- data.frame(G=values, i=i, t=t, a=a)
  stra[2:4] <- lapply(stra[2:4], as.factor)
  
  return (stra)
}


# Uncompress a (sparse) tree-ring array to its full form
unsparse_tra <- function(stra)
{
  # Sort levels of indexes to hopefully preserve temporal ordering
  indices <- lapply(stra[2:ncol(stra)], levels)
  indices[c("t", "a")] <- lapply(indices[c("t", "a")], as.numeric)
  indices <- lapply(indices, sort)
  indices[c("t", "a")] <- lapply(indices[c("t", "a")], as.character)
  
  # Generate an array of the appropriate size
  tra <- array(data=NA, dim=sapply(indices, length), dimnames=indices)
  
  # Fill in the recorded data  
  for (j in 1:nrow(stra))
  {
    stra_row <- stra[j,]
    
    i <- as.character(stra[j,"i"])
    t <- as.character(stra[j,"t"])
    a <- as.character(stra[j,"a"])
    G <- stra[j,"G"]
    
    tra[i,t,a] <- G
  }
  
  return(tra)
  
}

# Go directly from a rwl to stra
# Avoid unnecessary memory bottlenecks
# Increase chunk size as high as feasible for optimal performance
rwl.to.stra <- function (rwl, birth_years=NULL)
{
    if (is.null(birth_years)){
      # Determine birth for each tree
      birth_years <- foreach(i=colnames(rwl)) %do% {getEndpoints(rwl[i], side="start")}
      names (birth_years) <- names(rwl)
      
    }  else {
      # If rwl file is too short, add empty rows to it
      birth_limit <- min(sapply(birth_years, as.numeric))
      rwl_limit <- min(sapply(rownames(rwl), as.numeric))
      
      if (birth_limit < rwl_limit){
        num_missing <- rwl_limit-birth_limit
        
        empty_rows <- rbind(rwl[0,], matrix(NA, num_missing, ncol(rwl)))
        
        colnames(empty_rows) <- colnames(rwl)
        rownames(empty_rows) <- birth_limit:(rwl_limit-1)
        
        rwl <- rbind(empty_rows, rwl)
      }
    }
    
    # Find the indices for each year
    birth_index <- vector ()
    
    for (i in 1:length(birth_years)){
      birth_index[i] <- which(rownames(rwl)==toString(birth_years[i]))
    }
    
    names (birth_index) <- names (birth_years)
   
    # Find the full cells
    filled <- which(!is.na(rwl), arr.ind=TRUE)
    n_data_points <- nrow(filled)
     
    
    add_growth_data <- function (year, tree)
    {
      birth <- birth_index[tree]
      # Extract growth data
      G <- rwl [year, tree]
      
      # Find tree of growth data
      i <- names(rwl)[tree]
      
      # Find year of growth data
      t <- rownames(rwl)[year]
      
      # Find age of ring data
      a <- year - birth + 1
      
      out <- list(G, i, t, a)
      return(out)
    }
    
    # Find the growth data at each filled cell and record its position
    raw_stra <- mapply(add_growth_data, year=filled[,1], tree=filled[,2], SIMPLIFY=FALSE)
    stra <- data.frame(matrix(unlist(raw_stra),ncol=4, byrow=TRUE))
    names(stra) <- c("G", "i", "t", "a")
    stra[[1]] <- as.numeric(stra[[1]])
  
  return(stra)
  
}

# Geometric mean utility function ####
geomMean <- function(x){
  if (length(x)==0){
    return(NA)
  }
  val <- x[!is.na(x)]
  l_val <- log(val)
  out <- exp(mean(l_val))
  return (out)
}

# Get sample size along a dimension ####
sample_depth_tra <- function(tra, factor.dim=2, sparse=FALSE){ #1 is tree, 2 is time, 3 is age
  # Convert the tree-ring array to the appropriate form (sparse/full)
  if (sparse)
  {
    if (!is.data.frame(tra))
    {
      tra <- sparse_tra(tra)
    }
  } else
  {
    if (is.data.frame(tra))
    {
      tra <- unsparse_tra(tra)
    }
  }
  
  if(sparse)
  {
    positions <- tra[[factor.dim+1]]
    ids <- levels (positions)
    
    sample_depth <- sapply(ids, function(x){sum(positions==x)})
    
    if(factor.dim==1)
    {
      ordering <- sort(names(sample_depth))
    }  else
    {
      ordering <- sort(as.numeric(names(sample_depth)))
    }
   
    sample_depth <- sample_depth[sapply(ordering, as.character)]
  } else 
  {
    filled_cells <- !is.na(tra)
    
    sample_depth <- apply(filled_cells, factor.dim, sum)
  }
  
  return (sample_depth)
}
