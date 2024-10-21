
# number of unique Q matrices in the list
N_unique_Qs <- function(Qs){
  vector_list <-lapply(Qs, function(x) c(x))
  vector_matrix <- do.call(rbind, vector_list)
  unique_vectors <- unique(vector_matrix)
  nrow(unique_vectors)
}


# check if Q is trait independent, note the states should be sorted as: a0 b0 c0 d0 a1 b1 c1 d1
# works only for 8 states Qs
# is_trait_INdependent(Q_ehe8_C.t)
is_trait_INdependent <- function(Qin){
  Q <- initQ(c('a', 'b', 'c', 'd'), c(1:12))
  Q2<- initQ(c(0,1), c(13, 14))
  smm=amaSMM(Q2, Q)
  
  i=3
  # correlation according to Pagel
  out <- c()
  for (i in 1:12){
    cells <- which(smm==i)
    rates=Qin[cells]
    res <- all(rates==rates[1])
    out <- c(out, res)
  }
  
  # dual transitions
  #which((Q_ehe8_C.t + smm)==0)
  duals = c(7,  8, 15, 16, 21, 22, 29, 30, 35, 36, 43, 44, 49, 50, 57, 58)
  drates=Qin[duals]
  res <- all(drates==0)
  out <- c(out, res)
  
  return(all(out))
}


generateUpperTriangularMatrix <- function(element, N) {
  # Create an empty matrix
  matrix_result <- matrix(numeric(0), nrow = N, ncol = N)

  # Populate the matrix with indexed elements or 0 based on position
  # i=1
  # j=2
  for (i in 1:N) {
    for (j in 1:N) {
      if (i <= j) {
        if (class(element)=="character")
          matrix_result[i, j] <- paste0(element, i-1, j-1)
        else
          matrix_result[i, j] <- (element)
      } else {
        matrix_result[i, j] <- 0
      }
    }
  }

  return(matrix_result)
}
# Test
generateUpperTriangularMatrix('lam0', 3)
generateUpperTriangularMatrix(1, 3)



vector_to_upper_triangular_matrix <- function(N) {
  # Compute the side length of the matrix
  n <- (-1 + sqrt(1 + 8 * length(N))) / 2
  n <- as.integer(n)

  # Initialize a zero matrix of size n x n
  matrix <- matrix(0, n, n)

  counter <- 1
  for(i in 1:n) {
    for(j in i:n) {
      matrix[i, j] <- N[counter]
      counter <- counter + 1
    }
  }

  return(matrix)
}
# Test
vector_to_upper_triangular_matrix(c(1,2,3,4,5,6))


extract_upper_triangular <- function(matrix) {
  n <- nrow(matrix)
  upper_triangular_values <- c()

  for (i in 1:n) {
    for (j in i:n) {
      upper_triangular_values <- c(upper_triangular_values, matrix[i, j])
    }
  }

  return(upper_triangular_values)
}
# Test
u <- matrix(c(1,0,0,2,3,0,4,5,6), nrow=3, ncol=3)
extract_upper_triangular(u)
u <- generateUpperTriangularMatrix('lam0', 4)
u
u.ex <- extract_upper_triangular(u)
vector_to_upper_triangular_matrix(u.ex)



makeSymmetricMatrix <- function(upperTriangularMatrix) {
  # Add the upper triangular matrix to its transpose
  symmetricMatrix <- upperTriangularMatrix + t(upperTriangularMatrix)

  # Subtract the diagonal
  diag_values <- diag(upperTriangularMatrix)
  symmetricMatrix <- symmetricMatrix - diag(diag_values)

  return(symmetricMatrix)
}
# Test
u <- generateUpperTriangularMatrix(2, 3)
u
makeSymmetricMatrix(u)


make_matrix_diag2 <- function(N) {
  # Create an N x N matrix filled with 1s
  mat <- matrix(1, nrow=N, ncol=N)

  # Change the diagonal to 2
  diag(mat) <- 2

  return(mat)
}
# Test
make_matrix_diag2(3)


OffDiagonalHalve <- function(N) {
  # Create an N x N matrix filled with 1/2
  mat <- matrix(0.5, nrow=N, ncol=N)

  # Change the diagonal to 1
  diag(mat) <- 1

  return(mat)
}
# Test
OffDiagonalHalve(3)



create_mu_vector <- function(N) {
  mu_vector <- paste0("mu", 0:(N-1))
  return(mu_vector)
}
# Test
create_mu_vector(3)


generateRateMatrix <- function(N) {
  # Generate a matrix filled with symbols qij
  matrix <- matrix(0, nrow=N, ncol=N)
  for(i in 1:N) {
    for(j in 1:N) {
      if(i != j) {
        element_name <- paste0("q", i-1, j-1)
        matrix[i, j] <- element_name
      }
    }
  }

  # Set diagonal elements such that row sums are zero
  for(i in 1:N) {
    row_sum <- sum(ifelse(is.numeric(matrix[i, ]), matrix[i, ], 0))
    matrix[i, i] <- -row_sum
  }

  return(matrix)
}
generateRateMatrix(3)


extract_off_diagonal <- function(matrix) {
  N <- nrow(matrix) # Assuming it's a square matrix
  off_diagonal_elements <- c()

  for (i in 1:N) {
    for (j in 1:N) {
      if (i != j) {
        off_diagonal_elements <- c(off_diagonal_elements, matrix[i, j])
      }
    }
  }

  return(off_diagonal_elements)
}
# Test
QQ <- generateRateMatrix(3)
QQ
extract_off_diagonal(QQ)


fill_off_diagonal <- function(off_diagonal_elements, N) {
  # Determine N based on length of off_diagonal_elements
  #N <- (1 + sqrt(1 + 8 * length(off_diagonal_elements))) / 2
  #N <- as.integer(N)

  # Create an N x N zero matrix
  matrix <- matrix(0, nrow=N, ncol=N)

  counter <- 1

  for (i in 1:N) {
    for (j in 1:N) {
      if (i != j) {
        matrix[i, j] <- off_diagonal_elements[counter]
        counter <- counter + 1
      }
    }
  }

  return(matrix)
}
# Test
QQ <- generateRateMatrix(3)
QQ
extr <- extract_off_diagonal(QQ)
extr
fill_off_diagonal(extr, 3)


argnames_HiClaSSE <- function(Nstates){
  #Nstates <- 2
  # lam0 <- generateUpperTriangularMatrix('lam0', Nstates)
  # lam1 <- generateUpperTriangularMatrix('lam1', Nstates)
  #lam.tensor <- list(lam0, lam1)

  lam.tensor <- vector("list", length=Nstates)
  for (i in 1:Nstates){
    lamX <- generateUpperTriangularMatrix(paste0('lam', i-1), Nstates)
    lam.tensor[[i]] <- lamX
  }

  mu <- create_mu_vector(Nstates)
  Q <- generateRateMatrix(Nstates)

  lam.pars <- lapply(lam.tensor, function(x) extract_upper_triangular(x))
  lam.pars <- unlist(lam.pars)
  Q.pars <- extract_off_diagonal(Q)

  list(pars=c(lam.pars, mu, Q.pars), arrays=list(lam.tensor=lam.tensor, mu=mu, Q=Q))
}
# Test
argnames_HiClaSSE(3)



split_vector <- function(vec, K) {
  # Create a sequence repeating from 1 to (N/K) for each element of vec
  f <- rep(1:(length(vec)/K), each=K)

  # Split the vector based on the sequence
  list_of_vectors <- split(vec, f)

  return(list_of_vectors)
}
# Test
vec <- 1:20  # A vector with 20 elements
K <- 5       # We want to split it into vectors of size 5
result <- split_vector(vec, K)
result




pars_to_arrays <- function(pars, Nstates){
  # number of lambdas
  #0.5*(Nstates*Nstates-Nstates)+Nstates
  Nlambdas_one_state <- (0.5*(Nstates*Nstates+Nstates))
  Nlambdas <- Nlambdas_one_state*Nstates
  #Nqs <- Nstates*Nstates - Nstates

  lambdas <- pars[1:Nlambdas]
  mu <- pars[(Nlambdas+1):(Nlambdas+Nstates)]
  qs <- pars[-1:-(Nlambdas+Nstates)]

  lam.tensor <- split_vector(lambdas, Nlambdas_one_state)
  lam.tensor <- lapply(lam.tensor, function(x) vector_to_upper_triangular_matrix(x))
  Q <- fill_off_diagonal(qs, Nstates)
  diag(Q) <- -rowSums(Q)

  return(list(lam.tensor=lam.tensor, mu=mu, Q=Q))
}
# Test
Nstates = 3
argsHiClaSSE2 <- argnames_HiClaSSE(Nstates)
argsHiClaSSE2$pars
argsHiClaSSE2$arrays
pars <- c(1:27)
names(pars) <- argsHiClaSSE2$pars
pars
pars_to_arrays(pars, Nstates)



convert_to_array <- function(list_of_matrices) {
  # Check if the list is empty
  if (length(list_of_matrices) == 0) {
    stop("The list is empty.")
  }

  # Get the dimension of the first matrix
  matrix_dim <- dim(list_of_matrices[[1]])

  # Check if all matrices have the same dimension
  for (matrix in list_of_matrices) {
    if (!all(dim(matrix) == matrix_dim)) {
      stop("All matrices in the list must have the same dimension.")
    }
  }

  # Combine matrices into an array
  array_dim <- c(matrix_dim, length(list_of_matrices))
  result_array <- array(unlist(list_of_matrices), dim = array_dim)

  return(result_array)
}
# Example usage:
mat1 <- matrix(1:4, 2, 2)
mat2 <- matrix(5:8, 2, 2)
list_of_mats <- list(mat1, mat2)
convert_to_array(list_of_mats)


# reorder paramaters given partitioning scheme
reoder_lambdas <- function(arrays, v){
  tensor.neworder <- lapply(arrays$lam.tensor, function(x) x[v,v])
  tensor.neworder <- tensor.neworder[v]
  list(lam.tensor=tensor.neworder, mu=arrays$mu[v], Q=arrays$Q[v,v])
}
reoder_lambdas(argnames_HiClaSSE(4)$arrays, c(1,3, 2,4))



make_Qcol <- function(N, value_vector){
  if (length(value_vector) != N) {
    stop("Length of value_vector must be equal to N")
  }

  mat <- matrix(0, nrow = N, ncol = N)
  for (i in 1:N) {
    mat[, i] <- rep(value_vector[i], N)

  }

  diag(mat) <- 0
  diag(mat) <- -rowSums(mat)
  return(mat)
}
# Example usage
N <- 3
value_vector <- c(1, 2, 3)
make_Qcol(N, value_vector)


Lik <- function(y){
  lapply(y , function(x) x$loglik) %>% unlist
}

Aik <- function(y){
  lapply(y , function(x) x$AIC) %>% unlist
}


convert2ratesHisse <- function(lam, mu){
  turnover <- lam+mu
  eps <- mu/lam
  list(turnover=turnover, eps=eps)
}

convertHisse2Bisse <- function(turnover, eps) {
  # Calculate lambda
  lam <- turnover / (1 + eps)
    # Calculate mu
  mu <- eps * lam
  list(lam = lam, mu = mu)
}
#hi <- convert2ratesHisse(c(1,3), c(2,4))
#convertHisse2Bisse(hi$turnover, hi$eps)

# Function to set right diagonal elements to zero
set_right_diagonal_zero <- function(mat) {
  # Check if the matrix is square
  if (nrow(mat) != ncol(mat)) {
    stop("The matrix must be square.")
  }

  # Reverse the rows
  reversed_mat <- mat[nrow(mat):1,]

  # Get the indices of the diagonal elements in the reversed matrix
  diag_indices <- seq_len(nrow(mat)) + (seq_len(ncol(mat)) - 1) * nrow(mat)

  # Set the diagonal elements of the reversed matrix to zero
  reversed_mat[diag_indices] <- 0

  # Reverse the rows back to get the original matrix with right diagonal set to zero
  mat <- reversed_mat[nrow(reversed_mat):1,]

  return(mat)
}

# Test the function
my_matrix <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9), nrow = 3)
print("Original matrix:")
print(my_matrix)

new_matrix <- set_right_diagonal_zero(my_matrix)
print("Matrix with right diagonal set to zero:")
print(new_matrix)


# make arguments and global variables for HiClaSSE
makeArgs <- function(Args) {
  newArgs <- list(
    name = paste0('hiclasse', Args$Nstates),
    Nstates = as.integer(Args$Nstates),
    ny = as.integer(Args$Nstates*2L),
    idx.e = as.integer(c(1:Args$Nstates)),
    idx.d = as.integer(c(((Args$Nstates+1):(2*Args$Nstates)))),
    y = Args$y
  )

  ## Global vars to speed up likelihood computation
  # they are used internally by HiClasse
  MY_Nstates <<- as.integer(Args$Nstates)
  matrix_times_2    <<- make_matrix_diag2(as.integer(Args$Nstates))
  matrix_times_half <<- OffDiagonalHalve(as.integer(Args$Nstates))
  flat_vec          <<- rep(1, as.integer(Args$Nstates))

  return(newArgs)
}

# print global variables for HiClaSSE
printArgsGlobal <- function(){
  print('MY_Nstates:')
  print(MY_Nstates)
  print('matrix_times_2:')
  print(matrix_times_2)
  print('matrix_times_half:')
  print(matrix_times_half)
  print('flat_vec:')
  print(flat_vec)
}

# Args <- list(
#   Nstates = 4L,
#   y = list(
#     c(0,0,0,0, 1,1,0,0),
#     c(0,0,0,0, 0,0,1,1)
#   )
# )
# Args <- makeArgs(Args)
# printArgsGlobal()


reoder_pars <- function(pars.array, v){
  lam.tensor <- pars.array$lam.tensor
  mu <- pars.array$mu
  Q <- pars.array$Q

  lam.tensor <- lapply(lam.tensor, function(x) x[v,v])
  lam.tensor <- lam.tensor[v]
  mu <- mu[v]
  Q <-Q[v,v]

  list(lam.tensor=lam.tensor, mu=mu, Q=Q)
}

names2array <- function(pars.array, names){
  lam.tensor <- pars.array$lam.tensor
  mu <- pars.array$mu
  Q <- pars.array$Q

  lam.tensor <- lapply(lam.tensor, function(x) {rownames(x) <- colnames(x) <- names; x})
  names(lam.tensor) <- names
  names(mu) <- names
  rownames(Q) <- colnames(Q) <- names

  list(lam.tensor=lam.tensor, mu=mu, Q=Q)
}


get_item <- function(list, item, unlist=TRUE){
  list <- lapply(list, function(x) x[[item]])
  if (unlist) list <- unlist(list)
  list
}

get_aic <- function(vec, npar){
  2*npar - 2*vec
}

get_bic <- function(vec, npar, N){
  npar*log(N) - 2*vec
}

#X=pars.hc
formulas_zero_pars <- function(X) {
  # Find the names of the entities that have a value of 0
  zero_entities <- names(X[X == 0])

  # Create the output strings in the format 'entity_name ~ 0'
  output_strings <- paste0(zero_entities, " ~ 0")

  f <- lapply(output_strings, function(x) as.formula(x))
  return(f)
  # # Combine the output strings into a single line
  # output_line <- paste(output_strings, collapse = ", ")
  #
  # # Print the output line
  # print(output_line)
}


perturb_parameters <- function(params, mean = 0, sd = 0.1) {
  # Generate random perturbations from a normal distribution
  perturbations <- rnorm(length(params), mean, sd)

  # Add the perturbations to the original parameters
  perturbed_params <- params + perturbations

  # Ensure all values are >= 0
  perturbed_params <- pmax(perturbed_params, 0)
  perturbed_params <- perturbed_params + 0.00001

  return(perturbed_params)
}

# Example usage
original_params <- c(1.0, 2.0, 3.0, 4.0, 5.0)
perturb_parameters(original_params)


assign_classes_pairwise <- function(A, B) {
  # Initialize an empty list to store the classes
  classes <- list()

  # Initialize an empty vector to store zero-value pairs
  zero_pairs <- c()

  # Get the dimensions of the matrices
  n <- dim(A)[1]

  # Loop through the matrix A to find unique values and corresponding names in B
  for (i in 1:n) {
    for (j in 1:n) {
      # Skip NA and diagonal elements
      if (is.na(A[i, j]) || i == j) {
        next
      }

      # Handle zero-value pairs
      if (A[i, j] == 0) {
        zero_pairs <- c(zero_pairs, paste(B[i, j], "~ 0"))
        next
      }

      # Get the value from A
      value <- A[i, j]

      # Get the name from B
      name <- B[i, j]

      # Check if the value already exists in the classes list
      if (value %in% names(classes)) {
        # Append the name to the existing class
        classes[[as.character(value)]] <- c(classes[[as.character(value)]], name)
      } else {
        # Create a new class with the name
        classes[[as.character(value)]] <- c(name)
      }
    }
  }

  # Generate the output in formula notation
  output <- c()
  for (value in names(classes)) {
    #value="0.05"
    class_names <- classes[[value]]

    # Generate all pair-wise combinations for each class
    #pairs <- combn(class_names, 2, simplify = TRUE)
    if (length(class_names)>1){
      base <- class_names[1]
      formula <- paste(class_names[-1], paste0(" ~ ", base))
      output <- c(output, formula)
    }
  }

  # Add zero-value pairs to the output
  #output <- c(output, zero_pairs)
  output <- lapply(output, function(x) as.formula(x))

  return(output)
}

# Example usage
A <- matrix(c(NA, 0.05, 0.05, 0.00,
              0.05, NA, 0.00, 0.05,
              0.20, 0.00, NA, 0.05,
              0.00, 0.20, 0.05, NA), nrow = 4, byrow = TRUE)

B <- matrix(c("0", "q01", "q02", "q03",
              "q10", "0", "q12", "q13",
              "q20", "q21", "0", "q23",
              "q30", "q31", "q32", "0"), nrow = 4, byrow = TRUE)

assign_classes_pairwise(A, B)




# Function to get the right diagonal of a matrix
diagR <- function(mat) {
  if (nrow(mat) != ncol(mat)) {
    stop("The matrix must be square.")
  }

  right_diagonal <- numeric(nrow(mat))

  for (i in 1:nrow(mat)) {
    right_diagonal[i] <- mat[i, ncol(mat) - i + 1]
  }

  return(right_diagonal)
}

# Function to set the right diagonal of a matrix
`diagR<-` <- function(mat, value) {
  if (nrow(mat) != ncol(mat)) {
    stop("The matrix must be square.")
  }

  if (length(value) != nrow(mat)) {
    stop("The length of the value vector must match the dimensions of the square matrix.")
  }

  for (i in 1:nrow(mat)) {
    mat[i, ncol(mat) - i + 1] <- value[i]
  }

  return(mat)
}

# check matrix commutativity
Mcom <- function(A,B){
  A%*%B - B%*%A
}


#derivs.HiClasse_cpp_List(t=1, y=y, pars=pars.hc)
derivs.HiClasse_cpp_List <- function(t, y, pars){
  list(derivs.HiClasse_cpp(t, y, pars))
}


make_diag <- function(Q){
  diag(Q) <- 0
  diag(Q) <- - rowSums(Q)
  return(Q)
}

table_to_markdown <- function(vec) {
  headers <- names(vec)
  values <- vec
  header_line <- paste("|", paste(headers, collapse = " | "), "|")
  separator_line <- paste("|", paste(rep("---", length(headers)), collapse = " | "), "|")
  values_line <- paste("|", paste(values, collapse = " | "), "|")

  markdown_table <- paste(header_line, separator_line, values_line, sep = "\n")
  return(markdown_table)
}



library(igraph)

createGraph <- function(adj_matrix) {
  # Make sure the diagonal of the adjacency matrix is 0
  diag(adj_matrix) <- 0
  adj_matrix[adj_matrix>0] <- 1

  # Create a graph from the adjacency matrix
  graph <- graph.adjacency(adj_matrix, mode = "undirected")

  return(graph)
}

# Example usage:
# Assuming Qcon is your adjacency matrix
# graph <- createGraph(Qcon)


# counts number of extinct and alive species from diversitree simultation

count_species <- function(phy){
  sp <- phy$orig %>%   filter(!is.na(name))
  out <- c(sum(sp$extinct),  sum(!sp$extinct))
  names(out) <- c('extinct', 'alive')
  return(out)
}


library(ggplot2)
GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, draw_group = function(self, data, ..., draw_quantiles = NULL){
  data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
  grp <- data[1,'group']
  newdata <- plyr::arrange(transform(data, x = if(grp%%2==1) xminv else xmaxv), if(grp%%2==1) y else -y)
  newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
  newdata[c(1,nrow(newdata)-1,nrow(newdata)), 'x'] <- round(newdata[1, 'x']) 
  if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
    stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <= 
                                              1))
    quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
    aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
    aesthetics$alpha <- rep(1, nrow(quantiles))
    both <- cbind(quantiles, aesthetics)
    quantile_grob <- GeomPath$draw_panel(both, ...)
    ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
  }
  else {
    ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
  }
})

geom_split_violin <- function (mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, position = position, show.legend = show.legend, inherit.aes = inherit.aes, params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}

# make various parametrizations of Q with 4 states
ger_partitionsNK <- function(N, K_expression='<=2'){
  parts=partitions::setparts(N)
  n.parts=apply(parts, 2, function(x) length(unique(x)))
  #index=which(n.parts<=2)
  str=paste('which(n.parts', K_expression, ')')
  index <- eval(parse(text=str[1]))
  parts[,index]
}
make_Q4paramatrizations <- function(parts){
  Q1=initQ(c(0,1), c(1, 1))
  Q2 <- initQ(c('A','B'), c(1,1))
  tmp=amaSMM(Q2,Q1)
  #parts=setparts(8)
  index=which(tmp==1)
  
  Qout <- list()
  for (i in 1:ncol(parts)){
    Qnew <- tmp
    Qnew[index] <- parts[,i]
    Qout[[i]] <- Qnew
  }
  Qout
}
#set2 <- ger_partitionsNK(N=8, K_expression='==2')
#make_Q4paramatrizations(set2)




# Define the function to create a graph from the rate matrix
create_graph_from_matrix <- function(rate_matrix) {
  # Convert the rate matrix to a dataframe to easily handle row/column names
  rate_df <- as.data.frame(rate_matrix)
  
  # Get the row and column names
  nodes <- rownames(rate_df)
  
  # Create an empty graph
  g <- graph.empty(directed = TRUE)
  
  # Add nodes to the graph
  g <- add_vertices(g, nv = length(nodes), name = nodes)
  
  # Add edges to the graph (excluding main diagonal)
  for (i in seq_along(nodes)) {
    for (j in seq_along(nodes)) {
      if (i != j && rate_df[i, j] != 0) {
        g <- add_edges(g, c(nodes[i], nodes[j]), weight = rate_df[i, j])
      }
    }
  }
  
  return(g)
}

# Identify unidirectional edges
is_unidirectional <- function(g) {
  edges <- E(g)
  unidirectional <- logical(length(edges))
  
  for (edge in seq_along(edges)) {
    from <- ends(g, edge)[1]
    to <- ends(g, edge)[2]
    # Check if there is an edge in the opposite direction
    if (!are.connected(g, to, from)) {
      unidirectional[edge] <- TRUE
    }
  }
  
  return(unidirectional)
}

# g <- create_graph_from_matrix(Q_cid8.t)
# plot(g, vertex.label = V(g)$name, vertex.size = 30)
# g <- create_graph_from_matrix(Qc)
# unidirectional_edges <- is_unidirectional(g)
# plot(g, vertex.label = V(g)$name, vertex.size = 30, edge.color = ifelse(unidirectional_edges, "red", "black"))

# kronecker product for lambda tensor
kronecker_3d <- function(array1, array2) {
  result <- list()
  # Compute the Kronecker product for each combination of 2D matrices
  for (i in seq_along(array1)) {
    for (j in seq_along(array2)) {
      result[[length(result) + 1]] <- kronecker(array1[[i]], array2[[j]])
    }
  }
  return(result)
}




# Define the function
find_element <- function(list_of_matrices, element) {
  for (sublist_name in names(list_of_matrices)) {
    mat <- list_of_matrices[[sublist_name]]
    for (i in seq_len(nrow(mat))) {
      for (j in seq_len(ncol(mat))) {
        if (mat[i, j] == element) {
          return(c(sublist_name, rownames(mat)[i], colnames(mat)[j]))
        }
      }
    }
  }
  return(NULL)
}

#list_of_matrices=ar.re$lam.tensor
find_elements <- function(list_of_matrices, element) {
  res=matrix(NA, length(element), 4)
  #i=1
  for (i in 1:length(element)){
    el=element[i]
    nom=find_element(list_of_matrices, el)
    res[i,1] <- el
    res[i, 2:4] <- nom
  }
  return(res)
}

# element <- c("lam000", "lam001")
# find_elements(ar.re$lam.tensor, element)