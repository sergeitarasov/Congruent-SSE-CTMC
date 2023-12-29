#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// Helper function to fill off-diagonal elements of a square matrix
NumericMatrix fill_off_diagonal(NumericVector off_diagonal_elements, int N) {
  NumericMatrix matrix(N, N);
  int counter = 0;

  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      if (i != j) {
        matrix(i, j) = off_diagonal_elements[counter];
        counter++;
      }
    }
  }
  return matrix;
}

// Helper function to fill an upper triangular matrix from a vector
arma::mat vector_to_upper_triangular_matrix(NumericVector x, int Nstates) {
  arma::mat out(Nstates, Nstates, arma::fill::zeros);
  int k = 0;
  for (int i = 0; i < Nstates; ++i) {
    for (int j = i; j < Nstates; ++j) {
      out(i, j) = x[k++];
    }
  }
  return out;
}

// Function to split lambdas into Nstates chunks
std::vector<NumericVector> split_vector(NumericVector lambdas, int Nstates) {
  std::vector<NumericVector> out;
  int len = lambdas.size() / Nstates;
  for (int i = 0; i < Nstates; ++i) {
    out.push_back(lambdas[Range(i * len, (i + 1) * len - 1)]);
  }
  return out;
}

List cpp_lambdas_to_matrices(NumericVector lambdas, int Nstates) {
  std::vector<NumericVector> split_lambdas = split_vector(lambdas, Nstates);
  List matrices(Nstates);

  for (int i = 0; i < Nstates; ++i) {
    arma::mat mat = vector_to_upper_triangular_matrix(split_lambdas[i], Nstates);
    matrices[i] = wrap(mat);
  }

  return List::create(Named("matrices") = matrices);
}

// [[Rcpp::export]]
List cpp_pars_to_arrays(NumericVector pars, int Nstates) {
  int Nlambdas_one_state = 0.5 * (Nstates * Nstates + Nstates);
  int Nlambdas = Nlambdas_one_state * Nstates;

  NumericVector lambdas = pars[Range(0, Nlambdas - 1)];
  NumericVector mu = pars[Range(Nlambdas, Nlambdas + Nstates - 1)];
  NumericVector qs = pars[Range(Nlambdas + Nstates, pars.size() - 1)];

  // lam.tensor
  List lam_tensor = cpp_lambdas_to_matrices(lambdas, Nstates);

  // make Q matrix
  NumericMatrix Q = fill_off_diagonal(qs, Nstates);
  for (int i = 0; i < Nstates; ++i) {
    double row_sum = 0;
    for (int j = 0; j < Nstates; ++j) {
      row_sum += Q(i, j);
    }
    Q(i, i) = -row_sum;
  }

  return List::create(Named("lam.tensor") = lam_tensor["matrices"],
                      Named("mu") = mu,
                      Named("Q") = Q);
}



