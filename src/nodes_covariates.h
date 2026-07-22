inline mat softmax(const mat &Z)
{
    mat Zshift = Z;
    vec maxZ = max(Z, 1); // row max

    Zshift.each_col() -= maxZ; // stability

    mat P = exp(Zshift);
    vec rowsum = sum(P, 1);

    P.each_col() /= rowsum;

    return P;
}

inline mat log_softmax(const mat &Z)
{
    vec maxZ = max(Z, 1);
    mat Zshift = Z;
    Zshift.each_col() -= maxZ;

    vec logsumexp = log(sum(exp(Zshift), 1)) + maxZ;

    mat logP = Z;
    logP.each_col() -= logsumexp;

    return logP;
}

inline mat build_full_B(const mat &Btilde)
{
    mat B(Btilde.n_rows, Btilde.n_cols + 1, fill::zeros);
    B.cols(0, Btilde.n_cols - 1) = Btilde;
    return B;
}

double objective_gradient(
    const mat &X,
    const mat &T,
    mat &Btilde,
    mat &grad)
{
    mat B = build_full_B(Btilde);

    mat Z = X * B;

    mat logP = log_softmax(Z);
    double L = accu(T % logP);

    mat P = softmax(Z);
    // To verify
    vec s = sum(T, 1);
    mat SP = P;
    SP.each_col() %= s;

    mat Gfull = X.t() * (T - SP);

    // keep only identifiable parameters
    grad = Gfull.cols(0, Btilde.n_cols - 1);

    return L;
}
double line_search_bfgs(
    const mat &X,
    const mat &T,
    const mat &Btilde,
    const mat &grad,
    const mat &direction,
    double L)
{
    double alpha = 1.0;
    const double c1 = 1e-4;
    const double slope = accu(grad % direction);

    if (!std::isfinite(slope) || slope <= 0.0 || !std::isfinite(L))
        return 0.0;

    for (int i = 0; i < 25; i++)
    {
        mat Bnew = Btilde + alpha * direction;
        mat gtmp;
        double Lnew = objective_gradient(X, T, Bnew, gtmp);

        if (!std::isfinite(Lnew) || !gtmp.is_finite()) {
            alpha *= 0.5;
            continue;
        }

        // Armijo condition for maximization
        if (Lnew >= L + c1 * alpha * slope)
            return alpha;

        alpha *= 0.5;
    }

    return 0.0;
}

mat approx_B(const mat &X, const mat &T) {

    vec ref_tau = T.col(T.n_cols-1);
    mat log_unscaled_T = log(T.each_col() / ref_tau);
    mat pseudo_inv = pinv(X);

    return pseudo_inv * log_unscaled_T;
};

mat optimize_softmax(
    const mat &X,
    const mat &T)
{
    int p = X.n_cols;
    int R = T.n_cols;
    int q = R - 1;
    int nparam = p * q;
    int last_iter = 0;
    #ifdef DEBUG_M
    int no_h_iter = 0;
    #endif

    // only R-1 columns optimized
    mat Btilde(p, q);
    mat grad;

    Btilde = approx_B(X, T).cols(1,q);

    double L = objective_gradient(X, T, Btilde, grad);

    if (!std::isfinite(L) || !grad.is_finite()) {
        Rcpp::warning("Optimization for nodes covariates started from non-finite objective/gradient.");

        return build_full_B(Btilde);
    }

    mat H = eye<mat>(nparam, nparam); // inverse Hessian approximation

    for (int iter = 0; iter < BFGS_ITER_MAX; iter++)
    {
        vec g = vectorise(grad);
        if (!g.is_finite()) {
            Rcpp::warning("Optimization for nodes covariates encountered non-finite gradient.");
            break;
        }



        if (norm(g, 2) < TOL_NODE_COV)
            break;

        vec d = H * g; // ascent direction (maximize)

        if (!d.is_finite()) {
            d = g;
            H.eye();
        }

        if (dot(g, d) <= 0.0)
        {
            d = g;
            H.eye();
        }

        mat direction = reshape(d, p, q);
        double step = line_search_bfgs(X, T, Btilde, grad, direction, L);

        if (!std::isfinite(step) || step <= 0.0)
            break;

        mat Bold = Btilde;
        vec gold = g;

        Btilde += step * direction;

        double Lnew = objective_gradient(X, T, Btilde, grad);
        vec gnew = vectorise(grad);

        if (!std::isfinite(Lnew) || !gnew.is_finite()) {
            Btilde = Bold;
            grad = reshape(gold, p, q);
            H.eye();
            break;
        }

        double rel_obj_diff = std::abs(Lnew - L) / (1.0 + std::abs(L));

        vec s = vectorise(Btilde - Bold);
        vec y = gold - gnew; // Because we maximize
        double ys = dot(y, s);

        if (std::isfinite(ys) && ys > 1e-12 && s.is_finite() && y.is_finite())
        {
            double rho = 1.0 / ys;
            mat I = eye<mat>(nparam, nparam);
            mat syT = s * y.t();
            mat ysT = y * s.t();
            H = (I - rho * syT) * H * (I - rho * ysT) + rho * (s * s.t());
            if (!H.is_finite()) {
                H.eye();
            }
        }
        else
        {
        #ifdef DEBUG_M
            no_h_iter += 1;
        #endif
            H.eye();
        }

        L = Lnew;
        last_iter = iter;

        if (rel_obj_diff < TOL_NODE_COV)
            break;

        if (norm(gnew, 2) < TOL_NODE_COV)
            break;
        #ifdef DEBUG_M
        if (iter <= 10  || iter == BFGS_ITER_MAX - 1) {
        Rcpp::Rcout << "Iteration " << iter + 1 << "/" << BFGS_ITER_MAX << ": L = " << Lnew << endl;
        Rcpp::Rcout << "<yk, sk> = " << ys << endl;
        Rcpp::Rcout << "grad(xk+1) - grad(xk) = " << y << endl;
        Rcpp::Rcout << "xk+1 - xk = alphak*pk = " << s << endl;
        Rcpp::Rcout << no_h_iter << " steps without H update on " << iter << " steps" << endl;
        Rcpp::Rcout << "H:" << H << endl;
        }
        #endif
    }

    if (last_iter >= BFGS_ITER_MAX) {
        Rcpp::warning("Optimization for nodes covariates did not converge after %i steps.", last_iter);
    }

    return build_full_B(Btilde);
}
