
struct LBM
{
    mat Z1;
    mat Z2;

    // In case there are no covariates on nodes
    rowvec alpha1;
    rowvec alpha2;

    // For covariates on nodes
    bool has_row_covariates{false};
    bool has_col_covariates{false};

    mat row_covariates;
    mat col_covariates;

    // Parameters for the covariates
    mat B;
    mat G;

    mat alpha1mat;
    mat alpha2mat;

    LBM(Rcpp::List &membership_from_R)
    {
        mat origZ1 = membership_from_R["Z1"];
        mat origZ2 = membership_from_R["Z2"];
        Z1 = origZ1;
        Z2 = origZ2;
        double tol1 = TOL_P / Z1.n_rows;
        double tol2 = TOL_P / Z2.n_rows;
        boundaries(Z1, tol1, 1 - tol1);
        boundaries(Z2, tol2, 1 - tol2);
        Z1 /= repmat(sum(Z1, 1), 1, Z1.n_cols);
        Z2 /= repmat(sum(Z2, 1), 1, Z2.n_cols);
        alpha1 = sum(Z1, 0) / Z1.n_rows;
        alpha2 = sum(Z2, 0) / Z2.n_rows;

        // Init of RNG
        // arma::arma_rng::set_seed_random();

        if (membership_from_R.containsElementNamed("row_covariates"))
        {
            // TODO Ask JBL about why this is needed and the simpler way wont compile?
            mat orig_row_covariates = membership_from_R["row_covariates"];
            row_covariates = orig_row_covariates;
        }
        if (row_covariates.n_elem > 0) {
            // Rcpp::Rcout << "\nI'm a well behaved row covariates matrix" << endl;
            has_row_covariates = true;
            alpha1mat = repmat(alpha1, Z1.n_rows, 1);
        }

        if (membership_from_R.containsElementNamed("col_covariates"))
        {
            mat orig_col_covariates = membership_from_R["col_covariates"];
            col_covariates = orig_col_covariates;
        }
        if (col_covariates.n_elem > 0) {
            // Rcpp::Rcout << "\nI'm a well behaved col covariates matrix" << endl;
            has_col_covariates = true;
            alpha2mat = repmat(alpha2, Z2.n_rows, 1);
        }
    }

    LBM &operator=(const LBM &orig)
    {
        Z1 = orig.Z1;
        Z2 = orig.Z2;
        alpha1 = orig.alpha1;
        alpha2 = orig.alpha2;
        alpha1mat = orig.alpha1mat;
        alpha2mat = orig.alpha2mat;
        has_row_covariates = orig.has_row_covariates;
        has_col_covariates = orig.has_col_covariates;
        B = orig.B;
        G = orig.G;
        row_covariates = orig.row_covariates;
        col_covariates = orig.col_covariates;

        return *this;
    }

    inline double entropy()
    {
        return -accu(Z1 % log(Z1)) - accu(Z2 % log(Z2));
    }

    template <class model_type, class network_type>
    inline void e_step(model_type &model, network_type &net)
    {
        double tol1 = TOL_P / Z1.n_rows;
        double tol2 = TOL_P / Z2.n_rows;
        double step_size;
        unsigned int niter = 0;

        do
        {
            mat lZ1 = repmat(log(alpha1), Z1.n_rows, 1);
            mat lZ2 = repmat(log(alpha2), Z2.n_rows, 1);
            if (has_row_covariates)
            {
                lZ1 = log(alpha1mat);
            }
            if (has_col_covariates)
            {
                lZ2 = log(alpha2mat);
            }

            e_fixed_step(*this, model, net, lZ1, lZ2);

            lZ1 -= repmat(mean(lZ1, 1), 1, lZ1.n_cols);
            lZ2 -= repmat(mean(lZ2, 1), 1, lZ2.n_cols);

            lZ1 -= repmat(max(lZ1, 1), 1, lZ1.n_cols);
            lZ2 -= repmat(max(lZ2, 1), 1, lZ2.n_cols);

            lZ1 = exp(lZ1);
            lZ2 = exp(lZ2);

            lZ1 /= repmat(sum(lZ1, 1), 1, lZ1.n_cols);
            lZ2 /= repmat(sum(lZ2, 1), 1, lZ2.n_cols);

            boundaries(lZ1, tol1, 1 - tol1);
            boundaries(lZ2, tol2, 1 - tol2);
            lZ1 /= repmat(sum(lZ1, 1), 1, lZ1.n_cols);
            lZ2 /= repmat(sum(lZ2, 1), 1, lZ2.n_cols);

            double step_size1 = max(max(abs(Z1 - lZ1)));
            double step_size2 = max(max(abs(Z2 - lZ2)));

            step_size = (step_size1 > step_size2) ? step_size1 : step_size2;

#ifdef DEBUG_E
            fprintf(stderr, "E iteration: %i %f [%f %f]\n", niter, step_size, step_size1, step_size2);
#endif

            niter++;

            Z1 = lZ1;
            Z2 = lZ2;
        } while (step_size > TOL_F && niter < 10);
    }

    inline double m_step()
    {
        double dim1_out = 0;
        double dim2_out = 0;
        if (has_row_covariates && Z1.n_cols > 1)
        {
            B = optimize_softmax(row_covariates, Z1);
            alpha1mat = softmax(row_covariates * B);

            alpha1mat = clamp(alpha1mat, MIN_VAL, 1.0 - MIN_VAL);
            alpha1mat /= repmat(sum(alpha1mat, 1), 1, alpha1mat.n_cols);
            alpha1mat = clamp(alpha1mat, MIN_VAL, 1.0 - MIN_VAL);
            alpha1mat /= repmat(sum(alpha1mat, 1), 1, alpha1mat.n_cols);

            // What we want to sum is the diagonal of the below matrix
            // when we multiply pi_{i,q} and tau_{i,q} not tau_{j,q}
            mat cross_Z1_alpha1 = Z1 * log(alpha1mat).t();
            dim1_out = accu(cross_Z1_alpha1.diag());
        }
        else
        {
            alpha1 = sum(Z1, 0) / Z1.n_rows;
            dim1_out = accu(Z1 * log(alpha1).t());
        }
        if (has_col_covariates && Z2.n_cols > 1)
        {
            G = optimize_softmax(col_covariates, Z2);
            alpha2mat = softmax(col_covariates * G);

            alpha2mat = clamp(alpha2mat, MIN_VAL, 1.0 - MIN_VAL);
            alpha2mat /= repmat(sum(alpha2mat, 1), 1, alpha2mat.n_cols);
            alpha2mat = clamp(alpha2mat, MIN_VAL, 1.0 - MIN_VAL);
            alpha2mat /= repmat(sum(alpha2mat, 1), 1, alpha2mat.n_cols);

            mat cross_Z2_alpha2 = Z2 * log(alpha2mat).t();
            dim2_out = accu(cross_Z2_alpha2.diag());
        }
        else
        {
            alpha2 = sum(Z2, 0) / Z2.n_rows;
            dim2_out = accu(Z2 * log(alpha2).t());
        }
        return dim1_out + dim2_out;
    }

    inline Rcpp::List export_to_R()
    {
        Rcpp::List values;
        values["Z1"] = Z1;
        if (has_row_covariates && Z1.n_cols > 1)
        {
            values["alpha1"] = alpha1mat;
            values["B"] = B;
        }
        else
        {
            values["alpha1"] = alpha1;
        }
        values["Z2"] = Z2;
        if (has_col_covariates && Z2.n_cols > 1)
        {
            values["alpha2"] = alpha2mat;
            values["G"] = G;
        }
        else
        {
            values["alpha2"] = alpha2;
        }
        values["row_covariates"] = row_covariates;
        values["col_covariates"] = col_covariates;

        return values;
    }
};
