
struct SBM
{
    mat Z;
    rowvec alpha;

    // For nodes covariates
    bool has_nodes_covariates{false};
    mat nodes_covariates;
    mat B;
    mat alphamat;

    SBM(Rcpp::List & membership_from_R)
    {
        mat origZ = membership_from_R["Z"];
        Z=origZ;
        double tol = TOL_P/Z.n_rows;
        boundaries(Z,tol,1-tol);
        Z /= repmat( sum(Z,1), 1, Z.n_cols );
        alpha = sum(Z,0) / Z.n_rows;
        if (membership_from_R.containsElementNamed("nodes_covariates")) {
            has_nodes_covariates = true;
            mat orig_nodes_covariates = membership_from_R["nodes_covariates"];
            nodes_covariates = orig_nodes_covariates;
        }
    }

    SBM& operator=(const SBM& orig)
    {
        Z=orig.Z;
        alpha=orig.alpha;
        has_nodes_covariates = orig.has_nodes_covariates;
        nodes_covariates = orig.nodes_covariates;
        B = orig.B;
        alphamat = orig.alphamat;
        
        return *this;
    }
    
    inline
    double entropy()
    {
        return - accu( Z % log(Z) );
    }

    template<class model_type, class network_type>
    inline
    void e_step(model_type & model, network_type & net)
    {
        double tol = TOL_P/Z.n_rows;
        double step_size;
        unsigned int niter=0;
        do
        {
            // lZ the new log(Z) without renormalization
            mat lZ = repmat(log(alpha),Z.n_rows,1);
            if (has_nodes_covariates && Z.n_cols > 1)
            {
                lZ = log(alphamat);
            }
            
            // with a template, should be specialized by the model if possible
            e_fixed_step(*this, model, net, lZ);
            
            // This operation should change nothing theorically (but
            // pratically...)
            lZ -= repmat( max(lZ,1), 1, lZ.n_cols );

            // After this lZ is not the new log(Z) but the new Z without
            // normalization
            lZ = exp(lZ);

            // normalization
            lZ /= repmat( sum(lZ,1), 1, lZ.n_cols );

            boundaries(lZ,tol,1-tol);
            lZ /= repmat( sum(lZ,1), 1, lZ.n_cols );

            step_size = max(max(abs( Z-lZ )));

            Z=lZ;

            niter++;

            #ifdef DEBUG_E
                fprintf(stderr,"E iteration: %i %f\n",niter,step_size);
            #endif

        } while( step_size>TOL_F && niter<10);
    }

    inline
    double m_step()
    {
        if (has_nodes_covariates && Z.n_cols > 1)
        {
            B = optimize_softmax(nodes_covariates, Z); //compute_B
            alphamat = softmax(nodes_covariates * B);

            alphamat = clamp(alphamat, MIN_VAL, 1.0 - MIN_VAL);
            alphamat /= repmat(sum(alphamat, 1), 1, alphamat.n_cols);
            alphamat = clamp(alphamat, MIN_VAL, 1.0 - MIN_VAL);
            alphamat /= repmat(sum(alphamat, 1), 1, alphamat.n_cols);

            mat cross_Z_alpha = Z * log(alphamat).t();

            return accu(cross_Z_alpha.diag());
        }

        alpha = sum(Z,0) / Z.n_rows;
        return accu(Z * log(alpha).t());
    }

    inline
    Rcpp::List export_to_R()
    {
        Rcpp::List values;
        values["Z"] = Z;
        if (has_nodes_covariates && Z.n_cols > 1)
        {
            values["alpha"] = alphamat;
            values["B"] = B;
        }
        else
        {
            values["alpha"] = alpha;
        }
        values["nodes_covariates"] = nodes_covariates;

        return values;
    }
        
};

struct SBM_sym : public SBM
{
    SBM_sym(Rcpp::List & membership_from_R) : SBM(membership_from_R) {}
};



