
class bernoulli
{
    public:

    class network
    {
        public:
        mat adj; // adjacency matrix

        // For NA handling
        mat maskNA;
        mat maskNAt;

        // precalculated matrices for SBM
        mat adjZD; 
        mat ones_minus_adj_ZD;
        mat adjZDt;
        mat ones_minus_adj_ZDt;
        mat onesZD;

        // precalculated matrices for LBM
        mat ones_minus_adj;
        mat adjt;
        mat ones_minus_adjt;
        mat adj_ones;


        network(Rcpp::List & network_from_R)
        {
            mat adj_orig = network_from_R["adjacency"];
            double na_replace_value = 0;
            if(network_from_R.containsElementNamed("na_replace_value"))
            {
                na_replace_value = Rcpp::as<double>(network_from_R["na_replace_value"]);
            }

            maskNA = compute_mask(adj_orig);

            maskNAt = maskNA.t();

            adj = replace_missing_values(adj_orig, na_replace_value); // Replace NAs to perform
            adjZD = fill_diag(adj,0);
            ones_minus_adj_ZD = fill_diag(1-adj,0);
            adjZDt = adjZD.t();
            ones_minus_adj_ZDt = ones_minus_adj_ZD.t();
            onesZD = fill_diag(ones<mat>(adj.n_rows, adj.n_rows),0);

            ones_minus_adj = 1-adj;
            adjt = adj.t();
            ones_minus_adjt = ones_minus_adj.t();
            adj_ones = ones<mat>(adj.n_rows, adj.n_cols);

        }
    };

    /* parameters */
    unsigned int n_parameters;
    mat pi;

    bernoulli(SBM & membership, bernoulli::network & net)
    {
        n_parameters = membership.Z.n_cols * membership.Z.n_cols;
        pi.set_size(membership.Z.n_cols,membership.Z.n_cols);
    }
    
    bernoulli(SBM_sym & membership, bernoulli::network & net)
    {
        n_parameters = membership.Z.n_cols * (membership.Z.n_cols+1) / 2;
        pi.set_size(membership.Z.n_cols,membership.Z.n_cols);
    }
    
    bernoulli(LBM & membership, bernoulli::network & net)
    {
        n_parameters = membership.Z1.n_cols * membership.Z2.n_cols;
        pi.set_size(membership.Z1.n_cols,membership.Z2.n_cols);
    }

    inline
    Rcpp::List export_to_R()
    {
        Rcpp::List values;
        values["pi"] = pi;
        values["n_parameters"] = n_parameters;

        return values;
    }
};

template <>
inline void e_fixed_step(SBM &membership, bernoulli &model, bernoulli::network &net, mat &lZ)
{
    lZ += (net.adjZD % net.maskNA) * membership.Z * log(model.pi).t() + (net.ones_minus_adj_ZD % net.maskNA) * membership.Z * log(1 - model.pi).t() + (net.adjZDt % net.maskNAt) * membership.Z * log(model.pi) + (net.ones_minus_adj_ZDt % net.maskNAt) * membership.Z * log(1 - model.pi);
}

template<>
inline
void e_fixed_step(SBM_sym & membership, bernoulli & model, bernoulli::network & net, mat & lZ)
{
    lZ+= (net.adjZD%net.maskNA) * membership.Z * log(model.pi)
       + (net.ones_minus_adj_ZD%net.maskNA) * membership.Z * log(1-model.pi);
}

template<>
inline
void e_fixed_step(LBM & membership, bernoulli & model, bernoulli::network & net, mat & lZ1, mat & lZ2)
{
    lZ1 += (net.adj%net.maskNA) * membership.Z2 * log(model.pi).t()
         + (net.ones_minus_adj%net.maskNA) * membership.Z2 * log(1-model.pi).t();
    lZ2 += (net.adjt%net.maskNAt) * membership.Z1 * log(model.pi)
         + (net.ones_minus_adjt%net.maskNAt) * membership.Z1 * log(1-model.pi);
}

template<>
inline
double m_step(SBM & membership, bernoulli & model, bernoulli::network & net)
{
    model.pi = (membership.Z.t() * (net.adjZD%net.maskNA) * membership.Z)
                /
               (membership.Z.t() * (net.onesZD%net.maskNA) * membership.Z);

    return
        (
            accu(
                    (log(model.pi)-log(1-model.pi)) 
                    % 
                    (membership.Z.t() * (net.adjZD%net.maskNA) * membership.Z)
                )
            +
            accu(
                    log(1-model.pi)
                    %
                    (membership.Z.t() * (net.onesZD%net.maskNA) * membership.Z)
                )
        );
}

template<>
inline
double m_step(SBM_sym & membership, bernoulli & model, bernoulli::network & net)
{
    return m_step<SBM>(membership,model,net)/2;
}

template<>
inline
double m_step(LBM & membership, bernoulli & model, bernoulli::network & net)
{
    model.pi = (membership.Z1.t() * (net.adj%net.maskNA) * membership.Z2)
                /
               (membership.Z1.t() * (net.adj_ones%net.maskNA) * membership.Z2);

    return
        (
            accu(
                    log(model.pi) 
                    % 
                    (membership.Z1.t() * (net.adj%net.maskNA) * membership.Z2)
                )
            +
            accu(
                    log(1-model.pi)
                    %
                    (membership.Z1.t() * (net.ones_minus_adj%net.maskNA) * membership.Z2)
                )
        );
}


