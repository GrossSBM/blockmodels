/*
 * Multiplex binary network with independent layers.
 *
 * The network has d layers (adjacency matrices) on the same nodes. Given the
 * blocks q of i and l of j, the d binary variables X^k_{ij} are independent
 * and X^k_{ij} ~ Bernoulli(pi^k_{ql}): there is one connectivity matrix per
 * layer and no interaction between the layers.
 *
 * (By contrast, bernoulli_multiplex models the joint distribution of the
 * 2^d combinations of values over the layers.)
 */

class bernoulli_multiplex_independent
{
    public:

    class network
    {
        public:
        cube adj; // n x n x d (or n1 x n2 x d for LBM)

        // precomputed matrices (ZD: zero diagonal, for SBM and SBM_sym)
        cube adjZD;
        cube adjZDt;
        cube adjt;
        mat Mones;
        mat MonesZD;

        network(Rcpp::List & network_from_R)
        {
            Rcpp::List adj_list = network_from_R["adjacency"];
            mat first_mat = Rcpp::as<mat>( adj_list[0] );
            adj.set_size(first_mat.n_rows, first_mat.n_cols, adj_list.size());
            for(int k=0; k<adj_list.size(); k++)
                adj.slice(k) = Rcpp::as<mat>(adj_list[k]);

            Mones = ones<mat>(first_mat.n_rows, first_mat.n_cols);
            MonesZD = fill_diag(Mones,0);

            adjZD.set_size(adj.n_rows, adj.n_cols, adj.n_slices);
            adjZDt.set_size(adj.n_cols, adj.n_rows, adj.n_slices);
            adjt.set_size(adj.n_cols, adj.n_rows, adj.n_slices);
            for(unsigned int k=0; k<adj.n_slices; k++)
            {
                adjZD.slice(k) = fill_diag(adj.slice(k),0);
                adjZDt.slice(k) = adjZD.slice(k).t();
                adjt.slice(k) = adj.slice(k).t();
            }
        }
    };

    // parameters
    unsigned int n_parameters;
    cube pi; // pi.slice(k) is the connectivity matrix of the k-th layer

    bernoulli_multiplex_independent(SBM & membership, bernoulli_multiplex_independent::network & net)
    {
        n_parameters = membership.Z.n_cols * membership.Z.n_cols * net.adj.n_slices;
        pi.set_size(membership.Z.n_cols, membership.Z.n_cols, net.adj.n_slices);
    }

    bernoulli_multiplex_independent(SBM_sym & membership, bernoulli_multiplex_independent::network & net)
    {
        n_parameters = membership.Z.n_cols * (membership.Z.n_cols+1) / 2 * net.adj.n_slices;
        pi.set_size(membership.Z.n_cols, membership.Z.n_cols, net.adj.n_slices);
    }

    bernoulli_multiplex_independent(LBM & membership, bernoulli_multiplex_independent::network & net)
    {
        n_parameters = membership.Z1.n_cols * membership.Z2.n_cols * net.adj.n_slices;
        pi.set_size(membership.Z1.n_cols, membership.Z2.n_cols, net.adj.n_slices);
    }

    inline
    Rcpp::List export_to_R()
    {
        Rcpp::List values;
        values["n_parameters"] = n_parameters;

        Rcpp::List pi_as_list(pi.n_slices);
        for(unsigned int k=0; k<pi.n_slices; k++)
            pi_as_list[k] = pi.slice(k);
        values["pi"] = pi_as_list;

        return values;
    }

    /* Keep this. When this is not usefull, this is not used. */
    inline vec to_vector();
    bernoulli_multiplex_independent(SBM &, const vec &);
    bernoulli_multiplex_independent(SBM_sym &, const vec &);
    bernoulli_multiplex_independent(LBM &, const vec &);
};



/* log(p) where -Inf (p=0, a layer without any edge, or log(1-p) for a layer
 * full of edges) is replaced by a finite value, so that the terms multiplied
 * by a null count stay null instead of giving 0 * (-Inf) = NaN. */
inline
mat bmi_safe_log(const mat & p)
{
    mat l = log(p);
    const double lmin = log(std::numeric_limits<double>::min());
    l.transform( [lmin](double v) { return (v < lmin) ? lmin : v; } );
    return l;
}

/* sum of A log(p) + (O-A) log(1-p), with the convention 0 log(0) = 0 */
inline
double bmi_bernoulli_PL(const mat & p, const mat & A, const mat & O)
{
    double PL = 0;
    for(unsigned int i=0; i<p.n_elem; i++)
    {
        if(A(i) > 0)
            PL += A(i) * log(p(i));
        if(O(i) - A(i) > 0)
            PL += (O(i) - A(i)) * log(1 - p(i));
    }
    return PL;
}



/* E step
 *
 * For each layer k, with A = adjacency, O = ones (without diagonal for SBM):
 *   A Z log(pi)' + (O-A) Z log(1-pi)'
 *     = A Z (log(pi) - log(1-pi))' + O Z log(1-pi)'
 */

template<>
inline
void e_fixed_step(SBM & membership,
                  bernoulli_multiplex_independent & model,
                  bernoulli_multiplex_independent::network & net,
                  mat & lZ)
{
    mat OZ = net.MonesZD * membership.Z; // MonesZD is symmetric
    for(unsigned int k=0; k<net.adj.n_slices; k++)
    {
        mat lp = bmi_safe_log(model.pi.slice(k));
        mat l1mp = bmi_safe_log(1 - model.pi.slice(k));
        mat logit = lp - l1mp;
        lZ += net.adjZD.slice(k) * membership.Z * logit.t()
            + net.adjZDt.slice(k) * membership.Z * logit
            + OZ * (l1mp.t() + l1mp);
    }
}

template<>
inline
void e_fixed_step(SBM_sym & membership,
                  bernoulli_multiplex_independent & model,
                  bernoulli_multiplex_independent::network & net,
                  mat & lZ)
{
    mat OZ = net.MonesZD * membership.Z;
    for(unsigned int k=0; k<net.adj.n_slices; k++)
    {
        mat lp = bmi_safe_log(model.pi.slice(k));
        mat l1mp = bmi_safe_log(1 - model.pi.slice(k));
        lZ += net.adjZD.slice(k) * membership.Z * (lp - l1mp)
            + OZ * l1mp;
    }
}

template<>
inline
void e_fixed_step(LBM & membership,
                  bernoulli_multiplex_independent & model,
                  bernoulli_multiplex_independent::network & net,
                  mat & lZ1,
                  mat & lZ2)
{
    mat OZ2 = net.Mones * membership.Z2;
    mat OtZ1 = net.Mones.t() * membership.Z1;
    for(unsigned int k=0; k<net.adj.n_slices; k++)
    {
        mat lp = bmi_safe_log(model.pi.slice(k));
        mat l1mp = bmi_safe_log(1 - model.pi.slice(k));
        mat logit = lp - l1mp;
        lZ1 += net.adj.slice(k) * membership.Z2 * logit.t()
             + OZ2 * l1mp.t();
        lZ2 += net.adjt.slice(k) * membership.Z1 * logit
             + OtZ1 * l1mp;
    }
}



/* M step: explicit maximum, layer by layer */

template<>
inline
double m_step(SBM & membership,
              bernoulli_multiplex_independent & model,
              bernoulli_multiplex_independent::network & net)
{
    mat ZtOZ = membership.Z.t() * net.MonesZD * membership.Z;
    double PL = 0;
    for(unsigned int k=0; k<net.adj.n_slices; k++)
    {
        mat ZtAZ = membership.Z.t() * net.adjZD.slice(k) * membership.Z;
        model.pi.slice(k) = ZtAZ / ZtOZ;
        PL += bmi_bernoulli_PL(model.pi.slice(k), ZtAZ, ZtOZ);
    }
    return PL;
}

template<>
inline
double m_step(SBM_sym & membership,
              bernoulli_multiplex_independent & model,
              bernoulli_multiplex_independent::network & net)
{
    return m_step<SBM>(membership,model,net)/2;
}

template<>
inline
double m_step(LBM & membership,
              bernoulli_multiplex_independent & model,
              bernoulli_multiplex_independent::network & net)
{
    mat Z1tOZ2 = membership.Z1.t() * net.Mones * membership.Z2;
    double PL = 0;
    for(unsigned int k=0; k<net.adj.n_slices; k++)
    {
        mat Z1tAZ2 = membership.Z1.t() * net.adj.slice(k) * membership.Z2;
        model.pi.slice(k) = Z1tAZ2 / Z1tOZ2;
        PL += bmi_bernoulli_PL(model.pi.slice(k), Z1tAZ2, Z1tOZ2);
    }
    return PL;
}
