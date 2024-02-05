using ArgParse

s = ArgParseSettings()
@add_arg_table s begin
    "--rcut", "-r"
        help = "McLachlan distance thresh hold"
        arg_type = Float64
        default = 1e-3
    "--relative"
        help = "Whether to use relative McLachlan distance threshold"
        action => :store_true
    "--reduce"
        help = "Whether to reduce the density matrix to energies"
        action => :store_true
    "--lind", "-l"
        help = "Type of Lindblad operators"
        arg_type = String
        default = "z"
    "--gamma", "-g"
        help = "γ for Z Lindbladian"
        arg_type = Float64
        default = 0.0
    "--gamma_p"
        help = "γ₊ for σ₊ Lindbladian"
        arg_type = Float64
        default = 0.0
    "--gamma_m"
        help = "γ₋ for σ₋ Lindbladian"
        arg_type = Float64
        default = 0.0
    "--pool", "-p"
        help = "Type of ansatz pool"
        arg_type = String
        default = "all2"
    "--num_t"
        help = "Number of trajectories"
        arg_type = Int
        default = 1000
    "--no_state"
        help = "Whether to save the state of the quantum system"
        action => :store_true
    "--save_every"
        help = "Whether to save the ansatz parameters for every time step"
        action => :store_true
    "nqubit"
        help = "Number of qubits"
        arg_type = Int
        required = true
end

parsed_args = parse_args(ARGS, s)