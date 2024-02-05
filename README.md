# Adaptive varational open system simulation
This repo hosts the companion code for the paper [Adaptive variational simulation for open quantum systems](https://arxiv.org/abs/2305.06915).
Within, the AVQD directory houses the core package that includes the entirety of the adaptive variational algorithm.
The root directory serves as the primary project space, specifically configured for execution in a cluster computing environment. 
All dependencies for these two projects are readily available for installation via Julia’s package repository, except for DM. DM is a custom package I developed for convenient data storage and retrieval. It can be manually installed from this GitHub repository: https://github.com/neversakura/DM.jl.

## Quick start
1. Build the database index
```julia
    julia --project=. src/init_data_repo.jl
```
2. Exact solution
```julia
    julia --project=. src/exact.jl n -l "z" -g 1e-3
```
    The arguments: 
    * `n` : the total number of qubits;
    * `l` : the type of Lindblad operator `z` or `pm`, defaults to `z`;
    * `g` : the decay rate for ``\sigma_z`` Lindbladian, only works if `l` is set to `z`, defaults to `0`;
    * `gamma_p` : the decay rate for `\sigma_+` Lindbladian, only works if `l` is set to `pm`, defaults to `0`; 
    * `gamma_m` : the decay rate for `\sigma_-` Lindbladian, only works if `l` is set to `pm`, defaults to `0`; 

3. Trajectories
```julia
    julia -p 7 src/trajectory.jl n -r 1e-3 -l "z" -g 1e-2 --num_t 1000 --pool "all2" --relative --reduce --no_state
```
    The arguments:
    * `r` : the McLachlan distance threshold, defaults to `1e-3`;
    * `num_t` : number of trajectories, defaults to 1000;
    * `pool` : ansartz pool type, defaults to `"all2"`;
    * `relative` : whether to use unrestricted adaptive protocol;
    * `reduce` : whether to reduce the density matrix into observables before saving;
    * `no_state` : whether to save the state of the quantum system;
    * `save_every` : Whether to save the ansatz parameters for every time step.

4. Vectorization
```julia
julia --project=. src/vectorized.jl n -r 1e-5 -g 0.01 --relative
```

The arguments are the same as the previous cases.
