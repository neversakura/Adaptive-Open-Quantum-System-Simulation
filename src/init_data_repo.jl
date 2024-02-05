using DM 

data_root = "data";

# "test_benchmark_1" is for the dephasing model 
root = "test_benchmark_1"

params = Dict(
    "prob" => "%s",
    "type" => "%s",
    "nqubit" => "%d",
    "tf" => "%d",
    "dt" => "%.2e",
    "lind" => "%s",
    "gamma" => "%.2e"
)

files = [["prob", "type"], ["nqubit", "tf", "dt"]]
groups = [["gamma"]]

name = root * "/exact"
entry = DataEntry(name, root, params, files, groups)
save_to_index_file(data_root, entry);


params = Dict(
    "prob" => "%s",
    "type" => "%s",
    "nqubit" => "%d",
    "tf" => "%d",
    "dt" => "%.2e",
    "lind" => "%s",
    "gamma" => "%.2e",
    "rcut" => "%.2e",
    "pool" => "%s"
)

files = [["prob", "type"], ["nqubit", "tf", "rcut"], ["pool"]]
groups = [["dt", "gamma"]]

name = root * "/trajectory"
entry = DataEntry(name, root, params, files, groups)
save_to_index_file(data_root, entry);

name = root * "/vectorized"
entry = DataEntry(name, root, params, files, groups)
save_to_index_file(data_root, entry);

# "test_benchmark_2" is for the generalized amplitude damping model 
root = "test_benchmark_2"

params = Dict(
    "prob" => "%s",
    "type" => "%s",
    "nqubit" => "%d",
    "tf" => "%d",
    "dt" => "%.2e",
    "lind" => "%s",
    "gamma_p" => "%.2e",
    "gamma_m" => "%.2e"
)

files = [["prob", "type"], ["nqubit", "tf", "dt"]]
groups = [["gamma_p", "gamma_m"]]

name = root * "/exact"
entry = DataEntry(name, root, params, files, groups)
save_to_index_file(data_root, entry);

params = Dict(
    "prob" => "%s",
    "type" => "%s",
    "nqubit" => "%d",
    "tf" => "%d",
    "dt" => "%.2e",
    "lind" => "%s",
    "gamma_p" => "%.2e",
    "gamma_m" => "%.2e",
    "rcut" => "%.2e",
    "pool" => "%s"
)

files = [["prob", "type"], ["nqubit", "tf", "rcut"], ["pool"]]
groups = [["dt", "gamma_p", "gamma_m"]]

name = root * "/trajectory"
entry = DataEntry(name, root, params, files, groups)
save_to_index_file(data_root, entry);

name = root * "/vectorized"
entry = DataEntry(name, root, params, files, groups)
save_to_index_file(data_root, entry);