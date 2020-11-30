
topology_dir = "dn250_ic1"
output_file = "test_data.csv"
import FEM_module as fem
from csv_to_modelica import convert_to_modelica_package

topology = fem.loadTopologyDATA(topology_dir)
assembly = fem.SimulationAssembly(topology, idata=fem.default_test_data)
callback = fem.Sim_Callback(result_dir=topology_dir, CSVFile=output_file)
fem.simulate(assembly, callback=callback, dt_max=2000)
if fem.mpi_rank == 0:
    convert_to_modelica_package(f"{topology_dir}/{output_file}", dx=1000)

