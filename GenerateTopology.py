from TopologyGenerators import PipeData, TwoBuriedPipes, saveTopologyDATA

# Pipe Data
P = [PipeData(-0.25, -1, 0.250, name="Supply", t=[0.005, 0.059, 0.011],
              rho=[7850, 60, 960], cp=[525, 2000, 2300], lmbd=[50, 0.026, 0.43]),
     PipeData(0.25, -1, 0.250, name="Return", t=[0.005, 0.059, 0.011],
              rho=[7850, 60, 960], cp=[525, 2000, 2300], lmbd=[50, 0.026, 0.43])]

# Generate Topology
mesh_resolution = 300
topology = TwoBuriedPipes(P, mesh_resolution)
saveTopologyDATA(topology, "dn250_ic1")
