from fenics import *
import time
import csv
import numpy
import os
from input_data import default_training_data

parameters["std_out_all_processes"] = False
parameters["form_compiler"]["optimize"] = False
CW = MPI.comm_world
mpi_rank = MPI.rank(CW)

class SimulationAssembly():
    def __init__(self, td, idata=default_training_data):
        
        self.idata = idata
        self.T_ws = self.idata.T_ws
        self.T_wr = self.idata.T_wr
        self.T_a = self.idata.T_a
        self.alpha_s = self.idata.alpha_s
        self.alpha_r = self.idata.alpha_r
        self.alpha_a = self.idata.alpha_a

        # Definitinon of Function Space and facet normal vectors
        V = FunctionSpace(td.mesh, 'CG', 1)
        n = FacetNormal(td.mesh)

        # Dirichlet boundary conditions
        self.bcs = []
        if len(td.DirichletList) != 0:
            if type(td.DirichletList[0]) != list:
                td.DirichletList = [td.DirichletList]

        # Create list of Dirichlet BC according to TopologyData td
        for item in td.DirichletList:
            value = Constant(item[0])
            boundary_idx = int(item[1])
            self.bcs.append(DirichletBC(V, value, td.boundaries, boundary_idx))

        # Main Assembly
        u = TrialFunction(V)
        v = TestFunction(V)
        self.M = assemble((td.rho*td.cp*u*v)*dx)  # mass matrix
        self.K = assemble(inner(td.lmbd*grad(u), grad(v))*dx)  # conduc. matrix

        # Assembly for Robin Bcs
        ds = Measure('ds', domain=td.mesh, subdomain_data=td.boundaries)
        self.B_s = assemble((self.alpha_s*u*v)*ds(0))  # Γ_s
        self.B_r = assemble((self.alpha_r*u*v)*ds(1))  # Γ_r
        self.B_a = assemble((self.alpha_a*u*v)*ds(2))  # Γ_a
        self.BK = self.B_s + self.B_r + self.B_a + self.K  # convenient sum

        # Ghost values of 3rd order BCs are spatialy independent so we premultiply
        T_ones = interpolate(Constant(1), V)
        self.B_s = self.B_s*T_ones.vector()
        self.B_r = self.B_r*T_ones.vector()
        self.B_a = self.B_a*T_ones.vector()

        # Evaluate areas of the selected boundary surfaces
        self.Area_s = assemble(1.0*ds(0))  # area of Γ_s
        self.Area_r = assemble(1.0*ds(1))  # area of Γ_r
        self.Area_a = assemble(1.0*ds(2))  # area of Γ_a

        # Transformations for avg value evaluation
        # usage: (Avg_i*field.vector()).sum()
        self.Avg_s = assemble((u*v/self.Area_s)*ds(0))  # Γ_s
        self.Avg_r = assemble((u*v/self.Area_r)*ds(1))  # Γ_r
        self.Avg_a = assemble((u*v/self.Area_a)*ds(2))  # Γ_a

        # Transormation for flux evaluation through surfaces with Dirichlet BCs
        # usage: (flux_dbcs*field.vector()).sum()
        self.flux_dbcs = assemble((Constant(0)*u*v)*ds)  
        for item in td.DirichletList:
            boundary_idx = int(item[1])
            self.flux_dbcs += assemble((-td.lmbd*inner(grad(u), v*n))*ds(boundary_idx))

        self.T_prev = Function(V, name="Temperature", label="Celsius")  # previous step temperature field
        self.T = Function(V, name="Temperature", label="Celsius")  # new step temperature field

        self.t = 0.0  # current evaluated/considered time point
        self.t_prev = 0.0
        self.t_stop = 0.0  # maximum simulatin time
        self.dt = 0.0  # curent taken/considered time step
        self.Elapsed = 0.0  # elapsed real time after main loop entry
        self.E_rate = 0.0  # rate in which heat inside domain builds up
        self.T_avg_s = 0.0  # average temperature at Γ_s
        self.T_avg_r = 0.0  # average temperature at Γ_r
        self.T_avg_a = 0.0  # average temperature at Γ_a
        self.T_ws_value = 0.0  # average temperature of water in Supply pipe
        self.T_wr_value = 0.0  # average temperature of water in Return pipe
        self.T_a_value = 0.0  # average temperature of air
        self.Q_s = 0.0  # Total Heat Flow through Γ_s
        self.Q_r = 0.0  # Total Heat Flow through Γ_r
        self.Q_a = 0.0  # Total Heat Flow through Γ_a
        self.E_error = 0.0  # Total Energy Flow Error (W)
        self.callbackDATA = {}  # dictionary for possible arbitrary callBack related (temporary) data

class Sim_Callback():
    def __init__(self, result_dir="result", CSVFile="OutputData.csv", VTKFile=None):
        self.CSVFile = CSVFile
        self.VTKFile = VTKFile
        self.Dir = result_dir

    def FirstCall(self, sim):
        if mpi_rank == 0 and self.CSVFile is not None:

            # create folder and csv file
            if not os.path.isdir(self.Dir):
                os.mkdir(self.Dir)
            self.Log = open(self.Dir + f"/{self.CSVFile}", "w")
            self.writer = csv.writer(self.Log, delimiter=',', lineterminator='\n',)
            
            # write headers to csv
            self.writer.writerow(["time", "T_ws", "T_wr", "T_a", "T_avg_s",
                                  "T_avg_r", "T_avg_a", "Q_s",  "Q_r", "Q_a",
                                  "Q_dbcs", "E_error"])
            # write first data
            self.writer.writerow([sim.t, sim.T_ws_value, sim.T_wr_value,
                                  sim.T_a_value, sim.T_avg_s, sim.T_avg_r,
                                  sim.T_avg_a, sim.Q_s,  sim.Q_r, sim.Q_a,
                                  sim.E_error, sim.Q_dbcs])

        if self.VTKFile is not None:
            # create file for parallel writing of the tempe field
            self.Temperature = File(f"{self.Dir}/VTK/{self.VTKFile}")

    def __call__(self, sim):
        # print info state to Console and save selected to CSV (using ANSI escape chars.)
        if mpi_rank == 0:
            progress = 100*sim.t/sim.t_stop
            filled = int(0.5*progress)
            unfilled = 50 - filled
            RTime = (sim.t_stop-sim.t)/(sim.t/sim.Elapsed) if sim.t > 0 else 0.0
            print("Current Sim State:")
            print("\033[K[" + "#" * filled + "-" * unfilled + "]: %.1f %%" % progress)
            print("\033[KTime(s)// Remaining:%.1f, Elapsed:%.1f, Sim:%.1f, dt:%.3f" % (RTime, sim.Elapsed, sim.t, sim.dt))
            print("\033[KTemp(C)// Supply:%.3f, Return:%.3f, Ground:%.3f" % (sim.T_avg_s, sim.T_avg_r, sim.T_avg_a))
            print("\033[KFlux(W)// Supply:%.3f, Return:%.3f, Ground:%.3f, Dbcs:%.3f" % (sim.Q_s, sim.Q_r, sim.Q_a, sim.Q_dbcs))
            print("\033[KEnergy Rate(W)// Error Norm: %.3f" % sim.E_error, end="\033[F"*5)

            if self.CSVFile is not None:
                self.writer.writerow([sim.t, sim.T_ws_value, sim.T_wr_value,
                                      sim.T_a_value, sim.T_avg_s, sim.T_avg_r,
                                      sim.T_avg_a, sim.Q_s,  sim.Q_r, sim.Q_a,
                                      sim.E_error, sim.Q_dbcs])
        if self.VTKFile is not None:
            self.Temperature << sim.T_prev

    def close(self, sim):
        if mpi_rank == 0 and self.CSVFile is not None:
            self.Log.close()
        if mpi_rank == 0:
            print("\033[B"*5 + "\nSimulation Finished Sucesfully")

def generate_dt_sizes(dt_min, dt_max, quotient=2):
    StepSizeSpace = [dt_min]
    while StepSizeSpace[-1] < dt_max:
        StepSizeSpace.append(quotient*StepSizeSpace[-1])
    StepSizeSpace.append(dt_max)
    return StepSizeSpace

def eval_probes(sim, initial=False):
    sim.T_avg_s = (sim.Avg_s*sim.T.vector()).sum()  # avg temp Γ_s
    sim.T_avg_r = (sim.Avg_r*sim.T.vector()).sum()  # avg temp Γ_r
    sim.T_avg_a = (sim.Avg_a*sim.T.vector()).sum()  # avg temp Γ_a
    sim.T_ws_value = sim.T_ws(sim.t)                # avg temp Supply pipe
    sim.T_wr_value = sim.T_wr(sim.t)                # avg temp Return pipe
    sim.T_a_value = sim.T_a(sim.t)                  # avg temp Air
    sim.Q_s = sim.alpha_s*sim.Area_s*(sim.T_avg_s-sim.T_ws_value)  # heat flow Γ_s
    sim.Q_r = sim.alpha_r*sim.Area_r*(sim.T_avg_r-sim.T_wr_value)  # heat flow Γ_r
    sim.Q_a = sim.alpha_a*sim.Area_a*(sim.T_avg_a-sim.T_a_value)   # heat flow Γ_a
    sim.Q_dbcs = (sim.flux_dbcs*sim.T.vector()).sum()
    if not initial:
        Energy_prev = (sim.M*sim.T_prev.vector()).sum()
        Energy = (sim.M*sim.T.vector()).sum() 
        sim.E_rate = (Energy - Energy_prev)/sim.dt
        sim.E_error = abs(sim.E_rate+sim.Q_s+sim.Q_r+sim.Q_a+sim.Q_dbcs)

# Default solver definition
default_solver = KrylovSolver("gmres", "hypre_amg")
default_solver.parameters["absolute_tolerance"] = 1E-16
default_solver.parameters["relative_tolerance"] = 1E-12
default_solver.parameters["maximum_iterations"] = 1000
default_solver.parameters["nonzero_initial_guess"] = True

# Default callback
default_callback = Sim_Callback()

def simulate(sim, solver=default_solver, theta=0.5, dt_min=0.001, dt_max=100, t_start=0.0,
             t_stop=None, T_init=5, callback=default_callback):
     
    sim.t = t_start
    if t_stop is not None:
        sim.t_stop = t_stop
    else:
        sim.t_stop = sim.idata.t_stop
    if T_init == "steady-state":
        # Solve initial Steady-State problem
        A = sim.BK
        b = sim.T_ws(sim.t)*sim.B_s + sim.T_wr(sim.t)*sim.B_r + sim.T_a(sim.t)*sim.B_a
        for bc in sim.bcs:
            bc.apply(A, b)
        solver.solve(A, sim.T.vector(), b)
    else:
        sim.T.interpolate(Constant(T_init))

    def accept_solution():
        if callback is not None:
            sim.Elapsed = time.time() - start_time
            callback(sim)
        sim.T_prev.assign(sim.T)
        sim.t_prev = sim.t

    eval_probes(sim, initial=True)
    sim.T_prev.assign(sim.T)
    if callback is not None:
        callback.FirstCall(sim)

    # Solve TimeDependent problem
    A = []
    dt_sizes = generate_dt_sizes(dt_min, dt_max)
    for dt in dt_sizes:
        A.append(sim.M + theta*dt*sim.BK)
        for bc in sim.bcs:
            bc.apply(A[-1])

    start_time = time.time()
    sim.t_prev = t_start
    dtIDX = 0
    STOP = False
    while not STOP:
        sim.dt = dt_sizes[dtIDX]  # select dt according to adapation startegy
        sim.t = sim.t_prev + sim.dt     # new simulation time in k+1 step

        # Calculate new b
        b = sim.M*sim.T_prev.vector() \
            - sim.BK*sim.T_prev.vector()*(1-theta)*sim.dt \
            + (sim.T_ws(sim.t_prev)*sim.B_s + sim.T_wr(sim.t_prev)*sim.B_r \
               + sim.T_a(sim.t_prev)*sim.B_a)*(1-theta)*sim.dt \
            + (sim.T_ws(sim.t)*sim.B_s + sim.T_wr(sim.t)*sim.B_r \
               + sim.T_a(sim.t)*sim.B_a)*theta*sim.dt
        for bc in sim.bcs:
            bc.apply(b)

        solver.solve(A[dtIDX], sim.T.vector(), b)  # solve new T
        eval_probes(sim)

        # Adaptation strategy
        dtIDX_prev = dtIDX
        if sim.E_error > 1:
            dtIDX = max(0, dtIDX-1)
            if dtIDX_prev != dtIDX:
                continue
            else:
                accept_solution()
        elif sim.E_error < 0.5:
            dtIDX = min(len(dt_sizes)-1, dtIDX+1)
            accept_solution()
        else:
            accept_solution()

        # check termination condition
        STOP = sim.t >= sim.t_stop

    if callback is not None:
        callback.close(sim)
