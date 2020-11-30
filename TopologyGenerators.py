from fenics import *
from mshr import *
import numpy
import datetime
import os

class PipeData:
    def __init__(self, x0, y0, d,
                 name="Default Pipe Name",
                 t=[0.02, 0.1, 0.01],
                 rho=[7850, 60, 960],
                 cp=[525, 5, 1900],
                 lmbd=[50, 0.026, 0.43]):

        self.x0 = x0                        # x - coordinate of pipe center
        self.y0 = y0                        # y - coordinate of pipe center
        self.r = d/2                        # inner radius of pipe
        self.name = name                    # name of pipe
        self.t = t                          # thickness of pipe layers
        self.rho = rho                      # mass density of pipe layers
        self.cp = cp                        # heat capacity of pipe layers
        self.lmbd = lmbd                    # heat conductivity of pipe layers

class TopologyData():
    def __init__(self, mesh, domains, boundaries, DirichletList, rho, cp, lmbd, LogString):
        if MPI.rank(MPI.comm_world) == 0:
            print(LogString)
        self.mesh = mesh                    # Mesh
        self.domains = domains              # Mesh function dim()
        self.boundaries = boundaries        # Mesh function dim() - 1
        self.rho = rho                      # Mesh function DG0
        self.DirichletList = DirichletList  # List of tuples (value, boundary index)
        self.cp = cp                        # Function DG0
        self.lmbd = lmbd                    # Fesh function DG0
        self.LogString = LogString          # Text based info about topology

    def write_pvd_domains(self,path):
        f = File(path)
        f << self.domains

def saveTopologyDATA(td, DirName):
    File(f"{DirName}/mesh.xml") << td. mesh
    File(f"{DirName}/domains.xml") << td.domains
    File(f"{DirName}/boundaries.xml") << td.boundaries
    File(f"{DirName}/rho.xml") << td.rho
    File(f"{DirName}/cp.xml") << td.cp
    File(f"{DirName}/lmbd.xml") << td.lmbd
    numpy.savetxt(f"{DirName}/Dirichlet_list.csv", td.DirichletList, delimiter=',')

    # LogFile
    td.LogString += f"Topology data saved at: {datetime.datetime.now()}"
    Log_file = open(f"{DirName}/log.txt", "w")
    Log_file.write(td.LogString)
    Log_file.close()
    print("Topology data saved succesfully.")

def loadTopologyDATA(DirName):
    mesh = Mesh(f"{DirName}/mesh.xml")
    V0 = FunctionSpace(mesh, 'DG', 0)

    if os.stat(f"{DirName}/Dirichlet_list.csv").st_size == 0:
        DirichletList = []
    else:
        DirichletList = numpy.loadtxt(f"{DirName}/Dirichlet_list.csv", delimiter=',')

    # LogFile
    Log_file = open(f"{DirName}/log.txt", "r")
    LogString = Log_file.read()
    Log_file.close()

    return TopologyData(mesh,
                        MeshFunction("size_t", mesh, f"{DirName}/domains.xml"),
                        MeshFunction("size_t", mesh, f"{DirName}/boundaries.xml"),
                        DirichletList,
                        Function(V0, f"{DirName}/rho.xml"),
                        Function(V0, f"{DirName}/cp.xml"),
                        Function(V0, f"{DirName}/lmbd.xml"),
                        LogString)

# Anulus domain definition
def Anulus(x0, y0, r1, r2, s=100):
    return Circle(Point(x0, y0), r2, s) - Circle(Point(x0, y0), r1, s)

# ------------------------- TOPOLOGY GENERATOR -----------------------------
def TwoBuriedPipes(Pipes, Resolution=50, SoilProperty={"rho": 1300, "cp": 700, "lmbd": 0.34}):
    # HARDcoded soil dimesions
    Xmax = 4
    Xmin = -Xmax
    Ymin = -8

    LogString = f"Topology Instantiated at: {datetime.datetime.now()}\n"
    LogString += f"Soil: generated \n"
    Domain = Rectangle(Point(Xmin, Ymin), Point(Xmax, 0))
    rhoDATA = [SoilProperty["rho"]]
    cpDATA = [SoilProperty["cp"]]
    lmbdDATA = [SoilProperty["lmbd"]]

    # Cut holes representing fluid regions of the pipes
    for P in Pipes:
        LogString = LogString + f"{P.name} Interior: region subtracted \n"
        Domain = Domain - Circle(Point(P.x0, P.y0), P.r, 100)

    subDomain_idx = 0
    for P in Pipes:
        LogString += f"{P.name} Layers:\n"
        for j in range(len(P.t)):
            subDomain_idx += 1
            LogString += f"\tLayer {j}: domain_parts = {subDomain_idx}\n"
            Domain.set_subdomain(subDomain_idx, Anulus(P.x0, P.y0, P.r + sum(P.t[:j]), P.r + sum(P.t[:j+1])))
            rhoDATA.append(P.rho[j])
            cpDATA.append(P.rho[j])
            lmbdDATA.append(P.lmbd[j])

    mesh = generate_mesh(Domain, Resolution)
    LogString += f"Mesh: generated with \n\t{mesh.num_entities(2)} elements \n\t{mesh.num_entities(0)} nodes\n"
    domain_parts = MeshFunction('size_t', mesh, mesh.topology().dim(), mesh.domains())
    boundary_parts = MeshFunction('size_t', mesh, mesh.topology().dim()-1)
    boundary_parts.set_all(99)  # set all away from the ints it might actualy use

    # Function space and functions for Constant material properties
    V0 = FunctionSpace(mesh, 'DG', 0)
    rho = Function(V0)  # Spatial distribution of mass density
    cp = Function(V0)  # Spatial distribution of heat capacity
    lmbd = Function(V0)  # Spatial distribution of heat conductivity
    # choose material property based on subdomains into which individual elements belong
    help = numpy.asarray(domain_parts.array(), dtype=numpy.int32)
    rho.vector()[:] = numpy.choose(help, rhoDATA)
    cp.vector()[:] = numpy.choose(help, cpDATA)
    lmbd.vector()[:] = numpy.choose(help, lmbdDATA)

    # Mark all SubBoudaries
    SubIDX = -1
    for P in Pipes:
        # Mark inner pipe wall surfaces
        SubIDX += 1
        LogString += f"{P.name} Inner Wall: boundary_parts = {SubIDX}\n"
        InnerSurface = CompiledSubDomain(
                'on_boundary && pow(x[0]-x0, 2) + pow(x[1]-y0, 2) <= pow(r+tol, 2)', 
                x0=P.x0, y0=P.y0, r=P.r, tol=1E-14)
        InnerSurface.mark(boundary_parts, SubIDX)

    # Mark Top Boundary of the Soil
    OuterBoundary = CompiledSubDomain('on_boundary && near(x[1], y0)', y0=0)
    SubIDX += 1
    LogString += f"Rectangle TOP: boundary_parts = {SubIDX}\n"
    OuterBoundary.mark(boundary_parts, SubIDX)

    # Mark Bottom Boundary of the Soil
    OuterBoundary = CompiledSubDomain('on_boundary && near(x[1], y0)', y0=Ymin)
    SubIDX += 1
    LogString += f"Rectangle BOTTOM: boundary_parts = {SubIDX}\n"
    OuterBoundary.mark(boundary_parts, SubIDX)

    # Mark Left Boundary
    OuterBoundary = CompiledSubDomain('on_boundary && near(x[0], x0)', x0=Xmin)
    SubIDX += 1
    LogString += f"Rectangle LEFT: boundary_parts = {SubIDX}\n"
    OuterBoundary.mark(boundary_parts, SubIDX)

    # Mark Right Boundary
    OuterBoundary = CompiledSubDomain('on_boundary && near(x[0], x0)', x0=Xmax)
    SubIDX += 1
    LogString += f"Rectangle RIGHT: boundary_parts = {SubIDX}\n"
    OuterBoundary.mark(boundary_parts, SubIDX)

    # List all boundary dirichlet values with boundaries where the should apply
    DirichletList = numpy.array([5, 3])
    LogString += f"DirichletList: = {str(DirichletList)}\n"

    # Return Topology as TopologyData class
    return TopologyData(mesh, domain_parts, boundary_parts, DirichletList, 
                        rho, cp, lmbd, LogString)
