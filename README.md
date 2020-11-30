# Efficient integration of machine learning into district heating predictive models
Suplementary materials

Usage:
1) Generate mesh files etc.
```console
:~$ python3 GenerateTopology.py
```

2) Run FEM simulation to get training data
```console
:~$ mpirun -n 40 python3 SimulationTraining.py
```

3) Run FEM simulation to get test data
```console
:~$ mpirun -n 40 python3 SimulationTest.py
```
4) Train the Modelica model to fit training data
```console
:~$ mpirun -n 40 python3 optimize.py
```
