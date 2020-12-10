# Efficient integration of machine learning into district heating predictive models
---
### Published Article
DOI:  10.3390/en13236381
Link: https://www.mdpi.com/1996-1073/13/23/6381

---
Dependencies:
- FeniCS: https://fenicsproject.org/
- OpenModelica: https://www.openmodelica.org/
- lofi: https://github.com/LiborKudela/lofi

Usage:

1) Generate mesh files for the complex simulation etc.
```console
python3 GenerateTopology.py
```

2) Run FEM simulation to generate data for training
```console
mpirun -n 40 python3 SimulationTraining.py
```

3) Run FEM simulation to generate data for testing
```console
mpirun -n 40 python3 SimulationTest.py
```
4) Train the Modelica model to fit training data in a interactive mode
```console
mpirun -n 40 python3 optimize.py
```
