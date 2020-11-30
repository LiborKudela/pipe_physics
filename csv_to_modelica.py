import pandas as pd
import numpy as np

def write_lookup_function(name, x_key, y_key):
    code = f"""
function {name}
  input Real x;
  output Real y;
algorithm
  for i in 1:size({x_key}, 1) loop
    if {x_key}[i+1] > x then
      y := {y_key}[i] + (x - {x_key}[i]) * ({y_key}[i] - {y_key}[i + 1]) / ({x_key}[i] - {x_key}[i + 1]);
      break;
    end if;
  end for;
end {name};
"""
    return code

def write_direct_lookup_function(name, x_key, y_key):
    code = f"""
function {name}
  input Real x;
  output Real y;
protected
  Integer idx;
  Real delta_x;
algorithm
  idx := integer(x/dx);
  delta_x := x - idx*dx;
  y := {y_key}[idx+1] + delta_x*({y_key}[idx+2]-{y_key}[idx+1])/dx;
end {name};
"""
    return code

def write_output_connector():
    code = """
connector Output = output Real annotation(Icon(graphics = {Polygon(lineColor = {52, 101, 164},
                                                           fillColor = {144, 222, 236},
                                                           fillPattern = FillPattern.Solid,
                                                           lineThickness = 1,
                                                           points = {{-100, 100},
                                                                     {100, 0},
                                                                     {-100, -100},
                                                                     {-100, 100}})}));
"""
    return code

def write_source_base():
    code = """
model source_base
  Output y annotation(Placement(visible = true,
                      transformation(origin = {90, 0},
                                     extent = {{-10, -10}, {10, 10}}),
                      iconTransformation(origin = {90, 0},
                                         extent = {{-10, -10}, {10, 10}})));
equation
  annotation(Icon(graphics = {Rectangle(fillColor = {238, 238, 236},
                                        fillPattern = FillPattern.Solid,
                                        lineThickness = 1,
                                        extent = {{-80, 80}, {80, -80}}),
                              Text(origin = {-2, -1},
                                   extent = {{-74, 13}, {74, -13}},
                                   textString = "%name")}));
end source_base;
"""
    return code

def write_source(name, function):
    code = f"""
model {name}
  extends source_base;
equation
  y=Functions.{name}(time);
end {name};
"""
    return code

def convert_to_modelica_package(path, dx=None, make_callable=True, output=None):
    code = "package data\n"
    code += f"constant Real dx = {dx};\n" 
    # insert data
    df = pd.read_csv(path)
    keys = df.keys()
    x_key = keys[0]
    y_keys = keys[1:]
    x_min = df[x_key].iloc[0]
    x_max = df[x_key].iloc[-1]
    x_vals = np.arange(x_min, x_max, dx)
    if x_vals[-1] != x_max:
      np.append(x_vals,[x_max])
    for key in keys:
      vals = np.interp(x_vals, df[x_key], df[key]).astype(str)
      code += f"constant Real {key}_data[:] = " + "{" + ",".join(vals) + "};\n" 

    #insert 1D interpolations functions
    code += "package Functions\n"
    for y_key in y_keys:
        code += write_direct_lookup_function(y_key, f"{x_key}_data", f"{y_key}_data") 
    code += "end Functions;\n"

    # insert data sources blocks
    code += write_output_connector()
    code += write_source_base()
    for y_key in y_keys:
        code += write_source(y_key, f"Functions.{y_key}")

    code += "end data;"

    # save modelica file to disk
    if output is None:
        output = f"{path[:-4]}.mo" 
    f = open(output, "w")
    f.write(code)
    f.close()
