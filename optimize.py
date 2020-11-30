import lofi

model_file = 'Pipes.mo'
data_file = 'dn250_ic1/training_data.mo' 
model_name = 'Pipes.Tests.dn250_ic1_unsolved'

# create the API, JIT-compiles model etc.
API = lofi.APIs.open_modelica
model = API([model_file, data_file], model_name, visual_callback=None)

# instantiate optimizer and attach model to it
opt = lofi.optimizers.GRAPSO(model)

# sets some options of the optimizer and refresh them
opt.options.n = 39
opt.restart()

opt.train(500)
lofi.imode(globals())

# if the process stalls, restart the optimizer
# > opt.restart()
# > opt.train(200)
