import sys, os, time

# Add MIRTK path to Python environment
mirtk_path = 'D:/Dev/tests/mirtk-v2-201912/Lib/Python' # use forward slashs
sys.path.append(mirtk_path)
import mirtk

# Set verbose switch
mritk_verbose = 1

# Filenames for input and output
input_filename = 'D:/Data/test/mirtk/source.nii'
target_filename = 'D:/Data/test/mirtk/target.nii'
# the output transformation estimate : matrix for affine or 3D(2D) vector map for FFD
dofout_filename = 'D:/Data/test/mirtk/trans_estimate.dof.gz'

# Registration models & config files
model_affine = 'Rigid+Affine'
model_ffd = 'Rigid+Affine+FFD'

config_filename = 'ireg.cfg'
config_filename = os.path.join(os.path.dirname(__file__), 'cfg', config_filename)

# Run MIRTK register command
start_time = time.time()
# fill the command arguments
mritk_argv_register = []
mritk_argv_register.append('register')
mritk_argv_register.append(target_filename)
mritk_argv_register.append(input_filename)
mritk_argv_register.append('-model')
mritk_argv_register.append(model_affine)
mritk_argv_register.append('-parin')
mritk_argv_register.append(config_filename)
mritk_argv_register.append('-dofout')
mritk_argv_register.append(dofout_filename)
# run the command
mritk_exit_code = mirtk.call(mritk_argv_register, verbose=mritk_verbose)
# check the exit code
if mritk_exit_code != 0:
    sys.stderr.write('Error: command returned non-zero exit status ' + str(mritk_exit_code) + '\n')

print ("Registration took %s seconds" %(time.time() - start_time))

# Transform the input image on the target space using the transformation estimate
reg_output_filename = 'D:/Data/test/mirtk/reg_output.nii'
mritk_argv_transform = []
mritk_argv_transform.append('transform-image')
mritk_argv_transform.append(input_filename)  # source
mritk_argv_transform.append(reg_output_filename)  # output
if dofout_filename:
    mritk_argv_transform.append('-dofin')
    mritk_argv_transform.append(dofout_filename)
mritk_argv_transform.append('-target')
mritk_argv_transform.append(target_filename)
mritk_exit_code = mirtk.call(mritk_argv_transform, verbose=mritk_verbose)
if mritk_exit_code != 0:
    sys.stderr.write('Error: command returned non-zero exit status ' + str(mritk_exit_code) + '\n')