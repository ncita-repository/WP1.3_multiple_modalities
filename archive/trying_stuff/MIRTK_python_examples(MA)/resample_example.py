import sys

# Add MIRTK path to Python environment
mirtk_path = 'D:/Dev/tests/mirtk-v2-201912/Lib/Python' # use forward slashs
sys.path.append(mirtk_path)
import mirtk

# Set verbose switch
mritk_verbose = 1

# Filenames for input and output
input_filename = 'D:/Data/test/mirtk/source.nii'
output_filename = 'D:/Data/test/mirtk/source_voxilised.nii'

# Run MIRTK resample_image command to voxilise the source image
# fill the command arguments
mritk_argv_resample = []
mritk_argv_resample.append('resample-image')
mritk_argv_resample.append(input_filename)
mritk_argv_resample.append(output_filename)
mritk_argv_resample.append('-size')
mritk_argv_resample.append('1')
mritk_argv_resample.append('1')
mritk_argv_resample.append('1')
# run the command
mritk_exit_code = mirtk.call(mritk_argv_resample, verbose=mritk_verbose)
# check the exit code
if mritk_exit_code != 0:
    sys.stderr.write('Error: command returned non-zero exit status ' + str(mritk_exit_code) + '\n')