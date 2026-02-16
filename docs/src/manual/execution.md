

# !!! WORKING VERSION !!! DO NOT USE!!!





# Running CanRad

This guide covers different execution modes for running CanRad models (L2R, C2R, or T2R). The execution approach is the same regardless of which model you use.

## Execution Overview

CanRad can be run in three main modes:

| Mode | Best For | Parallelization | Typical Use Case |
|------|----------|-----------------|------------------|
| **Point-scale** | Single location or small domains | None | Testing, local site studies, radiometer comparisons |
| **Distributed** | Medium domains | Julia's `Distributed` | Multi-core processing of tiles |
| **HPC cluster** | Large domains | SLURM/PBS job arrays | > 100 km² scale applications |

## Prerequisites

Before running CanRad, ensure you have:
1. ✅ Preprocessed input data (see [Setup](setup.md))
2. ✅ Configured settings file (see [Settings](settings.md))
3. ✅ Input data paths correctly specified
4. ✅ Output directories created or writable


## Point-Scale Execution

### Single Point or Small Domain

For running CanRad at individual points or over small domains.
See the CanRad/testset folder for example settings files and scripts.

```julia


```

**When to use:**
- Testing settings on a subset
- Single-site studies
- Small domains (<1 km²)
- Development and debugging

### Running from Command Line
Create a simple run script (`run_canrad.jl`):

```julia
using CanRad, DelimitedFiles

# Load settings
include("C2R_settings.jl")

# Run model



```

Then execute:
```bash
julia run_canrad.jl
```

## Distributed Execution (Multi-Core)

For medium-sized domains, use Julia's `Distributed` package to process multiple tiles in parallel on a single machine.

### Tile-Based Processing

First, create a tile list file (`tiles.txt`) with coordinates for each tile:

```
2647000  1185000
2648000  1185000
2647000  1186000
2648000  1186000
```

### Distributed Wrapper Script

Create a wrapper script (`run_distributed.jl`):

```julia
using Pkg
Pkg.activate(".")  # Activate your project environment

using Distributed

# Add worker processes (adjust based on your CPU cores)
addprocs(5)  # Or use: addprocs(Sys.CPU_THREADS)

# Load packages on all workers
@everywhere begin
    using CanRad, DelimitedFiles, SpatialFileIO, NCDatasets
    
    tilesize = 1000      # Tile size in meters
    ptsampling = 25      # Point spacing for calculations
    
    # Read tile coordinates
    tiles = readdlm("tiles.txt")
    
    # Include your tile processing function
    include("tile_canrad.jl")
end

# Process tiles in parallel
pmap for tile_coords in eachrow(tiles)
    process_tile(tile_coords, tilesize, ptsampling)
end

println("All tiles complete!")
```


### Running Distributed Mode

```bash
julia run_distributed.jl
```

**Performance tips:**
- Set `addprocs(N)` where N ≤ physical CPU cores
- Each worker needs sufficient RAM for one tile
- Implement progress tracking to resume failed jobs

**Memory considerations:**
- L2R: ~2-8 GB per tile (depends on point density)
- C2R: ~0.5-2 GB per tile
- T2R: ~0.2-1 GB per tile

## HPC Cluster Execution

For large domains, use job arrays on HPC clusters with SLURM or PBS schedulers.

### Workflow Overview

1. Split domain into tiles
2. Create tile list file
3. Write job submission script
4. Submit job array
5. Monitor progress
6. Collate outputs (if needed)

### SLURM Job Array Script

Create `job_canrad.sh`:

```bash
#!/bin/bash
#SBATCH -J canrad           # Job name
#SBATCH -t 02:00:00         # Time limit (2 hours per tile)
#SBATCH --mem 4000          # Memory per task (4 GB)
#SBATCH --array=1-1000%50   # 1000 tiles, max 50 running simultaneously

# Load required modules (adjust for your cluster)
module load julia/1.9.0

# Path to your Julia script
myjlscript=/path/to/run_canrad_tile.jl

# Calculate tile index (process 4 tiles per job)
JOBNUM=$((((SLURM_ARRAY_TASK_ID-1)*4)+1))

# Run Julia
julia $myjlscript $JOBNUM
```

**Key SLURM directives:**
- `--array=1-1000`: Create 1000 array tasks (one per tile or group of tiles)
- `%50`: Run maximum 50 tasks simultaneously
- `-t 02:00:00`: 2-hour time limit per task
- `--mem 4000`: 4 GB RAM per task

### HPC Julia Run Script

Create `run_canrad_tile.jl`:

```julia
using Pkg
Pkg.activate(".")  # Activate your project environment

using CanRad, DelimitedFiles, SpatialFileIO, NCDatasets

# Get job number from command line argument
jobnum = parse(Int, ARGS[1])

# Configuration
basefolder = "/path/to/project"
datafolder = joinpath(basefolder, "datasets")
settingsfile = joinpath(basefolder, "C2R_settings.jl")
tilesfile = joinpath(basefolder, "tiles.txt")

# Output directories
outdir = joinpath(basefolder, "output_C2R_25m")
progdir = joinpath(basefolder, "progress")
hlmdir = joinpath(basefolder, "output_HLM")

# Create output directories if needed
mkpath(outdir)
mkpath(progdir)
mkpath(hlmdir)

# Read all tiles
all_tiles = readdlm(tilesfile)

# Process multiple tiles per job (e.g., 4 tiles per job)
tiles_per_job = 4
tile_indices = jobnum:min(jobnum+tiles_per_job-1, size(all_tiles, 1))

for idx in tile_indices
    
    xllcorner = Int(all_tiles[idx, 1])
    yllcorner = Int(all_tiles[idx, 2])
    tile_name = "$(xllcorner)_$(yllcorner)"
    
    # Check if already processed
    progfile = joinpath(progdir, "$(tile_name).txt")
    if isfile(progfile)
        println("Tile $(tile_name) already complete, skipping...")
        continue
    end
    
    try
        println("Processing tile $(idx): $(tile_name)")
        
        # Load settings
        include(settingsfile)
        
        # Configure for this tile
        CANRAD.xmin = xllcorner
        CANRAD.xmax = xllcorner + 1000  # 1 km tiles
        CANRAD.ymin = yllcorner
        CANRAD.ymax = yllcorner + 1000
        CANRAD.site = tile_name
        CANRAD.outdir = outdir
        CANRAD.progdir = progdir
        
        # Run model
        @time chm2rad!(CANRAD, SOLAR)
        
        # Mark complete
        open(progfile, "w") do f
            write(f, "Completed at $(now())\n")
            write(f, "Job: $(ENV["SLURM_ARRAY_TASK_ID"])\n")
        end
        
        println("Completed tile $(tile_name)")
        
    catch e
        @error "Error processing tile $(tile_name)" exception=(e, catch_backtrace())
        
        # Write error log
        open(joinpath(progdir, "$(tile_name)_ERROR.txt"), "w") do f
            write(f, "Error: $e\n")
            write(f, "At $(now())\n")
        end
    end
end

println("Job $(jobnum) complete!")
```

### Submitting HPC Jobs

```bash
# Submit job array
sbatch job_canrad.sh

# Check job status
squeue -u $USER

# Check specific job
squeue -j JOBID

# Cancel jobs if needed
scancel JOBID              # Cancel specific job
scancel -u $USER           # Cancel all your jobs
scancel --array=1-100      # Cancel specific array tasks
```

### Using Singularity/Apptainer

Many HPC systems use containers. Example with Singularity:

```bash
#!/bin/bash
#SBATCH -J canrad
#SBATCH -t 02:00:00
#SBATCH --mem 4000
#SBATCH --array=1-1000%50

module load singularityce

# Set Julia depot path (for package caching)
export JULIA_DEPOT_PATH=/scratch/$USER/julia

# Path to Julia container and script
julia_container=/path/to/julia.sif
myjlscript=/path/to/run_canrad_tile.jl

# Calculate job number
JOBNUM=$((((SLURM_ARRAY_TASK_ID-1)*4)+1))

# Run with Singularity
singularity exec $julia_container julia $myjlscript $JOBNUM
```

**First-time setup with containers:**
```bash
# Instantiate packages once before submitting jobs
singularity exec julia.sif julia -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
```

### PBS/Torque Job Script

For PBS-based systems:

```bash
#!/bin/bash
#PBS -N canrad
#PBS -l walltime=02:00:00
#PBS -l mem=4gb
#PBS -t 1-1000%50

cd $PBS_O_WORKDIR

module load julia/1.9.0

# Calculate job number
JOBNUM=$((((PBS_ARRAYID-1)*4)+1))

julia run_canrad_tile.jl $JOBNUM
```

## Progress Monitoring

### Track Completion

Check how many tiles are complete:

```bash
# Count completed tiles
ls output/progress/*.txt | wc -l

# Find failed tiles
ls output/progress/*ERROR.txt 2>/dev/null | wc -l

# Check specific tile
cat output/progress/2647000_1185000.txt
```

### Monitor Resource Usage (SLURM)

```bash
# Check efficiency of completed jobs
seff JOBID

# Monitor running jobs
sstat -j JOBID --format=JobID,MaxRSS,AveCPU

# Check account usage
sacct -u $USER --starttime=2024-01-01 --format=JobID,JobName,Elapsed,State,MaxRSS
```

### Resume Failed Jobs

To rerun only failed tiles, modify the tile list:

```julia
# Create list of incomplete tiles
using DelimitedFiles

all_tiles = readdlm("tiles.txt")
progdir = "output/progress"

incomplete = []
for i in 1:size(all_tiles, 1)
    tile_name = "$(Int(all_tiles[i,1]))_$(Int(all_tiles[i,2]))"
    if !isfile(joinpath(progdir, "$(tile_name).txt"))
        push!(incomplete, all_tiles[i, :])
    end
end

writedlm("tiles_incomplete.txt", incomplete)
println("$(length(incomplete)) tiles need processing")
```

Then resubmit with the incomplete tiles list.

## Output Collation

For tiled processing, you may want to merge outputs into a single file or mosaic.

### Collect Results

```julia
using NCDatasets, DelimitedFiles

tiles = readdlm("tiles.txt")
outdir = "output_C2R_25m"

# Example: Collect sky-view fractions
svf_data = Dict()

for tile in eachrow(tiles)
    tile_name = "$(Int(tile[1]))_$(Int(tile[2]))"
    ncfile = joinpath(outdir, "Output_$(tile_name).nc")
    
    if isfile(ncfile)
        ds = NCDataset(ncfile)
        svf_data[tile_name] = ds["svf_planar_evergreen"][:]
        close(ds)
    end
end

println("Collected data from $(length(svf_data)) tiles")
```

For a complete collation example, see `eurad/eurad25m/collatetiles_eurad.jl` in the repository.

## Performance Optimization

### Resource Allocation

**Memory requirements (per tile):**
- L2R: 2-8 GB (depends on point density)
- C2R: 0.5-2 GB (depends on resolution)
- T2R: 0.2-1 GB

**Time estimates (per 1 km² tile):**
- L2R: 10-60 minutes (depends on complexity)
- C2R: 2-15 minutes
- T2R: 1-5 minutes

**Recommended tile sizes:**
- Workstation: 1-2 km tiles
- HPC: 1 km tiles (good balance of granularity and overhead)

### Optimization Tips

1. **Pre-compute terrain** with T2R if running L2R or C2R multiple times
2. **Use smaller forest_peri** if possible (default 100m)
3. **Disable unnecessary outputs** (save_images, make_pngs)
4. **Process time subsets** for testing before full runs
5. **Use compressed formats** (.laz not .las) to reduce I/O
6. **Keep working directory on fast storage** (not network drives)

### Troubleshooting HPC Jobs

**Job fails immediately:**
- Check module loads work: `module load julia`
- Test script interactively: `srun --pty julia run_canrad_tile.jl 1`
- Verify file paths are absolute not relative

**Out of memory errors:**
- Increase `--mem` allocation
- Reduce tile size
- Check for memory leaks in long-running jobs

**Time limit exceeded:**
- Increase `-t` time limit
- Profile code to find bottlenecks
- Process fewer tiles per job

**Missing files errors:**
- Ensure input data accessible from compute nodes
- Use absolute paths
- Check file permissions

## Example Workflows

### Small Study (< 10 km²)

```bash
# Run directly on workstation
julia run_canrad.jl
```

### Medium Study (10-1000 km²)

```bash
# Use distributed processing
julia run_distributed.jl
```

### Large Study (> 1000 km²)

```bash
# Use HPC cluster
sbatch job_canrad.sh

# Monitor progress
watch -n 60 'ls output/progress/*.txt | wc -l'

# After completion, collate results
julia collate_outputs.jl
```

## Additional Resources

- **Distributed computing:** [Julia Distributed documentation](https://docs.julialang.org/en/v1/manual/distributed-computing/)
- **SLURM:** [SLURM documentation](https://slurm.schedmd.com/)
- **Example scripts:** See `eurad/` folder in repository for large-scale implementations

## Next Steps

After running CanRad:
1. Check progress logs for errors
2. Verify outputs exist for all tiles
3. Review output quality (see [Output](output.md))
4. Collate tiles if needed
5. Analyze results!