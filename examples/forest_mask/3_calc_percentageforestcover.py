
"""

Part of workflow to calculate forest mask for running CanRad-c2r over large domains

STEP 3:

Takes the buffered binary forest mask at canrad resolution and 
    calculates percentage forest cover at coarser resolutions for (e.g.) fsm-shd grids

Uses multi-processing to speed up processing of large datasets by breaking the input raster into blocks

- input: binary buffered forest mask
- calculations: calculates percentage forest cover at coarser resolutions for (e.g.) fsm-shd grids
- output: percentage forest cover raster at coarser resolutions

# .gitignore: note to CW: run on a SC-canrad machine in the environment gis-env in ~/windows/projects/eurad/eu_data


"""


import os
from concurrent.futures import ProcessPoolExecutor
import numpy as np

import rasterio
from rasterio.enums import Resampling
from rasterio.mask import mask
from math import gcd


# ---------------
# User configuration
# ---------------
# Multiprocessing configuration
use_multiprocessing = True
num_workers = 15
block_height = 2048

bin_mask_fn = "forest_data/ForestMask/calculations/eualps_ForestMask_010m_40mbuffer_agg025m.tif"

target_resolutions = [100, 250, 1000]  # in meters


# ---------------
# Utility functions
# ---------------


def format_resolution(res):
    num_digits = len(str(int(max(target_resolutions))))
    return f"{int(res):0{num_digits}d}m"

def process_percentage_block(args):
    """Process a single block for percentage calculation. Called by multiprocessing workers."""
    data_chunk, row_start, num_rows, scale_factor, trim_top = args
    
    # Calculate dimensions for aggregated block
    new_height = num_rows // scale_factor
    new_width = data_chunk.shape[1] // scale_factor
    
    # Crop to be divisible by scale_factor
    crop_height = new_height * scale_factor
    crop_width = new_width * scale_factor
    data_chunk = data_chunk[:crop_height, :crop_width]
    
    # Initialize output array for this block
    aggregated_chunk = np.zeros((new_height, new_width), dtype=np.float32)
    
    # Calculate percentage for each aggregated cell
    for i in range(new_height):
        for j in range(new_width):
            cell_data = data_chunk[
                i*scale_factor:(i+1)*scale_factor,
                j*scale_factor:(j+1)*scale_factor
            ]
            forest_cover = np.sum(cell_data == 1)
            total_pixels = cell_data.size
            aggregated_chunk[i, j] = round((forest_cover / total_pixels) * 100)
    
    # Trim overlapping rows from the top (if any)
    if trim_top > 0:
        trim_rows = trim_top // scale_factor
        aggregated_chunk = aggregated_chunk[trim_rows:, :]

    return row_start // scale_factor, aggregated_chunk

def aggregate_raster_percentage(src, output_raster, scale_factor):
    """Aggregate raster with percentage calculation using multiprocessing."""
    # Calculate the new shape
    new_height = src.height // scale_factor
    new_width = src.width // scale_factor
    
    print("Reading raster...")
    # Read entire raster
    data = src.read(1)
    
    # Create list of blocks to process
    block_list = []
    overlap = scale_factor - 1
    for row_start in range(0, data.shape[0], block_height):
        num_rows = min(block_height, data.shape[0] - row_start)
        # Ensure num_rows is divisible by scale_factor
        num_rows = (num_rows // scale_factor) * scale_factor
        if num_rows == 0:
            continue

        # Include overlap above (except first block)
        trim_top = 0
        read_start = row_start
        read_rows = num_rows
        if row_start > 0:
            trim_top = overlap
            read_start = row_start - overlap
            read_rows = num_rows + overlap

        # Ensure bounds
        read_rows = min(read_rows, data.shape[0] - read_start)

        data_chunk = data[read_start:read_start + read_rows, :]
        block_list.append((data_chunk, row_start, num_rows, scale_factor, trim_top))
    
    total_blocks = len(block_list)
    print(f"Processing {total_blocks} blocks...")
    
    # Process blocks in parallel
    if use_multiprocessing and len(block_list) > 1:
        with ProcessPoolExecutor(max_workers=num_workers) as executor:
            results = executor.map(process_percentage_block, block_list)
            block_results = []
            for idx, result in enumerate(results, 1):
                block_results.append(result)
                print(f"  [{idx}/{total_blocks}] Block processed")
    else:
        # Serial processing
        block_results = []
        for idx, block in enumerate(block_list, 1):
            result = process_percentage_block(block)
            block_results.append(result)
            print(f"  [{idx}/{total_blocks}] Block processed")
    
    # Scale the transform
    transform = src.transform * src.transform.scale(
        (src.width / new_width),
        (src.height / new_height)
    )
    
    # Update metadata
    profile = src.profile
    profile.update({
        'height': new_height,
        'width': new_width,
        'transform': transform,
        'dtype': 'float32',
        'count': 1,
        'crs': src.crs
    })
    
    # Sort blocks by row_start to ensure correct writing order
    block_results.sort(key=lambda x: x[0])
    
    # Write the aggregated data to a new file
    print("Writing aggregated raster...")
    with rasterio.open(output_raster, 'w', **profile) as dst:
        for row_start, aggregated_data in block_results:
            num_rows = aggregated_data.shape[0]
            window = rasterio.windows.Window(0, row_start, new_width, num_rows)
            dst.write(aggregated_data, 1, window=window)



def main():
    for target_res in target_resolutions:
        target_res_str = format_resolution(target_res)
        with rasterio.open(bin_mask_fn) as src:
            input_resolution = src.res[0]  # Assuming square pixels
            output_raster_fn = os.path.join(os.path.dirname(bin_mask_fn),f"ForestCoverFraction_{target_res_str}.tif")
            print(f"Calculating % forest cover at {target_res_str}...")
            scale_factor_res = int(target_res // input_resolution)
            aggregate_raster_percentage(
                src,
                output_raster_fn,
                scale_factor_res
            )
            print(f"  ✓ Calculated % forest cover at {target_res_str}")

if __name__ == '__main__':
    main()