
"""

Part of workflow to calculate forest mask for running CanRad-c2r over large domains

STEP 2:

Takes the buffered binary forest mask and aggregates it to desired coarser resolutions

Uses multi-processing to speed up processing of large datasets by breaking the input raster into blocks

- input: binary buffered forest mask
- calculations: aggregates the buffered forest mask to coarser resolutions using maximum value within aggregation window
- output: aggregated binary forest mask at coarser resolutions

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

bin_mask_fn = "forest_data/ForestMask/calculations/eualps_ForestMask_010m_40mbuffer.tif"

target_resolutions = [25, 100, 250, 1000]  # in meters

# ---------------
# Utility functions
# ---------------   
def format_resolution(res):
    num_digits = len(str(int(max(target_resolutions))))
    return f"{int(res):0{num_digits}d}m"

def common_resolution(input_res_m, target_res_m, scale=100):  # scale=100 for cm
    a = round(input_res_m * scale)
    b = round(target_res_m * scale)
    return gcd(a, b) / scale

def process_resample_block(args):
    """Process a single block for sub-sampling. Called by multiprocessing workers."""
    data_chunk, row_start, num_rows, scale_factor = args
    
    # Sub-sample using nearest neighbor repetition
    # Each pixel is repeated scale_factor times in each dimension
    resampled_chunk = np.repeat(np.repeat(data_chunk, scale_factor, axis=0), scale_factor, axis=1)
    
    return row_start * scale_factor, resampled_chunk

def resample_raster(src, new_res, output_raster_fn):
    res_str = format_resolution(new_res)
    input_resolution = src.res[0]
    scale_factor = int(input_resolution // new_res)
    
    print(f"Resampling to {res_str}...")
    
    # Read entire raster
    data = src.read(1)
    
    # Create list of blocks to process
    block_list = []
    for row_start in range(0, data.shape[0], block_height):
        num_rows = min(block_height, data.shape[0] - row_start)
        data_chunk = data[row_start:row_start + num_rows, :]
        block_list.append((data_chunk, row_start, num_rows, scale_factor))
    
    total_blocks = len(block_list)
    print(f"Processing {total_blocks} blocks...")
    
    # Process blocks in parallel
    if use_multiprocessing and len(block_list) > 1:
        with ProcessPoolExecutor(max_workers=num_workers) as executor:
            results = executor.map(process_resample_block, block_list)
            block_results = []
            for idx, result in enumerate(results, 1):
                block_results.append(result)
                print(f"  [{idx}/{total_blocks}] Block resampled")
    else:
        # Serial processing
        block_results = []
        for idx, block in enumerate(block_list, 1):
            result = process_resample_block(block)
            block_results.append(result)
            print(f"  [{idx}/{total_blocks}] Block resampled")
    
    # Calculate new dimensions
    new_height = data.shape[0] * scale_factor
    new_width = data.shape[1] * scale_factor
    
    # Update metadata
    profile = src.profile
    profile.update({
        'height': new_height,
        'width': new_width,
        'crs': src.crs,
        'dtype': 'uint8',
        'transform': src.transform * src.transform.scale(1/scale_factor, 1/scale_factor)
    })
    
    # Sort blocks by row_start to ensure correct writing order
    block_results.sort(key=lambda x: x[0])
    
    # Write all blocks to output
    print("Writing resampled raster...")
    with rasterio.open(output_raster_fn, 'w', **profile) as dst:
        for row_start, resampled_data in block_results:
            num_rows = resampled_data.shape[0]
            window = rasterio.windows.Window(0, row_start, new_width, num_rows)
            dst.write(resampled_data, 1, window=window)
    print(f"  ✓ Resampled raster written to {output_raster_fn}")


def process_aggregate_block(args):
    """Process a single block for aggregation with max. Called by multiprocessing workers."""
    row_start, num_rows, width, raster_path, scale_factor = args
    with rasterio.open(raster_path) as src:
        height = src.height

        # Determine the output rows this block is responsible for
        start_keep = ((row_start + scale_factor - 1) // scale_factor) * scale_factor
        end_keep = ((row_start + num_rows - 1) // scale_factor) * scale_factor

        max_output_rows = height // scale_factor
        if max_output_rows <= 0:
            return None
        last_full_start = (max_output_rows - 1) * scale_factor
        end_keep = min(end_keep, last_full_start)

        if start_keep >= row_start + num_rows or end_keep < start_keep:
            return None

        # Read aligned data range with overlap
        aligned_start = (row_start // scale_factor) * scale_factor
        read_start = aligned_start
        read_end = min(height, end_keep + scale_factor)
        read_rows = read_end - read_start

        window = rasterio.windows.Window(0, read_start, width, read_rows)
        data_chunk = src.read(1, window=window)

        # Ensure dimensions are divisible by scale_factor
        crop_height = (read_rows // scale_factor) * scale_factor
        crop_width = (width // scale_factor) * scale_factor
        data_chunk = data_chunk[:crop_height, :crop_width]

        # Reshape and take max over each scale_factor x scale_factor block
        new_height = crop_height // scale_factor
        new_width = crop_width // scale_factor

        reshaped = data_chunk.reshape(new_height, scale_factor, new_width, scale_factor)
        aggregated_chunk = reshaped.max(axis=1).max(axis=2)

        # Trim to only the output rows that belong to this block
        start_idx = (start_keep - read_start) // scale_factor
        end_idx = (end_keep - read_start) // scale_factor
        aggregated_chunk = aggregated_chunk[start_idx:end_idx + 1, :]

        return start_keep // scale_factor, aggregated_chunk


def aggregate_raster_max(src, output_raster, target_res):
    res_str = format_resolution(target_res)
    input_res = src.res[0]
    scale_factor = int(target_res // input_res)
    
    # Calculate the new shape
    new_height = src.height // scale_factor
    new_width = src.width // scale_factor

    # Create list of blocks to process
    block_list = []
    for row_start in range(0, src.height, block_height):
        num_rows = min(block_height, src.height - row_start)
        block_list.append((row_start, num_rows, src.width, src.name, scale_factor))

    total_blocks = len(block_list)
    print(f"Aggregating to {res_str}: Processing {total_blocks} blocks...")

    # Process blocks in parallel
    if use_multiprocessing and len(block_list) > 1:
        with ProcessPoolExecutor(max_workers=num_workers) as executor:
            results = executor.map(process_aggregate_block, block_list)
            block_results = []
            for idx, result in enumerate(results, 1):
                if result is not None:
                    block_results.append(result)
                print(f"  [{idx}/{total_blocks}] Block aggregated")
    else:
        # Serial processing
        block_results = []
        for idx, block in enumerate(block_list, 1):
            result = process_aggregate_block(block)
            if result is not None:
                block_results.append(result)
            print(f"  [{idx}/{total_blocks}] Block aggregated")

    # Update metadata
    profile = src.profile
    profile.update({
        'height': new_height,
        'width': new_width,
        'crs': src.crs,
        'dtype': 'uint8',
        'transform': src.transform * src.transform.scale(scale_factor, scale_factor)
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
    print(f"  ✓ Aggregated raster written to {output_raster}")

def main():
    for target_res in target_resolutions:
        target_res_str = format_resolution(target_res)
        with rasterio.open(bin_mask_fn) as src: 
            input_res = src.res[0]
            # check can achieve desired resolutions
            if int(target_res // input_res) != (target_res / input_res):
                print("cannot achieve desired resolution with integer scale factor from input resolution")
                print("resampling raster to a common resolution for output raster")
                common_res = common_resolution(input_res,target_res)
                common_res_str = format_resolution(common_res)
                output_raster_fn = os.path.join(os.path.dirname(bin_mask_fn),
                                    os.path.basename(bin_mask_fn[:-4]+f"_resample{common_res_str}.tif"))
                if not os.path.exists(output_raster_fn):
                    resample_raster(src,common_res, output_raster_fn)
                    print(f"  ✓ Resampled raster to {common_res_str}")
                else:
                    print(f"  ✓ Resampled raster to {common_res_str} already exists")
                input_raster_fn = output_raster_fn
            else:
                input_raster_fn = bin_mask_fn

        with rasterio.open(input_raster_fn) as src: 
            input_res = src.res[0]
            print(f"aggregating {input_res}m forest mask to {target_res_str} resolution")
            aggregated_raster_fn = os.path.join(os.path.dirname(input_raster_fn),
                        os.path.basename(bin_mask_fn[:-4]+f"_agg{target_res_str}.tif"))
            aggregate_raster_max(src,aggregated_raster_fn,target_res)
            print(f"  ✓ Aggregated forest mask to {target_res_str}")




if __name__ == '__main__':
    main()