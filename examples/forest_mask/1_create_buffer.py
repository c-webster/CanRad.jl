
"""

Part of workflow to calculate forest mask for running CanRad-c2r over large domains

STEP 1:

Takes a simple binary forest mask and adds a buffer of specified size around forested areas

Uses multi-processing to speed up processing of large datasets by breaking the input raster into blocks

- input: binary forest mask
- calculations: uses binary dilation to apply a buffer of <pixels> around forest pixels
- output: binary forest mask at same resolution as input

# .gitignore: note to CW: run on a SC-canrad machine in the environment gis-env in ~/windows/projects/eurad/eu_data


"""


import os
from concurrent.futures import ProcessPoolExecutor


import numpy as np
import rasterio
from scipy.ndimage import binary_dilation
import os

# ---------------
# User configuration
# ---------------


# Multiprocessing configuration
num_workers = 15 
use_multiprocessing = True

buffer_pixels = 4 # relative to pixel size of input raster
block_height = 1024

bin_mask_fn = "forest_data/ForestMask/calculations/eualps_ForestMask_010m_raw.tif"
buff_bin_mask_fn = "forest_data/ForestMask/calculations/eualps_ForestMask_010m_40mbuffer.tif"


# ---------------
# Utility functions
# ---------------


def process_buffering_block(args):
    """Process a single block for buffering. Called by multiprocessing workers."""
    row_start, num_rows, width, height, raster_path, buffer_pix = args
    with rasterio.open(raster_path) as src:
        # Read with overlap to avoid edge effects at block boundaries
        read_start = max(0, row_start - buffer_pix)
        read_end = min(height, row_start + num_rows + buffer_pix)
        read_rows = read_end - read_start
        window = rasterio.windows.Window(0, read_start, width, read_rows)
        dlt_chunk = src.read(1, window=window)

        forest_mask = np.isin(dlt_chunk, [1, 2])
        buffered_forest_mask = binary_dilation(
            forest_mask,
            iterations=buffer_pix
        )
        buffered_forest_data = buffered_forest_mask.astype(np.uint8)

        # Trim overlap to original block extent
        trim_top = row_start - read_start
        buffered_forest_data = buffered_forest_data[trim_top:trim_top + num_rows, :]

        return row_start, buffered_forest_data


def create_buffered_forest_mask(input_path, output_path):
    """Create buffered forest mask with multiprocessing for block processing."""
    with rasterio.open(input_path) as src:
        profile = src.profile.copy()
        profile.update({"dtype": "uint8", "nodata": 0, "compress": "lzw"})
        height = src.height
        width = src.width
        
        # Create list of blocks to process
        block_list = []
        for row_start in range(0, height, block_height):
            num_rows = min(block_height, height - row_start)
            block_list.append((row_start, num_rows, width, height, input_path, buffer_pixels))
        
        total_blocks = len(block_list)
        print(f"Processing {total_blocks} blocks...")
        
        # Process blocks in parallel
        if use_multiprocessing and len(block_list) > 1:
            with ProcessPoolExecutor(max_workers=num_workers) as executor:
                results = executor.map(process_buffering_block, block_list)
                block_results = []
                for idx, result in enumerate(results, 1):
                    block_results.append(result)
                    print(f"  [{idx}/{total_blocks}] Block processed")
        else:
            # Serial processing
            block_results = []
            for idx, block in enumerate(block_list, 1):
                result = process_buffering_block(block)
                block_results.append(result)
                print(f"  [{idx}/{total_blocks}] Block processed")
        
        # Write all blocks to output
        print("Writing output raster...")
        with rasterio.open(output_path, "w", **profile) as dst:
            for row_start, buffered_data in block_results:
                num_rows = buffered_data.shape[0]
                window = rasterio.windows.Window(0, row_start, width, num_rows)
                dst.write(buffered_data, 1, window=window)


def main():
    print(f"Using multiprocessing: {use_multiprocessing}")
    print(f"Number of workers: {num_workers}")
    print()
    
    if not os.path.exists(buff_bin_mask_fn):
        print("Creating buffered forest mask...")
        create_buffered_forest_mask(bin_mask_fn, buff_bin_mask_fn)
        print("✓ Buffered forest mask created")
    else:
        print("Buffered forest mask already exists. Skipping creation.")


if __name__ == "__main__":
    main()
