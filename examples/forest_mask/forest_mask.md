# Forest Mask Calculation Workflow

The following workflow was used to calculate the forest mask over the European Alps using ....

The use of the DLT dataset (Copernicus Land Monitoring Service Dominany Leaf Type) for the forest mask is specific to modelling using the CLMS data and a CHM calculating using Sentinel-2 images and the machine learning algorithm in Jiang et al., 2023. 

The use of `scipy.ndimage.binary_dilation()` is applicable to any forest mask for running `CanRad`. 

## Step 1: DLT Data Correction
- Corrected DLT data for where no data in CHM
- **Output:** `DLT_2018_010m_03035_V2_OSHD_2056_corrCHM.tif`

## Step 2: DLT Binarisation
- DLT binarised to 0/1 (using ArcGIS Raster Calculator)
- Resampled from 10m to 5m with nearest neighbour (using ArcGIS Resample)
- **Output:** `OSHD_nonCH_ForestMask_5m_corrCHM_raw.tif`

## Step 3: Buffer Added to Forest Mask
[`1_create_buffer.py`](1_create_buffer.py)
- Uses `scipy.ndimage.binary_dilation()` to add a buffer to the forest mask
- In the case of OSHD nonCH , used a dilation: 30m (6 × 5m pixels)
- **Output:** `OSHD_nonCH_ForestMask_5m_corrCHM_raw_40mbuffer.tif`

## Step 4: Aggregate from 5m to 25m
- Uses script: [`2_aggregate_binary_mask.py`](2_aggregate_binary_mask.py)
- Takes the buffered binary forest mask and aggregates it to 25m for CanRad/Canopy model resolution
- **Output:** `OSHD_nonCH_ForestMask_5m_corrCHM_raw_40mbuffer_agg25m.tif`

> **Note:** Above python scripts were run on Clare's UZH ScienceCloud VM in a GIS environment

## Step 5: Clip to OSHD Domain + set nan to 0
- Clip to `oshd_domain_poly.shp`
- DLT data has larger extent to the OSHD domain -> this doesn't matter because CanRad/CanopyParameters only load relevant data; but the forest mask should match the domain extent because CanRad/CP uses the mask to know where to run the calculations.
- the binary dilation and aggregation steps require that the forest mask is binary with 1 and nan, but for CanRad, 1 and 0 are better because it allows distinction with no forest (0) and out-of-domain (nan). 
- **Output:** `OSHD_nonCH_ForestMask_25m_corrCHM.tif`
