# Fynbos-Fire-Recovery-and-Vegetation-Recovery-in-Jonkershoek-Nature-Reserve
# Post-fire NDVI recovery, Jonkershoek (Feb 2021 fire)

Biomathematics final project.

Data: Sentinel-2 Level-2A surface reflectance (Harmonized collection,
`COPERNICUS/S2_SR_HARMONIZED`), accessed through Google Earth Engine.

## Files
- `ndvi_recovery.js` – Google Earth Engine script that extracts monthly NDVI
- `Jonkershoek_NDVI_quarterly.csv` – NDVI data exported from Earth Engine
- `ndvi_models.R` – fits the linear, exponential and seasonal recovery models

## Running the Earth Engine script
1. Open the `.js` file and click **Copy raw file** (top right).
2. Go to https://code.earthengine.google.com and sign in (an Earth Engine account is needed).
3. Paste the code into the editor and click **Run**.

## Running the R script
1. Download the `.R` file and the CSV into the same folder.
2. Open the `.R` file in RStudio and set the working directory to that folder
   (Session → Set Working Directory → To Source File Location).
3. Install packages if needed: `install.packages(c("tidyverse", "patchwork"))`
4. Run the script.
