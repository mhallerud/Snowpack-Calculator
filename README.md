NOTE: This script has been used in a limited number of circumstances and is not well-vetted. It was designed to work with ArcGIS Desktop Version 10.7.1.

# Snowpack-Calculator
Given a set of coordinates and associated dates, extracts SNODAS data and calculates overall and daily snow depth statistics at each set of locations and associated dates. Designed for camera-trap locations but could be applied with any location data, and could be easily adapted for daily calculations of other raster data.

# To Use
1. Install the numpy dependency, if not already installed.
2. Download the ExtractSnowData.py and BaseHeaderFile.hdr files.
3. Format locations and date ranges as seen in the Example_Input_CSV.csv file.
4. Set the user-defined file paths at the top of the script:
 - `temp_dir` = path to an empty folder where intermediate files will be stored (this can be deleted after the script is finished running)
 - `input_csv` = path to your input data CSV, which should be formatted as in the example CSV
 - `output_folder` = folder path where output calculations will be stored
 - `header_file` = path to the BaseHeaderFile.hdr downloaded from this page (this file is necessary for handling the SNODAS data)
5. Set the buffer distance that you want to use to calculate snow depth, in meters. For example, a `buffer_dist = 1000` would calculate snow depth statistics using a 1000-meter radius at each set of coordinates.

# Output
The output is 5 CSVs:
 - DailyMinsPerStation.csv: Minimum snow depth at each station and date.
 - DailyMaxPerStation.csv: Maximum snow depth at each station and date.
 - DailyMeansPerStation.csv: Mean snow depth at each station and date.
 - DailyStdPerStation.csv: Standard deviation of snow depth within the buffered area around each station and date.
 - OverallStatsPerStation.csv: Overall mean, standard deviation, minimum, maximum, and # of days with snow (Days > 0) at each station thoughout the range of dates. Calculations are made as the mean, minimum, and maximum of the daily values, the overall standard deviation is the standard deviation of the daily means, and the # of days with snow is based on the number of days where the mean snow depth is greater than 0.


