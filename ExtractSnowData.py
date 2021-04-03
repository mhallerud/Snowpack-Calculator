#-------------------------------------------------------------------------------------#
# Title: Snow Depth Data Extraction                                                   #
# Description: Given a set of coordinates and associated dates, extracts SNODAS data  #
#            and calculates overall and daily snow depth statistics for the           #
#            coordinates at the dates provided                                        #
#                                                                                     #
# Compatible with ArcGIS Desktop Version 10.7.1                                       #
#                                                                                     #
# Date Created: Aug 6 2020                                                            #
# Author: Maggie Hallerud (margaret.hallerud@maine.edu)                               #
#-------------------------------------------------------------------------------------#

## set user inputs
temp_dir = 'C:/Users/Maggie/Desktop/SNODAS Winter Snow Depth Data' # temporary folder that will be made to store temp files
buffer_dist = 1000 # buffer distance in meters
input_csv = 'C:/Users/Maggie/Desktop/Snow_Depth_Calculations/NECESSARY_CODE_FILES/InputCameraData.csv' # input CSV with coordinates and dates
output_folder = 'C:/Users/Maggie/Desktop/Snow_Depth_Calculations/Output_CSVs' # folder to store outputs in
header_file = 'C:/Users/Maggie/Desktop/Snow_Depth_Calculations/NECESSARY_CODE_FILES/BaseHeaderFile.hdr' # filepath to header file with specified raster attributes


# import dependencies
print "Setting up coding environment..."
import os
import csv
import arcpy
import numpy as np
import urllib2
import zipfile
import glob
from collections import defaultdict
import datetime
import tarfile
import gzip
import shutil
import csv


def main():
    # set up environment
    if not os.path.exists(temp_dir):
        os.mkdir(temp_dir)
    arcpy.CheckOutExtension("spatial")
    arcpy.env.overwriteOutput=True

    print "Reading in data..."
    # Make dictionary and list of station IDs from input CSV
    cam_list = []
    camDict = defaultdict(dict)
    with open(input_csv, "r") as infile:
        reader = csv.reader(infile)
        headers = next(reader)[1:]
        for row in reader:
            camDict[row[0]] = {key: value for key, value in zip(headers, row[1:])}
            cam_list.append(row[0])
    
    # Make SHP of camera stations (coords in WGS84)
    camera_lyr_wgs84 = arcpy.MakeXYEventLayer_management(input_csv, "Longitude", "Latitude", 'Camera_Stations_WGS84', arcpy.SpatialReference(4326))
    camera_shp_wgs84 = os.path.join(temp_dir, "Camera_Stations_WGS84.shp")
    if os.path.exists(camera_shp_wgs84):
        arcpy.Delete_management(camera_shp_wgs84)
    arcpy.FeatureClassToShapefile_conversion(camera_lyr_wgs84, temp_dir)
    
    # Set up empty lists for storing output values
    daily_means = []
    daily_mins = []
    daily_max = []
    daily_std = []
    station_stats = []
    
    # Loop through cameras and download/clip data to buffered points
    for cam in cam_list:
        try:
            print "Collecting data for camera station " + cam + "..."       
            # Make folder for camera in temp dir
            cam_folder = os.path.join(temp_dir, cam)
            if not os.path.exists(cam_folder):
                os.mkdir(cam_folder)
            
            # Pull start and end dates for camera from CSV
            str_start_date = camDict[cam]['StartDate']
            str_end_date = camDict[cam]['EndDate']
            start_date = datetime.datetime.strptime(str_start_date, "%m/%d/%Y")
            end_date = datetime.datetime.strptime(str_end_date, "%m/%d/%Y")
            
            print "....Downloading data"
            # Set up folder structure
            download_folder = os.path.join(cam_folder, "0_Downloaded_TAR_Files")
            if os.path.exists(download_folder):
                try:
                    shutil.rmtree(download_folder)
                except Exception as err:
                    pass
            os.mkdir(download_folder)
            i = 0 # set counter
            ndays = int((end_date - start_date).days + 1) # calculate # days 
            # Find and download data for each date associated with camera
            for n in range(ndays):
                i += 1
                date = start_date + datetime.timedelta(n)
                try:
                    # specify dates and number of days
                    i += 1
                    date = start_date + datetime.timedelta(n)
                    day_count = str(i).rjust(3,'0')
                    year = datetime.datetime.strftime(date, "%Y")
                    month = datetime.datetime.strftime(date, "%m")
                    day = datetime.datetime.strftime(date, "%d")
                    month_name = datetime.datetime.strftime(date, "%b")
                    # download data from FTP
                    url_filename = "SNODAS_" + year + month + day + ".tar"
                    download = "ftp://sidads.colorado.edu/DATASETS/NOAA/G02158/masked/" + year + "/" + month + '_' + month_name + "/" + url_filename
                    request = urllib2.urlopen(download)
                    # save data
                    out_filename = "DAY" + day_count + "_SNODAS_" + year + month + day + ".tar"
                    output = open(download_folder + '/' + out_filename, "wb")
                    output.write(request.read())
                    output.close()
                # error handling
                except Exception as err:
                    print "........WARNING: Could not download data for camera station " + cam + " at day " + str(date)
            
            print "....Extracting data"
            # Set up folder structure to hold extracted files
            gz_folder = os.path.join(cam_folder, "1_GZ_Files")
            extract_folder = os.path.join(cam_folder, "2_Extracted_SnowDepth_Files")
            if os.path.exists(gz_folder):
                try:
                    shutil.rmtree(gz_folder)
                except Exception as err:
                    pass
            os.mkdir(gz_folder)
            if os.path.exists(extract_folder):
                try:
                    shutil.rmtree(extract_folder)
                except Exception as err:
                    pass
            os.mkdir(extract_folder)
            
            # Extract GZ files from TAR files
            for tar in glob.glob(os.path.join(download_folder, "*.tar")):
                try:
                    my_tar = tarfile.open(tar)
                    my_tar.extractall(gz_folder)
                    my_tar.close()
                except Exception as err:
                    print "........WARNING: Data could not be extracted from file " + tar
                
            # Unzip GZ files associated with snow depth
            for gz in glob.glob(os.path.join(gz_folder, "*.gz")):
                gz_filepath = os.path.basename(gz)
                if gz_filepath[8:12] == "1036":
                    try:
                        with gzip.open(gz, 'rb') as f:
                            file_content = f.read()
                            f.close()
                            gz_extract_name = gz_filepath.split('.')[0] + '.' + gz_filepath.split('.')[1]
                            gz_extract_path = os.path.join(extract_folder, gz_extract_name)
                            gz_output = open(gz_extract_path, "wb")
                            gz_output.write(file_content)
                            gz_output.close()
                    except Exception as err:
                        print "........WARNING: Could not extract data from file " + gz

            # Add .hdr file and define coordinate system for all extracted snow depth files
            for dat in glob.glob(os.path.join(extract_folder, "*.dat")):
                try:
                    # copy over header file
                    shutil.copy2(header_file, dat.split('.')[0] + '.hdr')
                    # define projection
                    arcpy.DefineProjection_management(dat, arcpy.SpatialReference(4326))
                    # define "No Data" value
                    #nodata_vt = arcpy.ValueTable(2)
                    #nodata_vt.addRow('1 -9999')
                    #arcpy.SetRasterProperties_management(dat, data_type="SCIENTIFIC", nodata_vt=-9999)
                except Exception as err:
                    print "........WARNING: Could not define coordinate system for " + dat

            # Clean up non-extracted files
            print "....Deleting unnecessary files"
            try:
                shutil.rmtree(gz_folder)
                shutil.rmtree(download_folder)
            except Exception as err:
                print "........WARNING: Could not delete unnecessary files"
        except Exception as err:
            print "--------------------------------------------------------------------------"
            print "WARNING: There was a problem processing data for camera " + cam
            print "The exception thrown was:"
            print err
            print "--------------------------------------------------------------------------"

        print "--Starting calculations..."
        try:
            # Pull camera station point
            camera_site = os.path.join(cam_folder, 'Station_' + cam.replace('-', '_') + '.shp')
            arcpy.Select_analysis(camera_shp_wgs84, camera_site, """ "StationID" = '%s'""" % cam)

            # If buffer distance is larger than or equal to 500 meters (half cell resolution), then buffer camera stations and make calculations
            buffer_cam = os.path.join(cam_folder, "BufferedCam_"+str(buffer_dist)+"m.shp")
            if buffer_dist >= 500:
                print "....Clipping rasters to buffered camera point"
                arcpy.Buffer_analysis(camera_site, buffer_cam, "%d Meters" % buffer_dist)
                
                # Set up folder structure
                clipped_folder = os.path.join(cam_folder, "3_Clipped_Rasters")
                if os.path.exists(clipped_folder):
                    try:
                        shutil.rmtree(clipped_folder)
                        os.mkdir(clipped_folder)
                    except Exception as err:
                        pass
                else:
                    os.mkdir(clipped_folder)
                    
                # Clip all rasters to buffered camera station
                for dat in glob.glob(os.path.join(extract_folder, "*.dat")):
                    try:
                        dat_name = os.path.basename(dat).split('.')[0]
                        dat_clip = os.path.join(clipped_folder, dat_name+'_Clip.tif')
                        arcpy.Clip_management(dat, '', dat_clip, buffer_cam, "-9999", 'ClippingGeometry', 'NO_MAINTAIN_EXTENT')
                    except Exception as err:
                        print "........WARNING: There was an error clipping " + dat

                # Calculate stats per clipped raster (daily stats)
                # Set up arrays to save daily values for camera station
                # NOTE: Units are meters
                print "....Calculating daily statistics for station"
                daily_cam_means = [cam]
                daily_cam_mins = [cam]
                daily_cam_max = [cam]
                daily_cam_std = [cam]
                i=0 # set counter
                # Calculate statistics for each raster
                for tif in glob.glob(os.path.join(clipped_folder, "*.tif")):
                    i+=1
                    #day = os.path.basename(tif)[3:6]
                    #if int(day) != i:
                    #    print "WARNING! Day " + str(i) + " column for cam " + cam + " does not match - actual day is " + day
                    try:
                        # Calculate and pull raster stats
                        arcpy.CalculateStatistics_management(tif)
                        cam_day_mean = arcpy.GetRasterProperties_management(tif, "MEAN")[0] 
                        cam_day_min = arcpy.GetRasterProperties_management(tif, "MINIMUM")[0]
                        cam_day_max = arcpy.GetRasterProperties_management(tif, "MAXIMUM")[0]
                        cam_day_std = arcpy.GetRasterProperties_management(tif, "STD")[0]
                        # Add raster stats to list and convert to meters
                        daily_cam_means.append(float(cam_day_mean)/1000)
                        daily_cam_mins.append(float(cam_day_min)/1000)
                        daily_cam_max.append(float(cam_day_max)/1000)
                        daily_cam_std.append(float(cam_day_std)/1000)
                    # If calculations fail, adds "nodata" value for that day instead
                    except Exception as err:
                        print "WARNING! Statistics not calculated for " + tif
                        print "Exception thrown was:"
                        print err
                        daily_cam_means.append("NoData")
                        daily_cam_mins.append("NoData")
                        daily_cam_max.append("NoData")
                        daily_cam_std.append("NoData")

                # Add camera values to overall daily dataset
                daily_means.append(daily_cam_means)
                daily_mins.append(daily_cam_mins)
                daily_max.append(daily_cam_max)
                daily_std.append(daily_cam_std)

            # If buffer distance is less than 500 meters, extract data to exact points
            else:
                print "....Extracting daily values at camera point"
                # Make values folder
                vals_folder = os.path.join(cam_folder, "3_Exact_Values")
                if os.path.exists(vals_folder):
                    try:
                        shutil.rmtree(vals_folder)
                        os.mkdir(vals_folder)
                    except Exception as err:
                        print err
                else:
                    os.mkdir(vals_folder)

                # Emty list to hold camera output for all days
                cam_exact_values = [cam]
                
                # Loop through tifs and extract values to points
                for dat in glob.glob(os.path.join(extract_folder, "*.dat")):
                    # Avoiding errors
                    fields = [f.name for f in arcpy.ListFields(camera_site)]
                    if "RASTERVALU" in fields:
                        arcpy.DeleteField_management(camera_site, "RASTERVALU")
                    # extract values to camera point
                    base_name = os.path.basename(dat).split('.')[0]
                    out_val_pts = os.path.join(vals_folder, "Values_"+base_name+".shp")
                    arcpy.sa.ExtractValuesToPoints(camera_site, dat, out_val_pts)
                    # add to list of exact values
                    with arcpy.da.SearchCursor(out_val_pts, field_names=["RASTERVALU"]) as cursor:
                        for row in cursor:
                            cam_exact_values.append(float(row[0]))

                # Append camera values to lists
                daily_means.append(cam_exact_values)
                daily_mins.append(cam_exact_values)
                daily_max.append(cam_exact_values)
                daily_std.append([cam, "NA"])

                # Rename cam exact values for next step
                daily_cam_means = cam_exact_values
                daily_cam_mins = cam_exact_values
                daily_cam_max = cam_exact_values

            # Calculate overall stats for camera station
            daily_means_num = filter(lambda i:type(i) is float, daily_cam_means)
            if len(daily_means_num)==0:
                overall_mean = "NoData"
                overall_std = "NoData"
                days_over_zero = "NoData"
            else:
                overall_mean = np.mean(daily_means_num)
                overall_std = np.std(daily_means_num)
                np_cam_means = np.asarray(daily_cam_means)
                means_above_zero = np_cam_means > 0
                days_over_zero = len(np_cam_means[means_above_zero])
            if len(daily_cam_mins)==1:
                overall_min = "NoData"
            else:
                overall_min = min(daily_cam_mins)
            if len(daily_cam_max)==1:
                overall_max = "NoData"
            else:
                daily_max_num = filter(lambda i:type(i) is float, daily_cam_max)
                overall_max = max(daily_max_num)
            
            # append to overall array
            cam_stats = [cam, overall_mean, overall_min, overall_max, overall_std, days_over_zero]
            station_stats.append(cam_stats)

            # More clean up
            print "....Deleting unecessary camera files"
            try:
                arcpy.Delete_management(camera_site)
                if os.path.exists(buffer_cam):
                    arcpy.Delete_management(buffer_cam)
                shutil.rmtree(extract_folder)
            except Exception as err:
                print "........WARNING: Some files could not be deleted"

        # Messages thrown if calculations fail for a camera station
        except Exception as err:
            print "--------------------------------------------------------------------------"
            print "WARNING! There was a problem making calculations for camera station " + cam
            print "The Exception thrown was:"
            print err
            print "--------------------------------------------------------------------------"
            daily_means.append([cam, "Error collecting data"])
            daily_mins.append([cam, "Error collecting data"])
            daily_max.append([cam, "Error collecting data"])
            daily_std.append([cam, "Error collecting data"])
            station_stats.append([cam, "Error collecting data"])
            
    # Save output datasets: dailies, means, sd, min, max, # days > 0'
    print "Saving final outputs..."
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)
    if len(daily_means)>0:
        np.savetxt(os.path.join(output_folder, "DailyMeansPerStation.csv"), daily_means, fmt = '%s', delimiter=',')#header=
    if len(daily_mins)>0:
        np.savetxt(os.path.join(output_folder, "DailyMinsPerStation.csv"), daily_mins, fmt = '%s', delimiter=',')
    if len(daily_max)>0:
        np.savetxt(os.path.join(output_folder, "DailyMaxPerStation.csv"), daily_max, fmt = '%s', delimiter=',')
    if len(daily_std)>0:
        np.savetxt(os.path.join(output_folder, "DailyStdPerStation.csv"), daily_std, fmt = '%s', delimiter=',')
    if len(station_stats)>0:
        np.savetxt(os.path.join(output_folder, "OverallStatsPerStation.csv"), station_stats, fmt= '%s', delimiter=',', header="StationID, Mean, Min, Max, StDev, Days>0")


if __name__ == "__main__":
    main()
