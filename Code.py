print(''' 

Welcome to the land filtering module
this module provides filtering automation processes in order to filter the parcels of Stafford, Virginia.
The filtering processing will be:
1. Boundary Filter: Include the parcels that locate inside the development area inside Stafford.
2. Landuse Filter: Include just the parcels of the results of boundary filter that intersects with the landuses defined by the user.
3. Single Parcel Filter: Include the parcels of the results of landuse filter that meet the provided entries by the user.
4. Multi Parcel Assemblage Filter: Include the parcels of the results of landuse filter to include assembles of parcels that meets provided entries by the user.


''')

import geopandas as gpd
import pandas as pd
import os
import sys
from shapely.geometry import LineString
import uuid


def main():
    try:
        if len(sys.argv) != 1:
            raise ValueError('Usage: python script_name.py')
        # Input file paths
        input_parcels_path = r'inputs\parcels\StaffordZoning.parquet'
        input_boundary_path = r'inputs\dev_area'

        def read_shapefile(filepath):
            try:
                gdf = gpd.read_file(filepath)
                # print(f"Successfully read {filepath}")
                return gdf
            except FileNotFoundError:
                print(f"Error: The file {filepath} was not found.")
            except gpd.errors.DriverError:
                print(f"Error: The file {filepath} is not a valid shapefile or the driver is not supported.")
            except Exception as e:
                print(f"An unexpected error occurred while reading {filepath}: {e}")
            return None

        # Read the shapefiles into separate GeoDataFrames
        # print('Welcome')
        print('Reading Datasets .. ')
        print('                                          ')
        parcels = gpd.read_parquet(input_parcels_path)

        # Ensure there are no None geometries
        parcels = parcels[parcels['geometry'].notna()]

        # parcels = parcels.explode(index_parts=True)
        # parcels = parcels.reset_index(drop=True)

        # Initialize an empty list to hold the GeoDataFrames
        gdf_list = []

        # Loop through the files in the directory
        for filename in os.listdir(input_boundary_path):
            if filename.endswith('.geojson'):
                file_path = os.path.join(input_boundary_path, filename)
                gdf = gpd.read_file(file_path)
                gdf_list.append(gdf)

        # Concatenate all the GeoDataFrames into a single GeoDataFrame
        combined_gdf = pd.concat(gdf_list, ignore_index=True)

        # Optionally, ensure the combined GeoDataFrame has a valid geometry column
        boundary = gpd.GeoDataFrame(combined_gdf, geometry='geometry')

        landuseList = ['Residential', 'Mixed Use – Commercial/ Residential', 'Downtown Stafford - Mixed Use']
        landuse = boundary[boundary['AreaType'].isin(landuseList)]

        # Define the columns that will be included for parcels
        parcels_columns = ['PRCLID', 'LRSN', 'PRCLKEY', 'geometry']
        parcels = parcels[parcels_columns]

        # Rename the remain columns to match the same names of the table in the database
        parcels = parcels.copy()
        parcels.rename(columns={'PRCLKEY': 'parcel_key', 'PRCLID':'parcel_ID'}, inplace=True)

        # Modify the data types of columns that contains numbers
        parcels[['LRSN', 'parcel_key']] = parcels[['LRSN', 'parcel_key']].astype(int)
        parcels[['LRSN', 'parcel_key']] = parcels[['LRSN', 'parcel_key']].astype(str)

        # Recalculate area columns ['Sq_feet'] and ['Land_area']
        parcels['Sq_feet'] = parcels.geometry.area
        parcels['Land_area'] = parcels['Sq_feet'] / 43560
        parcels['Land_area'] = parcels['Land_area'].round(6)
        parcels = parcels.drop(columns='Sq_feet')

        # Use uuid module to fill the ['Id'] column with IDs with this format: 7FD82786-C9D5-41B5-861C-0599E14A85E5
        parcels['Id'] = [str(uuid.uuid4()).upper() for _ in range(len(parcels))]

        # Check if all geodataframes were read successfully
        if parcels is not None:
            print("Parcels shapefile was read successfully.")
        else:
            print("Failed to read parcels shapefile.")

        if boundary is not None:
            print("Boundary shapefile was read successfully.")
        else:
            print("Failed to read boundary shapefile.")

        if landuse is not None:
            print("Landuse shapefile was read successfully.")
        else:
            print("Failed to read landuse shapefile.")

        print('                                          ')

        # check for duplications in parcels and drop them if founded
        parcels_dubs_check = parcels.duplicated(subset='geometry')
        duplicated_parcels = parcels[parcels_dubs_check]
        print(f'Total parcels: {parcels.shape[0]} parcels')
        print('                                          ')
        if duplicated_parcels.shape[0] == 0:
            print('No duplicates found')
        elif duplicated_parcels.shape[0] >= 1:
            print(f'Count of duplicates founded: {duplicated_parcels.shape[0]}')
            print('Dropping duplicates ...')
            parcels = parcels.drop_duplicates(subset='geometry', keep='first')
            print(f'Count of parcels after excluding duplicates: {parcels.shape[0]} parcels')

        boundary = boundary.to_crs('EPSG:2924')
        landuse = landuse.to_crs('EPSG:2924')
        parcels = parcels.to_crs('epsg:2924')

        # remove spaces from the parcel_ID column
        parcels['parcel_ID'] = parcels['parcel_ID'].str.replace(' ', '')

        # ======================================================================================================================

        # Initiate the first filter (Boundary Filtering)

        print('======================================================================================================================')
        print('                                          ')
        print('                                          ')
        print('Boundary filtering initiated.')
        # print('This might take several minutes.')
        # print('Please, wait ...')

        parcels['interior_point'] = parcels['geometry'].apply(lambda geom: geom.representative_point())
        po_gdf = gpd.GeoDataFrame(parcels, geometry='interior_point')
        join = gpd.sjoin(po_gdf, boundary, how='inner', predicate='intersects')
        intersecting_parcels = parcels[parcels.index.isin(join.index)]
        intersecting_parcels = intersecting_parcels.drop(columns='interior_point')
        print('                                          ')
        parcels_number = intersecting_parcels.shape[0]

        ints_prcl_exp = intersecting_parcels.to_crs('EPSG:4326')
        ints_prcl_exp = ints_prcl_exp.drop(columns='Id')
        print('Boundary filtering completed.')
        print(f'Count of parcels that intersects with the development areas = {parcels_number} Parcels')
        print('                                          ')

        while True:
            export_option = str(input('Would you like to export the GeoJSON file of the filtered parcels?  Yes/No ?'))
            if export_option in ['Yes', 'yes', 'Y', 'y']:
                print(f'Your option: Yes')
                print('Exporting...')

                def boundary_filter_export(gdf, directory, filename):
                    # Ensure the directory exists
                    os.makedirs(directory, exist_ok=True)

                    # Define the full path for the output file
                    filepath = os.path.join(directory, filename)

                    # Export the GeoDataFrame to a GeoJSON file
                    gdf.to_file(filepath)

                # Specify the directory and filename
                output_path = r'outputs\GIS\boundary_filter'
                output_filename = 'boundary_filter_parcels.geojson'

                # Call the function to export the GeoDataFrame
                boundary_filter_export(ints_prcl_exp, output_path, output_filename)
                print(r'GeoJSON file exported to: outputs\GIS\boundary_filter')
                break
            if export_option in ['No', 'no', 'N', 'n']:
                print(f'Your option: No')
                break
            else:
                print('Invalid option')
                continue
    except ValueError as ve:
        # Handle value errors (e.g., incorrect number of arguments)
        print(f"Error: {ve}")
        print("Please provide the correct number of arguments.")
    except FileNotFoundError:
        # Handle file not found errors
        print(f"Error: The file '{sys.argv[0]}' was not found.")
    except Exception as e:
        # Handle any other unexpected exceptions
        print(f"An unexpected error occurred: {e}")

    # ======================================================================================================================

    # Initiate the second filter (designated future land use)

    print('======================================================================================================================')
    print('                                          ')
    print('Land use filter initiated.')
    print('''Please, insert the land use types
    1 = Residential
    2 = Mixed Use – Commercial/ Residential
    3 = Downtown Stafford - Mixed Use
    ''')

    landuse_dict = {'1': 'Residential', '2': 'Mixed Use – Commercial/ Residential', '3': 'Downtown Stafford - Mixed Use'}
    landuse_print = {'1': 'Residential', '2': 'Mixed Use – Commercial/ Residential', '3': 'Downtown Stafford - Mixed Use'}

    while True:
        try:
            land_use_list = []
            land_use_print = []
            if len(land_use_list) == 3 and len(land_use_print) == 3:
                break
            user_defined_landuse = input('Please, input the land use(s) (separate multiple entries with commas): ')
            if len(user_defined_landuse) > 1:
                input_num = [x.strip() for x in user_defined_landuse.split(',')]
                for n in input_num:
                    if int(n) > 3 or int(n) < 1:
                        print('Invalid entry')
                        continue
                for key in input_num:
                    land_use_list.append(landuse_dict[key])
                    land_use_print.append(landuse_print[key])
                print(f'Entered land use {land_use_print}')
                print('                                          ')
                # print('Land use filtering started...')

                f_landuse = landuse[landuse['AreaType'].isin(land_use_list)]

                # parcels_future_land_use = intersecting_parcels[intersecting_parcels.within(f_landuse.unary_union)]
                intersecting_parcels['interior_point'] = intersecting_parcels['geometry'].apply(lambda geom: geom.representative_point())
                points_gdf = gpd.GeoDataFrame(intersecting_parcels, geometry='interior_point')
                joined = gpd.sjoin(points_gdf, f_landuse, how='inner', predicate='intersects')
                parcels_future_land_use = intersecting_parcels[intersecting_parcels.index.isin(joined.index)]
                parcels_future_land_use = parcels_future_land_use.drop(columns='interior_point')
                parcels_future_land_use['Land_use'] = joined['AreaType']

                print('                                          ')
                print('Land use filtering completed.')
                print('                                          ')
                print(f'Land use filter results= {parcels_future_land_use.shape[0]} parcels(s)')
            else:
                if int(user_defined_landuse) > 3 or int(user_defined_landuse) < 1:
                    print('Invalid entry')
                    continue
                else:
                    land_use_list.append(str(landuse_dict[str(user_defined_landuse)]))
                    land_use_print.append(str(landuse_print[str(user_defined_landuse)]))
                    print(f'Entered land use {land_use_print}')
                    print('                                          ')
                    print('Land use filtering started...')

                    f_landuse = landuse[landuse['AreaType'] == land_use_list[0]]

                    # parcels_future_land_use = intersecting_parcels[intersecting_parcels.intersects(f_landuse.unary_union)]
                    intersecting_parcels['interior_point'] = intersecting_parcels['geometry'].apply(
                        lambda geom: geom.representative_point())
                    points_gdf = gpd.GeoDataFrame(intersecting_parcels, geometry='interior_point')
                    joined = gpd.sjoin(points_gdf, f_landuse, how='inner', predicate='intersects')
                    parcels_future_land_use = intersecting_parcels[intersecting_parcels.index.isin(joined.index)]
                    parcels_future_land_use = parcels_future_land_use.drop(columns='interior_point')
                    parcels_future_land_use['Land_use'] = joined['AreaType']

                    print('                                          ')
                    print('Land use filtering completed.')
                    print('                                          ')
                    print(f'Land use filter results= {parcels_future_land_use.shape[0]} parcels(s)')

            while True:
                repeat_landuse_filter = input('Would you like to repeat the filter of land use? Yes/No: ')
                if repeat_landuse_filter.lower() in ['no', 'n']:
                    exit_f = True  # Set a flag to exit the outer loop
                    print('Your option: No')
                    break  # Exit the inner loop
                elif repeat_landuse_filter.lower() in ['yes', 'y']:
                    print('Your option: Yes')
                    break  # Exit the inner loop and repeat the outer loop
                else:
                    print('Invalid option. Please enter Yes/No.')

            if 'exit_f' in locals() and exit_f:
                break  # Exit the outer loop if the user chose 'No'

        except ValueError:
            print("Please provide a valid number of land use.")
            continue
        except KeyError:
            print("Please provide a valid number of land use.")
            continue
        except Exception as lu_err:
            # Handle any other unexpected exceptions
            print(f"An unexpected error occurred: {lu_err}")
            continue

    # ======================================================================================================================

    # initiate the single parcel filter

    print('                                          ')
    print(
        '======================================================================================================================')
    print('                                          ')
    while True:
        try:
            print('Single parcel filter initiated')
            print('                                          ')
            user_defined_land_area = str(
                input('Please, input the minimum and the maximum land area (separate them with commas): '))
            input_land_area = [x.strip() for x in user_defined_land_area.split(',')]
            min_land_area = int(input_land_area[0])
            max_land_area = int(input_land_area[1])
            print(f'Entered minimum land area= {min_land_area} Acres')
            print(f'Entered maximum land area= {max_land_area} Acres')
            print('                                          ')
            print('Land area filtering started...')
            land_area_filter = parcels_future_land_use[(parcels_future_land_use['Land_area'] >= min_land_area) & (
                        parcels_future_land_use['Land_area'] <= max_land_area)]
            print('                                          ')
            # user_defined_threshold = int(input('Please, input the maximum assessment value: $'))
            # print(f'Entered maximum assessment value= {user_defined_threshold} $')
            # land_area = land_area_filter['Land_area']
            # bldg_assess = land_area_filter['Building_assessment']
            # threshold_filter = land_area_filter[(bldg_assess / land_area) <= user_defined_threshold]
            # print('Land area filtering completed.')
            count_single_filter_res = land_area_filter.shape[0]
            print('                                          ')
            print(f'Single parcel filter results= {count_single_filter_res} parcel(s)')
            print('                                          ')

            while True:
                repeat_single_filter = input('Would you like to repeat the single parcel filter? Yes/No: ')
                if repeat_single_filter.lower() in ['no', 'n']:
                    exit_flg = True  # Set a flag to exit the outer loop
                    print('Your option: No')
                    break  # Exit the inner loop
                elif repeat_single_filter.lower() in ['yes', 'y']:
                    print('Your option: Yes')
                    break  # Exit the inner loop and repeat the outer loop
                else:
                    print('Invalid option. Please enter Yes/No.')

            if 'exit_flg' in locals() and exit_flg:
                break  # Exit the outer loop if the user chose 'No'

        except IndexError:
            print("Please re-enter the land use numbers with a comma to separate between min and max like (min,max)")
            continue
        except Exception as sf_err:
            # Handle any other unexpected exceptions
            print(f"An unexpected error occurred: {sf_err}")
            continue

    print('                                          ')
    single_parcel = land_area_filter.copy()
    print('Single Parcel Filtering completed.')

    # ======================================================================================================================

    # initiate the multi parcel assemblage

    # Buffer the parcels
    input_parc_assemb = parcels_future_land_use.copy()
    buffer_size = 0.001
    input_parc_assemb['geometry'] = input_parc_assemb.buffer(buffer_size)

    # Parameters for finding assemblage opportunities
    max_parcels = 5
    print('                                          ')
    print(
        '======================================================================================================================')
    print('                                          ')
    print('Multi Parcels Assemblage filter initiated.')
    print('                                          ')
    while True:
        try:
            min_area = int(input('Please, input the minimum land area in acres: '))
            max_area = int(input('Please, input the maximum land area in acres: '))
            print(f'Entered minimum land area= {min_area} Acres')
            print(f'Entered maximum land area= {max_area} Acres')
            print('                                          ')
            print(f'Assemblage filter started with {input_parc_assemb.shape[0]} Parcels')
            print('                                          ')
            print('Filtering...')

            # Initialize an empty list to store assemblages and unions
            assemblages = []
            unions = []

            # Iterate through each parcel to find potential assemblages
            for parcel in input_parc_assemb.itertuples():
                nearby = input_parc_assemb[input_parc_assemb.geometry.intersects(parcel.geometry)]

                if (len(nearby) <= max_parcels) and (len(nearby) > 0):
                    # Convert the parcel to a GeoDataFrame
                    parcel_df = gpd.GeoDataFrame([parcel._asdict()], crs=input_parc_assemb.crs)

                    # Create an assemblage by combining the nearby parcels and the current parcel
                    assemblage = gpd.GeoDataFrame(pd.concat([nearby, parcel_df], ignore_index=True), crs=input_parc_assemb.crs)
                    assemblage = assemblage.drop_duplicates(subset='Id', keep='first')

                    # merge parcels/polygons of each assembly to create an outer boudaries for each assembly
                    assemblage_union = assemblage.unary_union
                    union_series = gpd.GeoSeries(assemblage_union)
                    union_gdf = gpd.GeoDataFrame(geometry=union_series, crs=input_parc_assemb.crs)

                    # Calculate the total area of the assemblage
                    total_area = assemblage['Land_area'].sum()

                    # Check if the total area meets the minimum area requirement
                    if (total_area >= min_area) and (total_area <= max_area):
                        assemblages.append(assemblage)
                        unions.append(union_gdf)

            if assemblages:
                # Combine all valid assemblages into a single GeoDataFrame
                assemblages_df = gpd.GeoDataFrame(pd.concat(assemblages, ignore_index=True), crs=input_parc_assemb.crs)

                # Combine unions of all valid assemblages into a single GeoDataFrame
                union_assemb = gpd.GeoDataFrame(pd.concat(unions, ignore_index=True), crs=input_parc_assemb.crs)

                # Preparing for change the geometries of the parcels of the assemblages to the actual geometries before create the buffer of 0.001
                # Assuming 'id' is the unique identifier column present in both GeoDataFrames
                unique_id_column = 'Id'

                # Set the unique identifier as the index to facilitate matching
                assemblages_df = assemblages_df.set_index(unique_id_column)
                parcels_future_land_use = parcels_future_land_use.set_index(unique_id_column)

                # Find the intersection of indices that exist in both GeoDataFrames
                common_indices = assemblages_df.index.intersection(parcels_future_land_use.index)

                # Filter out unmatched identifiers from parcels_future_land_use
                parcels_matched = parcels_future_land_use.loc[common_indices]

                # Update the geometry of the matched records in assemblages_df
                assemblages_df.loc[common_indices, 'geometry'] = parcels_matched['geometry']

                # Reset the index
                assemblages_df.reset_index(inplace=True)
                parcels_future_land_use.reset_index(inplace=True)

                # check for parcels duplications in the assemblages GeoDataFrame and drop them if found
                assemblages_df = assemblages_df.drop(columns='Index')
                assemb_dubs_check = assemblages_df.duplicated(subset='Id')
                duplicated_assembs = assemblages_df[assemb_dubs_check]
                assemblages_df = assemblages_df.drop_duplicates(subset='Id', keep='first')
                assemblages_df = assemblages_df.reset_index(drop=True)

                # check for parcels duplications in the unions GeoDataFrame and drop them if found
                union_dubs_check = union_assemb.duplicated(subset='geometry')
                duplicated_uinion = union_assemb[union_dubs_check]
                union_assemb = union_assemb.drop_duplicates(subset='geometry', keep='first')
                union_assemb = union_assemb.reset_index(drop=True)

                assemb_parcel_count = assemblages_df.shape[0]
                print(f'Assemblage filter results= {assemb_parcel_count} Parcel(s)')

                assemb_count = union_assemb.shape[0]
                print(f'Count of assemblages= {assemb_count} Assemblage(s)')
            else:
                assemblages_df = gpd.GeoDataFrame(columns=input_parc_assemb.columns, crs=input_parc_assemb.crs)
                assemb_parcel_count = assemblages_df.shape[0]
                print(f'Assemblage filter results= {assemb_parcel_count} Parcel(s)')

            # # Save the assemblages to a GeoJSON file
            # output_path = r'outputs\assemblage_filter\assemblages.geojson'
            # assemblages_df = assemblages_df.to_crs('epsg:4326')
            # assemblages_df.to_file(output_path, driver='GeoJSON')
            #
            # print(r'Assemblages exported to: outputs\assemblage_filter\assemblages.geojson')
            break
        except ValueError:
            print('Please, re-enter the number of min/max in number format like >>> (10)')
            print('                                          ')
            print('                                          ')
            continue
        except Exception as assm_err:
            # Handle any other unexpected exceptions
            print(f"An unexpected error occurred: {assm_err}")
            continue

    # convert the outer boundaries of unions of assemblages from polygons to linestring
    line_strings = union_assemb.geometry.apply(lambda x: LineString(x.exterior.coords))
    line_strings_gdf = gpd.GeoDataFrame(geometry=line_strings)

    # ======================================================================================================================

    # merge the GeoDataFrames of the single parcel filter and the multi parcels assemblages
    merged_gdf = gpd.GeoDataFrame(pd.concat([single_parcel, assemblages_df], ignore_index=True))

    # Remove duplicate rows based on all columns
    merged_dubs_check = merged_gdf['Id'].duplicated()
    duplicated_merge = merged_gdf[merged_dubs_check]
    merged_gdf = merged_gdf.drop_duplicates(subset='Id', keep='first')
    merged_gdf = merged_gdf.reset_index(drop=True)

    print('                                          ')
    print('                                          ')
    print('Results of single parcels filter and multi parcel assemblages have been checked for removing duplications if found. And have been merged in one GeoJSON file.')
    print('                                          ')
    print('                                          ')
    print(f'The final result contains: {merged_gdf.shape[0]} parcels')
    print('Exporting...')
    merged_gdf = merged_gdf.drop(columns='Id')
    merged_gdf = merged_gdf.to_crs('EPSG:4326')
    merged_gdf.to_file(r'outputs\GIS\single_and_assemblages\Parcels.geojson')
    line_strings_gdf = line_strings_gdf.to_crs('EPSG:4326')
    line_strings_gdf.to_file(r'outputs\GIS\single_and_assemblages\Boundaries_assemblages.geojson')
    print(r'GeoJSON files exported to: outputs\GIS\single_and_assemblages')
    print('                                          ')
    print('                                          ')
    exc = merged_gdf.drop(columns='geometry')
    exc = pd.DataFrame(exc)
    exc.to_csv(r'outputs\Excel\Parcels.csv')
    print(r'Excel file exported to: outputs\Excel\Parcels.csv')



if __name__ == "__main__":
    main()
