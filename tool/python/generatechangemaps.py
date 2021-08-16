# FOR CCDC PostProcessing
import os
import numpy as np
import pandas as pd
import gdal
import click
import datetime as datetime
import multiprocessing
from fixed_thread_pool_executor import FixedThreadPoolExecutor
import scipy.io


def singlerow_execution(reccg_path, filename, dt, cols, method, mode, targeted_year, t_c, bias, i_row, out_path):
    dat_pth = os.path.join(reccg_path, filename)
    ccd_scanline = np.fromfile(dat_pth, dtype=dt)
    results_row = np.full(cols, -9999, dtype=np.int32)

    # if it is not empty
    if len(ccd_scanline) == 0:
        print('the rec_cg file {} is missing'.format(dat_pth))
        return results_row

    # sorted ccd_scanline based on pos
    # pos_list = np.unique([x['pos'] for x in ccd_scanline])
    # sorted ccd_scanline based on pos
    ccd_scanline.sort(order='pos')
    last_pos = ccd_scanline[0]['pos']
    identical_pos_curve = []
    for count, curve in enumerate(ccd_scanline):
        if curve['pos'] != last_pos:
            if len(identical_pos_curve) > 0:
                identical_pos_curve_np = np.array(identical_pos_curve)

                # identical_pos_curve_np.sort(order='t_start')
                # identical_pos_curve_np = np.array([x for x in ccd_scanline if x['pos'] == pos])
                # identical_pos_curve_np.sort(order='t_start')
                if mode == 'r': # recent disturbance year, need sort in descending
                    identical_pos_curve_np = identical_pos_curve_np[::-1]
                i_col = int((identical_pos_curve_np[0]["pos"] - 1) % 5000)
                if i_col < 0:
                    print('Processing {} failed: i_row={}; i_col={}; pos = {}'.format(filename,
                                                                                      i_row, i_col,
                                                                                      np.ceil(identical_pos_curve_np[0]["pos"])))
                    return
                for n, element in enumerate(identical_pos_curve_np):
                    prob = element['change_prob']
                    if method is 'COLD' or method is 'obcold':
                        if prob < 100: # last segment
                            continue
                        green_direction = False
                        afforestation = False
                        if (element['magnitude'][3] > t_c and element['magnitude'][2] < -t_c and
                                element['magnitude'][4] < -t_c):
                            green_direction = True
                            if mode == 'r':
                                if n - 1 >= 0:
                                    if (identical_pos_curve_np[n - 1]['coefs'][3][1] > np.abs(
                                            element['coefs'][3][1])) and \
                                            (identical_pos_curve_np[n - 1]['coefs'][2][1] < -np.abs(
                                                element['coefs'][2][1])) and \
                                            (identical_pos_curve_np[n - 1]['coefs'][4][1] < -np.abs(
                                                element['coefs'][4][1])):
                                        afforestation = True
                                    else:
                                        afforestation = False
                            elif mode == 'f':
                                if n + 1 <= len(identical_pos_curve_np) - 1:
                                    if (identical_pos_curve_np[n+1]['coefs'][3][1] > np.abs(element['coefs'][3][1])) and \
                                        (identical_pos_curve_np[n+1]['coefs'][2][1] < -np.abs(element['coefs'][2][1])) and \
                                            (identical_pos_curve_np[n+1]['coefs'][4][1] < -np.abs(element['coefs'][4][1])):
                                        afforestation = True
                                    else:
                                        afforestation = False

                        if green_direction is not True or afforestation is True:
                            # print(element['t_break'])
                            try:
                                if prob == 100:
                                    year = pd.Timestamp.fromordinal(element['t_break'] - bias).year
                                    if targeted_year == 0:  #
                                        # array[i_row,i_col] = month
                                        results_row[i_col] = year
                                        break
                                    elif year == targeted_year:
                                        results_row[i_col] = element['t_break'] - \
                                                              (pd.Timestamp.toordinal(datetime.date(year,
                                                                                                    1, 1))
                                                               + bias) + 1  # convert doy
                                        break
                            except:
                                print("Fail at row {} and col {}; t_break = {}".format(i_row + 1, i_col + 1,
                                                                                       element['t_break']))
                    elif method is 'SCCD':
                        if element['category'] % 10 == 1 and element['category'] / 10 is not 2 and \
                                element['t_break'] > 0:
                            try:
                                if prob == 100:
                                    year = pd.Timestamp.fromordinal(element['t_break'] - bias).year
                                    if targeted_year == 0:  #
                                        results_row[i_col] = year
                                        break
                                    elif year == targeted_year:
                                        # array[i_row, i_col] = element['t_break']
                                        results_row[i_col] = element['t_break'] - \
                                                              (pd.Timestamp.toordinal(datetime.date(year,
                                                                                                    1, 1))
                                                               + bias) + 1 # convert doy
                                        break
                            except:
                                print("Fail at row {} and col {}; t_break = {}".format(i_row + 1, i_col + 1,
                                                                                       element['t_break']))
                                # clear list
            identical_pos_curve = []
            identical_pos_curve.append(curve)
            # go to the next pos
            last_pos = curve['pos']
        else:
            identical_pos_curve.append(curve)

        if count == (len(ccd_scanline) - 1):  # that means all curve have been collected for last_pos
            if len(identical_pos_curve) > 0:
                identical_pos_curve_np = np.array(identical_pos_curve)
                # identical_pos_curve_np.sort(order='t_start')
                # identical_pos_curve_np = np.array([x for x in ccd_scanline if x['pos'] == pos])
                # identical_pos_curve_np.sort(order='t_start')
                if mode == 'r': # recent disturbance year, need sort in descending
                    identical_pos_curve_np = identical_pos_curve_np[::-1]
                i_col = int((identical_pos_curve_np[0]["pos"] - 1) % 5000)
                if i_col < 0:
                    print('Processing {} failed: i_row={}; i_col={}; pos = {}'.format(filename,
                                                                                      i_row, i_col,
                                                                                      np.ceil(identical_pos_curve_np[0]["pos"])))
                    return
                for n, element in enumerate(identical_pos_curve_np):
                    if method is 'COLD' or method is 'obcold':
                        if prob < 100: # last segment
                            continue
                        green_direction = False
                        afforestation = False
                        if (element['magnitude'][3] > t_c and element['magnitude'][2] < -t_c and
                                element['magnitude'][4] < -t_c):
                            green_direction = True
                            if mode == 'r':
                                if n - 1 >= 0:
                                    if (identical_pos_curve_np[n - 1]['coefs'][3][1] > np.abs(
                                            element['coefs'][3][1])) and \
                                            (identical_pos_curve_np[n - 1]['coefs'][2][1] < -np.abs(
                                                element['coefs'][2][1])) and \
                                            (identical_pos_curve_np[n - 1]['coefs'][4][1] < -np.abs(
                                                element['coefs'][4][1])):
                                        afforestation = True
                                    else:
                                        afforestation = False
                            elif mode == 'f':
                                if n + 1 <= len(identical_pos_curve_np) - 1:
                                    if (identical_pos_curve_np[n+1]['coefs'][3][1] > np.abs(element['coefs'][3][1])) and \
                                        (identical_pos_curve_np[n+1]['coefs'][2][1] < -np.abs(element['coefs'][2][1])) and \
                                            (identical_pos_curve_np[n+1]['coefs'][4][1] < -np.abs(element['coefs'][4][1])):
                                        afforestation = True
                                    else:
                                        afforestation = False

                        if green_direction is not True or afforestation is True:
                            try:
                                if prob == 100:
                                    year = pd.Timestamp.fromordinal(element['t_break'] - bias).year
                                    if targeted_year == 0:  #
                                        # array[i_row,i_col] = month
                                        results_row[i_col] = year
                                        break
                                    elif year == targeted_year:
                                        results_row[i_col] = element['t_break'] - \
                                                              (pd.Timestamp.toordinal(datetime.date(year,
                                                                                                    1, 1))
                                                               + bias) + 1  # convert doy
                                        break
                            except:
                                print("Fail at row {} and col {}; t_break = {}".format(i_row + 1, i_col + 1,
                                                                                       element['t_break']))
                    elif method is 'SCCD':
                        if element['category'] % 10 == 1 and element['category'] / 10 is not 2 and \
                                element['t_break'] > 0:
                            try:
                                if prob == 100:
                                    year = pd.Timestamp.fromordinal(element['t_break'] - bias).year
                                    if targeted_year == 0:  #
                                        results_row[i_col] = year
                                        break
                                    elif year == targeted_year:
                                        # array[i_row, i_col] = element['t_break']
                                        results_row[i_col] = element['t_break'] - \
                                                              (pd.Timestamp.toordinal(datetime.date(year,
                                                                                                    1, 1))
                                                               + bias) + 1 # convert doy
                                        break
                            except:
                                print("Fail at row {} and col {}; t_break = {}".format(i_row + 1, i_col + 1,
                                                                                       element['t_break']))
                                # clear list
            identical_pos_curve = []

    outfile = os.path.join(out_path, 'tmp_map_row{}.npy'.format(str(i_row+1).zfill(4)))
    np.save(outfile, results_row)


@click.command()
@click.option('--reccg_path', type=str, help='rec_cg folder')
@click.option('--reference_path', type=str, help='image path used to provide georeference for output images')
@click.option('--out_path', type=str, help='output folder for saving image')
@click.option('--method', type=str, default='COLD', help='COLD or SCCD')
@click.option('--targeted_year', type=int, default=0, help='outputting a binary disturbance map for a specific year')
@click.option('--mode', type=str, default='f', help='recent disturbance year (r) or first disturbance year (f)')
@click.option('--threads', type=int, default=36, help='the number for parallelization')
def main(reccg_path, reference_path, out_path, method, targeted_year, mode, threads):
    # reference_path = '/Users/coloury/Dropbox/UCONN/spatial/test_results/h016v010/recentdist_map_COLD.tif'
    # method = 'COLD'
    # reccg_path ='/Users/coloury/Dropbox/UCONN/spatial/test_results/h030v006/july_version'
    # mode = 'r'
    # out_path = '/Users/coloury/tmp'
    # targeted_year = 0
    # results_row = np.full(cols, -9999, dtype=np.int32)
    # threads = 1
    # is_mat = False
    # singlerow_execution(reccg_path, filename, dt, cols, method, mode, targeted_year, t_c, bias, i_row, results_row)


    bias = 366
    if not os.path.exists(out_path):
        os.makedirs(out_path)

    if 'obcold' in reccg_path:
        method = 'obcold'

    # ref_file = [name for name in os.listdir(reference_path) if name.startswith('LC')][0]
    # ref_image = gdal.Open(os.path.join(reference_path, '{}/{}_MTLstack'.format(ref_file, ref_file)),
                          # gdal.GA_ReadOnly)
    # ref_file = [name for name in os.listdir(reference_path) if name.startswith('CM')][0]
    ref_image = gdal.Open(reference_path, gdal.GA_ReadOnly)
    trans = ref_image.GetGeoTransform()
    proj = ref_image.GetProjection()
    cols = ref_image.RasterXSize
    rows = ref_image.RasterYSize

    if method is 'COLD' or method is 'obcold':
        # outname'obcold':
        # outname = 'breakyear_ccd_h11v9_{}_{}_{}'.format(lower_year, upper_year, method)
        dt = np.dtype([('t_start', np.int32),
                       ('t_end', np.int32),
                       ('t_break', np.int32),
                       ('pos', np.int32),
                       ('num_obs', np.int32),
                       ('category', np.short),
                       ('change_prob', np.short),
                       ('coefs', np.double, (7, 8)),
                       ('rmse', np.double, 7),
                       ('magnitude', np.double, 7)])
    elif method is 'SCCD':
        # outname = 'breakyear_sccd_h11v9_{}_{}_{}'.format(lower_year, upper_year, method)
        dt = np.dtype([('t_start', np.int32),
                       ('t_end', np.int32),
                       ('t_break', np.int32),
                       ('pos', np.int32),
                       ('num_obs', np.int32),
                       ('category', np.int16),
                       ('land_type', np.int16),
                       ('t_confirmed', np.int32),
                       ('change_prob', np.int32),
                       ('coefs', np.double, (7, 6)),
                       ('rmse', np.double, (7, 1)),
                       ('magnitude', np.double, (7, 1))])


    t_c = -200

    # reccg_path = '/Users/coloury/Dropbox/UCONN/spatial/test_results/h015v008/singlepath'
    # cols = 5000
    # rows = 5000
    # threads = 1
    # mode = 'f'
    # targeted_year = 0
    # reccg_path = '/scratch/suy20004/h016v010_results'
    filelist = [f for f in os.listdir(reccg_path) if f.endswith(".dat") and f.startswith("record_change")]

    mapgeneration_executor = FixedThreadPoolExecutor(size=threads)

    # create an empty array
    # results = np.full((cols, rows), -9999, dtype=np.int32)
    for filename in filelist:
        if method == 'obcold':
            i_row = int(filename[17:filename.find('.dat')-7]) - 1
        else:
            i_row = int(filename[17:filename.find('.dat')-5]) - 1
        mapgeneration_executor.submit(singlerow_execution, reccg_path, filename, dt, cols, method, mode,
                                      targeted_year, t_c, bias, i_row, out_path)
        # results[i_row] = results_row

    # await all tile finished
    mapgeneration_executor.drain()

    # await thread pool to stop
    mapgeneration_executor.close()

    # assemble
    tmp_map_rows = [np.load(os.path.join(out_path, 'tmp_map_row{}.npy'.format(str(x+1).zfill(4))))
                    for x in range(rows)]
    results = np.vstack(tmp_map_rows)

    for x in range(rows):
        os.remove(os.path.join(out_path, 'tmp_map_row{}.npy'.format(str(x+1).zfill(4))))

    if targeted_year > 0:
        mode_string = str(targeted_year) + 'dist_map'
    elif mode == 'r':
        mode_string = 'recentdist_map'
    elif mode == 'f':
        mode_string = 'firstdist_map'

    outname = '{}_{}.tif'.format(mode_string, method)
    outfile = os.path.join(out_path, outname)
    outdriver1 = gdal.GetDriverByName("GTiff")
    outdata = outdriver1.Create(outfile, rows, cols, 1, gdal.GDT_Int16)
    outdata.GetRasterBand(1).WriteArray(results)
    outdata.FlushCache()
    outdata.SetGeoTransform(trans)
    outdata.FlushCache()
    outdata.SetProjection(proj)
    outdata.FlushCache()

if __name__ == '__main__':
    main()
# filename = 'record_change_row122.dat'
# dat_pth = os.path.join(in_path, filename)
# ccd_scanline = np.fromfile(dat_pth, dtype=dt)
# ccd_scanline
# for file in filelist:
#     os.remove(os.path.join(reccg_path, file))
# python generatechangemaps.py --reccg_path=/scratch/suy20004/h011v009_coldspatial2/obcold/ --reference_path=/scratch/suy20004/h011v009_stack/ --out_path=/scratch/suy20004/h011v009_coldspatial2/obcold_map
