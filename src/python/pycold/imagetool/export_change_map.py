# Author: Su Ye
# generating yearly, recent and first-disturbance maps from change records
import os
import numpy as np
import pandas as pd
from osgeo import gdal
import click
import datetime as datetime
from mpi4py import MPI
from osgeo import gdal_array
import pickle
import yaml
from collections import namedtuple

PACK_ITEM = 6
SccdOutput = namedtuple("SccdOutput", "position rec_cg min_rmse nrt_mode nrt_model nrt_queue")
output_sccd = np.dtype([('t_start', np.int32), ('t_break', np.int32), ('num_obs', np.int32),
                        ('coefs', np.float32, (6, 6)), ('rmse', np.float32, 6), ('magnitude', np.float32, 6)], align=True)

output_nrtqueue = np.dtype([('clry', np.short, 6), ('clrx_since1982', np.short)], align=True)
output_nrtmodel = np.dtype([('t_start_since1982', np.short), ('num_obs', np.short), ('obs', np.short, (6, 5)),
                            ('obs_date_since1982', np.short, 5), ('covariance', np.float32, (6, 36)),
                            ('nrt_coefs', np.float32, (6, 6)), ('H', np.float32, 6), ('rmse_sum', np.uint32, 6),
                            ('cm_outputs', np.short), ('cm_outputs_date', np.short)], align=True)


# copy from /pycold/src/python/pycold/pyclassifier.py because MPI has conflicts with the pycold package in UCONN HPC.
# Dirty approach!
def extract_features(cold_plot, band, ordinal_day_list, nan_val, feature_outputs=['a0', 'a1', 'b1']):
    """
    generate features for classification based on a plot-based rec_cg and a list of days to be predicted
    Parameters
    ----------
    cold_plot: nested array
        plot-based rec_cg
    band: integer
        the predicted band number range from 0 to 6
    ordinal_day_list: list
        a list of days that this function will predict every days as a list as output
    nan_val: integer
        NA value assigned to the output
    feature_outputs: a list of outputted feature name
        it must be within [a0, c1, a1, b1,a2, b2, a3, b3, rmse]
    Returns
    -------
        feature: a list (length = n_feature) of 1-array [len(ordinal_day_list)]
    """
    features = [np.full(len(ordinal_day_list), nan_val, dtype=np.double) for x in range(len(feature_outputs))]
    for index, ordinal_day in enumerate(ordinal_day_list):
        # print(index)
        for idx, cold_curve in enumerate(cold_plot):
            if idx == len(cold_plot) - 1:
                max_days = cold_plot[idx]['t_end']
            else:
                max_days = cold_plot[idx + 1]['t_start']

            if cold_curve['t_start'] <= ordinal_day < max_days:
                for n, feature in enumerate(feature_outputs):
                    if feature == 'a0':
                        features[n][index] = cold_curve['coefs'][band][0] + cold_curve['coefs'][band][1] * \
                                             ordinal_day / 100
                        if np.isnan(features[n][index]):
                            features[n][index] = 0
                    elif feature == 'c1':
                        features[n][index] = cold_curve['coefs'][band][1] / 100
                        if np.isnan(features[n][index]):
                            features[n][index] = 0
                    elif feature == 'a1':
                        features[n][index] = cold_curve['coefs'][band][2]
                        if np.isnan(features[n][index]):
                            features[n][index] = 0
                    elif feature == 'b1':
                        features[n][index] = cold_curve['coefs'][band][3]
                        if np.isnan(features[n][index]):
                            features[n][index] = 0
                    elif feature == 'a2':
                        features[n][index] = cold_curve['coefs'][band][4]
                        if np.isnan(features[n][index]):
                            features[n][index] = 0
                    elif feature == 'b2':
                        features[n][index] = cold_curve['coefs'][band][5]
                        if np.isnan(features[n][index]):
                            features[n][index] = 0
                    elif feature == 'a3':
                        features[n][index] = cold_curve['coefs'][band][6]
                        if np.isnan(features[n][index]):
                            features[n][index] = 0
                    elif feature == 'b3':
                        features[n][index] = cold_curve['coefs'][band][7]
                        if np.isnan(features[n][index]):
                            features[n][index] = 0
                    elif feature == 'rmse':
                        features[n][index] = cold_curve['rmse'][band]
                        if np.isnan(features[n][index]):
                            features[n][index] = 0
                    else:
                        raise Exception('the outputted feature must be in [a0, c1, a1, b1,a2, b2, a3, b3, rmse]')
                break
    return features


def index_sccdpack(sccd_pack_single):
    """
    convert list of sccdpack to namedtuple to facilitate parse,
    :param sccd_pack_single: a nested list
    :return: a namedtuple SccdOutput
    """
    if len(sccd_pack_single) != PACK_ITEM:
        raise Exception("the element number of sccd_pack_single must be {}".format(PACK_ITEM))

    # convert to named tuple
    sccd_pack_single = SccdOutput(*sccd_pack_single)

    # replace the element to structured array
    if len(sccd_pack_single.rec_cg) == 0:
        sccd_pack_single = sccd_pack_single._replace(rec_cg=np.asarray(sccd_pack_single.rec_cg,
                                                                       dtype=np.float64))
    else:
        sccd_pack_single = sccd_pack_single._replace(rec_cg=np.asarray(sccd_pack_single.rec_cg,
                                                                       dtype=output_sccd))
    if len(sccd_pack_single.nrt_model) > 0:
        sccd_pack_single = sccd_pack_single._replace(nrt_model=np.asarray(sccd_pack_single.nrt_model,
                                                                          dtype=output_nrtmodel))
    if len(sccd_pack_single.nrt_queue) > 0:
        sccd_pack_single = sccd_pack_single._replace(nrt_queue=np.asarray(sccd_pack_single.nrt_queue,
                                                                          dtype=output_nrtqueue))
    return sccd_pack_single


def getcategory_cold(cold_plot, i_curve):
    t_c = -200
    if cold_plot[i_curve]['magnitude'][3] > t_c and cold_plot[i_curve]['magnitude'][2] < -t_c and \
            cold_plot[i_curve]['magnitude'][4] < -t_c:
        if cold_plot[i_curve + 1]['coefs'][3, 1] > np.abs(cold_plot[i_curve]['coefs'][3, 1]) and \
                cold_plot[i_curve + 1]['coefs'][2, 1] < -np.abs(cold_plot[i_curve]['coefs'][2, 1]) and \
                cold_plot[i_curve + 1]['coefs'][4, 1] < -np.abs(cold_plot[i_curve]['coefs'][4, 1]):
            return 3  # aforestation
        else:
            return 2  # regrowth
    else:
        return 1  # land disturbance


def getcategory_obcold(cold_plot, i_curve, last_dist_type):
    t_c = -250
    if cold_plot[i_curve]['magnitude'][3] > t_c and cold_plot[i_curve]['magnitude'][2] < -t_c and cold_plot[i_curve]['magnitude'][4] < -t_c:
        if cold_plot[i_curve + 1]['coefs'][3, 1] > np.abs(cold_plot[i_curve]['coefs'][3, 1]) and \
                cold_plot[i_curve + 1]['coefs'][2, 1] < -np.abs(cold_plot[i_curve]['coefs'][2, 1]) and \
                cold_plot[i_curve + 1]['coefs'][4, 1] < -np.abs(cold_plot[i_curve]['coefs'][4, 1]):
            return 3  # aforestation
        else:
            return 2  # regrowth
    else:
        if i_curve > 0:
            if (cold_plot[i_curve]['t_break'] - cold_plot[i_curve - 1]['t_break'] > 365.25 * 5) or (last_dist_type != 1):
                return 1
            flip_count = 0
            for b in range(5):
                if cold_plot[i_curve]['magnitude'][b + 1] * cold_plot[i_curve - 1]['magnitude'][b + 1] < 0:
                    flip_count = flip_count + 1
            if flip_count >= 4:
                return 4
            else:
                return 1
        else:
            return 1  # land disturbance


# for sccd we won't consider afforestation
def getcategory_sccd(cold_plot, i_curve):
    t_c = -200
    if cold_plot[i_curve]['magnitude'][3] > t_c and cold_plot[i_curve]['magnitude'][2] < -t_c and \
            cold_plot[i_curve]['magnitude'][4] < -t_c:
        return 2  # regrowth
    else:
        return 1  # land disturbance


@click.command()
@click.option('--reccg_path', type=str, help='rec_cg folder')
@click.option('--reference_path', type=str, help='image path used to provide georeference for output images')
@click.option('--out_path', type=str, help='output folder for saving image')
@click.option('--method', type=click.Choice(['COLD', 'OBCOLD', 'SCCDOFFLINE']), default='COLD', help='the algorithm used for processing')
@click.option('--yaml_path', type=str, help='path for yaml file')
@click.option('--year_lowbound', type=int, default=1982, help='the starting year for exporting')
@click.option('--year_uppbound', type=int, default=2020, help='the ending year for exporting')
@click.option('--coefs', default=False, help='if output coefs layers')
@click.option('--coefs_bands', default=[0, 1, 2, 3, 4, 5], help='indicate the bands for output coefs_bands,'
                                                                ' only works when coefs is True; note that the band order is'
                                                                'b,g,r,n, s1,s2, t')
def main(reccg_path, reference_path, out_path, method, year_lowbound, year_uppbound, yaml_path, coefs, coefs_bands):
    # reference_path = '/Users/coloury/Dropbox/UCONN/spatial/test_results/h016v010/recentdist_map_COLD.tif'
    # method = 'SCCDOFFLINE'
    # yaml_path = '/home/coloury/Dropbox/Documents/PyCharmProjects/HLS_NRT/config_hls.yaml'
    # reccg_path ='/home/coloury'
    # out_path = '/home/coloury'
    # year_lowbound = 1982
    # year_uppbound = 2020

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    n_process = comm.Get_size()

    if method == 'OBCOLD':
        reccg_path = os.path.join(reccg_path, 'obcold')
        out_path = os.path.join(out_path, 'obcold_maps')
    elif method == 'COLD':
        out_path = os.path.join(out_path, 'cold_maps')
    elif method == 'SCCDOFFLINE':
        out_path = os.path.join(out_path, 'sccd_maps')

    # outname'obcold':
    # outname = 'breakyear_cold_h11v9_{}_{}_{}'.format(lower_year, upper_year, method)
    if method == 'SCCDOFFLINE':
        dt = np.dtype([('t_start', np.int32),
                       ('t_break', np.int32),
                       ('num_obs', np.int32),
                       ('coefs', np.float32, (6, 6)),  # note that the slope coefficient was scaled up by 10000
                       ('rmse', np.float32, 6),
                       ('magnitude', np.float32, 6)])
    else:
        dt = np.dtype([('t_start', np.int32),
                       ('t_end', np.int32),
                       ('t_break', np.int32),
                       ('pos', np.int32),
                       ('num_obs', np.int32),
                       ('category', np.short),
                       ('change_prob', np.short),
                       ('coefs', np.float32, (7, 8)),   # note that the slope coefficient was scaled up by 10000
                       ('rmse', np.float32, 7),
                       ('magnitude', np.float32, 7)])

    if rank == 0:
        if not os.path.exists(out_path):
            os.makedirs(out_path)

        ref_image = gdal.Open(reference_path, gdal.GA_ReadOnly)
        trans = ref_image.GetGeoTransform()
        proj = ref_image.GetProjection()
        cols = ref_image.RasterXSize
        rows = ref_image.RasterYSize

        with open(yaml_path, 'r') as yaml_obj:
            config = yaml.safe_load(yaml_obj)

        config['block_width'] = int(config['n_cols'] / config['n_block_x'])  # width of a block
        config['block_height'] = int(config['n_rows'] / config['n_block_y'])  # height of a block
        config['n_blocks'] = config['n_block_x'] * config['n_block_y']  # total number of blocks
    else:
        trans = None
        proj = None
        cols = None
        rows = None
        config = None

    trans = comm.bcast(trans, root=0)
    proj = comm.bcast(proj, root=0)
    cols = comm.bcast(cols, root=0)
    rows = comm.bcast(rows, root=0)
    config = comm.bcast(config, root=0)

    ranks_percore = int(np.ceil(config['n_blocks'] / n_process))
    for i in range(ranks_percore):
        iblock = n_process * i + rank
        if iblock >= config['n_blocks']:
            break
        current_block_y = int(np.floor(iblock / config['n_block_x'])) + 1
        current_block_x = iblock % config['n_block_x'] + 1
        if method == 'OBCOLD':
            filename = 'record_change_x{}_y{}_obcold.npy'.format(current_block_x, current_block_y)
        elif method == 'COLD':
            filename = 'record_change_x{}_y{}_cold.npy'.format(current_block_x, current_block_y)
        elif method == 'SCCDOFFLINE':
            filename = 'record_change_x{}_y{}_sccd.npy'.format(current_block_x, current_block_y)

        results_block = [np.full((config['block_height'], config['block_width']), -9999, dtype=np.int32)
                         for t in range(year_uppbound - year_lowbound + 1)]

        print('processing the rec_cg file {}'.format(os.path.join(reccg_path, filename)))
        if not os.path.exists(os.path.join(reccg_path, filename)):
            # reccg_path = '/home/coloury'
            print('the rec_cg file {} is missing'.format(os.path.join(reccg_path, filename)))
            for year in range(year_lowbound, year_uppbound + 1):
                outfile = os.path.join(out_path, 'tmp_map_block{}_{}.npy'.format(iblock + 1, year))
                np.save(outfile, results_block[year - year_lowbound])
            continue
        if method == 'SCCDOFFLINE':
            file = open(os.path.join(reccg_path, filename), 'rb')
            cold_block = []
            while True:
                try:
                    cold_block.append(index_sccdpack(pickle.load(file)))
                except EOFError:
                    break
            file.close()
        else:
            cold_block = np.array(np.load(os.path.join(reccg_path, filename)), dtype=dt)
        # cold_block = [np.array(element, dtype=dt) for element in cold_block]
        # if len(cold_block) == 0:
        #     print('the rec_cg file {} is missing'.format(dat_pth))
        #     for year in range(year_lowbound, year_uppbound+1):
        #         outfile = os.path.join(out_path, 'tmp_map_block{}_{}.npy'.format(iblock + 1, year))
        #         np.save(outfile, results_block[year - year_lowbound])
        #     continue
        if method == 'SCCDOFFLINE':
            for count, plot in enumerate(cold_block):
                for i_count, curve in enumerate(plot.rec_cg):
                    if curve['t_break'] == 0 or count == (len(cold_block) - 1):  # last segment
                        continue

                    i_col = int((plot.position  - 1) % config['n_cols']) - \
                            (current_block_x - 1) * config['block_width']
                    i_row = int((plot.position  - 1) / config['n_cols']) - \
                            (current_block_y - 1) * config['block_height']
                    if i_col < 0:
                        print('Processing {} failed: i_row={}; i_col={} for {}'.format(filename,
                                                                                       i_row, i_col, filename))
                        # return
                    current_dist_type = getcategory_sccd(plot.rec_cg, i_count)
                    break_year = pd.Timestamp.fromordinal(curve['t_break']).year
                    if break_year < year_lowbound or break_year > year_uppbound:
                        continue
                    results_block[break_year -
                                  year_lowbound][i_row][i_col] = current_dist_type * 1000 + curve['t_break'] - \
                                                                 (pd.Timestamp.toordinal(datetime.date(break_year, 1, 1))) + 1
        else:
            cold_block.sort(order='pos')
            current_processing_pos = cold_block[0]['pos']
            current_dist_type = 0
            for count, curve in enumerate(cold_block):
                if curve['pos'] != current_processing_pos:
                    current_processing_pos = curve['pos']
                    current_dist_type = 0

                if curve['change_prob'] < 100 or curve['t_break'] == 0 or count == (len(cold_block) - 1):  # last segment
                    continue

                i_col = int((curve["pos"] - 1) % config['n_cols']) - \
                        (current_block_x - 1) * config['block_width']
                i_row = int((curve["pos"] - 1) / config['n_cols']) - \
                        (current_block_y - 1) * config['block_height']
                if i_col < 0:
                    dat_pth = '?'
                    print('Processing {} failed: i_row={}; i_col={} for {}'.format(filename, i_row, i_col, dat_pth))
                    return

                if method == 'OBCOLD':
                    current_dist_type = getcategory_obcold(cold_block, count, current_dist_type)
                else:
                    current_dist_type = getcategory_cold(cold_block, count)
                break_year = pd.Timestamp.fromordinal(curve['t_break']).year
                if break_year < year_lowbound or break_year > year_uppbound:
                    continue
                results_block[break_year - year_lowbound][i_row][i_col] = current_dist_type * 1000 + curve['t_break'] - \
                    (pd.Timestamp.toordinal(datetime.date(break_year, 1, 1))) + 1
            # e.g., 1315 means that disturbance happens at doy of 315

        # save the temp dataset out
        for year in range(year_lowbound, year_uppbound + 1):
            outfile = os.path.join(out_path, 'tmp_map_block{}_{}.npy'.format(iblock + 1, year))
            np.save(outfile, results_block[year - year_lowbound])

    # wait for all processes
    comm.Barrier()

    if rank == 0:
        # assemble
        for year in range(year_lowbound, year_uppbound + 1):
            tmp_map_blocks = [np.load(os.path.join(out_path, 'tmp_map_block{}_{}.npy'.format(x + 1, year)))
                              for x in range(config['n_blocks'])]

            results = np.hstack(tmp_map_blocks)
            results = np.vstack(np.hsplit(results, config['n_block_x']))

            for x in range(config['n_blocks']):
                os.remove(os.path.join(out_path, 'tmp_map_block{}_{}.npy'.format(x + 1, year)))
            mode_string = str(year) + '_break_map'
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

        # output recent disturbance year
        recent_dist = np.full((config['n_rows'], config['n_cols']), 0, dtype=np.int16)
        for year in range(year_lowbound, year_uppbound + 1):
            mode_string = str(year) + '_break_map'
            outname = '{}_{}.tif'.format(mode_string, method)
            breakmap = gdal_array.LoadFile(os.path.join(out_path, outname))
            recent_dist[(breakmap / 1000).astype(np.byte) == 1] = year
        mode_string = 'recent_disturbance_map'
        outname = '{}_{}.tif'.format(mode_string, method)
        outfile = os.path.join(out_path, outname)
        outdriver1 = gdal.GetDriverByName("GTiff")
        outdata = outdriver1.Create(outfile, rows, cols, 1, gdal.GDT_Int16)
        outdata.GetRasterBand(1).WriteArray(recent_dist)
        outdata.FlushCache()
        outdata.SetGeoTransform(trans)
        outdata.FlushCache()
        outdata.SetProjection(proj)
        outdata.FlushCache()

        first_dist = np.full((config['n_rows'], config['n_cols']), 0, dtype=np.int16)
        for year in range(year_uppbound, year_lowbound - 1, -1):
            mode_string = str(year) + '_break_map'
            outname = '{}_{}.tif'.format(mode_string, method)
            breakmap = gdal_array.LoadFile(os.path.join(out_path, outname))
            first_dist[(breakmap / 1000).astype(np.byte) == 1] = year
        mode_string = 'first_disturbance_map'
        outname = '{}_{}.tif'.format(mode_string, method)
        outfile = os.path.join(out_path, outname)
        outdriver1 = gdal.GetDriverByName("GTiff")
        outdata = outdriver1.Create(outfile, rows, cols, 1, gdal.GDT_Int16)
        outdata.GetRasterBand(1).WriteArray(first_dist)
        outdata.FlushCache()
        outdata.SetGeoTransform(trans)
        outdata.FlushCache()
        outdata.SetProjection(proj)
        outdata.FlushCache()


if __name__ == '__main__':
    main()
