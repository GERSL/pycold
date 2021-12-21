# FOR CCDC PostProcessing
import os
import numpy as np
import pandas as pd
import gdal
import click
import datetime as datetime
import multiprocessing
import scipy.io
import yaml
from mpi4py import MPI
from osgeo import gdal_array


def getcategory_cold(ccd_plot, i_curve):
    t_c = -200
    if ccd_plot[i_curve]['magnitude'][3] > t_c and ccd_plot[i_curve]['magnitude'][2] < -t_c and \
            ccd_plot[i_curve]['magnitude'][4] < -t_c:
        if ccd_plot[i_curve + 1]['coefs'][3, 1] > np.abs(ccd_plot[i_curve]['coefs'][3, 1]) and \
                ccd_plot[i_curve + 1]['coefs'][2, 1] < -np.abs(ccd_plot[i_curve]['coefs'][2, 1]) and \
                ccd_plot[i_curve + 1]['coefs'][4, 1] < -np.abs(ccd_plot[i_curve]['coefs'][4, 1]):
            return 3  # aforestation
        else:
            return 2  # regrowth
    else:
        return 1  # land disturbance


def getcategory_obcold(ccd_plot, i_curve, last_dist_type):
    t_c = -250
    if ccd_plot[i_curve]['magnitude'][3] > t_c and ccd_plot[i_curve]['magnitude'][2] < -t_c and ccd_plot[i_curve]['magnitude'][4] < -t_c:
        if ccd_plot[i_curve + 1]['coefs'][3, 1] > np.abs(ccd_plot[i_curve]['coefs'][3, 1]) and \
                ccd_plot[i_curve + 1]['coefs'][2, 1] < -np.abs(ccd_plot[i_curve]['coefs'][2, 1]) and \
                ccd_plot[i_curve + 1]['coefs'][4, 1] < -np.abs(ccd_plot[i_curve]['coefs'][4, 1]):
            return 3  # aforestation
        else:
            return 2  # regrowth
    else:
        if i_curve > 0:
            if (ccd_plot[i_curve]['t_break'] - ccd_plot[i_curve-1]['t_break'] > 365.25 * 5) or (last_dist_type != 1):
                return 1
            flip_count = 0
            for b in range(5):
                if ccd_plot[i_curve]['magnitude'][b+1] * ccd_plot[i_curve-1]['magnitude'][b+1] < 0:
                    flip_count = flip_count + 1
            if flip_count >= 4:
                return 4
            else:
                return 1
        else:
            return 1  # land disturbance


@click.command()
@click.option('--reccg_path', type=str, help='rec_cg folder')
@click.option('--reference_path', type=str, help='image path used to provide georeference for output images')
@click.option('--out_path', type=str, help='output folder for saving image')
@click.option('--method', type=str, default='COLD', help='COLD, OBCOLD or SCCD')
@click.option('--year_lowbound', type=int, default=1982, help='the starting year for exporting')
@click.option('--year_uppbound', type=int, default=2020, help='the ending year for exporting')
def main(reccg_path, reference_path, out_path, method, year_lowbound, year_uppbound):
    # reference_path = '/Users/coloury/Dropbox/UCONN/spatial/test_results/h016v010/recentdist_map_COLD.tif'
    # method = 'COLD'
    # reccg_path ='/Users/coloury/Dropbox/UCONN/spatial/test_results/h030v006/july_version'
    # mode = 'r'
    # out_path = '/Users/coloury/tmp'
    # targeted_year = 0
    # results_block = np.full(cols, -9999, dtype=np.int32)
    # threads = 1
    # is_mat = False
    # singlerow_execution(reccg_path, filename, dt, cols, method, mode, targeted_year, t_c, bias, i_row, results_block)

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    n_process = comm.Get_size()
    bias = 366
    if method == 'OBCOLD':
        reccg_path = os.path.join(reccg_path, 'obcold')
        out_path = os.path.join(out_path, 'obcold_maps')
    elif method == 'COLD':
        out_path = os.path.join(out_path, 'cold_maps')


    dt = None
    if method == 'COLD' or method == 'OBCOLD':
        # outname'obcold':
        # outname = 'breakyear_ccd_h11v9_{}_{}_{}'.format(lower_year, upper_year, method)
        dt = np.dtype([('t_start', np.int32),
                       ('t_end', np.int32),
                       ('t_break', np.int32),
                       ('pos', np.int32),
                       ('num_obs', np.int32),
                       ('category', np.short),
                       ('change_prob', np.short),
                       ('coefs', np.float32, (7, 8)),
                       ('rmse', np.float32, 7),
                       ('magnitude', np.float32, 7)])
    elif method == 'SCCD':
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

    if rank == 0:
        if not os.path.exists(out_path):
            os.makedirs(out_path)

        ref_image = gdal.Open(reference_path, gdal.GA_ReadOnly)
        trans = ref_image.GetGeoTransform()
        proj = ref_image.GetProjection()
        cols = ref_image.RasterXSize
        rows = ref_image.RasterYSize

        with open('{}/spatial/parameters.yaml'.format(os.path.dirname(os.path.dirname(os.getcwd()))), 'r') as yaml_obj:
            parameters = yaml.safe_load(yaml_obj)

        obcold_params = parameters['obcold']
        obcold_params.update(parameters['common'])
        obcold_params['block_width'] = int(obcold_params['n_cols'] / obcold_params['n_block_h'])  # width of a block
        obcold_params['block_height'] = int(obcold_params['n_rows'] / obcold_params['n_block_v'])  # height of a block
        obcold_params['n_blocks'] = obcold_params['n_block_h'] * obcold_params['n_block_v']  # total number of blocks
    else:
        trans = None
        proj = None
        cols = None
        rows = None
        obcold_params = None

    trans = comm.bcast(trans, root=0)
    proj = comm.bcast(proj, root=0)
    cols = comm.bcast(cols, root=0)
    rows = comm.bcast(rows, root=0)
    obcold_params = comm.bcast(obcold_params, root=0)

    ranks_percore = int(np.ceil(obcold_params['n_blocks'] / n_process))
    print(rank, n_process)
    for i in range(ranks_percore):
        iblock = n_process * i + rank
        if iblock > obcold_params['n_blocks']:
            break
        current_block_y = int(np.floor(iblock / obcold_params['n_block_h'])) + 1
        current_block_x = iblock % obcold_params['n_block_h'] + 1
        if method == 'OBCOLD':
            filename = 'record_change_h{}_v{}_obcold.dat'.format(current_block_x, current_block_y)
        else:
            filename = 'record_change_h{}_v{}_cold.dat'.format(current_block_x, current_block_y)
        dat_pth = os.path.join(reccg_path, filename)
        ccd_block = np.fromfile(dat_pth, dtype=dt)
        results_block = [np.full((obcold_params['block_height'], obcold_params['block_width']), -9999, dtype=np.int32)
                         for t in range(year_uppbound - year_lowbound + 1)]
        if len(ccd_block) == 0:
            print('the rec_cg file {} is missing'.format(dat_pth))
            for year in range(year_lowbound, year_uppbound+1):
                outfile = os.path.join(out_path, 'tmp_map_block{}_{}.npy'.format(iblock + 1, year))
                np.save(outfile, results_block[year - year_lowbound])
            continue

        ccd_block.sort(order='pos')

        current_processing_pos = ccd_block[0]['pos']
        current_dist_type = 0
        for count, curve in enumerate(ccd_block):
            if curve['pos'] != current_processing_pos:
                current_processing_pos = curve['pos']
                current_dist_type = 0

            if curve['change_prob'] < 100 or curve['t_break'] == 0:  # last segment
                continue

            i_col = int((curve["pos"] - 1) % obcold_params['n_cols']) - \
                    (current_block_x - 1) * obcold_params['block_width']
            i_row = int((curve["pos"] - 1) / obcold_params['n_cols']) - \
                    (current_block_y - 1) * obcold_params['block_height']
            if i_col < 0:
                print('Processing {} failed: i_row={}; i_col={} for {}'.format(filename,
                                                                                  i_row, i_col, dat_pth))
                return

            if method == 'OBCOLD':
                current_dist_type = getcategory_obcold(ccd_block, count, current_dist_type)
            else:
                current_dist_type = getcategory_cold(ccd_block, count)
            break_year = pd.Timestamp.fromordinal(curve['t_break'] - 366).year
            if break_year < year_lowbound or break_year > year_uppbound:
                continue
            results_block[break_year - year_lowbound][i_row][i_col] = current_dist_type * 1000 + curve['t_break'] - \
                (pd.Timestamp.toordinal(datetime.date(break_year, 1, 1)) + bias) + 1
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
            tmp_map_blocks = [np.load(os.path.join(out_path, 'tmp_map_block{}_{}.npy'.format(x+1, year)))
                              for x in range(obcold_params['n_blocks'])]

            results = np.hstack(tmp_map_blocks)
            results = np.vstack(np.hsplit(results, obcold_params['n_block_v']))

            for x in range(obcold_params['n_blocks']):
                os.remove(os.path.join(out_path, 'tmp_map_block{}_{}.npy'.format(x+1, year)))
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
        recent_dist = np.full((obcold_params['n_rows'], obcold_params['n_cols']), 0, dtype=np.int16)
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




if __name__ == '__main__':
    main()
