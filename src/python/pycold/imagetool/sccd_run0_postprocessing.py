# this script is needed to be ran after run 0 is finished
import yaml
import os
import numpy as np
from datetime import datetime
from pytz import timezone
import click
import time
from pycold.utils import assemble_array
from pycold.pyclassifier import PyClassifierHPC
from pycold.app import defaults
from os.path import join
from pycold.imagetool.TileProcessing import phen_anchor_days

training_year = 2019


@click.command()
@click.option('--rank', type=int, default=0, help='the rank id')
@click.option('--n_cores', type=int, default=0, help='the total cores assigned')
@click.option('--sccdpack_path', type=str, default=None, help='the path for storing results')
@click.option('--yaml_path', type=str, default=None, help='YAML path')
@click.option('--seedmap_path', type=str, default=None, help='an existing label map path; '
                                                             'none means not using thematic info')
def main(rank, n_cores, sccdpack_path, yaml_path, seedmap_path):
    tz = timezone('US/Eastern')
    # Reading config
    with open(yaml_path, 'r') as yaml_obj:
        config = yaml.safe_load(yaml_obj)

    # set up some additional config
    block_width = int(config['n_cols'] / config['n_block_x'])  # width of a block
    block_height = int(config['n_rows'] / config['n_block_y'])  # height of a block
    nblock_eachcore = int(np.ceil(config['n_block_x'] * config['n_block_y'] * 1.0 / n_cores))
    if (config['n_cols'] % block_width != 0) or (config['n_rows'] % block_height != 0):
        print('n_cols, n_rows must be divisible respectively by block_width, block_height! Please double '
              'check your config yaml')
        exit()

    pyclassifier = PyClassifierHPC(config, record_path=sccdpack_path, n_features_perband=defaults['SCCD']['N_FEATURES'],
                                   band_num=defaults['SCCD']['NRT_BAND'], year_list_to_predict=[training_year],
                                   tmp_path=sccdpack_path,
                                   output_path=sccdpack_path,
                                   seedmap_path=seedmap_path)

    if not pyclassifier.is_finished_step4_assemble():
        if rank == 1:
            # 1) output rf.model
            pyclassifier.step2_train_rf(ref_year=training_year)
            for day in phen_anchor_days:
                pyclassifier.step2_train_rf(ref_year=day)

            print("Training rf ends: {}".format(datetime.now(tz).strftime('%Y-%m-%d %H:%M:%S')))

        for i in range(nblock_eachcore):
            if n_cores * i + rank > config['n_block_x'] * config['n_block_y']:
                break
            pyclassifier.step3_classification_sccd(block_id=n_cores * i + rank)

        if rank == 1:  # serial mode for assemble
            # 2) output yearly classification map
            pyclassifier.step4_assemble_sccd(clean=False)

    while not pyclassifier.is_finished_step4_assemble():
        time.sleep(15)

    if rank == 1:
        print("Assemble classification map ends: {}".format(datetime.now(tz).strftime('%Y-%m-%d %H:%M:%S')))
        # 3) output last change date map
        # lastchangedate_assemble = assemble_array(tmp_map_blocks, config['n_block_x'])
        lastchangedate_tile = np.full((config['n_rows'], config['n_cols']), 0)
        np.save(os.path.join(sccdpack_path, 'last_change_date.npy'), lastchangedate_tile)

        # 4) output object map
        object_map = np.full((config['n_rows'], config['n_cols']), 0)
        np.save(os.path.join(sccdpack_path, 'accumulated_oid.npy'), object_map)

        # clean
        tmp_filenames = [file for file in os.listdir(sccdpack_path)
                         if file.startswith('tmp_')]
        for file in tmp_filenames:
            os.remove(join(sccdpack_path, file))


if __name__ == '__main__':
    main()
