"""
This is a proof-of-concept for converting kwcoco files into the
expected data structures for pycold.

Relevant functions:
    * grab_demo_kwcoco_dataset - downloads a small kwcoco dataset for testing
    * stack_kwcoco - runs the stacking process on an entire kwcoco file
    * process_one_coco_image - runs the stacking for a single coco image.

Limitations:
    * Currently only handles Landsat-8

    * The quality bands are not exactly what I was expecting them to be,
      some of the quality filtering is stubbed out or disabled.

    * Not setup for an HPC environment yet, but that extension shouldn't be too
      hard.

    * Nodata values are currently not masked or handled

    * Configurations are hard-coded
"""
import kwcoco
import json
import numpy as np
import einops
import functools
import operator
import ubelt as ub
import itertools as it
import logging
from datetime import datetime
import numpy as geek
logger = logging.getLogger(__name__)


# TODO:
# For each sensor, register the specific bands we are interested in.
# This demo currently assumes landsat8
SENSOR_TO_INFO = {}
SENSOR_TO_INFO['L8'] = {
    'sensor_name': 'Landsat-8',
    'intensity_channels': 'blue|green|red|nir|swir16|swir22|lwir11',
    'quality_channels': 'cloudmask',
    'quality_interpretation': 'FMASK'
    } # The name of quality_channels for Drop 4 is 'cloudmask'.

# Register different quality bit standards.
QA_INTERPRETATIONS = {}

# These are specs for TA1 processed data
QA_BIT = {
    'clear'         : 1 << 0,
    'cloud'         : 1 << 1,
    'cloud_adj'     : 1 << 2,
    'shadow'        : 1 << 3,
    'snow'          : 1 << 4,
    'water'         : 1 << 5,
}

QA_INTERPRETATIONS['FMASK'] = {
    'clear'         : 0,
    'water'         : 1,
    'cloud_shadow'  : 2,
    'snow'          : 3,
    'cloud'         : 4,
    'no_obs'        : 255,
}
QUALITY_BIT_INTERPRETATIONS = {}

# function for decoding HLS qa band
def qa_decoding(qa_array):
    """
    This function is modified from qabitval_array_HLS function
    (https://github.com/GERSL/pycold/blob/c5b380eccc2916e5c3aec0bbd2b1982e114b75b1/src/python/pycold/imagetool/prepare_ard.py#L74)
    """
    unpacked = np.full(qa_array.shape, QA_INTERPRETATIONS['FMASK']['clear'])

    QA_CLOUD_unpacked = geek.bitwise_and(qa_array, QA_BIT['cloud'])
    QA_CLOUD_ADJ = geek.bitwise_and(qa_array, QA_BIT['cloud_adj'])
    QA_SHADOW_unpacked = geek.bitwise_and(qa_array, QA_BIT['shadow'])
    QA_SNOW_unpacked = geek.bitwise_and(qa_array, QA_BIT['snow'])
    QA_WATER_unpacked = geek.bitwise_and(qa_array, QA_BIT['water'])

    unpacked[QA_WATER_unpacked > 0] = QA_INTERPRETATIONS['FMASK']['water']
    unpacked[QA_SNOW_unpacked > 0] = QA_INTERPRETATIONS['FMASK']['snow']
    unpacked[QA_SHADOW_unpacked > 0] = QA_INTERPRETATIONS['FMASK']['cloud_shadow']
    unpacked[QA_CLOUD_ADJ > 0] = QA_INTERPRETATIONS['FMASK']['cloud']
    unpacked[QA_CLOUD_unpacked > 0] = QA_INTERPRETATIONS['FMASK']['cloud']
    unpacked[qa_array == QA_INTERPRETATIONS['FMASK']['no_obs']] = QA_INTERPRETATIONS['FMASK']['no_obs']

    return unpacked

def qa_decoding_no_boundary(qa_array):
    """
    This function is modified from qabitval_array_HLS function
    (https://github.com/GERSL/pycold/blob/c5b380eccc2916e5c3aec0bbd2b1982e114b75b1/src/python/pycold/imagetool/prepare_ard.py#L74)
    """
    unpacked = np.full(qa_array.shape, QA_INTERPRETATIONS['FMASK']['clear'])

    QA_CLOUD_unpacked = geek.bitwise_and(qa_array, QA_BIT['cloud'])
    QA_CLOUD_ADJ = geek.bitwise_and(qa_array, QA_BIT['cloud_adj'])
    QA_SHADOW_unpacked = geek.bitwise_and(qa_array, QA_BIT['shadow'])
    QA_SNOW_unpacked = geek.bitwise_and(qa_array, QA_BIT['snow'])
    QA_WATER_unpacked = geek.bitwise_and(qa_array, QA_BIT['water'])

    unpacked[QA_WATER_unpacked > 0] = QA_INTERPRETATIONS['FMASK']['water']
    unpacked[QA_SNOW_unpacked > 0] = QA_INTERPRETATIONS['FMASK']['snow']
    unpacked[QA_SHADOW_unpacked > 0] = QA_INTERPRETATIONS['FMASK']['cloud_shadow']
    unpacked[QA_CLOUD_unpacked > 0] = QA_INTERPRETATIONS['FMASK']['cloud']
    unpacked[QA_CLOUD_ADJ > 0] = QA_INTERPRETATIONS['FMASK']['clear']
    unpacked[qa_array == QA_INTERPRETATIONS['FMASK']['no_obs']] = QA_INTERPRETATIONS['FMASK']['no_obs']

    return unpacked

def setup_logging():
    # TODO: handle HPC things here in addition to stdout for doctests
    logging.basicConfig(level='INFO')

####################################################################################
#                   Functions for Artificial Surface Index (ASI)                   #
#  See original code: https://github.com/GERSL/ASI_py/blob/main/ASI_standalone.py  #
####################################################################################

def hist_cut(band, mask, fill_value=-9999, k=3, minmax='std'):
    if minmax == 'std':
        mean = band[mask].mean()
        std = band[mask].std()
        low_val = (mean - k * std)
        high_val = (mean + k * std)
    else:
        low_val, high_val = minmax # use specified value range.
    is_low = band < low_val
    is_high = band > high_val
    mask_invalid_index = is_low | is_high
    band[mask_invalid_index] = fill_value
    return band, ~mask_invalid_index


def minmax_norm(band, mask, fill_value=-9999):
    max_val = band[mask].max()
    min_val = band[mask].min()
    extent = max_val - min_val
    if extent != 0:
        shifted = band - min_val
        scaled = shifted / extent
        band[mask] = scaled[mask]    
    band[~mask] = fill_value
    return band

# Artificial Surface Index (ASI) is designed based the surface reflectance imagery of Landsat 8.
def artificial_surface_index(Blue, Green, Red, NIR, SWIR1, SWIR2, Scale, MaskValid_Obs, fillV):
    ##### The calculation chain.

    # Artificial surface Factor (AF).
    AF = (NIR - Blue) / (NIR + Blue) + 0.000001
    AF, MaskValid_AF = hist_cut(AF, MaskValid_Obs, fillV, 6, [-1, 1])
    MaskValid_AF_U = MaskValid_AF & MaskValid_Obs
    AF_Norm = minmax_norm(AF, MaskValid_AF_U, fillV)

    # Vegetation Suppressing Factor (VSF).
    MSAVI = ( (2*NIR+1*Scale) - np.sqrt((2*NIR+1*Scale)**2 - 8*(NIR-Red)) ) / 2 # Modified Soil Adjusted Vegetation Index (MSAVI).
    MSAVI, MaskValid_MSAVI = hist_cut( MSAVI, MaskValid_Obs, fillV, 6, [-1, 1])
    NDVI = (NIR - Red) / (NIR + Red) + 0.000001
    NDVI, MaskValid_NDVI  = hist_cut(NDVI, MaskValid_Obs, fillV, 6, [-1, 1])
    VSF = 1 - MSAVI*NDVI
    MaskValid_VSF = MaskValid_MSAVI & MaskValid_NDVI & MaskValid_Obs
    VSF_Norm = minmax_norm(VSF, MaskValid_VSF, fillV)

    # Soil Suppressing Factor (SSF).
    # Derive the Modified Bare soil Index (MBI).
    MBI = (SWIR1 - SWIR2 - NIR) / (SWIR1 + SWIR2 + NIR) + 0.5
    MBI, MaskValid_MBI = hist_cut(MBI, MaskValid_Obs, fillV, 6, [-0.5, 1.5])
    # Deriving Enhanced-MBI based on MBI and MNDWI.
    MNDWI = (Green - SWIR1) / (Green + SWIR1) + 0.000001
    MNDWI, MaskValid_MNDWI = hist_cut(MNDWI, MaskValid_Obs, fillV, 6, [-1, 1])
    EMBI = ((MBI+0.5) - (MNDWI+1)) / ((MBI+0.5) + (MNDWI+1))
    EMBI, MaskValid_EMBI = hist_cut(EMBI, MaskValid_Obs, fillV, 6, [-1, 1])
    # Derive SSF.
    SSF = (1 - EMBI)
    MaskValid_SSF = MaskValid_MBI & MaskValid_MNDWI & MaskValid_EMBI & MaskValid_Obs
    SSF_Norm = minmax_norm(SSF, MaskValid_SSF, fillV)

    # Modulation Factor (MF).
    MF = (Blue + Green - NIR - SWIR1) / (Blue + Green + NIR + SWIR1) + 0.000001
    MF, MaskValid_MF = hist_cut(MF, MaskValid_Obs, fillV, 6, [-1, 1])
    MaskValid_MF_U = MaskValid_MF & MaskValid_Obs
    MF_Norm = minmax_norm(MF, MaskValid_MF_U, fillV)

    # Derive Artificial Surface Index (ASI).
    ASI = AF_Norm * SSF_Norm * VSF_Norm * MF_Norm
    MaskValid_ASI = MaskValid_AF_U & MaskValid_VSF & MaskValid_SSF & MaskValid_MF_U & MaskValid_Obs    
    ASI[~MaskValid_ASI] = fillV
    
    return ASI
    
# def grab_demo_kwcoco_dataset():
#     """
#     Get a demo kwcoco dataset for use in testing
#
#     Returns:
#         Path: the path to the kwcoco dataset
#     """
#     # Register the name of the dataset, how to obtain it, and verify it.
#     dataset_info = {
#         'name': 'Aligned-DemoKHQ-2022-09-19-V7',
#         # The CID is the IPFS Content ID
#         'cid': 'bafybeigdkhphpa3n3rdv33w7g6tukmprdnch7g4bp4hc6ebmcr76y6yhwu',
#         # The sha512 is the hash of the zipfile we expect to grab.
#         'sha512': '6b98195f1d695c695d622ad9debeec93586e518de5934350a17001',
#     }
#     dpath = ub.Path.appdir('pycold/tests/demodata/kwcoco').ensuredir()
#
#     coco_fpath = dpath / dataset_info['name'] / 'data.kwcoco.json'
#     if not coco_fpath.exists():
#         # If the data does not already exist
#         # Use IPFS to download a demo kwcoco file with LandSat bands
#         filename = dataset_info['name'] + '.zip'
#         url = f'https://ipfs.io/ipfs/{dataset_info["cid"]}?filename={filename}'
#         zip_fpath = ub.download(
#             url=url,
#             fname=filename,
#             hash_prefix=dataset_info['sha512'],
#             hasher='sha512')
#         # Unzip the data
#         import zipfile
#         zfile = zipfile.ZipFile(zip_fpath)
#         zfile.extractall(dpath)
#
#         if __debug__:
#             # The first item should be the name of the bundle folder
#             expected_bundle_name = coco_fpath.parent.name
#             got_bundle_name = zfile.filelist[0].filename.strip('/')
#             if got_bundle_name != expected_bundle_name:
#                 print(f'expected_bundle_name={expected_bundle_name}')
#                 print(f'got_bundle_name={got_bundle_name}')
#
#         if not coco_fpath.exists():
#             raise AssertionError(
#                 'The zipfile did not contain the expected dataset')
#
#     return coco_fpath
#
#
# def _demo_kwcoco_bands():
#     """
#     This is a quick example illustrating how to dig into a single kwcoco image.
#
#     This does not belong in the main logic and should be moved to an examples
#     or tutorial folder. But it can live here for now.
#     """
#     from pycold.imagetool.prepare_kwcoco import grab_demo_kwcoco_dataset
#     coco_fpath = grab_demo_kwcoco_dataset()
#
#     import kwcoco
#     import kwimage
#     import kwplot
#     # Load the dataset
#     coco_dset = kwcoco.CocoDataset(coco_fpath)
#
#     # Grab one arbitrary image id
#     image_id = coco_dset.images()[0]
#
#     coco_img = coco_dset.coco_image(image_id)
#
#     # Concept: all important metadata lives in the coco_img.img dictionary
#     # But that's messsy to work with, so lets demo the API instead
#
#     # One thing to note is that all file paths will either be absolute
#     # or relative to the bundle "dpath".
#
#     print(f'coco_dset.bundle_dpath={coco_dset.bundle_dpath}')
#
#     # We can get a sense of what lives in the coco image by looking at its
#     # asset objects.
#     for asset_index, asset in enumerate(coco_img.iter_asset_objs()):
#         print('+ --- Asset {} --- '.format(asset_index))
#         print(f'  * file_name = {asset["file_name"]}')
#         print(f'  * channels = {asset["channels"]}')
#         print(f'  * parent_file_name= {asset["parent_file_name"]}')
#
#     # A concicse list of all channels is available here using the kwcoco
#     # channel spec.
#     print(f'coco_img.channels={coco_img.channels}')
#
#     # We can use the "delay" method to construct a delayed load object and
#     # manipluate how we will load the image for a particular set of bands
#     # before we actually do it. Lets load two items of data: the RGB and QA
#     # bands.
#
#     rgb_delay = coco_img.delay('red|green|blue')
#     qa_delay = coco_img.delay('qa_pixel')
#
#     # The finalize method tells the delayed operations to execute.
#     rgb_data = rgb_delay.finalize()
#     qa_data = qa_delay.finalize()
#
#     rgb_canvas = kwimage.normalize_intensity(rgb_data)
#
#     # Because the QA band is categorical, we should be able to make a short
#     # histogram that describes what is inside.
#     qabits_to_count = ub.dict_hist(qa_data.ravel())
#
#     # For the QA band lets assign a color to each category
#     colors = kwimage.Color.distinct(len(qabits_to_count))
#     qabits_to_color = dict(zip(qabits_to_count, colors))
#
#     # Colorize the QA bands
#     colorized = np.empty(qa_data.shape[0:2] + (3,), dtype=np.float32)
#     for qabit, color in qabits_to_color.items():
#         mask = qa_data[:, :, 0] == qabit
#         colorized[mask] = color
#
#     rgb_canvas = kwimage.normalize_intensity(rgb_data)
#
#     # Because the QA band is categorical, we should be able to make a short
#
#     qa_canvas = colorized
#     legend = kwplot.make_legend_img(qabits_to_color)  # Make a legend
#
#     # Stack things together into a nice single picture
#     qa_canvas = kwimage.stack_images([qa_canvas, legend], axis=1)
#     canvas = kwimage.stack_images([rgb_canvas, qa_canvas], axis=1)
#
#     ### Note, I'm not sure if the system the user is has an X server
#     ### or if this works in Jupyter notebooks. I will document the way I view
#     ### images on my machine in an IPython terminal, but I will also show how
#     ### to save these figures to disk so you can rsync them to a machine or do
#     ### whatever you normally do to look at an image.
#
#     ### HEADLESS METHOD
#     kwimage.imwrite('canvas.png', canvas)
#
#     ### IPYTHON METHOD
#     import kwplot
#     plt = kwplot.autoplt()
#
#     kwplot.imshow(canvas)
#
#     plt.show()


def stack_kwcoco(coco_fpath, out_dir):
    """
    Args:
        coco_fpath (str | PathLike | CocoDataset):
            the kwcoco dataset to convert

        out_dir (str | PathLike): path to write the data

    Returns:
        List[Dict]: a list of dictionary result objects

    Example:
        >>> from pycold.imagetool.prepare_kwcoco import *  # NOQA
        >>> setup_logging()
        >>> coco_fpath = grab_demo_kwcoco_dataset()
        >>> dpath = ub.Path.appdir('pycold/tests/stack_kwcoco').ensuredir()
        >>> out_dir = dpath / 'stacked'
        >>> results = stack_kwcoco(coco_fpath, out_dir)
    """

    # TODO: determine the block settings from the config
    config = {
        'n_block_x': 20,
        'n_block_y': 20,
        'adj_cloud': False,
        'mode'     : 'ASI' #None # 'ASI'
    }

    # TODO: configure
    out_dir = ub.Path(out_dir)

    # Load the kwcoco dataset
    dset = kwcoco.CocoDataset.coerce(coco_fpath)
    videos = dset.videos()
    results = []

    for video_id in videos:

        if video_id == 17: # testing for US_C000 site
            # Get the image ids of each image in this video seqeunce
            images = dset.images(video_id=video_id)

            for image_id in images:
                coco_image : kwcoco.CocoImage = dset.coco_image(image_id)
                coco_image = coco_image.detach()

                # For now, it supports only L8
                if coco_image.img['sensor_coarse'] == 'L8':
                    adj_cloud = False
                    # Transform the image data into the desired block structure.
                    result = process_one_coco_image(coco_image, config, out_dir)
                    results.append(result)

    return results

def process_one_coco_image(coco_image, config, out_dir):
    """
    Args:
        coco_image (kwcoco.CocoImage): the image to process
        out_dir (Path): path to write the image data

    Returns:
        Dict: result dictionary with keys:
            status (str) : either a string passed or failed
            fpaths (List[str]): a list of files that were written
    """
    n_block_x = config['n_block_x']
    n_block_y = config['n_block_y']
    adj_cloud = config['adj_cloud']
    mode      = config['mode']
    
    is_partition = True  # hard coded

    # Use the COCO name as a unique filename id.
    image_name = coco_image.img.get('name', None)
    video_name = coco_image.video.get('name', None)
    if image_name is None:
        image_name = 'img_{:06d}'.format(coco_image.img['id'])
    if video_name is None:
        video_name = 'vid_{:06d}'.format(coco_image.video['id'])

    video_dpath = (out_dir / video_name).ensuredir()

    # Other relevant coco metadata
    date_captured = coco_image.img['date_captured']
    ordinal_date = datetime.strptime(date_captured[:10], '%Y-%m-%d').toordinal()
    frame_index = coco_image.img['frame_index']
    n_cols = coco_image.img['width']
    n_rows = coco_image.img['height']
    # Determine what sensor the image is from.
    # Note: if kwcoco needs to register more fine-grained sensor
    # information we can do that.
    sensor = coco_image.img['sensor_coarse']
    assert sensor == 'L8', 'MWE only supports landsat-8 for now'

    # Given the sensor, determine what the intensity and quality band
    # we should request are.
    sensor_info = SENSOR_TO_INFO[sensor]
    intensity_channels = sensor_info['intensity_channels']
    quality_channels = sensor_info['quality_channels']
    quality_interpretation = sensor_info['quality_interpretation']
    quality_bits = QA_INTERPRETATIONS[quality_interpretation]
    # Specify how we are going to handle spatial resampling and nodata
    delay_kwargs = {
        'nodata_method': None,
        'space': 'video',
    }

    # Construct delayed images. These represent a tree of image
    # operations that will resample the image at the desired resolution
    # as well as align it with other images in the sequence.
    delayed_im = coco_image.delay(channels=intensity_channels, **delay_kwargs)
    delayed_qa = coco_image.delay(channels=quality_channels, **delay_kwargs)

    # Check what shape the data would be loaded with if we finalized right now.
    h, w = delayed_im.shape[0:2]
    # Determine if padding is necessary to properly break the data into blocks.
    padded_w = int(np.ceil(w / n_block_x) * n_block_x)
    padded_h = int(np.ceil(h / n_block_y) * n_block_y)

    if padded_w != h or padded_h != h:
        # cropping using an oversized slice with clip=False and wrap=False is
        # equivalent to padding. In the future a more efficient pad operation
        # where the padding value can be specified will be added, but this will
        # work well enough for now.
        slice_ = (slice(0, padded_h), slice(0, padded_w))
        delayed_im = delayed_im.crop(slice_, clip=False, wrap=False)
        delayed_qa = delayed_qa.crop(slice_, clip=False, wrap=False)

    # It is important that the categorical QA band is not interpolated or
    # antialiased, whereas the intensity bands should be.
    qa_data = delayed_qa.finalize(interpolation='nearest', antialias=False)
    # Decoding QA band
    adj_cloud = False
    if adj_cloud == True:
        qa_unpacked = qa_decoding(qa_data)
    else:
        qa_unpacked = qa_decoding_no_boundary(qa_data)

    # First check the quality bands before loading all of the image data.
    # FIXME: the quality bits in this example are wrong.
    # Setting the threshold to zero to bypass for now.
    clear_threshold = 0
    if clear_threshold > 0:
        clear_bits = functools.reduce(
            operator.or_, ub.take(quality_bits, ['clear_land', 'clear_water']))
        noobs_bits = functools.reduce(
            operator.or_, ub.take(quality_bits, ['no_observation']))
        is_clear = (qa_data & clear_bits) > 0
        is_noobs = (qa_data & noobs_bits) > 0
        is_obs = ~is_noobs
        is_obs_clear = is_clear & is_obs
        clear_ratio = is_obs_clear.sum() / is_obs.sum()
    else:
        clear_ratio = 1

    result = {
        'status': None
    }

    if clear_ratio <= clear_threshold:
        logger.warn('Not enough clear observations for {}/{}'.format(
            video_name, image_name))
        result['status'] = 'failed'
        return result

    im_data = delayed_im.finalize(interpolation='cubic', antialias=True)

    # NOTE: if we enable a nodata method, we will need to handle it here.
    # NOTE: if any intensity modification needs to be done handle it here.
    
    if mode == 'ASI':
        Scale = 10000
        fill_value = 0
        B1 = im_data[:, :, 0]
        B2 = im_data[:, :, 1]
        B3 = im_data[:, :, 2]
        B4 = im_data[:, :, 3]
        B5 = im_data[:, :, 4]
        B6 = im_data[:, :, 5]
        MaskValid_Obs = ((B1>0) & (B1<1*Scale) &
                        (B2>0) & (B2<1*Scale) &
                        (B3>0) & (B3<1*Scale) &
                        (B4>0) & (B4<1*Scale) &
                        (B5>0) & (B5<1*Scale) &
                        (B6>0) & (B6<1*Scale)
                        )
        
        # Calculating ASI
        ASI = artificial_surface_index(B1.astype(np.float32), B2.astype(np.float32), B3.astype(np.float32), B4.astype(np.float32), B5.astype(np.float32), B6.astype(np.float32), Scale, MaskValid_Obs, fill_value)
        # Get land mask.        
        MNDWI = (B2 - B5) / (B2 + B5)
        MNDWI, MaskValid_MNDWI = hist_cut(MNDWI, MaskValid_Obs, fill_value, 6, [-1, 1])
        Water_Th = 0; # Water threshold for MNDWI (may need to be adjusted for different study areas).
        MaskLand = (MNDWI<Water_Th)

        # Convert dtype from float32 to int16
        ASI = ASI * Scale
        ASI = ASI.astype('int16')
        ASI[ASI == 0] = fill_value     

        # Exclude water pixels.
        ASI[~MaskLand] = fill_value
        ASI = ASI.reshape(ASI.shape[0], ASI.shape[1], 1)
        false_band = np.full((ASI.shape[0], ASI.shape[1], 1), 0)
        # input for Hybrid-COLD (with ASI) = B2, B3, B4, B5, B6, ASI
        data = np.concatenate([im_data[:,:,1:6], ASI, false_band, qa_unpacked], axis=2)

    else:
        data = np.concatenate([im_data, qa_unpacked], axis=2)
    
    result_fpaths = []

    # TODO:
    # save necessary metadata alongside the npy file so we don't have
    # to rely on file names.
    metadata = {
        'image_name' : image_name,
        'date_captured': date_captured,
        'ordinal_date': ordinal_date,
        'n_cols': n_cols,
        'n_rows': n_rows,
        'padded_n_cols': padded_w,
        'padded_n_rows': padded_h,
        'n_block_x': n_block_x,
        'n_block_y': n_block_y,
        'adj_cloud': adj_cloud,
        'mode': mode
    }

    if is_partition:
        bw = int(padded_w / n_block_x)  # width of a block
        bh = int(padded_h / n_block_y)  # height of a block

        # Use einops to rearrange the data into blocks
        # Question: would using numpy strided tricks be faster?
        blocks = einops.rearrange(
            data, '(nby bh) (nbx bw) c -> nbx nby bh bw c', bw=bw, bh=bh)

        for i, j in it.product(range(n_block_y), range(n_block_x)):
            block = blocks[i, j]

            # FIXME: Disable skipping until QA bands are handled correctly
            SKIP_BLOCKS_WITH_QA = False
            if SKIP_BLOCKS_WITH_QA:
                # check if no valid pixels in the chip, then eliminate
                qa_unique = np.unique(block[..., -1])
                qa_unique
                # skip blocks are all cloud, shadow or filled values in DHTC,
                # we also don't need to save pixel that has qa value of
                # 'QA_CLOUD', 'QA_SHADOW', or FILLED value (255)
                if ... and False:
                    continue

            block_dname = 'block_x{}_y{}'.format(i + 1, j + 1)
            block_dpath = (video_dpath / block_dname).ensuredir()
            block_fpath = block_dpath / (image_name + '.npy')

            metadata.update({
                'x': i + 1,
                'y': j + 1,
                'total_pixels': int(np.prod(block.shape[0:2])),
                'total_bands': int(block.shape[-1]),
            })
            meta_fpath = block_dpath / (image_name + '.json')
            meta_fpath.write_text(json.dumps(metadata))
            np.save(block_fpath, block)
            result_fpaths.append(block_fpath)
            result_fpaths.append(meta_fpath)
        logger.info('Stacked blocked image {}/{}'.format(video_name, image_name))
    else:
        metadata.update({
            'total_pixels': int(np.prod(data.shape[0:2])),
            'total_bands': int(data.shape[-1]),
        })
        full_fpath = video_dpath / (image_name + '.npy')
        meta_fpath = video_dpath / (image_name + '.json')
        meta_fpath.write_text(json.dumps(metadata))
        np.save(full_fpath, data)
        result_fpaths.append(full_fpath)
        result_fpaths.append(meta_fpath)
        logger.info('Stacked full image {}/{}'.format(video_name, image_name))

    result['status'] = 'passed'
    result['fpaths'] = result_fpaths
    return result

# FIXME: I wasn't sure how to update grab_demo_kwcoco_dataset(). So I manually defined coco_fpath...
coco_fpath = '/home/jws18003/data/dvc-repos/smart_data_dvc/Aligned-Drop4-2022-08-08-TA1-S2-L8-ACC/data.kwcoco.json'
dpath = ub.Path.appdir('/gpfs/scratchfs1/zhz18039/jws18003/kwcoco').ensuredir()
out_dir = dpath / 'stacked'

stack_kwcoco(coco_fpath, out_dir)