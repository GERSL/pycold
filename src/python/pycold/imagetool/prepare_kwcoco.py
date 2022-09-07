import kwcoco
import numpy as np
import einops
import functools
import operator
import ubelt as ub
import itertools as it
import logging

logger = logging.getLogger(__name__)


# TODO:
# For each sensor, register the specific bands we are interested in.
# This demo currently assumes landsat8
SENSOR_TO_INFO = {}
SENSOR_TO_INFO['L8'] = {
    'sensor_name': 'Landsat-8',
    'intensity_channels': 'blue|green|red|nir|swir16|swir22|lwir11',
    'quality_channels': 'qa_pixel',
    # I'm not sure how to interpret the quality mask here. I'll need to dig in
    # more to figure out what it is. But it should be whatever is in the
    # landsat-c2ard-bt, landsat-c2ard-sr STAC collections.
    # But the idea is we should associate this data with a quality
    # interpretation to simplify the processing code.
    # The data is being read as 21824, which is "0101010101000000"
    'quality_interpretation': 'FMASK'  # I dont think this is right.
}


# Register different quality bit standards.
QUALITY_BIT_INTERPRETATIONS = {}

# These are specs for TA1 processed data
# https://smartgitlab.com/TE/standards/-/wikis/Data-Output-Specifications#quality-band
QUALITY_BIT_INTERPRETATIONS['TA1'] = {
    'TnE'           : 1 << 0,  # T&E binary mask
    'dilated_cloud' : 1 << 1,
    'cirrus'        : 1 << 2,
    'cloud'         : 1 << 3,
    'cloud_shadow'  : 1 << 4,
    'snow'          : 1 << 5,
    'clear'         : 1 << 6,
    'water'         : 1 << 7,
}

QUALITY_BIT_INTERPRETATIONS['FMASK'] = {
    'clear_land'     : 1 << 0,
    'clear_water'    : 1 << 1,
    'cloud_shadow'   : 1 << 2,
    'snow'           : 1 << 3,
    'cloud'          : 1 << 4,
    'no_observation' : 1 << 5,
}


def grab_demo_kwcoco_dataset():
    """
    Get a demo kwcoco dataset for use in testing

    Returns:
        Path: the path to the kwcoco dataset
    """
    dpath = ub.Path.appdir('pycold/tests/kwcoco/demo0').ensuredir()
    coco_fpath = dpath / 'Aligned-DemoKHQ/data.kwcoco.json'
    if not coco_fpath.exists():
        # If the data does not already exist
        # Use IPFS or Girder to download a demo kwcoco file with LandSat bands
        mirrors = [
            'https://ipfs.io/ipfs/bafybeihzcyzhacjplwmuygxapcwvrs6emfzkkztawex77y54dv2bbmufjq',  # Using IPFS is nice, but slow
            'https://data.kitware.com/api/v1/file/6318d19a11dab8142820733f/download',
        ]
        zip_fpath = ub.download(
            url=mirrors[1],
            fname='Aligned-DemoKHQ.zip',
            hash_prefix='d42b7f1b81940004d88227ae3a5520cef082388838122625e673efe4bd0d6824d4ed',
            hasher='sha512')
        # Unzip the data
        import zipfile
        zfile = zipfile.ZipFile(zip_fpath)
        zfile.extractall(dpath)
    return coco_fpath


def stack_kwcoco(coco_fpath, out_dpath):
    """
    Args:
        coco_fpath (str | PathLike): path to a kwcoco dataset
        out_dpath (str | PathLike): path to write the data

    Example:
        >>> from pycold.imagetool.prepare_kwcoco import *  # NOQA
        >>> coco_fpath = grab_demo_kwcoco_dataset()
        >>> out_dpath = dpath / 'stacked'
        >>> stack_kwcoco(coco_fpath, out_dpath)
    """

    # TODO: determine the block settings from the config
    config = {
        'n_block_x': 20,
        'n_block_y': 20,
    }

    # TODO: configure
    out_dir = ub.Path(out_dpath)

    # Load the kwcoco dataset
    dset = kwcoco.CocoDataset(coco_fpath)

    # A kwcoco dataset can point to multiple images and videos (which are
    # sequences of images). In the future we will likely want to split the main
    # kwcoco file into many smaller ones so we can take advantage of parallel
    # processing, but for now lets just loop over each video and each image in
    # that video.
    videos = dset.videos()

    for video_id in videos:

        # Get the image ids of each image in this video seqeunce
        images = dset.images(video_id=video_id)

        for image_id in images:

            # Construct a CocoImage object which provides access to all
            # underlying metadata and a concicse API to efficiently sample
            # image pixels.
            coco_image : kwcoco.CocoImage = dset.coco_image(image_id)

            # Using detach separates the image from the parent dataset which
            # will allow for easier parallel processing
            coco_image = coco_image.detach()

            # Transform the image data into the desired block structure.
            process_one_coco_image(coco_image, config, out_dir)


def process_one_coco_image(coco_image, config, out_dir):
    """
    Args:
        coco_image (kwcoco.CocoImage): the image to process
        out_dir (Path): path to write the image data
    """
    n_block_x = config['n_block_x']
    n_block_y = config['n_block_y']
    is_partition = True  # hard coded

    # Use the COCO name as a unique filename id.
    file_name = coco_image.img.get('name')

    # Other relevant coco metadata
    coco_image.img['date_captured']
    coco_image.img['frame_index']

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
    quality_bits = QUALITY_BIT_INTERPRETATIONS[quality_interpretation]

    # Specify how we are going to handle spatial resampling and nodata
    delay_kwargs = {
        # Note on the "nodata_method" argument to delay:
        # This argument tells kwcoco how to handle nodata values in the
        # underlying imagery.
        # Using None will tell kwcoco to return nodata value as-is.
        # Using 'float' tells kwcoco to convert all dtypes to float32 and
        # use nan to fill nodata values.
        # Using 'ma' will tell kwcoco to return numpy masked arrays
        # over nodata values.
        # Using ma would be best in this case, but it currently isn't
        # as tested as float.
        # 'nodata_method': 'ma',
        # 'nodata_method': 'float',
        'nodata_method': None,

        # Note on the "space" argument to delay:
        # By specifying "Video Space" kwcoco will ensure all images in
        # the sequence are pixel aligned
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

    if 0:
        # Developer note:
        # Try running this to get a feel for what the delayed image
        # representation. Also try removing the optimize call to see what
        # changes!
        delayed_im.optimize().write_network_text()
        delayed_qa.optimize().write_network_text()

    # It is important that the categorical QA band is not interpolated or
    # antialiased, whereas the intensity bands should be.
    qa_data = delayed_qa.finalize(interpolation='nearest', antialias=False)

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

    if clear_ratio <= clear_threshold:
        logger.warn('Not enough clear observations for {}'.format(file_name))
        return False

    im_data = delayed_im.finalize(interpolation='cubic', antialias=True)

    # NOTE: if we enable a nodata method, we will need to handle it here.
    # NOTE: if any intensity modification needs to be done handle it here.

    data = np.concatenate([im_data, qa_data], axis=2)

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

            block_dname = 'block_x{}_y{}'.format(j + 1, i + 1)
            block_dpath = (out_dir / block_dname).ensuredir()
            np.save(block_dpath / file_name, block)

    else:
        np.save(out_dir / file_name, blocks)
    return True
