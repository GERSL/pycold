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


# Register different quality bit standards. (This could/should be moved to a
# different module for for conciceness)
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


def setup_logging():
    # TODO: handle HPC things here in addition to stdout for doctests
    logging.basicConfig(level='INFO')


def grab_demo_kwcoco_dataset():
    """
    Get a demo kwcoco dataset for use in testing

    Returns:
        Path: the path to the kwcoco dataset
    """
    # Register the name of the dataset, how to obtain it, and verify it.
    dataset_info = {
        'name': 'Aligned-DemoKHQ-2022-09-19-V7',
        # The CID is the IPFS Content ID
        'cid': 'bafybeigdkhphpa3n3rdv33w7g6tukmprdnch7g4bp4hc6ebmcr76y6yhwu',
        # The sha512 is the hash of the zipfile we expect to grab.
        'sha512': '6b98195f1d695c695d622ad9debeec93586e518de5934350a17001',
    }
    dpath = ub.Path.appdir('pycold/tests/demodata/kwcoco').ensuredir()

    coco_fpath = dpath / dataset_info['name'] / 'data.kwcoco.json'
    if not coco_fpath.exists():
        # If the data does not already exist
        # Use IPFS to download a demo kwcoco file with LandSat bands
        filename = dataset_info['name'] + '.zip'
        url = f'https://ipfs.io/ipfs/{dataset_info["cid"]}?filename={filename}'
        zip_fpath = ub.download(
            url=url,
            fname=filename,
            hash_prefix=dataset_info['sha512'],
            hasher='sha512')
        # Unzip the data
        import zipfile
        zfile = zipfile.ZipFile(zip_fpath)
        zfile.extractall(dpath)

        if __debug__:
            # The first item should be the name of the bundle folder
            expected_bundle_name = coco_fpath.parent.name
            got_bundle_name = zfile.filelist[0].filename.strip('/')
            if got_bundle_name != expected_bundle_name:
                print(f'expected_bundle_name={expected_bundle_name}')
                print(f'got_bundle_name={got_bundle_name}')

        if not coco_fpath.exists():
            raise AssertionError(
                'The zipfile did not contain the expected dataset')

    return coco_fpath


def _demo_kwcoco_bands():
    """
    This is a quick example illustrating how to dig into a single kwcoco image.

    This does not belong in the main logic and should be moved to an examples
    or tutorial folder. But it can live here for now.
    """
    from pycold.imagetool.prepare_kwcoco import grab_demo_kwcoco_dataset
    coco_fpath = grab_demo_kwcoco_dataset()

    import kwcoco
    # Load the dataset
    coco_dset = kwcoco.CocoDataset(coco_fpath)

    # Grab one arbitrary image id
    image_id = coco_dset.images()[0]

    coco_img = coco_dset.coco_image(image_id)

    # Concept: all important metadata lives in the coco_img.img dictionary
    # But that's messsy to work with, so lets demo the API instead

    # One thing to note is that all file paths will either be absolute
    # or relative to the bundle "dpath".

    print(f'coco_dset.bundle_dpath={coco_dset.bundle_dpath}')

    # We can get a sense of what lives in the coco image by looking at its
    # asset objects.
    for asset_index, asset in enumerate(coco_img.iter_asset_objs()):
        print('+ --- Asset {} --- '.format(asset_index))
        print(f'  * file_name = file_name={asset["file_name"]}')
        print(f'  * channels = file_name={asset["channels"]}')

    # A concicse list of all channels is available here using the kwcoco
    # channel spec.
    print(f'coco_img.channels={coco_img.channels}')

    # We can use the "delay" method to construct a delayed load object and
    # manipluate how we will load the image for a particular set of bands
    # before we actually do it. Lets load two items of data: the RGB and QA
    # bands.

    rgb_delay = coco_img.delay('red|green|blue')
    qa_delay = coco_img.delay('qa_pixel')

    # The finalize method tells the delayed operations to execute.
    rgb_data = rgb_delay.finalize()
    qa_data = qa_delay.finalize()

    # Because the QA band is categorical, we should be able to make a short
    # histogram that describes what is inside.
    qabits_to_count = ub.dict_hist(qa_data.ravel())
    print('qabits_to_count = {}'.format(ub.repr2(qabits_to_count, nl=1)))

    # I'm not sure exactly what these are, so lets try to visualize them next
    # to the RGB data.  First we need to take the uint16 rgb data and make it
    # in a visual range. We can do this robustly with kwimage.normalize_intensity

    import kwimage
    import kwplot
    rgb_canvas = kwimage.normalize_intensity(rgb_data)

    # For the QA band lets assign a color to each category
    colors = kwimage.Color.distinct(len(qabits_to_count))
    qabits_to_color = dict(zip(qabits_to_count, colors))

    # Colorize the QA bands
    colorized = np.empty(qa_data.shape[0:2] + (3,), dtype=np.float32)
    for qabit, color in qabits_to_color.items():
        mask = qa_data[:, :, 0] == qabit
        colorized[mask] = color

    qa_canvas = colorized
    legend = kwplot.make_legend_img(qabits_to_color)  # Make a legend

    # Stack things together into a nice single picture
    qa_canvas = kwimage.stack_images([qa_canvas, legend], axis=1)
    canvas = kwimage.stack_images([rgb_canvas, qa_canvas], axis=1)

    ### Note, I'm not sure if the system the user is has an X server
    ### or if this works in Jupyter notebooks. I will document the way I view
    ### images on my machine in an IPython terminal, but I will also show how
    ### to save these figures to disk so you can rsync them to a machine or do
    ### whatever you normally do to look at an image.

    ### HEADLESS METHOD
    kwimage.imwrite('canvas.png', canvas)

    ### IPYTHON METHOD
    import kwplot
    plt = kwplot.autoplt()

    kwplot.imshow(canvas)

    plt.show()


def stack_kwcoco(coco_fpath, out_dir):
    """
    Args:
        coco_fpath (str | PathLike): path to a kwcoco dataset
        out_dir (str | PathLike): path to write the data

    Example:
        >>> from pycold.imagetool.prepare_kwcoco import *  # NOQA
        >>> setup_logging()
        >>> coco_fpath = grab_demo_kwcoco_dataset()
        >>> dpath = ub.Path.appdir('pycold/tests/stack_kwcoco').ensuredir()
        >>> out_dir = dpath / 'stacked'
        >>> stack_kwcoco(coco_fpath, out_dir)
    """

    # TODO: determine the block settings from the config
    config = {
        'n_block_x': 20,
        'n_block_y': 20,
    }

    # TODO: configure
    out_dir = ub.Path(out_dir)

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
    image_name = coco_image.img.get('name', None)
    video_name = coco_image.video.get('name', None)
    if image_name is None:
        image_name = 'img_{:06d}'.format(coco_image.img['id'])
    if video_name is None:
        video_name = 'vid_{:06d}'.format(coco_image.video['id'])

    video_dpath = (out_dir / video_name).ensuredir()

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
        logger.warn('Not enough clear observations for {}/{}'.format(
            video_name, image_name))
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
            block_dpath = (video_dpath / block_dname).ensuredir()
            block_fpath = block_dpath / (image_name + '.npy')
            np.save(block_fpath, block)
        logger.info('Stacked blocked image {}/{}'.format(video_name, image_name))
    else:
        full_fpath = video_dpath / (image_name + '.npy')
        np.save(full_fpath, data)
        logger.info('Stacked full image {}/{}'.format(video_name, image_name))
    return True
