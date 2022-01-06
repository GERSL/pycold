import numpy as np


def get_block_y(block_id, n_block_x):
    """
    Parameters
    ----------
    block_id: integer
    n_block_x: integer, number of blocks at x xis

    Returns
    -------
    current block id at y axis
    """
    return int(np.floor((block_id - 1) / n_block_x)) + 1


def get_block_x(block_id, n_block_x):
    """
    Parameters
    ----------
    block_id: integer
    n_block_x: integer, number of blocks at x xis

    Returns
    -------
    current block id at x axis
    """
    return (block_id - 1) % n_block_x + 1
