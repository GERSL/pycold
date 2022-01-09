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
    return int((block_id - 1) / n_block_x) + 1


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


def get_col_index(pos, n_cols, current_block_x, block_width):
    """
    Parameters
    ----------
    pos
    n_cols
    current_block_x
    block_width

    Returns
    -------

    """
    return int((pos - 1) % n_cols) - (current_block_x - 1) * block_width


def get_row_index(pos, n_cols, current_block_y, block_height):
    """
    Parameters
    ----------
    pos: start from 1
    n_cols
    current_block_y
    block_height

    Returns
    -------

    """
    return int((pos - 1) / n_cols) - (current_block_y - 1) * block_height

