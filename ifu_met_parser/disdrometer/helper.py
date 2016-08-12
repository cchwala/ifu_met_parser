#----------------------------------------------------------------------------
# Name:         
# Purpose:      
#
# Authors:      
#
# Created:      
# Copyright:    (c) Christian Chwala 2014
# Licence:      The MIT License
#----------------------------------------------------------------------------

import numpy as np


def generate_dsvd_mask(dsvd, D, v, precip_type):
    """Generates a mask for precipitation.

    The mask contains the valid range of a precipitation event.

    Parameters
    -----------

    dsvd : array_like
        array of measured dropsize velocity distribution.
    D : array_like
        array of drop siye classes.
    v : array_like
        array of velocity classes.
    precip_type : string
        type of precipitation for which the mask should be generated,
        e.g. 'rain', 'snow', 'hail', ...

    Returns
    -----------

    dsvd_mask : array_like
        array of boolean with True at indices that correspond  to
        the desired precipitation type

    """

    dsvd_mask = np.zeros_like(dsvd)

    D_grid, v_grid = np.meshgrid(D, v)
    v_t_grid = terminal_velocity(D_grid, approx_type='Maetzler')

    if precip_type == 'rain':
        dsvd_mask = ((v_grid < v_t_grid * 1.5 + 0.5) &
                     (v_grid > v_t_grid * 0.5 - 0.5) &
                     (D_grid < 5.5))

    return dsvd_mask


def terminal_velocity(D=None, approx_type='Maetzler', P=1013):
    """ Calculate the terminal fall velocity of rain drops

    """
    # TODO: Add adjustment for P
    if D is None:
        D = np.arange(0.1,5,0.1)
    if approx_type == 'Gossard':
        v = 9.65 * (1 - np.exp(-0.53 * D))
        return v
    elif approx_type == 'Maetzler':
        v = np.zeros_like(D)
        v[D <= 0.03] = 0
        v[(D > 0.03) & (D <= 0.6)] = 4.323 * (D[(D > 0.03) & (D <= 0.6)] - 0.03)
        v[D > 0.6] = 9.65 - 10.3 * np.exp(-0.6 * D[D > 0.6])
        return v
    elif approx_type == 'Gunn Kinzer data':
        D = np.array([0.0913, 1.687, 4.94, 5.76])
        v = np.array([0.25,   5.90,  9.08,  9.17])
        return v, D
    else:
        raise ValueError('No valid type for the approximation given')