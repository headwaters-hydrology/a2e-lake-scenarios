#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  5 09:40:11 2024

@author: mike
"""
import numpy as np
import math


##################################################
### Parameters

c_n_in_current = 4000
c_n_in_scenario = 800
c_p_in_current = 100
c_p_in_scenario = 40

z_maxes = [7, 20]
t = 10
u = 2
fetch = 20

r_n_in = c_n_in_scenario/c_n_in_current
r_p_in = c_p_in_scenario/c_p_in_current

#################################################
### Functions


def est_b(t):
    """

    """
    return 1 + 0.44*(t**0.13)


def est_log_c_lake_p(c_p_in, t, z_max):
    """

    """
    if z_max > 7.5:
        return np.log10(c_p_in)/est_b(t)
    else:
        return np.log10(c_p_in)


def est_log_c_lake_n(c_n_in, z_max):
    """

    """
    return 1.6 + 0.54*np.log10(c_n_in) - 0.41*np.log10(z_max)


def est_log_c_lake_chla(log_c_lake_p, log_c_lake_n):
    """

    """
    return -1.8 + 0.7*log_c_lake_n + 0.55*log_c_lake_p


def est_d_lake_secchi(log_c_lake_chla, z_max, u, fetch):
    """

    """
    if z_max >= 20:
        return (3.46 - 1.53*log_c_lake_chla)**2
    else:
        return (3.46 - 0.74*log_c_lake_chla - 0.35*np.log10((fetch*(u**2))/z_max))**2


def est_r_lake_p(r_p_in, t, z_max):
    """

    """
    if z_max > 7.5:
        return r_p_in**(1/est_b(t))
    else:
        return r_p_in


def est_r_lake_n(r_n_in):
    """

    """
    return r_n_in**0.54


def est_r_lake_chla(r_p_in, r_n_in):
    """

    """
    return (r_p_in**0.55) * (r_n_in**0.7)


def est_d_lake_secchi_scenario(r_lake_chla, d_secchi_lake_current, z_max):
    """

    """
    if z_max >= 20:
        return (np.log10(r_lake_chla**-1.53) + d_secchi_lake_current**0.5)**2
    else:
        return (np.log10(r_lake_chla**-0.74) + d_secchi_lake_current**0.5)**2



##################################################
### Tests

## TP


def test_tp_lake():
    """

    """
    for z_max in z_maxes:
        c_p_lake_current = 10**est_log_c_lake_p(c_p_in_current, t, z_max)
        c_p_lake_scenario = 10**est_log_c_lake_p(c_p_in_scenario, t, z_max)

        r_p_lake1 = round(c_p_lake_scenario/c_p_lake_current, 3)

        r_p_lake2 = round(est_r_lake_p(r_p_in, t, z_max), 3)

        assert r_p_lake1 == r_p_lake2


def test_tn_lake():
    """

    """
    for z_max in z_maxes:
        c_n_lake_current = 10**est_log_c_lake_n(c_n_in_current, z_max)
        c_n_lake_scenario = 10**est_log_c_lake_n(c_n_in_scenario, z_max)

        r_n_lake1 = round(c_n_lake_scenario/c_n_lake_current, 3)

        r_n_lake2 = round(est_r_lake_n(r_n_in), 3)

        assert r_n_lake1 == r_n_lake2


def test_chla_lake():
    """

    """
    for z_max in z_maxes:
        log_c_p_lake_current = est_log_c_lake_p(c_p_in_current, t, z_max)
        log_c_p_lake_scenario = est_log_c_lake_p(c_p_in_scenario, t, z_max)
        log_c_n_lake_current = est_log_c_lake_n(c_n_in_current, z_max)
        log_c_n_lake_scenario = est_log_c_lake_n(c_n_in_scenario, z_max)

        c_chla_lake_current = 10**est_log_c_lake_chla(log_c_p_lake_current, log_c_n_lake_current)
        c_chla_lake_scenario = 10**est_log_c_lake_chla(log_c_p_lake_scenario, log_c_n_lake_scenario)

        r_chla_lake1 = round(c_chla_lake_scenario/c_chla_lake_current, 3)

        r_p_lake2 = est_r_lake_p(r_p_in, t, z_max)
        r_n_lake2 = est_r_lake_n(r_n_in)
        r_chla_lake2 = round(est_r_lake_chla(r_p_lake2, r_n_lake2), 3)

        assert r_chla_lake1 == r_chla_lake2


def test_secchi_lake():
    """

    """
    for z_max in z_maxes:
        log_c_p_lake_current = est_log_c_lake_p(c_p_in_current, t, z_max)
        log_c_p_lake_scenario = est_log_c_lake_p(c_p_in_scenario, t, z_max)
        log_c_n_lake_current = est_log_c_lake_n(c_n_in_current, z_max)
        log_c_n_lake_scenario = est_log_c_lake_n(c_n_in_scenario, z_max)

        log_c_chla_lake_current = est_log_c_lake_chla(log_c_p_lake_current, log_c_n_lake_current)
        log_c_chla_lake_scenario = est_log_c_lake_chla(log_c_p_lake_scenario, log_c_n_lake_scenario)

        d_secchi_current = est_d_lake_secchi(log_c_chla_lake_current, z_max, u, fetch)
        d_secchi_scenario1 = round(est_d_lake_secchi(log_c_chla_lake_scenario, z_max, u, fetch), 3)

        r_p_lake2 = est_r_lake_p(r_p_in, t, z_max)
        r_n_lake2 = est_r_lake_n(r_n_in)
        r_chla_lake2 = est_r_lake_chla(r_p_lake2, r_n_lake2)

        d_secchi_scenario2 = round(est_d_lake_secchi_scenario(r_chla_lake2, d_secchi_current, z_max), 3)

        assert d_secchi_scenario1 == d_secchi_scenario2
















