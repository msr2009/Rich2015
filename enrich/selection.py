from __future__ import print_function
from enrich_error import EnrichError
from scipy import stats
from seqlib.basic import BasicSeqLib
from seqlib.barcodevariant import BarcodeVariantSeqLib, BarcodeMap
from seqlib.barcode import BarcodeSeqLib
from seqlib.overlap import OverlapSeqLib
from seqlib.variant import WILD_TYPE_VARIANT
from config_check import seqlib_type
from datacontainer import DataContainer
import os
import re
import math
import itertools
import time
import pandas as pd
import numpy as np
from collections import Counter
import logging
import copy

from sys import stderr

def nonsense_ns_carryover_apply_fn(row, position):
    """
    :py:meth:`pandas.DataFrame.apply` function for determining which rows 
    contribute counts to nonspecific carryover calculations. Returns ``True`` 
    if the variant has a change to stop at or before amino acid number 
    *position*.
    """
    m = re.search("p\.[A-Z][a-z][a-z](\d+)Ter", row.name)
    if m is not None:
        if int(m.group(1)) <= position:
            return True
        else:
            return False
    else:
        return False


def barcode_variation_apply_fn(row, barcode_data, mapping):
    """
    :py:meth:`pandas.DataFrame.apply` function for calculating the coefficient 
    of variation for a variant's barcodes.
    """
    bc_scores = barcode_data.ix[mapping.variants[row.name]]['score']
    bc_scores = bc_scores[~np.isnan(bc_scores)]
    if len(bc_scores) > 0:
        cv = stats.variation(bc_scores)
    else:
        cv = float("NaN")
    return pd.Series({'scored.unique.barcodes' : len(bc_scores), \
                      'barcode.cv' : cv})


def barcode_count_apply_fn(row, mapping):
    """
    :py:meth:`pandas.DataFrame.apply` function for counting the number of 
    unique barcodes for a variant.
    """
    return len(mapping.variants[row.name])


def barcode_varation_filter(row, cutoff):
    """
    Filtering function for barcode coefficient of variation.
    """
    if row['barcode.cv'] > cutoff:
        return False
    else:
        return True


def min_count_filter(row, cutoff):
    """
    Filtering function for minimum counts across all timepoints.
    """
    counts = row[[x for x in row.index if x.startswith("count")]].values
    if counts.min() < cutoff:
        return False
    else:
        return True


def min_input_count_filter(row, cutoff):
    """
    Filtering function for minimum count in input timepoint.
    """
    if row['count.0'] < cutoff:
        return False
    else:
        return True


def min_rsq_filter(row, cutoff):
    """
    Filtering function for minimum r-squared value. Entries with no 
    r-squared value are retained.
    """
    if np.isnan(row['r_sq']):
        return True
    elif row['r_sq'] < cutoff:
        return False
    else:
        return True


def linear_enrichment_apply_fn(row, timepoints):
    """
    :py:meth:`pandas.DataFrame.apply` apply function for calculating 
    enrichment scores and r-squared values using a linear regression.
    """
    if math.isnan(row[0]):
        # not present in input library
        score = float("NaN")
        r_sq = float("NaN")
        slope = float("NaN")
        intercept = float("NaN")
    else:
        row = row.values
	ratios = np.log2(row[~np.isnan(row)])
        # ratios = row[~np.isnan(row)]
        times = timepoints[~np.isnan(row)]
        if len(ratios) == 1:
            # only present in input library
            score = float("NaN")
            r_sq = float("NaN")
            slope = float("NaN")
            intercept = float("NaN")
        elif len(ratios) == 2:
            # rise over run
            score = (ratios[1] - ratios[0]) / (times[1] - times[0])
            r_sq = float("NaN")
            slope = float("NaN")
            intercept = float("NaN")
        else:
            slope, intercept, r, _, _ = stats.linregress(times, ratios)
            score = slope
            r_sq = r ** 2

    return pd.Series({'score' : score, 'r_sq' : r_sq, 'slope' : slope, 'intercept' : intercept})


class Selection(DataContainer):
    """
    Class for a single selection experiment, consisting of multiple 
    timepoints. This class coordinates :py:class:`~seqlib.seqlib.SeqLib` 
    objects. Creating a :py:class:`~selection.Selection` requires a valid 
    *config* object, usually from a ``.json`` configuration file.
    """
    def __init__(self, config):
        DataContainer.__init__(self, config)
        self.libraries = dict()
        self.timepoints = list()
        self.normalize_wt = False
        self.ns_carryover_fn = None
        self.ns_carryover_kwargs = None
        self.use_barcode_variation = False

        try:
            if 'barcodes' in config:
                if 'map file' in config['barcodes']:
                    self.barcode_map = BarcodeMap(config['barcodes']
                                                        ['map file'])
                else:
                    self.barcode_map = None
            else:
                self.barcode_map = None

            libnames = list()
            bcmfiles = list()
            for lib in config['libraries']:
                if 'output directory' not in lib:
                    lib['output directory'] = self.output_base
                libtype = seqlib_type(lib)
                if libtype is None:
                    raise EnrichError("Unrecognized SeqLib config", self.name)
                elif libtype == "BarcodeVariantSeqLib":
                    new = BarcodeVariantSeqLib(lib, barcode_map=self.barcode_map)
                    bcmfiles.append(new.barcode_map.filename)
                else:
                    new = globals()[libtype](lib)

                if new.output_base is None:
                    new.set_output_base(self.output_base)

                if new.timepoint not in self.libraries:
                    self.libraries[new.timepoint] = list()
                self.libraries[new.timepoint].append(new)
                libnames.append(new.name)
            self.timepoints = sorted(self.libraries.keys())

            if len(set(libnames)) != len(libnames):
                raise EnrichError("Non-unique library names", self.name)

            if len(bcmfiles) == len(libnames): # all BarcodeVariant
                if len(set(bcmfiles)) == 1:    # all the same BarcodeMap
                    self.use_barcode_variation = True
                    unify_barcode_maps = False
                    if self.barcode_map is None: # same BarcodeMap specified for all SeqLibs
                        unify_barcode_maps = True
                    elif bcmfiles[0] != self.barcode_map.filename: # all SeqLibs are overriding the Selection BarcodeMap
                        unify_barcode_maps = True
                    else: # this BarcodeMap is being used for all SeqLibs
                        pass
                    if unify_barcode_maps:
                        self.barcode_map = self.libraries[0][0].barcode_map
                        for lib in self.library_list():
                            lib.barcode_map = self.barcode_map

            if self.use_barcode_variation is False and 'max barcode variation' in config['filters']:
                raise EnrichError("Unable to filter on 'max barcode variation'", self.name)

            self.set_filters(config['filters'], {'min count' : 0,
                                      'min input count' : 0,
                                      'min rsquared' : 0.0,
                                      'max barcode variation' : None})

            if 'carryover correction' in config:
                if config['carrover correction']['method'] == "nonsense":
                    self.ns_carryover_fn = nonsense_ns_carryover_apply_fn
                    self.ns_carryover_kwargs = {'position' : int(config['carryover correction']['position'])}
                # add additional methods here using "elif" blocks
                else:
                    raise EnrichError("Unrecognized nonspecific carryover correction", self.name)

            if 'normalize wt' in config:
                if config['normalize wt'] is True:
                    self.normalize_wt = True

        except KeyError as key:
            raise EnrichError("Missing required config value {key}".format(key=key), 
                              self.name)
        except ValueError as value:
            raise EnrichError("Invalid parameter value {value}".format(value=value), self.name)

        if len(self.timepoints) < 2:
            raise EnrichError("Insufficient number of timepoints", 
                              self.name)

        if 0 not in self.timepoints:
            raise EnrichError("Missing timepoint 0", self.name)
        if self.timepoints[0] != 0:
            raise EnrichError("Invalid negative timepoint", self.name)

        # identify what kind of counts data is present in all timepoints
        dtype_counts = list()
        for lib in self.library_list():
            dtype_counts.extend(lib.df_dict.keys())
        dtype_counts = Counter(dtype_counts)
        for dtype in dtype_counts:
            if dtype_counts[dtype] == len(config['libraries']):
                self.df_dict[dtype] = True
        if 'barcodes_unmapped' in self.df_dict.keys(): # special case for BarcodeVariantSeqLib
            del self.df_dict['barcodes_unmapped']
        if 'barcodes_low_abuncande' in self.df_dict.keys(): # special case for BarcodeVariantSeqLib
            del self.df_dict['barcodes_low_abuncande']
        if len(self.df_dict.keys()) == 0:
            raise EnrichError("No count data present across all timepoints", 
                              self.name)


    def library_list(self):
        """
        Return the :py:class:`~seqlib.seqlib.SeqLib` objects as a list, sorted by timepoint.
        """
        libs = list()
        for tp in self.timepoints:
            libs.extend(self.libraries[tp])
        return libs


    def count_timepoints(self):
        """
        Combine :py:class:`~seqlib.seqlib.SeqLib` objects into individual timepoints and 
        tabulate counts for each timepoint. Counts are stored in the local
        :py:class:`pandas.DataFrame`. To tabulate counts for individual mutations 
        (not variants), see :py:meth:`count_mutations`.
        """
        # calculate counts for each SeqLib
        logging.info("Counting for each timepoint [{name}]".format(name=self.name))
        for lib in self.library_list():
            lib.calculate()
            lib.dump_data() # dump the data to save memory while calculating
                            # important for large barcode datasets

        # restore relevant count data
        for lib in self.library_list():
            lib.restore_data(keys=self.df_dict.keys())

        # perform the calculations
        for dtype in self.df_dict:
            self.calc_counts(dtype)


    def calc_counts(self, dtype):
        """
        Tabulate counts for each timepoint and create the :py:class:`pandas.DataFrame` indicated by 
        *dtype*. All :py:class:`~seqlib.seqlib.SeqLib` objects need to be counted before calling 
        this method.
        """
        # combine all libraries for a given timepoint
        tp_counts = dict()
        for tp in self.timepoints:
            c = self.libraries[tp][0].df_dict[dtype]
            for lib in self.libraries[tp][1:]:
                c = c.add(lib.df_dict[dtype], fill_value=0)
            tp_counts[tp] = c

        tp_frame = tp_counts[0]
        cnames = ["{cname}.0".format(cname=x) for x in tp_counts[0].columns]
        for tp in self.timepoints[1:]:
            tp_frame = tp_frame.join(tp_counts[tp], how="outer", rsuffix=str(tp))
            cnames += ["{cname}.{tp}".format(cname=x, tp=tp) for x in tp_counts[tp].columns]
        tp_frame.columns = cnames
        # remove data that are not in the initial timepoint
        try:
            tp_frame.dropna(axis=0, how="any", subset=['count.0'], inplace=True)
        except TypeError:
            # old versions of pandas don't support inplace
            tp_frame = tp_frame.dropna(axis=0, how="any", subset=['count.0'])
        self.df_dict[dtype] = tp_frame


    def count_mutations(self):
        """
        Creates and populates :py:class:`pandas.DataFrame` objects for individual mutations. This method 
        should be called after all filtering has been completed. The new 
        :py:class:`pandas.DataFrame` objects have dtype 'mutations_nt' and 'mutations_aa' (only if the 
        data set is coding).
        """
        # needs to happen all the filtering/exclusion of variants
        for lib in self.libraries:
            lib.count_mutations()

        if all(x.is_coding() for x in self.library_list()):
            mutation_dtypes = ('mutations_nt', 'mutations_aa')
        else:
            mutation_dtypes = ('mutations_nt')
        for dtype in mutation_dtypes:
            self.calc_counts(dtype)
            self.calc_frequencies(dtype)
            self.calc_ratios(dtype)
            self.calc_enrichments(dtype)
        self.sort_data('score', keys=mutation_dtypes)


    def calculate(self):
        """
        Wrapper method to calculate counts, frequencies, ratios, and enrichment scores 
        for all data in the :py:class:`Selection`. Enrichment scores are only calculated 
        if there are more than two timepoints.
        """
        self.count_timepoints()
        if self.ns_carryover_fn is not None:
            self.nonspecific_carryover(self.ns_carryover_fn, **self.ns_carryover_kwargs)
        for dtype in self.df_dict:
            self.calc_frequencies(dtype)
            self.calc_ratios(dtype)
            self.calc_enrichments(dtype)
        if self.use_barcode_variation: # all SeqLibs are BarcodeVariantSeqLibs and use the same BarcodeMap
            self.calc_barcode_variation()
            self.add_variants_to_barcodes()
        self.sort_data('score')

        if 'variants' in self.df_dict.keys():
            self.filter_variant_data()
            self.calc_frequencies('variants')
            self.calc_ratios('variants')
            self.calc_enrichments('variants')

        if self.normalize_wt:
            self.normalize_variants_to_wt()



    def calc_frequencies(self, dtype):
        """
        Calculate frequencies for each element in the :py:class:`pandas.DataFrame` indicated by 
        *dtype* ('variant' or 'barcode').
        """
        for tp in self.timepoints:
            self.df_dict[dtype]['frequency.{tp}'.format(tp=tp)] =  \
                self.df_dict[dtype]['count.{tp}'.format(tp=tp)] / \
                float(self.df_dict[dtype]['count.{tp}'.format(tp=tp)].sum())


    def calc_ratios(self, dtype):
        """
        Calculate ratios for each element in the :py:class:`pandas.DataFrame` indicated by 
        *dtype* ('variant' or 'barcode'). Assumes frequencies have been 
        calculated by :py:meth:`calc_frequencies`.
        """
        for tp in self.timepoints:
            if tp == 0: # input library
                self.df_dict[dtype]['ratio.{tp}'.format(tp=tp)] = 1.0
            else:
                self.df_dict[dtype]['ratio.{tp}'.format(tp=tp)] =  \
                        self.df_dict[dtype]['frequency.{tp}'.format(tp=tp)] / \
                        self.df_dict[dtype]['frequency.0']


    def calc_enrichments(self, dtype):
        """
        Calculate enrichment scores and r-squared values for each element in the :py:class:`pandas.DataFrame` indicated by 
        *dtype* ('variant' or 'barcode'). Assumes ratios have been 
        calculated by :py:meth:`calc_ratios`. Calculations performed using 
        :py:func:`linear_enrichment_apply_fn`.
        """
        # apply the enrichment-calculating function to a DataFrame
        # containing only ratio data
        ratio_df = self.df_dict[dtype][['ratio.{tp}'.format(tp=x) for x in self.timepoints]]
        enrichments = ratio_df.apply(linear_enrichment_apply_fn, axis=1, args=[np.asarray(self.timepoints)])
        self.df_dict[dtype] = pd.concat([self.df_dict[dtype], enrichments], axis=1)


    def calc_barcode_variation(self):
        """
        Calculate the `coefficient of variation <http://en.wikipedia.org/wiki/Coefficient_of_variation>`_ 
        for each variant's barcode enrichment scores. Requires both variant and barcode 
        data for all timepoints.
        """
        self.df_dict['variants']['barcode.count'] = \
                self.df_dict['variants'].apply(barcode_count_apply_fn, 
                axis=1, args=[self.barcode_map]).astype("int32")
        barcode_cv = self.df_dict['variants'].apply(\
                barcode_variation_apply_fn, axis=1,
                args=[self.df_dict['barcodes'], self.barcode_map])
        self.df_dict['variants']['scored.unique.barcodes'] = \
                barcode_cv['scored.unique.barcodes'].astype("int32")
        self.df_dict['variants']['barcode.cv'] = barcode_cv['barcode.cv']


    def add_variants_to_barcodes(self):
        """
        Add the associated variant information to each row of the barcode :py:class:`pandas.DataFrame`.
        """
        self.df_dict['barcodes']['variant'] = self.df_dict['barcodes'].apply(lambda x: self.barcode_map.bc_variant_strings[x.name], axis=1)


    def nonspecific_carryover(self, ns_apply_fn, **kwargs):
        """
        Correct the counts in the 'variants' :py:class:`pandas.DataFrame` for nonspecific carryover. 
        Nonspecific counts are defined by *ns_apply_fn* and its *kwargs*, which 
        takes a row as an argument and returns ``True`` if the row's counts 
        are nonspecific.
        """
        logging.info("Applying nonspecific carryover correction [{name}]".format(name=self.name))
        dtype = 'variants'
        ns_data = self.df_dict[dtype][self.df_dict[dtype].apply(ns_apply_fn, axis=1, **kwargs)]
        read_totals = [self.df_dict[dtype]['count.{tp}'.format(tp=tp)].sum() \
                       for tp in self.timepoints]
        read_totals = dict(zip(self.timepoints, read_totals))
        ns_frequencies = [ns_data['count.{tp}'.format(tp=tp)].sum() / \
                             read_totals[tp] for tp in self.timepoints]
        ns_frequencies = dict(zip(self.timepoints, ns_frequencies))
        for tp in self.timepoints[1:]: # don't modify time 0
            ns_mod = (ns_frequencies[tp] / ns_frequencies[0])
            self.df_dict[dtype]['count.{tp}'.format(tp=tp)] = self.df_dict[dtype]['count.{tp}'.format(tp=tp)] - self.df_dict[dtype]['count.{tp}'.format(tp=tp)] * ns_mod
            self.df_dict[dtype]['count.{tp}'.format(tp=tp)] = self.df_dict[dtype]['count.{tp}'.format(tp=tp)].astype("int32")


    def filter_variant_data(self):
        """
        Apply the filtering functions to the variant data, based on the filter 
        options present in the configuration object. Filtering is performed 
        using the appropriate apply function. Frequencies, ratios, and 
        enrichments must be recalculated after filtering.

        The data are written to the subdirectory ``"pre-filter"`` before filtering.
        """
        self.write_data(subdirectory="pre-filter")
        # for each filter that's specified
        # apply the filter
        if self.filters['max barcode variation']:
            nrows = len(self.df_dict['variants'])
            self.df_dict['variants'] = \
                    self.df_dict['variants'][self.df_dict['variants'].apply(\
                        barcode_varation_filter, axis=1, 
                        args=[self.filters['max barcode variation']])]
            self.filter_stats['max barcode variation'] = \
                    nrows - len(self.df_dict['variants'])
        if self.filters['min count'] > 0:
            nrows = len(self.df_dict['variants'])
            self.df_dict['variants'] = \
                    self.df_dict['variants'][self.df_dict['variants'].apply(\
                        min_count_filter, axis=1, 
                        args=[self.filters['min count']])]
            self.filter_stats['min count'] = \
                    nrows - len(self.df_dict['variants'])
        if self.filters['min input count'] > 0:
            nrows = len(self.df_dict['variants'])
            self.df_dict['variants'] = \
                    self.df_dict['variants'][self.df_dict['variants'].apply(\
                        min_input_count_filter, axis=1, 
                        args=[self.filters['min count']])]
            self.filter_stats['min input count'] = \
                    nrows - len(self.df_dict['variants'])
        if self.filters['min rsquared'] > 0.0:
            nrows = len(self.df_dict['variants'])
            self.df_dict['variants'] = \
                    self.df_dict['variants'][self.df_dict['variants'].apply(\
                        min_rsq_filter, axis=1, 
                        args=[self.filters['min rsquared']])]
            self.filter_stats['min rsquared'] = \
                    nrows - len(self.df_dict['variants'])

        self.filter_stats['total'] = sum(self.filter_stats.values())


    def write_all(self):
        self.write_data()
        for lib in self.library_list():
            lib.write_all()


    def normalize_variants_to_wt(self):
        """
        Normalizes variant scores (or ratios if only two timepoints are present) such that the wild type value is "neutral" (0 for scores, 1 for ratios). Does nothing for barcode-only data.

        The data are written to the subdirectory ``"pre-wtnorm"`` before normalization.
        """
        if 'variant' in self.df_dict.keys():
            self.write_data(subdirectory="pre-wtnorm")
            self.df_dict['variant']['score'] = self.df_dict['variant']['score'] - self.df_dict['variant']['score'][WILD_TYPE_VARIANT]
