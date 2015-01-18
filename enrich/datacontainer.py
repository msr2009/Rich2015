from __future__ import print_function
from enrich_error import EnrichError
import sys
import time
import os
import logging
import resource
import numpy as np
import pandas as pd


def fix_filename(s):
    """
    Clean up a filename *s* by removing invalid characters and converting 
    spaces to underscores. Returns the cleaned filename.
    """
    fname = "".join(c for c in s if c.isalnum() or c in (' ._~'))
    fname = fname.replace(' ', '_')
    return fname


class DataContainer(object):
    """
    Abstract class for all data-containing classes 
    (:py:class:`seqlib.seqlib.SeqLib`, :py:class:`selection.Selection`, and 
    :py:class:`experiment.Experiment`). Creating a 
    :py:class:`datacontainer.DataContainer`  requires a valid *config* object, 
    usually from a ``.json`` configuration file.

    .. note:: Example configuration files can be found in the documentation \
    for derived classes.

    Log file messages for filtered reads are defined here in the 
    ``_filter_messages`` dictionary. New read filtering options must have an 
    associated log file output message added to the dictionary.

    .. literalinclude:: ../datacontainer.py
        :lines: 41-56
    """

    # Note: the following block is referenced by line number above
    # When adding new messages, update the documentation line numbers also!
    _filter_messages = {
            # SeqLib messages
            'remove unresolvable' : "unresolvable mismatch",
            'min quality' : "single-base quality",
            'avg quality' : "average quality",
            'max mutations' : "excess mutations",
            'chastity' : "not chaste",
            'remove overlap indels' : "indel in read overlap",
            'merge failure' : "unable to merge reads",
            # Selection messages
            'min count' : "not enough variant reads",
            'min input count' : "not enough variant reads in input",
            'min rsquared' : "low r-squared",
            'max barcode variation' : "high barcode CV"
            # Experiment messages
        }


    def __init__(self, config):
        self.name = "Unnamed" + self.__class__.__name__
        self.df_dict = dict()
        self.df_files = dict()
        self.filters = None
        self.filter_stats = None
        self.output_base = None
        
        try:
            self.name = config['name']
        except KeyError as key:
            raise EnrichError("Missing required config value {key}".format(key=key), 
                              self.name)
        
        if 'output directory' in config:
            self.set_output_base(config['output directory'])


    def set_output_base(self, dirname):
        """
        Sets the object's base output directory (used for 
        :py:meth:`dump_data` and other class-specific methods) to *dirname* 
        and creates the directory if it doesn't exist.
        """
        try:
            if not os.path.exists(dirname):
                os.makedirs(dirname)
        except OSError:
            raise EnrichError("Failed to create output directory", self.name)
        self.output_base = dirname


    def dump_data(self, keys=None):
        """
        Save the :py:class:`pandas.DataFrame` objects as tab-separated files and 
        set the data to ``None`` to save memory. The 
        file names are stored for use by :py:meth:`restore_data`.
 
        The optional *keys* parameter is a list of types of data to be 
        dumped (variant, barcode, etc.). By default, all data are dumped.
        """
        self.log_memory_usage()
        if keys is None:
            keys = self.df_dict.keys()
        keys = [k for k in keys if self.df_dict[k] is not None]
        logging.info("Initiating data frame dump ({keys}) [{name}]".format(name=self.name, keys=", ".join(keys)))
        self.df_files.update(self.write_data(subdirectory="dump", keys=keys))
        for key in keys:
            self.df_dict[key] = None
        self.log_memory_usage()


    def write_data(self, subdirectory=None, keys=None):
        """
        Save the :py:class:`pandas.DataFrame` objects as tab-separated files 
        with the same name as the object.

        The optional *keys* parameter is a list of types of data to be 
        saved (variant, barcode, etc.). By default, all data are saved.

        Returns a dictionary with *keys* as the keys and corresponding filenames for the ``.tsv`` files as the values. This dictionary is required by :py:meth:`dump_data`
        """
        fname_dict = dict()
        if subdirectory is not None:
            directory = os.path.join(self.output_base, fix_filename(subdirectory), fix_filename(self.name))
        else:
            directory = os.path.join(self.output_base, fix_filename(self.name))
        
        try:
            if not os.path.exists(directory):
                os.makedirs(directory)
        except OSError:
            raise EnrichError("Failed to create output directory", self.name)

        if keys is None:
            keys = self.df_dict.keys()
        keys = [k for k in keys if self.df_dict[k] is not None]

        for key in keys:
            fname = os.path.join(directory, fix_filename(key + ".tsv"))
            self.df_dict[key].to_csv(fname, 
                    sep="\t", na_rep="NaN", float_format="%.4g", 
                    index_label="sequence")
            fname_dict[key] = fname
        logging.info("Successfully wrote data frames ({keys}) [{name}]".format(name=self.name, keys=", ".join(keys)))
        return fname_dict


    def restore_data(self, keys=None):
        """
        Load the data from the ``.tsv`` files written by :py:meth:`dump_data`.

        The optional *keys* parameter is a list of types of data to be 
        restored (variant, barcode, etc.). By default, all data are restored.
        """
        self.log_memory_usage()
        logging.info("Restoring data frames from dump ({keys}) [{name}]".format(name=self.name, keys=", ".join(keys)))
        if keys is None:
            keys = self.df_files.keys()
        for key in keys:
            self.df_dict[key] = pd.DataFrame.from_csv(self.df_files[key], sep="\t")
            self.df_files[key] = None
        self.log_memory_usage()


    def set_filters(self, config_filters, default_filters):
        """
        Sets the filtering options using the values from the *config_filters* 
        dictionary and *default_filters* dictionary. 

        .. note:: To help prevent user error, *config_filters* must be a \
        subset of *default_filters*.
        """
        self.filters = default_filters

        for key in self.filters:
            if key in config_filters:
                self.filters[key] = int(config_filters[key])

        unused = list()
        for key in config_filters:
            if key not in self.filters:
                unused.append(key)
        if len(unused) > 0:
            logging.warning("Unused filter parameters ({unused}) [{name}]".format(unused=', '.join(unused), name=self.name))

        self.filter_stats = dict()
        for key in self.filters:
            self.filter_stats[key] = 0
        self.filter_stats['total'] = 0


    def report_filter_stats(self):
        try:
            output_dir = os.path.join(self.output_base, fix_filename(self.name))
        except AttributeError:
            raise EnrichError("Invalid output directory specified for object", self.name)
        try:
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)
        except OSError:
            raise EnrichError("Failed to create output directory", self.name)
        with open(os.path.join(output_dir, "filter_stats.txt"), "w") as handle:
            elements = list()
            for key in sorted(self.filter_stats, key=self.filter_stats.__getitem__, reverse=True):
                if key != 'total':
                    print(DataContainer._filter_messages[key], self.filter_stats[key], sep="\t", file=handle)
            print('total', self.filter_stats['total'], sep="\t", file=handle)
            logging.info("Wrote filtering statistics [{name}]".format(name=self.name))


    def calculate(self):
        """
        Pure virtual method that defines how the data are calculated. 
        """
        raise NotImplementedError("must be implemented by subclass")


    def make_plots(self, subdirectory=None, keys=None):
        """
        Pure virtual method for creating plots.
        """
        raise NotImplementedError("must be implemented by subclass")


    def write_all(self):
        """
        Pure virtual method for writing :py:class:`pandas.DataFrame` contents to disk.
        """
        raise NotImplementedError("must be implemented by subclass")


    def sort_data(self, column, keys=None):
        """
        Sort the :py:class:`pandas.DataFrame` objects according the to the values in column *key*. The data are 
        sorted in descending order, with the ``NaN`` values at the end.

        The optional *keys* parameter is a list of types of data to be 
        sorted (variant, barcode, etc.). By default, all data are sorted.
        """
        if keys is None:
            keys = self.df_dict.keys()
        for key in keys:
            nan = self.df_dict[key][np.isnan(self.df_dict[key][column])]
            if len(nan) > 0:
                no_nan = self.df_dict[key][np.invert(np.isnan(self.df_dict[key][column]))]
                no_nan.sort(column, ascending=False, inplace=True)
                self.df_dict[key] = no_nan.append(nan)
            else:
                self.df_dict[key].sort(column, ascending=False, inplace=True)


    def log_memory_usage(self):
        """
        Write the current memory usage (for the whole process) to the log file.
        """
        logging.info("Current process memory usage is {size} [{name}]".format(size=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, name=self.name))

