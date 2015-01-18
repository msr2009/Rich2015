from __future__ import print_function
import time
import logging
from enrich_error import EnrichError
from datacontainer import DataContainer
import os.path
import enrich_plot


class SeqLib(DataContainer):
    """
    Abstract class for handling count data from a single sequencing library.
    Creating a :py:class:`seqlib.seqlib.SeqLib` requires a valid *config* 
    object, usually from a ``.json`` configuration file.

    .. note:: Example configuration files can be found in the documentation \
    for derived classes.
    """
    def __init__(self, config):
        DataContainer.__init__(self, config)

        try:
            self.timepoint = int(config['timepoint'])
        except KeyError as key:
            raise EnrichError("Missing required config value '{key}'".format(key=key), 
                              self.name)
        except ValueError as value:
            raise EnrichError("Invalid parameter value {value}".format(value=value), self.name)

        if 'report filtered reads' in config:
            self.report_filtered = config['report filtered reads']
        else:
            self.report_filtered = False


    def calculate(self):
        """
        Pure virtual method that defines how the data are counted.
        """
        raise NotImplementedError("must be implemented by subclass")


    def report_filtered_read(self, fq, filter_flags):
        """
        Write the :py:class:`~fqread.FQRead` object *fq* to the ``DEBUG``
        logging . The dictionary *filter_flags* contains ``True`` 
        values for each filtering option that applies to *fq*. Keys in 
        *filter_flags* are converted to messages using the 
        ``DataContainer._filter_messages`` dictionary.
        """
        logging.debug("Filtered read ({messages}) [{name}]\n{read!s}".format(
                      messages=', '.join(DataContainer._filter_messages[x] 
                                for x in filter_flags if filter_flags[x]), 
                      name=self.name, read=fq))


    def write_all(self):
        self.write_data()


    def make_plots(self, subdirectory=None, keys=None):
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
            fname = os.path.join(directory, fix_filename(key + ".pdf"))
            plot = enrich_plot.histogram(self.df_dict[key], "count")
            plot.savefig(fname)


