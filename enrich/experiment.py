from __future__ import print_function
from enrich_error import EnrichError
from datacontainer import DataContainer
import selection
import time
import sys
import pandas as pd


def condition_cv_apply_fn(row, condition):
    """
    :py:meth:`pandas.DataFrame.apply` function for calculating the 
    coefficient of variation for a variant's score (or ratio) in the 
    condition.
    """
    bc_scores = barcode_data.ix[mapping.variants[row.name]]['score']
    bc_scores = bc_scores[np.invert(np.isnan(bc_scores))]
    cv = stats.variation(bc_scores)
    return pd.Series({'scored.unique.barcodes' : len(bc_scores), \
                      'barcode.cv' : cv})


class Experiment(DataContainer):
    """
    Class for a coordinating multiple :py:class:`~.selection.Selection` 
    objects. Creating an 
    :py:class:`~experiment.Experiment` requires a valid *config* object, 
    usually from a ``.json`` configuration file.
    """
    def __init__(self, config):
        DataContainer.__init__(self, config)
        self.conditions = dict()
        self.control = None
        self.normalize_wt = False

        try:
            if 'normalize wt' in config:
                if config['normalize wt'] is True:
                    self.normalize_wt = True
            for cnd in config['conditions']:
                if not cnd['label'].isalnum():
                    raise EnrichError("Alphanumeric label required for condition '{label}'".format(label=cnd['label']), self.name)
                for sel_config in cnd['selections']: # assign output base if not present
                    if 'output directory' not in sel_config:
                        sel_config['output directory'] = self.output_base
                if cnd['label'] not in self.conditions:
                    self.conditions[cnd['label']] = [selection.Selection(x) for x in cnd['selections']]
                else:
                    raise EnrichError("Non-unique condition label '{label}'".format(label=cnd['label']), self.name)
                if 'control' in cnd:
                    if cnd['control']:
                        if self.control is None:
                            self.control = self.conditions[cnd['label']]
                        else:
                            raise EnrichError("Multiple control conditions", self.name)
        except KeyError as key:
            raise EnrichError("Missing required config value {key}".format(key=key), 
                              self.name)

        all_selections = self.selection_list()
        for dtype in all_selections[0].df_dict:
            if all(dtype in x.df_dict for x in all_selections):
                self.df_dict[dtype] = True
        if len(self.df_dict.keys()) == 0:
            raise EnrichError("No enrichment data present across all selections", 
                              self.name)

        # ensure consistency for wild type normalization
        for sel in all_selections:
            sel.normalize_wt = self.normalize_wt


    def selection_list(self):
        """
        Return the :py:class:`~selection.Selection` objects as a list.
        """
        selections = list()
        for key in self.conditions:
            selections.extend(self.conditions[key])
        return selections


    def calculate(self):
        """
        Calculate scores for all :py:class:`~selection.Selection` objects.
        """
        cnames = dict()
        for dtype in self.df_dict:
            self.df_dict[dtype] = pd.DataFrame()
            cnames[dtype] = list()

        for c in self.conditions:
            s_id = 1
            for s in self.conditions[c]:
                s_label = "{condition}.{id}".format(condition=c, id=s_id)
                s_id += 1
                s.calculate()
                for dtype in self.df_dict:
                    # score and r_sq columns
                    self.df_dict[dtype] = self.df_dict[dtype].join(s.df_dict[dtype][['score', 'r_sq']],
                        how="outer", rsuffix=s_label)
                    cnames[dtype].extend(["{cname}.{sel}".format(cname=x, sel=s_label) for x in ['score', 'r_sq']])
                s.dump_data()
        for dtype in self.df_dict:
            self.df_dict[dtype].columns = cnames[dtype]

        for s in self.selection_list():
            s.restore_data()


    def calc_variation(self):
        """
        Calculate the coefficient of variation for each variant's scores or ratios in each condition.
        """
        for dtype in self.df_dict:
            for c in self.conditions:
                c_columns = [x.startswith("score.{cnd}".format(cnd=c)) for x in self.df_dict[dtype].columns]
                c_values = self.df_dict[dtype][self.df_dict[dtype].columns[c_columns]]                
                self.df_dict[dtype]['{cnd}.cv'.format(cnd=c)] = c_values.apply(stats.variation, axis=1)


    def filter_data(self):
        """
        Apply the filtering functions to the data, based on the filter 
        options present in the configuration object. Filtering is performed 
        using the appropriate apply function.
        """
        self.write_data(os.path.join(self.output_base, "experiment_prefilter"))
        # for each filter that's specified
        # apply the filter
        self.filter_stats['total'] = sum(self.filter_stats.values())


    def write_all(self):
        self.write_data()
        for c in self.conditions:
            for sel in self.conditions[c]:
                sel.write_all()



