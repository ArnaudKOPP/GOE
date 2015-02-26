# coding=utf-8
"""
Go enrichement
"""

__author__ = "Arnaud KOPP"
__copyright__ = "© 2014-2015 KOPP Arnaud All Rights Reserved"
__credits__ = ["KOPP Arnaud"]
__license__ = "GNU GPL V2.0"
__version__ = "1.0"
__maintainer__ = "Arnaud KOPP"
__email__ = "kopp.arnaud@gmail.com"
__status__ = "Dev"

import collections
import numpy as np
import scipy.stats
import pandas as pd
import urllib.request
from HTSDataMining.utils import reporthook
from HTSDataMining.stat import adjustpvalues
import logging
log = logging.getLogger(__name__)

typedef_tag, term_tag = "[Typedef]", "[Term]"


class OBOreader(object):
    """
    Parse obo file
    """
    def __init__(self, obo_file="go.obo"):
        try:
            self._handle = open(obo_file)
        except IOError:
            log.info("Download ontologies")
            urllib.request.urlretrieve(
                "http://www.berkeleybop.org/ontologies/go/go.obo", "go.obo", reporthook)
            print('')  # add this for begin @ new line (it's ugly)
            self._handle = open(obo_file)

    def __iter__(self):
        line = self._handle.readline()
        if not line.startswith(term_tag):
            self.read_until(self._handle, term_tag)
        while 1:
            yield self.next()

    def next(self):
        """

        :return: :raise StopIteration:
        """
        lines = []
        line = self._handle.readline()
        if not line or line.startswith(typedef_tag):
            raise StopIteration

        # read until the next tag and save everything in between
        while 1:
            pos = self._handle.tell()  # save current postion for roll-back
            line = self._handle.readline()
            if not line or (line.startswith(typedef_tag)
                            or line.startswith(term_tag)):
                self._handle.seek(pos)  # roll-back
                break
            lines.append(line)

        rec = GOTerm()
        for line in lines:
            if line.startswith("id:"):
                rec.id = self.after_colon(line)
            if line.startswith("alt_id:"):
                rec.alt_ids.append(self.after_colon(line))
            elif line.startswith("name:"):
                rec.name = self.after_colon(line)
            elif line.startswith("namespace:"):
                rec.namespace = self.after_colon(line)
            elif line.startswith("is_a:"):
                rec.parents.append(self.after_colon(line).split()[0])
            elif line.startswith("is_obsolete:") and self.after_colon(line) == "true":
                rec.is_obsolete = True

        return rec

    @staticmethod
    def after_colon(line):
        """
        macro for getting anything after the :
        :param line:
        :return:
        """
        return line.split(":", 1)[1].strip()

    @staticmethod
    def read_until(handle, start):
        """
        read each line until it has a certain start, and then puts the start tag back
        :param handle:
        :param start:
        :return:
        """
        while 1:
            pos = handle.tell()
            line = handle.readline()
            if not line:
                break
            if line.startswith(start):
                handle.seek(pos)
                return
        raise EOFError("%s tag cannot be found" % start)


class GOTerm(object):
    """
    Go term, contain a lot more attribut than interfaced here
    """

    def __init__(self):
        self.id = ""  # GO:xxxxxx
        self.name = ""  # description
        self.namespace = ""  # BP, CC, MF
        self._parents = []  # is_a basestring of parents
        self.parents = []  # parent records
        self.children = []  # children records
        self.level = -1  # distance from root node
        self.is_obsolete = False  # is_obsolete
        self.alt_ids = []  # alternative identifiers

    def __str__(self):
        obsolete = "obsolete" if self.is_obsolete else ""
        return "%s\tlevel-%02d\t%s [%s] %s" % (self.id, self.level, self.name, self.namespace, obsolete)

    def __repr__(self):
        return "GOTerm('%s')" % self.id

    def has_parent(self, term):
        """
        Check if term has parents
        :param term:
        :return:
        """
        for p in self.parents:
            if p.id == term or p.has_parent(term):
                return True
        return False

    def has_child(self, term):
        """
        Check if term has children
        :param term:
        :return:
        """
        for p in self.children:
            if p.id == term or p.has_child(term):
                return True
        return False

    def get_all_parents(self):
        """
        Get all parents of term
        :return:
        """
        all_parents = set()
        for p in self.parents:
            all_parents.add(p.id)
            all_parents |= p.get_all_parents()
        return all_parents

    def get_all_children(self):
        """
        Get all children of term
        :return:
        """
        all_children = set()
        for p in self.children:
            all_children.add(p.id)
            all_children |= p.get_all_children()
        return all_children

    def get_all_parent_edges(self):
        """
        Get all edge parents
        :return:
        """
        all_parent_edges = set()
        for p in self.parents:
            all_parent_edges.add((self.id, p.id))
            all_parent_edges |= p.get_all_parent_edges()
        return all_parent_edges

    def get_all_child_edges(self):
        """
        Get all edge child
        :return:
        """
        all_child_edges = set()
        for p in self.children:
            all_child_edges.add((p.id, self.id))
            all_child_edges |= p.get_all_child_edges()
        return all_child_edges


class GOtree(object):
    """
    Class for construct a GO tree
    :param obo_file:
    """

    def __init__(self, obo_file="go.obo"):
        self.go_Term = {}
        self.load_obo_file(obo_file)

    def load_obo_file(self, obo_file):
        """
        load obo file into obo reader
        :param obo_file:
        """
        log.info("Load obo file ", obo_file)
        obo_reader = OBOreader(obo_file)
        for rec in obo_reader:
            self.go_Term[rec.id] = rec
            for alt in rec.alt_ids:
                self.go_Term[alt] = rec

        self.populate_terms()
        log.info("All GO nodes imported : ", len(self.go_Term))

    def populate_terms(self):
        """
        Construct go tree
        :return:
        """

        def depth(rec):
            """

            :param rec:
            :return:
            """
            if rec.level < 0:
                if not rec.parents:
                    rec.level = 0
                else:
                    rec.level = min(depth(rec) for rec in rec.parents) + 1
            return rec.level

        # make the parents references to the GO terms
        for rec in self.go_Term.items():
            rec = rec[1]
            rec.parents = [self.go_Term[x] for x in rec._parents]

        # populate children and levels
        for rec in self.go_Term.items():
            rec = rec[1]
            for p in rec.parents:
                p.children.append(rec)

            if rec.level < 0:
                depth(rec)

    def print_all_go_id(self):
        """
        Print all go ID
        """
        for rec_id, rec in sorted(self.go_Term.items()):
            print(rec)

    def query_term(self, term, verbose=True):
        """
        Search term
        :param term:
        :param verbose:
        :return:
        """
        if term not in self.go_Term:
            print("Term %s not found! ", term)
            return

        rec = self.go_Term[term]
        if verbose:
            print("all parents: ", rec.get_all_parents())
            print("all children: ", rec.get_all_children())

        return rec

    def paths_to_top(self, term):
        """
        search path to the top of tree
        :param term:
        :return:
        """
        if term not in self.go_Term:
            print("Term %s not found!", term)
            return

        def _paths_to_top_recursive(rec):
            if rec.level == 0:
                return [[rec]]
            paths = []
            for parent in rec.parents:
                top_paths = _paths_to_top_recursive(parent)
                for top_path in top_paths:
                    top_path.append(rec)
                    paths.append(top_path)
            return paths

        go_term = self.go_Term[term]
        return _paths_to_top_recursive(go_term)

    def _label_wrap(self, label):
        wrapped_label = r"%s\n%s" % (label, self.go_Term[label].name.replace(",", r"\n"))
        return wrapped_label

    def update_association(self, association):
        """
        Add parents in association dict
        :param association:
        """
        bad_terms = set()
        log.info("Update association")
        for key, terms in association.association.items():
            parents = set()
            for term in terms:
                try:
                    parents.update(self.go_Term[term].get_all_parents())
                except:
                    bad_terms.add(term)
            terms.update(parents)
        if bad_terms:
            print("terms not found: ", bad_terms)


class Association(object):
    """
    Association file in csv format
    first col = id
    second col = GO id reference
    id	go_id
    3804	GO:0003823
    3804	GO:0004872
    3804	GO:0005515
    3804	GO:0005886

    """

    def __init__(self, file):
        self.association = {}
        self._load_association_file(file)

    def _load_association_file(self, file):
        log.info("Load association file")
        try:
            assoc = pd.read_csv(file)
            assoc = assoc.dropna(axis=0)
            self._make_association(assoc)
        except:
            try:
                assoc = pd.read_csv(file, sep=',')
                assoc = assoc.dropna(axis=0)
                self._make_association(assoc)
            except Exception as e:
                print(e)

    def _make_association(self, assoc_data_frame):
        log.info("Making association ...")
        datagp = assoc_data_frame.groupby('id')
        for gene in assoc_data_frame.id.unique():
            go = datagp.get_group(gene)['go_id']
            self.association[gene] = set(go)

    def query(self, term):
        """
        search a term in association
        :param term:
        :return:
        """
        return self.association[term]


class GOEnrichmentRecord(object):
    """
    Represents one result (from a single GOTerm) in the GOEnrichmentStudy
    """
    _fields = "id ratio_in_study ratio_in_pop p_uncorrected description ".split()

    def __init__(self, id, ratio_in_study, ratio_in_pop, p_uncorrected):
        self.id = id
        self.ratio_in_study = ratio_in_study
        self.ratio_in_pop = ratio_in_pop
        self.p_uncorrected = p_uncorrected
        self.description = None
        self.goterm = None

    def __setattr__(self, name, value):
        self.__dict__[name] = value

    def __str__(self, indent=False):
        field_data = [self.__dict__[f] for f in self._fields]
        field_formatter = ["%s"] * 2 + ["%d/%d"] * 2 + ["%.3g/%.3g"] * 1 + ["%s"] * 1
        assert len(field_data) == len(field_formatter)
        return "\t".join(a % b for (a, b) in zip(field_formatter, field_data))

    def __repr__(self):
        return "GOEnrichmentRecord(%s)" % self.id

    def find_goterm(self, go):
        """
        Find go terminology
        :param go:
        """
        if self.id in go.go_Term:
            self.goterm = go.go_Term[self.id]
            self.description = self.goterm.name


class EnrichmentStudy(object):
    """
    Runs Fisher's exact test, as well as multiple corrections
    study file contain id
    pop file contain background id
    assoc file is csv format
    """

    def __init__(self, study, pop, assoc, compare=False):
        self.alpha = 0.05
        self.pval = 0.05
        self.compare = compare
        self.ration = 1
        self.indent = True
        self.min_ratio = self.ration
        if self.min_ratio is not None:
            assert 1 <= self.min_ratio <= 2
        assert 0 < self.alpha < 1, "Test-wise alpha must fall between (0, 1)"

        self.results = []
        self.study, self.pop = self.read_geneset(study, pop, compare=self.compare)
        self.association = Association(assoc)
        self.go_tree = GOtree()
        self.go_tree.update_association(self.association)

        self.term_study = self.count_terms(self.study, self.association, self.go_tree)
        self.term_pop = self.count_terms(self.pop, self.association, self.go_tree)

        self.pop_n, self.study_n = len(self.pop), len(self.study)

        self.run()

    def run(self):
        """
        run all
        :return:
        """
        for term, study_count in self.term_study.items():
            pop_count = self.term_pop[term]
            p = scipy.stats.fisher_exact(([[study_count, self.study_n], [pop_count, self.pop_n]]))

            one_record = GOEnrichmentRecord(id=term, p_uncorrected=p, ratio_in_study=(study_count, self.study_n),
                                            ratio_in_pop=(pop_count, self.pop_n))

            self.results.append(one_record)

        self.results.sort(key=lambda r: r.p_uncorrected[1])
        self.results = self.results

        for rec in self.results:
            # get go term for description and level
            rec.find_goterm(self.go_tree)

        return self.results

    def to_dataframe(self):
        """
        construct a numpy ndarray for storing result
        :return:
        """
        size = len(self.results)
        dataframe = np.zeros(size, dtype=[('ID', object), ('ratio_in_study', object), ('ratio_in_pop', object),
                                          ('Fischer Test OddsRation', object), ('Fischer Test P-value', object),
                                          ('FDR (BH)', object), ('BY', object), ('holm', object), ('hochberg', object),
                                          ('Description', object)])
        i = 0
        # # pretty ugly -> need to be more pythonic !!
        for record in self.results:
            assert isinstance(record, GOEnrichmentRecord)
            dataframe["ID"][i] = record.id
            dataframe["ratio_in_study"][i] = "%d/%d" % record.ratio_in_study
            dataframe["ratio_in_pop"][i] = "%d/%d" % record.ratio_in_pop
            dataframe["Fischer Test OddsRation"][i] = "%.3g" % record.p_uncorrected[0]
            dataframe["Fischer Test P-value"][i] = "%.3g" % record.p_uncorrected[1]
            dataframe["Description"][i] = record.description
            i += 1
        dataframe["holm"] = adjustpvalues(pvalues=dataframe["Fischer Test P-value"], method='holm')
        dataframe["hochberg"] = adjustpvalues(pvalues=dataframe["Fischer Test P-value"], method='hochberg')
        dataframe["FDR (BH)"] = adjustpvalues(pvalues=dataframe["Fischer Test P-value"])
        dataframe["BY"] = adjustpvalues(pvalues=dataframe["Fischer Test P-value"], method='BY')

        return dataframe

    @staticmethod
    def read_geneset(study_fn, pop_fn, compare=False):
        """
        Read study and population file
        :param study_fn:
        :param pop_fn:
        :param compare:
        :return:
        """
        log.info('Load study and population')
        pop = set(_.strip() for _ in open(pop_fn) if _.strip())
        study = frozenset(_.strip() for _ in open(study_fn) if _.strip())
        # some times the pop is a second group to compare, rather than the
        # population in that case, we need to make sure the overlapping terms
        # are removed first
        if compare:
            common = pop & study
            pop |= study
            pop -= common
            study -= common
            log.debug("removed ", len(common), " overlapping items")
            log.debug("Set 1: {0}, Set 2: {1}".format(len(study), len(pop)))
        return study, pop

    @staticmethod
    def count_terms(geneset, assoc, go_tree):
        """
        count the number of terms in the study group
        :param geneset:
        :param assoc:
        :param go_tree:
        """
        term_cnt = collections.defaultdict(int)
        for gene in geneset:
            try:
                for x in assoc.association[int(gene)]:
                    if x in go_tree.go_Term:
                        term_cnt[go_tree.go_Term[x].id] += 1
            except:
                continue
        return term_cnt