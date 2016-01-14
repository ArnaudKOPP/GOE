# coding=utf-8
"""
Go enrichement
"""

__author__ = "Arnaud KOPP"
__copyright__ = "© 2014-2016 KOPP Arnaud All Rights Reserved"
__credits__ = ["KOPP Arnaud"]
__license__ = "GPLv3"
__maintainer__ = "Arnaud KOPP"
__email__ = "kopp.arnaud@gmail.com"
__status__ = "Production"

import collections
import sys
import os
import numpy as np
import scipy.stats
import pandas as pd
import urllib.request
from GOE.utils import reporthook
from GOE.stat import adjustpvalues
import re
import logging

log = logging.getLogger(__name__)

class OBOreader(object):
    """
    Parse obo file
    """

    def __init__(self, obo_file="go.obo"):
        if os.path.isfile(obo_file):
            self._handle = obo_file
        else:
            log.info("Download ontologies")
            urllib.request.urlretrieve(
                "http://geneontology.org/ontology/go.obo", "go.obo", reporthook)
            print('')  # add this for begin @ new line (it's ugly)
            self._handle = obo_file

    def __iter__(self):
        """Return one GO Term record at a time from an obo file."""
        # Wait to open file until needed. Automatically close file when done.
        with open(self._handle) as fstream:
            rec_curr = None # Stores current GO Term
            for lnum, line in enumerate(fstream):
                # obo lines start with any of: [Term], [Typedef], /^\S+:/, or /^\s*/
                if line[0:6] == "[Term]":
                    rec_curr = self._init_goterm_ref(rec_curr, "Term", lnum)
                elif line[0:9] == "[Typedef]":
                    pass # Original OBOReader did not store these
                elif rec_curr is not None:
                    line = line.rstrip() # chomp
                    if ":" in line:
                        self._add_to_ref(rec_curr, line, lnum)
                    elif line == "":
                        if rec_curr is not None:
                            yield rec_curr
                            rec_curr = None
                    else:
                        self._die("UNEXPECTED LINE CONTENT: {L}".format(L=line), lnum)
            # Return last record, if necessary
            if rec_curr is not None:
                yield rec_curr

    def _init_goterm_ref(self, rec_curr, name, lnum):
        """Initialize new reference and perform checks."""
        if rec_curr is None:
            return GOTerm()
        msg = "PREVIOUS {REC} WAS NOT TERMINATED AS EXPECTED".format(REC=name)
        self._die(msg, lnum)

    def _add_to_ref(self, rec_curr, line, lnum):
        """Add new fields to the current reference."""
        # Examples of record lines containing ':' include:
        #   id: GO:0000002
        #   name: mitochondrial genome maintenance
        #   namespace: biological_process
        #   def: "The maintenance of ...
        #   is_a: GO:0007005 ! mitochondrion organization
        mtch = re.match(r'^(\S+):\s*(\S.*)$', line)
        if mtch:
            field_name = mtch.group(1)
            field_value = mtch.group(2)
            if field_name == "id":
                self._chk_none(rec_curr.id, lnum)
                rec_curr.id = field_value
            if field_name == "alt_id":
                rec_curr.alt_ids.append(field_value)
            elif field_name == "name":
                self._chk_none(rec_curr.name, lnum)
                rec_curr.name = field_value
            elif field_name == "namespace":
                self._chk_none(rec_curr.namespace, lnum)
                rec_curr.namespace = field_value
            elif field_name == "is_a":
                rec_curr._parents.append(field_value.split()[0])
            elif field_name == "is_obsolete" and field_value == "true":
                rec_curr.is_obsolete = True
        else:
            self._die("UNEXPECTED FIELD CONTENT: {L}\n".format(L=line), lnum)

    def _die(self, msg, lnum):
        """Raise an Exception if file read is unexpected."""
        raise Exception("**FATAL {FILE}({LNUM}): {MSG}\n".format(
            FILE=self._handle, LNUM=lnum, MSG=msg))

    def _chk_none(self, init_val, lnum):
        """Expect these lines to be uninitialized."""
        if init_val is None or init_val is "":
            return
        self._die("FIELD IS ALREADY INITIALIZED", lnum)


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
        self.depth = None  # longest distance from root node
        self.is_obsolete = False  # is_obsolete
        self.alt_ids = []  # alternative identifiers

    def __str__(self):
        obsolete = "obsolete" if self.is_obsolete else ""
        return "%s\tlevel-%02d\tdepth-%02d\t%s [%s] %s" % (self.id, self.level, self.depth, self.name, self.namespace,
                                                           obsolete)

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

    def write_hier_rec(self, gos_printed, out=sys.stdout, len_dash=1, max_depth=None, num_child=None, short_prt=False,
                       include_only=None, go_marks=None, depth=1, dp="-"):
        """Write hierarchy for a GO Term record."""
        GO_id = self.id
        # Shortens hierarchy report by only printing the hierarchy
        # for the sub-set of user-specified GO terms which are connected.
        if include_only is not None and GO_id not in include_only:
            return
        nrp = short_prt and GO_id in gos_printed
        if go_marks is not None:
            out.write('{} '.format('>' if GO_id in go_marks else ' '))
        if len_dash is not None:
            # Default character indicating hierarchy level is '-'.
            # '=' is used to indicate a hierarchical path printed in detail previously.
            letter = '-' if not nrp or not self.children else '='
            dp = ''.join([letter] * depth)
            out.write('{DASHES:{N}} '.format(DASHES=dp, N=len_dash))
        if num_child is not None:
            out.write('{N:>5} '.format(N=len(self.get_all_children())))
        out.write('{GO}\tL-{L:>02}\tD-{D:>02}\t{desc}\n'.format(
            GO=self.id, L=self.level, D=self.depth, desc=self.name))
        # Track GOs previously printed only if needed
        if short_prt:
            gos_printed.add(GO_id)
        # Do not print hierarchy below this turn if it has already been printed
        if nrp:
            return
        depth += 1
        if max_depth is not None and depth > max_depth:
            return
        for p in self.children:
            p.write_hier_rec(gos_printed, out, len_dash, max_depth, num_child, short_prt, include_only, go_marks,
                             depth, dp)


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
        log.info("Load obo file {}".format(obo_file))
        obo_reader = OBOreader(obo_file)
        for rec in obo_reader:
            self.go_Term[rec.id] = rec
            for alt in rec.alt_ids:
                self.go_Term[alt] = rec

        self.populate_terms()
        log.info("All GO nodes imported : {}".format(len(self.go_Term)))

    def populate_terms(self):
        """
        Construct go tree
        :return:
        """

        def __init_level(rec):
            """

            :param rec:
            :return:
            """
            if rec.level < 0:
                if not rec.parents:
                    rec.level = 0
                else:
                    rec.level = min(__init_level(rec) for rec in rec.parents) + 1
            return rec.level

        def _init_depth(rec):
            if rec.depth is None:
                if not rec.parents:
                    rec.depth = 0
                else:
                    rec.depth = max(_init_depth(rec) for rec in rec.parents) + 1
            return rec.depth

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
                __init_level(rec)

            if rec.depth is None:
                _init_depth(rec)

    def write_dag(self, out=sys.stdout):
        """
        Write info for all GO Terms in obo file, sorted numerically.
        """
        for rec_id, rec in sorted(self.go_Term.items()):
            print(rec, file=out)

    def write_hier_all(self, out=sys.stdout, len_dash=1, max_depth=None, num_child=None, short_prt=False):
        """
        Write hierarchy for all GO Terms in obo file.
        """
        # Print: [biological_process, molecular_function, and cellular_component]
        for go_id in ['GO:0008150', 'GO:0003674', 'GO:0005575']:
            self.write_hier(go_id, out, len_dash, max_depth, num_child, short_prt, None)

    def write_hier(self, GO_id, out=sys.stdout, len_dash=1, max_depth=None, num_child=None, short_prt=False,
                   include_only=None, go_marks=None):
        """
        Write hierarchy for a GO Term.
        """
        gos_printed = set()
        self.go_Term[GO_id].write_hier_rec(gos_printed, out, len_dash, max_depth, num_child, short_prt, include_only,
                                           go_marks)

    def write_summary_cnts(self, GO_ids, out=sys.stdout):
        """
        Write summary of level and depth counts for specific GO ids.
        """
        cnts = self.get_cnts_levels_depths_recs([self.go_Term[GO] for GO in GO_ids])
        self._write_summary_cnts(cnts, out)

    def write_summary_cnts_all(self, out=sys.stdout):
        """
        Write summary of level and depth counts for all active GO Terms.
        """
        cnts = self.get_cnts_levels_depths_recs(set(self.go_Term.values()))
        self._write_summary_cnts(cnts, out)

    def _write_summary_cnts(self, cnts, out=sys.stdout):
        """
        Write summary of level and depth counts for active GO Terms.
        """
        # Count level(shortest path to root) and depth(longest path to root)
        # values for all unique GO Terms.
        max_val = max(max(dep for dep in cnts['depth']),
                      max(lev for lev in cnts['level']))
        nss = ['biological_process', 'molecular_function', 'cellular_component']
        out.write('Dep <-Depth Counts->  <-Level Counts->\n')
        out.write('Lev   BP    MF    CC    BP    MF    CC\n')
        out.write('--- ----  ----  ----  ----  ----  ----\n')
        for i in range(max_val + 1):
            vals = ['{:>5}'.format(cnts[desc][i][ns]) for desc in cnts for ns in nss]
            out.write('{:>02} {}\n'.format(i, ' '.join(vals)))

    @staticmethod
    def get_cnts_levels_depths_recs(recs):
        """
        Collect counts of levels and depths in a Group of GO Terms.
        """
        cnts = collections.defaultdict(lambda: collections.defaultdict(collections.Counter))
        for rec in recs:
            if not rec.is_obsolete:
                cnts['level'][rec.level][rec.namespace] += 1
                cnts['depth'][rec.depth][rec.namespace] += 1
        return cnts

    @staticmethod
    def id2int(GO_id):
        return int(GO_id.replace("GO:", "", 1))

    def query_term(self, term, verbose=True):
        """
        Search term
        :param term:
        :param verbose:
        :return:
        """
        if term not in self.go_Term:
            log.error("Term {} not found!".format(term))
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
            log.error("Term {} not found!".format(term))
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
            log.info("terms not found: {}".format(bad_terms))

    def __str__(self):
        return repr(self.go_Term)

    def __repr__(self):
        return self.__str__()


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
            assoc = pd.read_csv(file, engine='c')
            assoc = assoc.dropna(axis=0)
            self._make_association(assoc)
        except Exception as e:
            log.error(e)

    def _make_association(self, assoc_data_frame):
        log.info("Making association ...")
        datagp = assoc_data_frame.groupby('id')
        for name, group in datagp:
            self.association[name] = set(group['go_id'].values)

    def query(self, term):
        """
        search a term in association
        :param term:
        :return:
        """
        return self.association[term]

    def __str__(self):
        return repr(self.association)

    def __repr__(self):
        return self.__str__()


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
            self.description = "{0} [{1}]".format(self.goterm.name, self.goterm.namespace)


class EnrichmentStudy(object):
    """
    Runs Fisher's exact test, as well as multiple corrections
    study file contain id
    pop file contain background id
    assoc file is csv format
    """

    def __init__(self, study, pop, assoc, compare=False, namespace_filter=None):
        self.compare = compare
        self.results = []
        self.study, self.pop = self.read_geneset(study, pop, compare=self.compare)
        self.association = Association(assoc)
        self.go_tree = GOtree()
        self.go_tree.update_association(self.association)

        if namespace_filter is not None:
            if namespace_filter not in ['molecular_function', 'biological_process', 'cellular_component']:
                namespace_filter = None

        self.term_study = self.count_terms(self.study, self.association, self.go_tree,
                                           namespace_filter=namespace_filter)
        self.term_pop = self.count_terms(self.pop, self.association, self.go_tree,
                                         namespace_filter=namespace_filter)
        self.pop_n, self.study_n = len(self.pop), len(self.study)

        self.run()

    def run(self):
        """
        run all
        :return:
        """
        log.info("Start go enrichement")
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

        log.info("Finished")
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
    def count_terms(geneset, assoc, go_tree, namespace_filter=None):
        """
        count the number of terms in the study group
        :param geneset:
        :param assoc:
        :param go_tree:
        :param namespace_filter: if want to filter by namespace
        """
        term_cnt = collections.defaultdict(int)
        for gene in geneset:
            try:
                for x in assoc.association[gene]:
                    if x in go_tree.go_Term:
                        if namespace_filter is not None:
                            if go_tree.go_Term[x].namespace == filter:
                                term_cnt[go_tree.go_Term[x].id] += 1
                        else:
                            term_cnt[go_tree.go_Term[x].id] += 1
            except:
                continue
        return term_cnt
