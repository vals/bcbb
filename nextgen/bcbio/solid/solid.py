import os
import sys
import glob
import shutil
from string import Template
__version__ = '0.2'

TEMPLATEDIR= os.path.join(os.path.dirname(__file__), "templates")

# Utility functions, not for export
def _make_dirs(d):
    for k in d.keys():
        if not os.path.exists(d[k]):
            os.makedirs(d[k])

def _write_template(fo, tmpl, write=True):
    """Write the template"""
    if write:
        if not os.path.exists(fo):
            fp = open(fo, "w")
            fp.write(tmpl)
            fp.close()
        return 
    else:
        return tmpl

def _clean_up(d, key):
    paths = []
    obj = os.listdir(d[key])
    for o in obj:
        if not o.endswith(".ini"):
            p = os.path.join(d[key], o)
            if os.path.isdir(p):
                shutil.rmtree(p)
            else:
                os.remove(p)
            paths.append(p)
    return "\n".join(paths)


class SOLiDProject(object):
    """Template class for SOLiD projects"""
    _key_map = {}
    def __init__(self, runname, samplename, reference, basedir):
        self.config = {} # For keys that can be None
        self.workflow = self.__class__.__name__
        self.basedirs = {'base' : basedir,
                         'work': os.path.join(basedir, "workdir"),
                         'output' : os.path.join(basedir, "output"),
                         'reads' : os.path.join(basedir, "reads"),
                         'log' : os.path.join(basedir, "log"),
                         'temp' : os.path.join(basedir, "temp"),
                         'globals' : os.path.join(basedir, "workdir", "globals"),
                         'intermediate' : os.path.join(basedir, "intermediate")
                         }
        # self.workdirs = {'dibayes' : os.path.join(basedir, "workdir", "dibayes"),
        #                  'positionErrors' : os.path.join(basedir, "workdir", "positionErrors")
        #                  }
        # self.outdirs = {'dibayes' : os.path.join(basedir, "output", "dibayes"),
        #                 'positionErrors' : os.path.join(basedir, "output", "positionErrors")
        #                 }
        self.template_path = os.path.join(TEMPLATEDIR, self.workflow)
        self.primersets = {}
        self.d = {'runname' :runname,
                  'samplename' : samplename,
                  'reference' : reference,
                  'basedir' : basedir,
                  'global_ini' : os.path.join(basedir, "workdir", "globals", "global.ini")
                  }
        _make_dirs(self.basedirs)
        
    def _set_d(self):
        d = {}
        for k in self.config.keys():
            d[k] = " = ".join([self._key_map[k], str(self.config[k])])
            if self.config[k] == None:
                d[k] = "# " + d[k]
        return d

    def ini_file(self, filename):
        inifile = os.path.join(self.template_path, filename)
        with open(inifile) as in_handle:
            tmpl = Template(in_handle.read())
        return tmpl.safe_substitute(self.d)

    def primerset_global(self):
        pass

    def global_ini(self, write=True):
        self.primerset_global()
        tmpl = self.ini_file('global.ini')
        return _write_template(os.path.join(self.basedirs['globals'], 'global.ini'), tmpl, write)
    
    def clean(self, verbose=False):
        pstr1 = _clean_up(self.basedirs, 'log')
        pstr2 = _clean_up(self.basedirs, 'temp')
        pstr3 = _clean_up(self.basedirs, 'intermediate') 
        pstr = pstr1 + pstr2 + pstr3
        if verbose and pstr != "":
            print "removing files for sample " + self.d['samplename']
            print pstr
        for k, p in self.primersets.items():
            p.clean(verbose)

class WT_SingleRead(SOLiDProject):
    def __init__(self, runname, samplename, reference, basedir, csfastafile, qualfile, filterref, exons_gtf, junction_ref, read_length=50):
        SOLiDProject.__init__(self, runname, samplename, reference, basedir)
        self.d.update({
                'read_length':read_length,
                'csfastafile':csfastafile,
                'qualfile':qualfile, 
                'filter_reference':filterref,
                'exons_gtf':exons_gtf,
                'junction_reference':junction_ref
                })

    def wt_single_read_ini(self):
        return self.ini_file('wt.single.read.workflow.ini')

# NOTE: the saet_target_file is a dummy file containing only one line
# indicating the target region size - otherwise bioscope crashes...
class TargetedFrag(SOLiDProject):
    def __init__(self, runname, samplename, reference, basedir, targetfile, saettargetfile, cmap, annotation_gtf_file=None, read_length=50, annotation_human_hg18=0):
        SOLiDProject.__init__(self, runname, samplename, reference, basedir)
        _key_map = self._key_map.update({'cmap':'cmap', 'annotation_gtf_file':'annotation.gtf.file'})
        self.config.update({
                'annotation_gtf_file':annotation_gtf_file
                })
        self.d.update( {
                'cmap' : cmap,
                'target_file' : targetfile,
                'saet_target_file' : saettargetfile,
                'annotation_human_hg18' : annotation_human_hg18
                } )
        self.d.update(self._set_d())
        self.primersets['F3'] = Primer("F3", read_length, self)

    def primerset_global(self):
        self.d.update({'read_length':self.primersets['F3'].d['read_length'],
                       'csfastafilebase':self.primersets['F3'].d['csfastafilebase']
                       })

    def init_project(self, saet=True, small_indel_frag=True, enrichment=True, targeted_workflow=True):
        analysis_plan = os.path.join(self.basedirs['base'], 'analysis.plan')
        ap = []
        self.global_ini()
        if saet:
            self.primersets['F3'].saet_ini()
            ap.append(os.path.join(self.primersets['F3'].dirs['work'], 'saet.ini'))
        if small_indel_frag:
            self.primersets['F3'].enrichment_ini()
            ap.append(os.path.join(self.primersets['F3'].dirs['work'], 'small.indel.frag.ini'))
        if enrichment:
            self.primersets['F3'].targeted_frag_workflow_ini()
            ap.append(os.path.join(self.primersets['F3'].dirs['work'], 'enrichment.ini'))
        if targeted_workflow:
            self.primersets['F3'].small_indel_frag_ini()
            ap.append(os.path.join(self.primersets['F3'].dirs['work'], 'targeted.frag.workflow.ini'))
        with open(analysis_plan, 'w') as apf:
            apf.write("\n".join(ap))
                          

class Primer(object):
    """Class for primer set"""
    def __init__(self, primer, readlength, project):
        self.project = project
        self.dirs = {'work': os.path.join(self.project.basedirs['work'],  primer + "_mapping"),
                     'output':os.path.join(self.project.basedirs['output'], primer + "_mapping"),
                     'reads':os.path.join(self.project.basedirs['reads'],  primer)
                     }
        _make_dirs(self.dirs)
        self.d = {'primer': primer,
                  'primerlabel': primer.lower(),
                  'read_length':readlength,
                  'csfastafilebase' : self.project.d['samplename'] + "_" + primer + ".csfasta",
                  'saet_input_csfastafile' : os.path.join(self.dirs['reads'], self.project.d['samplename'] + "_" + primer + ".csfasta"),
                  'saet_input_qualfile' : os.path.join(self.dirs['reads'], self.project.d['samplename'] + "_" + primer + "_QV.qual"),
                  'matobamqual' : self.project.d['samplename'] + "_" + primer + "_QV.qual" 
                  # As of yet I have no idea what this looks like
                  # 'small_indel_frag_qual' : self.project
                  }

    def saet_ini(self, write=True):
        tmpl = self.ini_file('saet.ini')
        return _write_template(os.path.join(self.dirs['work'], 'saet.ini'), tmpl, write)

    def enrichment_ini(self, write=True):
        tmpl = self.ini_file('enrichment.ini')
        return _write_template(os.path.join(self.dirs['work'], 'enrichment.ini'), tmpl, write)

    def targeted_frag_workflow_ini(self, write=True):
        tmpl = self.ini_file('targeted.frag.workflow.ini')
        return _write_template(os.path.join(self.dirs['work'], 'targeted.frag.workflow.ini'), tmpl, write)

    def small_indel_frag_ini(self, write=True):
        tmpl = self.ini_file('small.indel.frag.ini')
        return _write_template(os.path.join(self.dirs['work'], 'small.indel.frag.ini'), tmpl, write)
    
    def ini_file(self, filename):
        inifile = os.path.join(self.project.template_path, filename)
        with open(inifile) as in_handle:
            tmpl = Template(in_handle.read())
        # Global project dictionary
        d = self.project.d
        # Primer specific dictionary
        d.update(self.d)
        return tmpl.safe_substitute(d)

    def clean(self, verbose=False):
        pstr1 = _clean_up(self.dirs, 'work')
        pstr2 = _clean_up(self.dirs, 'output')
        pstr = pstr1 + pstr2
        if verbose and pstr != "":
            print "removing files for primer set " + self.d['primer']
            print pstr
        

class TargetedPE(SOLiDProject):
    def __init__(self):
        pass

class ReseqFrag(SOLiDProject):
    def __init__(self):
        pass



