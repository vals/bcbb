import os
import sys
import glob
import shutil
from string import Template
__version__ = '0.2.2'

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
        self.template_path = os.path.join(TEMPLATEDIR, self.workflow)
        self.primersets = {}
        self.d = {'runname' :runname,
                  'samplename' : samplename,
                  'reference' : reference,
                  'basedir' : basedir,
                  'global_ini' : os.path.join(basedir, "workdir", "globals", "global.ini"),
                  'file_base' : None
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

    def enrichment_ini(self, write=True):
        tmpl = self.ini_file('enrichment.ini')
        return _write_template(os.path.join(self.basedirs['work'], 'enrichment.ini'), tmpl, write)

    def targeted_paired_end_workflow_ini(self, write=True):
        tmpl = self.ini_file('targeted.paired.end.workflow.ini')
        return _write_template(os.path.join(self.basedirs['work'], 'targeted.paired.end.workflow.ini'), tmpl, write)

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

    def __str__(self):
        return "<%s>" % (self.__class__)

class WT_SingleRead(SOLiDProject):
    """Template class for WT_SingleRead pipeline
    Required input:

      runname      - name of instrument run
      samplename   - name of sample
      reference    - full path to reference file to which mapping is done
      basedir      - sample-based directory where analysis on particular sample takes place
      file_base    - base prefix for csfasta, qual files (should be just sample name?)
      filterref    - full path to file with sequences that should be filtered out
                     e.g. adapters, rRNA, tRNA, repeats
                     Usually delivered with bioscope
      exons_gtf    - full path to file with gene models and exon definitions
      PENDING: junction_ref - full path to file with extracted exon-exon junctions
    """
    def __init__(self, runname, samplename, reference, basedir, file_base, filterref, exons_gtf, junction_ref=None, read_length=50):
        SOLiDProject.__init__(self, runname, samplename, reference, basedir)
        self.d.update({
                'read_length':read_length,
                'filter_reference':filterref,
                'exons_gtf':exons_gtf,
                'file_base':file_base
                # 'junction_reference':junction_ref
                })
        self.primersets['F3'] = Primer("F3", read_length, self)

    def wt_single_read_ini(self):
        return self.ini_file('wt.single.read.workflow.ini')

    def init_project(self):
        analysis_plan = os.path.join(self.basedirs['base'], 'analysis.plan')
        ap = []
        self.global_ini()
        self.primersets['F3'].wt_single_read_workflow_ini()
        ap.append(os.path.join(self.primersets['F3'].dirs['work'], 'F3.ini'))
        with open(analysis_plan, 'w') as apf:
            apf.write("\n".join(ap))


# NOTE: the saet_target_file is a dummy file containing only one line
# indicating the target region size - otherwise bioscope crashes...
class TargetedFrag(SOLiDProject):
    def __init__(self, runname, samplename, reference, basedir, targetfile, saettargetfile, cmap, annotation_gtf_file=None, read_length=50, annotation_human_hg18=0, annotation_dbsnp_file_snpchrpos=None, annotation_dbsnp_file_snpcontigloc=None, annotation_dbsnp_file_snpcontiglocusid=None):
        SOLiDProject.__init__(self, runname, samplename, reference, basedir)
        _key_map = self._key_map.update({'cmap':'cmap', 'annotation_gtf_file':'annotation.gtf.file', 'annotation_dbsnp_file_snpchrpos':'annotation.dbsnp.file.snpchrpos', 'annotation_dbsnp_file_snpcontigloc':'annotation.dbsnp.file.snpcontigloc', 'annotation_dbsnp_file_snpcontiglocusid':'annotation.dbsnp.file.snpcontiglocusid'})
        self.config.update({
                'annotation_gtf_file':annotation_gtf_file,
                'annotation_dbsnp_file_snpchrpos':annotation_dbsnp_file_snpchrpos,
                'annotation_dbsnp_file_snpcontigloc':annotation_dbsnp_file_snpcontigloc,
                'annotation_dbsnp_file_snpcontiglocusid':annotation_dbsnp_file_snpcontiglocusid
                })
        self.d.update( {
                'cmap' : cmap,
                'target_file' : targetfile,
                'saet_target_file' : saettargetfile,
                'annotation_human_hg18' : annotation_human_hg18,
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
    def __init__(self, primer, readlength, project, small_indel_frag_run=1):
        self.project = project
        self.dirs = {'work': os.path.join(self.project.basedirs['work'],  primer + "_mapping"),
                     'output':os.path.join(self.project.basedirs['output'], primer),
                     'reads':os.path.join(self.project.basedirs['reads'],  primer)
                     }
        _make_dirs(self.dirs)
        self.d = {'primer': primer,
                  'primerlabel': primer.lower()[0:2],
                  'read_length':readlength,
                  'csfastafilebase' : self.project.d['samplename'] + "_" + primer + ".csfasta",
                  'saet_input_csfastafile' : os.path.join(self.dirs['reads'], str(self.project.d['file_base']) + "_" + primer + ".csfasta"),
                  'saet_input_qualfile' : os.path.join(self.dirs['reads'], str(self.project.d['file_base']) + "_" + primer + "_QV.qual"),
                  'csfastafile' : os.path.join(self.dirs['reads'], str(self.project.d['file_base']) + "_" + primer + ".csfasta"),
                  'qualfile' : os.path.join(self.dirs['reads'], str(self.project.d['file_base']) + "_" + primer + "_QV.qual"),
                  'matobamqual' : self.project.d['samplename'] + "_" + primer + "_QV.qual",
                  # 'small_indel_frag_run',
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

    def ma_to_bam_ini(self, write=True):
        tmpl = self.ini_file('ma.to.bam.ini')
        return _write_template(os.path.join(self.dirs['work'], 'ma.to.bam.ini'), tmpl, write)

    def primer_ini(self, write=True):
        tmpl = self.ini_file('primer.ini')
        return _write_template(os.path.join(self.dirs['work'], "%s.ini" % (self.d['primer'])), tmpl, write)

    def wt_single_read_workflow_ini(self, write=True):
        tmpl = self.ini_file('wt.single.read.workflow.ini')
        return _write_template(os.path.join(self.dirs['work'], "%s.ini" % (self.d['primer'])), tmpl, write)

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
    def __init__(self, runname, samplename, reference, basedir, targetfile, saettargetfile, cmap, file_base, annotation_gtf_file=None, read_length=[50, 35], annotation_human_hg18=0, primersetlabels=["F3", "F5-BC"], annotation_dbsnp_file_snpchrpos=None, annotation_dbsnp_file_snpcontigloc=None, annotation_dbsnp_file_snpcontiglocusid=None):
        SOLiDProject.__init__(self, runname, samplename, reference, basedir)
        _key_map = self._key_map.update({'cmap':'cmap', 'annotation_gtf_file':'annotation.gtf.file', 'annotation_dbsnp_file_snpchrpos':'annotation.dbsnp.file.snpchrpos', 'annotation_dbsnp_file_snpcontigloc':'annotation.dbsnp.file.snpcontigloc', 'annotation_dbsnp_file_snpcontiglocusid':'annotation.dbsnp.file.snpcontiglocusid'})
        self.config.update({
                'annotation_gtf_file':annotation_gtf_file,
                'annotation_dbsnp_file_snpchrpos':annotation_dbsnp_file_snpchrpos,
                'annotation_dbsnp_file_snpcontigloc':annotation_dbsnp_file_snpcontigloc,
                'annotation_dbsnp_file_snpcontiglocusid':annotation_dbsnp_file_snpcontiglocusid
                })
        self.d.update( {
                'cmap' : cmap,
                'target_file' : targetfile,
                'saet_target_file' : saettargetfile,
                'annotation_human_hg18' : annotation_human_hg18,
                'file_base' : file_base,
                'primer1' : primersetlabels[0],
                'primer2' : primersetlabels[1],
                'primerset' : ",".join(primersetlabels),
                'maximal_read_length' : max(read_length)
                } )
        self.d.update(self._set_d())
        self.primersets[primersetlabels[0]] = Primer(primersetlabels[0], read_length[0], self)
        self.primersets[primersetlabels[1]] = Primer(primersetlabels[1], read_length[1], self, 0)

    def primerset_global(self):
        self.d.update({
                       'csfastafilebase':self.primersets['F3'].d['csfastafilebase']
                       })

    def init_project(self, saet=True, ma_to_bam=True, enrichment=True, targeted_workflow=True):
        analysis_plan = os.path.join(self.basedirs['base'], 'analysis.plan')
        ap = []
        self.global_ini()
        primers = self.primersets.keys()
        if saet:
            self.primersets[primers[0]].saet_ini()
            self.primersets[primers[1]].saet_ini()
            ap.append("=%s" %(os.path.join(self.primersets[primers[0]].dirs['work'], 'saet.ini')))
            ap.append("=%s" %(os.path.join(self.primersets[primers[1]].dirs['work'], 'saet.ini')))
            ap.append("\n\n+\n")
        self.primersets[primers[0]].primer_ini()
        self.primersets[primers[1]].primer_ini()
        ap.append("=%s" %(os.path.join(self.primersets[primers[0]].dirs['work'], "%s.ini" %(primers[0]))))
        ap.append("=%s" %(os.path.join(self.primersets[primers[1]].dirs['work'], "%s.ini" %(primers[1]))))
        ap.append("\n\n")
        if ma_to_bam:
            self.primersets[primers[0]].ma_to_bam_ini()
            self.primersets[primers[1]].ma_to_bam_ini()
            ap.append("%s" %(os.path.join(self.primersets[primers[0]].dirs['work'], 'ma.to.bam.ini')))
            ap.append("%s" %(os.path.join(self.primersets[primers[1]].dirs['work'], 'ma.to.bam.ini')))
        if enrichment:
            self.enrichment_ini()
            ap.append(os.path.join(self.basedirs['work'], 'enrichment.ini'))
        if targeted_workflow:
            self.targeted_paired_end_workflow_ini()
            ap.append(os.path.join(self.basedirs['work'], 'targeted.paired.end.workflow.ini'))
        with open(analysis_plan, 'w') as apf:
            apf.write("\n".join(ap))


class ReseqFrag(SOLiDProject):
    def __init__(self):
        pass



