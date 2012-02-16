"""Manage and automate the Solexa read processing pipeline.
"""

# These index number-to-sequence mappings have been double-checked
# against the documentation from Illumina dated 2011-10-11.
# NOTE: index1-index27 are from the table
#       "TruSeq RNA and DNA Sample Prep Kits".
# NOTE: r1-r48 are from the table "TruSeq Small RNA Sample Prep Kits",
#       after reverse-complement conversion.
INDEX_LOOKUP = dict(index1='ATCACG',
                    index2='CGATGT',
                    index3='TTAGGC',
                    index4='TGACCA',
                    index5='ACAGTG',
                    index6='GCCAAT',
                    index7='CAGATC',
                    index8='ACTTGA',
                    index9='GATCAG',
                    index10='TAGCTT',
                    index11='GGCTAC',
                    index12='CTTGTA',
                    index13='AGTCAA',
                    index14='AGTTCC',
                    index15='ATGTCA',
                    index16='CCGTCC',
                    # index17 is "reserved" by Illumina
                    index18='GTCCGC',
                    index19='GTGAAA',
                    index20='GTGGCC',
                    index21='GTTTCG',
                    index22='CGTACG',
                    index23='GAGTGG',
                    # index24 is "reserved" by Illumina
                    index25='ACTGAT',
                    # index26 is "reserved" by Illumina
                    index27='ATTCCT',
                    # RPI indexes for "TruSeq Small RNA",
                    # These are reverse-complement of Illumina documentation
                    rpi1='ATCACG',
                    rpi2='CGATGT',
                    rpi3='TTAGGC',
                    rpi4='TGACCA',
                    rpi5='ACAGTG',
                    rpi6='GCCAAT',
                    rpi7='CAGATC',
                    rpi8='ACTTGA',
                    rpi9='GATCAG',
                    rpi10='TAGCTT',
                    rpi11='GGCTAC',
                    rpi12='CTTGTA',
                    rpi13='AGTCAA',
                    rpi14='AGTTCC',
                    rpi15='ATGTCA',
                    rpi16='CCGTCC',
                    rpi17='GTAGAG',
                    rpi18='GTCCGC',
                    rpi19='GTGAAA',
                    rpi20='GTGGCC',
                    rpi21='GTTTCG',
                    rpi22='CGTACG',
                    rpi23='GAGTGG',
                    rpi24='GGTAGC',
                    rpi25='ACTGAT',
                    rpi26='ATGAGC',
                    rpi27='ATTCCT',
                    rpi28='CAAAAG',
                    rpi29='CAACTA',
                    rpi30='CACCGG',
                    rpi31='CACGAT',
                    rpi32='CACTCA',
                    rpi33='CAGGCG',
                    rpi34='CATGGC',
                    rpi35='CATTTT',
                    rpi36='CAAACA',
                    rpi37='CGGAAT',
                    rpi38='CTAGCT',
                    rpi39='CTATAC',
                    rpi40='CTCAGA',
                    rpi41='GACGAC',
                    rpi42='TAATCG',
                    rpi43='TACAGC',
                    rpi44='TATAAT',
                    rpi45='TCATTC',
                    rpi46='TCCCGA',
                    rpi47='TCGAAG',
                    rpi48='TCGGCA')

INDEX_LOOKUP.update(dict([(k.replace('index', ''), v)
                          for k, v in INDEX_LOOKUP.items()]))
INDEX_LOOKUP.update(dict([(k.replace('index', 'idx'), v)
                          for k, v in INDEX_LOOKUP.items()]))
INDEX_LOOKUP.update(dict([(k.replace('index', 'in'), v)
                          for k, v in INDEX_LOOKUP.items()]))
INDEX_LOOKUP.update(dict([(k.replace('rpi', 'r'), v)
                          for k, v in INDEX_LOOKUP.items()]))
INDEX_LOOKUP.update(dict([(k.upper(), v)
                          for k, v in INDEX_LOOKUP.items()]))
