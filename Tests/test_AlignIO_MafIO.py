# Copyright 2016 - Tristan Bitard-Feildel.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for Bio.AlignIO.MafIO"""

import unittest

from Bio._py3k import StringIO

from Bio.AlignIO.MafIO import MafIterator, MafWriter


# from https://genome.ucsc.edu/FAQ/FAQformat.html#format5
maf_text = """track name=euArc visibility=pack
##maf version=1 scoring=tba.v8 
# tba.v8 (((human chimp) baboon) (mouse rat)) 
                   
a score=23262.0     
s hg18.chr7    27578828 38 + 158545518 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s panTro1.chr6 28741140 38 + 161576975 AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG
s baboon         116834 38 +   4622798 AAA-GGGAATGTTAACCAAATGA---GTTGTCTCTTATGGTG
s mm4.chr6     53215344 38 + 151104725 -AATGGGAATGTTAAGCAAACGA---ATTGTCTCTCAGTGTG
s rn3.chr4     81344243 40 + 187371129 -AA-GGGGATGCTAAGCCAATGAGTTGTTGTCTCTCAATGTG
                   
a score=5062.0                    
s hg18.chr7    27699739 6 + 158545518 TAAAGA
s panTro1.chr6 28862317 6 + 161576975 TAAAGA
s baboon         241163 6 +   4622798 TAAAGA 
s mm4.chr6     53303881 6 + 151104725 TAAAGA
s rn3.chr4     81444246 6 + 187371129 taagga

a score=6636.0
s hg18.chr7    27707221 13 + 158545518 gcagctgaaaaca
s panTro1.chr6 28869787 13 + 161576975 gcagctgaaaaca
s baboon         249182 13 +   4622798 gcagctgaaaaca
s mm4.chr6     53310102 13 + 151104725 ACAGCTGAAAATA

"""

# maf text2 from EPO ?


# maf text3 from UCSC Mouse chr10 first two alignment
maf_text2 = """
##maf version=1 scoring=autoMZ.v1                                                                                                        
a score=11761.000000                                                                                                                     
s hg19.chr10      60000 45 + 135534747 GAATTCCTTGAGGCCTAAATGCATCGGGGTGCTCTGGTTTTGTTG                                                     
s ponAbe2.chrUn 3569800 45 +  72422247 GAATTCCTTGAGGCCTAATTGCATCAGGGTGCTCTGGTTTTGTTG                                                     
q ponAbe2.chrUn                        999999999999999999999999999999999999999999999                                                     
i ponAbe2.chrUn N 0 C 0                                                                                                                  
s panTro2.chr18    5870 45 +  77261746 GAATTCCTTGAGGCCTAAATGCATCGGGGAGCTCTGGTTTTGTTG                                                     
q panTro2.chr18                        999999888999978999999999889979967847889999999                                                     
i panTro2.chr18 N 0 C 0

a score=45884.000000                                                                                                                     
s hg19.chr10           60045 108 + 135534747 TTGTTATTTCTGAATGACATTTACTTTGGTGCTCTTTATTTTGCGTATTTAAAACTATTAGATCGTGTGATTATATTTGACAGGTCTTAATTGACGCGCTGTTCAGCC                                                                                                                         
s panTro2.chr18         5915 108 +  77261746 TTGTTATTTCTGAATGACATTGACTTTGGTGCTTTTTATTTTGCATATTTAAAACTATTAGATCGTGTGATTATATTTGACAAGTCTTAATTGACGCGCTGTTCAGCG                                                                                                                         
q panTro2.chr18                              999999999999999999999999999995999999999999988999999999999999999963999999999999999999999999999999999999999999                                                                                                                         
i panTro2.chr18      C 0 C 0                                                                                                             
s ponAbe2.chrUn      3569845 108 +  72422247 TTGCTATTTCTGAATGACATTGACTTTGGTGCTCTTTATTTTGCATATTTAAAACTATTAGATCATGTGATTATATTTGACAGGTCTTAATTGATGCACTCTTCAGCG                                                                                                                         
q ponAbe2.chrUn                              999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999                                                                                                                         
i ponAbe2.chrUn      C 0 C 0                                                                                                             
s calJac1.Contig5564   47630 102 -    188978 TTGCTATTTCTGAGTGACATTGACTTTGGTGCACTTTATTTGGAATATTTAAAA-TATTAGATTGTG---TTGTATTTGAAAGGCCTTAATTGATGTGCTGTTCAG--                                                                                                                         
q calJac1.Contig5564                         999999999999999999999999999999999999999999999999999999-999999999999---999999999999999999999999999999999999--                                                                                                                         
i calJac1.Contig5564 N 0 I 14  

##eof maf
"""

class TestMafIO(unittest.TestCase):

    def test_one(self):
        test_ids = {0:['hg18.chr7', 'panTro1.chr6', 'baboon', 'mm4.chr6', 'rn3.chr4'],
                    1:['hg18.chr7', 'panTro1.chr6', 'baboon', 'mm4.chr6', 'rn3.chr4'],
                    2:['hg18.chr7', 'panTro1.chr6', 'baboon', 'mm4.chr6'],}

        test_seqs = {0:"-AA-GGGGATGCTAAGCCAATGAGTTGTTGTCTCTCAATGTG",
                    1:"taagga",
                    2:"ACAGCTGAAAATA"}
        
        handle = StringIO(maf_text)
        for i, alignment in enumerate(MafIterator(handle)):
            ids = []
            print(i, alignment._block_lines)
            for record in alignment:
                ids.append(record.id)
            self.assertEqual(ids, test_ids[i])

            expected = test_seqs[i].upper()
            self.assertEqual(str(record.seq).upper(), expected)

    def test_two(self):
        handle = StringIO(maf_text2)
        list2 = list(MafIterator(handle))
        handle.close()

        self.assertEqual(len(list2), 2)
        self.assertEqual(len(list2[0]), 3)

    def test_write_read(self):
        handle = StringIO(maf_text)
        list5 = list(MafIterator(handle))
        handle.close()

        handle = StringIO()
        MafWriter(handle).write_file(list5)
        handle.seek(0)
        list6 = list(MafIterator(handle))

        self.assertEqual(len(list5), len(list6))
        for a1, a2 in zip(list5, list6):
            self.assertEqual(len(a1), len(a2))
            for r1, r2 in zip(a1, a2):
                self.assertEqual(r1.id, r2.id)
                self.assertEqual(str(r1.seq), str(r2.seq))


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
