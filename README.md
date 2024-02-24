# AdvancedComputationalBiologyFinalPortfolio_PaytonMcAlphin
Final for Advanced Computational Biology, Payton McAlphin

## Sequence Objects (Parts 1-4)
```python
from Bio.Seq import Seq
```


```python
my_seq = Seq("GATCG")
```


```python
for index, letter in enumerate(my_seq):
    print("%i %s" % (index, letter))
```

    0 G
    1 A
    2 T
    3 C
    4 G



```python
# We can also print the length of each sequence
print(len(my_seq))
```

    5



```python
print(my_seq[0])
```

    G



```python
print(my_seq[4])
```

    G



```python
print(my_seq[2])
```

    T



```python
Seq("AAAA").count("AA")
```




    2




```python
my_seq = Seq("GATCGATGGGCCTATATAGGATCGAAAATCGC")
```


```python
len(my_seq)
```




    32




```python
my_seq.count("G")
```




    9




```python
100 * (my_seq.count("G") + my_seq.count("C")) / len(my_seq)
```




    46.875




```python
from Bio.SeqUtils import gc_fraction
```


```python
gc_fraction(my_seq)
```




    0.46875




```python
 my_seq[4:12]
```




    Seq('GATGGGCC')




```python
my_seq[0::3]
```




    Seq('GCTGTAGTAAG')




```python
my_seq[1::3]
```




    Seq('AGGCATGCATC')




```python
my_seq[2:3]
```




    Seq('T')




```python
my_seq[::-1]
```




    Seq('CGCTAAAAGCTAGGATATATCCGGGTAGCTAG')




```python
str(my_seq)
```




    'GATCGATGGGCCTATATAGGATCGAAAATCGC'




```python
fasta_format_string = ">Name\n%s\n" % my_seq
```


```python
print(fasta_format_string)
```

    >Name
    GATCGATGGGCCTATATAGGATCGAAAATCGC
    



```python
seq1 = Seq("ACGT")
seq2 = Seq("AACCGG")
```


```python
seq1 + seq2
```




    Seq('ACGTAACCGG')




```python
seq2 + seq1
```




    Seq('AACCGGACGT')




```python
contigs = [Seq("ATG"), Seq("ATCCCG"), Seq("TTGCA")]
```


```python
spacer = Seq("N" * 10)
```


```python
spacer.join(contigs)
```




    Seq('ATGNNNNNNNNNNATCCCGNNNNNNNNNNTTGCA')




```python
dna_seq = Seq("acgtACGT")
```


```python
dna_seq
```




    Seq('acgtACGT')




```python
dna_seq.upper()
```




    Seq('ACGTACGT')




```python
dna_seq.lower()
```




    Seq('acgtacgt')




```python
"gtac" in dna_seq
```




    False




```python
"GTAC" in dna_seq
```




    False




```python
dna_seq
```




    Seq('acgtACGT')




```python
dna_seq = dna_seq.upper()
```


```python
"GTAC" in dna_seq
```




    True




```python
my_seq = Seq("GATCGATGGGCCTATATAGGATCGAAAATCGC")
```


```python
my_seq.complement()
```




    Seq('CTAGCTACCCGGATATATCCTAGCTTTTAGCG')




```python
my_seq.reverse_complement()
```




    Seq('GCGATTTTCGATCCTATATAGGCCCATCGATC')




```python
protein_seq = Seq("EVRNAK")
protein_seq.complement()
```




    Seq('EBYNTM')




```python
coding_dna = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
```


```python
coding_dna
```




    Seq('ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG')




```python
template_dna = coding_dna.reverse_complement()
```


```python
template_dna
```




    Seq('CTATCGGGCACCCTTTCAGCGGCCCATTACAATGGCCAT')




```python
coding_dna
```




    Seq('ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG')




```python
messenger_rna = coding_dna.transcribe()
```


```python
messenger_rna
```




    Seq('AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG')




```python
template_dna.reverse_complement().transcribe()
```




    Seq('AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG')




```python
messenger_rna.back_transcribe()
```




    Seq('ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG')




```python
messenger_rna
```




    Seq('AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG')




```python
messenger_rna.translate()
```




    Seq('MAIVMGR*KGAR*')




```python
coding_dna.translate(table="Vertebrate Mitochondrial")
```




    Seq('MAIVMGRWKGAR*')




```python
coding_dna.translate(table = 2)
```




    Seq('MAIVMGRWKGAR*')




```python
coding_dna.translate(to_stop = True)
```




    Seq('MAIVMGR')




```python
coding_dna.translate(table =2, to_stop=True)
```




    Seq('MAIVMGRWKGAR')




```python
coding_dna.translate(table = 2, stop_symbol = "!")
```




    Seq('MAIVMGRWKGAR!')




```python
gene = Seq("GTGAAAAAGATGCAATCTATCGTACTCGCACTTTCCCTGGTTCTGGTCGCTCCCATGGCAGCACAGGCTGCGGAAATTACGTTAGTCCCGTCAGTAAAATTACAGATAGGCGATCGTGATAATCGTGGCTATTACTGGGATGGAGGTCACTGGCGCGACCACGGCTGGTGGAAACAACATTATGAATGGCGAGGCAATCGCTGGCACCTACACGGACCGCCGCCACCGCCGCGCCACCATAAGAAAGCTCCTCATGATCATCACGGCGGTCATGGTCCAGGCAAACATCACCGCTAA")
```


```python
gene.translate(table = "Bacterial")
```




    Seq('VKKMQSIVLALSLVLVAPMAAQAAEITLVPSVKLQIGDRDNRGYYWDGGHWRDH...HR*')




```python
gene.translate(table = "Bacterial", to_stop = True)
```




    Seq('VKKMQSIVLALSLVLVAPMAAQAAEITLVPSVKLQIGDRDNRGYYWDGGHWRDH...HHR')




```python
gene.translate(table = "Bacterial", cds = True)
```




    Seq('MKKMQSIVLALSLVLVAPMAAQAAEITLVPSVKLQIGDRDNRGYYWDGGHWRDH...HHR')




```python
from Bio.Data import CodonTable
```


```python
standard_table = CodonTable.unambiguous_dna_by_name["Standard"]
```


```python
mito_table = CodonTable.unambiguous_dna_by_name["Vertebrate Mitochondrial"]
```


```python
print(standard_table)
```

    Table 1 Standard, SGC0
    
      |  T      |  C      |  A      |  G      |
    --+---------+---------+---------+---------+--
    T | TTT F   | TCT S   | TAT Y   | TGT C   | T
    T | TTC F   | TCC S   | TAC Y   | TGC C   | C
    T | TTA L   | TCA S   | TAA Stop| TGA Stop| A
    T | TTG L(s)| TCG S   | TAG Stop| TGG W   | G
    --+---------+---------+---------+---------+--
    C | CTT L   | CCT P   | CAT H   | CGT R   | T
    C | CTC L   | CCC P   | CAC H   | CGC R   | C
    C | CTA L   | CCA P   | CAA Q   | CGA R   | A
    C | CTG L(s)| CCG P   | CAG Q   | CGG R   | G
    --+---------+---------+---------+---------+--
    A | ATT I   | ACT T   | AAT N   | AGT S   | T
    A | ATC I   | ACC T   | AAC N   | AGC S   | C
    A | ATA I   | ACA T   | AAA K   | AGA R   | A
    A | ATG M(s)| ACG T   | AAG K   | AGG R   | G
    --+---------+---------+---------+---------+--
    G | GTT V   | GCT A   | GAT D   | GGT G   | T
    G | GTC V   | GCC A   | GAC D   | GGC G   | C
    G | GTA V   | GCA A   | GAA E   | GGA G   | A
    G | GTG V   | GCG A   | GAG E   | GGG G   | G
    --+---------+---------+---------+---------+--



```python
print(mito_table)
```

    Table 2 Vertebrate Mitochondrial, SGC1
    
      |  T      |  C      |  A      |  G      |
    --+---------+---------+---------+---------+--
    T | TTT F   | TCT S   | TAT Y   | TGT C   | T
    T | TTC F   | TCC S   | TAC Y   | TGC C   | C
    T | TTA L   | TCA S   | TAA Stop| TGA W   | A
    T | TTG L   | TCG S   | TAG Stop| TGG W   | G
    --+---------+---------+---------+---------+--
    C | CTT L   | CCT P   | CAT H   | CGT R   | T
    C | CTC L   | CCC P   | CAC H   | CGC R   | C
    C | CTA L   | CCA P   | CAA Q   | CGA R   | A
    C | CTG L   | CCG P   | CAG Q   | CGG R   | G
    --+---------+---------+---------+---------+--
    A | ATT I(s)| ACT T   | AAT N   | AGT S   | T
    A | ATC I(s)| ACC T   | AAC N   | AGC S   | C
    A | ATA M(s)| ACA T   | AAA K   | AGA Stop| A
    A | ATG M(s)| ACG T   | AAG K   | AGG Stop| G
    --+---------+---------+---------+---------+--
    G | GTT V   | GCT A   | GAT D   | GGT G   | T
    G | GTC V   | GCC A   | GAC D   | GGC G   | C
    G | GTA V   | GCA A   | GAA E   | GGA G   | A
    G | GTG V(s)| GCG A   | GAG E   | GGG G   | G
    --+---------+---------+---------+---------+--



```python
mito_table.stop_codons
```




    ['TAA', 'TAG', 'AGA', 'AGG']




```python
mito_table.start_codons
```




    ['ATT', 'ATC', 'ATA', 'ATG', 'GTG']




```python
seq = Seq("ACGT")
```


```python
"ACGT" == seq1
```




    True




```python
seq1 == "ACGT"
```




    True




```python
unknown_seq = Seq(None, 10)
```


```python
unknown_seq
```




    Seq(None, length=10)




```python
len(unknown_seq)
```




    10




```python
seq = Seq({117512683: "TTGAAAACCTGAATGTGAGAGTCAGTCAAGGATAGT"}, length = 159345973)
```


```python
seq[1000:1020]
```




    Seq(None, length=20)




```python
seq[117512690:117512700]
```




    Seq('CCTGAATGTG')




```python
seq[117512670:]
```




    Seq({13: 'TTGAAAACCTGAATGTGAGAGTCAGTCAAGGATAGT'}, length=41833303)




```python
seq = Seq("ACGT")
```


```python
undefined_seq = Seq(None, length =10)
```


```python
seq + undefined_seq + seq
```




    Seq({0: 'ACGT', 14: 'ACGT'}, length=18)




```python
my_seq = Seq("GCCATTGTAATGGGCCGCTGAAAGGGTGCCCGA")
```


```python
from Bio.Seq import MutableSeq
```


```python
mutable_seq = MutableSeq(my_seq)
```


```python
mutable_seq
```




    MutableSeq('GCCATTGTAATGGGCCGCTGAAAGGGTGCCCGA')




```python
mutable_seq[5] = "C"
```


```python
mutable_seq
```




    MutableSeq('GCCATCGTAATGGGCCGCTGAAAGGGTGCCCGA')




```python
mutable_seq.remove("T")
```


```python
mutable_seq
```




    MutableSeq('GCCACGTAATGGGCCGCTGAAAGGGTGCCCGA')




```python
mutable_seq.reverse()
```


```python
mutable_seq
```




    MutableSeq('AGCCCGTGGGAAAGTCGCCGGGTAATGCACCG')




```python
new_seq = Seq(mutable_seq)
```


```python
new_seq
```




    Seq('AGCCCGTGGGAAAGTCGCCGGGTAATGCACCG')




```python
from Bio.Seq import reverse_complement, transcribe, back_transcribe, translate
```


```python
my_string = "GCTGTTATGGGTCGTTGGAAGGGTGGTCGTGCTGCTGGTTAG"
```


```python
reverse_complement(my_string)
```




    'CTAACCAGCAGCACGACCACCCTTCCAACGACCCATAACAGC'




```python
transcribe(my_string)
```




    'GCUGUUAUGGGUCGUUGGAAGGGUGGUCGUGCUGCUGGUUAG'




```python
back_transcribe(my_string)
```




    'GCTGTTATGGGTCGTTGGAAGGGTGGTCGTGCTGCTGGTTAG'




```python
translate(my_string)
```



    'AVMGRWKGGRAAG*'



## Sequence Annotations (Parts 1-4)
```python
from Bio.SeqRecord import SeqRecord
```


```python
from Bio.SeqRecord import SeqRecord
```


```python
from Bio.Seq import Seq
```


```python
simple_seq = Seq("GATC")
```


```python
simple_seq_r = SeqRecord(simple_seq)
```


```python
simple_seq_r
```




    SeqRecord(seq=Seq('GATC'), id='<unknown id>', name='<unknown name>', description='<unknown description>', dbxrefs=[])




```python
simple_seq_r.id = "AC12345"
```


```python
simple_seq_r.description = "Made up sequence for the VDB Computational Biology Class"
```


```python
print(simple_seq_r.description)
```

    Made up sequence for the VDB Computational Biology Class



```python
simple_seq_r.seq
```




    Seq('GATC')




```python
simple_seq_r
```




    SeqRecord(seq=Seq('GATC'), id='AC12345', name='<unknown name>', description='Made up sequence for the VDB Computational Biology Class', dbxrefs=[])




```python
simple_seq_r.letter_annotations["phred_quality"] = [40, 40, 38, 30]
```


```python
#https://raw.githubusercontent.com/biopython/biopython/master/Tests/GenBank/NC_005816.fna
```


```python
from Bio import SeqIO
```


```python
record = SeqIO.read("NC_005816.fna.txt", "fasta")
```


```python
record
```




    SeqRecord(seq=Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG'), id='gi|45478711|ref|NC_005816.1|', name='gi|45478711|ref|NC_005816.1|', description='gi|45478711|ref|NC_005816.1| Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=[])




```python
record.seq
```




    Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG')




```python
record.id
```




    'gi|45478711|ref|NC_005816.1|'




```python
record.description
```




    'gi|45478711|ref|NC_005816.1| Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence'




```python
record.dbxrefs
```




    []




```python
record.annotations
```




    {}




```python
record.features
```




    []




```python
record = SeqIO.read("NC_005816.gb.txt", "genbank")
```


```python
record
```




    SeqRecord(seq=Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG'), id='NC_005816.1', name='NC_005816', description='Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=['Project:58037'])




```python
record.seq
```




    Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG')




```python
record.id
```




    'NC_005816.1'




```python
record.name
```




    'NC_005816'




```python
record.description
```




    'Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence'




```python
record.letter_annotations
```




    {}




```python
len(record.annotations)
```




    13




```python
record.annotations["source"]
```




    'Yersinia pestis biovar Microtus str. 91001'




```python
record.dbxrefs
```




    ['Project:58037']




```python
len(record.features)
```




    41




```python
from Bio import SeqFeature
```


```python
start_pos = SeqFeature.AfterPosition(5)
```


```python
end_pos = SeqFeature.BetweenPosition(9, left=8, right=9)
```


```python
my_location = SeqFeature.SimpleLocation(start_pos, end_pos)
```


```python
print(my_location)
```

    [>5:(8^9)]



```python
my_location.start
```




    AfterPosition(5)




```python
my_location.end
```




    BetweenPosition(9, left=8, right=9)




```python
int(my_location.end)
```




    9




```python
int(my_location.start)
```




    5




```python
exact_location = SeqFeature.SimpleLocation(5,9)
```


```python
print(exact_location)
```

    [5:9]



```python
exact_location.start
```




    ExactPosition(5)




```python
from Bio.SeqRecord import SeqRecord
```


```python
record = SeqRecord(Seq("MMYQQGCFAGGTVLRLAKDLAENNRGARVLVVCSEITAVTFRGPSETHLDSMVGQALFGD"
                  "GAGAVIVGSDPDLSVERPLYELVWTGATLLPDSEGAIDGHLREVGLTFHLLKDVPGLISK"
                  "NIEKSLKEAFTPLGISDWNSTFWIAHPGGPAILDQVEAKLGLKEEKMRATREVLSEYGNM"
                  "SSAC"),
                        id="gi|14150838|gb|AAK54648.1|AF376133_1",
        description="chalcone synthase [Cucumis sativus]",
                  )
```


```python
print(record.format("fasta"))
```

    >gi|14150838|gb|AAK54648.1|AF376133_1 chalcone synthase [Cucumis sativus]
    MMYQQGCFAGGTVLRLAKDLAENNRGARVLVVCSEITAVTFRGPSETHLDSMVGQALFGD
    GAGAVIVGSDPDLSVERPLYELVWTGATLLPDSEGAIDGHLREVGLTFHLLKDVPGLISK
    NIEKSLKEAFTPLGISDWNSTFWIAHPGGPAILDQVEAKLGLKEEKMRATREVLSEYGNM
    SSAC
    



```python
print(record)
```

    ID: gi|14150838|gb|AAK54648.1|AF376133_1
    Name: <unknown name>
    Description: chalcone synthase [Cucumis sativus]
    Number of features: 0
    Seq('MMYQQGCFAGGTVLRLAKDLAENNRGARVLVVCSEITAVTFRGPSETHLDSMVG...SAC')



```python
from Bio import SeqIO
```


```python
record = SeqIO.read("NC_005816.gb.txt", "genbank")
```


```python
record
```




    SeqRecord(seq=Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG'), id='NC_005816.1', name='NC_005816', description='Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=['Project:58037'])




```python
len(record)
```




    9609




```python
len(record.features)
```




    41




```python
print(record.features[20])
```

    type: gene
    location: [4342:4780](+)
    qualifiers:
        Key: db_xref, Value: ['GeneID:2767712']
        Key: gene, Value: ['pim']
        Key: locus_tag, Value: ['YP_pPCP05']
    



```python
print(record.features[21])
```

    type: CDS
    location: [4342:4780](+)
    qualifiers:
        Key: codon_start, Value: ['1']
        Key: db_xref, Value: ['GI:45478716', 'GeneID:2767712']
        Key: gene, Value: ['pim']
        Key: locus_tag, Value: ['YP_pPCP05']
        Key: note, Value: ['similar to many previously sequenced pesticin immunity protein entries of Yersinia pestis plasmid pPCP, e.g. gi| 16082683|,ref|NP_395230.1| (NC_003132) , gi|1200166|emb|CAA90861.1| (Z54145 ) , gi|1488655| emb|CAA63439.1| (X92856) , gi|2996219|gb|AAC62543.1| (AF053945) , and gi|5763814|emb|CAB531 67.1| (AL109969)']
        Key: product, Value: ['pesticin immunity protein']
        Key: protein_id, Value: ['NP_995571.1']
        Key: transl_table, Value: ['11']
        Key: translation, Value: ['MGGGMISKLFCLALIFLSSSGLAEKNTYTAKDILQNLELNTFGNSLSHGIYGKQTTFKQTEFTNIKSNTKKHIALINKDNSWMISLKILGIKRDEYTVCFEDFSLIRPPTYVAIHPLLIKKVKSGNFIVVKEIKKSIPGCTVYYH']
    



```python
sub_record = record[4300:4800]
```


```python
len(sub_record)
```




    500




```python
len(sub_record.features)
```




    2




```python
sub_record.features[0]
```




    SeqFeature(SimpleLocation(ExactPosition(42), ExactPosition(480), strand=1), type='gene', qualifiers=...)




```python
sub_record.features[1]
```




    SeqFeature(SimpleLocation(ExactPosition(42), ExactPosition(480), strand=1), type='CDS', qualifiers=...)




```python
print(sub_record.features[0])
```

    type: gene
    location: [42:480](+)
    qualifiers:
        Key: db_xref, Value: ['GeneID:2767712']
        Key: gene, Value: ['pim']
        Key: locus_tag, Value: ['YP_pPCP05']
    



```python
print(sub_record.features[1])
```

    type: CDS
    location: [42:480](+)
    qualifiers:
        Key: codon_start, Value: ['1']
        Key: db_xref, Value: ['GI:45478716', 'GeneID:2767712']
        Key: gene, Value: ['pim']
        Key: locus_tag, Value: ['YP_pPCP05']
        Key: note, Value: ['similar to many previously sequenced pesticin immunity protein entries of Yersinia pestis plasmid pPCP, e.g. gi| 16082683|,ref|NP_395230.1| (NC_003132) , gi|1200166|emb|CAA90861.1| (Z54145 ) , gi|1488655| emb|CAA63439.1| (X92856) , gi|2996219|gb|AAC62543.1| (AF053945) , and gi|5763814|emb|CAB531 67.1| (AL109969)']
        Key: product, Value: ['pesticin immunity protein']
        Key: protein_id, Value: ['NP_995571.1']
        Key: transl_table, Value: ['11']
        Key: translation, Value: ['MGGGMISKLFCLALIFLSSSGLAEKNTYTAKDILQNLELNTFGNSLSHGIYGKQTTFKQTEFTNIKSNTKKHIALINKDNSWMISLKILGIKRDEYTVCFEDFSLIRPPTYVAIHPLLIKKVKSGNFIVVKEIKKSIPGCTVYYH']
    



```python
sub_record.annotations
```




    {'molecule_type': 'DNA'}




```python
sub_record.dbxrefs
```




    []




```python
sub_record.annotations["topology"] = "linear"
```


```python
sub_record.annotations
```




    {'molecule_type': 'DNA', 'topology': 'linear'}




```python
sub_record.id
```




    'NC_005816.1'




```python
sub_record.name
```




    'NC_005816'




```python
sub_record.description
```




    'Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence'




```python
sub_record.description = 'Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence'
```


```python
sub_record.description
```




    'Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence'




```python
print(sub_record.format("genbank")[:200] + "...")
```

    LOCUS       NC_005816                500 bp    DNA     linear   UNK 01-JAN-1980
    DEFINITION  Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete
                sequence.
    ACCESSION   NC_0058...



```python
record = SeqIO.read("NC_005816.gb.txt", "genbank")
```


```python
record
```




    SeqRecord(seq=Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG'), id='NC_005816.1', name='NC_005816', description='Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=['Project:58037'])




```python
len(record)
```




    9609




```python
len(record.features)
```




    41




```python
record.dbxrefs
```




    ['Project:58037']




```python
record.annotations.keys()
```




    dict_keys(['molecule_type', 'topology', 'data_file_division', 'date', 'accessions', 'sequence_version', 'gi', 'keywords', 'source', 'organism', 'taxonomy', 'references', 'comment'])




```python
shifted = record[2000:] + record[:2000]
```


```python
shifted
```




    SeqRecord(seq=Seq('GATACGCAGTCATATTTTTTACACAATTCTCTAATCCCGACAAGGTCGTAGGTC...GGA'), id='NC_005816.1', name='NC_005816', description='Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=[])




```python
len(shifted)
```




    9609




```python
len(shifted.features)
```




    40




```python
shifted.annotations.keys()
```




    dict_keys(['molecule_type'])




```python
shifted.dbxrefs
```




    []




```python
shifted.dbxrefs = record.dbxrefs[:]
```


```python
shifted.dbxrefs
```




    ['Project:58037']




```python
shifted.annotations = record.annotations.copy()
```


```python
shifted.annotations.keys()
```




    dict_keys(['molecule_type', 'topology', 'data_file_division', 'date', 'accessions', 'sequence_version', 'gi', 'keywords', 'source', 'organism', 'taxonomy', 'references', 'comment'])




```python
record
```




    SeqRecord(seq=Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG'), id='NC_005816.1', name='NC_005816', description='Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=['Project:58037'])




```python
print("%s %i %i %i %i" % (record.id, len(record), len(record.features), len(record.dbxrefs), len(record.annotations)))
```

    NC_005816.1 9609 41 1 13



```python
rc = record.reverse_complement(id ="Testing")
```


```python
rc
```




    SeqRecord(seq=Seq('CAGGGGTCGGGGTACGCATTCCCTCATGCGTCAATATTATCTGGCATTGCGATG...ACA'), id='Testing', name='<unknown name>', description='<unknown description>', dbxrefs=[])




```python
print("%s %i %i %i %i" % (rc.id, len(rc), len(rc.features), len(rc.dbxrefs), len(rc.annotations)))
```

    Testing 9609 41 0 0


## Sequence I/O (Parts 1-3)
```python
from Bio import SeqIO
```


```python
for seq_record in SeqIO.parse("ls_orchid.fasta.txt", "fasta"):
    print(seq_record.id)
    print(repr(seq_record.seq))
    print(len(seq_record))
```

    gi|2765658|emb|Z78533.1|CIZ78533
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGATGAGACCGTGG...CGC')
    740
    gi|2765657|emb|Z78532.1|CCZ78532
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAACAG...GGC')
    753
    gi|2765656|emb|Z78531.1|CFZ78531
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGCAG...TAA')
    748
    gi|2765655|emb|Z78530.1|CMZ78530
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAAACAACAT...CAT')
    744
    gi|2765654|emb|Z78529.1|CLZ78529
    Seq('ACGGCGAGCTGCCGAAGGACATTGTTGAGACAGCAGAATATACGATTGAGTGAA...AAA')
    733
    gi|2765652|emb|Z78527.1|CYZ78527
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAG...CCC')
    718
    gi|2765651|emb|Z78526.1|CGZ78526
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAG...TGT')
    730
    gi|2765650|emb|Z78525.1|CAZ78525
    Seq('TGTTGAGATAGCAGAATATACATCGAGTGAATCCGGAGGACCTGTGGTTATTCG...GCA')
    704
    gi|2765649|emb|Z78524.1|CFZ78524
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATAGTAG...AGC')
    740
    gi|2765648|emb|Z78523.1|CHZ78523
    Seq('CGTAACCAGGTTTCCGTAGGTGAACCTGCGGCAGGATCATTGTTGAGACAGCAG...AAG')
    709
    gi|2765647|emb|Z78522.1|CMZ78522
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGCAG...GAG')
    700
    gi|2765646|emb|Z78521.1|CCZ78521
    Seq('GTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAGAATATATGATCGAGT...ACC')
    726
    gi|2765645|emb|Z78520.1|CSZ78520
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGCAG...TTT')
    753
    gi|2765644|emb|Z78519.1|CPZ78519
    Seq('ATATGATCGAGTGAATCTGGTGGACTTGTGGTTACTCAGCTCGCCATAGGCTTT...TTA')
    699
    gi|2765643|emb|Z78518.1|CRZ78518
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGGAGGATCATTGTTGAGATAGTAG...TCC')
    658
    gi|2765642|emb|Z78517.1|CFZ78517
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAG...AGC')
    752
    gi|2765641|emb|Z78516.1|CPZ78516
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACAGTAT...TAA')
    726
    gi|2765640|emb|Z78515.1|MXZ78515
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGCTGAGACCGTAG...AGC')
    765
    gi|2765639|emb|Z78514.1|PSZ78514
    Seq('CGTAACAAGGTTTCCGTAGGTGGACCTTCGGGAGGATCATTTTTGAAGCCCCCA...CTA')
    755
    gi|2765638|emb|Z78513.1|PBZ78513
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...GAG')
    742
    gi|2765637|emb|Z78512.1|PWZ78512
    Seq('CGTAACAAGGTTTCCGTAGGTGGACCTTCGGGAGGATCATTTTTGAAGCCCCCA...AGC')
    762
    gi|2765636|emb|Z78511.1|PEZ78511
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTTCGGAAGGATCATTGTTGAGACCCCCA...GGA')
    745
    gi|2765635|emb|Z78510.1|PCZ78510
    Seq('CTAACCAGGGTTCCGAGGTGACCTTCGGGAGGATTCCTTTTTAAGCCCCCGAAA...TTA')
    750
    gi|2765634|emb|Z78509.1|PPZ78509
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...GGA')
    731
    gi|2765633|emb|Z78508.1|PLZ78508
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...TGA')
    741
    gi|2765632|emb|Z78507.1|PLZ78507
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCCCCA...TGA')
    740
    gi|2765631|emb|Z78506.1|PLZ78506
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCAA...TGA')
    727
    gi|2765630|emb|Z78505.1|PSZ78505
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCCA...TTT')
    711
    gi|2765629|emb|Z78504.1|PKZ78504
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTTCGGAAGGATCATTGTTGAGACCGCAA...TAA')
    743
    gi|2765628|emb|Z78503.1|PCZ78503
    Seq('CGTAACCAGGTTTCCGTAGGTGAACCTCCGGAAGGATCCTTGTTGAGACCGCCA...TAA')
    727
    gi|2765627|emb|Z78502.1|PBZ78502
    Seq('CGTAACCAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGACCGCCA...CGC')
    757
    gi|2765626|emb|Z78501.1|PCZ78501
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGACCGCAA...AGA')
    770
    gi|2765625|emb|Z78500.1|PWZ78500
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGCTCATTGTTGAGACCGCAA...AAG')
    767
    gi|2765624|emb|Z78499.1|PMZ78499
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAGGGATCATTGTTGAGATCGCAT...ACC')
    759
    gi|2765623|emb|Z78498.1|PMZ78498
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAAGGTCATTGTTGAGATCACAT...AGC')
    750
    gi|2765622|emb|Z78497.1|PDZ78497
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')
    788
    gi|2765621|emb|Z78496.1|PAZ78496
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...AGC')
    774
    gi|2765620|emb|Z78495.1|PEZ78495
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...GTG')
    789
    gi|2765619|emb|Z78494.1|PNZ78494
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGGTCGCAT...AAG')
    688
    gi|2765618|emb|Z78493.1|PGZ78493
    Seq('CGTAACAAGGATTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...CCC')
    719
    gi|2765617|emb|Z78492.1|PBZ78492
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...ATA')
    743
    gi|2765616|emb|Z78491.1|PCZ78491
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...AGC')
    737
    gi|2765615|emb|Z78490.1|PFZ78490
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')
    728
    gi|2765614|emb|Z78489.1|PDZ78489
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GGC')
    740
    gi|2765613|emb|Z78488.1|PTZ78488
    Seq('CTGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACGCAATAATTGATCGA...GCT')
    696
    gi|2765612|emb|Z78487.1|PHZ78487
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TAA')
    732
    gi|2765611|emb|Z78486.1|PBZ78486
    Seq('CGTCACGAGGTTTCCGTAGGTGAATCTGCGGGAGGATCATTGTTGAGATCACAT...TGA')
    731
    gi|2765610|emb|Z78485.1|PHZ78485
    Seq('CTGAACCTGGTGTCCGAAGGTGAATCTGCGGATGGATCATTGTTGAGATATCAT...GTA')
    735
    gi|2765609|emb|Z78484.1|PCZ78484
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGGGGAAGGATCATTGTTGAGATCACAT...TTT')
    720
    gi|2765608|emb|Z78483.1|PVZ78483
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GCA')
    740
    gi|2765607|emb|Z78482.1|PEZ78482
    Seq('TCTACTGCAGTGACCGAGATTTGCCATCGAGCCTCCTGGGAGCTTTCTTGCTGG...GCA')
    629
    gi|2765606|emb|Z78481.1|PIZ78481
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')
    572
    gi|2765605|emb|Z78480.1|PGZ78480
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')
    587
    gi|2765604|emb|Z78479.1|PPZ78479
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGT')
    700
    gi|2765603|emb|Z78478.1|PVZ78478
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCAGTGTTGAGATCACAT...GGC')
    636
    gi|2765602|emb|Z78477.1|PVZ78477
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGC')
    716
    gi|2765601|emb|Z78476.1|PGZ78476
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...CCC')
    592
    gi|2765600|emb|Z78475.1|PSZ78475
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GGT')
    716
    gi|2765599|emb|Z78474.1|PKZ78474
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACGT...CTT')
    733
    gi|2765598|emb|Z78473.1|PSZ78473
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGG')
    626
    gi|2765597|emb|Z78472.1|PLZ78472
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')
    737
    gi|2765596|emb|Z78471.1|PDZ78471
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')
    740
    gi|2765595|emb|Z78470.1|PPZ78470
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GTT')
    574
    gi|2765594|emb|Z78469.1|PHZ78469
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GTT')
    594
    gi|2765593|emb|Z78468.1|PAZ78468
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCGCAT...GTT')
    610
    gi|2765592|emb|Z78467.1|PSZ78467
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGA')
    730
    gi|2765591|emb|Z78466.1|PPZ78466
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...CCC')
    641
    gi|2765590|emb|Z78465.1|PRZ78465
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGC')
    702
    gi|2765589|emb|Z78464.1|PGZ78464
    Seq('CGTAACAAGGTTTCCGTAGGTGAGCGGAAGGGTCATTGTTGAGATCACATAATA...AGC')
    733
    gi|2765588|emb|Z78463.1|PGZ78463
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGTTCATTGTTGAGATCACAT...AGC')
    738
    gi|2765587|emb|Z78462.1|PSZ78462
    Seq('CGTCACGAGGTCTCCGGATGTGACCCTGCGGAAGGATCATTGTTGAGATCACAT...CAT')
    736
    gi|2765586|emb|Z78461.1|PWZ78461
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...TAA')
    732
    gi|2765585|emb|Z78460.1|PCZ78460
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...TTA')
    745
    gi|2765584|emb|Z78459.1|PDZ78459
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TTT')
    744
    gi|2765583|emb|Z78458.1|PHZ78458
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TTG')
    738
    gi|2765582|emb|Z78457.1|PCZ78457
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...GAG')
    739
    gi|2765581|emb|Z78456.1|PTZ78456
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGC')
    740
    gi|2765580|emb|Z78455.1|PJZ78455
    Seq('CGTAACCAGGTTTCCGTAGGTGGACCTTCGGGAGGATCATTTTTGAGATCACAT...GCA')
    745
    gi|2765579|emb|Z78454.1|PFZ78454
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AAC')
    695
    gi|2765578|emb|Z78453.1|PSZ78453
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GCA')
    745
    gi|2765577|emb|Z78452.1|PBZ78452
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...GCA')
    743
    gi|2765576|emb|Z78451.1|PHZ78451
    Seq('CGTAACAAGGTTTCCGTAGGTGTACCTCCGGAAGGATCATTGTTGAGATCACAT...AGC')
    730
    gi|2765575|emb|Z78450.1|PPZ78450
    Seq('GGAAGGATCATTGCTGATATCACATAATAATTGATCGAGTTAAGCTGGAGGATC...GAG')
    706
    gi|2765574|emb|Z78449.1|PMZ78449
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGC')
    744
    gi|2765573|emb|Z78448.1|PAZ78448
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGG')
    742
    gi|2765572|emb|Z78447.1|PVZ78447
    Seq('CGTAACAAGGATTCCGTAGGTGAACCTGCGGGAGGATCATTGTTGAGATCACAT...AGC')
    694
    gi|2765571|emb|Z78446.1|PAZ78446
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTCCGGAAGGATCATTGTTGAGATCACAT...CCC')
    712
    gi|2765570|emb|Z78445.1|PUZ78445
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...TGT')
    715
    gi|2765569|emb|Z78444.1|PAZ78444
    Seq('CGTAACAAGGTTTCCGTAGGGTGAACTGCGGAAGGATCATTGTTGAGATCACAT...ATT')
    688
    gi|2765568|emb|Z78443.1|PLZ78443
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACAT...AGG')
    784
    gi|2765567|emb|Z78442.1|PBZ78442
    Seq('GTAGGTGAACCTGCGGAAGGATCATTGTTGAGATCACATAATAATTGATCGAGT...AGT')
    721
    gi|2765566|emb|Z78441.1|PSZ78441
    Seq('GGAAGGTCATTGCCGATATCACATAATAATTGATCGAGTTAATCTGGAGGATCT...GAG')
    703
    gi|2765565|emb|Z78440.1|PPZ78440
    Seq('CGTAACAAGGTTTCCGTAGGTGGACCTCCGGGAGGATCATTGTTGAGATCACAT...GCA')
    744
    gi|2765564|emb|Z78439.1|PBZ78439
    Seq('CATTGTTGAGATCACATAATAATTGATCGAGTTAATCTGGAGGATCTGTTTACT...GCC')
    592



```python
identifiers = [seq_record.id for seq_record in SeqIO.parse("ls_orchid.gbk.txt", "genbank")]
```


```python
identifiers
```




    ['Z78533.1',
     'Z78532.1',
     'Z78531.1',
     'Z78530.1',
     'Z78529.1',
     'Z78527.1',
     'Z78526.1',
     'Z78525.1',
     'Z78524.1',
     'Z78523.1',
     'Z78522.1',
     'Z78521.1',
     'Z78520.1',
     'Z78519.1',
     'Z78518.1',
     'Z78517.1',
     'Z78516.1',
     'Z78515.1',
     'Z78514.1',
     'Z78513.1',
     'Z78512.1',
     'Z78511.1',
     'Z78510.1',
     'Z78509.1',
     'Z78508.1',
     'Z78507.1',
     'Z78506.1',
     'Z78505.1',
     'Z78504.1',
     'Z78503.1',
     'Z78502.1',
     'Z78501.1',
     'Z78500.1',
     'Z78499.1',
     'Z78498.1',
     'Z78497.1',
     'Z78496.1',
     'Z78495.1',
     'Z78494.1',
     'Z78493.1',
     'Z78492.1',
     'Z78491.1',
     'Z78490.1',
     'Z78489.1',
     'Z78488.1',
     'Z78487.1',
     'Z78486.1',
     'Z78485.1',
     'Z78484.1',
     'Z78483.1',
     'Z78482.1',
     'Z78481.1',
     'Z78480.1',
     'Z78479.1',
     'Z78478.1',
     'Z78477.1',
     'Z78476.1',
     'Z78475.1',
     'Z78474.1',
     'Z78473.1',
     'Z78472.1',
     'Z78471.1',
     'Z78470.1',
     'Z78469.1',
     'Z78468.1',
     'Z78467.1',
     'Z78466.1',
     'Z78465.1',
     'Z78464.1',
     'Z78463.1',
     'Z78462.1',
     'Z78461.1',
     'Z78460.1',
     'Z78459.1',
     'Z78458.1',
     'Z78457.1',
     'Z78456.1',
     'Z78455.1',
     'Z78454.1',
     'Z78453.1',
     'Z78452.1',
     'Z78451.1',
     'Z78450.1',
     'Z78449.1',
     'Z78448.1',
     'Z78447.1',
     'Z78446.1',
     'Z78445.1',
     'Z78444.1',
     'Z78443.1',
     'Z78442.1',
     'Z78441.1',
     'Z78440.1',
     'Z78439.1']




```python
record_iterator = SeqIO.parse("ls_orchid.fasta.txt", "fasta")
```


```python
first_record = next(record_iterator)
```


```python
print(first_record.id)
```

    gi|2765658|emb|Z78533.1|CIZ78533



```python
print(first_record.description)
```

    gi|2765658|emb|Z78533.1|CIZ78533 C.irapeanum 5.8S rRNA gene and ITS1 and ITS2 DNA



```python
second_record = next(record_iterator)
```


```python
print(second_record.id)
```

    gi|2765657|emb|Z78532.1|CCZ78532



```python
print(second_record.description)
```

    gi|2765657|emb|Z78532.1|CCZ78532 C.californicum 5.8S rRNA gene and ITS1 and ITS2 DNA



```python
records = list(SeqIO.parse("ls_orchid.gbk.txt", "genbank"))
```


```python
print("Found %i records" % len(records))
```

    Found 94 records



```python
print("The last record")
last_record = records[-1]
print(last_record.id)
print(repr(last_record.seq))
print(len(last_record))
```

    The last record
    Z78439.1
    Seq('CATTGTTGAGATCACATAATAATTGATCGAGTTAATCTGGAGGATCTGTTTACT...GCC')
    592



```python
print("The first record")
first_record = records[-3]
print(first_record.id)
print(repr(first_record.seq))
print(len(first_record))
```

    The first record
    Z78441.1
    Seq('GGAAGGTCATTGCCGATATCACATAATAATTGATCGAGTTAATCTGGAGGATCT...GAG')
    703



```python
record_iterator = SeqIO.parse("ls_orchid.gbk.txt", "genbank")
```


```python
first_record = next(record_iterator)
```


```python
print(first_record)
```

    ID: Z78533.1
    Name: Z78533
    Description: C.irapeanum 5.8S rRNA gene and ITS1 and ITS2 DNA
    Number of features: 5
    /molecule_type=DNA
    /topology=linear
    /data_file_division=PLN
    /date=30-NOV-2006
    /accessions=['Z78533']
    /sequence_version=1
    /gi=2765658
    /keywords=['5.8S ribosomal RNA', '5.8S rRNA gene', 'internal transcribed spacer', 'ITS1', 'ITS2']
    /source=Cypripedium irapeanum
    /organism=Cypripedium irapeanum
    /taxonomy=['Eukaryota', 'Viridiplantae', 'Streptophyta', 'Embryophyta', 'Tracheophyta', 'Spermatophyta', 'Magnoliophyta', 'Liliopsida', 'Asparagales', 'Orchidaceae', 'Cypripedioideae', 'Cypripedium']
    /references=[Reference(title='Phylogenetics of the slipper orchids (Cypripedioideae: Orchidaceae): nuclear rDNA ITS sequences', ...), Reference(title='Direct Submission', ...)]
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGATGAGACCGTGG...CGC')



```python
print(first_record.annotations)
```

    {'molecule_type': 'DNA', 'topology': 'linear', 'data_file_division': 'PLN', 'date': '30-NOV-2006', 'accessions': ['Z78533'], 'sequence_version': 1, 'gi': '2765658', 'keywords': ['5.8S ribosomal RNA', '5.8S rRNA gene', 'internal transcribed spacer', 'ITS1', 'ITS2'], 'source': 'Cypripedium irapeanum', 'organism': 'Cypripedium irapeanum', 'taxonomy': ['Eukaryota', 'Viridiplantae', 'Streptophyta', 'Embryophyta', 'Tracheophyta', 'Spermatophyta', 'Magnoliophyta', 'Liliopsida', 'Asparagales', 'Orchidaceae', 'Cypripedioideae', 'Cypripedium'], 'references': [Reference(title='Phylogenetics of the slipper orchids (Cypripedioideae: Orchidaceae): nuclear rDNA ITS sequences', ...), Reference(title='Direct Submission', ...)]}



```python
print(first_record.annotations.keys())
```

    dict_keys(['molecule_type', 'topology', 'data_file_division', 'date', 'accessions', 'sequence_version', 'gi', 'keywords', 'source', 'organism', 'taxonomy', 'references'])



```python
print(first_record.annotations.values())
```

    dict_values(['DNA', 'linear', 'PLN', '30-NOV-2006', ['Z78533'], 1, '2765658', ['5.8S ribosomal RNA', '5.8S rRNA gene', 'internal transcribed spacer', 'ITS1', 'ITS2'], 'Cypripedium irapeanum', 'Cypripedium irapeanum', ['Eukaryota', 'Viridiplantae', 'Streptophyta', 'Embryophyta', 'Tracheophyta', 'Spermatophyta', 'Magnoliophyta', 'Liliopsida', 'Asparagales', 'Orchidaceae', 'Cypripedioideae', 'Cypripedium'], [Reference(title='Phylogenetics of the slipper orchids (Cypripedioideae: Orchidaceae): nuclear rDNA ITS sequences', ...), Reference(title='Direct Submission', ...)]])



```python
print(first_record.annotations["source"])
```

    Cypripedium irapeanum



```python
print(first_record.annotations["organism"])
```

    Cypripedium irapeanum



```python
all_species = []
```


```python
for seq_record in SeqIO.parse("ls_orchid.gbk.txt", "genbank"):
    all_species.append(seq_record.annotations["organism"])
```


```python
print(all_species)
```

    ['Cypripedium irapeanum', 'Cypripedium californicum', 'Cypripedium fasciculatum', 'Cypripedium margaritaceum', 'Cypripedium lichiangense', 'Cypripedium yatabeanum', 'Cypripedium guttatum', 'Cypripedium acaule', 'Cypripedium formosanum', 'Cypripedium himalaicum', 'Cypripedium macranthon', 'Cypripedium calceolus', 'Cypripedium segawai', 'Cypripedium parviflorum var. pubescens', 'Cypripedium reginae', 'Cypripedium flavum', 'Cypripedium passerinum', 'Mexipedium xerophyticum', 'Phragmipedium schlimii', 'Phragmipedium besseae', 'Phragmipedium wallisii', 'Phragmipedium exstaminodium', 'Phragmipedium caricinum', 'Phragmipedium pearcei', 'Phragmipedium longifolium', 'Phragmipedium lindenii', 'Phragmipedium lindleyanum', 'Phragmipedium sargentianum', 'Phragmipedium kaiteurum', 'Phragmipedium czerwiakowianum', 'Phragmipedium boissierianum', 'Phragmipedium caudatum', 'Phragmipedium warszewiczianum', 'Paphiopedilum micranthum', 'Paphiopedilum malipoense', 'Paphiopedilum delenatii', 'Paphiopedilum armeniacum', 'Paphiopedilum emersonii', 'Paphiopedilum niveum', 'Paphiopedilum godefroyae', 'Paphiopedilum bellatulum', 'Paphiopedilum concolor', 'Paphiopedilum fairrieanum', 'Paphiopedilum druryi', 'Paphiopedilum tigrinum', 'Paphiopedilum hirsutissimum', 'Paphiopedilum barbigerum', 'Paphiopedilum henryanum', 'Paphiopedilum charlesworthii', 'Paphiopedilum villosum', 'Paphiopedilum exul', 'Paphiopedilum insigne', 'Paphiopedilum gratrixianum', 'Paphiopedilum primulinum', 'Paphiopedilum victoria', 'Paphiopedilum victoria', 'Paphiopedilum glaucophyllum', 'Paphiopedilum supardii', 'Paphiopedilum kolopakingii', 'Paphiopedilum sanderianum', 'Paphiopedilum lowii', 'Paphiopedilum dianthum', 'Paphiopedilum parishii', 'Paphiopedilum haynaldianum', 'Paphiopedilum adductum', 'Paphiopedilum stonei', 'Paphiopedilum philippinense', 'Paphiopedilum rothschildianum', 'Paphiopedilum glanduliferum', 'Paphiopedilum glanduliferum', 'Paphiopedilum sukhakulii', 'Paphiopedilum wardii', 'Paphiopedilum ciliolare', 'Paphiopedilum dayanum', 'Paphiopedilum hennisianum', 'Paphiopedilum callosum', 'Paphiopedilum tonsum', 'Paphiopedilum javanicum', 'Paphiopedilum fowliei', 'Paphiopedilum schoseri', 'Paphiopedilum bougainvilleanum', 'Paphiopedilum hookerae', 'Paphiopedilum papuanum', 'Paphiopedilum mastersianum', 'Paphiopedilum argus', 'Paphiopedilum venustum', 'Paphiopedilum acmodontum', 'Paphiopedilum urbanianum', 'Paphiopedilum appletonianum', 'Paphiopedilum lawrenceanum', 'Paphiopedilum bullenianum', 'Paphiopedilum superbiens', 'Paphiopedilum purpuratum', 'Paphiopedilum barbatum']



```python
all_species = [
    seq_record.annotations["organism"]
    for seq_record in SeqIO.parse("ls_orchid.gbk.txt", "genbank")
]
```


```python
print(all_species)
```

    ['Cypripedium irapeanum', 'Cypripedium californicum', 'Cypripedium fasciculatum', 'Cypripedium margaritaceum', 'Cypripedium lichiangense', 'Cypripedium yatabeanum', 'Cypripedium guttatum', 'Cypripedium acaule', 'Cypripedium formosanum', 'Cypripedium himalaicum', 'Cypripedium macranthon', 'Cypripedium calceolus', 'Cypripedium segawai', 'Cypripedium parviflorum var. pubescens', 'Cypripedium reginae', 'Cypripedium flavum', 'Cypripedium passerinum', 'Mexipedium xerophyticum', 'Phragmipedium schlimii', 'Phragmipedium besseae', 'Phragmipedium wallisii', 'Phragmipedium exstaminodium', 'Phragmipedium caricinum', 'Phragmipedium pearcei', 'Phragmipedium longifolium', 'Phragmipedium lindenii', 'Phragmipedium lindleyanum', 'Phragmipedium sargentianum', 'Phragmipedium kaiteurum', 'Phragmipedium czerwiakowianum', 'Phragmipedium boissierianum', 'Phragmipedium caudatum', 'Phragmipedium warszewiczianum', 'Paphiopedilum micranthum', 'Paphiopedilum malipoense', 'Paphiopedilum delenatii', 'Paphiopedilum armeniacum', 'Paphiopedilum emersonii', 'Paphiopedilum niveum', 'Paphiopedilum godefroyae', 'Paphiopedilum bellatulum', 'Paphiopedilum concolor', 'Paphiopedilum fairrieanum', 'Paphiopedilum druryi', 'Paphiopedilum tigrinum', 'Paphiopedilum hirsutissimum', 'Paphiopedilum barbigerum', 'Paphiopedilum henryanum', 'Paphiopedilum charlesworthii', 'Paphiopedilum villosum', 'Paphiopedilum exul', 'Paphiopedilum insigne', 'Paphiopedilum gratrixianum', 'Paphiopedilum primulinum', 'Paphiopedilum victoria', 'Paphiopedilum victoria', 'Paphiopedilum glaucophyllum', 'Paphiopedilum supardii', 'Paphiopedilum kolopakingii', 'Paphiopedilum sanderianum', 'Paphiopedilum lowii', 'Paphiopedilum dianthum', 'Paphiopedilum parishii', 'Paphiopedilum haynaldianum', 'Paphiopedilum adductum', 'Paphiopedilum stonei', 'Paphiopedilum philippinense', 'Paphiopedilum rothschildianum', 'Paphiopedilum glanduliferum', 'Paphiopedilum glanduliferum', 'Paphiopedilum sukhakulii', 'Paphiopedilum wardii', 'Paphiopedilum ciliolare', 'Paphiopedilum dayanum', 'Paphiopedilum hennisianum', 'Paphiopedilum callosum', 'Paphiopedilum tonsum', 'Paphiopedilum javanicum', 'Paphiopedilum fowliei', 'Paphiopedilum schoseri', 'Paphiopedilum bougainvilleanum', 'Paphiopedilum hookerae', 'Paphiopedilum papuanum', 'Paphiopedilum mastersianum', 'Paphiopedilum argus', 'Paphiopedilum venustum', 'Paphiopedilum acmodontum', 'Paphiopedilum urbanianum', 'Paphiopedilum appletonianum', 'Paphiopedilum lawrenceanum', 'Paphiopedilum bullenianum', 'Paphiopedilum superbiens', 'Paphiopedilum purpuratum', 'Paphiopedilum barbatum']



```python
all_species = []
```


```python
for seq_record in SeqIO.parse("ls_orchid.fasta.txt", "fasta"):
    all_species.append(seq_record.description.split()[1])
```


```python
print(all_species)
```

    ['C.irapeanum', 'C.californicum', 'C.fasciculatum', 'C.margaritaceum', 'C.lichiangense', 'C.yatabeanum', 'C.guttatum', 'C.acaule', 'C.formosanum', 'C.himalaicum', 'C.macranthum', 'C.calceolus', 'C.segawai', 'C.pubescens', 'C.reginae', 'C.flavum', 'C.passerinum', 'M.xerophyticum', 'P.schlimii', 'P.besseae', 'P.wallisii', 'P.exstaminodium', 'P.caricinum', 'P.pearcei', 'P.longifolium', 'P.lindenii', 'P.lindleyanum', 'P.sargentianum', 'P.kaiteurum', 'P.czerwiakowianum', 'P.boissierianum', 'P.caudatum', 'P.warszewiczianum', 'P.micranthum', 'P.malipoense', 'P.delenatii', 'P.armeniacum', 'P.emersonii', 'P.niveum', 'P.godefroyae', 'P.bellatulum', 'P.concolor', 'P.fairrieanum', 'P.druryi', 'P.tigrinum', 'P.hirsutissimum', 'P.barbigerum', 'P.henryanum', 'P.charlesworthii', 'P.villosum', 'P.exul', 'P.insigne', 'P.gratrixianum', 'P.primulinum', 'P.victoria', 'P.victoria', 'P.glaucophyllum', 'P.supardii', 'P.kolopakingii', 'P.sanderianum', 'P.lowii', 'P.dianthum', 'P.parishii', 'P.haynaldianum', 'P.adductum', 'P.stonei', 'P.philippinense', 'P.rothschildianum', 'P.glanduliferum', 'P.glanduliferum', 'P.sukhakulii', 'P.wardii', 'P.ciliolare', 'P.dayanum', 'P.hennisianum', 'P.callosum', 'P.tonsum', 'P.javanicum', 'P.fowliei', 'P.schoseri', 'P.bougainvilleanum', 'P.hookerae', 'P.papuanum', 'P.mastersianum', 'P.argus', 'P.venustum', 'P.acmodontum', 'P.urbanianum', 'P.appletonianum', 'P.lawrenceanum', 'P.bullenianum', 'P.superbiens', 'P.purpuratum', 'P.barbatum']



```python
record_iterator = SeqIO.parse("ls_orchid.fasta.txt", "fasta")
```


```python
first_record = next(record_iterator)
```


```python
first_record.id
```




    'gi|2765658|emb|Z78533.1|CIZ78533'




```python
first_record.id = "new_id"
```


```python
first_record.id
```




    'new_id'




```python
first_record
```




    SeqRecord(seq=Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGATGAGACCGTGG...CGC'), id='new_id', name='gi|2765658|emb|Z78533.1|CIZ78533', description='gi|2765658|emb|Z78533.1|CIZ78533 C.irapeanum 5.8S rRNA gene and ITS1 and ITS2 DNA', dbxrefs=[])




```python
first_record.description = first_record.id + " " + "desired new description"
```


```python
print(first_record.format("fasta"[:200]))
```

    >new_id desired new description
    CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGATGAGACCGTGGAATAAA
    CGATCGAGTGAATCCGGAGGACCGGTGTACTCAGCTCACCGGGGGCATTGCTCCCGTGGT
    GACCCTGATTTGTTGTTGGGCCGCCTCGGGAGCGTCCATGGCGGGTTTGAACCTCTAGCC
    CGGCGCAGTTTGGGCGCCAAGCCATATGAAAGCATCACCGGCGAATGGCATTGTCTTCCC
    CAAAACCCGGAGCGGCGGCGTGCTGTCGCGTGCCCAATGAATTTTGATGACTCTCGCAAA
    CGGGAATCTTGGCTCTTTGCATCGGATGGAAGGACGCAGCGAAATGCGATAAGTGGTGTG
    AATTGCAAGATCCCGTGAACCATCGAGTCTTTTGAACGCAAGTTGCGCCCGAGGCCATCA
    GGCTAAGGGCACGCCTGCTTGGGCGTCGCGCTTCGTCTCTCTCCTGCCAATGCTTGCCCG
    GCATACAGCCAGGCCGGCGTGGTGCGGATGTGAAAGATTGGCCCCTTGTGCCTAGGTGCG
    GCGGGTCCAAGAGCTGGTGTTTTGATGGCCCGGAACCCGGCAAGAGGTGGACGGATGCTG
    GCAGCAGCTGCCGTGCGAATCCCCCATGTTGTCGTGCTTGTCGGACAGGCAGGAGAACCC
    TTCCGAACCCCAATGGAGGGCGGTTGACCGCCATTCGGATGTGACCCCAGGTCAGGCGGG
    GGCACCCGCTGAGTTTACGC
    



```python
from Bio.Seq import Seq
```


```python
from Bio.SeqRecord import SeqRecord
```


```python
rec1 = SeqRecord(
    Seq("MMYQQGCFAGGTVLRLAKDLAENNRGARVLVVCSEITAVTFRGPSETHLDSMVGQALFGD"
    "GAGAVIVGSDPDLSVERPLYELVWTGATLLPDSEGAIDGHLREVGLTFHLLKDVPGLISK"
    "NIEKSLKEAFTPLGISDWNSTFWIAHPGGPAILDQVEAKLGLKEEKMRATREVLSEYGNM"
    "SSAC",
),
    id = "gi|14150838|gb|AAK54658.1|AF376133_1",
    description = "chalcone synthase [Cucumis sativus]",
)
```


```python
rec2 = SeqRecord(
    Seq("YPDYYFRITNREHKAELKEKFQRMCDKSMIKKRYMYLTEEILKENPSMCEYMAPSLDARQDMVVEIPKLGKEAAVKAIKEWGQ"),
    id="gi|13919613|gb|AAK33142.1|",
    description="chalcone synthase [Fragaria vesca subsp. bracteata]"
)
```


```python
rec3 = SeqRecord(
    Seq("MVTVEEFRRAQCAEGPATVMAIGTATPSNCVDQSTYPDYYFRITNSEHKVELKEKFKRMC"
        "EKSMIKKRYMHLTEEILKENPNICAYMAPSLDARQDIVVVEVPKLGKEAAQKAIKEWGQP"
        "KSKITHLVFCTTSGVDMPGCDYQLTKLLGLRPSVKRFMMYQQGCFAGGTVLRMAKDLAEN"
        "NKGARVLVVCSEITAVTFRGPNDTHLDSLVGQALFGDGAAAVIIGSDPIPEVERPLFELV"
        "SAAQTLLPDSEGAIDGHLREVGLTFHLLKDVPGLISKNIEKSLVEAFQPLGISDWNSLFW"
        "IAHPGGPAILDQVELKLGLKQEKLKATRKVLSNYGNMSSACVLFILDEMRKASAKEGLGT"
        "TGEGLEWGVLFGFGPGLTVETVVLHSVAT"),
    id="gi|13925890|gb|AAK49457.1|",
    description="chalcone synthase [Nicotiana tabacum]"
)
```


```python
my_records = [rec1, rec2, rec3]
```


```python
from Bio import SeqIO
SeqIO.write(my_records, "my_example.faa", "fasta")
```




    3




```python
records = SeqIO.parse("ls_orchid.gbk.txt", "genbank")
```


```python
count = SeqIO.write(records, "my_example.fasta", "fasta")
print("Converted %i records" % count)
```

    Converted 94 records



```python
for record in SeqIO.parse("ls_orchid.gbk.txt", "genbank"):
    print(record.id)
    print(record.seq.reverse_complement())
```

    Z78533.1
    GCGTAAACTCAGCGGGTGCCCCCGCCTGACCTGGGGTCACATCCGAATGGCGGTCAACCGCCCTCCATTGGGGTTCGGAAGGGTTCTCCTGCCTGTCCGACAAGCACGACAACATGGGGGATTCGCACGGCAGCTGCTGCCAGCATCCGTCCACCTCTTGCCGGGTTCCGGGCCATCAAAACACCAGCTCTTGGACCCGCCGCACCTAGGCACAAGGGGCCAATCTTTCACATCCGCACCACGCCGGCCTGGCTGTATGCCGGGCAAGCATTGGCAGGAGAGAGACGAAGCGCGACGCCCAAGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAAGACTCGATGGTTCACGGGATCTTGCAATTCACACCACTTATCGCATTTCGCTGCGTCCTTCCATCCGATGCAAAGAGCCAAGATTCCCGTTTGCGAGAGTCATCAAAATTCATTGGGCACGCGACAGCACGCCGCCGCTCCGGGTTTTGGGGAAGACAATGCCATTCGCCGGTGATGCTTTCATATGGCTTGGCGCCCAAACTGCGCCGGGCTAGAGGTTCAAACCCGCCATGGACGCTCCCGAGGCGGCCCAACAACAAATCAGGGTCACCACGGGAGCAATGCCCCCGGTGAGCTGAGTACACCGGTCCTCCGGATTCACTCGATCGTTTATTCCACGGTCTCATCAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78532.1
    GCCTCAACTCAGCGGGTGGCCCCGCCTGACCTGGGGTCGCATCTGAATGGAAATCAACTGCCCAATGGTTATTTTAGCTCCATTGGGGTTCAATTAGGTTCTTGTGTAGGTTCGAAAAAATACAACAACATGGGGGATTCAAATAGCAGCCTTATGACTGTTAGCATTCTCCACCTCGTGCCACATTCCTACCCATCAAAGCAACAATCCTTAGACCCACCGCACCTAGGCACAAGGGGCCAATCATTCACATCCGTATAATGCCAGCTTAGCGATATGCCAAGCAAGCATTGGTAGGAGAGACGCAACACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGAGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCGCATTTAGCTGCGTTCTTCATCGATGTTAGAGCCAAGATATCCGTTGCGAGAGTCATAAAAATTCACTGGCATGCTCAGTAGCATACTGCCCCTCTGATTTTTTCTGACAATAATGTCATTCATCAGTGATGCTTTATATGACTTGGCGCAAACTGCACCGTACTAAAGTTCAAACCTGCCATGAAAGCTCTTGAGGAGGCCCAACAACAAAGCAGGGTCACGACAAAAGCAGTGCCACGACGAGCTGAGTTACCACAGGTCCTCCAGATTCACTCGATCATATATTCTGTTGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78531.1
    TTACTCACGGGTGGCCCGCCTGACTGGGGTCGCATCTGAATGGAAATCACCGCCCGAGGGCGGTTTTGCGTCCACTGGGGGTTCAAACAGGTTCTTCTGTAGGCCTGACGAGTACGACAACACGGGGGATTCGAACGGCAGCCTTGCGGCTGCCAGCATCCTCCACCTACTGCCAAGTTCCGGGCCATCAGAACACCGATGCTTAGACCCGCCGCACCTAGGCGCAAGGGGCCAATCTTTCACGTCCTCACAATGACAGACTAGCCGCATGCCAATCAAGCATTATCAGGAGAGACGCAGCACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCGCGGGATTCTGCAATTCACACCACTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCAGATATCCCGTTGCCGAGAGTCGTCATAAATTCATTGGCACGCGCAACAGACGCCGCCCCTCCGAGTTTTTCTTGACAAAAATGCCATCCATCGGTGACGCTCCATATGACTTGGCGCAAACTGCGCCGCGCTAGACGTTCAAACCCGCCATGAAAGTTCCCGAGGCGGCCCGGTCGCAAACCGGGTTCACCACGAGAGCAAAGCCACGGTGAGCCGTGTAACCACGGGTCCTCCGGATTCACTCGATCGTATGTTCTGCTGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78530.1
    ATGGGGTGGCCCCCTGACTGGGGTCGCATCTGATTGGAAATCAACCACCCAAGGATTGTTTTGCCTTAATTAGGGTCCAAACAAGTTGTTCTATAGGCCCAACAAACATGACAACATGGGGGATTCAAACAATAGCCTTGCGGCTGCCAGCATCCTCCACCTCTTGCCAAGTTTCAGACCATCAAAACACATATCCTTAGACCCACCGCACCTAAGCACAAGGGGCCAATCTTTCACATCCACACAATGATGGCCTAGCTATATGCTGGACAAGCATTGGAAGGAGAGATAAAACATACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTTCGGTGCGTTCTTCATCGATGCAAGAGTCAAGATATCCGTCTGCCGAGAGTCATCAAAATTCATTGGAATGAACATTAACACACCACACCTCCGACGTTGTGTGGGGCAATAATGTCATTCATCGGTGATTCTTTATATACCTTGGTGCAAACTGGACCGTGCTAGAGGTTCAAACCCGCCATGAAAGCTCTCAATGAGGCCCAATGACAAATCATGGTCACCACAAAAGGATATCCCTAGCGAGCCAAATTACCACAAGTCCTCCAGATTCACTCAATCGTTTATTATGTTGTTTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78529.1
    TTTTAGCTCAGCTGGTGGCCCGCTGACTGGGGTCGCATCTGATTGGAAATCAACCACCCAAGGATTGTTTTGCCTTAATTAGGGTCCAAACAAGTTGTTCTATAGGCCCAACAAATATGACAACATGGGGGATTCAAACAATAGCCTTGCGGCTGCCAGCATCCTCCACCTCTTGCCAAGTTTCAGACCATCAAAACACATATCCTTAGACCCACCGCACCTAAGCACAAGGGGCCAATCTTTCACATCCACACAATGATGGCCTAGCTATATGCTGGACAAGCATTGGAAGGAGAGATAAAACATACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTTCGGTGGGATTCTTCATCGATGCAAGAGCCAAGATATCCCTCGGCGAGAGTCATCAACAATTCATTGGCATGACCAATAACACACCACACCTCCGACTTTTTTGACAATAATGTCATTCATCGGTGATTCTTTATATACCTTGGTGCAAACTGCACCGTGCTAGAGGTTCAAACCCGCCATGAAAGCTCTCAATGAGGCCCAATGACAAATCATGGTCACCACAAAAGGAAATCCCTAGCGAGCCAAATAACCACAAGTCCTCCAGATTCACTCAATCGTATATTCTGCTGTCTCAACAATGTCCTTCGGCAGCTCGCCGT
    Z78527.1
    GGGTGGCCCGCCTGACCTGGGGTCGCATCTAAATGGAAATCAACCGCCAAGGGTCATTTTACAATCCATTGGGGTTCAAGCATGTTCTTTTATAGGTTCGACAAATATGACAACATGGGGGATTCGAACGACAGCCTTGCGGCTTCCAGCATCCTCCACCTCCTGCCAGGTTTCGATCCATCAAAACATCAATCCTTAGACCCACCGCACCTAGGCACAAGGGGCCAATCTTTCACATCCACCCAATGCCAGCCTAGTTGTATGCCGGGTAAGCATTGGCAAGAGAGATGCGATACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCGCGGGATTCTGCAATTCACACCACTTATCGCATTTCGCTGCGTTCCTTCATCGATGCAAAGAGCCAGATATCCGTTGGCCGAGAGTCATCAAAAATCCATTAGCGCATGCAATAACACACTGCCACTCCAACTTTTGGTGACAGTAGTGCCATTTGCAGTTATGCTTCATATGACTTGGCGCAAACTGCACCGTGGTAGAGGTTCAAACCCGCCATGTAAGCTCCTAAGGAAGCCCAACAACAAATCAGGCGAGCCGAGTAACCACAGGTCATCCAGATTCACTCGATCATATATTCTACTGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78526.1
    ACATAAACTCAGCGGGTGGCCCCGCCTGACCTGGGGTCGCATCTAAATGGAAATCAACCGCCAAGGGTCATTTTACATCCATTGGGGTTCAAGCATGTTCTTTTATAGGTTCGACAAATATGACAACATGGGGGATTCGAACGACAGCCTTGCGGCTTCCAGCATCCTCCACCTCCTGCCAGGTTTCGATCCATCAAAACATCAATCCTTAGACCCACCGCACCTAGGCACAAGGGGCCAATCTTTCACATCCGCACAATGCCAGCCTAGTTGTATGCCGGGTAAGCATTGGCAAGAGAGATGCGATACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCGCGGGATTCTGCAATTCACACCACTTATCGCATTTCGCTGGGTTCTTCATCGATGCAAGAGCCAAGATATCCGTTGGCGAGAGTCATCAAAAATCCATTAGCGCATGCAATAGCACACTGCCACTCCAACTTTTGGTGACAATAATGCCATTGTTCAGTTATGCTTCATATGACTTGGCGCAAACTGCACCGTGGTAGAGGTTCAAACCCGGCATGTAAGCTCCTAAGGAAGCCCAACAACAAATCAGGGGAGCCGAGTAACCACAGGTCCTCCAGATTCACTCGATCATATATTCTACTGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78525.1
    TGCAGCGGGTGGCCCGCCTGACCTGGGGTCGCATATGAATGGCAGTCAACCATCCAAGGGTCATTTTGCCTCCACTGGGGTTCAAACAAGTTATTATATAGGTCCGATAAATACGATAACATGGGGGATTTAAATGACAACCTTGTGGCAGCCAATGTCCTCCACCTCCTGCCAAGTTTCGGGCCATCAAAACATTCCTTAGACCCAACGCGCCAAGGCACAAGGGGCCAATCTTTCGCATCCGCACAATGCAAGCCTAGATATATGGACAGGCATTGGCAAGAGAGACACAACACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAAGACTCGATGGTTCACGGGATCTTGCAATTCACACCACTTATCGCATTTCGCTGGGTTCTTCATCCGATGCAAGAGCCAAGATATCCGCTTGTCGAGAGTCATCAAAATTATATTCGGCATGCACAGGAGGACACCGCCCCTCCAACTTTTTTGACAATAATGTCATTTGTCGGTGATGCTTCATATGACTTGGCACAAACTGTGCCGTGCTAGAGTTTGAAACCTGCCATGAAAGCTCCCGAGGAGGCCCAACGACAAATTTGGGTCACCACAAAAGCAAAGCCCTCGGCAAGCCGAATAACCACAGGTCCTCCGGATTCACTCGATGTATATTCTGCTATCTCAACA
    Z78524.1
    GCTTAATCTCAGGGTGGCCACCTGACCTGGGGTTGCATCTGAATGAAAATCAACCGCCCAAGGGTCATTTTTCCTTCATTGGGGTTCAAACAAGTTCTTTTATAGGTTCAACAAATACGACAACATGGGGAATTCAAATGGCAGCCTTGCAGCTACCAACATTCTCCACCTCCTGCCGAGTTCCGAGCTATCAAAATACTAATCCTTAGACCCACCGCACCTAAGCACAAGGGGCCAATATTTCATCCGCACAATGCTAACCTAGCTATATGCTGGGCAAGCATTGGNAGGAGAGATGCAGCACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTTGGGCGCAACTTGCGTTCAAATACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCGCATTTCGCTGCGTTCTTCATCGGATGTAAAGAGCCAAGTTCCGTTGCGAGAGTCATCAAAAAATTCATTAGCATGCACAACAGCATGCCGCCCCTCCGACTTTTTTAACAACAATGCCATTCATCAGGGATGCTTATATGACTTGGCACAAACTGCGCCGTGCTAAAGGTTTGAACCCGCCATGAAAGCTCCTAAGGAGGCCCAATTAAGGTCACCACAAAAGCAAAGCACTGACGAGCCAAGTAACCACATGTCCTCCATATTCACTCAATCATATATTCTACTATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78523.1
    CTTAAACTCAGCGGTGGCCCCGCCTGACCTGGGTCGCAATATGAATGGCAATCAACCGCCCGAGGGTTGTTTGCCTCCATTGGGATTCAAACAGGTTCTTCTCTAGGTCCGACAAACACGACAATATGGGGGGTTCAAACGATAGCCTTATAGCTGCCAACATCCTCCACCTCCTGCCAGGTTTCGAACCATCAAAACACTAGTCCTTAGACCCACCGCACCTAGGCACAAGGGGCCAATCTTTCGCATCCACGCAATGCAAACCAATCTATATGATAGGAAAGCATTGATAGGAGAGACGCAGCACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCATCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACATATGTCGAGAGTCAATAAAAAATTCATTATGCATGCATGCAAGAGCACGCTTCCACTCTAGCTTTTGGGACAATAATGTCATTCCTCAGTGATGCTTCATATGACTTGGCGCATACTGCACCGTGCTAGAGGTTCACACCCACAAGGAAAGCTTGGGAGGAGGCCCAATGACAAGTTAGGGTCACCGCAAAAGCAAAGCCTATGTCGAGCTGAGTAACCACAAGTCCACCGGATTCACTCGATCATATATTCTGCTGTCTCAACAATGATCCTGCCGCAGGTTCACCTACGGAAACCTGGTTACG
    Z78522.1
    CTCAGCGGGTGGCCGCCTGACCTGGGGTCGCATATGAATGGCAATCAACCGCCCGAGGGTTGTTTGCCCCCATTGGGATTCAAACATGTTCTTCTCTAGGTCCGACAAACACGACAATATGGGGGATTCAAATGATAGCCTTATAGCTGCCAACATCCTCCACCTCCTGCCAGGTTTCGAACCATCAAAACACTAGTCCTTAGACCCACCGCACCTAGGCACAAGGGGCCAATCTTTCACATCCACGCAATGCAAACCTGTCTATATGACGTATAACCATTGACAGGAGAGACGCAGCACACGACGCCCAGGCAGGCGTTCCCTTAGCCTGATGGCATCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCGATCACACCACGATATGTCGAGAGTCATAAAAAATTCATTTGCATGTATGCAATAGCACGCTTCCACTCCAACTTTTTGGACAATAATGTCATTCATCAGTGATGCTTCATATGACTTGGCGCATACTGCACCGTGCTAGAGGTTCAAACCCACAAGGAAAGCTTGGGGGGAGGCCCAATGACAAATTAGGGTCACCGCAAAAGCAAAGCCTATGTCGAGCTGAGTAACCACAAGTCCACCGGATTCACTCGATCATATATTCTGCTGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78521.1
    GGTGGCCCCGCCTGACCTGGGGTCGCATATGAATGGCAATCAACCGCCCGAGAGTTGTTTGCCTCCATTGGGATTCAAACATGTTTTTCTCTAGGTCCGACAAACACGACAATATGGGGGATTCAAATGATAGCCTTATAGCTGCCAACATCCTCCACCTCCTGCCAGGTTTCGAACCATCAAAACACTAGTCCTTAGACCCACCGTACCTAGGCACAAGGGGCCATTCTTTCACATCCACGCAATGCAAACCTATCTATATGACATGAAAGCATTGACAGGAGAGACGCAGCACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCATCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCGCATTTTGCTGCGTTCTTCATCGATGCAAGAGCCAAGATATCCGTTGTCGAGAGTCATAAAAAATTCATTTGCATGCATGCAATAGCACGCTTCCACTCCAACTTTTTGACAATAATGTCATTCATCGGTGATGCTTCATATGACTTGGCGCATACTGCACCGTGCTAGAGGTTCAAACCCACAAGGAAAGCTTTGGAGGAGGCCCAATGGCAAATTAGGGTCACCGCCAAAGCCAAGCCTATGTCGAGCTGAGTAACCACAAGTCCACCGGATTCACTCGATCATATATTCTACTGTCTCAACAATGATCCTTCCGCAGGTTCACCTAC
    Z78520.1
    AAACTCAGCGGGTGGCCCGCCTGACCTGGGGTCGCATATGAATGGCAATCAACCGCCCGAGGGTTGTTTGCCTCCATTGGGATTCAAACAGGTTCTTCTCCAGGTCCGACAAACACGACAATATGGGGGATTCACATGACAACCTTATAGCTGCCAACATCCTCCACCTCCTGCCAGGTTTCGAACCATCAAAACACTAGTCCTTAGACCCACCGCACCTAGGCACAAGGGGCCAATCTTTCACATCCACACAATGCAAACCTATCGATATGACAGGAAAGCATTGACAGGAGAGACGCAGCACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCATCGGGCGCAACTTGCGTTCACATACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACAATTTCGCTGCGTTCTTCATCCGATGCAAGAGCCAAGATATCCCTTGTCGAGAGTCATAAAATATTCTATTGCATGCATGCAATAGCACGCTTCTACTCCAACTTTTTGGACAATAATGTCATTCATCGGTGATGATTCATATGACTTGGCGCATACTGCACCGTGCTAGAGGTTCAAACCCACAAGGAAAGCTTTGGAGGAGGCCCAATGACAAATTAGGGTCACCGCAAAAGCAAAGCCTATGGCGAGCTGAGTAACCACAAGTCCACCGGATTCACTCGATCATATATTCTGCTGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78519.1
    TAAACTCAGCGGGTGGCCCCGCCTGACCTGGGGTCGCATATGAATGGCAATCAACCGCCCGAGGGTTGTTTGCCTCCATTGGGATTCAAACATGTTCTTCTCTAGGTCCAACAAACGCGACAATATGGGGGATTCAAATGATAGCCTTATATAGCTGCCAACATCCTCCACCTCCTGCCAGGTTTCGAACCATCAAAACATTAAGTTCTTAGACCCACCGCACCTAGGCACAAGGGGCCAATCTTTCACATCCACACAATGCAAACCTATCTATATGATGGGAAAGCATTGACAGGAGAGACGCAGCACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCATCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCAAGATATCCGTTGTCGAGAGTCATAAAAAAATTCATTTGCATGCATGCAGTAGCACGCTTCCACTCCAACTTTTGGACAATAATGTCATTCATCGGCAATGCTTCATATGACTTGGTGCATACTGCACCGTGCTAGAGGTTCAAACCCACAAGGAAAGCTTGGGAGGAGGCCCAATGACAAATTAGGGTCACCGCAAAAGCAAAGCCTATGGCGAGCTGAGTAACCACAAGTCCACCAGATTCACTCGATCATAT
    Z78518.1
    GGAAAACTCAGCGGGTGGCCCCGCCTGACCTGGGGTCGTATCTGAATGGCAATCAATCGCTCAAGGGTTGTTTTGCCTCCATAGGGGTTCAAATAGGTTCTTCTGTGGGTCCGGCAAATACGACAACATGTGGGATTTGAACGACAGCCTTGCGACTATCAACATGCTCCACCTCCTGCCGGGGTTTGGGCCATCAAAACACCAATCCTTAGACCCACCGCACCTAAGCACAAGGGGTCAATCTTTCATATCTATGCAATACTGGCCTAATTGTATGTCGAGCAAGCATTGGCAAAAGAGATGCAGCACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTTGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGTATGCGTAACACACACCGTCCATCCGAGTTTGGGGGACAATAATGTAATTCGTCGGCGATGCTTCATATGACTTGGCGCAAACTGCGCCGTGATAGAGGCTCAAACCCGCCATAAAAGCTCTCGAGGATGCCCAACTACAAATCAGGGTCACCACAAAAGTTAAGCCTTCGACGAGCCGAGTAACCACAAGTCCTCCGGATTCACTCGATCGAATATTCTACTATCTCAACAATGATCCTCCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78517.1
    GCTTAAACTCAGCGGGTGGCCCCGCCTGACCTGGGGTCGTATCTGAATGGCAATCAATCGCTCGAGGGTTGTTTTGCCTCCATTGGGGTTCAAATAGGTTCTTCTGTGGGTCCGACAAATACGACAACATGTGGGATTTGAATGACAGTCTTGCGACTATCAACATGATCCACCTCCTGCCGGGGTTCGGGCCACCAAAACACCAATCCTTAGACCCACCGCACCTAAGCACAAGGGGTCAATCTTTCATATCTATGCAATACTGTCCTAATTGTATGTCGGGCAAGCATTGGCAAAAGAGATGCAGCACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTTGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCGCATTTCGCTGCGTTCTTCATCGATGGAAGAGCAAGATATCCGTTGGCCGAGAGTCATCAAAAATTTATTGGCATGCACAACAGCACACCGCCCATCCGACTTTTGTGACAATAATGTCATTCATCGGTGATGCTTCATATGACTTGGCGCAAACTGCACCGTGATAGAGGTTCAAACCCGCCATAAAATCTCCTGAGGATGCCCAACGACAAATCAGGGCCACAAAAGCTAAGACTTGGACGAGCCGAGTAACCACAAGTCCTCTGGATTCACTCGATCGTATATTCTACTGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78516.1
    TTAAACTCAGCGGGTGGCCCCGCCTGACCTGGGGTCATATCTGAATGGCAATCAATCACTCGAGGGTTGTTTTGCCTCCATTGGGGTTCAAATAGGTTCTTCTGTGGGTCCGACAAATACGACAACATGTGGGATTTGAACGATAGTCTTGCGACTATCAACATGCTCCACCTCCTGCCGGGGTTCGGGACATCAAAACACCAATCCTTAGACCCACCGCACCTAAGCACAAGGGGTCAATCTTTCATATCTATGCAATACTGGCCTAATTGTATGTCGGGCAAGCATTGGCAAAAGAGATGCAGCACACGACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTTGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCGCATTTCGCTAAGATATCCGTTGACGAGAGTCATCAAAAAATTATTGGTATGCACAATAGCACACCGCCCATCCGACTTTGGTGACAATAATGTCATTCGTCGGTGATGCTTCATATGACTTGGCGCAAACTGCGCCGTGATAGAGGTTCAAACCCGCCATAAAAGCTCCTGAGGATGCCCAACGACAAATCAGGGCCACAAAAGCTAAGACTTGGACGAGCCGAGTAACCACAAGTCCTCCGGATTCACTCGATCATATATTATACTGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78515.1
    GCTTAAACTCAGCGGGTGGCCTCACCTGACCTGGGGTTGCTTCCAAATGGCAATCAACCGCCCATGGGTTGGGTTGATGTGGCACCAATGGGGTTCCGAAGGGTCCTACGAATTATGATGGCCCATCAAATCGGACAACATACGGCATTCGCACAACAGCTTTGGCTCAGCCATTGCCCATCCACCTCTTGTCGAATTTTGGACCATCAAAAGCCCGAGGTCTTAGACCCATTGAACCGAAGAACAAGGGGCCAAACTCACATCCGCACCGACCGAACGGTATTCTTGGACAACCATTGGTGGGAGAGACACAGTACACAACGCCCAGGCAGGCGTGCCCTTAGCCTGACGGCCTCGAGCGCAACTTGCGTTCAAAGACTCGATGTTCACGGGATTCTGCAATTCACACCACTTATCGCATTTCGCTGCGTCCTCCATCGATGCAAGAGCCGAGATATCCGTTGCCGAGAGTTGTAAAAATGTAACTTGGGGGGCACGAATCAAGTGCCACCCCCACCGGCTACAAGGGAACACCCCTATGCCGATTCTTGTCTCAAAAGACTTGGCGCGAAGCTGCGCCGGGCTAGAAGTTTGATCGCCTTCAAGAACACTCCCAAGGAGGCCCGATGGCAGATCAAAGGCCACCACAGTAGCGAAAGCCCCTGGGAGATCGAAGTACACCGGTCCTCCAGATTAACTCGATCGTTTATTCTACGGTCTCAGCAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78514.1
    TAGCGGGTAGTCTCACCTGACCTGGGGTTGCATCCTAATGGCCGTCAACCGCCCATGGGTTGGGTCGACGTGCCTCCGACGGGGTTCGAAAGGGATCATTCCGGCCCATCAAATACGACAACATAGGGCATTCGCACAACAGCTTTGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTGTGGACCGTCACAAGCCCAAGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACTCTCATCCTCGCCGGCCGAACAGTATGGCTGTTCAACCTTAGCATTGGGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTTCGCTGGGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGACGAGAGTTGTCAGACGAATCTTAAGGGGGCACAGCAGCACGCCGGCCCTCCCCGGTCCATCCATGGAGGACGAAGACCCTTGTCGGTTCTATCTCATTTGACTTGGAGCGAAACTGCGCCGGGCCAAGGGTTCAATTGCCATAGAGAACTCTCCCAAGGAGGTTCGATGGGAAATCATGGTCACCACAGCGGCGGAAGCCCCTGGGAGACCAAAGTACACCGGTCCTCCGGATTATCTCGATCGTATATTTGGGGGCTTCAAAAATGATCCTCCCGAAGGTCCACCTACGGAAACCTTGTTACG
    Z78513.1
    CTCAGCGGTAGCTCACTGACTGGGTTGCATCCTAATGGCGTCACCGCCCATGGGTTGGGTCGACGTGCCTCCGAAGGGGTTCGAAAGGGATCGTTCTGGCACATCAAATACGACAACATAGGGCATTCGCACAAAAGCTTCGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTGTGGACCGTCACAAGCCCAAGTCTTGGACCCATCGCACAGAAGAACAAGGGGCCAAACTCTCATCCTCGCCGGCCGAACAGTATGGCTGTTCAACCTTAGCATTGGGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCCGAGAGTTGTCGAAAATTCTTTCGGGCACAGCAGCATGCCGGCCTCCCCGGCTCCATGGAGGACGATGCCCTTGCCGGTTCTATCTCATTTGACTTGGCGCGAAACTGCGCCGGGCCAAGGGTTCAATTGCCATGGAGAAATCTCCCAATGAGGCTCGATGGCAAATCATGGTCACCACAGCGGCGGAAGCCCCTGGGAGACCAAAGTACACCGGTCCTCCGGATTAACTCGATCGTGTATTTGGCGGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78512.1
    GCTTAAACTCAGCGGGTGGCCTCACCTGACCTGGGGTTGGATCCAAATGACCGTCGACCGCCCACGGGTTGACGTGCCTCCAATGGGGTTCAAAAGGGATTTATTCCGGCCCATCAAATACGACAACGCAGGGCATTCGCACGACAGCTTCGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTCCGGCCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACACTCATCCGCGCCGGCCGAACAGTATGCCTGTTCAGCCTTAGCAATGGGAGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCAGCTTATCACATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCGAGAGTTGTCAAAAAATAATTTCATTGGGGGAACGGAAGAACGCCGGCCCTCCCCGGTTCCATGGGGGACGATGCCCTCCGCCGGTTCCATCTCATTTGATTTGGGGCGAAACTGCGCCGGGCCAAGGGTTTCATTTGCCATCAAGAATTCCTCCCAGAGGGTTCGACGGAAATTCACGGTCACCACAGTAGCCAAGGCCCCTGGGAGACCAAACTACACCGGCCCTCCGGATTAATTCGATCGTTTTTTTGGGGGCTTCAAAAATGATCCTCCCGAAGGTCCACCTACGGAAACCTTGTTACG
    Z78511.1
    TCCTCACCTGACCTGGGGTTGCATCCAAATGGCCGTCGACCGCCCACGGGTTGACGTGCCTCCAATGGGGCTCAAAAGGGATTTATTCCGGCCCATCAAATACGACAGCGTAGGGCATTCGCACGACAGCTTCGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTCCGGCCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACACTCATCCGCGCCGGCCGAACAGTATGCCTGTTCAGCCTTAGCAATGGGAGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTTCGCTGGGTTCTTTCATCGGATGCAAAGAGCCGAGATATCCGTTGCCGAGAGTTGTCAAAAAAAAATTCATGGGGGGAACGGAAGAACGCCGGCCCTCCCCGGTTCCATGGAGGACGATGCCCTCCGCCGGTTCCATCTCATTTGACTTGGGGCGAAACTGCGCCGGGCCAAGGGTTCAATTGCCATCAAGAATTCTCCCAAGGAGGTTCGACGGAAATTCACGGCCACCACAGTAGCCAAAGCCCCTGGGAGACCAAACTACACCGGTCCTCCGGATTAACTCGATCGTTTTTTTGGGGGTCTCAACAATGATCCTTCCGAAGGTTCACCTACGGAAACCTTGTTACG
    Z78510.1
    TAAACTCAGCGGGTGGCTCACCTGACCTGGGGTTGCATCCAAATGGCCGTCAACCGCCCATGGGTTGACGTGCCTCCAATGGGGTTCAAAAGGGATTTATTCCGGCCCATCAAAAACGACAACGCAGGGCATTCGCACGACAGCTTTGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTTGGACCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAATCTCTCATCCGCGCCGGCCGAACAGTATGCCCGTTCAACCTTAGCAATGGGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTTCGCTGCGTTCTTCCATCCGATGCAAAGAGGCGAGATATCCGTTGCGAGAGTTGTCAAAAAATCCGAGGGGGGAACGGAAGAATGCCGGCCCTCCCCGGTTCCATGGGGGGGGACGCCCCCTGCCGGTCCTATCTCATTTGATTGGGGGGGAAATTGCGCCGGGCCAAGGGTCCATTGGCCACCAAGAATTCTCCCAGGGGGGTTCGATGGAAATTCACGGCCACCAAGGGGGGAAAGCCCCTGGGGAGACCAAATTAAACCGGTCCTCCGGTTTAATTCGATCGTTTTTTCGGGGGCTTAAAAAGGAATCCTCCCGAAGGTCACCTCGGAACCCTGGTTAG
    Z78509.1
    TCCTAAAATTAACGGGTGCCTTAACTTGCCGGGGTTAACCAAATGCCCGTCACCCCCCATTGGGTTGACGTGCCTCCAATGGGGTTCAAAAGGGATTTATTCTGGCCCATCAAAAACGACAACGTAGGGCATTCGCACGACAGCTCTGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTTGGACCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACTCTCATCCGCGCCGGCCGAACAGTATGCCCGTTCAACCTTAGCAATGGGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAACCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCTGCTTATCACATTTAAACCTTGTTACGTCCGTTGCCGAGAGTTGTCAGAAAAATTCATGGGGGGGCACGGAGCATGCCGGCCCTCCCCGCTCCATGGAGGACGACGCCCTCCGCCGGTTCTATCTCATATGACTTGGCGCGAAACTGCGCCGGGCCAAGGGTTCAATTGCCATCAAGAAATCTCCCAAGGAGGCTCGATGGCAAATCACGGTCACCACAGCGGCGAAAGCCCCTGGGAGACCAAACTACACCGGTCCTCCGGATTAACTCGATCGTATATTTGGCGGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78508.1
    TCAGCGGTGGCTCACTGACTGGGTTGCATCCAAGTGGCCGTCACCGCCCATGGGGTTGACGTGCCTCCAATGGGGTTCAAAGGGATTTATTCCGGCCCATCAAATACGACAACGTAGGGCATTCGCACGACAGCTCTGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTTGGACCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACTCTCATCCGCGCCGGCCGAACAGTATGCCCGTTCAACCTTAGCAATGGGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCCGAGAGTTGTCAAAAAATTCATTGGGGGCACGGCAGCATGCCGGCCCTCCCCGGCTCCATGGAGGACGACGCCCCCTGCCGGTTCTATCTCATATGACTTGGCGCGAAACTGCGCCGGGCCAAGGGTTCAATTGCCATCAAGAAATCTCCCAAGGAGGCTCGATGGCAAATCACGGTCACCACAGCGGCGAAAGCCCCTGGGAGACCAAACTACACCGGTCCTCCGGATTAACTCGATCGTATATTTGGCGGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78507.1
    TCAGCGGTGGCTCACTGACTGGGTTGCATCCAAATGGCCGTCGACCGCCACGGGTTGACGTGCCTCCAATGGGGTTCAAAGGGATTTATTCCGGCCCATCAAATACGACAACGCAGGGCATTCGCACGACAGCTTCGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTCCGGCCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACACTCATCCGCGCCGGCCGAACAGTATGCCTGTTCAGCCTTAGCAATGGGAGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTTCGCTGTGTTCTTCATCGATGCAAGAGCGAGATATCCGTTGCGAGAGTTGTCAAAAAAATTCATTGGGGGACGGAAGCACGCCGGCCCTCCCCGGTTCCATGGAGGACGATGCCCTCCGCCGGTTCCATCTCATTTGACTTGGGGCGAAACTGCGCCGGGCCAAGGGTTCAATTGCCATCAAGAAATCTCCCAAGGAGGCTCGACGGAAATTCACGGTCACCACAGTAGCCAAAGCCCCTGGGAGACCAAACTACACCGGTCCTCCGGATTAACTCGATCGTTTTTTTGGGGGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78506.1
    TCAGCGGGTGGCCTCACCTGACTGGGTTGGCATCCAATGGCCGTCAAACCGCCCATGGGTTGACGTGCCTCCAATGGGGTTCAAAAGGGATTTATTCTGGCCCATCAAATACGACAACGTAGGGCATTCGCACGACAGCTTTGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTTTTACCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACTCTCATCCGCGCCGGCCGAACAGTATATGCCTGTTCGACCTTAGCAACGGGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTATGAGCTGAGGGCTCCTTGAACGAGAGTTGTCAGCCCGATTCAAAGGGGGTACGGCGGCATGCCGGTCCTCCCCGGTTCCATGGAGGACGAAGCCCTCTGCCGGTTCTATCTCATACGACTTGGCGCGAAACTGCGCCGGGCCAAGGGTTCAATTGCCATCAAGAAATCTCCCGAGGAGGCTCGATGGAAAATCATGGTCACCGCAGTAGCCAACGCCCCTGGGAGACCAAACTACACCAGTCCTCCGGATTAACTCGATCGTATATTTTGCGGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78505.1
    AAAACTCAGCGGGTGGCTCACCTGACCTGGGGTTGCATCCAAATGGCCGTCAACCGCCCATGGGTTGACGTGCCTCCAATGGGGTTCAAAAGGGATTTATTCCGGCCCATCAAATACGACAACGCAGGGCATTCGCACGACAGCTTTGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTTTTACCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACTCTCATCCGCGCCGACCGAACAGTATGCCTGTTCGACCTTAGCAACGGGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTGACGAGAGTTGTCACAAAGATTCATTGGGGGTATGGTGGTATGCCGGGCCTCCCCGGTTCCATGGAGGACGAAGACCTCTGCCGGTTCTATCTCATACGACTTGGCGCGAAACTGCGCCGGGCCACGGGTTCAATTGCCATCAAGAAATCTCCCGAGGAGGCTCGATGGAATATCATGGTCACCGCAGTAGCCACCGCCCCTGGGAGACCAAACTACACCAGTCCTCCGGATTAACTCGATCGTATATTTGGCGGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78504.1
    TTACTCAGCGGTGGCTCACTGACTGGGGTTGCATCCAATGGCCGTCACCGCCCATGGGTTGACGTGCCTCCAATGGGGTTCAAAAGGGATTTATTCTGGCCCATCAAATACGACAACGTAGGGCATTCGCACGACAGTTTTGTCGCAGCCACCGTCCGTCCACCTCTTGTCGGATTTGGTACCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACTCTCATCCGCGCCGGCCGAACAGTATGTCTGTTCGACCTTAGCAACGGGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTTCGCTGCGTTCTTCATCGATGAAGAGCGAGATATCCGTTGCGAGAGTTGTCAAAAAGATTCATTGGGGGCACGGCGGGCATGCCGGCCCTCCCCGGCTCCATGGAGGGACGATGCCCTCTGACGGTTCTATCTCATACGACTTGGCGCGAAACTGCGCCGGGCCAAGGGTTCAATTGCCATCAAGAAATCTCCCGAGGAGGCTCGATGGAAAATCATGGTCACCGCAGTAGCCAACGCCCCTGGGAGACAAACTACACCAGTCCTCCGGATTAACTCGATCGTATATTTTGCGGTCTCAACAATGATCCTTCCGAAGGTTCACCTACGGAAACCTTGTTACG
    Z78503.1
    TTATCTCAGCGGGTGGCTCACCTGACCTGGGGTTGCATCCAAATGGCCGTCAACCGCCCATGGGTTGACGTGCCTCCAATGGGGTTCAAAAGGGATTTATTCTGGCCCATCAAATACGACAACGTAGGGCATTCGCACGACAGCTTTGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTTGGACCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACTCTCGTCCGCGCCGGCCGAACAGTATGCCCGTTCAACCTTAGCAATGGGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACAAAGAGCTGAGACTCCGTTGCCGAGAGTTGTCAAAATAATTCATTGGGGGTACGGTAGCATGCCGGCCCTCCCCGGTTCCATGGAGGACGACGCCCCCTGCCGGTTCTATCTCATATGACTTGGCGCGAAACTGCGCCGGGCCATGGGTTCAATTGCCATCAAGAAATCTCCCAAGGAGGCTCGATGGTAAATCACGGACACCACAGCGGCGAAAGCCCCTGGGAGACCAAACTACACCGGTCCTCCGGATTAACTCGATCGTATATTTGGCGGTCTCAACAAGGATCCTTCCGGAGGTTCACCTACGGAAACCTGGTTACG
    Z78502.1
    GCGTAAACTCACGGGGTGGCCTCACCTGACCTGGGGTTGCATCCAAATGGCCGTCAACCGCCCATGGGTTGACGTGCCTCCAATGGGGTTCAAAAGGGATTTATTCTGGCCCATCAAATACGACAACGTAGGGCATTCGCACGACAGCTTTGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTTGGACCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACTCTCGTCCGCGCCGGCCGAACAGTACGCCCGTTCAACCTTAGCAATGGGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACATTTTCGCTGCGTCCTTCCATCGATGCAAGAGCCGAGATTCCGGCCGAGAGTTGTCAATATAAAATTCATTGGGGGGACGGTAGTATGCCGGCCCTCCCCGGCTCCATGGAGGACGACGACCCCTGCCGGTTCTATCTCATATGACTTGGCGCGAAACTGCGCCGGGTCATGGGTTCAATTGCCATCAAGAAATCTCCCAAGGAGGCTCGATGGCAAATCACGGTCACCACAGGGGCGAAAGCCCCTGGGAGACCAAACTACACCGGTCCTCCGGATTAACTCGATCGTATATTTGGCGGTCTCAACAATGATCCTTCCGGAGGTTCACCTACGGAAACCTGGTTACG
    Z78501.1
    TCTTAAACTCAGCGGGTGGCCTCACCTGACCTGGGGTTGCATCCAAATGGCCGTCGACCGCCCACGGGTTGACGTGCCTCCAATGGGGTTCAAAAGGGATTTATTCCGGCCCATCAAATACGACAGCGTAGGGCATTCGCTCGACAGCTTCGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTCCGGCCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACACTCATCCGCGCCGGCCGAACAGTATGCCTGTTCAGCCTTAGCAATGGGAGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAAGACTCGATGGTTCACGGGATTCTGCAATTCACAACACTTATCACAATTTCGCTGCGTCCTTCCATCCGATGCAAGGAGCCGAGATTTCCCGTTTGGCCGAGAGTTGCCCAAAAAAAAATTCATTGGGGGCACGGCAACACGCCGGCCCTCCCCGGCTCCATGGAGGACGATTCCCTCCGCCGGTTCCATCTCATATGACTTGGCGCGAAACTGCGCCGGGCCAAGGGTTCAATTTCCATCAAGAAATCTCCCAAGGAGGCTCGACGGCAAAGTCACGGTCACCACAGTAACCAAAGCCCCTGGGAGACCAAACTACACCGGTCCTCCGGATTAACTCGATCGTATATTTTGCGGTCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78500.1
    CTTAAACTCAGCGGGTGGCCTCACCTGACCTGGGGTTGCATCCAAATGGCCGTCGACCGCCCACGGGTTGACGTGCCTCCAATGGGGTTCAAAAGGGATTTATTCCGGCCCATCAAATACGACAGCGTAGGGCATTCGCACGACAGCTTCGTCGCAGCCACCGTCCGTCCACCTCTTGCCGGATTTCCGGCCGTCAAAAGCCCAGGTCTTGGACCCATCGCACCGAAGAACAAGGGGCCAAACACTCATCCGCGCCGGCCGAACAGTATGCCTGTTCAGCCTTAGCAATGGGAGAGAGATGCAGCACGCAACGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCACAATTCGCTGCGTCCTTCCATCCGATGCAAAGAGCCGAGATTCCCGTTGGCCGAGAAGTTTTCCCAAAGAAAAATTCATTGGGGGCACGGCAGGACGCCGGCCCTCCCCGGCTCCATGGAGGGCGATGCCCTCCGCCGGTTCCATCTCATATGACTTGGCGCGAAACTGCGCCGGGCCAAGGGTTCAATTGCCATCAAGAAATCTCCCAAGGAGGCTCGACGGCAAAGTCACGGTCACCACAGTAACCAAAGCCCCTGGGAGACCAAACTACACCGGTCCTCCGGATTAACTCGATCAATTATTTTGCGGTCTCAACAATGAGCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78499.1
    GGTTGCTCACCTGACCTGGGGTCGCATCCAAAAGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCACAAGGGTCTCTAGATTATGATGGCCCATCTAATACGACAACCTGGAGCATTCGCCCAACAGTTTTGTTGTAGCATGATTCGTCCACCTCTTGCCGAATTTAGGACCATCAAAAGCCCCTGCAGCTCTTAGACCCCCAGCACCAAAGAACAAGGGGGCAAACTCACATCCGCACCGGCTGAACAGTATGCATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCTTGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATCTTGCAATTCACACACTTATCGCATTTCGCTGCGTCCTCCATCCGATGCAAAGAGCCGAGATATCCGTTGGCTGAGAGTTGTCAAAAAAATAATTGGGAGAGGTCAAGGCCACATGCCGTCCCCTCCTCGAATTCAACGAAGGGGGTATGTCCTTCACCATTTATGTGTCATATGACTTGGTGCAAAACTGCGACGGGCAAGGGTTAGATCGCCGGCAAGAAAGCTCCCAAGGAGGCTCCCTCCATGGAAAATCTAGGTCACCGCAATAGCAGGCGCCCAAGGGTGACCAAAGTAAACCGATCCTCTGGATTATCTCGATCAATTATTATGCGATCTCAACAATGATCCCTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78498.1
    GCTTAAACTCAGCGGGTTGCTCACTGACTGGGGTCGCAACAAATGGCCATCAACTGATCACGGGTTGATGTGCCTCCAGTGGGGTTCACAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCTGGAGCATTCGCACAACAGTTTTGTTGTAATTCGTCCACCTCTTGCCGAATTTAGGACCATCAAAAGCCCCTGCAGCTCTTAGACCCCCAGCACCGAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGCATGTACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTTGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATCCTGCAATTCACACCACTTATCGCATTTCGCTGGGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTCCAAAAATAATTTGGAGAGGGGGAGTGATTGTAGGGTCAAGGGCAAAATGCCGCCCCCTCCTTCGCCATTATTGTGTCATATGACTTGGGGCAAAACTGCCCCGGGCAAGGGTAAGATCGCCGGCAAGAAAACTCCCAAGGAGGCTCCATGGCAAATCTAGGTCACCGCAATAGCAAGCACCCATGGGTGACCGGAGTAAACTGATCCTCTGGAGTAACTCGATCAATTATTATGTGATCTCAACAATGACCTTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78497.1
    GCTTAAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCACAAGGGTGTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAGTTTTGTTGTAGCGCGTTGCGTCCACCTCTTGCCGAATTTAGGACCATCAAAAGCCCCTGCAGCTCTTAGACCCCCAGCACCAAAGAACAAGGGGCCAAACTCACACCCGCACCGGCTGAACAGTATGCCTGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTTGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCGCATTTCGCTGCGTTCCTTCATCGATGCAAGAGCTGAGATCTCCGTTGCTGAGAGTTGTTTCACAAAAAATAATTTGGAGAGGGGGAGGGAGTGTAGGCTCAAGGCCACATGCCCTTCACCTCCTTGAATTCACAAAGGGCCGTGCCCTTCACCATTTATGTGTCATGTGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTAGATCGCCGGCAAGAAAGCTCCCAAGGAGGCTCCATGGCAAATCTAGGTCACCGCAATAGCAAGCGCCCATGGGCGACCAAAGTAAACTGATCCTCTGGATTCACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78496.1
    GCTTAAACTCAGCTGGTTGCCTCACCTGACCTGGGGTCGCAACCCAAAAGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCACAAGGGTCTCTAGATTATGATGGCCCATCTAATACGACAACCTGGAGCATTCGCCCAACAGTTTTGTTGTAGCATGATTCGTCCACCTCTTGCCGAATTTAGGACCATCAAAAGCCCCTGCAGCTCTTAGACCCCCAGCACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGCATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTTGGGCGCAACTTGCGTTCAAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCGCATTTCGGTGCGTCCCTCCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTCAAAAAATTAATTGGAGAGGTCAAGGCCACATGCCGGCCCCTCCTCGAATTCACGAAGGGATATGCCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTAGATCGCCGGCAAGAAAGCTCCCAAGGAGGCTCCCTCCATGGCAAATCTAGGTCACCGCAATAGCAGGCGCCCAAGGGTGACCAAAGTAAACCGATCCTCTGGATTAACTCGATCAATTATTATGCGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78495.1
    CACTCAGCGGGTTGCTCACTGACCTGGGGTCGCAACCGAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCACAAGGGTCTCTGGATTATGCCGGCCCATCTAATACGACAACCTGGAGCATTCGCACAACAACAGTGTTGTTGTAGCCTGATTCGTCCACCTCTTGCCGAATTTAGGACCATCAAAAGCCCCTGCTGCTCTTAGACCCCCAGCACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGCATGGACAGCCTCATTGAGGGAGAGAGATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTCATGGCCTCGGGCGCAACTTGCGTTCAAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCACTTATCGGATTTCGCTGCGTTCCTCCATCCGATGCAAGAGCTGAGATATCCGTTGCTGAGAGTTGTTCAGAAAAATACTTTGGAGAGGGGGAGAGCGAGTGCAGTGTCAAGGCCACATGCCGCCCCCCCTCCTTGAATACACGAAGGGCTATGCCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACCGCGCCGGACAAGGGTTAGATCGCCGGCAAGAAAGCTCCCAAGGAGGCTCTCCATGGCAACTCTAGGTCACCGCAATAGCAAGCGCCCATGGGTGACCAAAGTAAACCGATCCTCTGGATTATCTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGGAGGTTCACCTACGGAAACCTTGTTACG
    Z78494.1
    CTTAAACTCAGCGGGTTGCCTCACCTGACCTGGGATCGCAACCAAATGGCCATTCAACTGATCATGGGTCGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCTGATTAAACACAACAACCCGGGGCATTCGCACAACAGTTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCGAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTACGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGACAGGCGTTCCCTTGGCCTGTTGGCCTCGGGCCCAACTTGCGTTCAAAGACTCGATGTTAACGGGATTCTGCAATTCACACCTGAGATATCCGAAGTTGAGAGTGTTACCTCAATAATGGGGAGTTGGTCAAGGCAGTATGCCGTCCCCCTTCATTAATTATGGGTCATTTGATTTGGCGCATTACTGCGCCGGGCAAGGGTATAGATCGCCAGCAGGAAAGTTCCCAAGGAGGCCCGATGGTAAATCTAGGTCACTGCAACAGTAAATGTCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGCGACCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78493.1
    GGGTTGCCTCACCTGACCTGGGATCGCAACCAAATGGCCATTCAACTGATCATGGGTCGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCTGATTAAACACAACAACCCGGGGCATTCGCACAACAGTTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCGAAGAACAAGGGGCCAAACTCACATCCACACCGGATGACCAGTACGTTTGGACAGCCTCATTAAGGGAGAGATATGAGTCGCAATGTCCAGGCAGGCGTCCCCTTGGGCTGATGGCCTCGCGCGCAACTTGCGTTCAAAGATTCGATGGTTCACGGGATTCTGCAAGGCACACCATTTATCGAATCTCGCTGCGTTCCTTCATCGATGCATGAGCCGAGATATCCGTTGCTGAGAGTTGTTACCTAAATACTTGGGGGCTGGTCAAGGAAGAATGCCGCCCCCCTTCACCACTTATGTAAAATTTGACTTGGCGCAAAACTGCGCCGGGCAAGGGTATAGATCGCCAGCAGGAAAGCTCCCAGGGAGGCCCGATGGAAATTCTAGGTCACTGCAACAGAAAATGCCCATGGGTGACCAAAGTAAACCGATCCTCCAGATTAACTCGATCAATTATTATGCGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAATCCTTGTTACG
    Z78492.1
    TATTAAACTCAGCGGTGTGCCTCACCTGACTGGGATCGCAACCAAATGGCCAATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTCCAAAAGGGTCTTCTGATTAAGCTGGCCCATGTAACACAACAACCCGGGGCATTCGCACAACAGTTTTGTTGTAGCATAATTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCGTGCAGCTCTTAGACCCCCCGTCCCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATTTGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTCCACGGGATTCTGCAAATCACACCAGTATCGCATTTCGCTGGGTTCTACATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGAGTGGGTCAAGGCAGTACGCTCCCCCCCTTCACCAATTATGTGACATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTATAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATTCTAGGTCACTGCAACAGCAAATGCCCGTGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGCGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78491.1
    GCTTATCTCAGCGGGTTGCTCACCTGACTGGGATCGCAACAAATGGCCATTCAACTGATCATGGGTTGATGGGCCTCTAATGGGGTCCAAAGGGTCTTCTGATTATGCTGGTCCATGTAACACAACAACCCGGGGCATTCGCACAACATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAATATAATTTGGAGCGGGTCAAGGCAGCACGCCTCCCCCCAACACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTATAGATCGCCAGCAAGAAAGCTCCTAAGGAGGCCCGATGGCAAATCTAGGTCACTGCAACAGCAAATGCCCGTGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGCGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78490.1
    TCAGCTGGTTGCTCACCTGACTGGGGTCGCAACCAAATGGCCGTCAACTGATCATGGGTTGACGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAAATTATGCTGGCCCATCTAATACGACAACGCGGGGCATTCGCACAACACTTGTTGCAGAACAGTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCGTTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTTTGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAATAATTTGGGGAGGGTCGAGGCAGCATGCCGCCCCCTTCACCAATTATGTGCCATATGACTTGGCGCAAAACTGCGCCGGACTAGGGTTCAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTCGGTCACTGCAACAGCAAATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78489.1
    GCCTAAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCATGGGGTTCAAAAGGGTCTTCAGATATGCTGGGCCCATCTAATACGACAACTCGGGGCATTCGCACAACAGTTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAATCCCATGCAGCTCTTAGACCCCCCGTACCAATGAACAAGGGGCCAAACCCACATCCCCACCGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAATAATCTGGGGAGGGTCAAGGCAGCATGCCGCCCCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCTCGGGTTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTCGGTCACTGCAATAGCAAATGCCCTCGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78488.1
    AGCCGGTTGCTCACCTGACTGGGGTCGCAACAAATGTCCATCAACTGATCATGGGTTGATGGGCCTCCGATGGGGTTCAAAAGGGTCTTCAGATTATGATGGCCCATCTAATACGACAACCCGGGGGATCCGCACAACAGTTTGGTTGTGGCGTTGTTCGTCCACCTCTTGCCGTTTTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAATGAACAAGGGGGCAAACTCACATCCGCACCGGTTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATTTGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGGGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTCCATCCCGTTGATGAGAGTTGTTAGAAAATAATTAGGGGAGGGTCAAGGCACCATGCCGCCCCCTTCACCAATTATGTGTCGTATGACTTGGCGCAAAACTGCGCCGGGCTAGGGTTAAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTCGGTCACTGCAATAACAAATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTGCGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACAG
    Z78487.1
    TTACTCAGCGGTTGCTCACCTGACTGGGGTCGCAACAAGTGGCCCCAACTGATCATGGGTTGATGGGCCTCTAATGGGGTTCAAAAGGGTCTTTAGATTACGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAGTTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGCATTTAGGACCATCGATAGCCCATGCTGCTCTTAGACCCCCCGTACCGAAGAACAAGGGGCCGAACTCACATCCGCACCGGCTGAACAATATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAATAATTTGGGGAGGGTAAAGGCAGCATGCCGCCCCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCTAGGGTTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCCCGGTCACTGCAATAGCAAATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGAGTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78486.1
    TCAGCGGGTTGCCTCACCTGACCTGGGGTCGCTACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAGTTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCACGCAGCTCTTAGACCCCCCGTACCAATGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCGAGATATCCGTTGCTGAGAGTTGTTAAGAAATAATTTGGGGAGGGTCAGGCACCATGCCGCACCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCTAGGGTTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTCGGGCACTGCAATAGCAAATGTCCATGGGTGACCAAAGTAAACTGATCCTCCAGATAATCTCGATCAATTATTATGTGATCTCAACAATGATCCTCCCGCAGATTCACCTACGGAAACCTCGTGACG
    Z78485.1
    TACTAAATCAGGGGTTGCTCAGCTGACTGGGGTCGCAACACATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAGGGGTCTTTAGATTGTGCTGGCCCGTCTAATACGACAACCCGGGGCATTCGCACAACAGTTTTGTTGTAGCATAGTTCGTCCACCTCTTGGCGTATTTAGGACCATCAAAAGCCCACGCAGCTCTTAGACCCCCCGTACCAATGAACAAGGGCCCAAACTCACATCCGCACCGGCTGCACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGTTCACGGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAGAAATAATTTGGGGAGGGTCAAGGCACCATGCCGCACCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCCCCGGGCTAGGGTTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTCGGTCACTGCAATAGCAAATGCCCATGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGATATCTCAACAATGATCCATCCGCAGATTCACCTTCGGACACCAGGTTCAG
    Z78484.1
    AAAACTCAGGGAGTTGCTTCACCTGACTTGGGGTCGCAACCCAATGGACCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTGTCAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAGTTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCACGCAGCTCTTAGACCCCCCGTACCAATGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCCAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCGAGATATCCGTTGCTGAGAGTTGGGAGGGTAGGAAACATGCCGCCCCTTCACAATTATGTGTCATATGACTTGGGCCAAAACTGCGCCGGGCTAGGGTTTAGATCGCCAGCAAGAAAGCTCCCATGGAGGCTCGATGGTAAACTCTCGGTCACTGCAATAGCCAAATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCCCAGGTTCACCTACGGAAACCTTGTTACG
    Z78483.1
    TGCTAAACTCAGGGGTTGCTTCACTTGACCTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAGTTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCAACAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAATGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGATATGTGTGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATTGTTCACGGGATTCTGCAATTCACACCATTTATGATTTCGGTTTGGTCTTCATCCATGCAAGAGGCCAGATTTCCGTTGCTGAGAGTTGTTAAAAAATAATTTGGGGGAGGGTCAAGGCAGCATGCCGCCCCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCTAGGGTTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTCGGTCACTGCAATAGCAAATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78482.1
    TGCTAAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAGTTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCACGCAGCTCTTAGACCCCCCGTACCAATGAACAAGGGGCCAAACGCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGNTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGAAGAGGCGAGATATCCGTTGCNGAGAGTTGTTAAAAAATAATTTGGGGAGGGTCAAGGCAGCATGCCGCCCCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCTACGGTTTAGATCGCCAGCAAGAAAGCTCCCAGGAGGCTCGATGGCAAATCTCGGTCACTGCAGTAGA
    Z78481.1
    TCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAGTTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCAACAAAAGCCCACGCAGCTCTTAGACCCCCCGTACCAATGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTGTGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTCCGTTGCTGAGAGTTGTTAAAAAAATGATTTGGGGAGGGTCAAGGCACCATGCCGCCCCCTTCGCCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCTAGGGTTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTCGGTCACTGCAATAGCAAATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78480.1
    TCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTATACGACAACCCGGGGCATTCGCACAACATTTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGTACCAACAAAAGCCCACGCAGCTCTTAGACCCCCCGTACCAATGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTGTGGACAGCCTCATTGAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTATCCGTTGCTGAGAGTTGTTAAAAAATAATTTGGGGAGGGTCAAGGCAGCATGCCGCCCCCTTCGCCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCTAGGGTTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTCGGTCACTGCAATAGCAAATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78479.1
    ACTCAGAGGGTTGCCTCACCTGACCTGGGGTCGCATCCAAATGGCCGTCAACTGATCTTGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAGTTCTGTTGTCGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACGACAGCCTCATTAAGGGAGAGATATGTCTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAGACTCGATGGTTCACGGGATCTCTGTATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGGCGAGATATCCGTTGCTGAGAGTTGTTAAAAAATAATTTGGGGGGCAAGGCAAGGCAGCATCCCTCCCCTTCCAGTAATGTCTCATATAAATTGGCGCAAATCTGCGCCGGGCAAGGGTTTAGTTCGGCTCGATGGCAAATCTAGGTCACTGAAACAGCAAATGCCCATGGGTGACCAAAGTAACCTGATCCTCCTGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78478.1
    GCCTAAACTCAGCGGGTGCCTCATCTGACCTGGGTCGCAACCAAATGTCCTCAACTGATAATGGGTTGATGGGCCTCCAATGGGGGTCAAAAGGGTCTTTAGATTATGCTGACCCATCTAAACGACAACCGGGGCATTCGCACAACAGTTCTGTTGTCACATAGCTCGTCCACCTCTTGCCGTATTTAGGACCATCCAACCCATACAGCTCTTAGACCCCCGTACCAAAGGACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGCAGCCTCATTAAGGGAGAGATATGTCTCACAATTTCCAGGCAGGCGTCCCCTTGACCTGATGGCCTCGGACGCAACATGCGTTCAAAACTCGATGTTCAGGGATCTGCTATCACGCCATTATCAATATTTGGGGGAGCCATGGCAGCATGCCGCCCCTTCCAGTTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCGGCTCGATGGCAAACTCTAGGCCACTGCAACAGCAAATGCCCATGGGTGACCAAAGTAAACTGATCCCCCAGATTAACTCGATCAATTATTATGTGATCTCAACACTGATCCTTCCGGAGGTTCACCTACGGAAACCTTGTTACG
    Z78477.1
    GCATAAACTCAGCGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAGTTCTGTTGTCGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTGTGGACAGCCTCATTAAGGGAGAGATATGTCTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGCCGCAACTTGCGTTCCAAAGACTCGATGGTTCAGGGATTCTTGCAATTCACACCATTTATCGCATTTCGCTGCGTCCTTCCATCCGATGCAAGAGCCGAGATACTCGTTGCTGAGAGTTGGTTAACAAATATTTGGGGAGGGCAAGGCAGCATGCCGCCCCTTCCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCGGCTCGATGGCAAATCTAGGTCACTGCAACAGCAAATGCCCATGGGTGACCAAAGCAAACTGATCCTCCAGATTAACTCGATCAATTATCATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78476.1
    GGGTCGCAACCAAATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGGATTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCGCAACAGTTCTGTTGTCGCATAGTTCGTCCACCTCTTGCCGTATTTAGGTGGATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACTGTATGTATGGACAGCCTCATTAAGGGAGAGATATGTCTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCGGATGGCCTAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGGGAGGGCCAAGGCAGCATGCCGCCCTTCCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCGGCTCGATGGCAAATCTAGGTCACTACAACAGCAAATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATCATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78475.1
    ACCTGACCTGGGTCGCAACCAAATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTTGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCATTGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNAAGAGCCGAGATATCCGTTGCTGAGAGTTGTCAAAAAAATAATTTGGGGAGGGTCAAGGCAGCATGCCACCCCCTTCCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGACGGGTTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCCAGGTCACTGCAAGAGCAGATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78474.1
    AAGCTCAACGGGTTGCTTCACCTGACCTGGGGTCGCAACCAAATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAGAAGGGTCTGTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTTGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGTGTCTTCATCGATGCAAGAGGCGAGATATCCGTTGCTGAGAGTTGTCAAAAATATAATTTGGGGAGGGTCAAGGCAGGATGCCGCCCTTCCAATTATGTGTCATATGATTGGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCCAGGTCACTGAAAGAGCAGATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTACGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78473.1
    CCTGGGTCGCAACCAAATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGTACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTAGGCAAGAACAAGGGGCCAAACTCACATCCGCACCGACTGAACAGTATGTATGGACAGCCTCATTGAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAACAAAATAATTTGGGGAGGGTCAAGGCAGCATGCCGCCCCCTTCCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCGCCACCAAGAAAGCTCCCATGGAGGCTCGATGGCAAATCCAGGTCACTACAAGAGCAGATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACCCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78472.1
    GCTTAAACTCAGCGGTTGCCTCACCTGACTTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATTGGCCTCTAAAGGGGTTCAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCTGGGCATTCGAACTACAATTTTGTTGTAGCATACTTCGTCCACCTCTTGCCGTATTTAAGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCTGTGGCTGAACAGTATGTATAGACAGCCTCATTAAAGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTCCTCCATCGATGCAAGAGCCGAGATTCCGTTGCTGAGAGTTATCAAAAAATAATTTGGGGAGGGTCTTAGACTGCATGCCGCCCCCTTCCAATTATGTGTGAAATGACTTGGCGCAAGACTGCGCCGGGCAACGATTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCCAGGTCACTGCAAGAGCAGATGCCCATGGGTGACTAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78471.1
    GCTTAAACTCAGCGGGTTGCCTCACCTGACTTGGGGTCGCAACCAAATGGCCATCAACTGATCAGGGTTGATGGGCCTCTAATGGGGTTCAAAAGGGTCTTTAAATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACTACAATTTTGTTGTAGCATACTTTGTCCACCTCTTGCCGTATTTAGGACCAGCAAAAGCCCATGCAGCTTTTAGACCCCCCGTACTAAAGAACAAGGGGCCAAACTCACATCCGCACTGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAAGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGATGTTGAGAGTTGTCAAAAAATTAATTTGGGGAAGGGTCAAGGCAGCATGCCGCCCCCTTCCATTATGTGTCGGGGTGACTTGACGCAAGGCTGCGTCGGGCAACGATTTAGATCGCCAGCAAGACAGCTCCCAAGGAGGCTCGATGGCAAATCCAGGTCACTGTAAGAGCAGATGCCCATGGGTGACCAAAGTAAACCGATCCTCTAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78470.1
    AACTGATCAGGGTTGATGGGCCTCTAATGGGGTTCAAAAGGGTCTTTAAATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACTACAATTTTGTTGTAGCATACTTTGTCCACCTCTTGCCGTATTTAGGACCAGCAAAAGCCCATGCAGCTTTTAGACCCCCCGTAGTAAAGAACAAGGGGCCAAACTCACATCCGCACTGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAAGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTAAAAAAATAATTTGGGAAGGGTCAAGGCAGCATGCCGCCCCCTTCCAATTATGTGTCATATGACTTGACGCAAGGCTGCGTCGGGCAACGATTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCCAGGTCACTGTAAGAGCAGATGCCCATGGGTGACCAAAGTAAACTGATCCTCTAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78469.1
    AACTGATCATGGGTTGATGGGCCTCTAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCTGGGCATTCGCACTATCATTTTGTTGTAGCATACTTCGTCCACCTCTTGCCGTATTTAGGTACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCTGTGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGAAGAGCCGAGATATCCGTTGCTGAGAGTTGTCAAAAAAATAATTTGGGGAGGGTCTAGACCGCATGCCGCCCCCTTCCAATTATGTGTGATATGACTTGGCGCAAGACTGCGCCGGGCAACGAATTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCCAGGTCACTGCAAGAGCAGATGCCCATGGGTGACTAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78468.1
    AACCAAATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTATACGACAACCCGGGGCATTCACACAATAATTTTGTTGTAGCATACTTCGTCCACCTCTTGCCGTATTTAGGTCCATCAAAAGCCCAGCTGCTCTTAGACCCCCCGTACCGAAGAACAAGGGGCCAAACTCACATCCGCACTGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTAGCCTGATGGCCTCGGGCGCAACTTGAGTTCCGGAGCTGAGAGTCGTTAAAAAAATGATATGGGGAGGGTCAAGGCAGCATGCTGCCCCTCCCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTCAGATCGCCATCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCCAGGTCACTACAAGAGCAGATTCCCATGAGTGACCAAAGTAAACTGATCCTCCAGATTAACTCAATCAATTATTATGCGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78467.1
    TCAGCGGGTTGCCTCACCTGACCTGGGTCGCAAACATATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAGAGGGTCTTTAGATTATGCTGGCCCAGCTAATACGACAACCCGGGTCATTCGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTTGGACCACCAAAGGCACATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTACCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTCAAAAAATAATTTGGGGAGGGTCAAGGCAGCATGCCACCCCCTTCCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCCAGGTCACTGCAAGAGCAGATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78466.1
    GGGTTGCCTCACCTGACCTGGGTCGCAACCAAATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACGCAAGAACAAGGGGCCAAACTCACATCCGCACAGGCTGAACTGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGCAAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTACAAAAATAATTTGGGGAGGGTCAAGGCAGCATGCCGCCCCCTTCCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCGCCAGCAAGAAAGCTCCCATGAAGGCTCGATGGCAAATCCAGGTCACTACAAGAGCAGATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78465.1
    GCAACCAAATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGCATTCGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACTCGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTTCCCTTGGCCTGATGGCCTCGGGCGCAACTTTCGTTCAAAGACTCGATGGTTTACGGGATTCTTCAATTCACACCAATTATCGCATTTCGCTGCGTTCCTCATCGATTCAAGAACCGAGAAATCCGTTGCTGAGAGTTGTCAAAAAAATAATTTGGGGAGGGTCAAGGCAGCATGCCGCCCCCTTCCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCGCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCCAGGTCACTGCAAGAGCAGATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78464.1
    GCTTAAACTCAGCGGGTTGCTCACCTGACCTGGGGTCGCACCAAATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTATTACGAAAGCCCGGGGAATCCGAACAGAAATTTGGTTGAAGAATAGTCCGCCCACCTCTTGCCGTTTTTAGGACCACCAAAAGCCCATGAAGTTCTTAGACCCCCCGCCCCAAAGAAAAAGGGGCCAAACTCACATCCGAACCGGTTGAACAGTATGTTTGGACAGCCTCATTAAGGGAGAGATTTGATTCGAAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGAACTTGGCGTTCAAAGACTCGATGGTTCACGGGATTCTGAAATTCACACCATTTATCGGATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAGAAAATAATTTGGGGAGGGTCAAGGCAGCATGCCGCCCCCTTCCAATTATGTGTCATATGACTTGGCGCAAAACTGCTCCGGGCAAGGGTTTAGATCTCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAATCCAGGTCACTACAAGAGCAGATGCCCATGGGTGACCAAAGTAAACTGATCCTCCAGACTCAACCGATCAATTATTATGTGATCTCAACAATGACCCTTCCGCTCACCTACGGAAACCTTGTTACG
    Z78463.1
    GCTTAAACTCAGCGGGTTGCTCACCTGATCTGGGGTCGCAACCAAATGGCCGTCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTATGCTGGCCCATCTAATACGACAACCCGGGGTATTCGCACAGCAATTTTGTTGCAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGCCCCAAAGAAAAAGGGGCCAAACTCACATCCGCACCGGTTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAGAAAATAATTTGGGGAGGGTCAAGGCAACATGCCGCCCCCTCCCAATTATGTGTCATATGACTTGGCGCAAAACTCCCCCGGGCGAGGGTTTAGATCCCCAGCAAGAAAGCTCCCAAGGAGGCTCGATGGCAATCCAGGTCACTACAAGAGCAGATGCCCATGGTGACCAAAGTAAACTGATCCTCCAGACTCAACCGAACAATGATTATGTGATCTCAACAATGAACCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78462.1
    ATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAATCTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGCAAATGCCCATGGGTGACCAAAGTAAACAGATCCTGCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGGTCACATCCGGAGACCTCGTGACG
    Z78461.1
    TTAACTCAGCGGCTTGCTCACTGACTGGGTCGCAACAATGGCATCAACTGATCATGGGTTGATGGGCTCCAATGGGGTTCAAAGGGTCTTCAGATTACGGTGGCCCATCTTATACGACAACCCGGGGGCTTTGGACAACAATTTTGGTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGGTCTTAGACCCCCCGTACCAAAGAACAAGGGGGCAAACTCACATCCGCACCGGCTGAACAGGATGTATGGACAGGCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTTCCCTTGGCCTGATGGCCTCGGGCGCAACTTGGGTTCAAAGACTCGATGGTTCACGGGATTCTGGAATTCACACCATTTATCGGATTTCGGTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTGAAAAAATTATTTGGGGGAGGGTCAGGGAAGAATGCCACCCCCTCCACCAATTATGTGTCATTTGACTGGGGGAAAAATTGCGCCGGGTAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGAAATTCTAGGTCACTCCAGCAGCAACTGCCCATGGGTGACCAAAGTAAACAGATCCTCCAGATTACCTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGGAGGTTCACCTACGGAAACCTTGTTACG
    Z78460.1
    TAAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTGGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCCTCCATCCGATGCAAGAGCCGAGATATCCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGGGAGGGTCAGGGCAGGATGCCACCCCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGCAAATGCTCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGGAGGTTCACCTACGGAAACCTTGTTACG
    Z78459.1
    AAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCAATTTATCGCAATTTCGCTGCGTTCCTCCATCCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGGGAGGGTCAGGGCAGGATGCCACCCCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGCAAATGCTCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78458.1
    CAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTCTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATCTGCAATTCACACCATTTATCGCATTTCGCTGCGTCCTCCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAATAATTTGGGGAGGGTCAGGGAAGGATGCCACCCCCTTCACCAATTATGTGTCGTACGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGCAAATGCTCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78457.1
    CTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTTAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACGATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGGTCCTTCATCCGATGAAAGAGCCGAGATATCCGTTGTTGAGAGTTGTTAAAAAAATATTTGGGGGAGGGTCAGGGAAGAATGCCACCCCCTCCACCAATTATGTGCCATTTGACTGGGGGAAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGAAATTCTAGGTCACTACAACAGAAAATGCCCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGGAGGTTCACCTACGGAAACCTTGTTACG
    Z78456.1
    GCTAAACTCAGGGTTGCCTCACCTGACCTGGGTCGCAACCCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTCCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGTTCACGGGATTCTGCAATCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCCAGACCGGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGGGAGGGTCAGGGCAGGATGCCACCCCCTCCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATTCTAGGTCACTTCAACAGCAAATGCTCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78455.1
    TGCTAAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGAAGCTCTTAGACCCCCCGTGCCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGGTTCTTCAACCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAATTATTTGGGGGAGGGTAAGGGAAGGATGCCACCCCCTCCACCATTTTTGTGTAATTTGATTGGGGGAAAAACTGCGCCGGGAAAGGGTTTAGTTCTCGCCAAAAAGAAAGTTCCCAAGGGGGTTCGATGGAAATTCTAGGTCACTCCAAAAGAAATTGTTCATGGGTGACCAAAGTAAAAAGTTCCTCCAGATTAACTCGATCAATTTTTATGTGATCTCAAAAATGATCCTCCCGAAGGTCCACCTACGGAAACCTGGTTACG
    Z78454.1
    GTTGCCTCACCTGACCTGGGGTCGCATCCAAATGGCCATCAACTGATCATGNGGTTGATGGGCCTCCAATGGGTTCAAAAGGGGCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAGAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGTCTGATGGCCTCGGCCGCAACTTGCGTTCACAGACTCGATGGTTCACGGGATTCTGCAAATCACACCATTTATCGCATGAGTTCCGTTGATGAGAGTTGTTAACAATAGTGGGGGGAGGGTCAGGGCAGGATGCCACCCCCTTCACCGATTATGTGTCAAATGACTTGGAGAAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGTAAATCTAGGTCACTTCAACAGCAAATGCTCATGGGTGACCAAAGTAAACAGATCCTCCAGATTATCTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78453.1
    TGCTAAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCTGAGATATCCGGTTGCTGAGAGTTGTTAAAAAATAATTTGGGGGAGGGTCAGGGAAGGATGCCACCCCCTTCACCAGTTATGTGTCAAATGAGTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAACTCTAGGTCACTTCAACAGCAAATGCTCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78452.1
    TGCTAAACTCAGCGGTTGCCTCACCTGACCTGGGGTCGCATCCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACAGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTGTTGTAGCATAGTTCATCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATTCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGGGAGGGTCAGGGCAGGATGCCACCCCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAAAAAGAAAGCTCCCAAGGAGGCTTGATGGCAAATCTAGGTCACTTCAACAGCAAATGCTCATGGGTGACCAAAGTAAATAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78451.1
    GCTCAGCTGGTTGCTCACTGACTGGTGTCGCATCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCTATAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGGCTTCGCACAACAATTTTATTGTAGGGATAGTTCGTCCACCTCTTGCCGTTTTTAGGACCATCAAAATCCCATGCAGGTCTTAGACCCCCCGTACCAAAGAACAAGGGGGCAAACTCACATCCGCACCGGGTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATTTGACTCGCGATGCCCAGGGAGGGGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGGGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGATTTCGCTGGTTCTTCATCGATGAAAGAGCGAGATTTCCGTTGTTGAGAGTTGCAAAAAATACTTTGGGGAGGGTCAGGGCAGGATGCCCGCCCCCTACACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCACGGGTTTAGATCTCGCCAGCAAGAAAGCTCCCAGGGGGGCTCCATGGAAATCTAGGTCACTACAACAGCAAATGCCCATGGGTGACCAAAGTAAACAGATCCTCCAGATTATCTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGGAGGTACACCTACGGAAACCTTGTTACG
    Z78450.1
    CTCAGGGGTTGCTCAGCTGATCTGGGGTCGCAACACATGGTCATCAACTGATCATGGGTTGATGGTCCTCCAATGGGGTTCAACAGGGTCTACAGGTTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTGGCACAACAATTTGGTTGTAGCATAGTTCGTCCACCTCTTGCCAGTTTTGAGGACCATCATATGCCCATGCAGTTCTTAGACCCCCCGGACCAAAGAACAAGGGGCCACACTCACATCCGCACCGGATGAACAGTATGTATGGACAGCCTCGTTAGGGGAGAGATATGACTCGGAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGGGAGGGTCAGGGCAGGATGCCACCCCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCTCGGCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGTAAATGCTCATGGATGACCAAAGTAAACAGATCCTCCAGCTTAACTCGATCAATTATTATGTGATATCAGCAATGATCCTTCC
    Z78449.1
    GCATAAACTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTAGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTTCGCTGCGTTCTTCATCGATGCAAGAGCTGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGGGGAGGGTCAGGGCAGGATGCCACCCTTCACCAATTATGTGTCAAATGAGTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGCAAATGCTCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78448.1
    CCTCAGCGGGTTGCCTCACCTGACCTGGGGTCGCANCCAAATGGCCATCAACTGATCATGGGCTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTNTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGGTGCGTTCTTCCATCCGATGCAAGAGCTGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGGGAGGGTCAAGGGAAGGATGCCACCCCCCTTCACCAATTATGTGTCATACGACTTGGCGCAAAACTGCGCCGGGAAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAATTCTAGGTCACTTCAACAGCAAATGCTCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78447.1
    GCTTAAACTCAGCGGGTTGCCTCACCTGACTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTGGTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGGGTTCTTCATCGATGCAAGAGCCGAGATAGGATGCCACCCCCTTCACCACTTATGTGTCATTTGACTGGGGGCAAAATTGCGCCGGGTAAGGGTTTAGTTCTCGCCAACAAGAAAGCTCCCAGGGAGGCTCGATGGCAACTCTAGGTCACTTCAACAGGAAATGCCCATGGGTGACCAAAGTAAACAGATCCTCCAGATTACTCGATCAATTATTATGTGATCTCAACAATGATCCTCCCGCAGGTTCACCTACGGAATCCTTGTTACG
    Z78446.1
    GGGTTGCCTCACCTGACCTGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGAGCCTTTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGGCCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCGAAGAACAAGGGGCCAAACTCACATCTGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGTTCACGGGATTCTGCAATCACACCATTATCGCAATTTCGCTGGTCCTCCATCGATGAAGAGCCGAGTATCCGTTGCTGAGAGTTGTTAAAAAATTATTTGGGGAGGATGCCACCCCCTTCACCAATTATGTGTCATATGACTTGGCGCAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGGAAATGCTCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGGAGGTTCACCTACGGAAACCTTGTTACG
    Z78445.1
    ACAGCTGTTGCTCACCTGACCTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTCTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTGCACGGGATTCTGCAATTCACACCATTTATCGCATTTGCCGAAGACTGAGTCTCCGTTGCTGAGAGTTGTTAAAAAAATAGTTGGGGGAGGGTCAGGGAAGGATGCCACCCCCTTCACCAATTATGTGTCAAACGACTTGGAGCAAAACTGCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAATTCTAGGTCACTTCAACAGCAAATGCTCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACCCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78444.1
    AATGGCCATCAACTGACATGGGTTGATGGGCTCCAATGGGGTCCAAAAGGGTCTTCAGATATCGGTGGCCCATCTTATACGACAACCCGGGGCCTTTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCATCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGGTGGGTTCTTCATCGATGCAAGAGCTGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATATTTGGGGGAGGGTAAGGGAAGAATGCCACCCCCTCCACCAATTTTGTGTAATTTGATTGGGGGAAAAATTGAGGGTTTAGATATCGCCAACAAGGATGCTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGCAAATGCCCATGGGTGACCCAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGTTCACCCTACGGAAACCTTGTTACG
    Z78443.1
    CCTCACCTGACTGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTATGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGGTGTGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTATCGCATTTCGCTGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTGTGGCCCAAAACTCTGCCGGGCAAGGGTTTAGATATCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGCAAATGCCCATGGGTGACCCAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTACGGAAACCTTGTTACG
    Z78442.1
    ACTCAGCGGGTTGCTCAGCTGACTGGGGTCGCAACAAATGGTCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGGTTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTGGGACAACAATTTTGTTGTGGCATAGTTCGTCCACCTCTTGCCGTTTTGAGGACCATCAAATGCCCATGCAGCTCTTAGACCCCCGGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGATGAACAGTATGTTTGGATAGCCTCGTTAGGGGAGAGATATGATTCCCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGGGAGGGTCAGGGCAGGATGCCACCCCCTTCACCAATTATGTGTCATATGACTTGGCGCCAAACTCCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGCAAATGCCCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATGATCCTTCCGCAGGTTCACCTAC
    Z78441.1
    CTCAGCGCGTTGCTCAGCTGACTGGGGTCGCAACACATGGTCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGGTTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTGGGACAACAATTTGGTTGTAGCATAGTTCGTCCACCTCTTGCCGTTTTGAGGACCATCAAATGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGATGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGGGAGGGTCAGGGCAGGATGCCACCCCCTTCACCAATTATGTGTCATATGACTTGGCGCACAACTCCGCCGGGCAAGGGTTTAGATCTCGCCAACAAGAAAACTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGCAAATGCTCATGGGTGACCAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATATCGGCAATGACCTTCC
    Z78440.1
    TGCTAAACTCAGGGGTTGCCTCACCTGACCTGGGGTCGCAACCAAATGGCCATCAACTGATCATGGGTTGATGGGCCTCCAATGGGGTTCAAAAGGGTCTTCAGATTACGGTGGCCCATCTAATACGACAACCCGGGGCCTTTGCACAACAATTTTGCTGTAGCATAGTTCGTCCACCTCTTGCCGTATTTAGGACCATCAAAAGCCCATGCAGCTCTTAGACCCCCCGTAGCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCGTTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGGTTCACGGGATTCTGCAATTCACACCATTTATCGCATTTCGGTGGGTTCTTCAACGATGCAAGAGCTGAGATATCCGTTGCTGAGAGTTGTTAAAAAATTATTTGGGGGAGGGTAAGGGAAGATATGCCACCCCCTCCACCATTTTTGTGTAATTTGATTGGGGGAAAAACTGCGCCGGGTAAGGGTTTAGTTCTCGCCAACAAGAATGCTCCCAAGGAGGGTCGATGGAAATTCTAGGTCACTACAACAGAAATTGCCCATGGGTGACCAAAATATGCAGATCCTCCAGATTACCTCGATCAATTATTATGTGATCTCAACAATGATCCTCCCGGAGGTCCACCTACGGAAACCTTGTTACG
    Z78439.1
    GGCCCAACTAAACGGCCAACCGGGGCATTTGCACAACAATTTTGTTGTAGCATAGTTCGTCCACCTCTTTCCGTATTTAGGACCATCCAAAGCCCATGCAGCTCTTAGACCCCCCGTACCAAAGAACAAGGGGCCAAACTCACATCCGCACCGGCTGAACAGTATGTATGGACAGCCTCATTAAGGGAGAGATATGACTCGCAATGCCCAGGCAGGCGTGCCCTTGGCCTGATGGCCTCGGGCGCAACTTGCGTTCAAAGACTCGATGTTCACGGGATTCTGCAATTCACACCATTATCGCATTTCGCTGCGTTCTTCATCGATGCAAGAGCCGAGATATCCGTTGCTGAGAGTTGTTAAAAAAATAATTTGGGGAGGGTCAGGGCAGGATTACCACCCCCTTCACCAATTATGTGTCATATGACTTGGCGCCCAACTCCGCCGGGCAGGGGTTTAGATCTCGCCAACAAGAAAGCTCCCAAGGAGGCTCGATGGCAAATCTAGGTCACTTCAACAGCAAATGCCCATGGGTGACCAAAGTAAACAGATCCTCCAGATTAACTCGATCAATTATTATGTGATCTCAACAATG



```python
records = [
    rec.reverse_complement(id = "rc_" + rec.id, description = "reverse complement")
    for rec in SeqIO.parse("ls_orchid.fasta.txt", "fasta")
]
```


```python
len(records)
```




    94




```python
records = [
    rec.reverse_complement(id = "rc" + rec.id, description = "reverse complement")
    for rec in SeqIO.parse("ls_orchid.fasta.txt", "fasta")
    if len(rec) < 700
]
```


```python
len(records)
```




    18




```python
records = (
rec.reverse_complement(id = "rc_" + rec.id, description = "reverse complement")
for rec in SeqIO.parse("ls_orchid.fasta.txt", "fasta")
    if len(rec) < 700
)

SeqIO.write(records, "rev_comp.fasta", "fasta")
```




    18

## Multiple Sequence Alignment (Parts 1-5)
```python
from Bio import AlignIO
```


```python
alignment = AlignIO.read("PF05371_seed.sth.txt", "stockholm")
```


```python
print(alignment)
```

    Alignment with 7 rows and 52 columns
    AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRL...SKA COATB_BPIKE/30-81
    AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKL...SRA Q9T0Q8_BPIKE/1-52
    DGTSTATSYATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRL...SKA COATB_BPI22/32-83
    AEGDDP---AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA COATB_BPM13/24-72
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA COATB_BPZJ2/1-49
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA Q9T0Q9_BPFD/1-49
    FAADDATSQAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKL...SRA COATB_BPIF1/22-73



```python
alignment = AlignIO.read("PF05371_seed.sth.txt", "stockholm")
```


```python
print("Alignment length %i" % alignment.get_alignment_length())
```

    Alignment length 52



```python
for record in alignment:
    print("%s - %s" % (record.seq, record.id))
```

    AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA - COATB_BPIKE/30-81
    AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA - Q9T0Q8_BPIKE/1-52
    DGTSTATSYATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA - COATB_BPI22/32-83
    AEGDDP---AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA - COATB_BPM13/24-72
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA - COATB_BPZJ2/1-49
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA - Q9T0Q9_BPFD/1-49
    FAADDATSQAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA - COATB_BPIF1/22-73



```python
for record in alignment:
    if record.dbxrefs:
        print("%s %s" % (record.id, record.dbxrefs))
```

    COATB_BPIKE/30-81 ['PDB; 1ifl ; 1-52;']
    COATB_BPM13/24-72 ['PDB; 2cpb ; 1-49;', 'PDB; 2cps ; 1-49;']
    Q9T0Q9_BPFD/1-49 ['PDB; 1nh4 A; 1-49;']
    COATB_BPIF1/22-73 ['PDB; 1ifk ; 1-50;']



```python
for record in alignment:
    print(record)
```

    ID: COATB_BPIKE/30-81
    Name: COATB_BPIKE
    Description: COATB_BPIKE/30-81
    Database cross-references: PDB; 1ifl ; 1-52;
    Number of features: 0
    /accession=P03620.1
    /start=30
    /end=81
    Per letter annotation for: secondary_structure
    Seq('AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA')
    ID: Q9T0Q8_BPIKE/1-52
    Name: Q9T0Q8_BPIKE
    Description: Q9T0Q8_BPIKE/1-52
    Number of features: 0
    /accession=Q9T0Q8.1
    /start=1
    /end=52
    Seq('AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA')
    ID: COATB_BPI22/32-83
    Name: COATB_BPI22
    Description: COATB_BPI22/32-83
    Number of features: 0
    /accession=P15416.1
    /start=32
    /end=83
    Seq('DGTSTATSYATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA')
    ID: COATB_BPM13/24-72
    Name: COATB_BPM13
    Description: COATB_BPM13/24-72
    Database cross-references: PDB; 2cpb ; 1-49;, PDB; 2cps ; 1-49;
    Number of features: 0
    /accession=P69541.1
    /start=24
    /end=72
    Per letter annotation for: secondary_structure
    Seq('AEGDDP---AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA')
    ID: COATB_BPZJ2/1-49
    Name: COATB_BPZJ2
    Description: COATB_BPZJ2/1-49
    Number of features: 0
    /accession=P03618.1
    /start=1
    /end=49
    Seq('AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA')
    ID: Q9T0Q9_BPFD/1-49
    Name: Q9T0Q9_BPFD
    Description: Q9T0Q9_BPFD/1-49
    Database cross-references: PDB; 1nh4 A; 1-49;
    Number of features: 0
    /accession=Q9T0Q9.1
    /start=1
    /end=49
    Per letter annotation for: secondary_structure
    Seq('AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA')
    ID: COATB_BPIF1/22-73
    Name: COATB_BPIF1
    Description: COATB_BPIF1/22-73
    Database cross-references: PDB; 1ifk ; 1-50;
    Number of features: 0
    /accession=P03619.2
    /start=22
    /end=73
    Per letter annotation for: secondary_structure
    Seq('FAADDATSQAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA')



```python
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
```


```python
align1 = MultipleSeqAlignment([
    SeqRecord(Seq("ACTGCTAGCTAG"), id="Alpha"),
    SeqRecord(Seq("ACT-CTAGCTAG"), id="Beta"),
    SeqRecord(Seq("ACTGCTAGDTAG"), id="Gamma")
])

align2 = MultipleSeqAlignment([
    SeqRecord(Seq("GTCAGC-AG"), id="Delta"),
    SeqRecord(Seq("GACAGCTAG"), id="Epsilon"),
    SeqRecord(Seq("GTCAGCTAG"), id="Zeta"),
])

align3 = MultipleSeqAlignment([
    SeqRecord(Seq("ACTAGTACAGCTG"), id="Eta"),
    SeqRecord(Seq("ACTAGTACAGCT-"), id="Theta"),
    SeqRecord(Seq("-CTACTACAGGTG"), id="Iota"),
])
```


```python
my_alignments = [align1, align2, align3]
```


```python
my_alignments
```




    [<<class 'Bio.Align.MultipleSeqAlignment'> instance (3 records of length 12) at 7f6f904730d0>,
     <<class 'Bio.Align.MultipleSeqAlignment'> instance (3 records of length 9) at 7f6f904731d0>,
     <<class 'Bio.Align.MultipleSeqAlignment'> instance (3 records of length 13) at 7f6f904732d0>]




```python
print(my_alignments)
```

    [<<class 'Bio.Align.MultipleSeqAlignment'> instance (3 records of length 12) at 7f6f904730d0>, <<class 'Bio.Align.MultipleSeqAlignment'> instance (3 records of length 9) at 7f6f904731d0>, <<class 'Bio.Align.MultipleSeqAlignment'> instance (3 records of length 13) at 7f6f904732d0>]



```python
from Bio import AlignIO
AlignIO.write(my_alignments, "my_example.phy", "phylip")
```




    3




```python
alignments = AlignIO.parse("my_example.phy", "phylip")
```


```python
for alignment in alignments:
    print(alignment)
    print()
```

    Alignment with 3 rows and 12 columns
    ACTGCTAGCTAG Alpha
    ACT-CTAGCTAG Beta
    ACTGCTAGDTAG Gamma
    
    Alignment with 3 rows and 9 columns
    GTCAGC-AG Delta
    GACAGCTAG Epsilon
    GTCAGCTAG Zeta
    
    Alignment with 3 rows and 13 columns
    ACTAGTACAGCTG Eta
    ACTAGTACAGCT- Theta
    -CTACTACAGGTG Iota
    



```python
alignments = list(AlignIO.parse("my_example.phy", "phylip"))
```


```python
last_align = alignments[-1]
```


```python
print(last_align)
```

    Alignment with 3 rows and 13 columns
    ACTAGTACAGCTG Eta
    ACTAGTACAGCT- Theta
    -CTACTACAGGTG Iota



```python
first_align = alignments[0]
```


```python
print(first_align)
```

    Alignment with 3 rows and 12 columns
    ACTGCTAGCTAG Alpha
    ACT-CTAGCTAG Beta
    ACTGCTAGDTAG Gamma



```python
from Bio import AlignIO
```


```python
count = AlignIO.convert("PF05371_seed.sth.txt", "stockholm", "PF05371_seed.aln", "clustal")
```


```python
print("Converted %i alignments" % count)
```

    Converted 1 alignments



```python
alignments = AlignIO.parse("PF05371_seed.sth.txt", "stockholm")
```


```python
count = AlignIO.write(alignments, "PF05371_seed.aln", "clustal")
```


```python
print("Converted %i alignments" % count)
```

    Converted 1 alignments



```python
AlignIO.convert("PF05371_seed.sth.txt", "stockholm", "PF05371_seed.phy", "phylip")
```




    1




```python

```
AlignIO.convert("PF05371_seed.sth.txt", "stockholm", "PF05371_seed.phy", "phylip-relaxed")

```python
alignment = AlignIO.read("PF05371_seed.sth.txt", "stockholm")
name_mapping = {}
for i, record in enumerate(alignment):
    name_mapping[i] = record.id
    record.id = "seq%i" % i
```


```python
print(name_mapping)
```

    {0: 'COATB_BPIKE/30-81', 1: 'Q9T0Q8_BPIKE/1-52', 2: 'COATB_BPI22/32-83', 3: 'COATB_BPM13/24-72', 4: 'COATB_BPZJ2/1-49', 5: 'Q9T0Q9_BPFD/1-49', 6: 'COATB_BPIF1/22-73'}



```python
AlignIO.write([alignment], "PF05371_seed.phy", "phylip")
```




    1




```python
alignment = AlignIO.read("PF05371_seed.sth.txt", "stockholm")
```


```python
print("Number of rows: %i" %len(alignment))
```

    Number of rows: 7



```python
for record in alignment:
    print("%s - %s" % (record.seq, record.id))
```

    AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA - COATB_BPIKE/30-81
    AEPNAATNYATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA - Q9T0Q8_BPIKE/1-52
    DGTSTATSYATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA - COATB_BPI22/32-83
    AEGDDP---AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA - COATB_BPM13/24-72
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA - COATB_BPZJ2/1-49
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA - Q9T0Q9_BPFD/1-49
    FAADDATSQAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA - COATB_BPIF1/22-73



```python
print(alignment[3:7])
```

    Alignment with 4 rows and 52 columns
    AEGDDP---AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA COATB_BPM13/24-72
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA COATB_BPZJ2/1-49
    AEGDDP---AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKL...SKA Q9T0Q9_BPFD/1-49
    FAADDATSQAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKL...SRA COATB_BPIF1/22-73



```python
print(alignment[2, 6])
```

    T



```python
print(alignment[2].seq[6])
```

    T



```python
print(alignment[:, 6])
```

    TTT---T



```python
print(alignment[3:6, :6])
```

    Alignment with 3 rows and 6 columns
    AEGDDP COATB_BPM13/24-72
    AEGDDP COATB_BPZJ2/1-49
    AEGDDP Q9T0Q9_BPFD/1-49



```python
print(alignment[:, 6:9])
```

    Alignment with 7 rows and 3 columns
    TNY COATB_BPIKE/30-81
    TNY Q9T0Q8_BPIKE/1-52
    TSY COATB_BPI22/32-83
    --- COATB_BPM13/24-72
    --- COATB_BPZJ2/1-49
    --- Q9T0Q9_BPFD/1-49
    TSQ COATB_BPIF1/22-73



```python
print(alignment[:, 9:])
```

    Alignment with 7 rows and 43 columns
    ATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA COATB_BPIKE/30-81
    ATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA Q9T0Q8_BPIKE/1-52
    ATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA COATB_BPI22/32-83
    AKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA COATB_BPM13/24-72
    AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA COATB_BPZJ2/1-49
    AKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA Q9T0Q9_BPFD/1-49
    AKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA COATB_BPIF1/22-73



```python
edited = alignment[:, :6] + alignment[:, 9:]
```


```python
print(edited)
```

    Alignment with 7 rows and 49 columns
    AEPNAAATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA COATB_BPIKE/30-81
    AEPNAAATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA Q9T0Q8_BPIKE/1-52
    DGTSTAATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA COATB_BPI22/32-83
    AEGDDPAKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA COATB_BPM13/24-72
    AEGDDPAKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA COATB_BPZJ2/1-49
    AEGDDPAKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA Q9T0Q9_BPFD/1-49
    FAADDAAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA COATB_BPIF1/22-73



```python
edited.sort()
```


```python
print(edited)
```

    Alignment with 7 rows and 49 columns
    DGTSTAATEAMNSLKTQATDLIDQTWPVVTSVAVAGLAIRLFKKFSSKA COATB_BPI22/32-83
    FAADDAAKAAFDSLTAQATEMSGYAWALVVLVVGATVGIKLFKKFVSRA COATB_BPIF1/22-73
    AEPNAAATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIRLFKKFSSKA COATB_BPIKE/30-81
    AEGDDPAKAAFNSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA COATB_BPM13/24-72
    AEGDDPAKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFASKA COATB_BPZJ2/1-49
    AEPNAAATEAMDSLKTQAIDLISQTWPVVTTVVVAGLVIKLFKKFVSRA Q9T0Q8_BPIKE/1-52
    AEGDDPAKAAFDSLQASATEYIGYAWAMVVVIVGATIGIKLFKKFTSKA Q9T0Q9_BPFD/1-49



```python
from Bio.Seq import Seq
```


```python
from Bio.SeqRecord import SeqRecord
```


```python
from Bio.Align import MultipleSeqAlignment
```


```python
alignment = MultipleSeqAlignment(
    [
        SeqRecord(Seq("ACTCCTA"), id="seq1"),
        SeqRecord(Seq("AAT-CTA"), id="seq2"),
        SeqRecord(Seq("CCTACT-"), id="seq3"),
        SeqRecord(Seq("TCTCCTC"), id="seq4"),
    ]
)
```


```python
print(alignment)
```

    Alignment with 4 rows and 7 columns
    ACTCCTA seq1
    AAT-CTA seq2
    CCTACT- seq3
    TCTCCTC seq4



```python
substitutions = alignment.substitutions
```


```python
print(substitutions)
```

        A    C    T
    A 2.0  4.5  1.0
    C 4.5 10.0  0.5
    T 1.0  0.5 12.0
    



```python
m = substitutions.select("ATCG")
```


```python
print(m)
```

        A    T    C   G
    A 2.0  1.0  4.5 0.0
    T 1.0 12.0  0.5 0.0
    C 4.5  0.5 10.0 0.0
    G 0.0  0.0  0.0 0.0
    



```python
m = substitutions.select("ACTG")
```


```python
print(m)
```

        A    C    T   G
    A 2.0  4.5  1.0 0.0
    C 4.5 10.0  0.5 0.0
    T 1.0  0.5 12.0 0.0
    G 0.0  0.0  0.0 0.0
    



```python
import Bio.Align.Applications
```


```python
dir(Bio.Align.Applications)
```




    ['ClustalOmegaCommandline',
     'ClustalwCommandline',
     'DialignCommandline',
     'MSAProbsCommandline',
     'MafftCommandline',
     'MuscleCommandline',
     'PrankCommandline',
     'ProbconsCommandline',
     'TCoffeeCommandline',
     '_ClustalOmega',
     '_Clustalw',
     '_Dialign',
     '_MSAProbs',
     '_Mafft',
     '_Muscle',
     '_Prank',
     '_Probcons',
     '_TCoffee',
     '__all__',
     '__builtins__',
     '__cached__',
     '__doc__',
     '__file__',
     '__loader__',
     '__name__',
     '__package__',
     '__path__',
     '__spec__']




```python
from Bio.Align.Applications import ClustalwCommandline
```


```python
from Bio import AlignIO
```


```python
align = AlignIO.read("opuntia.aln.txt", "clustal")
```


```python
print(align)
```

    Alignment with 7 rows and 906 columns
    TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273285|gb|AF191659.1|AF191
    TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273284|gb|AF191658.1|AF191
    TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273287|gb|AF191661.1|AF191
    TATACATAAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273286|gb|AF191660.1|AF191
    TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273290|gb|AF191664.1|AF191
    TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273289|gb|AF191663.1|AF191
    TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273291|gb|AF191665.1|AF191



```python
from Bio import Phylo
```


```python
tree = Phylo.read("opuntia.dnd.txt", "newick")
```


```python
Phylo.draw_ascii(tree)
```

                                 _______________ gi|6273291|gb|AF191665.1|AF191665
      __________________________|
     |                          |   ______ gi|6273290|gb|AF191664.1|AF191664
     |                          |__|
     |                             |_____ gi|6273289|gb|AF191663.1|AF191663
     |
    _|_________________ gi|6273287|gb|AF191661.1|AF191661
     |
     |__________ gi|6273286|gb|AF191660.1|AF191660
     |
     |    __ gi|6273285|gb|AF191659.1|AF191659
     |___|
         | gi|6273284|gb|AF191658.1|AF191658
    



```python
from Bio import Align
```


```python
aligner = Align.PairwiseAligner()
```


```python
aligner = Align.PairwiseAligner(match_score = 1.0)
```


```python
target = "GAACT"
```


```python
query = "GAT"
```


```python
score = aligner.score(target, query)
```


```python
score
```




    3.0




```python
alignments = aligner.align(target, query)
```


```python
for alignment in alignments:
    print(alignment)
```

    target            0 GAACT 5
                      0 ||--| 5
    query             0 GA--T 3
    
    target            0 GAACT 5
                      0 |-|-| 5
    query             0 G-A-T 3
    



```python
aligner.mode = "local"
```


```python
target = "AGAACTC"
```


```python
query = "GAACT"
```


```python
score = aligner.score(target, query)
```


```python
score
```




    5.0




```python
alignments = aligner.align(target, query)
```


```python
for alignment in alignments:
    print(alignment)
```

    target            1 GAACT 6
                      0 ||||| 5
    query             0 GAACT 5
    



```python
print(aligner)
```

    Pairwise sequence aligner with parameters
      wildcard: None
      match_score: 1.000000
      mismatch_score: 0.000000
      target_internal_open_gap_score: 0.000000
      target_internal_extend_gap_score: 0.000000
      target_left_open_gap_score: 0.000000
      target_left_extend_gap_score: 0.000000
      target_right_open_gap_score: 0.000000
      target_right_extend_gap_score: 0.000000
      query_internal_open_gap_score: 0.000000
      query_internal_extend_gap_score: 0.000000
      query_left_open_gap_score: 0.000000
      query_left_extend_gap_score: 0.000000
      query_right_open_gap_score: 0.000000
      query_right_extend_gap_score: 0.000000
      mode: local
    



```python
aligner.algorithm
```




    'Smith-Waterman'




```python
aligner.epsilon
```




    1e-06




```python
from Bio import Align
```


```python
aligner = Align.PairwiseAligner()
```


```python
target = "GAACT"
```


```python
query = "GAT"
```


```python
alignments = aligner.align(target, query)
```


```python
alignment = alignments[0]
```


```python
alignment
```




    <Alignment object (2 rows x 5 columns) at 0x7f6f90452610>




```python
alignment.score
```




    3.0




```python
alignment.target
```




    'GAACT'




```python
alignment.query
```




    'GAT'




```python
print(alignment)
```

    target            0 GAACT 5
                      0 ||--| 5
    query             0 GA--T 3
    



```python
alignment.coordinates
```




    array([[0, 2, 4, 5],
           [0, 2, 2, 3]])




```python
len(alignment)
```




    2




```python
alignment.shape
```




    (2, 5)




```python
aligner.mode = "local"
```


```python
local_alignments = aligner.align("TGAACT", "GAC")
```


```python
local_alignment = local_alignments[0]
```


```python
print(local_alignment)
```

    target            1 GAAC 5
                      0 ||-| 4
    query             0 GA-C 3
    



```python
local_alignment.shape
```




    (2, 4)




```python
aligner.mode = "global"
```


```python
aligner = Align.PairwiseAligner(match = 1.0, mismatch_score = -10)
```


```python
print(aligner)
```

    Pairwise sequence aligner with parameters
      wildcard: None
      match_score: 1.000000
      mismatch_score: -10.000000
      target_internal_open_gap_score: 0.000000
      target_internal_extend_gap_score: 0.000000
      target_left_open_gap_score: 0.000000
      target_left_extend_gap_score: 0.000000
      target_right_open_gap_score: 0.000000
      target_right_extend_gap_score: 0.000000
      query_internal_open_gap_score: 0.000000
      query_internal_extend_gap_score: 0.000000
      query_left_open_gap_score: 0.000000
      query_left_extend_gap_score: 0.000000
      query_right_open_gap_score: 0.000000
      query_right_extend_gap_score: 0.000000
      mode: global
    



```python
alignments = aligner.align("AAACAAA", "AAAGAAA")
```


```python
len(alignments)
```




    2




```python
print(alignments[0])
```

    target            0 AAAC-AAA 7
                      0 |||--||| 8
    query             0 AAA-GAAA 7
    



```python
print(alignments[1])
```

    target            0 AAA-CAAA 7
                      0 |||--||| 8
    query             0 AAAG-AAA 7
    



```python
print(local_alignment)
```

    target            1 GAAC 5
                      0 ||-| 4
    query             0 GA-C 3
    



```python
local_alignment.sort()
```


```python
print(local_alignment)
```

    target            0 GA-C 3
                      0 ||-| 4
    query             1 GAAC 5
    



```python
from Bio import Align
```


```python
from Bio.Seq import reverse_complement
```


```python
target = "AAACCC"
```


```python
query = "AACC"
```


```python
aligner = Align.PairwiseAligner(mismatch_score= -1, internal_gap_score = -1)
```


```python
aligner.score(target, query)
```




    4.0




```python
aligner.score(target, reverse_complement(query))
```




    0.0




```python
aligner.score(target, reverse_complement(query), strand = "-")
```




    4.0




```python
aligner.score(target, query, strand = "-" )
```




    0.0




```python
alignments = aligner.align(target, query)
```


```python
len(alignments)
```




    1




```python
print(alignments[0])
```

    target            0 AAACCC 6
                      0 -||||- 6
    query             0 -AACC- 4
    



```python
print(alignments[0].format("bed"))
```

    target	1	5	query	4.0	+	1	5	0	1	4,	0,
    



```python
alignments = aligner.align(target, reverse_complement(query), strand = "-")
```


```python
print(alignments[0].format("bed"))
```

    target	1	5	query	4.0	-	1	5	0	1	4,	0,
    



```python
alignments = aligner.align(target, query, strand = "-")
```


```python
len(alignments)
```




    2




```python
print(alignments[0])
```

    target            0 AAACCC----  6
                      0 ---------- 10
    query             4 ------GGTT  0
    



```python
print(alignments[1])
```

    target            0 ----AAACCC  6
                      0 ---------- 10
    query             4 GGTT------  0
    



```python
aligner.left_gap_score = -0.5
```


```python
aligner.right_gap_score = -0.2
```


```python
aligner.score(target, query)
```




    3.3




```python
alignments = aligner.align(target, query)
```


```python
len(alignments)
```




    1




```python
print(alignments[0])
```

    target            0 AAACCC 6
                      0 -||||- 6
    query             0 -AACC- 4
    



```python
alignments = aligner.align(target, reverse_complement(query), strand = "-")
```


```python
print(alignments)
```

    <Bio.Align.PairwiseAlignments object at 0x7f74644af110>



```python
aligner.score(target, reverse_complement(query), strand = "-")
```




    3.3




```python
print(alignments[0])
```

    target            0 AAACCC 6
                      0 -||||- 6
    query             4 -AACC- 0
    



```python
aligner.score(target, query, strand = "+")
```




    3.3


## Blast + Challenge
### Blast
```python
from Bio.Blast import NCBIWWW
```


```python
NCBIWWW.email = "paytonmcalphin@yahoo.com"
```


```python
result_handle = NCBIWWW.qblast("blastn", "nt", "8332116")
```


```python
from Bio import SeqIO
```


```python
record = SeqIO.read("m_cold.fasta.txt", format = "fasta")
```


```python
print(record)
```

    ID: gi|8332116|gb|BE037100.1|BE037100
    Name: gi|8332116|gb|BE037100.1|BE037100
    Description: gi|8332116|gb|BE037100.1|BE037100 MP14H09 MP Mesembryanthemum crystallinum cDNA 5' similar to cold acclimation protein, mRNA sequence
    Number of features: 0
    Seq('CACTAGTACTCGAGCGTNCTGCACCAATTCGGCACGAGCAAGTGACTACGTTNT...TTC')



```python
result_handle = NCBIWWW.qblast("blastn", "nt", record.seq)
```


```python
with open("m_cold.fasta", "w") as out_handle:
    out_handle.write(result_handle.read())
result_handle.close()
```


```python
from Bio.Blast import NCBIXML
```


```python
result_handle = open("my_blast.xml")
```


```python
blast_record = NCBIXML.read(result_handle)
```


```python
E_VALUE_THRESH = 0.04
```


```python
for alignment in blast_record.alignments:
    for hsp in alignment.hsps:
        if hsp.expect < E_VALUE_THRESH:
            print("****ALIGNMENT****")
            print("sequence:", alignment.title)
            print("length:", alignment.length)
            print("e value:", hsp.expect)
            print(hsp.query[0:75] + "...")
            print(hsp.match[0:75] + "...")
            print(hsp.sbjct[0:75] + "...")
```

    ****ALIGNMENT****
    sequence: gi|262205317|ref|NR_030195.1| Homo sapiens microRNA 520b (MIR520B), microRNA
    length: 61
    e value: 4.91307e-23
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171311|ref|NR_035856.1| Pan troglodytes microRNA mir-520b (MIR520B), microRNA
    length: 60
    e value: 1.71483e-22
    CCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    CCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|270133242|ref|NR_032573.1| Macaca mulatta microRNA mir-519a (MIR519A), microRNA
    length: 85
    e value: 2.54503e-20
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||||| ||||||||||||||||| |||||||||||||||||||||||||||||||||||...
    CCCTCTAGAGGGAAGCGCTTTCTGTGGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171322|ref|NR_035857.1| Pan troglodytes microRNA mir-520c (MIR520C), microRNA
    length: 86
    e value: 8.88303e-20
    CCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    |||||| |||||||| ||||||||||||||||||||||||||||||||||||||||||||...
    CCTCTAGAGGGAAGCACTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171322|ref|NR_035857.1| Pan troglodytes microRNA mir-520c (MIR520C), microRNA
    length: 86
    e value: 3.31332e-06
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||||| | |||||| |||||| || | | | | || ||||||||||||| | ||||||...
    CCCTCTAAAAGGAAGCACTTTCTTTTCTTTCAGACAACAGAAAGTGCTTCCCTCTAGAGG...
    ****ALIGNMENT****
    sequence: gi|301171267|ref|NR_035851.1| Pan troglodytes microRNA mir-519b (MIR519B), microRNA
    length: 80
    e value: 8.88303e-20
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||||| ||||||||||||||||||||||||||||||||||||||| ||||||||||||...
    CCCTCTAGAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCATCCTTTTAGAGG...
    ****ALIGNMENT****
    sequence: gi|262205330|ref|NR_030198.1| Homo sapiens microRNA 520c (MIR520C), microRNA
    length: 87
    e value: 8.88303e-20
    CCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    |||||| |||||||| ||||||||||||||||||||||||||||||||||||||||||||...
    CCTCTAGAGGGAAGCACTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205330|ref|NR_030198.1| Homo sapiens microRNA 520c (MIR520C), microRNA
    length: 87
    e value: 3.31332e-06
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||||| | |||||| |||||| || | | | | || ||||||||||||| | ||||||...
    CCCTCTAAAAGGAAGCACTTTCTTTTCTTTCAGACAACAGAAAGTGCTTCCCTCTAGAGG...
    ****ALIGNMENT****
    sequence: gi|262205302|ref|NR_030191.1| Homo sapiens microRNA 519b (MIR519B), microRNA
    length: 81
    e value: 8.88303e-20
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||||| ||||||||||||||||||||||||||||||||||||||| ||||||||||||...
    CCCTCTAGAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCATCCTTTTAGAGG...
    ****ALIGNMENT****
    sequence: gi|301171259|ref|NR_035850.1| Pan troglodytes microRNA mir-519a (MIR519A), microRNA
    length: 86
    e value: 3.10048e-19
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    |||||||||||||||||||||||||||||||||||||| |||||||| |||||||||||...
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAGGAAAGTGCATCCTTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|262205451|ref|NR_030222.1| Homo sapiens microRNA 519a-2 (MIR519A2), microRNA
    length: 87
    e value: 3.10048e-19
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    |||||||||||||||||||||||||||||||||||||| |||||||| |||||||||||...
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAGGAAAGTGCATCCTTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|301171447|ref|NR_035871.1| Pan troglodytes microRNA mir-526b (MIR526B), microRNA
    length: 82
    e value: 3.77716e-18
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||||  |||||||| |||||||||||||||||||||||||||||||||||||||||||...
    CCCTCTTGAGGGAAGCACTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ****ALIGNMENT****
    sequence: gi|301171447|ref|NR_035871.1| Pan troglodytes microRNA mir-526b (MIR526B), microRNA
    length: 82
    e value: 0.000140886
    CCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    |||||| | |||||| |||||| || | | | | || ||||||||||||| |  ||||||...
    CCTCTAAAAGGAAGCACTTTCTTTTCTTTCAGACAACAGAAAGTGCTTCCCTCAAGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171276|ref|NR_035852.1| Pan troglodytes microRNA mir-519c (MIR519C), microRNA
    length: 86
    e value: 3.77716e-18
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||||| ||||||||||||||||||||||||||||||||||||||| || |||||||||...
    CCCTCTAGAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCATCTTTTTAGAGG...
    ****ALIGNMENT****
    sequence: gi|262205290|ref|NR_030188.1| Homo sapiens microRNA 519c (MIR519C), microRNA
    length: 87
    e value: 3.77716e-18
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||||| ||||||||||||||||||||||||||||||||||||||| || |||||||||...
    CCCTCTAGAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCATCTTTTTAGAGG...
    ****ALIGNMENT****
    sequence: gi|301171354|ref|NR_035860.1| Pan troglodytes microRNA mir-520f (MIR520F), microRNA
    length: 86
    e value: 1.31836e-17
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||||| ||||||||||||||||| ||| |||||||||| ||||||||||||||||||||...
    CCCTCTAAAGGGAAGCGCTTTCTGTGGTCAGAAAGAAAAGCAAGTGCTTCCTTTTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205281|ref|NR_030186.1| Homo sapiens microRNA 520f (MIR520F), microRNA
    length: 87
    e value: 1.31836e-17
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||||| ||||||||||||||||| ||| |||||||||| ||||||||||||||||||||...
    CCCTCTAAAGGGAAGCGCTTTCTGTGGTCAGAAAGAAAAGCAAGTGCTTCCTTTTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205298|ref|NR_030190.1| Homo sapiens microRNA 526b (MIR526B), microRNA
    length: 83
    e value: 4.60152e-17
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||||  |||||||| |||||||||||||||||||| ||||||||||||||||||||||...
    CCCTCTTGAGGGAAGCACTTTCTGTTGTCTGAAAGAAGAGAAAGTGCTTCCTTTTAGAGG...
    ****ALIGNMENT****
    sequence: gi|262205298|ref|NR_030190.1| Homo sapiens microRNA 526b (MIR526B), microRNA
    length: 83
    e value: 0.000140886
    CCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    |||||| | |||||| |||||| || | | | | || ||||||||||||| |  ||||||...
    CCTCTAAAAGGAAGCACTTTCTCTTCTTTCAGACAACAGAAAGTGCTTCCCTCAAGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171394|ref|NR_035865.1| Pan troglodytes microRNA mir-522 (MIR522), microRNA
    length: 86
    e value: 1.60609e-16
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||| ||||||||||||||||||||||||||||||||||| || |||| |||||||...
    CCCTCTAGAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAATGGTTCCCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|262205429|ref|NR_030218.1| Homo sapiens microRNA 519a-1 (MIR519A1), microRNA
    length: 85
    e value: 1.60609e-16
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||| |||||||||||||||||||||||||||||| |||||||| |||||||||||...
    CTCTAGAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAGGAAAGTGCATCCTTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|262205423|ref|NR_030217.1| Homo sapiens microRNA 522 (MIR522), microRNA
    length: 87
    e value: 1.60609e-16
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||| ||||||||||||||||||||||||||||||||||| || |||| |||||||...
    CCCTCTAGAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAATGGTTCCCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|301171401|ref|NR_035866.1| Pan troglodytes microRNA mir-523 (MIR523), microRNA
    length: 79
    e value: 5.6058e-16
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||||| |||||||||||||||||||||||||||||||||| | |||||| | |||||||...
    CCCTCTAGAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAACGCGCTTCCCTATAGAGGG...
    ****ALIGNMENT****
    sequence: gi|270133247|ref|NR_032574.1| Macaca mulatta microRNA mir-519b (MIR519B), microRNA
    length: 81
    e value: 5.6058e-16
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||||  |||||||||||||||||||||||||||||||||| |||| ||| |||||||||...
    CCCTCTGGAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAACGTGCATCCCTTTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205309|ref|NR_030193.1| Homo sapiens microRNA 523 (MIR523), microRNA
    length: 87
    e value: 5.6058e-16
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||||| |||||||||||||||||||||||||||||||||| | |||||| | |||||||...
    CCCTCTAGAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAACGCGCTTCCCTATAGAGGG...
    ****ALIGNMENT****
    sequence: gi|270132717|ref|NR_032716.1| Macaca mulatta microRNA mir-526b (MIR526B), microRNA
    length: 83
    e value: 1.95662e-15
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||||  |||||||| ||||||||||||||||  |||||||||||||||||||||||||...
    CCCTCTTGAGGGAAGCACTTTCTGTTGTCTGAATAAAAAGAAAGTGCTTCCTTTTAGAGG...
    ****ALIGNMENT****
    sequence: gi|270132717|ref|NR_032716.1| Macaca mulatta microRNA mir-526b (MIR526B), microRNA
    length: 83
    e value: 0.00171634
    CCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    |||||| | |||||| |||||| ||   | | | || ||||||||||||| |  ||||||...
    CCTCTAAAAGGAAGCACTTTCTTTTTATTCAGACAACAGAAAGTGCTTCCCTCAAGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171437|ref|NR_035870.1| Pan troglodytes microRNA mir-526a-2 (MIR526A-2), microRNA
    length: 68
    e value: 6.82927e-15
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||||| |||||||| |||||||||||||||||||||||||  ||| |||||| ||||||...
    CCCTCTAGAGGGAAGCACTTTCTGTTGTCTGAAAGAAAAGAACATGCATCCTTTCAGAGGG...
    ****ALIGNMENT****
    sequence: gi|270133306|ref|NR_032587.1| Macaca mulatta microRNA mir-523a (MIR523A), microRNA
    length: 87
    e value: 6.82927e-15
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||||| |||||||| |||||||||||||| |||||||||| | |||||| |||||||||...
    CCCTCTAGAGGGAAGCACTTTCTGTTGTCTGGAAGAAAAGAATGCGCTTCCCTTTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|270133306|ref|NR_032587.1| Macaca mulatta microRNA mir-523a (MIR523A), microRNA
    length: 87
    e value: 9.49283e-07
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||||| |||||||||| |||| || |   | | || ||||||||||||| | |||||||...
    CCCTCTAAAGGGAAGCGCATTCTTTTCTTCCAGACAACAGAAAGTGCTTCCCTCTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171428|ref|NR_035869.1| Pan troglodytes microRNA mir-526a-1 (MIR526A-1), microRNA
    length: 84
    e value: 2.38365e-14
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||||| |||||||| ||||||||||  |||||||| |||||| |||||||||||||||...
    CCCTCTAGAGGGAAGCACTTTCTGTTGCTTGAAAGAAGAGAAAGCGCTTCCTTTTAGAGG...
    ****ALIGNMENT****
    sequence: gi|301171428|ref|NR_035869.1| Pan troglodytes microRNA mir-526a-1 (MIR526A-1), microRNA
    length: 84
    e value: 2.71974e-07
    CCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    |||||| | ||||||||||||| || | | ||  || ||||||||||||| | |||||||...
    CCTCTAAAAGGAAGCGCTTTCTCTTCTTTCAAGCAACAGAAAGTGCTTCCCTCTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171211|ref|NR_035845.1| Pan troglodytes microRNA mir-518b (MIR518B), microRNA
    length: 82
    e value: 2.38365e-14
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||| | ||||||||||||||||||||||||||||||| |||| ||| || ||||||||...
    CCCTCCAGAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAACAAAGCGCTCCCCTTTAGAGG...
    ****ALIGNMENT****
    sequence: gi|301171153|ref|NR_035838.1| Pan troglodytes microRNA mir-516a-2 (MIR516A-2), microRNA
    length: 89
    e value: 2.38365e-14
    GAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||| ||||||||||||||||||||||||||||||||||||| ||||||...
    GAAGCACTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTCAGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171153|ref|NR_035838.1| Pan troglodytes microRNA mir-516a-2 (MIR516A-2), microRNA
    length: 89
    e value: 0.0209094
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTT...
    ||||||  | |||||| |||||| || | | | | || |||||||||||| ||...
    CCCTCTGAAAGGAAGCACTTTCTTTTCTTTCAGACAACAGAAAGTGCTTCTTT...
    ****ALIGNMENT****
    sequence: gi|301171146|ref|NR_035837.1| Pan troglodytes microRNA mir-516a-1 (MIR516A-1), microRNA
    length: 89
    e value: 2.38365e-14
    GAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||| ||||||||||||||||||||||||||||||||||||| ||||||...
    GAAGCACTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTCAGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171146|ref|NR_035837.1| Pan troglodytes microRNA mir-516a-1 (MIR516A-1), microRNA
    length: 89
    e value: 0.0209094
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTT...
    ||||||  | |||||| |||||| || | | | | || |||||||||||| ||...
    CCCTCTGAAAGGAAGCACTTTCTTTTCTTTCAGACAACAGAAAGTGCTTCTTT...
    ****ALIGNMENT****
    sequence: gi|270133254|ref|NR_032575.1| Macaca mulatta microRNA mir-519c (MIR519C), microRNA
    length: 87
    e value: 2.38365e-14
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||||| | |||||| ||||||||||| |||||||||||||||||| || |||||||||...
    CCCTCTAGAAGGAAGCACTTTCTGTTGTTTGAAAGAAAAGAAAGTGCATCATTTTAGAGG...
    ****ALIGNMENT****
    sequence: gi|270133254|ref|NR_032575.1| Macaca mulatta microRNA mir-519c (MIR519C), microRNA
    length: 87
    e value: 3.31332e-06
    CCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    |||||| |  || || |||||| || | | ||| || ||||||||||||||| |||||||...
    CCTCTAAAATGATGCACTTTCTTTTCTTTCAAACAACAGAAAGTGCTTCCTTCTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205445|ref|NR_030221.1| Homo sapiens microRNA 516a-2 (MIR516A2), microRNA
    length: 90
    e value: 2.38365e-14
    GAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||| ||||||||||||||||||||||||||||||||||||| ||||||...
    GAAGCACTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTCAGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205445|ref|NR_030221.1| Homo sapiens microRNA 516a-2 (MIR516A2), microRNA
    length: 90
    e value: 0.0209094
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTT...
    ||||||  | |||||| |||||| || | | | | || |||||||||||| ||...
    CCCTCTGAAAGGAAGCACTTTCTTTTCTTTCAGACAACAGAAAGTGCTTCTTT...
    ****ALIGNMENT****
    sequence: gi|262205441|ref|NR_030220.1| Homo sapiens microRNA 516a-1 (MIR516A1), microRNA
    length: 90
    e value: 2.38365e-14
    GAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||| ||||||||||||||||||||||||||||||||||||| ||||||...
    GAAGCACTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTCAGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205441|ref|NR_030220.1| Homo sapiens microRNA 516a-1 (MIR516A1), microRNA
    length: 90
    e value: 0.0209094
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTT...
    ||||||  | |||||| |||||| || | | | | || |||||||||||| ||...
    CCCTCTGAAAGGAAGCACTTTCTTTTCTTTCAGACAACAGAAAGTGCTTCTTT...
    ****ALIGNMENT****
    sequence: gi|262205325|ref|NR_030197.1| Homo sapiens microRNA 526a-1 (MIR526A1), microRNA
    length: 85
    e value: 2.38365e-14
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||||| |||||||| ||||||||||  |||||||| |||||| |||||||||||||||...
    CCCTCTAGAGGGAAGCACTTTCTGTTGCTTGAAAGAAGAGAAAGCGCTTCCTTTTAGAGG...
    ****ALIGNMENT****
    sequence: gi|262205325|ref|NR_030197.1| Homo sapiens microRNA 526a-1 (MIR526A1), microRNA
    length: 85
    e value: 2.71974e-07
    CCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    |||||| | ||||||||||||| || | | ||  || ||||||||||||| | |||||||...
    CCTCTAAAAGGAAGCGCTTTCTCTTCTTTCAAGCAACAGAAAGTGCTTCCCTCTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205321|ref|NR_030196.1| Homo sapiens microRNA 518b (MIR518B), microRNA
    length: 83
    e value: 2.38365e-14
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||| | ||||||||||||||||||||||||||||||| |||| ||| || ||||||||...
    CCCTCCAGAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAACAAAGCGCTCCCCTTTAGAGG...
    ****ALIGNMENT****
    sequence: gi|301171236|ref|NR_035848.1| Pan troglodytes microRNA mir-518e (MIR518E), microRNA
    length: 87
    e value: 8.31975e-14
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||| ||||||||||||||||||| || ||||||||||||| |||||| || ||||...
    CCCTCTAGAGGGAAGCGCTTTCTGTTGGCTAAAAGAAAAGAAAGCGCTTCCCTTCAGAG...
    ****ALIGNMENT****
    sequence: gi|301171236|ref|NR_035848.1| Pan troglodytes microRNA mir-518e (MIR518E), microRNA
    length: 87
    e value: 0.000491741
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||  ||||||||||||||| || | | |   || |||||| |||||| | |||||||...
    CTCTGAAGGGAAGCGCTTTCTTTTCTTTTAGCCAACAGAAAGCGCTTCCCTCTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205383|ref|NR_030209.1| Homo sapiens microRNA 518e (MIR518E), microRNA
    length: 88
    e value: 8.31975e-14
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||| ||||||||||||||||||| || ||||||||||||| |||||| || ||||...
    CCCTCTAGAGGGAAGCGCTTTCTGTTGGCTAAAAGAAAAGAAAGCGCTTCCCTTCAGAG...
    ****ALIGNMENT****
    sequence: gi|262205383|ref|NR_030209.1| Homo sapiens microRNA 518e (MIR518E), microRNA
    length: 88
    e value: 0.000491741
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||  ||||||||||||||| || | | |   || |||||| |||||| | |||||||...
    CTCTGAAGGGAAGCGCTTTCTTTTCTTTTAGCCAACAGAAAGCGCTTCCCTCTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171342|ref|NR_035859.1| Pan troglodytes microRNA mir-520e (MIR520E), microRNA
    length: 86
    e value: 2.90388e-13
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||| | | ||||||  ||||||||||||||||| ||||||||||||||||||| |||||...
    CCCTCAAGATGGAAGCAGTTTCTGTTGTCTGAAAGGAAAGAAAGTGCTTCCTTTTTGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205265|ref|NR_030183.1| Homo sapiens microRNA 520e (MIR520E), microRNA
    length: 87
    e value: 2.90388e-13
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||| | | ||||||  ||||||||||||||||| ||||||||||||||||||| |||||...
    CCCTCAAGATGGAAGCAGTTTCTGTTGTCTGAAAGGAAAGAAAGTGCTTCCTTTTTGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171247|ref|NR_035849.1| Pan troglodytes microRNA mir-518f (MIR518F), microRNA
    length: 86
    e value: 1.01355e-12
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||||| |||||||| |||||| |||||| ||||||||||||| |||||  ||||||||...
    CCCTCTAGAGGGAAGCACTTTCTCTTGTCTAAAAGAAAAGAAAGCGCTTCTCTTTAGAGG...
    ****ALIGNMENT****
    sequence: gi|301171247|ref|NR_035849.1| Pan troglodytes microRNA mir-518f (MIR518F), microRNA
    length: 86
    e value: 2.71974e-07
    CCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    |||||| || |||||||||||| || | | | | || ||||||||||||| | |||||||...
    CCTCTAAAGAGAAGCGCTTTCTTTTCTTTTAGACAAGAGAAAGTGCTTCCCTCTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|270133223|ref|NR_032569.1| Macaca mulatta microRNA mir-518c (MIR518C), microRNA
    length: 101
    e value: 1.01355e-12
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||  ||||||||||||||||||||||||||||||| |||| |||||  |||||||...
    CTCTGGAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAACAAAGCGCTTCTCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|262205313|ref|NR_030194.1| Homo sapiens microRNA 518f (MIR518F), microRNA
    length: 87
    e value: 1.01355e-12
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||||| |||||||| |||||| |||||| ||||||||||||| |||||  ||||||||...
    CCCTCTAGAGGGAAGCACTTTCTCTTGTCTAAAAGAAAAGAAAGCGCTTCTCTTTAGAGG...
    ****ALIGNMENT****
    sequence: gi|262205313|ref|NR_030194.1| Homo sapiens microRNA 518f (MIR518F), microRNA
    length: 87
    e value: 2.71974e-07
    CCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    |||||| || |||||||||||| || | | | | || ||||||||||||| | |||||||...
    CCTCTAAAGAGAAGCGCTTTCTTTTCTTTTAGACAAGAGAAAGTGCTTCCCTCTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171225|ref|NR_035847.1| Pan troglodytes microRNA mir-518d (MIR518D), microRNA
    length: 86
    e value: 3.53765e-12
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||| |||||||| |||||||||||||||||||||  |||| |||||| ||| |||...
    CCCTCTAGAGGGAAGCACTTTCTGTTGTCTGAAAGAAACCAAAGCGCTTCCCTTTGGAG...
    ****ALIGNMENT****
    sequence: gi|301171225|ref|NR_035847.1| Pan troglodytes microRNA mir-518d (MIR518D), microRNA
    length: 86
    e value: 0.000491741
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||| | |||||||||||||   || | | | | || ||||||||||||| | |||||||...
    CTCCAAAGGGAAGCGCTTTGGTTTCTTTCAGACAACAGAAAGTGCTTCCCTCTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|270133232|ref|NR_032571.1| Macaca mulatta microRNA mir-518e (MIR518E), microRNA
    length: 87
    e value: 3.53765e-12
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||| ||||||||| |||||||  |||||||||||||||| || |||| |||||||...
    CCCTCTAGAGGGAAGCGATTTCTGTGATCTGAAAGAAAAGAAAATGGTTCCCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|262205392|ref|NR_030211.1| Homo sapiens microRNA 518d (MIR518D), microRNA
    length: 87
    e value: 3.53765e-12
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||| |||||||| |||||||||||||||||||||  |||| |||||| ||| |||...
    CCCTCTAGAGGGAAGCACTTTCTGTTGTCTGAAAGAAACCAAAGCGCTTCCCTTTGGAG...
    ****ALIGNMENT****
    sequence: gi|262205392|ref|NR_030211.1| Homo sapiens microRNA 518d (MIR518D), microRNA
    length: 87
    e value: 0.000491741
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||| | |||||||||||||   || | | | | || ||||||||||||| | |||||||...
    CTCCAAAGGGAAGCGCTTTGGTTTCTTTCAGACAACAGAAAGTGCTTCCCTCTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171300|ref|NR_035855.1| Pan troglodytes microRNA mir-520a (MIR520A), microRNA
    length: 84
    e value: 1.23476e-11
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGA...
    ||||| | |||||||  ||||||||||||||| |||||||||||||||||| ||| ||...
    CCCTCCAGAGGGAAGTACTTTCTGTTGTCTGAGAGAAAAGAAAGTGCTTCCCTTTGGA...
    ****ALIGNMENT****
    sequence: gi|301171161|ref|NR_035839.1| Pan troglodytes microRNA mir-516b-1 (MIR516B-1), microRNA
    length: 89
    e value: 1.23476e-11
    GAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||| ||||||||| | ||||||||||||||||||||||||| ||||||...
    GAAGCACTTTCTGTTTTGTGAAAGAAAAGAAAGTGCTTCCTTTCAGAGGG...
    ****ALIGNMENT****
    sequence: gi|270133259|ref|NR_032576.1| Macaca mulatta microRNA mir-519d (MIR519D), microRNA
    length: 91
    e value: 1.23476e-11
    CCCTCTACAGGGAAGCGCTTTCTGTT----GTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||| | |||||||| |||||||||    |||||| |||||| |||||||||||||||||||...
    CCCTCCAAAGGGAAGCACTTTCTGTTTGTTGTCTGAGAGAAAACAAAGTGCTTCCTTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|270133259|ref|NR_032576.1| Macaca mulatta microRNA mir-519d (MIR519D), microRNA
    length: 91
    e value: 0.00171634
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGA---AAGAAA-AGAAAGTGCTTCCTTTTAGAGGG...
    ||||| | |||||| |||| | || ||| |   || ||| ||||||||||||| ||| |||||...
    CTCTAAAAGGAAGCACTTTGTTTTCTCTCAGACAACAAACAGAAAGTGCTTCCCTTTGGAGGG...
    ****ALIGNMENT****
    sequence: gi|270133218|ref|NR_032568.1| Macaca mulatta microRNA mir-518b (MIR518B), microRNA
    length: 83
    e value: 1.23476e-11
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||| | |||||||| |||||||||||||||||||||  |||| ||| || ||||||||...
    CCCTCCAGAGGGAAGCACTTTCTGTTGTCTGAAAGAAAGCAAAGCGCTCCCCTTTAGAGG...
    ****ALIGNMENT****
    sequence: gi|270133218|ref|NR_032568.1| Macaca mulatta microRNA mir-518b (MIR518B), microRNA
    length: 83
    e value: 0.00171634
    CCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    |||||| |||| ||||||||   || | | | | || ||||||||||||| | | |||||...
    CCTCTAAAGGGGAGCGCTTTGCTTTCTTTCAGACAACAGAAAGTGCTTCCCTCTGGAGGG...
    ****ALIGNMENT****
    sequence: gi|270132725|ref|NR_032718.1| Macaca mulatta microRNA mir-516 (MIR516), microRNA
    length: 90
    e value: 1.23476e-11
    GAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||| ||||||||| | ||||||||||||||||||||||||| ||||||...
    GAAGCACTTTCTGTTTTGTGAAAGAAAAGAAAGTGCTTCCTTTCAGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205396|ref|NR_030212.1| Homo sapiens microRNA 516b-1 (MIR516B1), microRNA
    length: 90
    e value: 1.23476e-11
    GAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||| ||||||||| | ||||||||||||||||||||||||| ||||||...
    GAAGCACTTTCTGTTTTGTGAAAGAAAAGAAAGTGCTTCCTTTCAGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205294|ref|NR_030189.1| Homo sapiens microRNA 520a (MIR520A), microRNA
    length: 85
    e value: 1.23476e-11
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGA...
    ||||| | |||||||  ||||||||||||||| |||||||||||||||||| ||| ||...
    CCCTCCAGAGGGAAGTACTTTCTGTTGTCTGAGAGAAAAGAAAGTGCTTCCCTTTGGA...
    ****ALIGNMENT****
    sequence: gi|301171602|ref|NR_035634.1| Pan troglodytes microRNA mir-1283 (MIR1283), microRNA
    length: 86
    e value: 4.30974e-11
    CTCTACA--GGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    |||||||  || ||||||||||||||||| ||||||| |||||| |||||| ||| |||||...
    CTCTACAAAGGAAAGCGCTTTCTGTTGTCAGAAAGAAGAGAAAGCGCTTCCCTTTTGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171602|ref|NR_035634.1| Pan troglodytes microRNA mir-1283 (MIR1283), microRNA
    length: 86
    e value: 0.00599063
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGC-TTCCTTT-TAGAG...
    ||||| | ||||||||||||||| || | |   | || |||||| || ||||||| |||||...
    CCCTCAAAAGGGAAGCGCTTTCTCTTCTTTCTGACAACAGAAAGCGCTTTCCTTTGTAGAG...
    ****ALIGNMENT****
    sequence: gi|301171454|ref|NR_035872.1| Pan troglodytes microRNA mir-527 (MIR527), microRNA
    length: 86
    e value: 4.30974e-11
    AGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTT...
    |||||||| ||||||||||||| |||||||||||||||||||| |||...
    AGGGAAGCCCTTTCTGTTGTCTAAAAGAAAAGAAAGTGCTTCCCTTT...
    ****ALIGNMENT****
    sequence: gi|301171333|ref|NR_035858.1| Pan troglodytes microRNA mir-520d (MIR520D), microRNA
    length: 86
    e value: 4.30974e-11
    CTCTACA--GGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    |||||||  ||||||| ||||||||||||| |||||||||||||||||||  ||| | |||...
    CTCTACAAAGGGAAGCCCTTTCTGTTGTCTAAAAGAAAAGAAAGTGCTTCTCTTTGGTGGG...
    ****ALIGNMENT****
    sequence: gi|301171291|ref|NR_035854.1| Pan troglodytes microRNA mir-519e (MIR519E), microRNA
    length: 86
    e value: 4.30974e-11
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||| | | || ||| |||||||||||||||||||||| ||||||| |||||||||||...
    CTCCAAAAGGGAGCACTTTCTGTTGTCTGAAAGAAAACAAAGTGCCTCCTTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|301171219|ref|NR_035846.1| Pan troglodytes microRNA mir-518c (MIR518C), microRNA
    length: 100
    e value: 4.30974e-11
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||  |||||||| |||||||||||||||||||||| |||| |||||  |||||||...
    CTCTGGAGGGAAGCACTTTCTGTTGTCTGAAAGAAAACAAAGCGCTTCTCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|301171219|ref|NR_035846.1| Pan troglodytes microRNA mir-518c (MIR518C), microRNA
    length: 100
    e value: 0.00599063
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||| || |||||||||| | || | | | | || ||||||||||||| |  ||||...
    CTCTAAAGAGAAGCGCTTTGTTTTCTTTCAGACAACAGAAAGTGCTTCCCTCCAGAG...
    ****ALIGNMENT****
    sequence: gi|269847477|ref|NR_031696.1| Homo sapiens microRNA 1283-2 (MIR1283-2), microRNA
    length: 87
    e value: 4.30974e-11
    TCTACA--GGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||  || |||||||||||||||||||||||||||||||  |||||| ||| |||...
    TCTACAAAGGAAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAATCGCTTCCCTTTGGAG...
    ****ALIGNMENT****
    sequence: gi|262205435|ref|NR_030219.1| Homo sapiens microRNA 527 (MIR527), microRNA
    length: 85
    e value: 4.30974e-11
    AGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTT...
    |||||||| ||||||||||||| |||||||||||||||||||| |||...
    AGGGAAGCCCTTTCTGTTGTCTAAAAGAAAAGAAAGTGCTTCCCTTT...
    ****ALIGNMENT****
    sequence: gi|262205336|ref|NR_030199.1| Homo sapiens microRNA 518c (MIR518C), microRNA
    length: 101
    e value: 4.30974e-11
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||  |||||||| |||||||||||||||||||||| |||| |||||  |||||||...
    CTCTGGAGGGAAGCACTTTCTGTTGTCTGAAAGAAAACAAAGCGCTTCTCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|262205336|ref|NR_030199.1| Homo sapiens microRNA 518c (MIR518C), microRNA
    length: 101
    e value: 0.00599063
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||| || |||||||||| | || | | | | || ||||||||||||| |  ||||...
    CTCTAAAGAGAAGCGCTTTGTTTTCTTTCAGACAACAGAAAGTGCTTCCCTCCAGAG...
    ****ALIGNMENT****
    sequence: gi|270133264|ref|NR_032577.1| Macaca mulatta microRNA mir-520a (MIR520A), microRNA
    length: 85
    e value: 1.50425e-10
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGA...
    ||||| | |||||||   ||||||||||||||| ||||||||||||||||| ||| ||...
    CCCTCCAGAGGGAAGTATTTTCTGTTGTCTGAAGGAAAAGAAAGTGCTTCCCTTTGGA...
    ****ALIGNMENT****
    sequence: gi|269846946|ref|NR_031573.1| Homo sapiens microRNA 1283-1 (MIR1283-1), microRNA
    length: 87
    e value: 1.50425e-10
    TCTACA--GGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||||  || ||||||||||||||||| ||||||| |||||| |||||| ||| |||||...
    TCTACAAAGGAAAGCGCTTTCTGTTGTCAGAAAGAAGAGAAAGCGCTTCCCTTTTGAGGG...
    ****ALIGNMENT****
    sequence: gi|269846946|ref|NR_031573.1| Homo sapiens microRNA 1283-1 (MIR1283-1), microRNA
    length: 87
    e value: 0.0209094
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGC-TTCCTTT-TAGA...
    ||||| | ||||||||||||||| || | |   | || |||||| || ||||||| ||||...
    CCCTCAAAAGGGAAGCGCTTTCTCTTCTTTCTGACAACAGAAAGCGCTTTCCTTTGTAGA...
    ****ALIGNMENT****
    sequence: gi|262205358|ref|NR_030204.1| Homo sapiens microRNA 520d (MIR520D), microRNA
    length: 87
    e value: 1.50425e-10
    TCTACA--GGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||||  ||||||| ||||||||||||| |||||||||||||||||||  ||| | |||...
    TCTACAAAGGGAAGCCCTTTCTGTTGTCTAAAAGAAAAGAAAGTGCTTCTCTTTGGTGGG...
    ****ALIGNMENT****
    sequence: gi|301171206|ref|NR_035844.1| Pan troglodytes microRNA mir-518a (MIR518A), microRNA
    length: 86
    e value: 5.25034e-10
    AGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTT...
    |||||||| ||||||||||||| ||||||||||||| |||||| |||...
    AGGGAAGCCCTTTCTGTTGTCTAAAAGAAAAGAAAGCGCTTCCCTTT...
    ****ALIGNMENT****
    sequence: gi|301171206|ref|NR_035844.1| Pan troglodytes microRNA mir-518a (MIR518A), microRNA
    length: 86
    e value: 0.00599063
    AGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTT...
    ||||||||||||||| || | | | | || |||||| |||||| |||...
    AGGGAAGCGCTTTCTTTTCTTTTAGACAACAGAAAGGGCTTCCCTTT...
    ****ALIGNMENT****
    sequence: gi|301171169|ref|NR_035840.1| Pan troglodytes microRNA mir-516b-2 (MIR516B-2), microRNA
    length: 89
    e value: 5.25034e-10
    GAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||| |||| |||| | ||||||||||||||||||||||||| ||||||...
    GAAGCACTTTGTGTTTTGTGAAAGAAAAGAAAGTGCTTCCTTTCAGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205387|ref|NR_030210.1| Homo sapiens microRNA 518a-1 (MIR518A1), microRNA
    length: 85
    e value: 5.25034e-10
    AGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTT...
    |||||||| |||||||||||||||||||| |||||| |||||| |||...
    AGGGAAGCCCTTTCTGTTGTCTGAAAGAAGAGAAAGCGCTTCCCTTT...
    ****ALIGNMENT****
    sequence: gi|262205387|ref|NR_030210.1| Homo sapiens microRNA 518a-1 (MIR518A1), microRNA
    length: 85
    e value: 0.00599063
    AGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTT...
    ||||||||||||||| || | | | | || |||||| |||||| |||...
    AGGGAAGCGCTTTCTCTTCTTTCAGACAACAGAAAGGGCTTCCCTTT...
    ****ALIGNMENT****
    sequence: gi|301171381|ref|NR_035863.1| Pan troglodytes microRNA mir-521-1 (MIR521-1), microRNA
    length: 86
    e value: 1.83255e-09
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||| | |||||||  ||||||||||||| ||||||||||| |  ||||| |||||||...
    CCCTCCAAAGGGAAGAACTTTCTGTTGTCTAAAAGAAAAGAACGCACTTCCCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|270133201|ref|NR_032565.1| Macaca mulatta microRNA mir-517 (MIR517), microRNA
    length: 86
    e value: 1.83255e-09
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||| | |||||| || ||||| |||| ||||||||||  |||| |||||||||||...
    CCCTCTAGATGGAAGCACTGTCTGTGGTCTAAAAGAAAAGATCGTGCATCCTTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|262205417|ref|NR_030216.1| Homo sapiens microRNA 521-1 (MIR521-1), microRNA
    length: 87
    e value: 1.83255e-09
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||| | |||||||  ||||||||||||| ||||||||||| |  ||||| |||||||...
    CCCTCCAAAGGGAAGAACTTTCTGTTGTCTAAAAGAAAAGAACGCACTTCCCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|270133214|ref|NR_032567.1| Macaca mulatta microRNA mir-518a (MIR518A), microRNA
    length: 87
    e value: 6.39622e-09
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTC-CTTT...
    ||||  | |||||||| ||||||||||||| || |||||||||||||||| ||||...
    CCCTACAAAGGGAAGCCCTTTCTGTTGTCTAAACGAAAAGAAAGTGCTTCTCTTT...
    ****ALIGNMENT****
    sequence: gi|262205407|ref|NR_030214.1| Homo sapiens microRNA 517c (MIR517C), microRNA
    length: 95
    e value: 6.39622e-09
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||| | |||||| || ||||||||||  |||||||||  |||| |||||||||||...
    CCCTCTAGATGGAAGCACTGTCTGTTGTCT--AAGAAAAGATCGTGCATCCTTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|262205377|ref|NR_030208.1| Homo sapiens microRNA 526a-2 (MIR526A2), microRNA
    length: 65
    e value: 6.39622e-09
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||||| |||||||| ||||||||||    |||||||||||  ||| |||||| ||||||...
    CCCTCTAGAGGGAAGCACTTTCTGTTG----AAAGAAAAGAACATGCATCCTTTCAGAGGG...
    ****ALIGNMENT****
    sequence: gi|270133311|ref|NR_032588.1| Macaca mulatta microRNA mir-523b (MIR523B), microRNA
    length: 89
    e value: 2.2325e-08
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCT--GAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||| || |||||||||||||||| ||   |||||| || ||| |||||| |||||||...
    CCCTCTAGAGCGAAGCGCTTTCTGTTGGCTAGAAAAGAATAGGAAGCGCTTCCCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|270133297|ref|NR_032585.1| Macaca mulatta microRNA mir-521 (MIR521), microRNA
    length: 87
    e value: 2.2325e-08
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||| | |||||||  ||||||||||||| ||||||||||| |  ||||| ||| |||...
    CCCTCCAAAGGGAAGTACTTTCTGTTGTCTAAAAGAAAAGAACGCACTTCCCTTTGGAG...
    ****ALIGNMENT****
    sequence: gi|262205401|ref|NR_030213.1| Homo sapiens microRNA 518a-2 (MIR518A2), microRNA
    length: 87
    e value: 2.2325e-08
    AGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTT...
    |||||||| ||||||||||||| |||||| |||||| |||||| |||...
    AGGGAAGCCCTTTCTGTTGTCTAAAAGAAGAGAAAGCGCTTCCCTTT...
    ****ALIGNMENT****
    sequence: gi|262205401|ref|NR_030213.1| Homo sapiens microRNA 518a-2 (MIR518A2), microRNA
    length: 87
    e value: 0.00599063
    AGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTT...
    ||||||||||||||| || | | | | || |||||| |||||| |||...
    AGGGAAGCGCTTTCTCTTCTTTTAGACAACAGAAAGGGCTTCCCTTT...
    ****ALIGNMENT****
    sequence: gi|301171372|ref|NR_035862.1| Pan troglodytes microRNA mir-520h (MIR520H), microRNA
    length: 89
    e value: 7.79219e-08
    CCCTCTACAGGGAAGCGCTTTCTGTT----GTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||| ||| |||| |||||||||    |||||| | |||| ||||||||||| |||||||...
    CCCTCTAGAGG-AAGCACTTTCTGTTTGTTGTCTGAGAAAAAACAAAGTGCTTCCCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|301171372|ref|NR_035862.1| Pan troglodytes microRNA mir-520h (MIR520H), microRNA
    length: 89
    e value: 0.00171634
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGA---AAGAAA-AGAAAGTGCTTCCTTTTAGAGGG...
    ||||| |||||||| |||| | || ||| |   || ||| ||||||||||||| | |||||||...
    CTCTAAAGGGAAGCACTTTGTTTTTTCTCAGACAACAAACAGAAAGTGCTTCC-TCTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205412|ref|NR_030215.1| Homo sapiens microRNA 520h (MIR520H), microRNA
    length: 88
    e value: 7.79219e-08
    CCCTCTACAGGGAAGCGCTTTCTGTT----GTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||| ||| |||| |||||||||    |||||| | |||| ||||||||||| |||||||...
    CCCTCTAGAGG-AAGCACTTTCTGTTTGTTGTCTGAGAAAAAACAAAGTGCTTCCCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|262205412|ref|NR_030215.1| Homo sapiens microRNA 520h (MIR520H), microRNA
    length: 88
    e value: 0.00171634
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGA---AAGAAA-AGAAAGTGCTTCCTTTTAGAGGG...
    ||||| |||||||| |||| | || ||| |   || ||| ||||||||||||| | |||||||...
    CTCTAAAGGGAAGCACTTTGTTTTTTCTCAGACAACAAACAGAAAGTGCTTCC-TCTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205366|ref|NR_030206.1| Homo sapiens microRNA 520g (MIR520G), microRNA
    length: 90
    e value: 7.79219e-08
    CCCTCTACAGGGAAGCGCTTTCTGTT----GTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||| ||| |||| |||||||||    |||||| | |||| ||||||||||| |||||||...
    CCCTCTAGAGG-AAGCACTTTCTGTTTGTTGTCTGAGAAAAAACAAAGTGCTTCCCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|262205366|ref|NR_030206.1| Homo sapiens microRNA 520g (MIR520G), microRNA
    length: 90
    e value: 0.00171634
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGA---AAGAAA-AGAAAGTGCTTCCTTTTAGAGGG...
    ||||| |||||||| |||| | || ||| |   || ||| ||||||||||||| | |||||||...
    CTCTAAAGGGAAGCACTTTGTTTTTTCTCAGACAACAAACAGAAAGTGCTTCC-TCTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|270132721|ref|NR_032717.1| Macaca mulatta microRNA mir-524 (MIR524), microRNA
    length: 85
    e value: 2.71974e-07
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTT...
    ||||  | ||| |||| |||||| |||||| ||||||||||| |||||||| |||...
    CCCTACAAAGGCAAGCACTTTCTCTTGTCTAAAAGAAAAGAAGGTGCTTCCCTTT...
    ****ALIGNMENT****
    sequence: gi|262205276|ref|NR_030185.1| Homo sapiens microRNA 519e (MIR519E), microRNA
    length: 84
    e value: 9.49283e-07
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||| | | || ||| |||||||||   |||||||||| ||||||| |||||||||||...
    CTCCAAAAGGGAGCACTTTCTGTT---TGAAAGAAAACAAAGTGCCTCCTTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|301171363|ref|NR_035861.1| Pan troglodytes microRNA mir-520g (MIR520G), microRNA
    length: 89
    e value: 3.31332e-06
    CCCTCTACAGGGAAGCGCTTTCTGTT----GTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||| ||| | || |||||||||    |||||| | |||| ||||||||||| |||||||...
    CCCTCTAGAGG-ACGCACTTTCTGTTTGTTGTCTGAGAAAAAACAAAGTGCTTCCCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|270133196|ref|NR_032564.1| Macaca mulatta microRNA mir-516a-2 (MIR516A-2), microRNA
    length: 89
    e value: 3.31332e-06
    GAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||| ||||||||||||| ||||||||| ||||| ||||||   |||||...
    GAAGCACTTTCTGTTGTCT-AAAGAAAAGGAAGTGTTTCCTTCCCGAGGG...
    ****ALIGNMENT****
    sequence: gi|270133191|ref|NR_032563.1| Macaca mulatta microRNA mir-516a-1 (MIR516A-1), microRNA
    length: 89
    e value: 3.31332e-06
    GAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||| ||||||||||||| ||||||||| ||||| ||||||   |||||...
    GAAGCACTTTCTGTTGTCT-AAAGAAAAGGAAGTGTTTCCTTCCCGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171387|ref|NR_035864.1| Pan troglodytes microRNA mir-521-2 (MIR521-2), microRNA
    length: 86
    e value: 1.15646e-05
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||| | |||||||   ||||| |||||| ||||||||||| |  ||||| |||||||...
    CTCCAAAGGGAAGAATTTTCTCTTGTCTAAAAGAAAAGAACGCACTTCCCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|262205353|ref|NR_030203.1| Homo sapiens microRNA 521-2 (MIR521-2), microRNA
    length: 87
    e value: 1.15646e-05
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||| | |||||||   ||||| |||||| ||||||||||| |  ||||| |||||||...
    CTCCAAAGGGAAGAATTTTCTCTTGTCTAAAAGAAAAGAACGCACTTCCCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|270133237|ref|NR_032572.1| Macaca mulatta microRNA mir-518f (MIR518F), microRNA
    length: 87
    e value: 0.000140886
    CCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    |||||  ||||||||||||||| |  | | | | || | ||||||||||| | |||||||...
    CCTCTGAAGGGAAGCGCTTTCTTTCCTTTCACACAAGATAAAGTGCTTCCCTCTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171283|ref|NR_035853.1| Pan troglodytes microRNA mir-519d (MIR519D), microRNA
    length: 87
    e value: 0.00171634
    CCCTCTACAGGGAAGCGCTTTCTG-TTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||| | |||||||||||||||| |||| |      ||| ||||||| ||| |||||||...
    CCCTCCAAAGGGAAGCGCTTTCTGTTTGTTTTCTCTCAAACAAAGTGCCTCCCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|270133288|ref|NR_032583.1| Macaca mulatta microRNA mir-520g (MIR520G), microRNA
    length: 90
    e value: 0.00171634
    CCCTCTACAGGGAAGCGCTTTCTGTT----GTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||| ||  |||| |||||||||    ||||||   |||| ||||||||||| || ||||...
    CCCTCTAGAGA-AAGCACTTTCTGTTTGTTGTCTGAGGAAAAACAAAGTGCTTCCCTTCAGAG...
    ****ALIGNMENT****
    sequence: gi|262205371|ref|NR_030207.1| Homo sapiens microRNA 516b-2 (MIR516B2), microRNA
    length: 85
    e value: 0.00171634
    GAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||| |||| |||| | ||||||     |||||||||||||| ||||||...
    GAAGCACTTTGTGTTTTGTGAAAG-----AAAGTGCTTCCTTTCAGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205349|ref|NR_030202.1| Homo sapiens microRNA 519d (MIR519D), microRNA
    length: 88
    e value: 0.00171634
    CCCTCTACAGGGAAGCGCTTTCTG-TTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||| | |||||||||||||||| |||| |      ||| ||||||| ||| |||||||...
    CCCTCCAAAGGGAAGCGCTTTCTGTTTGTTTTCTCTTAAACAAAGTGCCTCCCTTTAGAG...

### Challenge 1
```python
from Bio.Blast import NCBIWWW
```


```python
NCBIWWW.email = "paytonmcalphin@yahoo.com"
```


```python
result_handle = NCBIWWW.qblast("blastn", "nt", "12567")
```


```python
from Bio import SeqIO
```


```python
record = SeqIO.read("CDK4.fasta.txt", format = "fasta")
```


```python
print(record)
```

    ID: NC_000076.7:126899404-126903157
    Name: NC_000076.7:126899404-126903157
    Description: NC_000076.7:126899404-126903157 Cdk4 [organism=Mus musculus] [GeneID=12567] [chromosome=10]
    Number of features: 0
    Seq('AACGTCCGGCGCCCGCCCCGCCCCCCGGGTTTCCGCGCGCCTCTCTGGCAGCTG...GTA')



```python
result_handle = NCBIWWW.qblast("blastn", "nt", record.seq)
```


```python
with open("CDK4.fasta.txt", "w") as out_handle:
    out_handle.write(result_handle.read())
result_handle.close()
```


```python
from Bio.Blast import NCBIXML
```


```python
result_handle = open("my_blast.xml")
```


```python
blast_record = NCBIXML.read(result_handle)
```


```python
E_VALUE_THRESH = 0.04
```


```python
for alignment in blast_record.alignments:
    for hsp in alignment.hsps:
        if hsp.expect < E_VALUE_THRESH:
            print("****ALIGNMENT****")
            print("sequence:", alignment.title)
            print("length:", alignment.length)
            print("e value:", hsp.expect)
            print(hsp.query[0:75] + "...")
            print(hsp.match[0:75] + "...")
            print(hsp.sbjct[0:75] + "...")
```

    ****ALIGNMENT****
    sequence: gi|262205317|ref|NR_030195.1| Homo sapiens microRNA 520b (MIR520B), microRNA
    length: 61
    e value: 4.91307e-23
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171311|ref|NR_035856.1| Pan troglodytes microRNA mir-520b (MIR520B), microRNA
    length: 60
    e value: 1.71483e-22
    CCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||...
    CCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|270133242|ref|NR_032573.1| Macaca mulatta microRNA mir-519a (MIR519A), microRNA
    length: 85
    e value: 2.54503e-20
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||||| ||||||||||||||||| |||||||||||||||||||||||||||||||||||...
    CCCTCTAGAGGGAAGCGCTTTCTGTGGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171322|ref|NR_035857.1| Pan troglodytes microRNA mir-520c (MIR520C), microRNA
    length: 86
    e value: 8.88303e-20
    CCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    |||||| |||||||| ||||||||||||||||||||||||||||||||||||||||||||...
    CCTCTAGAGGGAAGCACTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171322|ref|NR_035857.1| Pan troglodytes microRNA mir-520c (MIR520C), microRNA
    length: 86
    e value: 3.31332e-06
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||||| | |||||| |||||| || | | | | || ||||||||||||| | ||||||...
    CCCTCTAAAAGGAAGCACTTTCTTTTCTTTCAGACAACAGAAAGTGCTTCCCTCTAGAGG...
    ****ALIGNMENT****
    sequence: gi|301171267|ref|NR_035851.1| Pan troglodytes microRNA mir-519b (MIR519B), microRNA
    length: 80
    e value: 8.88303e-20
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||||| ||||||||||||||||||||||||||||||||||||||| ||||||||||||...
    CCCTCTAGAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCATCCTTTTAGAGG...
    ****ALIGNMENT****
    sequence: gi|262205330|ref|NR_030198.1| Homo sapiens microRNA 520c (MIR520C), microRNA
    length: 87
    e value: 8.88303e-20
    CCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    |||||| |||||||| ||||||||||||||||||||||||||||||||||||||||||||...
    CCTCTAGAGGGAAGCACTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205330|ref|NR_030198.1| Homo sapiens microRNA 520c (MIR520C), microRNA
    length: 87
    e value: 3.31332e-06
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||||| | |||||| |||||| || | | | | || ||||||||||||| | ||||||...
    CCCTCTAAAAGGAAGCACTTTCTTTTCTTTCAGACAACAGAAAGTGCTTCCCTCTAGAGG...
    ****ALIGNMENT****
    sequence: gi|262205302|ref|NR_030191.1| Homo sapiens microRNA 519b (MIR519B), microRNA
    length: 81
    e value: 8.88303e-20
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||||| ||||||||||||||||||||||||||||||||||||||| ||||||||||||...
    CCCTCTAGAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCATCCTTTTAGAGG...
    ****ALIGNMENT****
    sequence: gi|301171259|ref|NR_035850.1| Pan troglodytes microRNA mir-519a (MIR519A), microRNA
    length: 86
    e value: 3.10048e-19
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    |||||||||||||||||||||||||||||||||||||| |||||||| |||||||||||...
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAGGAAAGTGCATCCTTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|262205451|ref|NR_030222.1| Homo sapiens microRNA 519a-2 (MIR519A2), microRNA
    length: 87
    e value: 3.10048e-19
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    |||||||||||||||||||||||||||||||||||||| |||||||| |||||||||||...
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAGGAAAGTGCATCCTTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|301171447|ref|NR_035871.1| Pan troglodytes microRNA mir-526b (MIR526B), microRNA
    length: 82
    e value: 3.77716e-18
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||||  |||||||| |||||||||||||||||||||||||||||||||||||||||||...
    CCCTCTTGAGGGAAGCACTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ****ALIGNMENT****
    sequence: gi|301171447|ref|NR_035871.1| Pan troglodytes microRNA mir-526b (MIR526B), microRNA
    length: 82
    e value: 0.000140886
    CCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    |||||| | |||||| |||||| || | | | | || ||||||||||||| |  ||||||...
    CCTCTAAAAGGAAGCACTTTCTTTTCTTTCAGACAACAGAAAGTGCTTCCCTCAAGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171276|ref|NR_035852.1| Pan troglodytes microRNA mir-519c (MIR519C), microRNA
    length: 86
    e value: 3.77716e-18
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||||| ||||||||||||||||||||||||||||||||||||||| || |||||||||...
    CCCTCTAGAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCATCTTTTTAGAGG...
    ****ALIGNMENT****
    sequence: gi|262205290|ref|NR_030188.1| Homo sapiens microRNA 519c (MIR519C), microRNA
    length: 87
    e value: 3.77716e-18
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||||| ||||||||||||||||||||||||||||||||||||||| || |||||||||...
    CCCTCTAGAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCATCTTTTTAGAGG...
    ****ALIGNMENT****
    sequence: gi|301171354|ref|NR_035860.1| Pan troglodytes microRNA mir-520f (MIR520F), microRNA
    length: 86
    e value: 1.31836e-17
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||||| ||||||||||||||||| ||| |||||||||| ||||||||||||||||||||...
    CCCTCTAAAGGGAAGCGCTTTCTGTGGTCAGAAAGAAAAGCAAGTGCTTCCTTTTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205281|ref|NR_030186.1| Homo sapiens microRNA 520f (MIR520F), microRNA
    length: 87
    e value: 1.31836e-17
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||||| ||||||||||||||||| ||| |||||||||| ||||||||||||||||||||...
    CCCTCTAAAGGGAAGCGCTTTCTGTGGTCAGAAAGAAAAGCAAGTGCTTCCTTTTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205298|ref|NR_030190.1| Homo sapiens microRNA 526b (MIR526B), microRNA
    length: 83
    e value: 4.60152e-17
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||||  |||||||| |||||||||||||||||||| ||||||||||||||||||||||...
    CCCTCTTGAGGGAAGCACTTTCTGTTGTCTGAAAGAAGAGAAAGTGCTTCCTTTTAGAGG...
    ****ALIGNMENT****
    sequence: gi|262205298|ref|NR_030190.1| Homo sapiens microRNA 526b (MIR526B), microRNA
    length: 83
    e value: 0.000140886
    CCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    |||||| | |||||| |||||| || | | | | || ||||||||||||| |  ||||||...
    CCTCTAAAAGGAAGCACTTTCTCTTCTTTCAGACAACAGAAAGTGCTTCCCTCAAGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171394|ref|NR_035865.1| Pan troglodytes microRNA mir-522 (MIR522), microRNA
    length: 86
    e value: 1.60609e-16
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||| ||||||||||||||||||||||||||||||||||| || |||| |||||||...
    CCCTCTAGAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAATGGTTCCCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|262205429|ref|NR_030218.1| Homo sapiens microRNA 519a-1 (MIR519A1), microRNA
    length: 85
    e value: 1.60609e-16
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||| |||||||||||||||||||||||||||||| |||||||| |||||||||||...
    CTCTAGAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAGGAAAGTGCATCCTTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|262205423|ref|NR_030217.1| Homo sapiens microRNA 522 (MIR522), microRNA
    length: 87
    e value: 1.60609e-16
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||| ||||||||||||||||||||||||||||||||||| || |||| |||||||...
    CCCTCTAGAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAATGGTTCCCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|301171401|ref|NR_035866.1| Pan troglodytes microRNA mir-523 (MIR523), microRNA
    length: 79
    e value: 5.6058e-16
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||||| |||||||||||||||||||||||||||||||||| | |||||| | |||||||...
    CCCTCTAGAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAACGCGCTTCCCTATAGAGGG...
    ****ALIGNMENT****
    sequence: gi|270133247|ref|NR_032574.1| Macaca mulatta microRNA mir-519b (MIR519B), microRNA
    length: 81
    e value: 5.6058e-16
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||||  |||||||||||||||||||||||||||||||||| |||| ||| |||||||||...
    CCCTCTGGAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAACGTGCATCCCTTTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205309|ref|NR_030193.1| Homo sapiens microRNA 523 (MIR523), microRNA
    length: 87
    e value: 5.6058e-16
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||||| |||||||||||||||||||||||||||||||||| | |||||| | |||||||...
    CCCTCTAGAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAACGCGCTTCCCTATAGAGGG...
    ****ALIGNMENT****
    sequence: gi|270132717|ref|NR_032716.1| Macaca mulatta microRNA mir-526b (MIR526B), microRNA
    length: 83
    e value: 1.95662e-15
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||||  |||||||| ||||||||||||||||  |||||||||||||||||||||||||...
    CCCTCTTGAGGGAAGCACTTTCTGTTGTCTGAATAAAAAGAAAGTGCTTCCTTTTAGAGG...
    ****ALIGNMENT****
    sequence: gi|270132717|ref|NR_032716.1| Macaca mulatta microRNA mir-526b (MIR526B), microRNA
    length: 83
    e value: 0.00171634
    CCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    |||||| | |||||| |||||| ||   | | | || ||||||||||||| |  ||||||...
    CCTCTAAAAGGAAGCACTTTCTTTTTATTCAGACAACAGAAAGTGCTTCCCTCAAGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171437|ref|NR_035870.1| Pan troglodytes microRNA mir-526a-2 (MIR526A-2), microRNA
    length: 68
    e value: 6.82927e-15
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||||| |||||||| |||||||||||||||||||||||||  ||| |||||| ||||||...
    CCCTCTAGAGGGAAGCACTTTCTGTTGTCTGAAAGAAAAGAACATGCATCCTTTCAGAGGG...
    ****ALIGNMENT****
    sequence: gi|270133306|ref|NR_032587.1| Macaca mulatta microRNA mir-523a (MIR523A), microRNA
    length: 87
    e value: 6.82927e-15
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||||| |||||||| |||||||||||||| |||||||||| | |||||| |||||||||...
    CCCTCTAGAGGGAAGCACTTTCTGTTGTCTGGAAGAAAAGAATGCGCTTCCCTTTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|270133306|ref|NR_032587.1| Macaca mulatta microRNA mir-523a (MIR523A), microRNA
    length: 87
    e value: 9.49283e-07
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||||| |||||||||| |||| || |   | | || ||||||||||||| | |||||||...
    CCCTCTAAAGGGAAGCGCATTCTTTTCTTCCAGACAACAGAAAGTGCTTCCCTCTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171428|ref|NR_035869.1| Pan troglodytes microRNA mir-526a-1 (MIR526A-1), microRNA
    length: 84
    e value: 2.38365e-14
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||||| |||||||| ||||||||||  |||||||| |||||| |||||||||||||||...
    CCCTCTAGAGGGAAGCACTTTCTGTTGCTTGAAAGAAGAGAAAGCGCTTCCTTTTAGAGG...
    ****ALIGNMENT****
    sequence: gi|301171428|ref|NR_035869.1| Pan troglodytes microRNA mir-526a-1 (MIR526A-1), microRNA
    length: 84
    e value: 2.71974e-07
    CCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    |||||| | ||||||||||||| || | | ||  || ||||||||||||| | |||||||...
    CCTCTAAAAGGAAGCGCTTTCTCTTCTTTCAAGCAACAGAAAGTGCTTCCCTCTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171211|ref|NR_035845.1| Pan troglodytes microRNA mir-518b (MIR518B), microRNA
    length: 82
    e value: 2.38365e-14
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||| | ||||||||||||||||||||||||||||||| |||| ||| || ||||||||...
    CCCTCCAGAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAACAAAGCGCTCCCCTTTAGAGG...
    ****ALIGNMENT****
    sequence: gi|301171153|ref|NR_035838.1| Pan troglodytes microRNA mir-516a-2 (MIR516A-2), microRNA
    length: 89
    e value: 2.38365e-14
    GAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||| ||||||||||||||||||||||||||||||||||||| ||||||...
    GAAGCACTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTCAGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171153|ref|NR_035838.1| Pan troglodytes microRNA mir-516a-2 (MIR516A-2), microRNA
    length: 89
    e value: 0.0209094
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTT...
    ||||||  | |||||| |||||| || | | | | || |||||||||||| ||...
    CCCTCTGAAAGGAAGCACTTTCTTTTCTTTCAGACAACAGAAAGTGCTTCTTT...
    ****ALIGNMENT****
    sequence: gi|301171146|ref|NR_035837.1| Pan troglodytes microRNA mir-516a-1 (MIR516A-1), microRNA
    length: 89
    e value: 2.38365e-14
    GAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||| ||||||||||||||||||||||||||||||||||||| ||||||...
    GAAGCACTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTCAGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171146|ref|NR_035837.1| Pan troglodytes microRNA mir-516a-1 (MIR516A-1), microRNA
    length: 89
    e value: 0.0209094
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTT...
    ||||||  | |||||| |||||| || | | | | || |||||||||||| ||...
    CCCTCTGAAAGGAAGCACTTTCTTTTCTTTCAGACAACAGAAAGTGCTTCTTT...
    ****ALIGNMENT****
    sequence: gi|270133254|ref|NR_032575.1| Macaca mulatta microRNA mir-519c (MIR519C), microRNA
    length: 87
    e value: 2.38365e-14
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||||| | |||||| ||||||||||| |||||||||||||||||| || |||||||||...
    CCCTCTAGAAGGAAGCACTTTCTGTTGTTTGAAAGAAAAGAAAGTGCATCATTTTAGAGG...
    ****ALIGNMENT****
    sequence: gi|270133254|ref|NR_032575.1| Macaca mulatta microRNA mir-519c (MIR519C), microRNA
    length: 87
    e value: 3.31332e-06
    CCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    |||||| |  || || |||||| || | | ||| || ||||||||||||||| |||||||...
    CCTCTAAAATGATGCACTTTCTTTTCTTTCAAACAACAGAAAGTGCTTCCTTCTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205445|ref|NR_030221.1| Homo sapiens microRNA 516a-2 (MIR516A2), microRNA
    length: 90
    e value: 2.38365e-14
    GAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||| ||||||||||||||||||||||||||||||||||||| ||||||...
    GAAGCACTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTCAGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205445|ref|NR_030221.1| Homo sapiens microRNA 516a-2 (MIR516A2), microRNA
    length: 90
    e value: 0.0209094
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTT...
    ||||||  | |||||| |||||| || | | | | || |||||||||||| ||...
    CCCTCTGAAAGGAAGCACTTTCTTTTCTTTCAGACAACAGAAAGTGCTTCTTT...
    ****ALIGNMENT****
    sequence: gi|262205441|ref|NR_030220.1| Homo sapiens microRNA 516a-1 (MIR516A1), microRNA
    length: 90
    e value: 2.38365e-14
    GAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||| ||||||||||||||||||||||||||||||||||||| ||||||...
    GAAGCACTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTCAGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205441|ref|NR_030220.1| Homo sapiens microRNA 516a-1 (MIR516A1), microRNA
    length: 90
    e value: 0.0209094
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTT...
    ||||||  | |||||| |||||| || | | | | || |||||||||||| ||...
    CCCTCTGAAAGGAAGCACTTTCTTTTCTTTCAGACAACAGAAAGTGCTTCTTT...
    ****ALIGNMENT****
    sequence: gi|262205325|ref|NR_030197.1| Homo sapiens microRNA 526a-1 (MIR526A1), microRNA
    length: 85
    e value: 2.38365e-14
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||||| |||||||| ||||||||||  |||||||| |||||| |||||||||||||||...
    CCCTCTAGAGGGAAGCACTTTCTGTTGCTTGAAAGAAGAGAAAGCGCTTCCTTTTAGAGG...
    ****ALIGNMENT****
    sequence: gi|262205325|ref|NR_030197.1| Homo sapiens microRNA 526a-1 (MIR526A1), microRNA
    length: 85
    e value: 2.71974e-07
    CCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    |||||| | ||||||||||||| || | | ||  || ||||||||||||| | |||||||...
    CCTCTAAAAGGAAGCGCTTTCTCTTCTTTCAAGCAACAGAAAGTGCTTCCCTCTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205321|ref|NR_030196.1| Homo sapiens microRNA 518b (MIR518B), microRNA
    length: 83
    e value: 2.38365e-14
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||| | ||||||||||||||||||||||||||||||| |||| ||| || ||||||||...
    CCCTCCAGAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAACAAAGCGCTCCCCTTTAGAGG...
    ****ALIGNMENT****
    sequence: gi|301171236|ref|NR_035848.1| Pan troglodytes microRNA mir-518e (MIR518E), microRNA
    length: 87
    e value: 8.31975e-14
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||| ||||||||||||||||||| || ||||||||||||| |||||| || ||||...
    CCCTCTAGAGGGAAGCGCTTTCTGTTGGCTAAAAGAAAAGAAAGCGCTTCCCTTCAGAG...
    ****ALIGNMENT****
    sequence: gi|301171236|ref|NR_035848.1| Pan troglodytes microRNA mir-518e (MIR518E), microRNA
    length: 87
    e value: 0.000491741
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||  ||||||||||||||| || | | |   || |||||| |||||| | |||||||...
    CTCTGAAGGGAAGCGCTTTCTTTTCTTTTAGCCAACAGAAAGCGCTTCCCTCTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205383|ref|NR_030209.1| Homo sapiens microRNA 518e (MIR518E), microRNA
    length: 88
    e value: 8.31975e-14
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||| ||||||||||||||||||| || ||||||||||||| |||||| || ||||...
    CCCTCTAGAGGGAAGCGCTTTCTGTTGGCTAAAAGAAAAGAAAGCGCTTCCCTTCAGAG...
    ****ALIGNMENT****
    sequence: gi|262205383|ref|NR_030209.1| Homo sapiens microRNA 518e (MIR518E), microRNA
    length: 88
    e value: 0.000491741
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||  ||||||||||||||| || | | |   || |||||| |||||| | |||||||...
    CTCTGAAGGGAAGCGCTTTCTTTTCTTTTAGCCAACAGAAAGCGCTTCCCTCTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171342|ref|NR_035859.1| Pan troglodytes microRNA mir-520e (MIR520E), microRNA
    length: 86
    e value: 2.90388e-13
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||| | | ||||||  ||||||||||||||||| ||||||||||||||||||| |||||...
    CCCTCAAGATGGAAGCAGTTTCTGTTGTCTGAAAGGAAAGAAAGTGCTTCCTTTTTGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205265|ref|NR_030183.1| Homo sapiens microRNA 520e (MIR520E), microRNA
    length: 87
    e value: 2.90388e-13
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||| | | ||||||  ||||||||||||||||| ||||||||||||||||||| |||||...
    CCCTCAAGATGGAAGCAGTTTCTGTTGTCTGAAAGGAAAGAAAGTGCTTCCTTTTTGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171247|ref|NR_035849.1| Pan troglodytes microRNA mir-518f (MIR518F), microRNA
    length: 86
    e value: 1.01355e-12
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||||| |||||||| |||||| |||||| ||||||||||||| |||||  ||||||||...
    CCCTCTAGAGGGAAGCACTTTCTCTTGTCTAAAAGAAAAGAAAGCGCTTCTCTTTAGAGG...
    ****ALIGNMENT****
    sequence: gi|301171247|ref|NR_035849.1| Pan troglodytes microRNA mir-518f (MIR518F), microRNA
    length: 86
    e value: 2.71974e-07
    CCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    |||||| || |||||||||||| || | | | | || ||||||||||||| | |||||||...
    CCTCTAAAGAGAAGCGCTTTCTTTTCTTTTAGACAAGAGAAAGTGCTTCCCTCTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|270133223|ref|NR_032569.1| Macaca mulatta microRNA mir-518c (MIR518C), microRNA
    length: 101
    e value: 1.01355e-12
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||  ||||||||||||||||||||||||||||||| |||| |||||  |||||||...
    CTCTGGAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAACAAAGCGCTTCTCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|262205313|ref|NR_030194.1| Homo sapiens microRNA 518f (MIR518F), microRNA
    length: 87
    e value: 1.01355e-12
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||||| |||||||| |||||| |||||| ||||||||||||| |||||  ||||||||...
    CCCTCTAGAGGGAAGCACTTTCTCTTGTCTAAAAGAAAAGAAAGCGCTTCTCTTTAGAGG...
    ****ALIGNMENT****
    sequence: gi|262205313|ref|NR_030194.1| Homo sapiens microRNA 518f (MIR518F), microRNA
    length: 87
    e value: 2.71974e-07
    CCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    |||||| || |||||||||||| || | | | | || ||||||||||||| | |||||||...
    CCTCTAAAGAGAAGCGCTTTCTTTTCTTTTAGACAAGAGAAAGTGCTTCCCTCTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171225|ref|NR_035847.1| Pan troglodytes microRNA mir-518d (MIR518D), microRNA
    length: 86
    e value: 3.53765e-12
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||| |||||||| |||||||||||||||||||||  |||| |||||| ||| |||...
    CCCTCTAGAGGGAAGCACTTTCTGTTGTCTGAAAGAAACCAAAGCGCTTCCCTTTGGAG...
    ****ALIGNMENT****
    sequence: gi|301171225|ref|NR_035847.1| Pan troglodytes microRNA mir-518d (MIR518D), microRNA
    length: 86
    e value: 0.000491741
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||| | |||||||||||||   || | | | | || ||||||||||||| | |||||||...
    CTCCAAAGGGAAGCGCTTTGGTTTCTTTCAGACAACAGAAAGTGCTTCCCTCTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|270133232|ref|NR_032571.1| Macaca mulatta microRNA mir-518e (MIR518E), microRNA
    length: 87
    e value: 3.53765e-12
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||| ||||||||| |||||||  |||||||||||||||| || |||| |||||||...
    CCCTCTAGAGGGAAGCGATTTCTGTGATCTGAAAGAAAAGAAAATGGTTCCCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|262205392|ref|NR_030211.1| Homo sapiens microRNA 518d (MIR518D), microRNA
    length: 87
    e value: 3.53765e-12
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||| |||||||| |||||||||||||||||||||  |||| |||||| ||| |||...
    CCCTCTAGAGGGAAGCACTTTCTGTTGTCTGAAAGAAACCAAAGCGCTTCCCTTTGGAG...
    ****ALIGNMENT****
    sequence: gi|262205392|ref|NR_030211.1| Homo sapiens microRNA 518d (MIR518D), microRNA
    length: 87
    e value: 0.000491741
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||| | |||||||||||||   || | | | | || ||||||||||||| | |||||||...
    CTCCAAAGGGAAGCGCTTTGGTTTCTTTCAGACAACAGAAAGTGCTTCCCTCTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171300|ref|NR_035855.1| Pan troglodytes microRNA mir-520a (MIR520A), microRNA
    length: 84
    e value: 1.23476e-11
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGA...
    ||||| | |||||||  ||||||||||||||| |||||||||||||||||| ||| ||...
    CCCTCCAGAGGGAAGTACTTTCTGTTGTCTGAGAGAAAAGAAAGTGCTTCCCTTTGGA...
    ****ALIGNMENT****
    sequence: gi|301171161|ref|NR_035839.1| Pan troglodytes microRNA mir-516b-1 (MIR516B-1), microRNA
    length: 89
    e value: 1.23476e-11
    GAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||| ||||||||| | ||||||||||||||||||||||||| ||||||...
    GAAGCACTTTCTGTTTTGTGAAAGAAAAGAAAGTGCTTCCTTTCAGAGGG...
    ****ALIGNMENT****
    sequence: gi|270133259|ref|NR_032576.1| Macaca mulatta microRNA mir-519d (MIR519D), microRNA
    length: 91
    e value: 1.23476e-11
    CCCTCTACAGGGAAGCGCTTTCTGTT----GTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||| | |||||||| |||||||||    |||||| |||||| |||||||||||||||||||...
    CCCTCCAAAGGGAAGCACTTTCTGTTTGTTGTCTGAGAGAAAACAAAGTGCTTCCTTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|270133259|ref|NR_032576.1| Macaca mulatta microRNA mir-519d (MIR519D), microRNA
    length: 91
    e value: 0.00171634
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGA---AAGAAA-AGAAAGTGCTTCCTTTTAGAGGG...
    ||||| | |||||| |||| | || ||| |   || ||| ||||||||||||| ||| |||||...
    CTCTAAAAGGAAGCACTTTGTTTTCTCTCAGACAACAAACAGAAAGTGCTTCCCTTTGGAGGG...
    ****ALIGNMENT****
    sequence: gi|270133218|ref|NR_032568.1| Macaca mulatta microRNA mir-518b (MIR518B), microRNA
    length: 83
    e value: 1.23476e-11
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGG...
    ||||| | |||||||| |||||||||||||||||||||  |||| ||| || ||||||||...
    CCCTCCAGAGGGAAGCACTTTCTGTTGTCTGAAAGAAAGCAAAGCGCTCCCCTTTAGAGG...
    ****ALIGNMENT****
    sequence: gi|270133218|ref|NR_032568.1| Macaca mulatta microRNA mir-518b (MIR518B), microRNA
    length: 83
    e value: 0.00171634
    CCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    |||||| |||| ||||||||   || | | | | || ||||||||||||| | | |||||...
    CCTCTAAAGGGGAGCGCTTTGCTTTCTTTCAGACAACAGAAAGTGCTTCCCTCTGGAGGG...
    ****ALIGNMENT****
    sequence: gi|270132725|ref|NR_032718.1| Macaca mulatta microRNA mir-516 (MIR516), microRNA
    length: 90
    e value: 1.23476e-11
    GAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||| ||||||||| | ||||||||||||||||||||||||| ||||||...
    GAAGCACTTTCTGTTTTGTGAAAGAAAAGAAAGTGCTTCCTTTCAGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205396|ref|NR_030212.1| Homo sapiens microRNA 516b-1 (MIR516B1), microRNA
    length: 90
    e value: 1.23476e-11
    GAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||| ||||||||| | ||||||||||||||||||||||||| ||||||...
    GAAGCACTTTCTGTTTTGTGAAAGAAAAGAAAGTGCTTCCTTTCAGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205294|ref|NR_030189.1| Homo sapiens microRNA 520a (MIR520A), microRNA
    length: 85
    e value: 1.23476e-11
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGA...
    ||||| | |||||||  ||||||||||||||| |||||||||||||||||| ||| ||...
    CCCTCCAGAGGGAAGTACTTTCTGTTGTCTGAGAGAAAAGAAAGTGCTTCCCTTTGGA...
    ****ALIGNMENT****
    sequence: gi|301171602|ref|NR_035634.1| Pan troglodytes microRNA mir-1283 (MIR1283), microRNA
    length: 86
    e value: 4.30974e-11
    CTCTACA--GGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    |||||||  || ||||||||||||||||| ||||||| |||||| |||||| ||| |||||...
    CTCTACAAAGGAAAGCGCTTTCTGTTGTCAGAAAGAAGAGAAAGCGCTTCCCTTTTGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171602|ref|NR_035634.1| Pan troglodytes microRNA mir-1283 (MIR1283), microRNA
    length: 86
    e value: 0.00599063
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGC-TTCCTTT-TAGAG...
    ||||| | ||||||||||||||| || | |   | || |||||| || ||||||| |||||...
    CCCTCAAAAGGGAAGCGCTTTCTCTTCTTTCTGACAACAGAAAGCGCTTTCCTTTGTAGAG...
    ****ALIGNMENT****
    sequence: gi|301171454|ref|NR_035872.1| Pan troglodytes microRNA mir-527 (MIR527), microRNA
    length: 86
    e value: 4.30974e-11
    AGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTT...
    |||||||| ||||||||||||| |||||||||||||||||||| |||...
    AGGGAAGCCCTTTCTGTTGTCTAAAAGAAAAGAAAGTGCTTCCCTTT...
    ****ALIGNMENT****
    sequence: gi|301171333|ref|NR_035858.1| Pan troglodytes microRNA mir-520d (MIR520D), microRNA
    length: 86
    e value: 4.30974e-11
    CTCTACA--GGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    |||||||  ||||||| ||||||||||||| |||||||||||||||||||  ||| | |||...
    CTCTACAAAGGGAAGCCCTTTCTGTTGTCTAAAAGAAAAGAAAGTGCTTCTCTTTGGTGGG...
    ****ALIGNMENT****
    sequence: gi|301171291|ref|NR_035854.1| Pan troglodytes microRNA mir-519e (MIR519E), microRNA
    length: 86
    e value: 4.30974e-11
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||| | | || ||| |||||||||||||||||||||| ||||||| |||||||||||...
    CTCCAAAAGGGAGCACTTTCTGTTGTCTGAAAGAAAACAAAGTGCCTCCTTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|301171219|ref|NR_035846.1| Pan troglodytes microRNA mir-518c (MIR518C), microRNA
    length: 100
    e value: 4.30974e-11
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||  |||||||| |||||||||||||||||||||| |||| |||||  |||||||...
    CTCTGGAGGGAAGCACTTTCTGTTGTCTGAAAGAAAACAAAGCGCTTCTCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|301171219|ref|NR_035846.1| Pan troglodytes microRNA mir-518c (MIR518C), microRNA
    length: 100
    e value: 0.00599063
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||| || |||||||||| | || | | | | || ||||||||||||| |  ||||...
    CTCTAAAGAGAAGCGCTTTGTTTTCTTTCAGACAACAGAAAGTGCTTCCCTCCAGAG...
    ****ALIGNMENT****
    sequence: gi|269847477|ref|NR_031696.1| Homo sapiens microRNA 1283-2 (MIR1283-2), microRNA
    length: 87
    e value: 4.30974e-11
    TCTACA--GGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||  || |||||||||||||||||||||||||||||||  |||||| ||| |||...
    TCTACAAAGGAAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAATCGCTTCCCTTTGGAG...
    ****ALIGNMENT****
    sequence: gi|262205435|ref|NR_030219.1| Homo sapiens microRNA 527 (MIR527), microRNA
    length: 85
    e value: 4.30974e-11
    AGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTT...
    |||||||| ||||||||||||| |||||||||||||||||||| |||...
    AGGGAAGCCCTTTCTGTTGTCTAAAAGAAAAGAAAGTGCTTCCCTTT...
    ****ALIGNMENT****
    sequence: gi|262205336|ref|NR_030199.1| Homo sapiens microRNA 518c (MIR518C), microRNA
    length: 101
    e value: 4.30974e-11
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||  |||||||| |||||||||||||||||||||| |||| |||||  |||||||...
    CTCTGGAGGGAAGCACTTTCTGTTGTCTGAAAGAAAACAAAGCGCTTCTCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|262205336|ref|NR_030199.1| Homo sapiens microRNA 518c (MIR518C), microRNA
    length: 101
    e value: 0.00599063
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||| || |||||||||| | || | | | | || ||||||||||||| |  ||||...
    CTCTAAAGAGAAGCGCTTTGTTTTCTTTCAGACAACAGAAAGTGCTTCCCTCCAGAG...
    ****ALIGNMENT****
    sequence: gi|270133264|ref|NR_032577.1| Macaca mulatta microRNA mir-520a (MIR520A), microRNA
    length: 85
    e value: 1.50425e-10
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGA...
    ||||| | |||||||   ||||||||||||||| ||||||||||||||||| ||| ||...
    CCCTCCAGAGGGAAGTATTTTCTGTTGTCTGAAGGAAAAGAAAGTGCTTCCCTTTGGA...
    ****ALIGNMENT****
    sequence: gi|269846946|ref|NR_031573.1| Homo sapiens microRNA 1283-1 (MIR1283-1), microRNA
    length: 87
    e value: 1.50425e-10
    TCTACA--GGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||||  || ||||||||||||||||| ||||||| |||||| |||||| ||| |||||...
    TCTACAAAGGAAAGCGCTTTCTGTTGTCAGAAAGAAGAGAAAGCGCTTCCCTTTTGAGGG...
    ****ALIGNMENT****
    sequence: gi|269846946|ref|NR_031573.1| Homo sapiens microRNA 1283-1 (MIR1283-1), microRNA
    length: 87
    e value: 0.0209094
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGC-TTCCTTT-TAGA...
    ||||| | ||||||||||||||| || | |   | || |||||| || ||||||| ||||...
    CCCTCAAAAGGGAAGCGCTTTCTCTTCTTTCTGACAACAGAAAGCGCTTTCCTTTGTAGA...
    ****ALIGNMENT****
    sequence: gi|262205358|ref|NR_030204.1| Homo sapiens microRNA 520d (MIR520D), microRNA
    length: 87
    e value: 1.50425e-10
    TCTACA--GGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||||  ||||||| ||||||||||||| |||||||||||||||||||  ||| | |||...
    TCTACAAAGGGAAGCCCTTTCTGTTGTCTAAAAGAAAAGAAAGTGCTTCTCTTTGGTGGG...
    ****ALIGNMENT****
    sequence: gi|301171206|ref|NR_035844.1| Pan troglodytes microRNA mir-518a (MIR518A), microRNA
    length: 86
    e value: 5.25034e-10
    AGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTT...
    |||||||| ||||||||||||| ||||||||||||| |||||| |||...
    AGGGAAGCCCTTTCTGTTGTCTAAAAGAAAAGAAAGCGCTTCCCTTT...
    ****ALIGNMENT****
    sequence: gi|301171206|ref|NR_035844.1| Pan troglodytes microRNA mir-518a (MIR518A), microRNA
    length: 86
    e value: 0.00599063
    AGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTT...
    ||||||||||||||| || | | | | || |||||| |||||| |||...
    AGGGAAGCGCTTTCTTTTCTTTTAGACAACAGAAAGGGCTTCCCTTT...
    ****ALIGNMENT****
    sequence: gi|301171169|ref|NR_035840.1| Pan troglodytes microRNA mir-516b-2 (MIR516B-2), microRNA
    length: 89
    e value: 5.25034e-10
    GAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||| |||| |||| | ||||||||||||||||||||||||| ||||||...
    GAAGCACTTTGTGTTTTGTGAAAGAAAAGAAAGTGCTTCCTTTCAGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205387|ref|NR_030210.1| Homo sapiens microRNA 518a-1 (MIR518A1), microRNA
    length: 85
    e value: 5.25034e-10
    AGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTT...
    |||||||| |||||||||||||||||||| |||||| |||||| |||...
    AGGGAAGCCCTTTCTGTTGTCTGAAAGAAGAGAAAGCGCTTCCCTTT...
    ****ALIGNMENT****
    sequence: gi|262205387|ref|NR_030210.1| Homo sapiens microRNA 518a-1 (MIR518A1), microRNA
    length: 85
    e value: 0.00599063
    AGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTT...
    ||||||||||||||| || | | | | || |||||| |||||| |||...
    AGGGAAGCGCTTTCTCTTCTTTCAGACAACAGAAAGGGCTTCCCTTT...
    ****ALIGNMENT****
    sequence: gi|301171381|ref|NR_035863.1| Pan troglodytes microRNA mir-521-1 (MIR521-1), microRNA
    length: 86
    e value: 1.83255e-09
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||| | |||||||  ||||||||||||| ||||||||||| |  ||||| |||||||...
    CCCTCCAAAGGGAAGAACTTTCTGTTGTCTAAAAGAAAAGAACGCACTTCCCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|270133201|ref|NR_032565.1| Macaca mulatta microRNA mir-517 (MIR517), microRNA
    length: 86
    e value: 1.83255e-09
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||| | |||||| || ||||| |||| ||||||||||  |||| |||||||||||...
    CCCTCTAGATGGAAGCACTGTCTGTGGTCTAAAAGAAAAGATCGTGCATCCTTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|262205417|ref|NR_030216.1| Homo sapiens microRNA 521-1 (MIR521-1), microRNA
    length: 87
    e value: 1.83255e-09
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||| | |||||||  ||||||||||||| ||||||||||| |  ||||| |||||||...
    CCCTCCAAAGGGAAGAACTTTCTGTTGTCTAAAAGAAAAGAACGCACTTCCCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|270133214|ref|NR_032567.1| Macaca mulatta microRNA mir-518a (MIR518A), microRNA
    length: 87
    e value: 6.39622e-09
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTC-CTTT...
    ||||  | |||||||| ||||||||||||| || |||||||||||||||| ||||...
    CCCTACAAAGGGAAGCCCTTTCTGTTGTCTAAACGAAAAGAAAGTGCTTCTCTTT...
    ****ALIGNMENT****
    sequence: gi|262205407|ref|NR_030214.1| Homo sapiens microRNA 517c (MIR517C), microRNA
    length: 95
    e value: 6.39622e-09
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||| | |||||| || ||||||||||  |||||||||  |||| |||||||||||...
    CCCTCTAGATGGAAGCACTGTCTGTTGTCT--AAGAAAAGATCGTGCATCCTTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|262205377|ref|NR_030208.1| Homo sapiens microRNA 526a-2 (MIR526A2), microRNA
    length: 65
    e value: 6.39622e-09
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||||| |||||||| ||||||||||    |||||||||||  ||| |||||| ||||||...
    CCCTCTAGAGGGAAGCACTTTCTGTTG----AAAGAAAAGAACATGCATCCTTTCAGAGGG...
    ****ALIGNMENT****
    sequence: gi|270133311|ref|NR_032588.1| Macaca mulatta microRNA mir-523b (MIR523B), microRNA
    length: 89
    e value: 2.2325e-08
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCT--GAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||| || |||||||||||||||| ||   |||||| || ||| |||||| |||||||...
    CCCTCTAGAGCGAAGCGCTTTCTGTTGGCTAGAAAAGAATAGGAAGCGCTTCCCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|270133297|ref|NR_032585.1| Macaca mulatta microRNA mir-521 (MIR521), microRNA
    length: 87
    e value: 2.2325e-08
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||| | |||||||  ||||||||||||| ||||||||||| |  ||||| ||| |||...
    CCCTCCAAAGGGAAGTACTTTCTGTTGTCTAAAAGAAAAGAACGCACTTCCCTTTGGAG...
    ****ALIGNMENT****
    sequence: gi|262205401|ref|NR_030213.1| Homo sapiens microRNA 518a-2 (MIR518A2), microRNA
    length: 87
    e value: 2.2325e-08
    AGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTT...
    |||||||| ||||||||||||| |||||| |||||| |||||| |||...
    AGGGAAGCCCTTTCTGTTGTCTAAAAGAAGAGAAAGCGCTTCCCTTT...
    ****ALIGNMENT****
    sequence: gi|262205401|ref|NR_030213.1| Homo sapiens microRNA 518a-2 (MIR518A2), microRNA
    length: 87
    e value: 0.00599063
    AGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTT...
    ||||||||||||||| || | | | | || |||||| |||||| |||...
    AGGGAAGCGCTTTCTCTTCTTTTAGACAACAGAAAGGGCTTCCCTTT...
    ****ALIGNMENT****
    sequence: gi|301171372|ref|NR_035862.1| Pan troglodytes microRNA mir-520h (MIR520H), microRNA
    length: 89
    e value: 7.79219e-08
    CCCTCTACAGGGAAGCGCTTTCTGTT----GTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||| ||| |||| |||||||||    |||||| | |||| ||||||||||| |||||||...
    CCCTCTAGAGG-AAGCACTTTCTGTTTGTTGTCTGAGAAAAAACAAAGTGCTTCCCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|301171372|ref|NR_035862.1| Pan troglodytes microRNA mir-520h (MIR520H), microRNA
    length: 89
    e value: 0.00171634
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGA---AAGAAA-AGAAAGTGCTTCCTTTTAGAGGG...
    ||||| |||||||| |||| | || ||| |   || ||| ||||||||||||| | |||||||...
    CTCTAAAGGGAAGCACTTTGTTTTTTCTCAGACAACAAACAGAAAGTGCTTCC-TCTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205412|ref|NR_030215.1| Homo sapiens microRNA 520h (MIR520H), microRNA
    length: 88
    e value: 7.79219e-08
    CCCTCTACAGGGAAGCGCTTTCTGTT----GTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||| ||| |||| |||||||||    |||||| | |||| ||||||||||| |||||||...
    CCCTCTAGAGG-AAGCACTTTCTGTTTGTTGTCTGAGAAAAAACAAAGTGCTTCCCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|262205412|ref|NR_030215.1| Homo sapiens microRNA 520h (MIR520H), microRNA
    length: 88
    e value: 0.00171634
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGA---AAGAAA-AGAAAGTGCTTCCTTTTAGAGGG...
    ||||| |||||||| |||| | || ||| |   || ||| ||||||||||||| | |||||||...
    CTCTAAAGGGAAGCACTTTGTTTTTTCTCAGACAACAAACAGAAAGTGCTTCC-TCTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205366|ref|NR_030206.1| Homo sapiens microRNA 520g (MIR520G), microRNA
    length: 90
    e value: 7.79219e-08
    CCCTCTACAGGGAAGCGCTTTCTGTT----GTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||| ||| |||| |||||||||    |||||| | |||| ||||||||||| |||||||...
    CCCTCTAGAGG-AAGCACTTTCTGTTTGTTGTCTGAGAAAAAACAAAGTGCTTCCCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|262205366|ref|NR_030206.1| Homo sapiens microRNA 520g (MIR520G), microRNA
    length: 90
    e value: 0.00171634
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGA---AAGAAA-AGAAAGTGCTTCCTTTTAGAGGG...
    ||||| |||||||| |||| | || ||| |   || ||| ||||||||||||| | |||||||...
    CTCTAAAGGGAAGCACTTTGTTTTTTCTCAGACAACAAACAGAAAGTGCTTCC-TCTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|270132721|ref|NR_032717.1| Macaca mulatta microRNA mir-524 (MIR524), microRNA
    length: 85
    e value: 2.71974e-07
    CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTT...
    ||||  | ||| |||| |||||| |||||| ||||||||||| |||||||| |||...
    CCCTACAAAGGCAAGCACTTTCTCTTGTCTAAAAGAAAAGAAGGTGCTTCCCTTT...
    ****ALIGNMENT****
    sequence: gi|262205276|ref|NR_030185.1| Homo sapiens microRNA 519e (MIR519E), microRNA
    length: 84
    e value: 9.49283e-07
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||| | | || ||| |||||||||   |||||||||| ||||||| |||||||||||...
    CTCCAAAAGGGAGCACTTTCTGTT---TGAAAGAAAACAAAGTGCCTCCTTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|301171363|ref|NR_035861.1| Pan troglodytes microRNA mir-520g (MIR520G), microRNA
    length: 89
    e value: 3.31332e-06
    CCCTCTACAGGGAAGCGCTTTCTGTT----GTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||| ||| | || |||||||||    |||||| | |||| ||||||||||| |||||||...
    CCCTCTAGAGG-ACGCACTTTCTGTTTGTTGTCTGAGAAAAAACAAAGTGCTTCCCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|270133196|ref|NR_032564.1| Macaca mulatta microRNA mir-516a-2 (MIR516A-2), microRNA
    length: 89
    e value: 3.31332e-06
    GAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||| ||||||||||||| ||||||||| ||||| ||||||   |||||...
    GAAGCACTTTCTGTTGTCT-AAAGAAAAGGAAGTGTTTCCTTCCCGAGGG...
    ****ALIGNMENT****
    sequence: gi|270133191|ref|NR_032563.1| Macaca mulatta microRNA mir-516a-1 (MIR516A-1), microRNA
    length: 89
    e value: 3.31332e-06
    GAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||| ||||||||||||| ||||||||| ||||| ||||||   |||||...
    GAAGCACTTTCTGTTGTCT-AAAGAAAAGGAAGTGTTTCCTTCCCGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171387|ref|NR_035864.1| Pan troglodytes microRNA mir-521-2 (MIR521-2), microRNA
    length: 86
    e value: 1.15646e-05
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||| | |||||||   ||||| |||||| ||||||||||| |  ||||| |||||||...
    CTCCAAAGGGAAGAATTTTCTCTTGTCTAAAAGAAAAGAACGCACTTCCCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|262205353|ref|NR_030203.1| Homo sapiens microRNA 521-2 (MIR521-2), microRNA
    length: 87
    e value: 1.15646e-05
    CTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||| | |||||||   ||||| |||||| ||||||||||| |  ||||| |||||||...
    CTCCAAAGGGAAGAATTTTCTCTTGTCTAAAAGAAAAGAACGCACTTCCCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|270133237|ref|NR_032572.1| Macaca mulatta microRNA mir-518f (MIR518F), microRNA
    length: 87
    e value: 0.000140886
    CCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    |||||  ||||||||||||||| |  | | | | || | ||||||||||| | |||||||...
    CCTCTGAAGGGAAGCGCTTTCTTTCCTTTCACACAAGATAAAGTGCTTCCCTCTAGAGGG...
    ****ALIGNMENT****
    sequence: gi|301171283|ref|NR_035853.1| Pan troglodytes microRNA mir-519d (MIR519D), microRNA
    length: 87
    e value: 0.00171634
    CCCTCTACAGGGAAGCGCTTTCTG-TTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||| | |||||||||||||||| |||| |      ||| ||||||| ||| |||||||...
    CCCTCCAAAGGGAAGCGCTTTCTGTTTGTTTTCTCTCAAACAAAGTGCCTCCCTTTAGAG...
    ****ALIGNMENT****
    sequence: gi|270133288|ref|NR_032583.1| Macaca mulatta microRNA mir-520g (MIR520G), microRNA
    length: 90
    e value: 0.00171634
    CCCTCTACAGGGAAGCGCTTTCTGTT----GTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||||| ||  |||| |||||||||    ||||||   |||| ||||||||||| || ||||...
    CCCTCTAGAGA-AAGCACTTTCTGTTTGTTGTCTGAGGAAAAACAAAGTGCTTCCCTTCAGAG...
    ****ALIGNMENT****
    sequence: gi|262205371|ref|NR_030207.1| Homo sapiens microRNA 516b-2 (MIR516B2), microRNA
    length: 85
    e value: 0.00171634
    GAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG...
    ||||| |||| |||| | ||||||     |||||||||||||| ||||||...
    GAAGCACTTTGTGTTTTGTGAAAG-----AAAGTGCTTCCTTTCAGAGGG...
    ****ALIGNMENT****
    sequence: gi|262205349|ref|NR_030202.1| Homo sapiens microRNA 519d (MIR519D), microRNA
    length: 88
    e value: 0.00171634
    CCCTCTACAGGGAAGCGCTTTCTG-TTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAG...
    ||||| | |||||||||||||||| |||| |      ||| ||||||| ||| |||||||...
    CCCTCCAAAGGGAAGCGCTTTCTGTTTGTTTTCTCTTAAACAAAGTGCCTCCCTTTAGAG...

## Open CV (Parts 1-3)
```python
import numpy as np
import matplotlib.pyplot as plt
%matplotlib inline
```


```python
import cv2
```


```python
img = cv2.imread("mushroom.jpg")
```


```python
type(img)
```




    numpy.ndarray




```python
img_wrong = cv2.imread('wrong/path/doesnot/abcdegh.jpg')
```


```python
type(img_wrong)
```




    NoneType




```python
plt.imshow(img)
```




    <matplotlib.image.AxesImage at 0x7fc3588cac50>




![png](output_6_1.png)



```python
fix_img = cv2.cvtColor(img, cv2.COLOR_BGR2RGB)
```


```python
plt.imshow(fix_img)
```




    <matplotlib.image.AxesImage at 0x7fc3587e9e10>




![png](output_8_1.png)



```python
img_gray = cv2.imread("mushroom.jpg", cv2.IMREAD_GRAYSCALE)
img_gray.shape
```




    (1200, 1200)




```python
plt.imshow(img_gray)
```




    <matplotlib.image.AxesImage at 0x7fc35875bfd0>




![png](output_10_1.png)



```python
plt.imshow(img_gray, cmap = "gray")
```




    <matplotlib.image.AxesImage at 0x7fc355316810>




![png](output_11_1.png)



```python
fix_img.shape
```




    (1200, 1200, 3)




```python
new_img = cv2.resize(fix_img,(1000,400))
plt.imshow(new_img)
```




    <matplotlib.image.AxesImage at 0x7fc3552864d0>




![png](output_13_1.png)



```python
new_img.shape
```




    (400, 1000, 3)




```python
w_ratio = 0.5
h_ratio = 0.5


new_img = cv2.resize(fix_img, (0,0), fix_img, w_ratio, h_ratio)
```


```python
plt.imshow(new_img)
```




    <matplotlib.image.AxesImage at 0x7fc3551e9290>




![png](output_16_1.png)



```python
new_img.shape
```




    (600, 600, 3)




```python
flip_img = cv2.flip(fix_img, 0)
plt.imshow(flip_img)
```




    <matplotlib.image.AxesImage at 0x7fc3551d7090>




![png](output_18_1.png)



```python
flip_img2 = cv2.flip(fix_img, -1)
plt.imshow(flip_img2)
```




    <matplotlib.image.AxesImage at 0x7fc355138650>




![png](output_19_1.png)



```python
type(fix_img)
```




    numpy.ndarray




```python
cv2.imwrite('mushroom_fixed_image.jpg', fix_img)
```




    True




```python
import cv2
import matplotlib.pyplot as plt
%matplotlib inline
```


```python
img = cv2.imread("mushroom.jpg")
```


```python
plt.imshow(img)
```




    <matplotlib.image.AxesImage at 0x7fc35511a310>




![png](output_24_1.png)



```python
img1 = cv2.cvtColor(img, cv2.COLOR_BGR2RGB)
```


```python
plt.imshow(img1)
```




    <matplotlib.image.AxesImage at 0x7fc3586ef9d0>




![png](output_26_1.png)



```python
img2 = cv2.cvtColor(img, cv2.COLOR_BGR2HSV)
```


```python
plt.imshow(img2)
```




    <matplotlib.image.AxesImage at 0x7fc358667250>




![png](output_28_1.png)



```python
img3 = cv2.cvtColor(img, cv2.COLOR_BGR2HLS)
```


```python
plt.imshow(img3)
```




    <matplotlib.image.AxesImage at 0x7fc35864c190>




![png](output_30_1.png)



```python
img1 = cv2.imread('donotcopy.jpg')
img2 = cv2.imread('mushroom.jpg')
```


```python
plt.imshow(img1)
```




    <matplotlib.image.AxesImage at 0x7fc3585a5390>




![png](output_32_1.png)



```python
img1 = cv2.cvtColor(img1, cv2.COLOR_BGR2RGB)
img2 = cv2.cvtColor(img2, cv2.COLOR_BGR2RGB)
```


```python
plt.imshow(img2)
```




    <matplotlib.image.AxesImage at 0x7fc35857f950>




![png](output_34_1.png)



```python
img1 = cv2.resize(img1, (1200,1200))
img2 = cv2.resize(img2, (1200, 1200))
```


```python
alpha = 0.5
beta = 0.5
```


```python
blended = cv2.addWeighted(img1, alpha, img2, beta, gamma=0)
```


```python
plt.imshow(blended)
```




    <matplotlib.image.AxesImage at 0x7fc3584f2cd0>




![png](output_38_1.png)



```python
alpha = 0.8
beta = 0.2

blended1 = cv2.addWeighted(img1, alpha, img2, beta, 0)
plt.imshow(blended1)
```




    <matplotlib.image.AxesImage at 0x7fc35845b390>




![png](output_39_1.png)



```python
img1 = cv2.imread('donotcopy.jpg')
iimg2 = cv2.imread('mushroom.jpg')

img1 = cv2.cvtColor(img1, cv2.COLOR_BGR2RGB)
img1 = cv2.cvtColor(img2, cv2.COLOR_BGR2RGB)

img1 = cv2.resize(img1, (200,200))
```


```python
large_img = img2
small_img = img1

x_offset = 0
y_offset = 0

x_end = x_offset + small_img.shape[1]
y_end = y_offset + small_img.shape[0]

large_img[y_offset:y_end, x_offset:x_end] = small_img

plt.imshow(large_img)
```




    <matplotlib.image.AxesImage at 0x7fc358438e90>




![png](output_41_1.png)



```python
import cv2
import matplotlib.pyplot as plt
%matplotlib inline
```


```python
img = cv2.imread('rainbow.jpg')
```


```python
plt.imshow(img)
```




    <matplotlib.image.AxesImage at 0x7fc35839f650>




![png](output_44_1.png)



```python
img = cv2.imread('rainbow.jpg', 0)
```


```python
plt.imshow(img, cmap = 'gray')
```




    <matplotlib.image.AxesImage at 0x7fc3583872d0>




![png](output_46_1.png)



```python
ret1, thresh1 = cv2.threshold(img, 127, 255, cv2.THRESH_BINARY)
```


```python
ret1
```




    127.0




```python
plt.imshow(thresh1, cmap = "gray")
```




    <matplotlib.image.AxesImage at 0x7fc35506efd0>




![png](output_49_1.png)



```python
img2 = cv2.imread('rainbow.jpg', 0)
ret1, thresh1 = cv2.threshold(img2, 127, 255, cv2.THRESH_TRUNC)
plt.imshow(thresh1, cmap = "gray")
```




    <matplotlib.image.AxesImage at 0x7fc354fef210>




![png](output_50_1.png)



```python

```


```python
img3 = cv2.imread('rainbow.jpg', 0)
ret1, thresh1 = cv2.threshold(img3, 127, 255, cv2.THRESH_TOZERO)
plt.imshow(thresh1, cmap = "gray")
```




    <matplotlib.image.AxesImage at 0x7fc354fbb710>




![png](output_52_1.png)



```python
img_r = cv2.imread('crossword.jpg', 0)
plt.imshow(img_r, cmap = 'gray')
```




    <matplotlib.image.AxesImage at 0x7fc354f9e0d0>




![png](output_53_1.png)



```python
def show_pic(img):
    fig = plt.figure(figsize = (15,15))
    ax = fig.add_subplot(111)
    ax.imshow(img, cmap = 'gray')
```


```python
show_pic(img_r)
```


![png](output_55_0.png)



```python
ret, th1 = cv2.threshold(img_r, 127, 255, cv2.THRESH_BINARY)
show_pic(th1)
```


![png](output_56_0.png)



```python
ret, th1 = cv2.threshold(img_r, 200, 255, cv2.THRESH_BINARY)
show_pic(th1)
```


![png](output_57_0.png)



```python
th2 = cv2.adaptiveThreshold(img_r, 255, cv2.ADAPTIVE_THRESH_MEAN_C, cv2.THRESH_BINARY, 11, 8)
```


```python
show_pic(th2)
```


![png](output_59_0.png)



```python
blended = cv2.addWeighted(src1 = th1, alpha = 0.6,
                         src2 = th2, beta = 0.4, gamma = 0)

show_pic(blended)
```


![png](output_60_0.png)



```python
th3 = cv2.adaptiveThreshold(img_r, 255, cv2.ADAPTIVE_THRESH_MEAN_C, cv2.THRESH_BINARY, 11, 8)

blended = cv2.addWeighted(src1=th1, alpha=0.6,
                          src2=th2, beta=0.4, gamma=0)

show_pic(blended)
```


![png](output_61_0.png)


## Aspect Detection (Corner and Edge Detection)
### Corner Detection
```python
import cv2
import numpy as np
import matplotlib.pyplot as plt
%matplotlib inline
```


```python
real_chess = cv2.imread('fakechessboard.jpg')
flat_chess = cv2.cvtColor(real_chess, cv2.COLOR_BGR2RGB) 
plt.imshow(flat_chess)
plt.show()
```


![png](output_1_0.png)



```python
gray_flat_chess = cv2.cvtColor(flat_chess, cv2.COLOR_BGR2GRAY)
plt.imshow(gray_flat_chess, cmap = "gray")
```




    <matplotlib.image.AxesImage at 0x7fef743a6dd0>




![png](output_2_1.png)



```python
real_chess = cv2.imread("realchessboard.jpeg")
real_chess = cv2.cvtColor(real_chess, cv2.COLOR_BGR2RGB)
```


```python
plt.imshow(real_chess)
```




    <matplotlib.image.AxesImage at 0x7fef75b4c310>




![png](output_4_1.png)



```python
gray_real_chess = cv2.cvtColor(real_chess, cv2.COLOR_BGR2GRAY)
plt.imshow(gray_real_chess, cmap = 'gray')
```




    <matplotlib.image.AxesImage at 0x7fef74d3f110>




![png](output_5_1.png)



```python
gray = np.float32(gray_flat_chess)
dst = cv2.cornerHarris(src = gray, blockSize = 2, ksize = 3, k =0.04)

dst = cv2.dilate(dst, None)
```


```python
flat_chess[dst>0.01*dst.max()] = [255,0,0]

plt.imshow(flat_chess)
```




    <matplotlib.image.AxesImage at 0x7fef742fb950>




![png](output_7_1.png)



```python
gray = np.float32(gray_real_chess)
dst = cv2.cornerHarris(src = gray, blockSize =2, ksize=3, k=0.04)
dst = cv2.dilate(dst, None)

real_chess[dst>0.01*dst.max()] = [255, 0, 0]

plt.imshow(real_chess)
```




    <matplotlib.image.AxesImage at 0x7fef7425c150>




![png](output_8_1.png)



```python
#Shi-Tomasi Corner Detection

corners = cv2.goodFeaturesToTrack(gray_flat_chess, 64, 0.01, 10)
```


```python
corners = np.int0(corners)

for i in corners:
    x,y = i.ravel()
    cv2.circle(flat_chess, (x,y),3,(255,0,0), -1)

plt.imshow(flat_chess)
    
```




    <matplotlib.image.AxesImage at 0x7fef7423e550>




![png](output_10_1.png)



```python
corners = cv2.goodFeaturesToTrack(gray_real_chess, 100, 0.01, 10)

corners = np.int0(corners)

for i in corners:
    x, y = i.ravel()
    cv2.circle(real_chess, (x,y), 3, (0,255,0), -1)

plt.imshow(real_chess)
```




    <matplotlib.image.AxesImage at 0x7fef7419a690>




![png](output_11_1.png)



### Edge Detection
```python
import cv2
```


```python
import numpy as np
```


```python
import matplotlib.pyplot as plt
%matplotlib inline
```


```python
img = cv2.imread("mushroom.jpg")
plt.imshow(img)
```




    <matplotlib.image.AxesImage at 0x7fd9bebd2450>




![png](output_3_1.png)



```python
edges = cv2.Canny(image = img, threshold1 = 127, threshold2 = 127)

plt.imshow(edges)
```




    <matplotlib.image.AxesImage at 0x7fd9be315110>




![png](output_4_1.png)



```python
med_value = np.median(img)
med_value
```




    70.0




```python
lower = int(max(0, 0.7*med_value))
upper = int(min(255,1.3*med_value))

edges = cv2.Canny(img, threshold1 = lower, threshold2 = upper)

plt.imshow(edges)
```




    <matplotlib.image.AxesImage at 0x7fd9be27fdd0>




![png](output_6_1.png)



```python
edges = cv2.Canny(image = img, threshold1 = lower, threshold2 = upper +100)

plt.imshow(edges)
```




    <matplotlib.image.AxesImage at 0x7fd9bc5b8b10>




![png](output_7_1.png)



```python
blurred_img = cv2.blur(img, ksize = (5,5))

edges = cv2.Canny(image=blurred_img,
                  threshold1 = lower,
                  threshold2 = upper)
plt.imshow(edges)
```




    <matplotlib.image.AxesImage at 0x7fd9bc5a0ed0>




![png](output_8_1.png)



```python
blurred_img = cv2.blur(img, ksize = (7,7))

edges = cv2.Canny(image=blurred_img,
                 threshold1 = lower,
                 threshold2 = upper)
plt.imshow(edges)
```




    <matplotlib.image.AxesImage at 0x7fd9bc492890>




![png](output_9_1.png)



```python
blurred_img = cv2.blur(img, ksize = (7,7))

edges = cv2.Canny(image=blurred_img,
                 threshold1 = lower,
                 threshold2 = upper + 100)
plt.imshow(edges)
```




    <matplotlib.image.AxesImage at 0x7fd9bc473b50>




![png](output_10_1.png)


## Feature Detection

### Feature Matches
```python
import cv2
import numpy as np
import matplotlib.pyplot as plt
%matplotlib inline
```


```python
def display(img, cmap = 'gray'):
    fig = plt.figure(figsize = (12, 10))
    ax = fig.add_subplot(111)
    ax.imshow(img, cmap = 'gray')
```


```python
apple_jacks = cv2.imread("applejacks.jpg")
display(apple_jacks)
```


![png](output_2_0.png)



```python
cereals = cv2.imread('manycereal.jpg', 0)
display(cereals)
```


![png](output_3_0.png)



```python
orb = cv2.ORB_create()

kp1,des1 = orb.detectAndCompute(apple_jacks, mask=None)
kp2,des2 = orb.detectAndCompute(cereals, mask=None)
```


```python
bf = cv2.BFMatcher(cv2.NORM_HAMMING, crossCheck = True)
matches= bf.match(des1, des2)
```


```python
matches = sorted(matches, key = lambda X:X.distance)
```


```python
apple_jacks_matches = cv2.drawMatches(apple_jacks, kp1, cereals,  kp2, matches[:25], None, flags = 2)
```


```python
display(apple_jacks_matches)
```


![png](output_8_0.png)



```python
sift = cv2.SIFT_create()
```


```python
kp1, des1 = sift.detectAndCompute(apple_jacks, None)
kp2, des2 = sift.detectAndCompute(cereals, None)
```


```python
bf = cv2.BFMatcher()
matches = bf.knnMatch(des1, des2, k=2)
```


```python

```




### Object Detection
```python
import cv2
```


```python
import numpy as np
```


```python
import matplotlib.pyplot as plt
```


```python
%matplotlib inline
```


```python
full = cv2.imread('sunflower.jpg')
```


```python
full = cv2.cvtColor(full, cv2.COLOR_BGR2RGB)
```


```python
plt.imshow(full)
```




    <matplotlib.image.AxesImage at 0x7f78c5eb3690>




![png](output_6_1.png)



```python
test = cv2.imread('sunflower.jpg')
```


```python
test = cv2.cvtColor(test, cv2.COLOR_BGR2RGB)
```


```python
plt.imshow(test)
```




    <matplotlib.image.AxesImage at 0x7f78cc0ca850>




![png](output_9_1.png)



```python
print('Test image shape:', full.shape)
print('Training image shape:', test.shape)
```

    Test image shape: (994, 661, 3)
    Training image shape: (994, 661, 3)



```python
methods = ['cv2.TM_CCOEFF', 'cv2.TM_CCOEFF_NORMED', 'cv2.TM_CCORR', 'cv2.TM_SQDIFF', 'cv2.TM_SQDIFF_NORMED']
```


```python
for m in methods:
    test_copy = test.copy()
    method = eval(m)
    
    res = cv2.matchTemplate(test_copy, full, method)
    
    min_val, max_val, min_loc, max_loc = cv2.minMaxLoc(res)
    
    if method in [cv2.TM_SQDIFF_NORMED]:
        top_left = min_loc
    else:
        top_left = max_loc
    
    height, width, channels = full.shape
    bottom_right = (top_left[0] + width, top_left[1] + height)
    
    cv2.rectangle(test_copy, top_left, bottom_right, (255,0,0), 10)  # Fixed this line
    
    plt.subplot(121)
    plt.imshow(res)
    plt.title("Heatmap of template matching")
    
    plt.subplot(122)
    plt.imshow(test_copy)
    plt.title('Detection of template')
    
    plt.suptitle(m)  # Corrected to plt.suptitle
    
    plt.show()
    print('\n')
    print('\n')
```


![png](output_12_0.png)


    
    
    
    



![png](output_12_2.png)


    
    
    
    



![png](output_12_4.png)


    
    
    
    



![png](output_12_6.png)


    
    
    
    



![png](output_12_8.png)


    





















