print("-"*8,end="")
print(" DNA SEQUENCE ANALYZER ",end="")
print("-"*8)
s = input("Enter DNA Sequence: ").upper().strip()
seq= set(s)
valid={"A","T","G","C"}
if not seq.issubset(valid):
    print("="*40)
    print("\tDNA ANALYSIS REPORT")
    print("="*40)
    print("\nStatus : INVALID DNA SEQUENCE\n")
    
elif not s:
    print("="*40)
    print("\tDNA ANALYSIS REPORT")
    print("="*40)
    print("\nStatus : NO SEQUENCE ENTERED\n")
    
else:
    print("="*40)
    print("\tDNA ANALYSIS REPORT")
    print("="*40)

    total_len = len(s)
    print(f"Sequence Length : {total_len}")

    print("\nNucleotide Compositiion")
    print("-"*15)
    freq = {"A":0,"T":0,"G":0,"C":0}
    for char in s:
        freq[char] += 1

    for char in freq:
        frequency=(freq[char] / total_len) * 100
        print(f"{char} : {freq[char]} ({round(frequency, 1)} %)")

    gc= freq["G"] + freq["C"]
    at= freq["A"] + freq["T"]

    print(f"\nGC Content : {round((gc / total_len) * 100, 2)} %")
    print(f"AT Content : {round((at / total_len) * 100, 2)} %")

    print("\nDNA Analysis")
    print("-"*15)
    print(f"Original : {s}")
    pairs = {"A":"T","T":"A","G":"C","C":"G"}
    complement=""
    for char in s:
        complement += pairs[char]
    print(f"Complement : {complement}")

    rev=s[::-1]
    comp_rev=complement[::-1]
    print(f"Reverse Sequence : {rev}")
    print(f"Reverse Complement : {comp_rev}")



    print("\nTranscription")
    print("-"*15)

    mrna = s.replace("T","U")
    print(f"mRNA (from coding strand) : {mrna}")

    RNA_CODON_TABLE = {
    'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',
    'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
    'UAU': 'Y', 'UAC': 'Y', 'UAA': '*', 'UAG': '*',  
    'UGU': 'C', 'UGC': 'C', 'UGA': '*', 'UGG': 'W',
    'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
    'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',  
    'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
    'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'}

    def translate_mrna (mrna):
        protein_seq = ""
      
        for i in range(0,len(mrna)-2,3):
            codon = mrna[i:i+3]
           
            amino_acid = RNA_CODON_TABLE.get(codon,"X")
            if amino_acid == "*":
                break
            protein_seq += amino_acid
        return protein_seq


    print("\nORF and Protein Analysis")
    print("-"*15)
    
    def orf_finder(s):
        orf=[]
        for strand in (s,comp_rev):
            for frame in range(3):
                sub=strand[frame:]
                i=0
                while i<len(sub)-2:
                    if sub[i:i+3]=="ATG":
                        j=i+3
                        found_stop = False
                        while j<len(sub)-2:
                            codon= sub[j:j+3]
                            if codon in ("TAA","TGA","TAG"):
                                found_stop = True
                                break
                            j+=3
                            
                        if found_stop:
                            orf.append(sub[i:j+3])
                            i = j+3
                        else:
                            i+=3
                    else:
                        i+=1
        return set(orf)
    
    lst = orf_finder(s)
    for i,seq in enumerate(lst, start=1):
        valid_seq=seq.replace("T","U")
        protein= translate_mrna(valid_seq)
        print(f"ORF{i} : {seq}")
        print(f"Protein : {protein}\n")

    print("Motif Search")
    print("-"*15)
    def find_motif(s,motif):
        positions=[]
        for i in range(len(s)-len(motif)+1):
            if s[i:i+len(motif)] == motif:
                positions.append(i+1)
        return positions    
    
    motif=input("Enter Motif to search : ").upper().strip()
    found_positions = find_motif(s,motif)
    if not found_positions:
        print("Motif not found")
    else:
        print(f"Motif found at position: {found_positions}")