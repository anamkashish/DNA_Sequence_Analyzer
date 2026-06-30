print("DNA Sequence Analyzer")
seq = input("Enter DNA Sequence: ")
if not seq:
    print("="*40)
    print("\tDNA ANALYSIS REPORT")
    print("="*40)
    print("\nStatus : NO SEQUENCE ENTERED\n")
    print("="*40)
else:
    s= seq.upper()
    s1= set(s)
    s2= {"A","T","G","C"}
    if not(s1-s2):
        print("="*40)
        print("\tDNA ANALYSIS REPORT")
        print("="*40)
        print(f"Length = {len(s)}")

        print("\nNucleotide Counts")
        print("-"*15)
        print(f"A Count : {s.count("A")}")
        print(f"T Count : {s.count("T")}")
        print(f"G Count : {s.count("G")}")
        print(f"C Count : {s.count("C")}")

        print("\nSequence Operations")
        print("-"*15)
        print(f"RNA Sequence : {s.replace("T","U")}")
        print(f"Reverse Sequence : {s[::-1]}")

        print("\nAnalysis")
        print("-"*15)
        content = (100*(s.count("G") + s.count("C")))/len(s)
        print(f"GC Content : {round(content,2)}%")

        print("\nStatus : VALID DNA SEQUENCE\n")
        print("="*40)

    else:
        print("="*40)
        print("\tDNA ANALYSIS REPORT")
        print("="*40)
        print("\nStatus : INVALID DNA SEQUENCE\n")
        print(f"Invalid charcter found : {s1-s2}\n")
        print("="*40)
