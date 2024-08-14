import os
from collections import defaultdict

# Updated dictionary with rsIDs (SNPs) and their associated genotypes (allele pairs) for each haplogroup
haplogroup_snps = {
    # Broad Haplogroups and their subclades
    "R1a": {
        "R1a": {
            "rs17222573": ["A;A", "A;G", "G;G"],  
            "rs17307677": ["C;C", "C;T", "T;T"],
            "rs17306692": ["A;A", "A;C", "C;C"],
            "rs17250535": ["A;T", "T;T"],
        },
        "R1a1": {
            "rs17307105": ["A;A", "A;G", "G;G"],
            "rs9786587": ["G;G", "G;T", "T;T"],
            "rs17316227": ["A;A", "A;G", "G;G"],
            "rs2534636": ["A;A"],
        },
        "R1a1a": {
            "rs3908": [],
            "rs2020857": ["T;T"],
            "rs1722146": ["T;T"],
            "rs17315926": ["C;C","C;T","T;T"],
            "rs17221601": ["A;A","A;T","T;T"],
        },
        "R1a1a1": {
            "rs17316771": ["A;A","A;G","G;G"]
        },
        "R1a1a1a": {
            "rs2032622": ["T;T"]
        },
        "R1a1a1c": {
            "rs2032626": ["G;G"],
            "rs2032644": ["C;C"],
            "rs2032655": ["G;G"]
        },
        "R1a1a1h1": {
            "rs7892981": ["A;A","A;G","G;G"]
        },
        "R1a1a1i1": {
            "rs7067418": ["A;A","A;G","G;G"]
        }
    },
    "R1b": {
        "R1b": {
            "rs9786184": ["A;A"],
        },
        "R1b1": {
            "rs17249854": ["C;C","C;T","T;T"],
            "rs9786194": ["A;A","A;C","C;C"],
            "rs150173": ["A;A"],
        },
        "R1b1a": {
            "rs9785702": ["C;C"],
            "rs2917400": ["C;C","C;T","T;T"],
        },
        "R1b1a1": {
            "rs2032634": ["G;T"],
        },
        "R1b1a2": {
            "rs9786882": ["A;A","A;G","G;G"],
            "rs9786153": ["T;T"],
            "rs877756": ["C;C","C;T","T;T"],
            "rs2058276": ["A;A", "A;G", "G;G"],
        },
        "R1b1a2a": {
            "rs9785971": ["A;A","A;G","G;G"],
            "rs9786142": ["A;A","A;T","T;T"],
        },
        "R1b1b2a1": {
            "rs9785831": ["C;C","C;T","T;T"],
        },
        "R1b1a2a1a": {
            "rs9786140": ["A;A","A;G","G;G"],
        },
        "R1b1a2a1a1": {
            "rs9786076": ["C;C","C;T","T;T"],
            "rs13304168": ["C;C","C;T","T;T"],
            "rs2082033": ["C;C","C;T","T;T"],
            "rs9786283": ["A;A","A;C","C;C"],
            "rs9785659": ["A;A","A;G","G;G"],
        },
        "R1b1a2a1a1a": {
            "rs16981293": ["T;T"],
        },
        "R1b1a2a1a1a5": {
            "rs34001725": ["A;A","A;G","G;G"],
        },
        "R1b1a2a1a1a5a": {
            "rs17222279": ["A;A"],
        },
        "R1b1a2a1a1a5c1a1": {
            "rs13304625": ["A;A"],
        },
        "R1b1a2a1a1a5c1a1a": {
            "rs35760092": ["G;G","G;T","T;T"],
            "rs13305517": ["A;A","A;G","G;G"],
            "rs13305070": ["A;A","A;C","C;C"],
        },
        "R1b1a2a1a1b": {
            "rs34276300": ["A;A"],
        },
        "R1b1a2a1a1b2b1": {
            "rs1800865": ["T;T"],
        },
        "R1b1a2a1a1b3": {
            "rs1236440": ["T;T"],
        },
        "R1b1a2a1a1b3a": {
            "rs2032615": [],
        },
        "R1b1a2a1a1b3c": {
            "rs2566671": ["C;C"],
        },
        "R1b1a2a1a1b3c1": {
            "rs7067305": ["A;A","A;G","G;G"],
        },
        "R1b1a2a1a1b3c1a": {
            "rs9341273": ["C;C","C;T","T;T"],
        },
        "R1b1a2a1a1b4": {
            "rs11799226": ["G;G"],
        },
        "R1b1a2a1a1b4b": {
            "rs20321": ["A;A"],
        },
        "R1b1a2a1a1b4e": {
            "rs9306842": ["A;A","A;T","T;T"],
        },
        "R1b1a2a1a1b5": {
            "rs35199432": ["A;A","A;G","G;G"],
        },
        "R1b1a2a1b": {
            "rs9786602": ["A;A","A;C","C;C"],
        },
        "R1b1c1": {
            "rs3909": ["A;A"],
        },
    },

    # I1
    "I1": {
        "I1": {
            "rs17307252": ["A;A", "A;G", "G;G"],
            "rs871626": ["A;A", "A;G", "G;G"],
            "rs17221531": ["C;C","C;T","T;T"],
            "rs111697819": ["A;A","A;T","T;T"],
            "rs35960273": ["A;A","A;G","G;G"],
            "rs3906451": ["A;A", "A;C", "C;C"],
            "rs17307007": ["A;A","A;C","C;C"],
            "rs17222167": ["A;A", "A;G", "G;G"],
            "rs17307315": ["C;C","C;T","T;T"],
            "rs17307586": ["C;C","C;T","T;T"],
            "rs17249791": ["C;C", "C;T", "T;T"],
            "rs17222657": ["A;A","A;T","T;T"],
            "rs17250114": ["C;C","C;G","G;G"],
            "rs17316639": ["G;G","G;T","T;T"],
            "rs9341296": ["T;T"],
            "rs13447354": ["A;A"],
            "rs17316597": ["A;A"],
            "rs112707890": ["A;A","A;G","G;G"],
            "rs113686221": ["C;C","C;T","T;T"],
        },
        "I1a1a": {
            "rs9341274": ["G;G"],
        },
        "I1a1a1": {
            "rs2032637": ["G;G"],
        },
        "I1a1b": {
            "rs34626372": ["C;C"],
        },
    },

    # I2
    "I2": {
        "I2": {
            "rs35547782": ["C;C","C;T","T;T"],
            "rs17307294": ["G;G"],
        },
    },
}

def parse_dna_file(file_path):
    detected_snps = {}
    with open(file_path, 'r') as f:
        for line in f:
            if not line.startswith("#") and "rsid" not in line:
                parts = line.strip().split()
                if len(parts) >= 4:
                    rsid = parts[0]
                    allele1 = parts[3]
                    allele2 = parts[4] if len(parts) > 4 else allele1  # For 23andMe format compatibility
                    genotype = f"{allele1};{allele2}"
                    detected_snps[rsid] = genotype
    return detected_snps

def calculate_haplogroups(detected_snps):
    broad_haplogroup_counts = defaultdict(int)
    specific_haplogroup_counts = defaultdict(int)
    snps_detected_in_haplogroups = defaultdict(list)

    for broad_haplogroup, subclades in haplogroup_snps.items():
        for subclade, snps_list in subclades.items():
            detected_snps_in_subclade = [snp for snp in snps_list if snp in detected_snps and detected_snps[snp] in snps_list[snp]]
            if detected_snps_in_subclade:
                specific_haplogroup_counts[subclade] += len(detected_snps_in_subclade)
                broad_haplogroup_counts[broad_haplogroup] += len(detected_snps_in_subclade)
                snps_detected_in_haplogroups[subclade].extend(detected_snps_in_subclade)

    return broad_haplogroup_counts, specific_haplogroup_counts, snps_detected_in_haplogroups

def predict_haplogroup(broad_haplogroup_counts):
    if not broad_haplogroup_counts:
        return "Unknown"
    predicted_haplogroup = max(broad_haplogroup_counts, key=broad_haplogroup_counts.get)
    return predicted_haplogroup

def print_results(predicted_haplogroup, broad_haplogroup_counts, specific_haplogroup_counts, snps_detected_in_haplogroups):
    print(f"\033[1mPredicted Y Haplogroup:\033[0m {predicted_haplogroup}\n")
    
    print("\033[1mHaplogroup Detection:\033[0m")
    for broad_haplogroup, count in broad_haplogroup_counts.items():
        print(f"{broad_haplogroup}: {count} SNPs detected")

    print("\n\033[1mSubclade Detection:\033[0m")
    for subclade, count in specific_haplogroup_counts.items():
        print(f"{subclade}: {count} SNPs detected")

    print("\n\033[1mDetected SNPs by Haplogroup:\033[0m")
    for subclade, snps in snps_detected_in_haplogroups.items():
        print(f"{subclade}: {', '.join(snps)}")

def main():
    dna_file_path = input("Enter the path to your DNA file: ")
    if not os.path.exists(dna_file_path):
        print("Error: The specified file does not exist.")
        return

    detected_snps = parse_dna_file(dna_file_path)
    broad_haplogroup_counts, specific_haplogroup_counts, snps_detected_in_haplogroups = calculate_haplogroups(detected_snps)
    predicted_haplogroup = predict_haplogroup(broad_haplogroup_counts)
    print_results(predicted_haplogroup, broad_haplogroup_counts, specific_haplogroup_counts, snps_detected_in_haplogroups)

if __name__ == "__main__":
    main()
