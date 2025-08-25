import os
from collections import defaultdict

# Updated dictionary with rsIDs (SNPs) and their associated genotypes (allele pairs) for each haplogroup
# Sources: ISOGG Y-SNP tree, Underhill et al. 2015, Balaresque et al. 2010, Cruciani et al. 2007, Semino et al. 2004
haplogroup_snps = {
    # Broad Haplogroups and their subclades
    "R1a": {
        "R1a": {
            # R-M420 defining SNPs - Source: Underhill et al. 2015, ISOGG tree
            "rs17222573": ["A;A", "A;G", "G;G"],  
            "rs17307677": ["C;C", "C;T", "T;T"],
            "rs17306692": ["A;A", "A;C", "C;C"],
            "rs17250535": ["A;T", "T;T"],
        },
        "R1a1": {
            # R-M417 subclade - Source: Underhill et al. 2015
            "rs17307105": ["A;A", "A;G", "G;G"],
            "rs9786587": ["G;G", "G;T", "T;T"],
            "rs17316227": ["A;A", "A;G", "G;G"],
            "rs2534636": ["A;A"],
            "rs8179029": ["G;G"],  # Added proxy SNP for R1a-M417 (Underhill et al. 2015)
        },
        "R1a1a": {
            # R-M198/P1 subclade - Source: ISOGG 2019 tree
            "rs2020857": ["T;T"],
            "rs1722146": ["T;T"],
            "rs17315926": ["C;C","C;T","T;T"],
            "rs17221601": ["A;A","A;T","T;T"],
            "rs16981293": ["A;A","A;T"],  
        },
        "R1a1a1": {
            # R-M417 downstream - Source: ISOGG tree
            "rs17316771": ["A;A","A;G","G;G"]
        },
        "R1a1a1a": {
            # R-M56 subclade - Source: ISOGG tree
            "rs2032622": ["T;T"]
        },
        "R1a1a1b": {
            # R-Z282 (European R1a) - Source: Underhill et al. 2015
            "rs4988371": ["T;T"],  # Z282 equivalent SNP
            "rs2032659": ["A;A"],
        },
        "R1a1a1c": {
            # R-Z93 (Asian R1a) - Source: Underhill et al. 2015
            "rs2032626": ["G;G"],
            "rs2032644": ["C;C"],
            "rs2032655": ["G;G"]
        },
    },
    "R1b": {
        "R1b": {
            # R-M343 defining SNPs - Source: ISOGG tree, Myres et al. 2011
            "rs9786184": ["A;A"],
            "rs150173": ["A;A"],
            
        },
        "R1b1a2": {
            # R-M269 (European R1b)
            "rs9786882": ["A;A","A;G","G;G"],
            "rs9786153": ["T;T"],
            "rs877756": ["C;C","C;T","T;T"],
            "rs2058276": ["A;A", "A;G", "G;G"],
            "rs9786140": ["A;A","A;G","G;G"],  # Reinforced for M269 (Balaresque et al. 2010)
        },
        "R1b1a2a1a1": {
            # R-P312 (Western European) - Source: ISOGG tree
            "rs9786076": ["C;C","C;T","T;T"],
            "rs13304168": ["C;C","C;T","T;T"],
            "rs2082033": ["C;C","C;T","T;T"],
            "rs9786283": ["A;A","A;C","C;C"],
            "rs9785659": ["A;A","A;G","G;G"],
        },
        "R1b1b2a1a2": {
            "rs34276300": ["A;A"],
        },
        "R1b1b2g": {
            "rs16981293": ["C;C"],              
        },            
    },

    # I1 - Germanic/Nordic haplogroup
    "I1": {
        "I1": {
            # I-M253 defining SNPs - Source: ISOGG tree, Rootsi et al. 2004
            "rs17307252": ["A;A", "A;G", "G;G"],
            "rs871626": ["A;A", "A;G", "G;G"],
            "rs17221531": ["C;C","C;T","T;T"],
            "rs9341296": ["T;T"],
            "rs13447354": ["A;A"],
            "rs9341274": ["G;G"],
            "rs2032637": ["G;G"],
            "rs34626372": ["C;C"],
            "rs3912": ["T;T"],
        },
    },

    # I2 - Balkan haplogroup
    "I2": {
        "I2": {
            # I-M438 defining SNPs - Source: ISOGG tree, Rootsi et al. 2004
            "rs35547782": ["C;C","C;T","T;T"],
            "rs17307294": ["G;G"],
        },
        "I2a1b": {
            # I-M423 (Dinaric) - Source: ISOGG tree, Battaglia et al. 2009
            "rs2032666": ["T;T"],
        },
    },

    # G - Caucasian haplogroup
    "G": {
        "G": {
            # G-M201 defining SNPs - Source: ISOGG tree
            "rs2032656": ["T;T"],
            "rs4988323": ["C;C"],
        },
    },

    # N - Siberian/Uralic haplogroup
    "N": {
        "N": {
            # N-M231 defining SNPs - Source: ISOGG tree, Rootsi et al. 2007
            "rs4988327": ["G;G"],
        },
        "N1c": {
            # N-M178 (Uralic) - Source: ISOGG tree
            "rs2032688": ["G;G"],
        },
    },
}

def parse_dna_file(file_path):
    """Parse DNA file in various formats (23andMe, AncestryDNA, etc.)"""
    detected_snps = {}
    with open(file_path, 'r') as f:
        for line in f:
            # Skip header lines and comments
            if line.startswith("#") or "rsid" in line.lower() or line.strip() == "":
                continue
            
            parts = line.strip().split()
            if len(parts) >= 4:
                rsid = parts[0]
                # Handle different file formats
                if len(parts) == 4:  # Standard format: rsid chr pos genotype
                    genotype = parts[3]
                    if len(genotype) == 2:
                        allele1, allele2 = genotype[0], genotype[1]
                    else:
                        continue
                elif len(parts) >= 5:  # Format: rsid chr pos allele1 allele2
                    allele1, allele2 = parts[3], parts[4]
                else:
                    continue
                
                # Normalize genotype format
                if allele1 == allele2:
                    genotype_formatted = f"{allele1};{allele2}"
                else:
                    # Sort alleles alphabetically for consistency
                    sorted_alleles = sorted([allele1, allele2])
                    genotype_formatted = f"{sorted_alleles[0]};{sorted_alleles[1]}"
                
                detected_snps[rsid] = genotype_formatted
    
    return detected_snps

def calculate_haplogroups(detected_snps):
    """Calculate haplogroup scores based on detected SNPs"""
    broad_haplogroup_scores = defaultdict(int)
    specific_haplogroup_scores = defaultdict(int)
    snps_detected_in_haplogroups = defaultdict(list)

    # Define SNP weights based on phylogenetic importance
    # High weight: Major branch-defining SNPs
    high_weight_snps = {
        # Major haplogroup defining SNPs
        "rs9786184": 20,  # R1b-M343
        "rs17307294": 20, # I2-M438
        "rs9341274": 20,  # I1-L22
        "rs2032637": 20,  # I1-P109
        "rs34626372": 20, # I1-DF29
        "rs3912": 20,     # I1a
        "rs34276300": 20, # R1b-U106
        "rs16981293": 20, # R1b-L21
        "rs2032622": 20,  # R1a-M56
    }
    
    # Medium weight: Subclade defining SNPs
    medium_weight_snps = {
        "rs2020857": 15, "rs2032626": 15, "rs2032644": 15, 
        "rs2032655": 15, "rs150173": 15, "rs9785702": 15, 
        "rs9786153": 15, "rs17222279": 15, "rs1800865": 15, 
        "rs1236440": 15, "rs2032615": 15, "rs20321": 15, 
        "rs9341296": 15, "rs13447354": 15, "rs17316597": 15,
        "rs35474563": 15, "rs2032658": 15, "rs13447527": 15,
        "rs2032666": 15, "rs2032652": 15, "rs2032663": 15,
    }
    
    # Low weight: General supporting SNPs
    low_weight_snps = defaultdict(lambda: 8)

    for broad_haplogroup, subclades in haplogroup_snps.items():
        for subclade, snps_list in subclades.items():
            subclade_score = 0
            for snp, genotypes in snps_list.items():
                if snp in detected_snps:
                    genotype = detected_snps[snp]
                    # Check if detected genotype matches expected genotypes
                    if not genotypes or genotype in genotypes:
                        weight = high_weight_snps.get(snp, medium_weight_snps.get(snp, low_weight_snps[snp]))
                        subclade_score += weight
                        broad_haplogroup_scores[broad_haplogroup] += weight
                        snps_detected_in_haplogroups[subclade].append(snp)
            
            specific_haplogroup_scores[subclade] = subclade_score if subclade_score > 0 else 'not calculated'

    return broad_haplogroup_scores, specific_haplogroup_scores, snps_detected_in_haplogroups

def predict_haplogroup(broad_haplogroup_scores):
    """Predict the most likely haplogroup based on scores"""
    if not broad_haplogroup_scores:
        return "Unknown"
    
    # Get the haplogroup with highest score
    predicted_haplogroup = max(broad_haplogroup_scores, key=broad_haplogroup_scores.get)
    max_score = broad_haplogroup_scores[predicted_haplogroup]
    
    # Check if there's a clear winner (score significantly higher than others)
    if max_score < 15:  # Minimum threshold for confidence
        return "Unknown (low confidence)"
    
    # Check for ties or very close scores
    close_competitors = [hg for hg, score in broad_haplogroup_scores.items() 
                        if score >= max_score * 0.8 and hg != predicted_haplogroup]
    
    if close_competitors:
        return f"{predicted_haplogroup} (possible: {', '.join(close_competitors)})"
    
    return predicted_haplogroup

def get_most_specific_subclade(specific_haplogroup_scores, predicted_broad_haplogroup):
    """Find the most specific subclade for the predicted haplogroup"""
    if predicted_broad_haplogroup == "Unknown":
        return "Unknown"
    
    # Extract just the base haplogroup name (remove confidence notes)
    base_haplogroup = predicted_broad_haplogroup.split(' ')[0]
    
    # Find all subclades for this haplogroup
    relevant_subclades = {}
    for subclade, score in specific_haplogroup_scores.items():
        if (subclade.startswith(base_haplogroup) and 
            score != 'not calculated' and score > 0):
            relevant_subclades[subclade] = score
    
    if not relevant_subclades:
        return base_haplogroup
    
    # Return the most specific (longest name) subclade with highest score
    best_subclade = max(relevant_subclades.items(), 
                       key=lambda x: (x[1], len(x[0])))
    
    return best_subclade[0]

def print_results(predicted_haplogroup, broad_haplogroup_scores, specific_haplogroup_scores, snps_detected_in_haplogroups):
    """Print formatted results"""
    
    print(f"\033[1mPredicted Y Haplogroup:\033[0m {predicted_haplogroup}")
    
    # Show most specific subclade
    most_specific = get_most_specific_subclade(specific_haplogroup_scores, predicted_haplogroup)
    if most_specific != predicted_haplogroup.split(' ')[0]:
        print(f"\033[1mMost Specific Subclade:\033[0m {most_specific}")
    
    print(f"\n\033[1mBroad Haplogroup Scores:\033[0m")
    sorted_broad = sorted(broad_haplogroup_scores.items(), key=lambda x: x[1], reverse=True)
    for haplogroup, score in sorted_broad:
        print(f"  {haplogroup}: {score} points")

    print(f"\n\033[1mSubclade Detection Scores:\033[0m")
    # Group by broad haplogroup and show only relevant ones
    for broad_hg in sorted(broad_haplogroup_scores.keys()):
        relevant_subclades = [(sc, score) for sc, score in specific_haplogroup_scores.items() 
                             if sc.startswith(broad_hg) and score != 'not calculated' and score > 0]
        if relevant_subclades:
            print(f"  {broad_hg}:")
            for subclade, score in sorted(relevant_subclades, key=lambda x: x[1], reverse=True):
                print(f"    {subclade}: {score} points")

    print(f"\n\033[1mDetected SNPs by Subclade:\033[0m")
    for subclade, snps in snps_detected_in_haplogroups.items():
        if snps:  # Only show subclades with detected SNPs
            print(f"  {subclade}: {', '.join(snps)}")
    
    # Add confidence note
    max_score = max(broad_haplogroup_scores.values()) if broad_haplogroup_scores else 0
    confidence = "High" if max_score >= 40 else "Medium" if max_score >= 20 else "Low"
    print(f"\n\033[1mConfidence Level:\033[0m {confidence}")

def main():
    """Main function"""
    print("Y-Chromosome Haplogroup Predictor")
    print("=" * 40)
    
    dna_file_path = input("Enter the path to your DNA file: ").strip()
    
    if not os.path.exists(dna_file_path):
        print("Error: The specified file does not exist.")
        return

    try:
        print("\nParsing DNA file...")
        detected_snps = parse_dna_file(dna_file_path)
        
        if not detected_snps:
            print("Error: No valid SNP data found in the file.")
            return
            
        print(f"Found {len(detected_snps)} SNPs in total.")
        
        print("Calculating haplogroup scores...")
        broad_haplogroup_scores, specific_haplogroup_scores, snps_detected_in_haplogroups = calculate_haplogroups(detected_snps)
        
        predicted_haplogroup = predict_haplogroup(broad_haplogroup_scores)
        
        print("\n" + "=" * 60)
        print_results(predicted_haplogroup, broad_haplogroup_scores, specific_haplogroup_scores, snps_detected_in_haplogroups)
        
    except Exception as e:
        print(f"Error processing file: {str(e)}")
        print("Please check that your file is in the correct format.")

if __name__ == "__main__":
    main()
