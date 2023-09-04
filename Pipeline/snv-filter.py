#snv annotation
python snpEff.py

#!snpEff.py:
import csv
import subprocess
i_j_mapping = {}
with open('info.csv', newline='') as csvfile:
    reader = csv.reader(csvfile)
    for row in reader:
        i = row[0]
        j = row[1]
        if j not in i_j_mapping:
            i_j_mapping[j] = []
        i_j_mapping[j].append(i)

for j, i_list in i_j_mapping.items():
    i_values = ' '.join(i_list)
    cmd = f"java -jar snpEff.jar {i_values} /home/ec2-user/yjw/yjw/snv/vcf/{j}.vcf > /efs/yjw/snpeffresult/{j}.vcf"
    subprocess.run(cmd, shell=True)

#filter synonymous SNVs

import os
import glob

def parse_vcf_file(vcf_file):
    total_snps = 0
    missense_variant = 0
    synonymous_snps = 0

    with open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            info_field = fields[7]

            if 'ANN=' in info_field:
                annotations = info_field.split('ANN=')[1].split(',')
                for annotation in annotations:
                    if '|synonymous_variant|' in annotation:
                        synonymous_snps += 1
                    if '|missense_variant|' in annotation:
                        missense_variant += 1

            total_snps += 1

    return total_snps, synonymous_snps, missense_variant

def calculate_ratios(vcf_file):
    total_snps, synonymous_snps, missense_variant = parse_vcf_file(vcf_file)

    if total_snps == 0:
        return 0.0, 0.0, 0.0

    synonymous_ratio = synonymous_snps / total_snps
    missense_ratio = missense_variant / total_snps

    return synonymous_ratio, missense_ratio, total_snps

def main():
    input_folder = '/efs/yjw/snpeffresult' 
    output_file = '/efs/yjw/mutation_ratios.txt' 

    with open(output_file, 'w') as out:
        out.write("Sample\tSynonymous_Ratio\tMissense_Ratio\tTotal_SNPs\n")

        for vcf_file in glob.glob(os.path.join(input_folder, '*.vcf')):
            sample_name = os.path.basename(vcf_file).rsplit('.', 1)[0] #以第二个点进行分割
            synonymous_ratio, missense_ratio, total_snps = calculate_ratios(vcf_file)
            out.write(f"{sample_name}\t{synonymous_ratio:.4f}\t{missense_ratio:.4f}\t{total_snps}\n")

if __name__ == "__main__":
    main()
