# config_and_print.py

path = '../projects/GSE63525/GM12878/'
resolutions = [1000000]
chromosomes = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22']
data_types = ['oe']
genomeID = "hg19"
hic_url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525%5FGM12878%5Finsitu%5Fprimary%2Breplicate%5Fcombined%5F30%2Ehic"
bigwig_file = '/home/dwk681/workspace/hypermatrix_test/hypermatrix/projects/softwarefiles/ENCFF000EHJ_hg19_wgEncodeCrgMapabilityAlign36mer.bigWig'
mappability_threshold = 0.5

print(f"path='{path}'")
print(f"resolutions=({', '.join(map(str, resolutions))})")
print(f"chromosomes=({' '.join(chromosomes)})")
print(f"data_types=({' '.join(data_types)})")
print(f"genomeID='{genomeID}'")
print(f"hic_url='{hic_url}'")
print(f"bigwig_file='{bigwig_file}'")
print(f"mappability_threshold={mappability_threshold}")

